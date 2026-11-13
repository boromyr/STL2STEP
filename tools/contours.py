"""contours.py MESH.stl [ORIG.step] [--png out.png cx cy cz dx dy dz size] [--save labels.npy]
The CAD's edges read from the mesh alone: every mesh edge scored, the scored ones closed into a network,
the triangles between them grouped into faces. With ORIG.step (the true CAD, aligned by ICP) the result is measured:
edges caught / false, and how many true faces come out whole (one predicted face = one CAD face)."""
import sys, os, math
import numpy as np
from collections import defaultdict, Counter
sys.path.insert(0, "D:/Users/PC/Desktop/refit")
import refit as R
R.Log.level = 100


def tri_feat(V):
    e = np.roll(V, -1, axis=0) - V
    L = np.linalg.norm(e, axis=1)
    a = 0.5 * float(np.linalg.norm(np.cross(V[1] - V[0], V[2] - V[0])))
    lmax = float(L.max())
    h = 2.0 * a / max(lmax, 1e-12)
    j = int(np.argmax(L))
    return float(np.median(L)), lmax / max(h, 1e-12), e[j] / max(lmax, 1e-12)


CREASE_MIN, CREASE_K = 1e9, 2.5
STEP_MIN, STEP_REL = 0.3, 0.3


class Mesh:
    def __init__(self, topo):
        self.topo = topo
        nF = topo.nF
        self.N = np.array([topo.norms[i] for i in range(nF)])
        self.size = np.zeros(nF); self.elong = np.zeros(nF); self.dir = np.zeros((nF, 3))
        for i in range(nF):
            self.size[i], self.elong[i], self.dir[i] = tri_feat(topo.verts[i])
        self.E = [k for k, fs in enumerate(topo.e_faces) if len(fs) == 2]
        self.ef = np.array([topo.e_faces[k] for k in self.E])
        self.ev = np.array([topo.e_verts[k] for k in self.E])
        a, b = self.ef[:, 0], self.ef[:, 1]
        self.ang = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", self.N[a], self.N[b]), -1, 1)))
        self.dsize = np.abs(np.log2(self.size[a] / self.size[b]))
        self.cosd = np.abs(np.einsum("ij,ij->i", self.dir[a], self.dir[b]))
        self.el_a, self.el_b = self.elong[a], self.elong[b]
        # how much each triangle turns from its OTHER neighbours (the largest of them, sharp ones apart)
        nb = defaultdict(list)
        for k in range(len(a)):
            nb[int(a[k])].append(k); nb[int(b[k])].append(k)
        def other(i, k):
            return max([self.ang[k2] for k2 in nb[i] if k2 != k and self.ang[k2] < 40.0] or [0.0])
        # the strip a triangle is half of (its one flat neighbour), and the angle at the strip's far side:
        # inside a face the strips turn by one step, so an edge whose angle isn't the far side's is where the step changes
        tv = [set(int(x) for x in topo.f_verts[i]) if hasattr(topo, "f_verts") else None for i in range(topo.nF)]
        def far(i, k):
            flat = [k2 for k2 in nb[i] if k2 != k and self.ang[k2] < 0.05]
            if len(flat) != 1:
                return -1.0
            p = int(a[flat[0]]) if int(b[flat[0]]) == i else int(b[flat[0]])
            uv = {int(self.ev[k, 0]), int(self.ev[k, 1])}
            c = [k2 for k2 in nb[p] if k2 != flat[0] and not ({int(self.ev[k2, 0]), int(self.ev[k2, 1])} & uv)]
            return float(self.ang[c[0]]) if len(c) == 1 else -1.0
        self.far_a = np.array([far(int(a[k]), k) for k in range(len(a))])
        self.far_b = np.array([far(int(b[k]), k) for k in range(len(a))])
        self.bend_a = np.array([other(int(a[k]), k) for k in range(len(a))])
        self.bend_b = np.array([other(int(b[k]), k) for k in range(len(a))])


def score(m, v=1):
    """strong / weak marks per interior mesh edge"""
    strong = (m.ang > 20.0) | ((m.dsize > 1.5) & (m.ang > 0.2))
    # a crease: the two sides turn from each other more than either turns from anything else beside it
    hi = np.maximum(m.bend_a, m.bend_b)
    strong |= (m.ang > CREASE_MIN) & (m.ang > CREASE_K * (hi + 0.1))
    def step(far):
        top = np.maximum(m.ang, far)
        return (far >= 0.0) & (far < 40.0) & (top > STEP_MIN) & (np.abs(m.ang - far) > STEP_REL * top)
    strong |= (m.ang < 20.0) & (step(m.far_a) | step(m.far_b)) & (m.ang > 0.05)
    both = (m.el_a >= 3.0) & (m.el_b >= 3.0)
    mixed = ((m.el_a >= 4.0) & (m.el_b < 2.0)) | ((m.el_b >= 4.0) & (m.el_a < 2.0))
    weak = (both & (m.cosd < 0.9)) | mixed | (m.dsize > 1.0) | (m.ang > 10.0)
    return strong, weak & ~strong


def close_network(m, strong, weak):
    """weak edges stay only where they are part of a network with no loose ends (pruned until none is left); strong ones stay"""
    keep = strong | weak
    deg = Counter()
    for k in np.nonzero(keep)[0]:
        deg[int(m.ev[k, 0])] += 1; deg[int(m.ev[k, 1])] += 1
    at = defaultdict(list)
    for k in np.nonzero(keep)[0]:
        at[int(m.ev[k, 0])].append(k); at[int(m.ev[k, 1])].append(k)
    stack = [k for k in np.nonzero(weak)[0] if deg[int(m.ev[k, 0])] == 1 or deg[int(m.ev[k, 1])] == 1]
    while stack:
        k = stack.pop()
        if not keep[k] or strong[k]:
            continue
        v0, v1 = int(m.ev[k, 0]), int(m.ev[k, 1])
        if deg[v0] != 1 and deg[v1] != 1:
            continue
        keep[k] = False
        for v in (v0, v1):
            deg[v] -= 1
            if deg[v] == 1:
                stack.extend(j for j in at[v] if keep[j] and not strong[j])
    return keep


def faces_of(m, cut):
    """triangles grouped across every interior edge that isn't cut"""
    nF = m.topo.nF
    par = list(range(nF))
    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for k in np.nonzero(~cut)[0]:
        a, b = find(int(m.ef[k, 0])), find(int(m.ef[k, 1]))
        if a != b:
            par[a] = b
    return np.array([find(i) for i in range(nF)])


def absorb_small(m, lab, cut, min_tris=3):
    """a face of one or two triangles goes to the neighbour it shares the flattest side with"""
    cnt = Counter(lab.tolist())
    nb = defaultdict(list)
    for k in range(len(m.ef)):
        a, b = int(m.ef[k, 0]), int(m.ef[k, 1])
        nb[a].append((b, k)); nb[b].append((a, k))
    lab = lab.copy()
    for i in np.argsort(m.size):
        if cnt[int(lab[i])] >= min_tris:
            continue
        best = min(((m.ang[k], j) for j, k in nb[i] if lab[j] != lab[i]), default=None)
        if best is not None and best[0] < 20.0:
            cnt[int(lab[i])] -= 1
            lab[i] = lab[best[1]]
            cnt[int(lab[i])] += 1
    return lab


def extract(m):
    strong, weak = score(m)
    cut = close_network(m, strong, weak)
    lab = faces_of(m, cut)
    lab = absorb_small(m, lab, cut)
    return lab, strong, weak, cut


def evaluate(m, lab, of, areas):
    a, b = m.ef[:, 0], m.ef[:, 1]
    true = of[a] != of[b]
    pred = lab[a] != lab[b]
    smooth = m.ang < 10.0
    out = dict(edges=len(a), cad=int(true.sum()), caught=float((pred & true).sum() / max(true.sum(), 1)), false=float((pred & ~true).sum() / max((~true).sum(), 1)),
               smooth_cad=int((true & smooth).sum()), smooth_caught=float((pred & true & smooth).sum() / max((true & smooth).sum(), 1)))
    # faces: a CAD face is whole if 95% of its area is in one predicted face, and that face is 95% its own
    A_tp = defaultdict(float); A_t = defaultdict(float); A_p = defaultdict(float)
    for i in range(len(lab)):
        A_tp[(int(of[i]), int(lab[i]))] += areas[i]; A_t[int(of[i])] += areas[i]; A_p[int(lab[i])] += areas[i]
    whole = split = merged = 0
    best = {}
    for (t, p), ar in A_tp.items():
        if ar > best.get(t, (0, None))[0]:
            best[t] = (ar, p)
    for t, (ar, p) in best.items():
        comp = ar / A_t[t]; pur = ar / A_p[p]
        if comp >= 0.95 and pur >= 0.95:
            whole += 1
        elif comp < 0.95:
            split += 1
        else:
            merged += 1
    out.update(cad_faces=len(A_t), pred_faces=len(A_p), whole=whole, split=split, merged=merged)
    return out


def surf_sig(shA):
    """one id per underlying surface: two CAD faces on the same plane / cylinder / cone / sphere / torus are one surface"""
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    sig, ids = [], {}
    for n, f in enumerate(R.explore(shA, R.TopAbs_FACE)):
        s_ = BRepAdaptor_Surface(R.td_Face(f)); t = str(s_.GetType()).split("_")[-1]
        r = lambda x, d=3: round(float(x), d) + 0.0
        try:
            if t == "Plane":
                p = s_.Plane(); n_ = p.Axis().Direction(); o = p.Location()
                v = np.array([n_.X(), n_.Y(), n_.Z()]); d = float(v @ np.array([o.X(), o.Y(), o.Z()]))
                if v[np.argmax(np.abs(v))] < 0: v, d = -v, -d
                key = (t, r(v[0], 4), r(v[1], 4), r(v[2], 4), r(d))
            elif t in ("Cylinder", "Cone", "Sphere", "Torus"):
                g = getattr(s_, t)(); ax = g.Position() if t != "Sphere" else None
                if t == "Sphere":
                    o = g.Location(); key = (t, r(o.X()), r(o.Y()), r(o.Z()), r(g.Radius()))
                else:
                    dd = ax.Direction(); o = ax.Location(); v = np.array([dd.X(), dd.Y(), dd.Z()]); p0 = np.array([o.X(), o.Y(), o.Z()])
                    if v[np.argmax(np.abs(v))] < 0: v = -v
                    q = p0 - (p0 @ v) * v if t == "Cylinder" else p0
                    rad = (g.Radius(),) if t == "Cylinder" else ((g.RefRadius(), g.SemiAngle()) if t == "Cone" else (g.MajorRadius(), g.MinorRadius()))
                    key = (t, r(v[0], 4), r(v[1], 4), r(v[2], 4), r(q[0]), r(q[1]), r(q[2])) + tuple(r(x) for x in rad)
            else:
                key = ("face", n)
        except Exception:
            key = ("face", n)
        sig.append(ids.setdefault(key, len(ids)))
    return np.array(sig)


def ground_truth(step, shB, topo, m):
    """per triangle: the CAD face it came from - faces on one surface that touch counted as one"""
    from icp import align, closest_on_tris
    shA = R.read_input(step)[0]
    toA, med, mx, tA, fA, tipA, arA = align(shA, shB)
    cent = np.array([topo.cents[i] for i in range(topo.nF)])
    dd, _, jj = closest_on_tris(toA(cent), tA)
    sg = surf_sig(shA)[fA[jj]]
    par = list(range(topo.nF))
    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for k in range(len(m.ef)):
        a_, b_ = int(m.ef[k, 0]), int(m.ef[k, 1])
        if sg[a_] == sg[b_]:
            ra, rb = find(a_), find(b_)
            if ra != rb: par[ra] = rb
    return np.array([find(i) for i in range(topo.nF)])


if __name__ == "__main__":
    args = sys.argv[1:]
    shB = R.read_input(args[0])[0]
    topo = R.Topo(shB)
    m = Mesh(topo)
    lab, strong, weak, cut = extract(m)
    areas = np.array([topo.areas[i] for i in range(topo.nF)])
    print(f"{os.path.basename(args[0])}: {topo.nF} triangles, {len(m.E)} edges; strong {int(strong.sum())}, weak {int(weak.sum())}, kept {int(cut.sum())}; faces {len(set(lab.tolist()))}")
    if len(args) > 1 and args[1].lower().endswith((".step", ".stp")):
        of = ground_truth(args[1], shB, topo, m)
        r = evaluate(m, lab, of, areas)
        print(f"  CAD edges caught {100*r['caught']:.1f}%  false {100*r['false']:.2f}%  | smooth CAD edges ({r['smooth_cad']}) caught {100*r['smooth_caught']:.1f}%")
        print(f"  CAD faces {r['cad_faces']}: whole {r['whole']}, split {r['split']}, merged with another {r['merged']}  (predicted faces {r['pred_faces']})")
    if "--save" in args:
        cents = np.array([topo.cents[i] for i in range(topo.nF)])
        np.save(args[args.index("--save") + 1], np.column_stack([cents, lab]))
