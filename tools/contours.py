"""
contours.py - the CAD's edges read from the mesh alone.

Every interior edge of the mesh is scored from the tessellation itself - no surface is fitted:
  strong  a crease (over 20 degrees, or turning more than either side turns from anything else
          beside it), a jump in facet size, a change of the turning step across a strip;
  weak    a change of tessellation pattern (long facets meeting at an angle, long against round,
          a size change, over 10 degrees).
The weak ones stay only where they close a network with no loose ends; the triangles between the
kept edges are grouped into faces (one or two triangles go to the neighbour they share the
flattest side with). What is left are the mesh's contours: the lines where the exporter's
tessellation says one CAD face ends and another begins - also between tangent faces, where no
angle shows it.

  python tools/contours.py MESH.stl                       the network and its faces
  python tools/contours.py MESH.stl --cad ORIG.step       measured against the true CAD (aligned by
                                                          ICP): edges caught / false, faces whole /
                                                          split / merged with another
  python tools/contours.py MESH.stl --step OUT.step       does a conversion's every edge lie on a mesh
                                                          contour? An edge off them is a face boundary
                                                          in the wrong place (a row of facets given to
                                                          the neighbour: test17's arches, 13/11); a
                                                          strong contour no edge covers is a face
                                                          across a crease. Worst ones listed with
                                                          where they are and the faces they part.
  --png FILE [--at X,Y,Z] [--dir DX,DY,DZ] [--span MM]    picture: the contours red (strong thick,
                                                          weak thin), and with --step the output's
                                                          faces by kind (plane grey, cylinder blue,
                                                          cone orange, sphere violet, torus red,
                                                          B-spline green), its edges black
  --near X,Y,Z,R                                          reports restricted to a ball
  --on MM                                                 how close an edge must lie to a contour to
                                                          be on it (default 5 x refit's automatic
                                                          tolerance for the mesh's size)
  --save FILE.npy                                         per triangle: centroid x, y, z, face label
  --json FILE                                             everything measured, for scripts
  (the second positional argument, if a .step, is --cad: the old usage)

As a module:
  import contours as C
  c = C.Contours.from_file("part.stl")      # or C.Contours(topo) on a refit.Topo
  c.labels                                  # face label per triangle
  c.polylines()                             # the contours as chains of mesh vertices
  c.faces()                                 # per face: triangles, area, flat or curved
  c.follow(shape)                           # a B-Rep's edges against the contours
"""
import sys, os, math, json, argparse
import numpy as np
from collections import defaultdict, Counter

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
import refit as R

R.Log.level = 100


# --- the scores ---------------------------------------------------------------------------------

def tri_feat(V):
    """median side, elongation (longest side over its height), direction of the longest side"""
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
    """per-triangle and per-edge features of a refit.Topo (triangles only)"""

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


def score(m):
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
    """(labels, strong, weak, kept) - the whole pipeline on a Mesh"""
    strong, weak = score(m)
    cut = close_network(m, strong, weak)
    lab = faces_of(m, cut)
    lab = absorb_small(m, lab, cut)
    return lab, strong, weak, cut


# --- the result -----------------------------------------------------------------------------------

def _seg_dist(P, A, B):
    """distance of each point P[i] from the segment A[i]-B[i]"""
    ab = B - A
    L2 = np.maximum(np.einsum("ij,ij->i", ab, ab), 1e-300)
    t = np.clip(np.einsum("ij,ij->i", P - A, ab) / L2, 0.0, 1.0)
    return np.linalg.norm(P - (A + t[:, None] * ab), axis=1)


class Contours:
    """the contours of a triangle mesh (a refit.Topo of it)"""

    def __init__(self, topo):
        self.topo = topo
        self.m = Mesh(topo)
        self.labels, self.strong, self.weak, self.kept = extract(self.m)
        m = self.m
        # the contours are where the faces part (after the small ones went to their neighbours)
        self.border = self.labels[m.ef[:, 0]] != self.labels[m.ef[:, 1]]
        self.P0 = topo.vpos[m.ev[:, 0]]
        self.P1 = topo.vpos[m.ev[:, 1]]
        self.areas = np.array([topo.areas[i] for i in range(topo.nF)])
        self.cents = np.array([topo.cents[i] for i in range(topo.nF)])
        V = topo.vpos
        self.diag = float(np.linalg.norm(V.max(axis=0) - V.min(axis=0)))
        self._kd = None

    @classmethod
    def from_file(cls, path):
        sh = R.read_input(path)[0]
        return cls(R.Topo(sh))

    # (what it found)
    def summary(self):
        b = self.border
        L = np.linalg.norm(self.P1 - self.P0, axis=1)
        return dict(triangles=int(self.topo.nF), edges=int(len(self.m.E)), strong=int(self.strong.sum()), weak=int(self.weak.sum()),
                    kept=int(self.kept.sum()), faces=int(len(set(self.labels.tolist()))), contour_edges=int(b.sum()),
                    contour_length=float(L[b].sum()), strong_length=float(L[b & self.strong].sum()))

    def polylines(self, near=None):
        """the contours as chains of mesh edges between junctions (a vertex where other than two contour edges meet):
        [{"verts": [...], "points": array, "length", "strong": fraction of the length on strong edges,
          "angle": the mean dihedral angle (degrees), "closed"}]"""
        m = self.m
        ks = np.nonzero(self.border)[0]
        if near is not None:
            c, r = np.asarray(near[:3], float), float(near[3])
            mid = 0.5 * (self.P0[ks] + self.P1[ks])
            ks = ks[np.linalg.norm(mid - c, axis=1) <= r]
        at = defaultdict(list)
        for k in ks:
            at[int(m.ev[k, 0])].append(int(k)); at[int(m.ev[k, 1])].append(int(k))
        used = set()
        out = []

        def walk(v, k):
            vs, es = [v], []
            while True:
                used.add(k); es.append(k)
                a, b = int(m.ev[k, 0]), int(m.ev[k, 1])
                v = b if a == v else a
                vs.append(v)
                if len(at[v]) != 2:
                    return vs, es, False
                nxt = [k2 for k2 in at[v] if k2 != k]
                if not nxt or nxt[0] in used:
                    return vs, es, vs[0] == vs[-1]
                k = nxt[0]

        starts = [v for v in at if len(at[v]) != 2]
        for v in starts:
            for k in at[v]:
                if k not in used:
                    out.append(walk(v, k))
        for k in ks:  # rings with no junction
            if int(k) not in used:
                out.append(walk(int(m.ev[k, 0]), int(k)))
        res = []
        for vs, es, closed in out:
            P = self.topo.vpos[vs]
            Ls = np.linalg.norm(self.P1[es] - self.P0[es], axis=1)
            res.append(dict(verts=vs, points=P, length=float(Ls.sum()), closed=bool(closed),
                            strong=float(Ls[self.strong[es]].sum() / max(Ls.sum(), 1e-300)),
                            angle=float(np.average(m.ang[es], weights=Ls) if Ls.sum() > 0 else 0.0)))
        return res

    def faces(self):
        """per face label: triangles, area, centroid, normal spread (degrees: under 1, flat)"""
        out = {}
        groups = defaultdict(list)
        for i, l in enumerate(self.labels.tolist()):
            groups[l].append(i)
        for l, ix in groups.items():
            ix = np.array(ix)
            w = self.areas[ix]
            n = (self.m.N[ix] * w[:, None]).sum(axis=0)
            n /= max(np.linalg.norm(n), 1e-300)
            spread = float(np.degrees(np.arccos(np.clip(self.m.N[ix] @ n, -1, 1))).max())
            out[int(l)] = dict(triangles=int(len(ix)), area=float(w.sum()), centroid=((self.cents[ix] * w[:, None]).sum(axis=0) / max(w.sum(), 1e-300)).tolist(),
                               spread=spread, flat=spread < 1.0)
        return out

    # (against a B-Rep)
    def _dense(self, step):
        """the contour segments sampled every `step` mm: points, and the segment each came from"""
        ks = np.nonzero(self.border)[0]
        A, B = self.P0[ks], self.P1[ks]
        L = np.linalg.norm(B - A, axis=1)
        n = np.maximum(1, np.ceil(L / step).astype(int))
        seg = np.repeat(ks, n + 1)
        t = np.concatenate([np.linspace(0.0, 1.0, k + 1) for k in n])
        P = np.repeat(A, n + 1, axis=0) + t[:, None] * np.repeat(B - A, n + 1, axis=0)
        return P, seg

    def contour_dist(self, P, step, with_seg=False):
        """the distance of each point from the nearest contour segment (and, with_seg, that segment's index)"""
        from scipy.spatial import cKDTree
        if self._kd is None or self._kd[0] != step:
            D, seg = self._dense(step)
            self._kd = (step, cKDTree(D), seg)
        _, kd, seg = self._kd
        if not len(P):
            return (np.zeros(0), np.zeros(0, dtype=int)) if with_seg else np.zeros(0)
        _, j = kd.query(P, k=min(4, kd.n))
        j = np.atleast_2d(j)
        if j.shape[0] != len(P):
            j = j.T
        d = np.full(len(P), np.inf)
        best = np.zeros(len(P), dtype=int)
        for c in range(j.shape[1]):
            k = seg[j[:, c]]
            dk = _seg_dist(P, self.P0[k], self.P1[k])
            better = dk < d
            d = np.where(better, dk, d)
            best = np.where(better, k, best)
        return (d, best) if with_seg else d

    def follow(self, shape, on=None, near=None, n_worst=15):
        """
        Every edge of the B-Rep sampled and measured from the contours (on: within `on` mm), and every strong contour
        edge from the B-Rep's edges. Returns a dict: per edge (length, max distance, fraction on, midpoint, the kinds of
        the faces it parts), the totals, the worst edges, the strong contours no edge covers.
        """
        from scipy.spatial import cKDTree
        from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
        on = float(on) if on else 5.0 * R.auto_tol(self.diag)
        step = 0.5 * on
        emap = R.edge_face_map(shape)
        fidx = {}
        for n_, f in enumerate(R.explore(shape, R.TopAbs_FACE)):
            fidx.setdefault(f.__hash__(), []).append((f, n_ + 1))

        def face_no(f):
            for f2, n_ in fidx.get(f.__hash__(), ()):
                if f2.IsSame(f):
                    return n_
            return 0

        kinds = {}

        def kind(f):
            n_ = face_no(f)
            if n_ not in kinds:
                kinds[n_] = str(BRepAdaptor_Surface(R.td_Face(f)).GetType()).split("_")[-1].replace("Surface", "")
            return kinds[n_]

        c_near = None if near is None else (np.asarray(near[:3], float), float(near[3]))
        edges, samples = [], []
        for i in range(1, R._size(emap) + 1):
            e = R.td_Edge(emap.FindKey(i))
            if R.bt_Degenerated(e):
                continue
            try:
                cv = BRepAdaptor_Curve(e)
                t0, t1 = cv.FirstParameter(), cv.LastParameter()
                P = np.array([[q.X(), q.Y(), q.Z()] for q in (cv.Value(float(t)) for t in np.linspace(t0, t1, 9))])
                L = float(np.linalg.norm(np.diff(P, axis=0), axis=1).sum())
                n = int(min(400, max(9, math.ceil(L / step) + 1)))
                if n > 9:
                    P = np.array([[q.X(), q.Y(), q.Z()] for q in (cv.Value(float(t)) for t in np.linspace(t0, t1, n))])
                ctype = str(cv.GetType()).split("_")[-1]
            except Exception:
                continue
            if c_near is not None and np.linalg.norm(P - c_near[0], axis=1).min() > c_near[1]:
                continue
            fl = emap.FindFromIndex(i)
            fs = []
            it = R.TopTools_ListIteratorOfListOfShape(fl) if hasattr(R, "TopTools_ListIteratorOfListOfShape") else None
            try:
                fs = [R.td_Face(x) for x in fl]
            except Exception:
                fs = []
                if it is not None:
                    while it.More():
                        fs.append(R.td_Face(it.Value())); it.Next()
            seams = len(fs) == 2 and fs[0].IsSame(fs[1])
            edges.append(dict(curve=ctype, length=L, faces=[face_no(f) for f in fs], kinds=[kind(f) for f in fs], seam=bool(seams), mid=P[len(P) // 2].tolist()))
            samples.append(P)
        # (off a contour by less than a facet: a contour runs alongside and the edge isn't on it - a boundary moved by a
        # row of facets; further: no contour there at all, which may be the detector's miss - long facets meeting at a
        # tangent edge often show nothing)
        facet = float(np.median(np.linalg.norm(self.P1 - self.P0, axis=1)))
        for E, P in zip(edges, samples):
            # ⚠️ THE CONTOUR'S VERTICES ON THE EDGE, NOT ITS CHORDS. An exact curve runs through the mesh's vertices and
            # between two of them stands off their chord by its sag (test17's r 21.5 circles, 0.047 off the 47-gon's
            # sides; a thread's helix more): measured from the chords, every exact edge looked moved. So each sample is
            # matched to its nearest contour segment, and that segment's ends are measured from the edge - ends beyond
            # the edge's own ends left out (the next contour along).
            d, k = self.contour_dist(P, step, with_seg=True)
            Pd = P if len(P) >= 2 else np.vstack([P, P])
            A_, B_ = Pd[:-1], Pd[1:]

            def to_edge(Q):
                """distance of each point from the sampled edge, and whether its nearest point is one of the edge's ends"""
                dd = np.full(len(Q), np.inf)
                tt = np.zeros(len(Q))
                for s_ in range(len(A_)):
                    ab = B_[s_] - A_[s_]
                    L2 = max(float(ab @ ab), 1e-300)
                    t = np.clip(((Q - A_[s_]) @ ab) / L2, 0.0, 1.0)
                    ds = np.linalg.norm(Q - (A_[s_] + t[:, None] * ab), axis=1)
                    better = ds < dd
                    dd = np.where(better, ds, dd)
                    tt = np.where(better, s_ + t, tt)
                end = (tt <= 1e-9) | (tt >= len(A_) - 1e-9)
                return dd, end

            ku = np.unique(k)
            e0, x0 = to_edge(self.P0[ku])
            e1, x1 = to_edge(self.P1[ku])
            # (an end past the edge's own ends doesn't count; one that is a contour going off sideways - its other end
            # on the edge, it is far - is another contour meeting the edge, not this one's)
            e0 = np.where(x0 & (e0 > on), 0.0, e0)
            e1 = np.where(x1 & (e1 > on), 0.0, e1)
            side = (np.minimum(e0, e1) <= on) & (np.maximum(e0, e1) > on)
            off_k = dict(zip(ku.tolist(), np.where(side, 0.0, np.maximum(e0, e1)).tolist()))
            dv = np.array([off_k[int(x)] for x in k])
            dv = np.where(d > facet, np.maximum(dv, d), dv)
            seg_l = np.linalg.norm(np.diff(P, axis=0), axis=1)
            ok = (dv[:-1] <= on) & (dv[1:] <= on)
            d = dv
            E["max"] = float(d.max())
            E["on"] = float(seg_l[ok].sum() / max(seg_l.sum(), 1e-300))
            E["class"] = "on" if E["on"] >= 0.999 else ("moved" if E["max"] <= facet else "unseen")
        real = [E for E in edges if not E["seam"]]
        L_all = sum(E["length"] for E in real)
        L_on = sum(E["length"] * E["on"] for E in real)
        L_moved = sum(E["length"] * (1.0 - E["on"]) for E in real if E["class"] == "moved")
        L_unseen = sum(E["length"] * (1.0 - E["on"]) for E in real if E["class"] == "unseen")
        # (the moved ones first, the longest stretch off first; a closed curve's seam is apart)
        worst = sorted([E for E in real if E["on"] < 0.999], key=lambda E: (E["class"] != "moved", -(1.0 - E["on"]) * E["length"]))[:n_worst]
        # the strong contours no edge of the B-Rep covers
        ks = np.nonzero(self.border & self.strong)[0]
        if c_near is not None:
            mid = 0.5 * (self.P0[ks] + self.P1[ks])
            ks = ks[np.linalg.norm(mid - c_near[0], axis=1) <= c_near[1]]
        unc_len = 0.0
        unc = []
        if len(ks) and samples:
            kd = cKDTree(np.vstack([np.vstack([P[:-1], 0.5 * (P[:-1] + P[1:])]) for P in samples]))
            for k in ks:
                Q = self.P0[k] + np.linspace(0.0, 1.0, 5)[:, None] * (self.P1[k] - self.P0[k])
                dd, _ = kd.query(Q)
                if float(dd.max()) > on + step:
                    Lk = float(np.linalg.norm(self.P1[k] - self.P0[k]))
                    unc_len += Lk
                    unc.append(dict(mid=(0.5 * (self.P0[k] + self.P1[k])).tolist(), length=Lk, angle=float(self.m.ang[k])))
        strong_len = float(np.linalg.norm(self.P1[ks] - self.P0[ks], axis=1).sum()) if len(ks) else 0.0
        unc.sort(key=lambda u: -u["length"])
        return dict(on=on, facet=facet, edges=len(real), seams=len(edges) - len(real), length=L_all, length_on=L_on,
                    length_moved=L_moved, length_unseen=L_unseen,
                    edges_moved=sum(1 for E in real if E["class"] == "moved"), edges_unseen=sum(1 for E in real if E["class"] == "unseen"),
                    edges_off=sum(1 for E in real if E["on"] < 0.999), worst=worst, all=edges,
                    strong_length=strong_len, strong_uncovered=unc_len, uncovered=unc[:n_worst])


# --- against the true CAD -------------------------------------------------------------------------

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


# --- the picture ----------------------------------------------------------------------------------

def render(c, png, at=None, view=(1.0, 1.0, 1.0), span=None, shape=None, size=(1400, 1100)):
    """the contours red (strong thick, weak thin) on the mesh - or on the B-Rep's faces by kind with its edges black"""
    import pyvista as pv
    topo = c.topo
    at = np.asarray(at, float) if at is not None else 0.5 * (topo.vpos.max(axis=0) + topo.vpos.min(axis=0))
    span = float(span) if span else 0.5 * c.diag
    d = np.asarray(view, float); d /= np.linalg.norm(d)
    pl = pv.Plotter(off_screen=True, window_size=size); pl.set_background("white")
    near = lambda P: np.linalg.norm(P - at, axis=1).min() <= 2.5 * span
    if shape is None:
        sel = np.nonzero(np.linalg.norm(c.cents - at, axis=1) <= 2.5 * span)[0]
        V = np.vstack([topo.verts[i] for i in sel]); F = np.hstack([[3, 3 * j, 3 * j + 1, 3 * j + 2] for j in range(len(sel))])
        pl.add_mesh(pv.PolyData(V, F), color="#c8c8c8", show_edges=True, edge_color="#999999", line_width=0.5)
    else:
        from OCP.BRepMesh import BRepMesh_IncrementalMesh
        from OCP.BRep import BRep_Tool
        from OCP.TopLoc import TopLoc_Location
        from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
        col = {"Plane": "#c8c8c8", "Cylinder": "#6aa0e0", "Cone": "#e0b050", "Sphere": "#d070d0", "Torus": "#e06070"}
        BRepMesh_IncrementalMesh(shape, max(1e-3, 1e-3 * span), False, 0.2, True)
        for f in R.explore(shape, R.TopAbs_FACE):
            loc = TopLoc_Location(); tr = R._st(BRep_Tool, "Triangulation")(R.td_Face(f), loc)
            if tr is None:
                continue
            V = np.array([[q.X(), q.Y(), q.Z()] for q in (tr.Node(i).Transformed(loc.Transformation()) for i in range(1, tr.NbNodes() + 1))])
            if not near(V):
                continue
            F = np.array([[3] + [x - 1 for x in tr.Triangle(i).Get()] for i in range(1, tr.NbTriangles() + 1)])
            k = str(BRepAdaptor_Surface(R.td_Face(f)).GetType()).split("_")[-1]
            pl.add_mesh(pv.PolyData(V, F.ravel()), color=col.get(k, "#70c070"), smooth_shading=True)
        for e in R.explore(shape, R.TopAbs_EDGE):
            try:
                a = BRepAdaptor_Curve(R.td_Edge(e))
                P = np.array([[q.X(), q.Y(), q.Z()] for q in (a.Value(float(t)) for t in np.linspace(a.FirstParameter(), a.LastParameter(), 40))])
            except Exception:
                continue
            if near(P):
                pl.add_mesh(pv.lines_from_points(P), color="black", line_width=2)
    for strong, w in ((False, 1.5), (True, 3.5)):
        ks = np.nonzero(c.border & (c.strong == strong))[0]
        if not len(ks):
            continue
        mid = 0.5 * (c.P0[ks] + c.P1[ks])
        ks = ks[np.linalg.norm(mid - at, axis=1) <= 2.5 * span]
        if not len(ks):
            continue
        P = np.empty((2 * len(ks), 3)); P[0::2] = c.P0[ks]; P[1::2] = c.P1[ks]
        L = pv.PolyData(P); L.lines = np.hstack([[2, 2 * i, 2 * i + 1] for i in range(len(ks))])
        pl.add_mesh(L, color="red", line_width=w)
    up = (0, 0, 1) if abs(d[2]) < 0.9 else (0, 1, 0)
    pl.camera_position = [tuple(at + d * 3 * span), tuple(at), up]
    pl.camera.view_angle = 2.0 * math.degrees(math.atan(1.0 / 3.0))
    pl.screenshot(png)


# --- command line ---------------------------------------------------------------------------------

def _vec(s, n):
    v = [float(x) for x in s.replace(";", ",").split(",")]
    if len(v) != n:
        raise argparse.ArgumentTypeError(f"{n} numbers, comma-separated")
    return v


def main(argv=None):
    ap = argparse.ArgumentParser(description="The CAD's edges read from the mesh alone (see the docstring).", formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mesh", help=".stl, or a .step born from a mesh")
    ap.add_argument("orig", nargs="?", help="(old usage) the true CAD, as --cad")
    ap.add_argument("--cad", help="the true CAD (.step): edges and faces measured against it")
    ap.add_argument("--step", help="a conversion (.step): do its edges lie on the contours?")
    ap.add_argument("--on", type=float, help="mm: an edge within this of a contour is on it (default 5 x refit's automatic tolerance)")
    ap.add_argument("--near", type=lambda s: _vec(s, 4), help="X,Y,Z,R: reports restricted to a ball")
    ap.add_argument("--png", help="picture of the contours (and of the --step's faces)")
    ap.add_argument("--at", type=lambda s: _vec(s, 3), help="X,Y,Z: centre of the picture (default: the --near ball's, or the part's)")
    ap.add_argument("--dir", type=lambda s: _vec(s, 3), default=[1.0, 1.0, 1.0], help="DX,DY,DZ: looking from there")
    ap.add_argument("--span", type=float, help="mm: half the picture's height (default: the --near radius, or half the part)")
    ap.add_argument("--save", help=".npy: per triangle centroid x, y, z and face label")
    ap.add_argument("--json", help="everything measured")
    ap.add_argument("--worst", type=int, default=15, help="how many edges off the contours to list")
    a = ap.parse_args(argv)
    if a.orig and not a.cad:
        if not a.orig.lower().endswith((".step", ".stp")):
            ap.error("the second argument is the true CAD, a .step")
        a.cad = a.orig
    shB = R.read_input(a.mesh)[0]
    topo = R.Topo(shB)
    c = Contours(topo)
    s = c.summary()
    rep = dict(mesh=os.path.abspath(a.mesh), summary=s)
    print(f"{os.path.basename(a.mesh)}: {s['triangles']} triangles, {s['edges']} edges; strong {s['strong']}, weak {s['weak']}, kept {s['kept']}; "
          f"faces {s['faces']}; contours {s['contour_length']:.1f} mm ({s['strong_length']:.1f} strong)")
    if a.near:
        pls = c.polylines(near=a.near)
        print(f"  near {a.near[:3]} r {a.near[3]}: {len(pls)} contour lines, {sum(p['length'] for p in pls):.2f} mm")
        rep["near_polylines"] = [dict(length=p["length"], strong=p["strong"], angle=p["angle"], closed=p["closed"], ends=[p["points"][0].tolist(), p["points"][-1].tolist()]) for p in pls]
    if a.cad:
        of = ground_truth(a.cad, shB, topo, c.m)
        r = evaluate(c.m, c.labels, of, c.areas)
        print(f"  CAD edges caught {100*r['caught']:.1f}%  false {100*r['false']:.2f}%  | smooth CAD edges ({r['smooth_cad']}) caught {100*r['smooth_caught']:.1f}%")
        print(f"  CAD faces {r['cad_faces']}: whole {r['whole']}, split {r['split']}, merged with another {r['merged']}  (predicted faces {r['pred_faces']})")
        rep["cad"] = r
    shape = None
    if a.step:
        shape = R.read_input(a.step)[0]
        f = c.follow(shape, on=a.on, near=a.near, n_worst=a.worst)
        Lt = max(f['length'], 1e-300)
        print(f"  {os.path.basename(a.step)}: {f['edges']} edges ({f['seams']} seams apart), {f['length']:.1f} mm; on the contours (within {f['on']:.4f} mm) "
              f"{100 * f['length_on'] / Lt:.1f}% of the length")
        print(f"    moved (a contour within a facet, {f['facet']:.3f} mm, and the edge not on it): {f['edges_moved']} edges, {f['length_moved']:.2f} mm ({100 * f['length_moved'] / Lt:.2f}%)")
        print(f"    unseen (no contour there - the detector's miss, or a boundary where the mesh shows none): {f['edges_unseen']} edges, {f['length_unseen']:.2f} mm ({100 * f['length_unseen'] / Lt:.2f}%)")
        if f["strong_length"] > 0:
            print(f"  strong contours no edge covers: {f['strong_uncovered']:.2f} of {f['strong_length']:.1f} mm ({100 * f['strong_uncovered'] / f['strong_length']:.2f}%)")
        if f["worst"]:
            print("  edges off the contours (moved first, the longest stretch off first):")
            for E in f["worst"]:
                print(f"    {E['class']:6s} {E['curve']:12s} L {E['length']:7.3f}  off {100 * (1 - E['on']):5.1f}%  max {E['max']:.4f}  at {np.round(E['mid'], 2).tolist()}  "
                      f"faces {E['faces']} {'/'.join(E['kinds'])}")
        if f["uncovered"]:
            print("  strong contours no edge covers (longest first):")
            for u in f["uncovered"][:a.worst]:
                print(f"    L {u['length']:.3f}  angle {u['angle']:5.1f}  at {np.round(u['mid'], 2).tolist()}")
        rep["follow"] = {k: v for k, v in f.items() if k != "all"}
    if a.png:
        at = a.at or (a.near[:3] if a.near else None)
        span = a.span or (a.near[3] if a.near else None)
        render(c, a.png, at=at, view=a.dir, span=span, shape=shape)
        print(f"  picture: {a.png}")
    if a.save:
        np.save(a.save, np.column_stack([c.cents, c.labels]))
    if a.json:
        with open(a.json, "w", encoding="utf-8") as fh:
            json.dump(rep, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    return 0


if __name__ == "__main__":
    sys.exit(main())
