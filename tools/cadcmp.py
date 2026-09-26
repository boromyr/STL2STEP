"""Face-by-face comparison of an output B-Rep with the original CAD.
  python cadcmp.py original.step ours.step [-v]
Same frame assumed when the bounding boxes agree, otherwise ICP (icp.align).

Every face of ours is sampled and each sample is assigned to the nearest face
of the original. Reported:
  - clean:   original faces covered by ONE face of ours, of the same type
  - split:   original faces covered by several of ours (with the pieces)
  - wrong type, and our faces spanning several originals (merged faces)
  - planar pieces of ours on curved originals: the mesh left behind
  - analytic parameters (radius, cone angle) of ours against the CAD's
  - our edges lying on no original edge, and our straight segments along
    curved original edges (polylines where the CAD has a curve)
-v lists the details."""
import sys, os, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import refit as R
from collections import defaultdict, Counter
from scipy.spatial import cKDTree
from OCP.BRepAdaptor import BRepAdaptor_Surface, BRepAdaptor_Curve
R.Log.level = 40


def tri_of(shape, defl):
    from OCP.BRepMesh import BRepMesh_IncrementalMesh
    BRepMesh_IncrementalMesh(shape, defl, False, 0.2, True)
    tris, fid, faces = [], [], []
    for k, f in enumerate(R.explore(shape, R.TopAbs_FACE)):
        f = R.td_Face(f)
        faces.append(f)
        loc = R.TopLoc_Location()
        t = R.bt_Triangulation(f, loc)
        for d2 in (defl * 4, defl * 20, defl * 100):
            if t is not None:
                break
            BRepMesh_IncrementalMesh(f, d2, False, 0.5, True)
            t = R.bt_Triangulation(f, loc)
        if t is None:
            print("no triangulation for face", k)
            continue
        tr = loc.Transformation()
        P = np.array([[q.X(), q.Y(), q.Z()] for q in (t.Node(i).Transformed(tr) for i in range(1, t.NbNodes() + 1))])
        I = np.array([t.Triangle(i).Get() for i in range(1, t.NbTriangles() + 1)]) - 1
        tris.append(P[I]); fid.append(np.full(len(I), k))
    return np.concatenate(tris), np.concatenate(fid), faces


def stype(f):
    return str(BRepAdaptor_Surface(f, True).GetType()).split("_")[-1].replace("Surface", "")


def sparams(f):
    s = BRepAdaptor_Surface(f, True)
    t = stype(f)
    if t == "Cylinder":
        c = s.Cylinder(); return dict(r=c.Radius())
    if t == "Cone":
        c = s.Cone(); return dict(a=math.degrees(c.SemiAngle()))
    if t == "Sphere":
        return dict(r=s.Sphere().Radius())
    if t == "Torus":
        c = s.Torus(); return dict(R=c.MajorRadius(), r=c.MinorRadius())
    return {}


def main():
    A = R.read_input(sys.argv[1])[0]
    B = R.read_input(sys.argv[2])[0]
    verb = "-v" in sys.argv
    tA, fA, facesA = tri_of(A, 0.005)
    tB, fB, facesB = tri_of(B, 0.005)
    PA, PB = tA.reshape(-1, 3), tB.reshape(-1, 3)
    toA = None
    if np.abs(PA.min(0) - PB.min(0)).max() > 0.05 or np.abs(PA.max(0) - PB.max(0)).max() > 0.05:
        from icp import align
        toA = align(A, B)[0]
        tB = toA(tB.reshape(-1, 3)).reshape(-1, 3, 3)
        print("frames differ: ICP alignment")
    # dense samples of A, tagged with the face
    from qa import _split
    diag = float(np.linalg.norm(PA.max(0) - PA.min(0)))
    h = diag / 600
    S, SF = [], []
    for k in np.unique(fA):
        sub = _split(tA[fA == k], h)
        S.append(np.vstack([sub.mean(1), sub.reshape(-1, 3)])); SF.append(np.full(len(sub) * 4, k))
    S, SF = np.vstack(S), np.concatenate(SF)
    tree = cKDTree(S)
    # our triangles, split too, so that a big triangle spanning two CAD faces is shared out
    arB = np.zeros(len(facesB))
    w = defaultdict(float)
    dmax = np.zeros(len(facesB))
    for k in np.unique(fB):
        sub = _split(tB[fB == k], h)
        ar = 0.5 * np.linalg.norm(np.cross(sub[:, 1] - sub[:, 0], sub[:, 2] - sub[:, 0]), axis=1)
        d, j = tree.query(sub.mean(1))
        o = SF[j]
        arB[k] = ar.sum()
        dmax[k] = d.max()
        arr = np.zeros(len(facesA))
        np.add.at(arr, o, ar)
        for oo in np.nonzero(arr)[0]:
            w[(int(oo), int(k))] = float(arr[oo])
    arA = np.array([R.face_area(f) for f in facesA])
    tipA = [stype(f) for f in facesA]
    tipB = [stype(f) for f in facesB]
    on = defaultdict(list)     # original -> [(area, ours)]
    of = defaultdict(list)     # ours -> [(area, original)]
    for (o, b), a in w.items():
        on[o].append((a, b)); of[b].append((a, o))
    # an original face is "clean" when one of our faces covers >= 97% of it with the same type
    clean, split, wrongtype, miss = 0, [], [], []
    for o in range(len(facesA)):
        lst = sorted(on.get(o, []), reverse=True)
        sig = [(a, b) for a, b in lst if a > max(1e-3, 0.01 * arA[o]) and (a > 0.02 * arB[b] or b == lst[0][1])]
        if not sig:
            miss.append(o); continue
        a0, b0 = sig[0]
        if tipB[b0] != tipA[o]:
            wrongtype.append((o, b0, a0))
        if len(sig) > 1:
            split.append((o, sig))
        elif tipB[b0] == tipA[o]:
            clean += 1
    # our faces spanning several originals (under-segmentation)
    merged = []
    for b, lst in of.items():
        sig = sorted([(a, o) for a, o in lst if a > max(1e-3, 0.03 * arB[b])], reverse=True)
        if len(sig) > 1:
            merged.append((b, sig))
    print(f"original {len(facesA)} faces {dict(Counter(tipA))}")
    print(f"ours     {len(facesB)} faces {dict(Counter(tipB))}")
    print(f"clean (1:1 same type): {clean}/{len(facesA)}   split: {len(split)}   wrong type: {len(wrongtype)}   "
          f"our faces spanning >1 original: {len(merged)}")
    # residual tessellation: our planar faces on non-planar originals
    tess = defaultdict(lambda: [0, 0.0])
    for (o, b), a in w.items():
        if tipB[b] == "Plane" and tipA[o] != "Plane":
            tess[o][0] += 1; tess[o][1] += a
    ta = sum(v[1] for v in tess.values())
    print(f"planar pieces of ours on curved originals: area {ta:.3f} mm2 over {len(tess)} originals")
    # analytic parameter check
    perr = []
    for o in range(len(facesA)):
        lst = sorted(on.get(o, []), reverse=True)
        if not lst:
            continue
        b = lst[0][1]
        if tipA[o] == tipB[b] and tipA[o] in ("Cylinder", "Cone", "Sphere", "Torus"):
            pa, pb = sparams(facesA[o]), sparams(facesB[b])
            e = max(abs(pa[k] - pb[k]) for k in pa)
            perr.append((e, o, b, pa, pb))
    perr.sort(key=lambda r: -r[0])
    if perr:
        print(f"analytic params, worst: " + ", ".join(f"{tipA[o]} {pa}->{ {k: round(v, 4) for k, v in pb.items()} }" for e, o, b, pa, pb in perr[:4]))
    def lab(o):
        return f"{tipA[o]}#{o}({arA[o]:.2f})"
    if split:
        print("SPLIT originals (largest first):")
        for o, sig in sorted(split, key=lambda r: -arA[r[0]])[:25 if verb else 12]:
            print(f"   {lab(o)} -> " + ", ".join(f"{tipB[b]}#{b} {a:.3f}" for a, b in sig[:6]) + (" ..." if len(sig) > 6 else ""))
    if wrongtype:
        print("WRONG TYPE (dominant face):")
        for o, b, a in sorted(wrongtype, key=lambda r: -arA[r[0]])[:25 if verb else 12]:
            print(f"   {lab(o)} -> {tipB[b]}#{b} {a:.3f}")
    if merged:
        print("OUR faces over several originals:")
        for b, sig in sorted(merged, key=lambda r: -arB[r[0]])[:15]:
            print(f"   {tipB[b]}#{b}({arB[b]:.2f}) <- " + ", ".join(f"{lab(o)} {a:.3f}" for a, o in sig[:5]))
    if miss:
        print(f"originals with nothing on them: {[lab(o) for o in miss[:10]]}")
    # edges
    def ecount(sh):
        c = Counter()
        for e in R.explore(sh, R.TopAbs_EDGE):
            try:
                c[str(BRepAdaptor_Curve(R.td_Edge(e)).GetType()).split("_")[-1]] += 1
            except Exception:
                c["?"] += 1
        return dict(c)
    print(f"edges original {ecount(A)}\n      ours     {ecount(B)}")
    if verb and tess:
        print("planar residue on curved originals:")
        for o, (n, a) in sorted(tess.items(), key=lambda r: -r[1][1])[:20]:
            print(f"   {lab(o)}: {n} planar faces of ours, {a:.4f} mm2")
    # edge mapping: every edge of ours -> original edge under it
    from OCP.GCPnts import GCPnts_UniformAbscissa
    def edge_samples(sh, n=24):
        out, eid, etyp, lens = [], [], [], []
        for k, e in enumerate(R.explore(sh, R.TopAbs_EDGE)):
            e = R.td_Edge(e)
            try:
                c = BRepAdaptor_Curve(e)
            except Exception:
                continue
            t = str(c.GetType()).split("_")[-1]
            u0, u1 = c.FirstParameter(), c.LastParameter()
            m = int(min(20000, max(50, R.edge_length(e) / 0.004))) if n > 50 else n
            P = np.array([[p.X(), p.Y(), p.Z()] for p in (c.Value(u0 + (u1 - u0) * (i + 0.5) / m) for i in range(m))])
            out.append(P); eid.append(np.full(m, k)); etyp.append(t)
            lens.append(float(R.edge_length(e)))
        return np.vstack(out), np.concatenate(eid), etyp, np.array(lens)
    EA, eA, tyA, lA = edge_samples(A, 200)
    EB, eB, tyB, lB = edge_samples(B, 8)
    if toA is not None:
        EB = toA(EB)
    tr = cKDTree(EA)
    d, j = tr.query(EB)
    # per our edge: majority original edge, and max distance
    per = defaultdict(list)
    for k, dd, jj in zip(eB, d, j):
        per[int(k)].append((float(dd), int(eA[jj])))
    bad = defaultdict(lambda: [0, 0.0])   # original edge -> (our line pieces, length) where the original is curved
    for k, lst in per.items():
        dd = max(x[0] for x in lst)
        o = Counter(x[1] for x in lst).most_common(1)[0][0]
        if dd > 0.02:
            continue      # not on an original edge (tangent seams etc.)
        if tyB[k] == "Line" and tyA[o] != "Line":
            bad[o][0] += 1; bad[o][1] += lB[k]
    ph = Counter(); phl = Counter(); phpts = []
    for k, lst in per.items():
        if max(x[0] for x in lst) > 0.02:
            ph[tyB[k]] += 1; phl[tyB[k]] += lB[k]
            phpts.append((lB[k], tyB[k], EB[eB == k].mean(0)))
    print(f"our edges NOT on any original edge: {dict(ph)}  length { {k: round(v, 2) for k, v in phl.items()} }")
    if verb:
        for L, t, c in sorted(phpts, key=lambda r: -r[0])[:12]:
            print(f"   {t} len {L:.3f} at {c.round(2)}")
    nb = sum(v[0] for v in bad.values())
    print(f"our straight segments along curved original edges: {nb} over {len(bad)} original edges")
    if verb:
        for o, (n, L) in sorted(bad.items(), key=lambda r: -r[1][0])[:25]:
            P = EA[eA == o]
            print(f"   {tyA[o]}#{o} len {lA[o]:.2f} -> {n} segments ({L:.2f})  mid {P[len(P)//2].round(2)}")


if __name__ == "__main__":
    main()
