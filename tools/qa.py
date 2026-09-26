"""
QA of a refit output against its input mesh.
  python qa.py input.stl output.step [--json out.json]
Two-sided deviation (output->mesh and mesh->output), B-Rep statistics:
face types, edge curve types, tolerances, validity, free edges, volume.
"""
import sys, json, math, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.spatial import cKDTree
import refit as rf

rf.Log.level = 40


def mesh_tris_from_shape(shape, defl):
    from OCP.BRepMesh import BRepMesh_IncrementalMesh
    BRepMesh_IncrementalMesh(shape, defl, False, 0.1, True)
    T = []
    for f in rf.explore(shape, rf.TopAbs_FACE):
        f = rf.td_Face(f)
        loc = rf.TopLoc_Location()
        tri = rf._st(rf.BRep_Tool, "Triangulation")(f, loc)
        if tri is None:
            continue
        tr = loc.Transformation()
        nodes = np.array([[q.X(), q.Y(), q.Z()] for q in (tri.Node(k).Transformed(tr) for k in range(1, tri.NbNodes() + 1))])
        for k in range(1, tri.NbTriangles() + 1):
            a, b, c = tri.Triangle(k).Get()
            T.append((nodes[a - 1], nodes[b - 1], nodes[c - 1]))
    return np.array(T)


def stl_tris(path):
    RWStl = rf._m("RWStl").RWStl
    tri = rf._st(RWStl, "ReadFile")(path)
    nodes = np.array([[tri.Node(k).X(), tri.Node(k).Y(), tri.Node(k).Z()] for k in range(1, tri.NbNodes() + 1)])
    idx = np.array([tri.Triangle(k).Get() for k in range(1, tri.NbTriangles() + 1)]) - 1
    return nodes[idx]


def sample(T, n, rng):
    A, B, C = T[:, 0], T[:, 1], T[:, 2]
    ar = 0.5 * np.linalg.norm(np.cross(B - A, C - A), axis=1)
    p = ar / ar.sum()
    k = rng.choice(len(T), n, p=p)
    r1, r2 = rng.random(n), rng.random(n)
    s = np.sqrt(r1)
    return (1 - s)[:, None] * A[k] + (s * (1 - r2))[:, None] * B[k] + (s * r2)[:, None] * C[k]


def pt_tri(P, A, B, C):
    """exact distance point-triangle, row-wise (P,A,B,C: m x 3)."""
    ab, ac, ap = B - A, C - A, P - A
    d1, d2 = np.einsum("ij,ij->i", ab, ap), np.einsum("ij,ij->i", ac, ap)
    bp = P - B
    d3, d4 = np.einsum("ij,ij->i", ab, bp), np.einsum("ij,ij->i", ac, bp)
    cp = P - C
    d5, d6 = np.einsum("ij,ij->i", ab, cp), np.einsum("ij,ij->i", ac, cp)
    va = d3 * d6 - d5 * d4
    vb = d5 * d2 - d1 * d6
    vc = d1 * d4 - d3 * d2
    den = va + vb + vc
    den = np.where(np.abs(den) < 1e-300, 1e-300, den)
    v = vb / den
    w = vc / den
    Q = A + v[:, None] * ab + w[:, None] * ac
    # regions
    m = (d1 <= 0) & (d2 <= 0); Q[m] = A[m]
    m2 = (d3 >= 0) & (d4 <= d3); Q[m2] = B[m2]
    m3 = (vc <= 0) & (d1 >= 0) & (d3 <= 0) & ~m & ~m2
    t = d1 / np.where(np.abs(d1 - d3) < 1e-300, 1e-300, d1 - d3)
    Q[m3] = (A + t[:, None] * ab)[m3]
    m4 = (d6 >= 0) & (d5 <= d6) & ~m & ~m2
    Q[m4] = C[m4]
    m5 = (vb <= 0) & (d2 >= 0) & (d6 <= 0) & ~m & ~m2 & ~m4
    t = d2 / np.where(np.abs(d2 - d6) < 1e-300, 1e-300, d2 - d6)
    Q[m5] = (A + t[:, None] * ac)[m5]
    m6 = (va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0) & ~m & ~m2 & ~m4
    t = (d4 - d3) / np.where(np.abs((d4 - d3) + (d5 - d6)) < 1e-300, 1e-300, (d4 - d3) + (d5 - d6))
    Q[m6] = (B + t[:, None] * (C - B))[m6]
    return np.linalg.norm(P - Q, axis=1)


def _split(T, h):
    """longest-edge bisection until every edge <= h (exact cover of T)."""
    out = []
    while len(T):
        e = np.stack([np.linalg.norm(T[:, 1] - T[:, 0], axis=1), np.linalg.norm(T[:, 2] - T[:, 1], axis=1), np.linalg.norm(T[:, 0] - T[:, 2], axis=1)], 1)
        ok = e.max(1) <= h
        out.append(T[ok])
        T = T[~ok]
        if not len(T):
            break
        k = e[~ok].argmax(1)
        # rotate so the longest edge is (0,1)
        idx = np.stack([k, (k + 1) % 3, (k + 2) % 3], 1)
        T = np.take_along_axis(T, idx[:, :, None], 1)
        M = 0.5 * (T[:, 0] + T[:, 1])
        T = np.concatenate([np.stack([T[:, 0], M, T[:, 2]], 1), np.stack([M, T[:, 1], T[:, 2]], 1)])
    return np.concatenate(out)


def dist_to_tris(P, T, k=8):
    """EXACT point-to-mesh distance: small sub-triangles, k nearest centroids
    as a first bound, then every sub-triangle whose centroid lies within
    bound + circumradius is checked."""
    diag = float(np.linalg.norm(T.reshape(-1, 3).max(0) - T.reshape(-1, 3).min(0)))
    h = diag / 300.0
    S = _split(T, h)
    A, B, C = S[:, 0], S[:, 1], S[:, 2]
    cen = S.mean(1)
    rad = np.max(np.linalg.norm(S - cen[:, None, :], axis=2), axis=1)
    rmax = float(rad.max())
    tree = cKDTree(cen)
    _, nn = tree.query(P, k=k)
    best = np.full(len(P), np.inf)
    for j in range(k):
        t = nn[:, j]
        best = np.minimum(best, pt_tri(P, A[t], B[t], C[t]))
    # exactness: any triangle closer than best has its centroid within best + rmax
    cand = tree.query_ball_point(P, best + rmax)
    qi = np.repeat(np.arange(len(P)), [len(c) for c in cand])
    ti = np.concatenate([np.asarray(c, dtype=int) for c in cand]) if len(qi) else np.zeros(0, int)
    for s in range(0, len(qi), 2_000_000):
        q, t = qi[s:s + 2_000_000], ti[s:s + 2_000_000]
        d = pt_tri(P[q], A[t], B[t], C[t])
        np.minimum.at(best, q, d)
    return best


def brep_stats(shape):
    from collections import Counter
    fc = Counter()
    for f in rf.explore(shape, rf.TopAbs_FACE):
        fc[str(rf.face_surface_type(f)).split("_")[-1]] += 1
    ec = Counter()
    many = 0
    for e in rf.explore(shape, rf.TopAbs_EDGE):
        try:
            ec[str(rf.BRepAdaptor_Curve(rf.td_Edge(e)).GetType()).split("_")[-1]] += 1
        except Exception:
            ec["?"] += 1
    for f in rf.explore(shape, rf.TopAbs_FACE):
        if rf.count_sub(f, rf.TopAbs_EDGE) > 8:
            many += 1
    etol = [float(rf.bt_Tolerance(rf.td_Edge(e))) for e in rf.explore(shape, rf.TopAbs_EDGE)]
    vtol = [float(rf.bt_Tolerance(rf.td_Vertex(v))) for v in rf.explore(shape, rf.TopAbs_VERTEX)]
    return dict(faces=sum(fc.values()), face_types=dict(fc), edges=sum(ec.values()), edge_types=dict(ec), faces_gt8=many,
                etol_max=max(etol), etol_p99=float(np.percentile(etol, 99)), etol_mean=float(np.mean(etol)),
                vtol_max=max(vtol), vtol_p99=float(np.percentile(vtol, 99)))


def main():
    inp, out = sys.argv[1], sys.argv[2]
    jout = sys.argv[sys.argv.index("--json") + 1] if "--json" in sys.argv else None
    t0 = time.time()
    Tm = stl_tris(inp)
    sh, _ = rf.read_step(out)
    st = brep_stats(sh)
    diag = float(np.linalg.norm(Tm.reshape(-1, 3).max(0) - Tm.reshape(-1, 3).min(0)))
    To = mesh_tris_from_shape(sh, 2e-5 * diag)
    rng = np.random.default_rng(0)
    n = 100000
    Po = sample(To, n, rng)
    d_om = dist_to_tris(Po, Tm)
    Pm = sample(Tm, n, rng)
    d_mo = dist_to_tris(Pm, To)
    res = dict(file=os.path.basename(inp), valid=rf.is_valid(sh), free=rf.count_free_edges(sh), volume=rf.shape_volume(sh),
               out2mesh_mean=float(d_om.mean()), out2mesh_p99=float(np.percentile(d_om, 99)), out2mesh_max=float(d_om.max()),
               mesh2out_mean=float(d_mo.mean()), mesh2out_p99=float(np.percentile(d_mo, 99)), mesh2out_max=float(d_mo.max()),
               out2mesh_at=[round(float(x), 3) for x in Po[int(d_om.argmax())]],
               mesh2out_at=[round(float(x), 3) for x in Pm[int(d_mo.argmax())]], **st)
    print(json.dumps(res, indent=1))
    if jout:
        with open(jout, "w") as fh:
            json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
