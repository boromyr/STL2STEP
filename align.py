"""Aligns two B-Reps of the same part (inertia axes) and compares the faces.

usage:  python align.py original.step ours.step
"""
import sys, itertools
import numpy as np
sys.path.insert(0, "D:/Users/PC/Desktop/refit")
import refit as R
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.BRep import BRep_Tool
from OCP.TopLoc import TopLoc_Location
from OCP.TopAbs import TopAbs_FACE
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps


def tri_of(shape, defl):
    """(triangles (n,3,3), face id per triangle, surface types, areas)."""
    BRepMesh_IncrementalMesh(shape, defl, False, 0.3, True)
    tris, fid, tipi, aree = [], [], [], []
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    k = 0
    while ex.More():
        f = TopoDS.Face_s(ex.Current()); ex.Next()
        tipi.append(str(BRepAdaptor_Surface(f, True).GetType())[16:].replace("ype.", ""))
        aree.append(R.face_area(f))
        loc = TopLoc_Location()
        t = BRep_Tool.Triangulation_s(f, loc)
        if t is not None:
            tr = loc.Transformation()
            P = np.array([[t.Node(i).Transformed(tr).X(), t.Node(i).Transformed(tr).Y(),
                           t.Node(i).Transformed(tr).Z()] for i in range(1, t.NbNodes() + 1)])
            for i in range(1, t.NbTriangles() + 1):
                a, b, c = t.Triangle(i).Get()
                tris.append(P[[a - 1, b - 1, c - 1]])
                fid.append(k)
        k += 1
    return np.array(tris), np.array(fid), tipi, np.array(aree)


def inertia_frame(shape):
    g = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, g)
    cm = g.CentreOfMass()
    c = np.array([cm.X(), cm.Y(), cm.Z()])
    M = g.MatrixOfInertia()
    I = np.array([[M.Value(i, j) for j in (1, 2, 3)] for i in (1, 2, 3)])
    w, Vx = np.linalg.eigh(I)
    return c, Vx, w


def best_transform(shA, shB, PA, PB):
    """Rigid transform taking A onto B, chosen among the 24 axis combinations."""
    cA, VA, wA = inertia_frame(shA)
    cB, VB, wB = inertia_frame(shB)
    best = None
    for perm in itertools.permutations(range(3)):
        for sg in itertools.product((1, -1), repeat=3):
            Q = VA[:, list(perm)] * np.array(sg)
            if np.linalg.det(Q) < 0:
                continue
            Rm = VB @ Q.T                      # x_B = cB + Rm @ (x_A - cA)
            if abs(np.linalg.det(Rm) - 1) > 1e-6:
                continue
            Pt = cB + (PA - cA) @ Rm.T
            d = chamfer(Pt, PB)
            if best is None or d < best[0]:
                best = (d, cA, cB, Rm)
    return best


def chamfer(P, Q, cap=1500):
    a = P[np.linspace(0, len(P) - 1, min(cap, len(P))).astype(int)]
    b = Q[np.linspace(0, len(Q) - 1, min(cap, len(Q))).astype(int)]
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    return float(d.min(axis=1).mean() + d.min(axis=0).mean())


def point_tri_dist(P, T, chunk=400):
    """Distance of each point from the closest triangle: (dist, index)."""
    A, B, C = T[:, 0], T[:, 1], T[:, 2]
    out_d = np.full(len(P), np.inf)
    out_i = np.zeros(len(P), dtype=int)
    for s in range(0, len(T), chunk):
        a, b, c = A[s:s + chunk], B[s:s + chunk], C[s:s + chunk]
        ab, ac = b - a, c - a
        n = np.cross(ab, ac)
        nn = np.maximum(np.linalg.norm(n, axis=1), 1e-12)
        nu = n / nn[:, None]
        w = P[:, None, :] - a[None, :, :]
        h = np.einsum("pij,ij->pi", w, nu)
        proj = P[:, None, :] - h[:, :, None] * nu[None, :, :]
        v0 = proj - a[None, :, :]
        d00 = np.einsum("ij,ij->i", ab, ab)
        d01 = np.einsum("ij,ij->i", ab, ac)
        d11 = np.einsum("ij,ij->i", ac, ac)
        d20 = np.einsum("pij,ij->pi", v0, ab)
        d21 = np.einsum("pij,ij->pi", v0, ac)
        den = np.maximum(d00 * d11 - d01 * d01, 1e-20)
        u = (d11 * d20 - d01 * d21) / den
        v = (d00 * d21 - d01 * d20) / den
        inside = (u >= 0) & (v >= 0) & (u + v <= 1)
        d = np.where(inside, np.abs(h), np.inf)
        # outside the triangle: distance to the three sides
        for p0, p1 in ((a, b), (b, c), (c, a)):
            e = p1 - p0
            ee = np.maximum(np.einsum("ij,ij->i", e, e), 1e-20)
            tt = np.clip(np.einsum("pij,ij->pi", P[:, None, :] - p0[None, :, :], e) / ee, 0, 1)
            q = p0[None, :, :] + tt[:, :, None] * e[None, :, :]
            d = np.minimum(d, np.linalg.norm(P[:, None, :] - q, axis=2))
        j = np.argmin(d, axis=1)
        dm = d[np.arange(len(P)), j]
        m = dm < out_d
        out_d[m] = dm[m]
        out_i[m] = j[m] + s
    return out_d, out_i


def main():
    shA = R.read_input(sys.argv[1])[0]          # original
    shB = R.read_input(sys.argv[2])[0]          # ours
    diag = 1.0
    tA, fA, tipA, arA = tri_of(shA, 0.05)
    tB, fB, tipB, arB = tri_of(shB, 0.05)
    PA = tA.reshape(-1, 3)
    PB = tB.reshape(-1, 3)
    d, cA, cB, Rm = best_transform(shA, shB, PA, PB)
    print(f"alignment: mean residual {d:.4f} mm")
    # bring OURS into the original's frame
    Rinv = Rm.T
    def toA(P):
        return cA + (P - cB) @ Rinv.T
    cen = tB.mean(axis=1)
    cenA = toA(cen)
    dd, jj = point_tri_dist(cenA, tA)
    print(f"distance of our triangles' centroids from the original: "
          f"mean {dd.mean():.4f}  p99 {np.percentile(dd,99):.4f}  max {dd.max():.4f} mm")
    # for each original face: how many OUR faces sit on it
    origf = fA[jj]                    # original face under each of our triangles
    from collections import defaultdict
    sotto = defaultdict(set)
    for tri_i, of in enumerate(origf):
        sotto[of].add(fB[tri_i])
    print(f"\noriginal faces {len(tipA)}, ours {len(tipB)}")
    print("original face                          -> our faces")
    righe = []
    for of in range(len(tipA)):
        ns = sotto.get(of, set())
        tt = sorted({tipB[n] for n in ns})
        righe.append((len(ns), arA[of], tipA[of], tt))
    righe.sort(key=lambda r: -r[0])
    for n, a, t, tt in righe[:20]:
        print(f"  {t:16s} area {a:8.3f}  -> {n:4d} faces  {','.join(tt)}")
    tot = sum(r[0] for r in righe)
    print(f"  ... total assigned faces of ours: {tot}")


if __name__ == "__main__":
    main()
