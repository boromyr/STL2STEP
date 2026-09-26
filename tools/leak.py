"""How closely do our boundaries follow the CAD edges?

For each REGION (and each SECTION) we look at which original face its
facets belong to: the ones that don't sit in the majority face are
"bleed", i.e. the boundary zigzags across a CAD edge instead of
following it.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import refit as R
from align import tri_of, best_transform, point_tri_dist
from icp import align, closest_on_tris
from trama3 import load_brep
from collections import defaultdict


def mappa(orig, topo):
    shA = R.read_input(orig)[0]
    tA, fA, tipA, arA = tri_of(shA, 0.05)
    PA = tA.reshape(-1, 3)
    d, cA, cB, Rm = best_transform(shA, R.Topo, PA, topo.vpos) if False else \
        best_transform_shapes(shA, topo, PA)
    return d, cA, cB, Rm, tA, fA, tipA, arA


def best_transform_shapes(shA, topo, PA):
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    import itertools

    def frame_pts(P):
        c = P.mean(axis=0)
        H = (P - c).T @ (P - c)
        w, V = np.linalg.eigh(H)
        return c, V

    cA, VA = frame_pts(PA)
    cB, VB = frame_pts(topo.vpos)
    best = None
    for perm in itertools.permutations(range(3)):
        for sg in itertools.product((1, -1), repeat=3):
            Q = VA[:, list(perm)] * np.array(sg)
            if np.linalg.det(Q) < 0:
                continue
            Rm = VB @ Q.T
            Pt = cB + (PA - cA) @ Rm.T
            a = Pt[np.linspace(0, len(Pt) - 1, 900).astype(int)]
            b = topo.vpos[np.linspace(0, len(topo.vpos) - 1, 900).astype(int)]
            dd = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
            sc = float(dd.min(axis=1).mean() + dd.min(axis=0).mean())
            if best is None or sc < best[0]:
                best = (sc, cA, cB, Rm)
    return best


def main():
    orig, brep = sys.argv[1], sys.argv[2]
    shB = load_brep(brep)
    topo = R.Topo(shB)
    shA = R.read_input(orig)[0]
    toA, med, mx, tA, fA, tipA, arA = align(shA, shB)
    cent = np.array([topo.cents[i] for i in range(topo.nF)])
    dd, _, jj = closest_on_tris(toA(cent), tA)
    of = fA[jj]
    print(f"ICP alignment: mean {med:.5f} max {mx:.5f} mm · "
          f"centroids: mean {dd.mean():.5f} max {dd.max():.5f} mm")

    diag = float(np.linalg.norm(topo.vpos.max(axis=0) - topo.vpos.min(axis=0)))
    tf = max(2e-4, 1e-5 * diag)
    regs = R.segment_curved(topo, tf, min(10 * tf, 1e-3 * diag), diag, min_faces=4, threads=1)
    areas = np.array([topo.areas[i] for i in range(topo.nF)])
    sg = R.Segmenter(topo, tf, min(10 * tf, 1e-3 * diag), diag)
    lab, nc = sg.sections()

    def sbava(gruppi, nome):
        tot_bad, tot_f, righe = 0, 0, []
        for g in gruppi:
            if len(g) < 4:
                continue
            c = defaultdict(int)
            for i in g:
                c[of[i]] += 1
            top = max(c.values())
            bad = len(g) - top
            tot_bad += bad
            tot_f += len(g)
            if bad:
                righe.append((bad, len(g), sorted(c.values(), reverse=True)[:4]))
        righe.sort(reverse=True)
        print(f"{nome}: {len(gruppi)} groups, facets {tot_f}, bleed {tot_bad} "
              f"({100.0*tot_bad/max(tot_f,1):.1f}%)")
        for bad, n, det in righe[:8]:
            print(f"    group of {n:5d} facets: {bad:4d} outside  (breakdown {det})")

    sbava([list(Rg.faces) for Rg in regs], "REGIONS")
    sbava([list(np.where(lab == c)[0]) for c in range(nc)], "SECTIONS")


if __name__ == "__main__":
    main()
