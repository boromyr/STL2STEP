"""For each face of the original CAD: do the mesh facets that belong to it
really sit on ONE primitive, within our tolerance?"""
import sys
import numpy as np
sys.path.insert(0, "D:/Users/PC/Desktop/refit")
import refit as R
from align import tri_of
from icp import align, closest_on_tris
from trama3 import load_brep

shA = R.read_input(sys.argv[1])[0]
shB = load_brep(sys.argv[2])
topo = R.Topo(shB)
toA, med, mx, tA, fA, tipA, arA = align(shA, shB)
print(f"ICP alignment: mean {med:.5f} max {mx:.5f} mm")
cent = np.array([topo.cents[i] for i in range(topo.nF)])
dd, _, jj = closest_on_tris(toA(cent), tA)
of = fA[jj]
diag = float(np.linalg.norm(topo.vpos.max(axis=0) - topo.vpos.min(axis=0)))
tf = max(2e-4, 1e-5 * diag)
sg = R.Segmenter(topo, tf, min(10 * tf, 1e-3 * diag), diag)
print(f"tol_fit {tf:.2e}")
for o in range(len(tipA)):
    m = [i for i in range(topo.nF) if of[i] == o and dd[i] < 0.05]
    if len(m) < 8 or "Plane" in tipA[o]:
        continue
    P = np.vstack([topo.verts[i] for i in m])
    p = R.refit_best(m, topo.verts, topo.norms, topo.areas)
    if p is None:
        print(f"  {tipA[o].replace('GeomAbs_',''):16s} area {arA[o]:8.3f} {len(m):5d} faces -> fit None")
        continue
    p2 = sg._snap_cylinder(p, m)
    d = np.abs(p2.dist(P))
    fuori = [i for i in m if not sg.within(p2, i)]
    if not fuori or "BSpline" in tipA[o]:
        print(f"  {tipA[o].replace('GeomAbs_',''):16s} area {arA[o]:8.3f} {len(m):5d} faces -> "
              f"{p2.label():7s} max residual {d.max():.2e}  outside {len(fuori)}")
        continue
    print(f"  {tipA[o].replace('GeomAbs_',''):16s} area {arA[o]:8.3f} {len(m):5d} faces -> "
          f"{p2.label():7s} max residual {d.max():.2e}  outside {len(fuori)}")
    for i in fuori:
        V = topo.verts[i]
        res = float(np.abs(p2.dist(V)).max())
        print(f"       facet area {topo.areas[i]:.5f}  residual {res:.2e}  "
              f"threshold {sg.face_tol(p2, i):.2e}  dist from CAD {dd[i]:.2e}")
