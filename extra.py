"""Le facce di TROPPO: dove il CAD ha una faccia sola e noi ne mettiamo due."""
import sys
import numpy as np
sys.path.insert(0, "D:/Users/PC/Desktop/refit")
import refit as R
from align import tri_of
from icp import align, closest_on_tris
from collections import defaultdict

shA = R.read_input(sys.argv[1])[0]
shB = R.read_input(sys.argv[2])[0]
toA, med, mx, tA, fA, tipA, arA = align(shA, shB)
tB, fB, tipB, arB = tri_of(shB, 0.05)
print(f"allineamento ICP: medio {med:.5f} max {mx:.5f} mm")
# ogni triangolo NOSTRO -> faccia originale
cen = tB.mean(axis=1)
dd, _, jj = closest_on_tris(toA(cen), tA)
of = fA[jj]
# area di ogni nostra faccia su ogni faccia originale
peso = defaultdict(float)
ar_tri = 0.5 * np.linalg.norm(np.cross(tB[:, 1] - tB[:, 0], tB[:, 2] - tB[:, 0]), axis=1)
for k in range(len(tB)):
    peso[(int(of[k]), int(fB[k]))] += float(ar_tri[k])
sotto = defaultdict(list)
for (o, b), a in peso.items():
    sotto[o].append((a, b))
print(f"facce originali {len(tipA)}, nostre {len(tipB)}")
extra = 0
for o in sorted(sotto, key=lambda q: -len(sotto[q])):
    lst = sorted(sotto[o], reverse=True)
    if len(lst) < 2:
        continue
    extra += len(lst) - 1
    print(f"  {tipA[o].replace('GeomAbs_',''):16s} area {arA[o]:8.3f} -> {len(lst)} nostre:")
    for a, b in lst[:6]:
        print(f"        {tipB[b].replace('GeomAbs_',''):16s} area totale {arB[b]:8.3f}  "
              f"qui {a:8.3f}")
print(f"facce nostre in eccesso su facce originali condivise: {extra}")
