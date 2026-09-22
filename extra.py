"""The EXTRA faces: where the CAD has a single face and we put two."""
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
print(f"ICP alignment: mean {med:.5f} max {mx:.5f} mm")
# every OUR triangle -> original face
cen = tB.mean(axis=1)
dd, _, jj = closest_on_tris(toA(cen), tA)
of = fA[jj]
# area of each of our faces on each original face
peso = defaultdict(float)
ar_tri = 0.5 * np.linalg.norm(np.cross(tB[:, 1] - tB[:, 0], tB[:, 2] - tB[:, 0]), axis=1)
for k in range(len(tB)):
    peso[(int(of[k]), int(fB[k]))] += float(ar_tri[k])
sotto = defaultdict(list)
for (o, b), a in peso.items():
    sotto[o].append((a, b))
print(f"original faces {len(tipA)}, ours {len(tipB)}")
extra = 0
for o in sorted(sotto, key=lambda q: -len(sotto[q])):
    lst = sorted(sotto[o], reverse=True)
    if len(lst) < 2:
        continue
    extra += len(lst) - 1
    print(f"  {tipA[o].replace('GeomAbs_',''):16s} area {arA[o]:8.3f} -> {len(lst)} of ours:")
    for a, b in lst[:6]:
        print(f"        {tipB[b].replace('GeomAbs_',''):16s} total area {arB[b]:8.3f}  "
              f"here {a:8.3f}")
print(f"our faces in excess over shared original faces: {extra}")
