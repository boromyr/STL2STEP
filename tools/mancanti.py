"""Which BOUNDARIES between CAD faces escape the detector, and why."""
import sys, os, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import refit as R
from align import tri_of
from icp import align, closest_on_tris
from segnali import tri_metrics
from collections import defaultdict

shA = R.read_input(sys.argv[1])[0]
shB = R.read_input(sys.argv[2])[0]
topo = R.Topo(shB)
toA, med, mx, tA, fA, tipA, arA = align(shA, shB)
cent = np.array([topo.cents[i] for i in range(topo.nF)])
dd, _, jj = closest_on_tris(toA(cent), tA)
of = fA[jj]
nF = topo.nF
scala = np.zeros(nF); dirs = np.zeros((nF, 3))
for i in range(nF):
    scala[i], _, dirs[i], _ = tri_metrics(topo.verts[i])
N = np.array([topo.norms[i] for i in range(nF)])

coppie = defaultdict(lambda: [0, 0, [], [], []])   # (a,b) -> n, taken, angles, textures, cos
for i in range(nF):
    for j in topo.adj[i]:
        if j <= i or of[i] == of[j]:
            continue
        key = (int(min(of[i], of[j])), int(max(of[i], of[j])))
        ang = math.degrees(math.acos(float(np.clip(N[i] @ N[j], -1, 1))))
        ds = abs(math.log2(max(scala[i], 1e-12) / max(scala[j], 1e-12)))
        cd = abs(float(dirs[i] @ dirs[j]))
        rec = coppie[key]
        rec[0] += 1
        rec[1] += 1 if (ang > 30 or ds > 1.5) else 0
        rec[2].append(ang); rec[3].append(ds); rec[4].append(cd)
print(f"boundaries between original faces: {len(coppie)}")
righe = []
for (a, b), (n, presi, angs, dss, cds) in coppie.items():
    righe.append((presi / n, n, a, b, float(np.median(angs)), float(np.median(dss)),
                  float(np.median(cds))))
righe.sort()
persi3 = [r for r in righe if r[0] < 0.5]
persi3.sort(key=lambda r: -r[1])
print("  the MISSED boundaries, longest first:")
for f, n, a, b, ang, ds, cd in persi3[:14]:
    print(f"   {tipA[a].replace('GeomAbs_',''):14s} area {arA[a]:7.2f} | "
          f"{tipA[b].replace('GeomAbs_',''):14s} area {arA[b]:7.2f} : "
          f"{n:4d} edges, taken {100*f:5.1f}%  (dihedral {ang:6.2f}, texture {ds:5.2f}, cos {cd:5.2f})")
tot = sum(r[1] for r in righe)
ok = sum(r[0] * r[1] for r in righe)
persi = [r for r in righe if r[0] < 0.5]
print(f"  missed boundaries (under 50%): {len(persi)} out of {len(righe)}")
