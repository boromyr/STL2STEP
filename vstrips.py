# -*- coding: utf-8 -*-
import sys, numpy as np, collections
sys.path.insert(0, "D:/Users/PC/Desktop/refit")
import refit as R
from OCP.BRepAdaptor import BRepAdaptor_Surface

sh, _ = R.read_input(sys.argv[1])
topo = R.Topo(sh)
out = []
for i in range(topo.nF):
    V = topo.verts[i]
    if V.size == 0: continue
    d = V.max(axis=0) - V.min(axis=0)
    if str(BRepAdaptor_Surface(topo.faces[i]).GetType()).split("_")[-1] != "Plane": continue
    if d[2] > 4 and d[0] < 1.2 and d[1] < 1.2:
        out.append((i, V.mean(axis=0), topo.norms[i], topo.areas[i]))
print(len(out), "strisce verticali")
# raggruppa per xy
g = collections.defaultdict(list)
for i, c, n, a in out:
    g[(round(c[0]), round(c[1]))].append((i, np.round(c,2), np.round(n,3), round(a,2)))
for k in sorted(g, key=lambda k: -len(g[k]))[:10]:
    print("gruppo", k, len(g[k]))
    for r in g[k]: print("   ", r)
