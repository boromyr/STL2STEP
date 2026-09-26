# -*- coding: utf-8 -*-
"""How much CURVATURE is left tessellated in the output: planar faces bordering
other planar faces with a dihedral angle between 1 and 30 degrees = a curved band cut into pieces."""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import refit as R
from OCP.BRepAdaptor import BRepAdaptor_Surface

for path in sys.argv[1:]:
    sh = R.read_input(path)[0]
    t = R.Topo(sh)
    pian = np.array([str(BRepAdaptor_Surface(t.faces[i], True).GetType()).endswith("Plane")
                     for i in range(t.nF)])
    N = np.array([t.norms[i] for i in range(t.nF)])
    A = np.array([t.areas[i] for i in range(t.nF)])
    curva = np.zeros(t.nF, bool)
    for i in range(t.nF):
        if not pian[i]:
            continue
        for j in t.adj[i]:
            if not pian[j]:
                continue
            c = float(np.clip(N[i] @ N[j], -1, 1))
            a = math.degrees(math.acos(c))
            if 0.8 <= a <= 30.0:
                curva[i] = True
                break
    tot = A.sum()
    print(f"--- {path}")
    print(f"    faces {t.nF}  planar {int(pian.sum())}  analytic {int((~pian).sum())}")
    print(f"    planar faces that are part of a tessellated curved band: {int(curva.sum())} "
          f"({100*A[curva].sum()/tot:.1f}% of the area)")
    print(f"    area covered by analytic faces: {100*A[~pian].sum()/tot:.1f}%")
