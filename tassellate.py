# -*- coding: utf-8 -*-
"""Quanta CURVA e' rimasta tassellata nell'uscita: facce piane che confinano
con altre facce piane con un diedro fra 1 e 30 gradi = una banda curva spezzata."""
import sys, math
sys.path.insert(0, "D:/Users/PC/Desktop/refit")
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
    print(f"    facce {t.nF}  piane {int(pian.sum())}  analitiche {int((~pian).sum())}")
    print(f"    piane che fanno parte di una banda curva tassellata: {int(curva.sum())} "
          f"({100*A[curva].sum()/tot:.1f}% dell'area)")
    print(f"    area coperta da facce analitiche: {100*A[~pian].sum()/tot:.1f}%")
