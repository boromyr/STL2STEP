"""Anatomia di un B-Rep: facce per tipo di superficie, spigoli per tipo di curva."""
import sys
import numpy as np
sys.path.insert(0, "D:/Users/PC/Desktop/refit")
import refit as R
from collections import Counter
from OCP.BRepAdaptor import BRepAdaptor_Surface, BRepAdaptor_Curve


def info(path):
    sh = R.read_input(path)[0]
    topo = R.Topo(sh)
    sf = Counter()
    lati = Counter()
    for i in range(topo.nF):
        sf[str(BRepAdaptor_Surface(topo.faces[i], True).GetType())[16:].replace("GeomAbs_", "")] += 1
        lati[len(topo.f_edges[i])] += 1
    cu = Counter()
    for e in topo.edges:
        try:
            cu[str(BRepAdaptor_Curve(e).GetType())[14:].replace("GeomAbs_", "")] += 1
        except Exception:
            cu["?"] += 1
    V = topo.vpos
    print(f"--- {path}")
    print(f"    facce {topo.nF}   spigoli {len(topo.edges)}   vertici {len(V)}")
    print(f"    superfici: {dict(sf.most_common())}")
    print(f"    curve    : {dict(cu.most_common())}")
    n8 = sum(v for k, v in lati.items() if k > 8)
    print(f"    lati per faccia: max {max(lati)}  con >8 lati {n8}  "
          f"distribuzione {dict(sorted(lati.items())[:10])}")
    print(f"    volume {R.shape_volume(sh):.4f} mm3   bbox {(V.max(axis=0)-V.min(axis=0)).round(2)}")
    return sh, topo


if __name__ == "__main__":
    for p in sys.argv[1:]:
        info(p)
