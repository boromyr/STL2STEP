"""Anatomy of a B-Rep: faces by surface type, edges by curve type."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
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
    print(f"    faces {topo.nF}   edges {len(topo.edges)}   vertices {len(V)}")
    print(f"    surfaces: {dict(sf.most_common())}")
    print(f"    curves  : {dict(cu.most_common())}")
    n8 = sum(v for k, v in lati.items() if k > 8)
    print(f"    sides per face: max {max(lati)}  with >8 sides {n8}  "
          f"distribution {dict(sorted(lati.items())[:10])}")
    print(f"    volume {R.shape_volume(sh):.4f} mm3   bbox {(V.max(axis=0)-V.min(axis=0)).round(2)}")
    return sh, topo


if __name__ == "__main__":
    for p in sys.argv[1:]:
        info(p)
