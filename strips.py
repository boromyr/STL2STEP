# -*- coding: utf-8 -*-
import sys, numpy as np, collections
sys.path.insert(0, "D:/Users/PC/Desktop/refit")
import refit as R

def main(path):
    sh, _ = R.read_input(path)
    topo = R.Topo(sh)
    kinds = collections.Counter()
    strips = []
    for i in range(topo.nF):
        s = R.surf_of(topo.faces[i]) if hasattr(R, "surf_of") else None
        V = topo.verts[i]
        if V.size == 0: continue
        c = V.mean(axis=0); Q = V - c
        _,sv,_ = np.linalg.svd(Q, full_matrices=False)
        L = sv[0]; W = sv[1]
        if topo.planar[i] and W > 1e-12 and L / max(W,1e-12) > 8.0 and topo.areas[i] > 0.2:
            strips.append((L/W, topo.areas[i], np.round(c,1)))
        kinds["planar" if topo.planar[i] else "curved"] += 1
    print(path, dict(kinds), "  thin strips:", len(strips))
    for s in sorted(strips, key=lambda t: -t[1])[:10]:
        print(f"    ratio {s[0]:6.1f} area {s[1]:7.2f} center {s[2]}")
if __name__ == "__main__":
    for p in sys.argv[1:]:
        main(p)
