"""Which signals tell a CAD EDGE apart from internal tessellation?

Ground truth: the original STEP file, aligned via ICP. A RAW mesh edge is
"CAD" if the two facets sharing it belong to different original faces.
"""
import sys, os, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import refit as R
from align import tri_of
from icp import align, closest_on_tris


def tri_metrics(V):
    """(scale, elongation, dominant direction) of a triangle."""
    e = V - np.roll(V, -1, axis=0)
    L = np.linalg.norm(e, axis=1)
    s = float(np.median(L))
    lmax = float(L.max())
    # minimum height = 2*area/longest_side
    a = 0.5 * float(np.linalg.norm(np.cross(V[1] - V[0], V[2] - V[0])))
    h = 2.0 * a / max(lmax, 1e-12)
    allung = lmax / max(h, 1e-12)
    d = e[int(np.argmax(L))]
    d = d / max(np.linalg.norm(d), 1e-12)
    return s, allung, d, a


def main():
    shA = R.read_input(sys.argv[1])[0]
    shB = R.read_input(sys.argv[2])[0]
    topo = R.Topo(shB)
    toA, med, mx, tA, fA, tipA, arA = align(shA, shB)
    cent = np.array([topo.cents[i] for i in range(topo.nF)])
    dd, _, jj = closest_on_tris(toA(cent), tA)
    of = fA[jj]
    print(f"ICP alignment mean {med:.5f} max {mx:.5f} · {topo.nF} facets")

    nF = topo.nF
    scala = np.zeros(nF); allung = np.zeros(nF); dirs = np.zeros((nF, 3)); ar = np.zeros(nF)
    for i in range(nF):
        scala[i], allung[i], dirs[i], ar[i] = tri_metrics(topo.verts[i])
    N = np.array([topo.norms[i] for i in range(nF)])

    righe = []
    for i in range(nF):
        for j in topo.adj[i]:
            if j <= i:
                continue
            ang = math.degrees(math.acos(float(np.clip(N[i] @ N[j], -1, 1))))
            ds = abs(math.log2(max(scala[i], 1e-12) / max(scala[j], 1e-12)))
            da = abs(math.log2(max(allung[i], 1e-12) / max(allung[j], 1e-12)))
            dA = abs(math.log2(max(ar[i], 1e-12) / max(ar[j], 1e-12)))
            cd = abs(float(dirs[i] @ dirs[j]))
            cad = 1 if of[i] != of[j] else 0
            righe.append((cad, ang, ds, da, dA, cd))
    A = np.array(righe)
    cad = A[:, 0] > 0.5
    print(f"edges: {len(A)}, CAD {cad.sum()} ({100*cad.mean():.1f}%)")
    nomi = ("dihedral", "texture(sides)", "elongation", "area", "cos directions")
    for k, nome in enumerate(nomi, start=1):
        b, d = A[cad, k], A[~cad, k]
        print(f"  {nome:14s} CAD  p10 {np.percentile(b,10):7.3f} p50 {np.median(b):7.3f} "
              f"p90 {np.percentile(b,90):7.3f}")
        print(f"  {'':14s} int. p50 {np.median(d):7.3f} p90 {np.percentile(d,90):7.3f} "
              f"p99 {np.percentile(d,99):7.3f}")
    # how much each criterion catches, at ~1% false positives
    print("\n  criterion                      caught/CAD   false/internal")
    def prova(nome, mask):
        print(f"  {nome:28s} {mask[cad].mean()*100:6.1f}%      {mask[~cad].mean()*100:6.2f}%")
    prova("dihedral > 30", A[:, 1] > 30)
    prova("dihedral > 10", A[:, 1] > 10)
    prova("texture > 1.5", A[:, 2] > 1.5)
    prova("texture > 1.0", A[:, 2] > 1.0)
    prova("elongation > 1.5", A[:, 3] > 1.5)
    prova("elongation > 1.0", A[:, 3] > 1.0)
    prova("area > 2.0", A[:, 4] > 2.0)
    prova("texture>1 or elong>1.5", (A[:, 2] > 1.0) | (A[:, 3] > 1.5))
    prova("dihedral>30 or texture>1.5", (A[:, 1] > 30) | (A[:, 2] > 1.5))
    prova("dihedral>30 or texture>1 or elong>1.5",
          (A[:, 1] > 30) | (A[:, 2] > 1.0) | (A[:, 3] > 1.5))
    prova("cos directions < 0.90", A[:, 5] < 0.90)
    prova("cos directions < 0.50", A[:, 5] < 0.50)
    prova("cos directions < 0.99", A[:, 5] < 0.99)
    prova("dihedral>30 or cos<0.9", (A[:, 1] > 30) | (A[:, 5] < 0.9))
    prova("dihedral>30 or texture>1.5 or cos<0.9",
          (A[:, 1] > 30) | (A[:, 2] > 1.5) | (A[:, 5] < 0.9))
    # the HARD cases: smooth CAD edges
    hard = cad & (A[:, 1] < 10)
    print("")
    print(f"  CAD edges with dihedral < 10 degrees: {hard.sum()} "
          f"({100*hard.sum()/max(cad.sum(),1):.1f}% of the CAD ones)")
    soft = (~cad) & (A[:, 1] < 10)
    def prova2(nome, mask):
        print(f"  {nome:28s} {mask[hard].mean()*100:6.1f}%      {mask[soft].mean()*100:6.2f}%")
    print("  on those (and false positives among smooth internal edges):")
    prova2("texture > 1.5", A[:, 2] > 1.5)
    prova2("texture > 1.0", A[:, 2] > 1.0)
    prova2("cos < 0.90", A[:, 5] < 0.90)
    prova2("cos < 0.50", A[:, 5] < 0.50)
    prova2("texture>1.5 or cos<0.9", (A[:, 2] > 1.5) | (A[:, 5] < 0.9))
    prova2("texture>1.0 or cos<0.9", (A[:, 2] > 1.0) | (A[:, 5] < 0.9))


if __name__ == "__main__":
    main()
