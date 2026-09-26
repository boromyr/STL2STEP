# -*- coding: utf-8 -*-
"""
How many POLYLINES in the output are actually an ARC (or a segment) from the CAD.

For every wire, the maximal runs of consecutive straight edges are taken and
we try to fit: a line (meaning Phase A failed to merge it) or a circle. If it
fits within tolerance, that run is a single curve in disguise as a polyline.
"""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import refit as R
from collections import Counter
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCP.TopAbs import TopAbs_WIRE, TopAbs_EDGE
from OCP.BRepTools import BRepTools_WireExplorer
from OCP.TopoDS import TopoDS


def curve_type(e):
    try:
        return str(BRepAdaptor_Curve(e).GetType()).split("_")[-1]
    except Exception:
        return "?"


def fit_circle_3d(P):
    """(center, radius, max deviation) of the least-squares circle for P (n x 3)."""
    c0 = P.mean(axis=0)
    Q = P - c0
    u, s, vt = np.linalg.svd(Q, full_matrices=False)
    n = vt[2]
    piano = float(np.abs(Q @ n).max())
    X, Y = vt[0], vt[1]
    xy = np.c_[Q @ X, Q @ Y]
    A = np.c_[2 * xy, np.ones(len(xy))]
    b = (xy ** 2).sum(axis=1)
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    cen = sol[:2]
    rad2 = sol[2] + cen @ cen
    if not np.isfinite(rad2) or rad2 <= 0:
        return None
    rad = math.sqrt(rad2)
    d = np.abs(np.linalg.norm(xy - cen, axis=1) - rad)
    return c0 + cen[0] * X + cen[1] * Y, rad, max(float(d.max()), piano)


def main(path, tol, ang_vivo=20.0):
    sh = R.read_input(path)[0]
    topo = R.Topo(sh)
    diag = float(np.linalg.norm(topo.vpos.max(axis=0) - topo.vpos.min(axis=0)))
    print(f"--- {path}   diagonal {diag:.1f} mm   judgment tolerance {tol:g} mm")
    tipi = Counter()
    for e in topo.edges:
        tipi[curve_type(e)] += 1
    print(f"    edges {len(topo.edges)}: {dict(tipi.most_common())}")

    visti = set()
    restanti = set()
    archi, rette, corse = [], [], 0
    for f in topo.faces:
        for w in R.explore(f, TopAbs_WIRE):
            ex = BRepTools_WireExplorer(R.td_Wire(w))
            seq = []
            while ex.More():
                e = ex.Current()
                seq.append((e, R.vpos(ex.CurrentVertex())))
                ex.Next()
            if not seq:
                continue
            key = tuple(sorted(topo.emap.FindIndex(e) for e, _ in seq))
            if key in visti:
                continue
            visti.add(key)
            n = len(seq)
            chiuso = n > 2
            lin = [curve_type(e) == "Line" for e, _ in seq]
            # maximal runs of consecutive straight edges (on the wire, which is cyclic)
            if not any(lin):
                continue
            start = 0
            if all(lin):
                runs = [list(range(n))]
            else:
                while lin[start]:
                    start -= 1
                    if start < -n:
                        break
                runs, cur = [], []
                for k in range(n):
                    j = (start + k) % n
                    if lin[j]:
                        cur.append(j)
                    elif cur:
                        runs.append(cur); cur = []
                if cur:
                    runs.append(cur)
            # ⚠️ A RUN MUST ALSO BE SPLIT AT SHARP CORNERS. A ring that is
            # entirely straight is not a single curve: it's a chain of arcs
            # and lines meeting at corners. Without this cut you try to fit
            # ONE circle through the whole ring, it doesn't fit, and you
            # wrongly conclude that there was no arc there at all.
            spezzate = []
            for run in runs:
                if len(run) < 3:
                    continue
                Pr = np.array([seq[j][1] for j in run] + [seq[(run[-1] + 1) % n][1]])
                dd = np.diff(Pr, axis=0)
                ll = np.linalg.norm(dd, axis=1)
                uu = dd / np.maximum(ll, 1e-30)[:, None]
                cc = np.einsum("ij,ij->i", uu[:-1], uu[1:])
                vivo = np.degrees(np.arccos(np.clip(cc, -1, 1))) > ang_vivo
                cur = [run[0]]
                for q in range(1, len(run)):
                    if vivo[q - 1]:
                        spezzate.append(cur); cur = []
                    cur.append(run[q])
                spezzate.append(cur)
            runs = spezzate
            for run in runs:
                if len(run) < 3:
                    continue
                corse += 1
                P = np.array([seq[j][1] for j in run] +
                             [seq[(run[-1] + 1) % n][1]])
                # is it a line?
                d0 = P - P[0]
                dr = P[-1] - P[0]
                L = float(np.linalg.norm(dr))
                if L > 1e-12:
                    dr = dr / L
                    dev = float(np.linalg.norm(d0 - np.outer(d0 @ dr, dr), axis=1).max())
                    if dev <= tol:
                        rette.append((len(run), L))
                        continue
                else:
                    continue
                dev_retta = dev
                c = fit_circle_3d(P)
                if c is None:
                    continue
                cen, rad, dev = c
                # ⚠️ the circle has to explain the run MUCH better than the
                # line, otherwise a nearly-straight polyline fits any old
                # circle and the count inflates for nothing.
                if dev <= tol and rad < 0.5 * diag and dev < 0.2 * dev_retta:
                    V = P - cen
                    ang = 0.0
                    for q in range(len(V) - 1):
                        ang += math.degrees(math.acos(float(np.clip(
                            V[q] @ V[q + 1] / max(np.linalg.norm(V[q]) *
                                                  np.linalg.norm(V[q + 1]), 1e-30), -1, 1))))
                    if ang < 5.0:
                        continue
                    restanti.discard(id(run))
                    archi.append((len(run), rad, ang, dev))
    print(f"    runs of >=3 consecutive segments: {corse}")
    print(f"    MISSED ARCS : {len(archi)}  ({sum(a[0] for a in archi)} edges that "
          f"would be {len(archi)} circles)")
    if archi:
        r = np.array([a[1] for a in archi]); g = np.array([a[2] for a in archi])
        nn = np.array([a[0] for a in archi])
        print(f"       radius  min {r.min():.3f}  median {np.median(r):.3f}  max {r.max():.3f}")
        print(f"       median span {np.median(g):.1f} degrees   segments per arc "
              f"median {int(np.median(nn))}  max {nn.max()}")
        for k, v in sorted(Counter(np.round(r, 3)).most_common(12)):
            pass
        print("       most frequent radii:",
              dict(Counter(np.round(r, 2)).most_common(8)))
    print(f"    MISSED LINES : {len(rette)}  ({sum(a[0] for a in rette)} edges that "
          f"would be {len(rette)} segments)")


if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else 1e-3)
