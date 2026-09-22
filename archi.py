# -*- coding: utf-8 -*-
"""
Quante SPEZZATE dell'uscita sono in realta' un ARCO (o un segmento) del CAD.

Per ogni wire si prendono le corse massimali di spigoli rettilinei consecutivi
e si prova a farci passare: una retta (allora la Fase A non ha unito) oppure
un cerchio. Se passa entro tolleranza, quella corsa e' una curva sola travestita
da spezzata.
"""
import sys, math
sys.path.insert(0, "D:/Users/PC/Desktop/refit")
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
    """(centro, raggio, scarto max) del cerchio ai minimi quadrati per P (n x 3)."""
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
    print(f"--- {path}   diagonale {diag:.1f} mm   tolleranza di giudizio {tol:g} mm")
    tipi = Counter()
    for e in topo.edges:
        tipi[curve_type(e)] += 1
    print(f"    spigoli {len(topo.edges)}: {dict(tipi.most_common())}")

    visti = set()
    restanti = set()
    archi, rette, corse = [], [], 0
    for f in topo.faces:
        for w in R.explore(f, TopAbs_WIRE):
            ex = BRepTools_WireExplorer(TopoDS.Wire_s(w))
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
            # corse massimali di rettilinei consecutivi (sul wire, che e' ciclico)
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
            # ⚠️ UNA CORSA VA SPEZZATA ANCHE AGLI SPIGOLI VIVI. Un anello tutto
            # rettilineo non e' una curva sola: e' una catena di archi e rette
            # attaccati agli angoli. Senza questo taglio si prova a far passare
            # UN cerchio per tutto l'anello, non passa, e si conclude - a torto -
            # che li' non c'era nessun arco.
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
                # retta?
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
                # ⚠️ il cerchio deve spiegare la corsa MOLTO meglio della retta,
                # altrimenti una spezzata quasi dritta passa per un cerchio
                # qualunque e il conto si gonfia di niente.
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
    print(f"    corse di >=3 segmenti consecutivi: {corse}")
    print(f"    ARCHI mancati : {len(archi)}  ({sum(a[0] for a in archi)} spigoli che "
          f"sarebbero {len(archi)} cerchi)")
    if archi:
        r = np.array([a[1] for a in archi]); g = np.array([a[2] for a in archi])
        nn = np.array([a[0] for a in archi])
        print(f"       raggio  min {r.min():.3f}  mediana {np.median(r):.3f}  max {r.max():.3f}")
        print(f"       apertura mediana {np.median(g):.1f} gradi   segmenti per arco "
              f"mediana {int(np.median(nn))}  max {nn.max()}")
        for k, v in sorted(Counter(np.round(r, 3)).most_common(12)):
            pass
        print("       raggi piu' frequenti:",
              dict(Counter(np.round(r, 2)).most_common(8)))
    print(f"    RETTE mancate : {len(rette)}  ({sum(a[0] for a in rette)} spigoli che "
          f"sarebbero {len(rette)} segmenti)")


if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else 1e-3)
