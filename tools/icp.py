"""Precise alignment between two B-Reps of the same part: inertia axes + ICP."""
import sys, os, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import refit as R
from align import tri_of, best_transform


def closest_on_tris(P, T, chunk=300):
    """(distance, closest point, triangle index) for each point."""
    A, B, C = T[:, 0], T[:, 1], T[:, 2]
    n = len(P)
    out_d = np.full(n, np.inf)
    out_q = np.zeros((n, 3))
    out_i = np.zeros(n, dtype=int)
    for s in range(0, len(T), chunk):
        a, b, c = A[s:s + chunk], B[s:s + chunk], C[s:s + chunk]
        ab, ac = b - a, c - a
        nrm = np.cross(ab, ac)
        nn = np.maximum(np.linalg.norm(nrm, axis=1), 1e-12)
        nu = nrm / nn[:, None]
        w = P[:, None, :] - a[None, :, :]
        h = np.einsum("pij,ij->pi", w, nu)
        proj = P[:, None, :] - h[:, :, None] * nu[None, :, :]
        v0 = proj - a[None, :, :]
        d00 = np.einsum("ij,ij->i", ab, ab)
        d01 = np.einsum("ij,ij->i", ab, ac)
        d11 = np.einsum("ij,ij->i", ac, ac)
        d20 = np.einsum("pij,ij->pi", v0, ab)
        d21 = np.einsum("pij,ij->pi", v0, ac)
        den = np.maximum(d00 * d11 - d01 * d01, 1e-20)
        u = (d11 * d20 - d01 * d21) / den
        v = (d00 * d21 - d01 * d20) / den
        inside = (u >= 0) & (v >= 0) & (u + v <= 1)
        d = np.where(inside, np.abs(h), np.inf)
        Q = np.where(inside[:, :, None], proj, 0.0)
        for p0, p1 in ((a, b), (b, c), (c, a)):
            e = p1 - p0
            ee = np.maximum(np.einsum("ij,ij->i", e, e), 1e-20)
            tt = np.clip(np.einsum("pij,ij->pi", P[:, None, :] - p0[None, :, :], e) / ee, 0, 1)
            q = p0[None, :, :] + tt[:, :, None] * e[None, :, :]
            dq = np.linalg.norm(P[:, None, :] - q, axis=2)
            m = dq < d
            d = np.where(m, dq, d)
            Q = np.where(m[:, :, None], q, Q)
        j = np.argmin(d, axis=1)
        ar = np.arange(len(P))
        dm = d[ar, j]
        m = dm < out_d
        out_d[m] = dm[m]
        out_q[m] = Q[ar, j][m]
        out_i[m] = j[m] + s
    return out_d, out_q, out_i


def kabsch(P, Q):
    """Rigid transform that takes P onto Q (rows = points)."""
    cp, cq = P.mean(axis=0), Q.mean(axis=0)
    H = (P - cp).T @ (Q - cq)
    U, S, Vt = np.linalg.svd(H)
    D = np.diag([1.0, 1.0, np.sign(np.linalg.det(Vt.T @ U.T))])
    Rm = Vt.T @ D @ U.T
    return cp, cq, Rm


def align(shA, shB, defl=0.05, iters=25, cap=2500):
    """Returns a function that maps B's points into A's frame, and the residual."""
    tA, fA, tipA, arA = tri_of(shA, defl)
    tB, fB, tipB, arB = tri_of(shB, defl)
    PA = tA.reshape(-1, 3)
    PB = tB.reshape(-1, 3)
    d0, cA, cB, Rm = best_transform(shA, shB, PA, PB)
    # B's points to align: the triangulation vertices
    S = np.unique(np.round(PB, 6), axis=0)
    if len(S) > cap:
        S = S[np.linspace(0, len(S) - 1, cap).astype(int)]
    # current transform: x_A = cA + (x_B - cB) @ Rm
    off = cA - cB @ Rm
    M = Rm.copy()
    for _ in range(iters):
        X = S @ M + off
        dd, Q, _ = closest_on_tris(X, tA)
        keep = dd <= max(3.0 * np.median(dd), 1e-6)
        cp, cq, Rk = kabsch(X[keep], Q[keep])
        # x' = (x - cp) @ Rk.T + cq
        M2 = M @ Rk.T
        off2 = (off - cp) @ Rk.T + cq
        if np.allclose(M2, M, atol=1e-12) and np.allclose(off2, off, atol=1e-12):
            M, off = M2, off2
            break
        M, off = M2, off2
    X = S @ M + off
    dd, _, _ = closest_on_tris(X, tA)
    return (lambda P: np.atleast_2d(P) @ M + off), float(dd.mean()), float(dd.max()), tA, fA, tipA, arA


if __name__ == "__main__":
    shA = R.read_input(sys.argv[1])[0]
    shB = R.read_input(sys.argv[2])[0]
    f, med, mx, tA, fA, tipA, arA = align(shA, shB)
    print(f"ICP: mean residual {med:.5f} mm, max {mx:.5f} mm")
