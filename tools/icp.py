"""Precise alignment between two B-Reps of the same part: inertia axes + ICP."""
import sys, os, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import refit as R
from align import tri_of, best_transform


def _pt_tri(P, a, b, c):
    """(distance, closest point) of each point P[i] from the triangle a[i] b[i] c[i]"""
    ab, ac = b - a, c - a
    nrm = np.cross(ab, ac)
    n_len = np.linalg.norm(nrm, axis=1)
    nu = nrm / np.maximum(n_len, 1e-300)[:, None]
    w = P - a
    h = np.einsum("ij,ij->i", w, nu)
    proj = P - h[:, None] * nu
    v0 = proj - a
    d00 = np.einsum("ij,ij->i", ab, ab); d01 = np.einsum("ij,ij->i", ab, ac); d11 = np.einsum("ij,ij->i", ac, ac)
    d20 = np.einsum("ij,ij->i", v0, ab); d21 = np.einsum("ij,ij->i", v0, ac)
    den = np.maximum(d00 * d11 - d01 * d01, 1e-300)
    u = (d11 * d20 - d01 * d21) / den
    v = (d00 * d21 - d01 * d20) / den
    inside = (u >= 0) & (v >= 0) & (u + v <= 1) & (n_len > 1e-12 * np.maximum(d00, d11))
    d = np.where(inside, np.abs(h), np.inf)
    Q = np.where(inside[:, None], proj, 0.0)
    for p0, p1 in ((a, b), (b, c), (c, a)):
        e = p1 - p0
        ee = np.maximum(np.einsum("ij,ij->i", e, e), 1e-300)
        tt = np.clip(np.einsum("ij,ij->i", P - p0, e) / ee, 0, 1)
        q = p0 + tt[:, None] * e
        dq = np.linalg.norm(P - q, axis=1)
        m = dq < d
        d = np.where(m, dq, d)
        Q = np.where(m[:, None], q, Q)
    return d, Q


class TriIndex:
    """
    Exact closest-triangle queries on many triangles: the triangles are sorted by size into classes (each twice the one
    before), a k-d tree of the centroids per class. The closest triangle to a point is no farther by its centroid than
    the best distance known plus its class's size - so a ball per class finds it, and nothing is approximate.
    """

    def __init__(self, T):
        from scipy.spatial import cKDTree
        self.T = T
        self.c = T.mean(axis=1)
        self.r = np.linalg.norm(T - self.c[:, None, :], axis=2).max(axis=1)
        r0 = max(float(np.percentile(self.r, 5)), 1e-9) if len(T) else 1.0
        cls = np.maximum(0, np.ceil(np.log2(np.maximum(self.r, 1e-300) / r0))).astype(int)
        self.classes = []
        for k in np.unique(cls):
            ix = np.nonzero(cls == k)[0]
            self.classes.append((r0 * 2.0 ** k, ix, cKDTree(self.c[ix])))

    def _exact(self, P, pi, ti):
        d = np.empty(len(pi)); Q = np.empty((len(pi), 3))
        for s in range(0, len(pi), 200000):
            sl = slice(s, s + 200000)
            T = self.T[ti[sl]]
            d[sl], Q[sl] = _pt_tri(P[pi[sl]], T[:, 0], T[:, 1], T[:, 2])
        return d, Q

    def near(self, P, margin=0.0):
        """every (point, triangle, distance, closest point) with the triangle within the point's closest distance + margin"""
        P = np.atleast_2d(np.asarray(P, float))
        n = len(P)
        # an upper bound: the triangle of the nearest centroid, per class
        d0 = np.full(n, np.inf)
        ar = np.arange(n)
        for R_, ix, kd in self.classes:
            _, j = kd.query(P)
            d, _ = self._exact(P, ar, ix[j])
            d0 = np.minimum(d0, d)
        pis, tis = [], []
        for R_, ix, kd in self.classes:
            for s in range(0, n, 20000):
                ball = kd.query_ball_point(P[s:s + 20000], d0[s:s + 20000] + margin + R_ * 1.0000001)
                cnt = np.fromiter((len(b) for b in ball), dtype=int, count=len(ball))
                if cnt.sum():
                    pis.append(np.repeat(np.arange(s, s + len(ball)), cnt))
                    tis.append(ix[np.concatenate([np.asarray(b, dtype=int) for b in ball if len(b)])])
        pi = np.concatenate(pis); ti = np.concatenate(tis)
        d, Q = self._exact(P, pi, ti)
        # (d0 was only an upper bound: the margin is from the true closest distance)
        dmin = np.full(n, np.inf)
        np.minimum.at(dmin, pi, d)
        keep = d <= dmin[pi] + margin
        return pi[keep], ti[keep], d[keep], Q[keep]

    def closest(self, P):
        """(distance, closest point, triangle index) for each point"""
        P = np.atleast_2d(np.asarray(P, float))
        pi, ti, d, Q = self.near(P)
        o = np.lexsort((d, pi))
        pi, ti, d, Q = pi[o], ti[o], d[o], Q[o]
        first = np.r_[True, pi[1:] != pi[:-1]]
        out_d = np.full(len(P), np.inf); out_q = np.zeros((len(P), 3)); out_i = np.zeros(len(P), dtype=int)
        out_d[pi[first]] = d[first]; out_q[pi[first]] = Q[first]; out_i[pi[first]] = ti[first]
        return out_d, out_q, out_i


_INDEX = {}


def tri_index(T):
    """the TriIndex of an array of triangles, built once"""
    k = id(T)
    if k not in _INDEX or _INDEX[k][0] is not T:
        if len(_INDEX) > 4:
            _INDEX.clear()
        _INDEX[k] = (T, TriIndex(T))
    return _INDEX[k][1]


def closest_on_tris(P, T, chunk=300):
    """(distance, closest point, triangle index) for each point."""
    if len(T) > 2000 and len(P):
        return tri_index(T).closest(P)
    A, B, C = T[:, 0], T[:, 1], T[:, 2]
    n = len(P)
    out_d = np.full(n, np.inf)
    out_q = np.zeros((n, 3))
    out_i = np.zeros(n, dtype=int)
    for s in range(0, len(T), chunk):
        a, b, c = A[s:s + chunk], B[s:s + chunk], C[s:s + chunk]
        ab, ac = b - a, c - a
        nrm = np.cross(ab, ac)
        n_len = np.linalg.norm(nrm, axis=1)
        nn = np.maximum(n_len, 1e-12)
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
        # (a zero-area triangle has no plane: its "inside" put every point
        # at distance 0 and the ICP "converged" on test14 0.5 mm off)
        inside = (u >= 0) & (v >= 0) & (u + v <= 1) & (n_len > 1e-12 * np.maximum(d00, d11))[None, :]
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
