"""
contours.py - the CAD's edges read from the mesh alone.

Every interior edge of the mesh is scored from the tessellation itself - no surface is fitted:
  strong  a crease (over 20 degrees), a jump in facet size (4x, turning half a degree at least),
          a flat patch much wider than what is across the edge (a plane beside a round's strip:
          3x for 4+ triangles, 4x for 3, 5x and 1 degree for a quad), two flat patches turning 2
          degrees where the strips beside don't turn alike (two planes - not a cylinder's quads);
  weak    a change of tessellation pattern (long facets meeting at an angle, long against round,
          a size change, over 10 degrees where neither side's next strip turns alike);
  hint    a weak edge, or a strip turning less for its width than the next one across it (one in
          five is a CAD edge: only where a loose end may go on).
  flat    under 0.05 degrees, or coplanar to float32: never a contour.
A strong contour that stops in the middle goes on along hint edges (over one short unmarked one at
most) to the next contour; then the weak ones join, an edge where nothing changes (under 3 degrees,
alike triangles on both sides: Mesh.tess) is dropped, and the network is pruned of its loose ends
(creases stay). The triangles between the kept edges are grouped into faces; a face of one or two
triangles and under a quarter of its neighbour's area joins the neighbour it shares the flattest
side with. Last, faces next to each other that are one plane, or one cylinder or cone, together
(every vertex within 5 microns of the surface fitted to both) become one: the exporter's vertices lie
on the CAD's surfaces, whatever its slivers between them do. And the other way round: a face no plane,
cylinder or cone fits is split where one grown from its triangles stops, if every part then fits one
(a cylinder and the narrow plane tangent to it, two tangent cylinders), or fits and the rest is
free-form with the triangles changing across the border (a narrow round at a B-spline wall's corner;
cylinders and cones refined by least squares on the vertices); or where single mesh edges from
border to border, the triangles changing across, cut off a piece that fits one. What is left are the mesh's
contours: the lines where the exporter's tessellation says one CAD face ends and another begins - also
between tangent faces, where no angle shows it.
Measured (--cad, 18/11/2026), CAD faces whole by area / CAD edges' length caught / false contours:
test4 97/98%/3 mm, test0 97/98%/16 mm, test10 97/99%/10 mm, test13 91/84%/79 mm, test12 87/88%/0.95 m,
test14 89/90%/3.2 m (the version before: 28-35% whole, 47-92% caught, 6.4 and 19.8 m false on 12 and 14).
Every threshold was tried on all six: what helped four and hurt test12/14 was left out.

  python tools/contours.py MESH.stl                       the network and its faces
  python tools/contours.py MESH.stl --cad ORIG.step       measured against the true CAD (aligned by
                                                          ICP): edges caught / false, faces whole /
                                                          split / merged with another
  python tools/contours.py MESH.stl --step OUT.step       does a conversion's every edge lie on a mesh
                                                          contour? An edge off them is a face boundary
                                                          in the wrong place (a row of facets given to
                                                          the neighbour: test17's arches, 13/11); a
                                                          strong contour no edge covers is a face
                                                          across a crease. Worst ones listed with
                                                          where they are and the faces they part.
  python tools/contours.py MESH.stl --step OUT.step --hints H.json
                                                          the faces found here whole that the
                                                          conversion cut into pieces, as refit's
                                                          hints (refit.py --hints H.json; all three
                                                          steps: tools/refit_contours.py)
  --png FILE [--at X,Y,Z] [--dir DX,DY,DZ] [--span MM]    picture: the contours red (strong thick,
                                                          weak thin), and with --step the output's
                                                          faces by kind (plane grey, cylinder blue,
                                                          cone orange, sphere violet, torus red,
                                                          B-spline green), its edges black
  --near X,Y,Z,R                                          reports restricted to a ball
  --on MM                                                 how close an edge must lie to a contour to
                                                          be on it (default 5 x refit's automatic
                                                          tolerance for the mesh's size)
  --save FILE.npy                                         per triangle: centroid x, y, z, face label
  --json FILE                                             everything measured, for scripts
  (the second positional argument, if a .step, is --cad: the old usage)

As a module:
  import contours as C
  c = C.Contours.from_file("part.stl")      # or C.Contours(topo) on a refit.Topo
  c.labels                                  # face label per triangle
  c.polylines()                             # the contours as chains of mesh vertices
  c.faces()                                 # per face: triangles, area, flat or curved
  c.follow(shape)                           # a B-Rep's edges against the contours
  c.hints(shape)                            # refit's hints where the B-Rep cut a face found here
"""
import sys, os, math, json, argparse
import numpy as np
from collections import defaultdict, Counter

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
import refit as R

R.Log.level = 100


# --- the scores ---------------------------------------------------------------------------------

SHARP = 20.0                    # degrees: an edge turning this much is a contour by itself
FLAT = 0.05                     # degrees: under this two triangles are one plane - never a contour
DSIZE, DSIZE_ANG = 2.0, 0.5     # triangles 2^DSIZE times bigger across the edge, turning by DSIZE_ANG degrees at least
STEP_K, STEP_MIN = 0.7, 0.1     # an edge turning, for its width, less than STEP_K x the next edge across the strip (itself STEP_MIN degrees at least)
PAR, PAR_LEN = 0.94, 0.5        # the edge across the strip: within 20 degrees of this one, and no shorter than half of it (or longer than twice)
CREASE_ALIKE = 0.45             # a weak crease (over 10 degrees): neither side's next edge across turns within 45% of as much
WIDE = 3.0                      # one side's flat patch (4+ triangles) this many times wider than the other's, from the edge
WIDE3 = 4.0                     # ... for a patch of just 3 triangles (3: a strip of test13's hole parted it, 76% -> 90% whole)
PLANE2_RATIO, PLANE2_ANG = 5.0, 1.0  # ... of 2 triangles (a quad): this much wider, and turning this many degrees
FF_ANG, FF_ALIKE = 2.0, 0.6     # two flat patches turning this much, neither side's next edge across within 60% of it
GAP_EDGES, GAP_TURN = 1, 0.5    # a loose end goes on over this many unmarked mesh edges at most, each within 60 degrees of the last
                                # (2: 0.7 m more false contours on the six measured, nothing more caught; 0: test13's rounds lost)
GAP_COST = 4.0                  # an unmarked edge costs this many times a weak one of its length
GAP_LEN = 2.0                   # ... and is no longer than this many times the edge whose end goes on
SAME_TESS, SAME_ANG = 0.2, 3.0  # an edge whose triangles change less than this (Mesh.tess), turning less, is never kept
ABSORB_TRIS, ABSORB_FRAC = 3, 0.25  # a face of fewer triangles joins its neighbour - if this much smaller than it
FIT_TOL = 5e-3                  # mm: two faces next to each other are one surface if all their vertices are this close to it
PLANE_NANG = 0.1                # ... a plane: and each face's mean normal within this many degrees of the plane's
AXIAL_TRIS = 4                  # ... a cylinder or cone: and each face of this many triangles at least
GROW_ANG, GROW_OFF = 10.0, 0.05  # grow_axial(): a neighbour's triangles turn from the surface by this many degrees at
                                #   most - all but this share of its area
SMALL_TOL, SMALL_FOLD, SMALL_VERTS = 1e-3, 5.0, 8   # ... small faces together on one cylinder or cone: mm, degrees, vertices
SCAN_VERTS = 400                # _axial_fit(): the cone's axis scanned (the normals can't tell it) on this many vertices at most
SPLIT_TOL = 1.5e-3              # mm: a face no one surface fits is split if each part fits one this closely
SPLIT_TRIS, SPLIT_SEEDS = 4, 24  # ... a face of this many triangles at least; surfaces grown from this many of its triangles
CHORD_TESS, CHORD_NEAR, CHORD_TRIS = 3.0, 2.0, 4  # split_chords: tess across a chord; its ends not this x its length apart
                                # along the border; the piece cut off of this many triangles at least
SPLIT_TRY = 1                   # ... grown parts tried, biggest first (2: 19 faces more whole, 85 mm more false; 4: 25, 192 mm)
FREE_TESS, FREE_TOL = 1.0, 1e-3  # ... or, the rest fitting nothing (free-form), the part this closely - and the triangles changing
                                # across its border (Mesh.tess, mean by length) this much at least
SPLIT_TRIVIAL, SPLIT_NEAR = 6, 2.0  # ... not if every part has this many vertices or fewer (two strips: three rulings, a
                                # cylinder through any) and the whole face is on one cylinder or cone within this x FIT_TOL
SPLIT_PART = 0.05              # ... the grown part between this fraction of the face's area and 1 minus it
HINT_COVER = 0.95               # hints(): a conversion's face holding this share of ours leaves it alone
HINT_FRAG, HINT_LOOSE = 15, 0.03  # ... a free-form / torus / sphere face only if cut into this many (mm: torus, sphere)
HINT_TRIS = 6
HINT_AXIAL = 3.0                # ... a cylinder or cone within this x FIT_TOL before a torus or sphere within HINT_LOOSE
HINT_EXACT = 1e-3               # ... less the conversion's cylinders, cones, tori its facets' vertices lie on within this (mm)
HINT_DEV = 2.0                  # refit's "max_dev": edge tolerance ceilings the hint's surface may stand off vertices and facets


class Mesh:
    """per-triangle and per-edge features of a refit.Topo (triangles only)"""

    def __init__(self, topo):
        self.topo = topo
        nF = topo.nF
        P = topo.vpos
        self.N = np.array([topo.norms[i] for i in range(nF)])
        # (the triangles' vertices by number, from their edges)
        tv = np.zeros((nF, 3), dtype=int)
        for i in range(nF):
            vs = []
            for k in topo.f_edges[i]:
                for v in topo.e_verts[k]:
                    if v not in vs:
                        vs.append(v)
            tv[i] = (vs + vs[:1] * 3)[:3]
        self.tv = tv
        V = P[tv]
        e = np.roll(V, -1, axis=1) - V
        L = np.linalg.norm(e, axis=2)
        area = 0.5 * np.linalg.norm(np.cross(V[:, 1] - V[:, 0], V[:, 2] - V[:, 0]), axis=1)
        self.area = area
        lmax = L.max(axis=1)
        j = L.argmax(axis=1)
        self.size = np.median(L, axis=1)
        self.elong = lmax / np.maximum(2.0 * area / np.maximum(lmax, 1e-12), 1e-12)
        self.dir = e[np.arange(nF), j] / np.maximum(lmax, 1e-12)[:, None]
        self.E = [k for k, fs in enumerate(topo.e_faces) if len(fs) == 2]
        self.ef = np.array([topo.e_faces[k] for k in self.E]).reshape(-1, 2)
        self.ev = np.array([topo.e_verts[k] for k in self.E]).reshape(-1, 2)
        a, b = self.ef[:, 0], self.ef[:, 1]
        self.ang = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", self.N[a], self.N[b]), -1, 1)))
        self.dsize = np.abs(np.log2(np.maximum(self.size[a], 1e-300) / np.maximum(self.size[b], 1e-300)))
        self.cosd = np.abs(np.einsum("ij,ij->i", self.dir[a], self.dir[b]))
        # how much the triangles themselves change across the edge: size, elongation, smallest angle, and the way the
        # longest side points. Along one CAD face the exporter cuts alike triangles; between two, often not
        cosv = np.einsum("tij,tij->ti", e, -np.roll(e, 1, axis=1)) / np.maximum(L * np.roll(L, 1, axis=1), 1e-300)
        amin = np.arccos(np.clip(cosv, -1, 1)).min(axis=1)
        lel = np.log2(np.maximum(self.elong, 1.0))
        self.tess = (self.dsize + np.abs(lel[a] - lel[b]) + np.abs(amin[a] - amin[b]) / 0.5 + (1.0 - self.cosd))
        self.el_a, self.el_b = self.elong[a], self.elong[b]
        d = P[self.ev[:, 1]] - P[self.ev[:, 0]]
        self.len = np.linalg.norm(d, axis=1)
        self.edir = d / np.maximum(self.len, 1e-300)[:, None]

        # the vertex of each side's triangle that isn't on the edge
        def third(t):
            T = tv[t]
            on = (T == self.ev[:, :1]) | (T == self.ev[:, 1:2])
            return T[np.arange(len(t)), np.argmin(on, axis=1)]

        self.opp = np.column_stack([third(a), third(b)])
        # how wide each side's triangle is from the edge (its third vertex's distance from the edge's line)
        def height(o):
            w = P[o] - P[self.ev[:, 0]]
            return np.linalg.norm(w - np.einsum("ij,ij->i", w, self.edir)[:, None] * self.edir, axis=1)

        self.wid = np.column_stack([height(self.opp[:, 0]), height(self.opp[:, 1])])
        # the turning per unit of width: a strip's chord turns by its width over the radius, half at each side - so
        # inside a face this is the curvature across the edge whatever the strips' widths, and at a tangent edge less
        self.curv = self.ang / np.maximum(self.wid.sum(axis=1), 1e-300)
        # one plane: no angle to speak of, or neither far vertex off the other's plane by more than float32 can tell
        eps = 4.0 * float(np.spacing(np.float32(np.abs(P).max()))) if len(P) else 0.0
        ha = np.abs(np.einsum("ij,ij->i", P[self.opp[:, 1]] - P[self.ev[:, 0]], self.N[a]))
        hb = np.abs(np.einsum("ij,ij->i", P[self.opp[:, 0]] - P[self.ev[:, 0]], self.N[b]))
        self.flat = (self.ang < FLAT) | ((np.maximum(ha, hb) < eps) & (self.ang < SHARP))
        # vertices on the mesh's own border (an open or non-manifold edge): a contour may end there
        self.rim = np.zeros(topo.nV, dtype=bool)
        for k, fs in enumerate(topo.e_faces):
            if len(fs) != 2:
                self.rim[list(topo.e_verts[k])] = True
        self.t_edges = [[] for _ in range(nF)]
        for k in range(len(a)):
            self.t_edges[int(a[k])].append(k); self.t_edges[int(b[k])].append(k)
        self.by_verts = {(min(int(u), int(v)), max(int(u), int(v))): k for k, (u, v) in enumerate(self.ev)}
        (self.far_a, self.fcurv_a), (self.far_b, self.fcurv_b) = self._far(0), self._far(1)
        self._patches()

    def _patches(self):
        """
        The flat patches (triangles in one plane, joined) and, per edge and side, how far that side's patch reaches from
        the edge's line: a plane's face is one patch, as wide as the face; a strip of a curved face is as wide as the
        strip - whatever the triangles each is cut into.
        """
        nF = self.topo.nF
        P = self.topo.vpos
        par = list(range(nF))

        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]; x = par[x]
            return x

        for k in np.nonzero(self.flat)[0]:
            a, b = find(int(self.ef[k, 0])), find(int(self.ef[k, 1]))
            if a != b:
                par[a] = b
        self.patch = np.array([find(i) for i in range(nF)])
        cnt = np.bincount(self.patch, minlength=nF)
        self.pn = np.column_stack([cnt[self.patch[self.ef[:, 0]]], cnt[self.patch[self.ef[:, 1]]]])
        self.pwid = self.wid.copy()
        tris = defaultdict(list)
        for i in np.nonzero(cnt[self.patch] >= 2)[0]:
            tris[int(self.patch[i])].append(int(i))
        sides = defaultdict(list)
        for k in np.nonzero(~self.flat)[0]:
            for s_ in (0, 1):
                if self.pn[k, s_] >= 2:
                    sides[int(self.patch[self.ef[k, s_]])].append((int(k), s_))
        for p, es in sides.items():
            Q = P[np.unique(self.tv[tris[p]])]
            ks = np.array([k for k, _ in es]); ss = np.array([s_ for _, s_ in es])
            n = max(1, 2000000 // max(len(Q), 1))
            for c0 in range(0, len(ks), n):
                kk = ks[c0:c0 + n]
                w = Q[None, :, :] - P[self.ev[kk, 0]][:, None, :]
                t = np.einsum("enj,ej->en", w, self.edir[kk])
                self.pwid[kk, ss[c0:c0 + n]] = np.linalg.norm(w - t[:, :, None] * self.edir[kk][:, None, :], axis=2).max(axis=1)

    def _far(self, side):
        """
        The angle at the edge ACROSS the strip, per edge, from one of its two triangles: the triangle and a neighbour of
        it make a quadrilateral, whose side opposite the edge - the one most parallel to it - is where the surface turns
        next (a cylinder's next ruling, a grid's next line). -1 where there is none.
        """
        out = np.full(len(self.ef), -1.0)
        cur = np.zeros(len(self.ef))
        ang, ev, ef, opp = self.ang, self.ev, self.ef, self.opp
        for k in range(len(ef)):
            t = int(ef[k, side])
            u, v = int(ev[k, 0]), int(ev[k, 1])
            best = PAR
            for d in self.t_edges[t]:
                if d == k or ang[d] >= SHARP:
                    continue
                s2 = 1 if int(ef[d, 0]) == t else 0          # the neighbour's side of d
                w = int(opp[d, s2])
                q = int(ev[d, 0]) if int(ev[d, 0]) not in (u, v) else int(ev[d, 1])
                if q in (u, v) or w in (u, v):
                    continue
                c = self.by_verts.get((min(q, w), max(q, w)))
                if c is None:
                    continue
                par = abs(float(self.edir[c] @ self.edir[k]))
                if par > best and PAR_LEN <= self.len[c] / max(self.len[k], 1e-300) <= 1.0 / PAR_LEN:
                    best = par; out[k] = ang[c]; cur[k] = self.curv[c]
        return out, cur


def score(m):
    """strong / weak marks per interior mesh edge (and m.hint: where a loose end may go on)"""
    # (a crease is a contour whatever turns beside it: "not if the strips on both sides turn as much" - a round in a
    # few strips - took faces away on test10/12/14 and no false contour: test14's false creases at 90-179 degrees are
    # the truth's, slivers along a crease given to the face across)
    sharp = m.ang >= SHARP
    m.sharp = sharp
    strong = sharp | ((m.dsize > DSIZE) & (m.ang > DSIZE_ANG))

    def alike(far, tol):
        """the edge across the strip (on that side) turns within tol x this one's angle: one more step of a run"""
        return (far >= 0) & (np.abs(far - m.ang) <= tol * m.ang)

    # a plane's face beside a strip, or a grid's row: one flat patch much wider than what is across the edge
    lo, hi = m.pwid.min(axis=1), m.pwid.max(axis=1)
    i = np.argmax(m.pwid, axis=1)
    pn_big = m.pn[np.arange(len(i)), i]
    pn_max = m.pn.max(axis=1)
    strong |= (hi >= WIDE * np.maximum(lo, 1e-300)) & (pn_max >= 4)
    strong |= (hi >= WIDE3 * np.maximum(lo, 1e-300)) & (pn_max == 3)
    # ... a plane of just two triangles too (a quad), if wider still and turning from what is across: a round's first
    # strip beside it - test12/14's tangent planes (the six measured: 2762 -> 2839 faces whole, caught 86.4 -> 87.6%,
    # 100 mm more false. The strip's elongation, tried too, told nothing more)
    strong |= (pn_big == 2) & (hi >= PLANE2_RATIO * np.maximum(lo, 1e-300)) & (m.ang > PLANE2_ANG)
    # two flat patches meeting at a few degrees - two planes - unless either side's next edge turns about as much: a
    # cylinder's quads, or a chord wall's planes, turn alike at every strip (the six measured: 614 mm caught, 82 false,
    # mostly test12/14's planes 3-10 degrees apart, under the crease rules)
    strong |= (m.pn.min(axis=1) >= 2) & (m.ang > FF_ANG) & ~(alike(m.far_a, FF_ALIKE) | alike(m.far_b, FF_ALIKE))
    both = (m.el_a >= 3.0) & (m.el_b >= 3.0)
    mixed = ((m.el_a >= 4.0) & (m.el_b < 2.0)) | ((m.el_b >= 4.0) & (m.el_a < 2.0))
    # (over 10 degrees - unless the edge across the strip, on either side, turns about as much: a coarse round turns
    # alike at every strip. Measured on test4/0/10/12/13/14: 21 -> 5.6 m of false contours, 3000 faces whole of 5438
    # (2700 without). "Over 1.5x the next" lost more of test10's chord walls; "alike on both sides" made 4.3 m more
    # false contours on test12/14)
    crease = (m.ang > 10.0) & ~(alike(m.far_a, CREASE_ALIKE) | alike(m.far_b, CREASE_ALIKE))
    weak = (both & (m.cosd < 0.9)) | mixed | (m.dsize > 1.0) | crease

    # (a face's first strip turns from its tangent neighbour by half a step - or by half of each one's step: less, for
    # its width, than the strip turns on its other side. ONLY A HINT: where the strips are cut into triangles between
    # rails that don't match, the angles inside a face are as uneven - measured against the true CAD of test4, test0,
    # test10 and test13, one such edge in five is a CAD edge. It shows a loose end where to go on, and no more)
    def step(far, fcurv):
        return (far >= STEP_MIN) & (far < SHARP) & (m.curv < STEP_K * fcurv)

    strong &= ~m.flat
    weak &= ~strong & ~m.flat
    m.hint = (weak | step(m.far_a, m.fcurv_a) | step(m.far_b, m.fcurv_b)) & ~m.flat
    return strong, weak


def _degrees(m, keep):
    deg = np.zeros(m.topo.nV, dtype=int)
    np.add.at(deg, m.ev[keep, 0], 1); np.add.at(deg, m.ev[keep, 1], 1)
    return deg


def close_gaps(m, keep, lead, hint):
    """
    A contour that stops in the middle of the mesh goes on, along the mesh's edges, to the nearest contour ahead: along
    `hint` edges, over a couple of unmarked ones at most, none of them across a plane, each within 60 degrees of the
    one before. Only the ends of `lead` edges are continued. Returns the edges added.
    """
    import heapq
    P = m.topo.vpos
    keep = keep.copy()
    added = np.zeros(len(keep), dtype=bool)
    at = defaultdict(list)
    for k in np.nonzero(~m.flat)[0]:
        at[int(m.ev[k, 0])].append(int(k)); at[int(m.ev[k, 1])].append(int(k))
    deg = _degrees(m, keep)
    for v0 in [int(v) for v in np.nonzero((deg == 1) & ~m.rim)[0]]:
        if deg[v0] != 1:
            continue
        k0 = next((k for k in at[v0] if keep[k]), None)
        if k0 is None or not lead[k0]:
            continue
        o = int(m.ev[k0, 0]) if int(m.ev[k0, 1]) == v0 else int(m.ev[k0, 1])
        d0 = P[v0] - P[o]; L0 = float(np.linalg.norm(d0)); d0 = d0 / max(L0, 1e-300)
        # (Dijkstra over the vertices; the direction kept at each is the one of the cheapest way in)
        best = {v0: 0.0}
        heap = [(0.0, 0, v0, tuple(d0))]
        prev = {}
        hit = None
        while heap:
            c, n, v, dv = heapq.heappop(heap)
            if c > best.get(v, np.inf):
                continue
            if v != v0 and (deg[v] >= 1 or m.rim[v]):
                hit = v
                break
            for k in at[v]:
                if keep[k]:
                    continue
                w = int(m.ev[k, 0]) if int(m.ev[k, 1]) == v else int(m.ev[k, 1])
                if w == o:
                    continue
                d = P[w] - P[v]; L = float(np.linalg.norm(d)); d = d / max(L, 1e-300)
                turn = float(d @ np.asarray(dv))
                if turn < GAP_TURN:
                    continue
                n2 = n + (0 if hint[k] else 1)
                # (an unmarked step is a short one: not a cylinder's whole ruling)
                if n2 > GAP_EDGES or (not hint[k] and L > GAP_LEN * L0):
                    continue
                c2 = c + L * (2.0 - turn) * (1.0 if hint[k] else GAP_COST)
                if c2 < best.get(w, np.inf):
                    best[w] = c2; prev[w] = (v, k)
                    heapq.heappush(heap, (c2, n2, w, tuple(d)))
        if hit is None:
            continue
        v = hit
        while v != v0:
            v, k = prev[v]
            keep[k] = True; added[k] = True
            deg[int(m.ev[k, 0])] += 1; deg[int(m.ev[k, 1])] += 1
    return added


def close_network(m, keep, fixed):
    """edges stay only where they are part of a network with no loose ends (pruned until none is left); `fixed` ones stay"""
    keep = keep.copy()
    deg = _degrees(m, keep) + 2 * m.rim
    at = defaultdict(list)
    for k in np.nonzero(keep)[0]:
        at[int(m.ev[k, 0])].append(int(k)); at[int(m.ev[k, 1])].append(int(k))
    stack = [int(k) for k in np.nonzero(keep & ~fixed)[0] if deg[m.ev[k, 0]] == 1 or deg[m.ev[k, 1]] == 1]
    while stack:
        k = stack.pop()
        if not keep[k] or fixed[k]:
            continue
        v0, v1 = int(m.ev[k, 0]), int(m.ev[k, 1])
        if deg[v0] != 1 and deg[v1] != 1:
            continue
        keep[k] = False
        for v in (v0, v1):
            deg[v] -= 1
            if deg[v] == 1:
                stack.extend(j for j in at[v] if keep[j] and not fixed[j])
    return keep


def faces_of(m, cut):
    """triangles grouped across every interior edge that isn't cut"""
    nF = m.topo.nF
    par = list(range(nF))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x

    for k in np.nonzero(~cut)[0]:
        a, b = find(int(m.ef[k, 0])), find(int(m.ef[k, 1]))
        if a != b:
            par[a] = b
    return np.array([find(i) for i in range(nF)])


def absorb_small(m, lab, cut, min_tris=ABSORB_TRIS, max_frac=ABSORB_FRAC):
    """
    A face of a few triangles goes, whole, to the neighbour it shares the flattest side with. (Triangle by triangle, the
    two of one face went to two neighbours - and the diagonal between them, in a plane, was a contour)
    Not a big one: a plane of two triangles is a face (test13's 441 mm2 side went into its round) - only one under
    `max_frac` of the neighbour's area goes (56 faces of test4/0/10/13 back whole). ⚠️ min_tris 8 took 100 mm of false
    contours away there (the rings where a mesh's density changes) - and lost 64 faces and 2.5% of the CAD edges' length
    of test12, which it wasn't tuned on: 3 stays.
    """
    lab = lab.copy()
    members = defaultdict(list)
    for i, l in enumerate(lab.tolist()):
        members[l].append(i)
    area = m.area
    farea = defaultdict(float)
    for l, ts in members.items():
        farea[l] = float(area[ts].sum())
    small = sorted((l for l, ts in members.items() if len(ts) < min_tris), key=lambda l: m.size[members[l]].sum())
    for l in small:
        ts = members.get(l)
        if not ts or len(ts) >= min_tris:
            continue
        best = None
        for i in ts:
            for k in m.t_edges[i]:
                j = int(m.ef[k, 1]) if int(m.ef[k, 0]) == i else int(m.ef[k, 0])
                if lab[j] != l and (best is None or m.ang[k] < best[0]):
                    best = (m.ang[k], int(lab[j]))
        if best is None or best[0] >= SHARP or farea[l] > max_frac * farea[best[1]]:
            continue
        to = best[1]
        lab[ts] = to
        members[to].extend(ts)
        farea[to] += farea.pop(l)
        del members[l]
    return lab


def merge_fitted(m, lab, kind, tol=FIT_TOL):
    """
    Faces next to each other that are one plane ("plane") or one cylinder or cone ("axial") together - every vertex of
    both within `tol` of the surface fitted to them - become one face: the exporter's vertices lie ON the CAD's surface,
    what the triangles between them do doesn't matter. A plane cut into flat patches with slivers between (their far
    vertex a micron off the plane, turning a degree) was parted at every sliver: test12, a 558 mm2 plane in 8 faces.
    A plane: each face's own mean normal within PLANE_NANG of the plane's too (without it a narrow round's strip, which
    fits a plane, went into its neighbour: faces lost). A cylinder or cone: both faces of AXIAL_TRIS triangles at least
    (any one or two facets fit a cylinder). The six measured, from 4941 mm of false contours: planes -355 mm, nothing
    caught lost, 10 faces more whole (test10 79 -> 84% of the area, test12 73 -> 76%); cylinders and cones -323 mm
    more, 0.2% of the CAD edges' length lost, test10 99 -> 106 faces whole. (Spheres: 1 face less, 170 mm false less;
    a looser tol, or no minimum of triangles, lost 20-60 faces. A height-field polynomial for the free-form faces of
    test12/14, degree 2-4: 180 mm false less, but test12 0.7% of the area whole less, test4 3 faces with degree 4)
    """
    P = m.topo.vpos
    a, b = m.ef[:, 0], m.ef[:, 1]
    tris = defaultdict(list)
    for i, l in enumerate(lab.tolist()):
        tris[l].append(i)
    verts = {l: np.unique(m.tv[ts]) for l, ts in tris.items()}

    def tilt(x, n):
        ts = np.array(tris[x])
        nm = (m.N[ts] * m.area[ts][:, None]).sum(axis=0)
        nm /= max(np.linalg.norm(nm), 1e-300)
        return math.degrees(math.acos(min(1.0, abs(float(nm @ n)))))

    def fit(la, lb, vs):
        """the largest distance of the vertices from the surface fitted to both faces (inf: not one surface)"""
        V = P[vs]
        if kind == "plane":
            ctr = V.mean(axis=0)
            n = np.linalg.svd(V - ctr, full_matrices=False)[2][2]
            r = float(np.abs((V - ctr) @ n).max())
            return r if r < tol and tilt(la, n) <= PLANE_NANG and tilt(lb, n) <= PLANE_NANG else math.inf
        ts = np.array(tris[la] + tris[lb])
        try:
            p = R.fit_axial(V, m.N[ts], m.area[ts], shallow=True)
            if p is not None:
                p = R.lm_refine(p, V, stall=True)
        except Exception:
            p = None
        return float(np.abs(p.dist(V)).max()) if p is not None else math.inf

    small = 1 if kind == "plane" else AXIAL_TRIS
    pairs = set()
    for k in np.nonzero((lab[a] != lab[b]) & ~m.sharp)[0]:
        la, lb = int(lab[a[k]]), int(lab[b[k]])
        if min(len(tris[la]), len(tris[lb])) >= small:
            pairs.add((min(la, lb), max(la, lb)))
    cand = []
    for la, lb in pairs:
        r = fit(la, lb, np.union1d(verts[la], verts[lb]))
        if r < tol:
            cand.append((r, la, lb))
    if not cand:
        return lab
    cand.sort()     # (the closest fits first; each merge fitted again on the faces as they are by then)
    par = {}

    def find(x):
        while par.get(x, x) != x:
            x = par[x]
        return x

    for _, la, lb in cand:
        ra, rb = find(la), find(lb)
        if ra == rb:
            continue
        vs = np.union1d(verts[ra], verts[rb])
        if fit(ra, rb, vs) >= tol:
            continue
        par[ra] = rb
        verts[rb] = vs
        tris[rb] = tris[rb] + tris[ra]
    if not par:
        return lab
    return np.array([find(int(l)) for l in lab.tolist()])


def _join_small(m, lab, tol=SMALL_TOL, fold=SMALL_FOLD, min_verts=SMALL_VERTS):
    """
    Faces too small to fit on their own (under AXIAL_TRIS triangles), next to each other across folds of `fold`
    degrees at most, become one face where together they are one cylinder or cone - `min_verts` vertices at least,
    all within `tol` - and no plane: a coarse band of quads and fans, each its own flat (test9: 3 faces of 3
    triangles, 192 mm2, on a cylinder within 5e-7 mm).
    """
    P = m.topo.vpos
    a, b = m.ef[:, 0], m.ef[:, 1]
    cnt = defaultdict(int)
    for l in lab.tolist():
        cnt[l] += 1
    small = {l for l, n in cnt.items() if n < AXIAL_TRIS}
    k = np.nonzero((lab[a] != lab[b]) & ~m.sharp & (m.ang < fold))[0]
    par = {l: l for l in small}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for x, y in zip(lab[a[k]].tolist(), lab[b[k]].tolist()):
        if x in small and y in small:
            par[find(x)] = find(y)
    groups = defaultdict(list)
    for l in small:
        groups[find(l)].append(l)
    lab = lab.copy()
    for r, ls in groups.items():
        if len(ls) < 2:
            continue
        g = np.nonzero(np.isin(lab, ls))[0]
        if len(np.unique(m.tv[g])) < min_verts:
            continue
        s = _surfaces(m, g, scan=True)
        if s[0][0] >= FIT_TOL and len(s) > 1 and s[1][0] < tol:
            lab[g] = r
    return lab


def grow_axial(m, lab, tol=FIT_TOL, ang=GROW_ANG):
    """
    A curved face grows over its neighbours on ITS OWN cylinder or cone: the face of AXIAL_TRIS triangles or more
    that one holds within `tol` (and no plane), the largest first, takes each neighbour - of any size - whose
    vertices all lie on that surface within `tol` and whose triangles turn from it by `ang` degrees at most; the
    surface fitted again on what it grew to, from where it was. ⚠️ NOT merge_fitted's pairs: that fits each two
    faces again, the axis read off their normals - on a coarse wall the facets' normals span a few degrees, and the
    axis came out wrong (test9's lettered wall, r 271 mm: one face of 437 mm2 in 31, single fans of 1-90 triangles
    between the letters, each pair "a" cylinder r 444 tilted 18 degrees, or none). The surface of the largest one
    holds them all within 2e-5 mm.
    """
    P = m.topo.vpos
    a, b = m.ef[:, 0], m.ef[:, 1]
    lab = _join_small(m, lab)
    tris = defaultdict(list)
    for i, l in enumerate(lab.tolist()):
        tris[l].append(i)
    nb = defaultdict(set)
    k = np.nonzero((lab[a] != lab[b]) & ~m.sharp)[0]
    for x, y in zip(lab[a[k]].tolist(), lab[b[k]].tolist()):
        nb[x].add(y)
        nb[y].add(x)
    cos_ang = math.cos(math.radians(ang))

    def surface(ts, start=None):
        V = P[np.unique(m.tv[ts])]
        try:
            p = start if start is not None else R.fit_axial(V, m.N[ts], m.area[ts], shallow=True)
            p = R.lm_refine(p, V, stall=True) if p is not None else None
        except Exception:
            return None
        if p is None or getattr(p, "kind", None) != "axial":
            return None
        d = np.abs(p.dist(V))
        return p if np.isfinite(d).all() and float(d.max()) < tol else None

    def lies_on(p, ts):
        ts = np.asarray(ts)
        if float(np.abs(p.dist(P[np.unique(m.tv[ts])])).max()) >= tol:
            return False
        # (by area: a sliver between a letter's outline and a far vertex has its three vertices on the surface and
        # any plane at all around its long side - test9's fans, hundreds of them, kept two faces of one r 404 wall)
        n = p.normal_at(P[m.tv[ts]].mean(axis=1))
        off = np.abs(np.einsum("ij,ij->i", n, m.N[ts])) < cos_ang
        return float(m.area[ts][off].sum()) <= GROW_OFF * float(m.area[ts].sum())

    tnb = defaultdict(list)          # triangle -> its neighbours across edges no crease marks
    for x, y in zip(a[~m.sharp].tolist(), b[~m.sharp].tolist()):
        tnb[x].append(y)
        tnb[y].append(x)

    def take_part(l, p, mine, x):
        """the triangles of a free-form neighbour x (no plane, cylinder or cone holds it within SPLIT_TOL) that lie on
        the surface p within FREE_TOL, joined to the face l, if the triangles change across where they part from the
        rest of x (FREE_TESS, as split_fitted's free-form rest): test9's rim over the lettering, six coarse facets of
        the r 404 wall in one face with the lip's chord wall beside them, 5.4 degrees apart"""
        ts = tris[x]
        if len(ts) < SPLIT_TRIS or fits_any(ts, SPLIT_TOL):
            return None
        ts = np.asarray(ts)
        vs = np.unique(m.tv[ts])
        ok = np.zeros(len(P), bool)
        ok[vs[np.abs(p.dist(P[vs])) < FREE_TOL]] = True
        n = p.normal_at(P[m.tv[ts]].mean(axis=1))
        on = ok[m.tv[ts]].all(axis=1) & (np.abs(np.einsum("ij,ij->i", n, m.N[ts])) >= cos_ang)
        cand = set(ts[on].tolist())
        if not cand or len(cand) == len(ts):
            return None
        mineset = set(mine)
        st = [t for t in cand if any(y in mineset for y in tnb[t])]
        part = set(st)
        while st:
            for y in tnb[st.pop()]:
                if y in cand and y not in part:
                    part.add(y)
                    st.append(y)
        if not part or 2 * len(part) >= len(ts):
            return None
        inP = np.zeros(len(lab), bool)
        inP[list(part)] = True
        inX = lab == x
        bd = (inP[a] != inP[b]) & inX[a] & inX[b]
        if not bd.any() or float((m.tess[bd] * m.len[bd]).sum() / m.len[bd].sum()) < FREE_TESS:
            return None
        return sorted(part)

    def holds(ts):
        """a fit of its own holds the triangles ts too - what _surfaces will read, the hints' kind: test10's face of
        14 triangles took 3 more within 5e-3 of its cylinder, read on its own at 5.5e-3 it went a sphere"""
        u = np.asarray(ts)
        V = P[np.unique(m.tv[u])]
        q2 = _axial_fit(m, u, V)
        return q2 is not None and float(np.abs(q2.dist(V)).max()) < tol

    def fits_any(ts, t):
        return any(r < t for r, _ in _surfaces(m, np.asarray(ts)))

    seeds = []
    for l, ts in tris.items():
        if len(ts) < AXIAL_TRIS or not nb[l]:
            continue
        s = _surfaces(m, ts)
        if s[0][0] < tol:
            continue        # (a plane: merge_fitted's)
        p = surface(ts)
        if p is not None:
            seeds.append((-float(m.area[ts].sum()), l, p))
    seeds.sort(key=lambda q: q[:2])
    par = {}

    def find(x):
        while par.get(x, x) != x:
            x = par[x]
        return x

    for _, l, p in seeds:
        if find(l) != l:
            continue
        mine = list(tris[l])
        todo = set(nb[l])
        seen = {l}
        late, retried = [], False
        while todo or (late and not retried):
            if not todo:
                # ⚠️ ONCE MORE, THE ONES ITS OWN FIT REFUSED, ON ALL IT GREW TO: on part of an arc the axis read off the
                # normals can be anywhere (test9's r 271 wall: a 34 mm2 facet, vertices on it to 1e-5, refused at
                # 1669 triangles - the fit 6.4 mm off - and fine at 2102). Which neighbour came first decided it.
                todo, retried = {x for x in late if find(x) == x}, True
                seen -= todo
            x = todo.pop()
            if x in seen or find(x) != x:
                continue
            seen.add(x)
            if not lies_on(p, tris[x]):
                part = take_part(l, p, mine, x)
                if part and not holds(mine + part):
                    late.append(x)      # (the same partial-arc axis as below: once more at the end)
                elif part:
                    mine += part
                    pset = set(part)
                    tris[x] = [t for t in tris[x] if t not in pset]
                    lab[part] = l
                    # (and the faces the part touches, on the far side of x: theirs to try too)
                    todo |= {find(int(lab[y])) for t in part for y in tnb[t]} - seen
                continue
            q = surface(mine + tris[x], start=p)
            if q is None:
                continue
            if not holds(mine + tris[x]):
                late.append(x)
                continue
            p = q
            par[x] = l
            mine += tris[x]
            todo |= nb[x] - seen
    if not par:
        return lab
    return np.array([find(int(x)) for x in lab.tolist()])


def _axial_fit(m, g, V, tol=FIT_TOL):
    """
    The cylinder or cone through the triangles g's vertices V (refit's fit_axial, refined on V), or None.
    ⚠️ THE NORMALS WEIGHED BY THE TRIANGLE'S SHAPE TOO: a sliver between a letter's outline and a far vertex has
    any plane at all around its long side - test9's r 404 wall, 770 triangles mostly slivers, read as nothing
    (0.22 mm) though its vertices lie on the cylinder within 2e-4. ⚠️ AND WHERE THE NORMALS CAN'T TELL THE AXIS,
    SCANNED: a cone band of a few facets has its normals on a short arc, and the cone's axis is anywhere in the
    plane of their mean and of the axis across their spread (the cylinder's); a few dozen axes in that plane, the
    closest refined (test9's coarse 45-degree chamfers, 19 facets of 147 mm2: 2.5e-2 mm read, 4e-7 scanned).
    Only on a few hundred vertices: the coarse faces it is for.
    """
    T = m.topo.vpos[m.tv[g]]
    e2 = sum(np.einsum("ij,ij->i", T[:, i] - T[:, j], T[:, i] - T[:, j]) for i, j in ((0, 1), (1, 2), (2, 0)))
    q = np.clip(4.0 * math.sqrt(3.0) * m.area[g] / np.maximum(e2, 1e-300), 0.0, 1.0)
    W = m.area[g] * q * q
    if W.sum() <= 1e-12 * max(m.area[g].sum(), 1e-300):
        W = m.area[g]

    def done(p):
        try:
            p = R.lm_refine(p, V, stall=True) if p is not None else None
        except Exception:
            return None, math.inf
        if p is None or getattr(p, "kind", None) != "axial":
            return None, math.inf
        d = np.abs(p.dist(V))
        return (p, float(d.max())) if np.isfinite(d).all() else (None, math.inf)
    try:
        best, r = done(R.fit_axial(V, m.N[g], W, shallow=True))
    except Exception:
        best, r = None, math.inf
    if r < tol or len(V) < 6 or len(V) > SCAN_VERTS:
        return best
    w = W / max(W.sum(), 1e-300)
    N = m.N[g]
    nb = (N * w[:, None]).sum(axis=0)
    D = N - nb
    sp = np.linalg.eigh((D * w[:, None]).T @ D)[1][:, 2]
    a0 = np.cross(nb, sp)
    if np.linalg.norm(nb) < 0.5 or np.linalg.norm(a0) < 1e-9:
        return best
    a0 /= np.linalg.norm(a0)
    u = nb / np.linalg.norm(nb)
    scan = []
    for beta in np.radians(np.linspace(-80.0, 80.0, 33)):
        ax = math.cos(beta) * a0 + math.sin(beta) * u
        try:
            p = R._axial_on_axis(V, ax / np.linalg.norm(ax), False, 3.0)
        except Exception:
            p = None
        if p is not None:
            scan.append((p.rms, p))
    for _, p in sorted(scan, key=lambda s: s[0])[:2]:
        p, rr = done(p)
        if rr < r:
            best, r = p, rr
    return best


def _surfaces(m, g, scan=False):
    """[(largest distance of the triangles g's vertices, distance function)]: the plane and the cylinder or cone
    fitted to them (`scan`: by _axial_fit - slivers weighed down, the cone's axis scanned where the normals can't
    tell it; for the coarse faces and the hints only. split_fitted's seeds read so parted test10's blend
    otherwise, a "sphere" hint of 17 triangles came out of it)"""
    P = m.topo.vpos
    g = np.asarray(g)
    V = P[np.unique(m.tv[g])]
    ctr = V.mean(axis=0)
    n = np.linalg.svd(V - ctr, full_matrices=False)[2][2]
    out = [(float(np.abs((V - ctr) @ n).max()), lambda X, c=ctr, n=n: np.abs((X - c) @ n))]
    if len(g) >= 3:
        if scan:
            p = _axial_fit(m, g, V)
        else:
            try:
                p = R.fit_axial(V, m.N[g], m.area[g], shallow=True)
                if p is not None:     # (fit_axial reads the axis off the normals: on a narrow strip microns off)
                    p = R.lm_refine(p, V, stall=True)
            except Exception:
                p = None
        if p is not None:
            d = np.abs(p.dist(V))
            if np.isfinite(d).all():
                f = lambda X, p=p: np.abs(p.dist(X))
                f.prim = p
                out.append((float(d.max()), f))
    return out


def split_fitted(m, lab, tol=FIT_TOL, part_tol=SPLIT_TOL):
    """
    The other way round from merge_fitted: a face that no plane, cylinder or cone fits (within part_tol) holds two
    CAD faces with no edge marked between them - a cylinder and the narrow plane tangent to it, two tangent cylinders
    of different radius, a cone and its plane. From SPLIT_SEEDS of its triangles a plane and a cylinder or cone are
    fitted to the seed's two rings of neighbours, and grown over every connected triangle whose vertices are within
    `tol` of them; the biggest, fitted again on all it took until it stops growing. If the rest's connected pieces
    each fit a surface too, every triangle goes to the part whose surface its vertices are closest to (the grown
    part's border was a row off where the other surface is close to it: a cone's last facets went to the plane
    beside), the parts fitted again, three times; the face is split if each part then fits its surface within
    part_tol, and the parts are tried again. The six measured, from 4263 mm of false contours: 84 faces more whole,
    0.6% more of the CAD edges' length caught, 22 mm false more (test0 86 -> 96% of the area whole, test10 86 -> 96%,
    test12 76 -> 79%). part_tol = tol: 16 faces more, but 237 mm false (test14's B-splines, flat on one side, fit a
    plane and a cone within 5 microns); the parts not settled by the closest surface: 19 faces less, 80 mm more
    false; the rest not required to fit: 2.6 m false.
    A free-form rest (a B-spline wall and the narrow round at its corner: test12's hexagonal hole, the round one row
    of facets as long as the wall is deep) is split off too if the grown part fits within FREE_TOL and the triangles
    change across its border (FREE_TESS) - unsettled, the rest's surface fits nothing to settle by. With the cylinders
    and cones refined by least squares on the vertices (lm_refine: fit_axial reads the axis off the normals, microns
    off on a strip - the true cylinders of test12 fitted to 1.2 microns median, refined to 0.09), here and in
    merge_fitted: 3087 -> 3144 faces whole, 88.8 -> 89.6% caught, 4285 -> 4224 mm false (test12 79 -> 85% of the
    area whole). (Seeds on the most elongated triangles as well: 7 faces less)
    """
    P = m.topo.vpos
    a, b = m.ef[:, 0], m.ef[:, 1]
    nb = defaultdict(list)          # triangle -> its neighbours in the same face
    for x, y in zip(a.tolist(), b.tolist()):
        if lab[x] == lab[y]:
            nb[x].append(y)
            nb[y].append(x)
    tris = defaultdict(list)
    for i, l in enumerate(lab.tolist()):
        tris[l].append(i)
    lab = lab.copy()
    nxt = int(lab.max()) + 1

    def fits(g, t):
        return any(r < t for r, _ in _surfaces(m, g))

    def trivial(g):
        """a part any surface fits: a strip or two of a coarse wall, as few vertices as the surface has unknowns"""
        return len(np.unique(m.tv[np.asarray(list(g))])) <= SPLIT_TRIVIAL

    def one_wall(ts, parts):
        """the whole face on one cylinder or cone within SPLIT_NEAR x tol (its axis scanned), and nothing changes across
        where the parts meet: alike triangles (Mesh.tess under SAME_TESS, mean by length)"""
        pid = np.full(len(lab), -1)
        for j, c in enumerate(parts):
            pid[list(c)] = j
        bd = (pid[a] >= 0) & (pid[b] >= 0) & (pid[a] != pid[b])
        if not bd.any() or float((m.tess[bd] * m.len[bd]).sum() / m.len[bd].sum()) >= SAME_TESS:
            return False
        s = _surfaces(m, ts, scan=True)
        return len(s) > 1 and s[1][0] < SPLIT_NEAR * tol

    def trivial_split(ts, parts):
        # ⚠️ A FEW STRIPS OF A COARSE WALL FIT A CYLINDER EACH, ALWAYS: test9's rim, one facet tall in strips 5.7 mm
        # wide, on a cylinder r 130 within 5.6e-3 (past tol), came apart in two and two strips - three rulings, three
        # unknowns. Only where the triangles don't change across: test12's tangent planes and cylinders of 5-6 coarse
        # triangles, as few vertices, part where they do (tess 0.65-2.1) - kept split, 4 of its 6 such faces right
        return all(trivial(c) for c in parts) and one_wall(ts, parts)

    def comps(sel):
        """the connected pieces of the triangles sel"""
        sel = set(sel)
        out = []
        while sel:
            s = sel.pop()
            st, c = [s], [s]
            while st:
                for y in nb[st.pop()]:
                    if y in sel:
                        sel.remove(y)
                        st.append(y)
                        c.append(y)
            out.append(c)
        return out

    def piece(s, sel):
        """the connected piece of the triangles sel holding s"""
        sel.discard(s)
        st, c = [s], [s]
        while st:
            for y in nb[st.pop()]:
                if y in sel:
                    sel.remove(y)
                    st.append(y)
                    c.append(y)
        return c

    def grow(seed, ts, vs, reg, only=None):
        """(area, triangles, which surface) - the connected triangles on the seed's side of a surface fitted to
        reg, within tol; the biggest of the plane's and the axial's"""
        best = None
        for j, (r, f) in enumerate(_surfaces(m, reg)):
            if r >= tol or (only is not None and j != only):
                continue
            ok = np.zeros(len(P), bool)
            ok[vs[f(P[vs]) < tol]] = True
            if not ok[m.tv[seed]].all():
                continue
            c = piece(seed, set(ts[ok[m.tv[ts]].all(axis=1)].tolist()))
            ar = float(m.area[c].sum())
            if best is None or ar > best[0]:
                best = (ar, c, j)
        return best

    def settle(ts, groups, iters=3):
        """each triangle to the group whose surface its vertices are closest to, the groups fitted again; their
        connected pieces"""
        ts = np.array(ts)
        cur = np.zeros(len(ts), int)
        pos = {t: i for i, t in enumerate(ts.tolist())}
        for j, g in enumerate(groups):
            for t in g:
                cur[pos[t]] = j
        X = P[m.tv[ts]].reshape(-1, 3)
        for _ in range(iters):
            D = np.full((len(ts), len(groups)), np.inf)
            for j in range(len(groups)):
                if (cur == j).any():
                    f = min(_surfaces(m, ts[cur == j]), key=lambda x: x[0])[1]
                    D[:, j] = f(X).reshape(-1, 3).max(axis=1)
            nw = D.argmin(axis=1)
            if (nw == cur).all():
                break
            cur = nw
        out = [c for j in range(len(groups)) for c in comps(ts[cur == j].tolist())]
        return sorted(out, key=lambda c: -float(m.area[c].sum()))

    def attempt(l, ts, vs, A, best):
        """the face's parts if the part grown `best` parts it (else None)"""
        for _ in range(3):      # fitted again on all it took, until it stops growing
            g = grow(best[3], ts, vs, np.array(best[1]), best[2])
            if g is None or g[0] <= best[0] + 1e-12:
                break
            best = g + (best[3],)
        if not SPLIT_PART * A <= best[0] <= (1 - SPLIT_PART) * A:
            return None
        S = set(best[1])
        rest = comps([t for t in ts.tolist() if t not in S])
        if all(fits(c, tol) for c in rest):
            groups = settle(ts, [sorted(S)] + rest)
            if trivial_split(ts, groups):
                return None
            return groups if all(fits(c, part_tol) for c in groups) else None
        if not fits(np.array(sorted(S)), FREE_TOL):
            return None
        inS = np.zeros(len(lab), bool)
        inS[list(S)] = True
        inF = lab == l
        bd = (inS[a] != inS[b]) & inF[a] & inF[b]
        if not bd.any() or float((m.tess[bd] * m.len[bd]).sum() / m.len[bd].sum()) < FREE_TESS:
            return None
        return None if trivial_split(ts, [sorted(S)] + rest) else [sorted(S)] + rest

    work = [l for l, ts in tris.items() if len(ts) >= SPLIT_TRIS]
    while work:
        l = work.pop()
        ts = np.array(tris[l])
        if fits(ts, part_tol):
            continue
        vs = np.unique(m.tv[ts])
        A = float(m.area[ts].sum())
        seeds = ts if len(ts) <= SPLIT_SEEDS else np.random.default_rng(len(ts)).choice(ts, SPLIT_SEEDS, replace=False)
        cands, seen = [], set()
        for s in seeds:
            s = int(s)
            reg = {s}
            fr = [s]
            for _ in range(2):
                fr = [y for x in fr for y in nb[x] if y not in reg and not reg.add(y)]
            if len(reg) < 3:
                continue
            g = grow(s, ts, vs, np.array(sorted(reg)))
            if g is not None and frozenset(g[1]) not in seen:
                seen.add(frozenset(g[1]))
                cands.append(g + (s,))
        cands.sort(key=lambda x: -x[0])
        groups = None
        for best in cands[:SPLIT_TRY]:  # the biggest part first; the next if it doesn't part the face
            groups = attempt(l, ts, vs, A, best)
            if groups is not None:
                break
        if groups is None:
            continue
        for i, c in enumerate(groups):
            k = l if i == 0 else nxt
            if i:
                nxt += 1
            lab[c] = k
            tris[k] = c
            if i and len(c) >= SPLIT_TRIS:
                work.append(k)
        for t in ts.tolist():
            nb[t] = [y for y in nb[t] if lab[y] == lab[t]]
    return lab


def split_chords(m, lab, tol=FIT_TOL):
    """
    A face's own mesh edges that join two vertices of its border (not CHORD_NEAR x their length apart along it) with
    the triangles changing across (Mesh.tess at least CHORD_TESS) are cut; a piece they cut off that fits a plane,
    cylinder or cone within `tol` - CHORD_TRIS triangles at least, and not one surface with a piece across a chord -
    becomes a face, and what is left its connected pieces. A narrow round tessellated as one row of facets as long as
    the wall is deep, tangent to B-spline walls on both sides (test12's hexagonal hole): its two long sides are each
    ONE mesh edge from the ring of contours in front to the one behind, unmarked (0.5-3 degrees), and its vertices -
    all on the boundary curves, approximated - lie on the cylinder only within 1.3-3.8 microns, too loose for
    split_fitted to grow it apart from the walls' slivers. The six measured: 14 faces more whole, 6 mm more false
    (test12 85 -> 87% of the area whole, its hole's six rounds all whole). (The chords as contours of their own: right
    one time in four, 0.3-0.6 m false; the piece not checked against its neighbour: a cylinder cut in strips by chords
    that all fit it, 0.3 m false)
    """
    import heapq
    a, b = m.ef[:, 0], m.ef[:, 1]
    bd = lab[a] != lab[b]
    on = np.zeros(m.topo.nV, bool)
    on[m.ev[bd].ravel()] = True
    on |= m.rim
    cand = ~bd & on[m.ev[:, 0]] & on[m.ev[:, 1]] & (m.tess >= CHORD_TESS) & ~m.flat
    if not cand.any():
        return lab
    adj = defaultdict(list)             # the border, vertex -> (vertex, length)
    for k in np.nonzero(bd)[0]:
        x, y = int(m.ev[k, 0]), int(m.ev[k, 1])
        adj[x].append((y, float(m.len[k])))
        adj[y].append((x, float(m.len[k])))

    def near(u, v, lim):
        """u reaches v along the border within lim"""
        d = {u: 0.0}
        h = [(0.0, u)]
        while h:
            du, x = heapq.heappop(h)
            if x == v:
                return True
            if du > d.get(x, math.inf):
                continue
            for y, L in adj[x]:
                if du + L <= lim and du + L < d.get(y, math.inf):
                    d[y] = du + L
                    heapq.heappush(h, (du + L, y))
        return False

    byface = defaultdict(list)
    for k in np.nonzero(cand)[0]:
        if not near(int(m.ev[k, 0]), int(m.ev[k, 1]), CHORD_NEAR * m.len[k]):
            byface[int(lab[a[k]])].append(int(k))
    if not byface:
        return lab
    lab = lab.copy()
    nxt = int(lab.max()) + 1

    def pieces(ts, links):
        nb = defaultdict(list)
        for x, y in links:
            nb[x].append(y)
            nb[y].append(x)
        left, out = set(ts), []
        while left:
            s = left.pop()
            st, c = [s], [s]
            while st:
                for y in nb[st.pop()]:
                    if y in left:
                        left.remove(y)
                        st.append(y)
                        c.append(y)
            out.append(c)
        return out

    for l, ks in byface.items():
        ts = np.nonzero(lab == l)[0].tolist()
        inF = np.nonzero((lab[a] == l) & (lab[b] == l))[0]
        cut = set(ks)
        comps = pieces(ts, [(int(a[k]), int(b[k])) for k in inF if k not in cut])
        if len(comps) < 2:
            continue
        cid = {t: i for i, c in enumerate(comps) for t in c}
        across = defaultdict(set)
        for k in ks:
            i, j = cid[int(a[k])], cid[int(b[k])]
            if i != j:
                across[i].add(j)
                across[j].add(i)

        def fits(g):
            return any(r < tol for r, _ in _surfaces(m, np.array(g)))

        off = [i for i, c in enumerate(comps) if len(c) >= CHORD_TRIS and fits(c)
               and not any(fits(c + comps[j]) for j in across[i])]
        if not off or len(off) == len(comps):
            continue
        gone = set()
        for i in off:
            lab[comps[i]] = nxt
            nxt += 1
            gone.update(comps[i])
        rest = [t for t in ts if t not in gone]
        left = pieces(rest, [(int(a[k]), int(b[k])) for k in inF if int(a[k]) not in gone and int(b[k]) not in gone])
        for c in sorted(left, key=lambda c: -float(m.area[c].sum()))[1:]:
            lab[c] = nxt
            nxt += 1
    return lab


def extract(m):
    """(labels, strong, weak, kept) - the whole pipeline on a Mesh"""
    strong, weak = score(m)
    keep = strong.copy()
    keep |= close_gaps(m, keep, strong, m.hint)
    keep |= weak
    # nothing changes across: hardly an angle, and alike triangles on both sides - not a contour, whatever marked it
    # (the six measured: 0.66 m of false contours less - mostly in test12/14's free-form faces - for 0.3% of the CAD
    # edges' length. As a mark of its own "the triangles change here" was right one time in five: the exporter's own
    # density changes inside a B-spline face look the same)
    keep &= ~((m.tess < SAME_TESS) & (m.ang < SAME_ANG) & ~m.sharp)
    cut = close_network(m, keep, m.sharp)
    lab = faces_of(m, cut)
    lab = absorb_small(m, lab, cut)
    lab = merge_fitted(m, lab, "plane")
    lab = merge_fitted(m, lab, "axial")
    lab = grow_axial(m, lab)
    lab = split_fitted(m, lab)
    lab = split_chords(m, lab)
    # (again: a piece split_fitted parted from a face no surface held is on its neighbour's own surface - test9's
    # r 404 wall, 79 triangles of a mixed face of 117 beside the 691 of the rest)
    lab = grow_axial(m, lab)
    return lab, strong, weak, cut


# --- the result -----------------------------------------------------------------------------------

def _seg_dist(P, A, B):
    """distance of each point P[i] from the segment A[i]-B[i]"""
    ab = B - A
    L2 = np.maximum(np.einsum("ij,ij->i", ab, ab), 1e-300)
    t = np.clip(np.einsum("ij,ij->i", P - A, ab) / L2, 0.0, 1.0)
    return np.linalg.norm(P - (A + t[:, None] * ab), axis=1)


class Contours:
    """the contours of a triangle mesh (a refit.Topo of it)"""

    def __init__(self, topo):
        self.topo = topo
        self.m = Mesh(topo)
        self.labels, self.strong, self.weak, self.kept = extract(self.m)
        m = self.m
        # the contours are where the faces part (after the small ones went to their neighbours)
        self.border = self.labels[m.ef[:, 0]] != self.labels[m.ef[:, 1]]
        self.P0 = topo.vpos[m.ev[:, 0]]
        self.P1 = topo.vpos[m.ev[:, 1]]
        self.areas = m.area      # (the facets' own: OCC's per-face area took a third of the time)
        self.cents = np.array([topo.cents[i] for i in range(topo.nF)])
        V = topo.vpos
        self.diag = float(np.linalg.norm(V.max(axis=0) - V.min(axis=0)))
        self._kd = None

    @classmethod
    def from_file(cls, path):
        sh = R.read_input(path)[0]
        return cls(R.Topo(sh))

    # (what it found)
    def summary(self):
        b = self.border
        L = np.linalg.norm(self.P1 - self.P0, axis=1)
        return dict(triangles=int(self.topo.nF), edges=int(len(self.m.E)), strong=int(self.strong.sum()), weak=int(self.weak.sum()),
                    kept=int(self.kept.sum()), faces=int(len(set(self.labels.tolist()))), contour_edges=int(b.sum()),
                    contour_length=float(L[b].sum()), strong_length=float(L[b & self.strong].sum()))

    def polylines(self, near=None):
        """the contours as chains of mesh edges between junctions (a vertex where other than two contour edges meet):
        [{"verts": [...], "points": array, "length", "strong": fraction of the length on strong edges,
          "angle": the mean dihedral angle (degrees), "closed"}]"""
        m = self.m
        ks = np.nonzero(self.border)[0]
        if near is not None:
            c, r = np.asarray(near[:3], float), float(near[3])
            mid = 0.5 * (self.P0[ks] + self.P1[ks])
            ks = ks[np.linalg.norm(mid - c, axis=1) <= r]
        at = defaultdict(list)
        for k in ks:
            at[int(m.ev[k, 0])].append(int(k)); at[int(m.ev[k, 1])].append(int(k))
        used = set()
        out = []

        def walk(v, k):
            vs, es = [v], []
            while True:
                used.add(k); es.append(k)
                a, b = int(m.ev[k, 0]), int(m.ev[k, 1])
                v = b if a == v else a
                vs.append(v)
                if len(at[v]) != 2:
                    return vs, es, False
                nxt = [k2 for k2 in at[v] if k2 != k]
                if not nxt or nxt[0] in used:
                    return vs, es, vs[0] == vs[-1]
                k = nxt[0]

        starts = [v for v in at if len(at[v]) != 2]
        for v in starts:
            for k in at[v]:
                if k not in used:
                    out.append(walk(v, k))
        for k in ks:  # rings with no junction
            if int(k) not in used:
                out.append(walk(int(m.ev[k, 0]), int(k)))
        res = []
        for vs, es, closed in out:
            P = self.topo.vpos[vs]
            Ls = np.linalg.norm(self.P1[es] - self.P0[es], axis=1)
            res.append(dict(verts=vs, points=P, length=float(Ls.sum()), closed=bool(closed),
                            strong=float(Ls[self.strong[es]].sum() / max(Ls.sum(), 1e-300)),
                            angle=float(np.average(m.ang[es], weights=Ls) if Ls.sum() > 0 else 0.0)))
        return res

    def faces(self):
        """per face label: triangles, area, centroid, normal spread (degrees: under 1, flat)"""
        out = {}
        groups = defaultdict(list)
        for i, l in enumerate(self.labels.tolist()):
            groups[l].append(i)
        for l, ix in groups.items():
            ix = np.array(ix)
            w = self.areas[ix]
            n = (self.m.N[ix] * w[:, None]).sum(axis=0)
            n /= max(np.linalg.norm(n), 1e-300)
            spread = float(np.degrees(np.arccos(np.clip(self.m.N[ix] @ n, -1, 1))).max())
            out[int(l)] = dict(triangles=int(len(ix)), area=float(w.sum()), centroid=((self.cents[ix] * w[:, None]).sum(axis=0) / max(w.sum(), 1e-300)).tolist(),
                               spread=spread, flat=spread < 1.0)
        return out

    # (against a B-Rep)
    def _dense(self, step):
        """the contour segments sampled every `step` mm: points, and the segment each came from"""
        ks = np.nonzero(self.border)[0]
        A, B = self.P0[ks], self.P1[ks]
        L = np.linalg.norm(B - A, axis=1)
        n = np.maximum(1, np.ceil(L / step).astype(int))
        seg = np.repeat(ks, n + 1)
        t = np.concatenate([np.linspace(0.0, 1.0, k + 1) for k in n])
        P = np.repeat(A, n + 1, axis=0) + t[:, None] * np.repeat(B - A, n + 1, axis=0)
        return P, seg

    def contour_dist(self, P, step, with_seg=False):
        """the distance of each point from the nearest contour segment (and, with_seg, that segment's index)"""
        from scipy.spatial import cKDTree
        if self._kd is None or self._kd[0] != step:
            D, seg = self._dense(step)
            self._kd = (step, cKDTree(D), seg)
        _, kd, seg = self._kd
        if not len(P):
            return (np.zeros(0), np.zeros(0, dtype=int)) if with_seg else np.zeros(0)
        _, j = kd.query(P, k=min(4, kd.n))
        j = np.atleast_2d(j)
        if j.shape[0] != len(P):
            j = j.T
        d = np.full(len(P), np.inf)
        best = np.zeros(len(P), dtype=int)
        for c in range(j.shape[1]):
            k = seg[j[:, c]]
            dk = _seg_dist(P, self.P0[k], self.P1[k])
            better = dk < d
            d = np.where(better, dk, d)
            best = np.where(better, k, best)
        return (d, best) if with_seg else d

    def follow(self, shape, on=None, near=None, n_worst=15):
        """
        Every edge of the B-Rep sampled and measured from the contours (on: within `on` mm), and every strong contour
        edge from the B-Rep's edges. Returns a dict: per edge (length, max distance, fraction on, midpoint, the kinds of
        the faces it parts), the totals, the worst edges, the strong contours no edge covers.
        """
        from scipy.spatial import cKDTree
        from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
        on = float(on) if on else 5.0 * R.auto_tol(self.diag)
        step = 0.5 * on
        emap = R.edge_face_map(shape)
        fidx = {}
        for n_, f in enumerate(R.explore(shape, R.TopAbs_FACE)):
            fidx.setdefault(f.__hash__(), []).append((f, n_ + 1))

        def face_no(f):
            for f2, n_ in fidx.get(f.__hash__(), ()):
                if f2.IsSame(f):
                    return n_
            return 0

        kinds = {}

        def kind(f):
            n_ = face_no(f)
            if n_ not in kinds:
                kinds[n_] = str(BRepAdaptor_Surface(R.td_Face(f)).GetType()).split("_")[-1].replace("Surface", "")
            return kinds[n_]

        c_near = None if near is None else (np.asarray(near[:3], float), float(near[3]))
        edges, samples = [], []
        for i in range(1, R._size(emap) + 1):
            e = R.td_Edge(emap.FindKey(i))
            if R.bt_Degenerated(e):
                continue
            try:
                cv = BRepAdaptor_Curve(e)
                t0, t1 = cv.FirstParameter(), cv.LastParameter()
                P = np.array([[q.X(), q.Y(), q.Z()] for q in (cv.Value(float(t)) for t in np.linspace(t0, t1, 9))])
                L = float(np.linalg.norm(np.diff(P, axis=0), axis=1).sum())
                n = int(min(400, max(9, math.ceil(L / step) + 1)))
                if n > 9:
                    P = np.array([[q.X(), q.Y(), q.Z()] for q in (cv.Value(float(t)) for t in np.linspace(t0, t1, n))])
                ctype = str(cv.GetType()).split("_")[-1]
            except Exception:
                continue
            if c_near is not None and np.linalg.norm(P - c_near[0], axis=1).min() > c_near[1]:
                continue
            fl = emap.FindFromIndex(i)
            fs = []
            it = R.TopTools_ListIteratorOfListOfShape(fl) if hasattr(R, "TopTools_ListIteratorOfListOfShape") else None
            try:
                fs = [R.td_Face(x) for x in fl]
            except Exception:
                fs = []
                if it is not None:
                    while it.More():
                        fs.append(R.td_Face(it.Value())); it.Next()
            seams = len(fs) == 2 and fs[0].IsSame(fs[1])
            edges.append(dict(curve=ctype, length=L, faces=[face_no(f) for f in fs], kinds=[kind(f) for f in fs], seam=bool(seams), mid=P[len(P) // 2].tolist()))
            samples.append(P)
        # (off a contour by less than a facet: a contour runs alongside and the edge isn't on it - a boundary moved by a
        # row of facets; further: no contour there at all, which may be the detector's miss - long facets meeting at a
        # tangent edge often show nothing)
        facet = float(np.median(np.linalg.norm(self.P1 - self.P0, axis=1)))
        for E, P in zip(edges, samples):
            # ⚠️ THE CONTOUR'S VERTICES ON THE EDGE, NOT ITS CHORDS. An exact curve runs through the mesh's vertices and
            # between two of them stands off their chord by its sag (test17's r 21.5 circles, 0.047 off the 47-gon's
            # sides; a thread's helix more): measured from the chords, every exact edge looked moved. So each sample is
            # matched to its nearest contour segment, and that segment's ends are measured from the edge - ends beyond
            # the edge's own ends left out (the next contour along).
            d, k = self.contour_dist(P, step, with_seg=True)
            Pd = P if len(P) >= 2 else np.vstack([P, P])
            A_, B_ = Pd[:-1], Pd[1:]

            def to_edge(Q):
                """distance of each point from the sampled edge, and whether its nearest point is one of the edge's ends"""
                dd = np.full(len(Q), np.inf)
                tt = np.zeros(len(Q))
                for s_ in range(len(A_)):
                    ab = B_[s_] - A_[s_]
                    L2 = max(float(ab @ ab), 1e-300)
                    t = np.clip(((Q - A_[s_]) @ ab) / L2, 0.0, 1.0)
                    ds = np.linalg.norm(Q - (A_[s_] + t[:, None] * ab), axis=1)
                    better = ds < dd
                    dd = np.where(better, ds, dd)
                    tt = np.where(better, s_ + t, tt)
                end = (tt <= 1e-9) | (tt >= len(A_) - 1e-9)
                return dd, end

            ku = np.unique(k)
            e0, x0 = to_edge(self.P0[ku])
            e1, x1 = to_edge(self.P1[ku])
            # (an end past the edge's own ends doesn't count; one that is a contour going off sideways - its other end
            # on the edge, it is far - is another contour meeting the edge, not this one's)
            e0 = np.where(x0 & (e0 > on), 0.0, e0)
            e1 = np.where(x1 & (e1 > on), 0.0, e1)
            side = (np.minimum(e0, e1) <= on) & (np.maximum(e0, e1) > on)
            off_k = dict(zip(ku.tolist(), np.where(side, 0.0, np.maximum(e0, e1)).tolist()))
            dv = np.array([off_k[int(x)] for x in k])
            dv = np.where(d > facet, np.maximum(dv, d), dv)
            seg_l = np.linalg.norm(np.diff(P, axis=0), axis=1)
            ok = (dv[:-1] <= on) & (dv[1:] <= on)
            d = dv
            E["max"] = float(d.max())
            E["on"] = float(seg_l[ok].sum() / max(seg_l.sum(), 1e-300))
            E["class"] = "on" if E["on"] >= 0.999 else ("moved" if E["max"] <= facet else "unseen")
        real = [E for E in edges if not E["seam"]]
        L_all = sum(E["length"] for E in real)
        L_on = sum(E["length"] * E["on"] for E in real)
        L_moved = sum(E["length"] * (1.0 - E["on"]) for E in real if E["class"] == "moved")
        L_unseen = sum(E["length"] * (1.0 - E["on"]) for E in real if E["class"] == "unseen")
        # (the moved ones first, the longest stretch off first; a closed curve's seam is apart)
        worst = sorted([E for E in real if E["on"] < 0.999], key=lambda E: (E["class"] != "moved", -(1.0 - E["on"]) * E["length"]))[:n_worst]
        # the strong contours no edge of the B-Rep covers
        ks = np.nonzero(self.border & self.strong)[0]
        if c_near is not None:
            mid = 0.5 * (self.P0[ks] + self.P1[ks])
            ks = ks[np.linalg.norm(mid - c_near[0], axis=1) <= c_near[1]]
        unc_len = 0.0
        unc = []
        if len(ks) and samples:
            kd = cKDTree(np.vstack([np.vstack([P[:-1], 0.5 * (P[:-1] + P[1:])]) for P in samples]))
            for k in ks:
                Q = self.P0[k] + np.linspace(0.0, 1.0, 5)[:, None] * (self.P1[k] - self.P0[k])
                dd, _ = kd.query(Q)
                if float(dd.max()) > on + step:
                    Lk = float(np.linalg.norm(self.P1[k] - self.P0[k]))
                    unc_len += Lk
                    unc.append(dict(mid=(0.5 * (self.P0[k] + self.P1[k])).tolist(), length=Lk, angle=float(self.m.ang[k])))
        strong_len = float(np.linalg.norm(self.P1[ks] - self.P0[ks], axis=1).sum()) if len(ks) else 0.0
        unc.sort(key=lambda u: -u["length"])
        return dict(on=on, facet=facet, edges=len(real), seams=len(edges) - len(real), length=L_all, length_on=L_on,
                    length_moved=L_moved, length_unseen=L_unseen,
                    edges_moved=sum(1 for E in real if E["class"] == "moved"), edges_unseen=sum(1 for E in real if E["class"] == "unseen"),
                    edges_off=sum(1 for E in real if E["on"] < 0.999), worst=worst, all=edges,
                    strong_length=strong_len, strong_uncovered=unc_len, uncovered=unc[:n_worst])

    def step_faces(self, shape, defl=0.01, margin=0.03):
        """per triangle: the face of `shape` (a conversion of this mesh, in its frame) it lies on - of the faces within
        `margin` the one whose normal is the triangle's own, then the nearest; -1 off every face"""
        from align import tri_of
        from icp import tri_index
        tA, fA, _, _ = tri_of(shape, defl)
        idx = tri_index(tA)
        nA = np.cross(tA[:, 1] - tA[:, 0], tA[:, 2] - tA[:, 0])
        nA /= np.maximum(np.linalg.norm(nA, axis=1), 1e-300)[:, None]
        pi, ti, d, _ = idx.near(self.cents, margin)
        f = fA[ti]
        o = np.lexsort((d, f, pi)); pi, ti, d, f = pi[o], ti[o], d[o], f[o]
        first = np.r_[True, (pi[1:] != pi[:-1]) | (f[1:] != f[:-1])]; pi, ti, d, f = pi[first], ti[first], d[first], f[first]
        ang = np.degrees(np.arccos(np.clip(np.abs(np.einsum("ij,ij->i", self.m.N[pi], nA[ti])), 0.0, 1.0)))
        o = np.lexsort((d, np.round(ang / 0.2), pi)); pi, f = pi[o], f[o]
        first = np.r_[True, pi[1:] != pi[:-1]]
        out = np.full(len(self.labels), -1); out[pi[first]] = f[first]
        return out

    def hints(self, shape, frag=HINT_FRAG, loose=HINT_LOOSE, min_tris=HINT_TRIS):
        """
        refit's hints (refit.py --hints) where the conversion `shape` cut one of these faces into pieces: the face that
        one plane, cylinder or cone holds within FIT_TOL always (refit fits it again: a wrong one goes back); one only
        a torus or sphere holds within `loose` mm, or none (one free-form patch), only if the conversion cut it into
        `frag` faces or more. ⚠️ ONLY WHERE THE CONVERSION CUT: where its face is larger than ours, ours is the one cut
        (test0: 14 of 14 such faces were halves of a whole CAD face). ⚠️ NOT ON A FEW PIECES: a tangent torus, sphere
        and cylinder this file sees as one face, the conversion parts right into 2-12 analytic faces (test13: 9 free-
        form hints, 57 -> 48 whole faces); a free-form wall it leaves as mesh comes in 15-76 (test12_foro). ⚠️ LESS ITS
        EXACT CURVED PIECES: a free-form patch over the whole face took test18's tube bends - 27 of 104 cylinders and 7
        of 21 tori gone into B-splines, the largest tolerance doubled. The conversion's cylinders, cones, tori whose
        facets' vertices lie on them within HINT_EXACT stay; the hint is what is left, each connected piece still cut
        in two or more. (Not the curved share as a test: test12's B-spline walls, one CAD face each, come in strips of
        cylinders 12-78% of the area, and without those hints its false contours went 2960 -> 4417 mm.)
        Returns [{"type", "label", "points" (centroids), "areas"}].
        """
        sf = self.step_faces(shape)
        from OCP.BRepAdaptor import BRepAdaptor_Surface
        stype = [str(BRepAdaptor_Surface(R.td_Face(f)).GetType()).split("_")[-1] for f in R.explore(shape, R.TopAbs_FACE)]
        curved = [t in ("Cylinder", "Cone", "Torus", "Sphere") for t in stype]
        spline = [t in ("BSplineSurface", "BezierSurface") for t in stype]
        m = self.m
        P = m.topo.vpos
        exact = {}

        def is_exact(p):
            """the conversion's face p is a cylinder, cone, torus or sphere its facets' vertices lie on"""
            if p < 0 or not curved[p]:
                return False
            if p not in exact:
                g = np.nonzero(sf == p)[0]
                ok = False
                if len(g) >= 4:
                    V = P[np.unique(m.tv[g])]
                    ok = min(r for r, _ in _surfaces(m, g)) < HINT_EXACT
                    if not ok:
                        try:
                            q = R.fit_torus(self.cents[g], m.N[g], m.area[g])
                            q = R.lm_refine(q, V, stall=True) if q is not None else None
                            ok = q is not None and float(np.abs(q.dist(V)).max()) < HINT_EXACT
                        except Exception:
                            ok = False
                exact[p] = ok
            return exact[p]

        def kind_of(g):
            """a cylinder or cone within HINT_AXIAL x FIT_TOL, a torus or sphere within `loose`, else freeform"""
            V = P[np.unique(m.tv[g])]
            # ⚠️ THE CYLINDER BEFORE THE SPHERE: a wall one facet tall, a few strips on a cylinder within 5.6e-3 (test9's
            # rim) is a "sphere" within `loose` too - and refit made it one
            ax = _surfaces(m, g, scan=True)
            if len(ax) > 1 and ax[1][0] < HINT_AXIAL * FIT_TOL:
                return "cylinder" if abs(float(ax[1][1].prim.slope)) < 1e-3 else "cone"
            for k, fit in (("torus", lambda: R.fit_torus(self.cents[g], m.N[g], m.area[g])), ("sphere", lambda: R.fit_sphere(V))):
                try:
                    q = fit()
                    q = R.lm_refine(q, V, stall=True) if q is not None else None
                except Exception:
                    q = None
                if q is not None and float(np.abs(q.dist(V)).max()) < loose:
                    return k
            return "freeform"

        def pieces(sel):
            """the connected pieces of the triangles sel (across the mesh's edges)"""
            on = np.zeros(len(self.labels), bool); on[sel] = True
            a, b = m.ef[:, 0], m.ef[:, 1]
            k = on[a] & on[b]
            par = {int(t): int(t) for t in sel}

            def find(x):
                while par[x] != x:
                    par[x] = par[par[x]]; x = par[x]
                return x
            for x, y in zip(a[k].tolist(), b[k].tolist()):
                rx, ry = find(x), find(y)
                if rx != ry:
                    par[rx] = ry
            comp = defaultdict(list)
            for t in sel:
                comp[find(int(t))].append(int(t))
            return [np.array(c) for c in comp.values()]
        def coarse_exact(ts, tol=SMALL_TOL):
            """...but a curved face of a few large facets, SMALL_VERTS vertices on a cylinder or cone within
            `tol` and on no plane, is one surface whatever its count (test9: a wall r 403 mm, 236 mm2 in 5
            triangles, cut in 5 flats)"""
            if len(ts) < 3 or len(np.unique(m.tv[ts])) < SMALL_VERTS - 1:
                return False
            s = _surfaces(m, ts, scan=True)
            return s[0][0] >= FIT_TOL and len(s) > 1 and s[1][0] < tol

        def auto_hint(l, g, s):
            h = dict(type="auto", label=f"contours face {l}", points=self.cents[g].tolist(), areas=self.areas[g].tolist(), max_dev=HINT_DEV)
            if len(s) > 1 and getattr(s[1][1], "prim", None) is not None:
                p = s[1][1].prim      # (refit's start where the facets' normals can't tell the axis - see fit_hint)
                h["surface"] = dict(kind="axial", center=np.asarray(p.center, float).tolist(), axis=np.asarray(p.axis, float).tolist(),
                                    r0=float(p.r0), slope=float(p.slope))
            return h
        out = []
        groups = defaultdict(list)
        for i, l in enumerate(self.labels.tolist()):
            groups[l].append(i)
        for l, ts in groups.items():
            ts = np.asarray(ts)
            # (within FIT_TOL as any hint: test9's rim, 5 triangles on a cone within 2.3e-3 cut in 4 flats)
            if len(ts) < min_tris and not coarse_exact(ts, FIT_TOL):
                continue
            A = float(self.areas[ts].sum())
            pf = np.unique(sf[ts])
            share = np.array([self.areas[ts[sf[ts] == p]].sum() for p in pf])
            if share.max() >= HINT_COVER * A:
                continue
            s = _surfaces(m, ts, scan=True)
            kind = "plane" if s[0][0] < FIT_TOL else "auto" if len(s) > 1 and s[1][0] < FIT_TOL else None
            if kind is None:
                # ⚠️ ...BUT A FREE-FORM FACE IN FREE-FORM PIECES IS CUT, HOWEVER FEW: test9's chamfer and blend
                # patches came in 2-12 B-spline pieces, its coarse blends (6-60 triangles, 25-90 mm2) in planes
                # of a facet or two - 840 mm2, no hint. The conversion's own patches and flats, not the analytic
                # faces test13's tangent fillets part into: those few pieces only if all
                # B-splines or planes, less the exact ones.
                few = len(pf) < frag
                flat_or_free = lambda p: p >= 0 and (spline[int(p)] or stype[int(p)] == "Plane")
                keep = ts[~np.array([is_exact(int(p)) for p in sf[ts]])]
                ps = [int(p) for p in np.unique(sf[keep])]
                if few and not all(flat_or_free(p) for p in ps):
                    # ⚠️ ...LESS THOSE ANALYTIC PIECES, NOT THE WHOLE FACE, WHERE FLATS OF A FACET OR THREE ARE LEFT:
                    # test9's top rims (12-47 mm2) in B-splines and single-facet flats beside one cone or cylinder
                    # the vertices miss by 2.6e-3 - no hint, the flats left. With no such flat the analytic pieces
                    # are the face's (test13's fillet: B-splines, a torus and cylinders - hinted, 59 -> 67 faces)
                    if not any(p >= 0 and stype[p] == "Plane" and int((sf == p).sum()) <= 3 for p in ps):
                        continue
                    big = {p for p in ps if not flat_or_free(p) and int((sf == p).sum()) >= 4}
                    keep = keep[~np.isin(sf[keep], list(big))] if big else keep
                for g in pieces(keep):
                    if len(np.unique(sf[g])) < 2:
                        continue
                    if len(g) >= min_tris:
                        out.append(dict(type=kind_of(g), label=f"contours face {l}", points=self.cents[g].tolist(), areas=self.areas[g].tolist(), max_dev=HINT_DEV))
                    elif coarse_exact(g, FIT_TOL):
                        # ⚠️ ...AND A FEW COARSE FLATS ONE CYLINDER OR CONE HOLDS: test9's top rim, five single-
                        # triangle flats (7 vertices on a cone within 2.1e-3) an exact cylinder parted from its
                        # B-spline - too few for a patch, left as triangles
                        out.append(auto_hint(l, g, _surfaces(m, g, scan=True)))
                continue
            out.append(auto_hint(l, ts, s) if kind == "auto" else
                       dict(type=kind, label=f"contours face {l}", points=self.cents[ts].tolist(), areas=self.areas[ts].tolist(), max_dev=HINT_DEV))
        return out


# --- against the true CAD -------------------------------------------------------------------------

def evaluate(m, lab, of, areas):
    a, b = m.ef[:, 0], m.ef[:, 1]
    true = of[a] != of[b]
    pred = lab[a] != lab[b]
    smooth = m.ang < 10.0
    L = m.len
    out = dict(edges=len(a), cad=int(true.sum()), caught=float((pred & true).sum() / max(true.sum(), 1)), false=float((pred & ~true).sum() / max((~true).sum(), 1)),
               smooth_cad=int((true & smooth).sum()), smooth_caught=float((pred & true & smooth).sum() / max((true & smooth).sum(), 1)))
    # faces: a CAD face is whole if 95% of its area is in one predicted face, and that face is 95% its own
    A_tp = defaultdict(float); A_t = defaultdict(float); A_p = defaultdict(float)
    for i in range(len(lab)):
        A_tp[(int(of[i]), int(lab[i]))] += areas[i]; A_t[int(of[i])] += areas[i]; A_p[int(lab[i])] += areas[i]
    whole = split = merged = 0
    whole_area = 0.0
    best = {}
    for (t, p), ar in A_tp.items():
        if ar > best.get(t, (0, None))[0]:
            best[t] = (ar, p)
    for t, (ar, p) in best.items():
        comp = ar / A_t[t]; pur = ar / A_p[p]
        if comp >= 0.95 and pur >= 0.95:
            whole += 1; whole_area += A_t[t]
        elif comp < 0.95:
            split += 1
        else:
            merged += 1
    out.update(cad_faces=len(A_t), pred_faces=len(A_p), whole=whole, split=split, merged=merged, area=float(sum(A_t.values())), whole_area=whole_area,
               cad_len=float(L[true].sum()), caught_len=float(L[pred & true].sum()), false_len=float(L[pred & ~true].sum()),
               smooth_len=float(L[true & smooth].sum()), smooth_caught_len=float(L[pred & true & smooth].sum()))
    return out


def surf_sig(shA):
    """one id per underlying surface: two CAD faces on the same plane / cylinder / cone / sphere / torus are one surface"""
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    sig, ids = [], {}
    for n, f in enumerate(R.explore(shA, R.TopAbs_FACE)):
        s_ = BRepAdaptor_Surface(R.td_Face(f)); t = str(s_.GetType()).split("_")[-1]
        r = lambda x, d=3: round(float(x), d) + 0.0
        try:
            if t == "Plane":
                p = s_.Plane(); n_ = p.Axis().Direction(); o = p.Location()
                v = np.array([n_.X(), n_.Y(), n_.Z()]); d = float(v @ np.array([o.X(), o.Y(), o.Z()]))
                if v[np.argmax(np.abs(v))] < 0: v, d = -v, -d
                key = (t, r(v[0], 4), r(v[1], 4), r(v[2], 4), r(d))
            elif t in ("Cylinder", "Cone", "Sphere", "Torus"):
                g = getattr(s_, t)(); ax = g.Position() if t != "Sphere" else None
                if t == "Sphere":
                    o = g.Location(); key = (t, r(o.X()), r(o.Y()), r(o.Z()), r(g.Radius()))
                else:
                    dd = ax.Direction(); o = ax.Location(); v = np.array([dd.X(), dd.Y(), dd.Z()]); p0 = np.array([o.X(), o.Y(), o.Z()])
                    if v[np.argmax(np.abs(v))] < 0: v = -v
                    q = p0 - (p0 @ v) * v if t == "Cylinder" else p0
                    rad = (g.Radius(),) if t == "Cylinder" else ((g.RefRadius(), g.SemiAngle()) if t == "Cone" else (g.MajorRadius(), g.MinorRadius()))
                    key = (t, r(v[0], 4), r(v[1], 4), r(v[2], 4), r(q[0]), r(q[1]), r(q[2])) + tuple(r(x) for x in rad)
            else:
                key = ("face", n)
        except Exception:
            key = ("face", n)
        sig.append(ids.setdefault(key, len(ids)))
    return np.array(sig)


def cad_faces(step, shB, topo, m, defl=0.01, margin=0.03):
    """
    Per triangle of the mesh: the face of the true CAD it was cut from (the CAD aligned by ICP), and how far its
    centroid is from it. ⚠️ NOT THE NEAREST FACE: beside a tangent edge both faces are within a chord's sag of the
    centroid, and the nearest one is a coin's toss (the truth of 12/11 was noisy exactly where it mattered). Of the faces
    within `margin` of the nearest, the one whose normal there is the triangle's own - a facet in a plane has the
    plane's normal, the first strip of the round beside it is half a step away; by distance only where the normals
    agree (two patches of one smooth surface).
    Returns (face index per triangle, distance per triangle, the CAD shape).
    """
    from icp import align, tri_index
    shA = R.read_input(step)[0]
    toA, med, mx, tA, fA, tipA, arA = align(shA, shB, defl=defl)
    idx = tri_index(tA)
    nA = np.cross(tA[:, 1] - tA[:, 0], tA[:, 2] - tA[:, 0])
    nA /= np.maximum(np.linalg.norm(nA, axis=1), 1e-300)[:, None]
    C0 = toA(np.array([topo.cents[i] for i in range(topo.nF)]))
    V = toA(topo.vpos)[m.tv]
    nM = np.cross(V[:, 1] - V[:, 0], V[:, 2] - V[:, 0])
    nM /= np.maximum(np.linalg.norm(nM, axis=1), 1e-300)[:, None]
    pi, ti, d, _ = idx.near(C0, margin)
    f = fA[ti]
    # per (triangle, CAD face): the closest of the face's triangles
    o = np.lexsort((d, f, pi))
    pi, ti, d, f = pi[o], ti[o], d[o], f[o]
    first = np.r_[True, (pi[1:] != pi[:-1]) | (f[1:] != f[:-1])]
    pi, ti, d, f = pi[first], ti[first], d[first], f[first]
    ang = np.degrees(np.arccos(np.clip(np.abs(np.einsum("ij,ij->i", nM[pi], nA[ti])), 0.0, 1.0)))
    # (the normals to within 0.2 degrees, then the distance)
    o = np.lexsort((d, np.round(ang / 0.2), pi))
    pi, d, f = pi[o], d[o], f[o]
    first = np.r_[True, pi[1:] != pi[:-1]]
    face = np.zeros(topo.nF, dtype=int); dist = np.full(topo.nF, np.inf)
    face[pi[first]] = f[first]; dist[pi[first]] = d[first]
    return face, dist, shA


def ground_truth(step, shB, topo, m):
    """per triangle: the CAD face it came from - faces on one surface that touch counted as one"""
    face, _, shA = cad_faces(step, shB, topo, m)
    sg = surf_sig(shA)[face]
    par = list(range(topo.nF))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x

    for k in range(len(m.ef)):
        a_, b_ = int(m.ef[k, 0]), int(m.ef[k, 1])
        if sg[a_] == sg[b_]:
            ra, rb = find(a_), find(b_)
            if ra != rb: par[ra] = rb
    return np.array([find(i) for i in range(topo.nF)])


# --- the picture ----------------------------------------------------------------------------------

def render(c, png, at=None, view=(1.0, 1.0, 1.0), span=None, shape=None, size=(1400, 1100)):
    """the contours red (strong thick, weak thin) on the mesh - or on the B-Rep's faces by kind with its edges black"""
    import pyvista as pv
    topo = c.topo
    at = np.asarray(at, float) if at is not None else 0.5 * (topo.vpos.max(axis=0) + topo.vpos.min(axis=0))
    span = float(span) if span else 0.5 * c.diag
    d = np.asarray(view, float); d /= np.linalg.norm(d)
    pl = pv.Plotter(off_screen=True, window_size=size); pl.set_background("white")
    near = lambda P: np.linalg.norm(P - at, axis=1).min() <= 2.5 * span
    if shape is None:
        sel = np.nonzero(np.linalg.norm(c.cents - at, axis=1) <= 2.5 * span)[0]
        V = np.vstack([topo.verts[i] for i in sel]); F = np.hstack([[3, 3 * j, 3 * j + 1, 3 * j + 2] for j in range(len(sel))])
        pl.add_mesh(pv.PolyData(V, F), color="#c8c8c8", show_edges=True, edge_color="#999999", line_width=0.5)
    else:
        from OCP.BRepMesh import BRepMesh_IncrementalMesh
        from OCP.BRep import BRep_Tool
        from OCP.TopLoc import TopLoc_Location
        from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
        col = {"Plane": "#c8c8c8", "Cylinder": "#6aa0e0", "Cone": "#e0b050", "Sphere": "#d070d0", "Torus": "#e06070"}
        BRepMesh_IncrementalMesh(shape, max(1e-3, 1e-3 * span), False, 0.2, True)
        for f in R.explore(shape, R.TopAbs_FACE):
            loc = TopLoc_Location(); tr = R._st(BRep_Tool, "Triangulation")(R.td_Face(f), loc)
            if tr is None:
                continue
            V = np.array([[q.X(), q.Y(), q.Z()] for q in (tr.Node(i).Transformed(loc.Transformation()) for i in range(1, tr.NbNodes() + 1))])
            if not near(V):
                continue
            F = np.array([[3] + [x - 1 for x in tr.Triangle(i).Get()] for i in range(1, tr.NbTriangles() + 1)])
            k = str(BRepAdaptor_Surface(R.td_Face(f)).GetType()).split("_")[-1]
            pl.add_mesh(pv.PolyData(V, F.ravel()), color=col.get(k, "#70c070"), smooth_shading=True)
        for e in R.explore(shape, R.TopAbs_EDGE):
            try:
                a = BRepAdaptor_Curve(R.td_Edge(e))
                P = np.array([[q.X(), q.Y(), q.Z()] for q in (a.Value(float(t)) for t in np.linspace(a.FirstParameter(), a.LastParameter(), 40))])
            except Exception:
                continue
            if near(P):
                pl.add_mesh(pv.lines_from_points(P), color="black", line_width=2)
    for strong, w in ((False, 1.5), (True, 3.5)):
        ks = np.nonzero(c.border & (c.strong == strong))[0]
        if not len(ks):
            continue
        mid = 0.5 * (c.P0[ks] + c.P1[ks])
        ks = ks[np.linalg.norm(mid - at, axis=1) <= 2.5 * span]
        if not len(ks):
            continue
        P = np.empty((2 * len(ks), 3)); P[0::2] = c.P0[ks]; P[1::2] = c.P1[ks]
        L = pv.PolyData(P); L.lines = np.hstack([[2, 2 * i, 2 * i + 1] for i in range(len(ks))])
        pl.add_mesh(L, color="red", line_width=w)
    up = (0, 0, 1) if abs(d[2]) < 0.9 else (0, 1, 0)
    pl.camera_position = [tuple(at + d * 3 * span), tuple(at), up]
    pl.camera.view_angle = 2.0 * math.degrees(math.atan(1.0 / 3.0))
    pl.screenshot(png)


# --- command line ---------------------------------------------------------------------------------

def _vec(s, n):
    v = [float(x) for x in s.replace(";", ",").split(",")]
    if len(v) != n:
        raise argparse.ArgumentTypeError(f"{n} numbers, comma-separated")
    return v


def main(argv=None):
    ap = argparse.ArgumentParser(description="The CAD's edges read from the mesh alone (see the docstring).", formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mesh", help=".stl, or a .step born from a mesh")
    ap.add_argument("orig", nargs="?", help="(old usage) the true CAD, as --cad")
    ap.add_argument("--cad", help="the true CAD (.step): edges and faces measured against it")
    ap.add_argument("--step", help="a conversion (.step): do its edges lie on the contours?")
    ap.add_argument("--on", type=float, help="mm: an edge within this of a contour is on it (default 5 x refit's automatic tolerance)")
    ap.add_argument("--near", type=lambda s: _vec(s, 4), help="X,Y,Z,R: reports restricted to a ball")
    ap.add_argument("--png", help="picture of the contours (and of the --step's faces)")
    ap.add_argument("--at", type=lambda s: _vec(s, 3), help="X,Y,Z: centre of the picture (default: the --near ball's, or the part's)")
    ap.add_argument("--dir", type=lambda s: _vec(s, 3), default=[1.0, 1.0, 1.0], help="DX,DY,DZ: looking from there")
    ap.add_argument("--span", type=float, help="mm: half the picture's height (default: the --near radius, or half the part)")
    ap.add_argument("--hints", help="with --step: refit's hints file (refit.py --hints) for the faces the conversion cut into pieces")
    ap.add_argument("--save", help=".npy: per triangle centroid x, y, z and face label")
    ap.add_argument("--json", help="everything measured")
    ap.add_argument("--worst", type=int, default=15, help="how many edges off the contours to list")
    a = ap.parse_args(argv)
    if a.orig and not a.cad:
        if not a.orig.lower().endswith((".step", ".stp")):
            ap.error("the second argument is the true CAD, a .step")
        a.cad = a.orig
    shB = R.read_input(a.mesh)[0]
    topo = R.Topo(shB)
    c = Contours(topo)
    s = c.summary()
    rep = dict(mesh=os.path.abspath(a.mesh), summary=s)
    print(f"{os.path.basename(a.mesh)}: {s['triangles']} triangles, {s['edges']} edges; strong {s['strong']}, weak {s['weak']}, kept {s['kept']}; "
          f"faces {s['faces']}; contours {s['contour_length']:.1f} mm ({s['strong_length']:.1f} strong)")
    if a.near:
        pls = c.polylines(near=a.near)
        print(f"  near {a.near[:3]} r {a.near[3]}: {len(pls)} contour lines, {sum(p['length'] for p in pls):.2f} mm")
        rep["near_polylines"] = [dict(length=p["length"], strong=p["strong"], angle=p["angle"], closed=p["closed"], ends=[p["points"][0].tolist(), p["points"][-1].tolist()]) for p in pls]
    if a.cad:
        of = ground_truth(a.cad, shB, topo, c.m)
        r = evaluate(c.m, c.labels, of, c.areas)
        print(f"  CAD edges caught {100*r['caught']:.1f}%  false {100*r['false']:.2f}%  | smooth CAD edges ({r['smooth_cad']}) caught {100*r['smooth_caught']:.1f}%")
        print(f"  CAD faces {r['cad_faces']}: whole {r['whole']}, split {r['split']}, merged with another {r['merged']}  (predicted faces {r['pred_faces']})")
        rep["cad"] = r
    shape = None
    if a.step:
        shape = R.read_input(a.step)[0]
        f = c.follow(shape, on=a.on, near=a.near, n_worst=a.worst)
        Lt = max(f['length'], 1e-300)
        print(f"  {os.path.basename(a.step)}: {f['edges']} edges ({f['seams']} seams apart), {f['length']:.1f} mm; on the contours (within {f['on']:.4f} mm) "
              f"{100 * f['length_on'] / Lt:.1f}% of the length")
        print(f"    moved (a contour within a facet, {f['facet']:.3f} mm, and the edge not on it): {f['edges_moved']} edges, {f['length_moved']:.2f} mm ({100 * f['length_moved'] / Lt:.2f}%)")
        print(f"    unseen (no contour there - the detector's miss, or a boundary where the mesh shows none): {f['edges_unseen']} edges, {f['length_unseen']:.2f} mm ({100 * f['length_unseen'] / Lt:.2f}%)")
        if f["strong_length"] > 0:
            print(f"  strong contours no edge covers: {f['strong_uncovered']:.2f} of {f['strong_length']:.1f} mm ({100 * f['strong_uncovered'] / f['strong_length']:.2f}%)")
        if f["worst"]:
            print("  edges off the contours (moved first, the longest stretch off first):")
            for E in f["worst"]:
                print(f"    {E['class']:6s} {E['curve']:12s} L {E['length']:7.3f}  off {100 * (1 - E['on']):5.1f}%  max {E['max']:.4f}  at {np.round(E['mid'], 2).tolist()}  "
                      f"faces {E['faces']} {'/'.join(E['kinds'])}")
        if f["uncovered"]:
            print("  strong contours no edge covers (longest first):")
            for u in f["uncovered"][:a.worst]:
                print(f"    L {u['length']:.3f}  angle {u['angle']:5.1f}  at {np.round(u['mid'], 2).tolist()}")
        rep["follow"] = {k: v for k, v in f.items() if k != "all"}
    if a.hints:
        if shape is None:
            ap.error("--hints needs --step (the conversion to compare)")
        hs = c.hints(shape)
        with open(a.hints, "w", encoding="utf-8") as fh:
            json.dump({"refit_hints": 1, "groups": hs}, fh)
        print(f"  hints: {len(hs)} faces the conversion cut ({', '.join(f'{k} {n}' for k, n in Counter(h['type'] for h in hs).items()) or 'none'}) -> {a.hints}")
    if a.png:
        at = a.at or (a.near[:3] if a.near else None)
        span = a.span or (a.near[3] if a.near else None)
        render(c, a.png, at=at, view=a.dir, span=span, shape=shape)
        print(f"  picture: {a.png}")
    if a.save:
        np.save(a.save, np.column_stack([c.cents, c.labels]))
    if a.json:
        with open(a.json, "w", encoding="utf-8") as fh:
            json.dump(rep, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    return 0


if __name__ == "__main__":
    sys.exit(main())
