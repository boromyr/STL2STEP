#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
refit.py
=============================================================================
"Soft" reconversion of a mesh (STL, or STEP born from a mesh) into an analytic B-Rep.

The principle is a single one: EVERY change to the geometry is local, is
verified immediately (solid still closed, faces valid, consistent
orientation, area and volume consistent with the mesh) and if it fails the
check it is DISCARDED. Whatever can't be converted stays tessellated as it
was: the output file is always consistent with the input one, at worst it's
just less "clean". Nothing is sewn (Sewing) and nothing is rebuilt globally.

PHASE A (-a) : merging of coplanar faces and collinear edges.
               Two faces merge only if ALL the vertices of the combined
               group lie within the linear tolerance of the common plane.
PHASE B (-b) : CIRCULAR HOLES. Cylindrical walls closed 360 degrees, concave,
               opening onto two planar faces orthogonal to the axis, with
               both boundaries being exact circles -> an analytic cylinder
               with two circles. Nothing else: it's the cautious phase.
PHASE C (-c) : fillets, chamfers, countersinks, counterbores, corner spheres,
               bosses: every region of facets lying on a cylinder / cone /
               sphere / torus is replaced by the analytic surface, one
               region at a time, with the same check and the same tolerance.
               What Phase B rejects (counterbored holes, holes opening onto
               a curved face) passes here, with the boundaries left
               polygonal where no exact curve describes them.
               Where NO quadric describes the surface - a chamfer running
               along a curved edge, a three-way fillet at a corner - the
               region is rebuilt with a B-spline (see "CURVED FILLETS"
               below). Disabled with --no-free.

THE PHASES RUN IN THE ORDER THEY'RE WRITTEN, and none pulls another one in.
Without flags the sequence is A B C A (the last one re-merges the planar
faces split up by the replacements). You can write whatever you need:
  refit.py x.stl -c                 only Phase C, on the raw mesh
  refit.py x.stl -a -b -c           without the final re-merge
  refit.py x.stl -a -b -c -a 0.005  wide final re-merge (the cleanest)
It's also useful for finding where a defect comes from: if '-c' alone shows
it, it's not Phase A's fault.
⚠️ The final 0.005 mm re-merge doesn't move anything: measured on all the
test files, the maximum distance from the starting part is the same as the
default, while the face count collapses (test4 333 -> 54, test5 217 -> 116).

HOW TO READ THE REPORT (--report)
----------------------------------
One line per region, with the primitive type, the angular coverage, the
facet count, the vertex-to-surface deviation (rms and max) and the outcome:

  [N analytic edges . N reused . N polygonal]
      analytic = boundaries rebuilt with the exact curve (circle, ellipse,
                 intersection of the two surfaces);
      reused   = boundaries already good, taken from the mesh untouched;
      polygonal= boundaries left as the mesh's polyline because no exact
                 curve describes them within tolerance.
  "boundaries with the neighboring planes left polygonal" = second attempt:
      the face is analytic but the boundaries with the neighboring planes
      stay the mesh's.
  "contour left polygonal"               = third attempt: no new edge at
      all, the contour stays identical to the mesh (needed when the
      neighbor is an already-closed analytic face, with its own seam, that
      must not be touched).
  "rejected: ..." explains why: the region stays tessellated. In
      particular "small region isolated among tessellated facets" means it
      was a primitive fitted on noise (random radii in a fillet zone): it's
      left to the mesh, to be cleaned up by hand.

THE CAD's EDGES ARE ALREADY IN THE MESH

A tessellator works one CAD face at a time: inside a face the texture is
uniform, across a CAD edge it changes abruptly, because the two sides were
tessellated by two independent passes with two different curvatures (the
plane in big triangles, the fillet in thin strips). It's visible to the eye
looking at the mesh, and it's measurable: on test4 the facet-size jump is
0.00 (median, in log2) inside a planar face and 1.48 (p90) across different
faces.
Those edges are recognized and used twice.
1. PROTECTED: Phase A never merges across one of them. The rule is
   "texture jump > 1.5 (i.e. almost threefold) AND angle > 0.2 degrees":
   the second condition is what makes it safe, because inside a real
   planar face the facets are coplanar to within 0.04 degrees (p99
   measured) and a texture jump alone, there, is just Delaunay
   triangulation. On test4: zero false cuts, 573 CAD edges protected.
2. AS WORKING BOUNDARIES: cutting the mesh where there's a sharp edge
   (dihedral > 30 degrees) or a texture jump yields the SECTIONS, which
   are the original CAD's faces. Checked against test4's original STEP
   file: the sections cut out exactly the countersink's cone (26.99 mm2),
   the cylinders (105.52 and 60.31) and the free-form surfaces (5.21 and
   5.21) that the quadric fit was covering with five spheres and a torus.
   It's on those sections that the free-form surface is attempted.
   WHAT IT'S WORTH, in numbers, still on test4 against the real CAD: of
   8,523 mesh edges, 1,323 are CAD edges; the rule finds 1,285 of them
   (97.1%) with 18 false positives out of 7,200 internal ones (0.25%).
   Tried and DISCARDED: the triangle-elongation jump (41% caught but 2.8%
   false) and the cosine between dominant directions, which looked very
   strong (inside a face it's 1.000 up to the p99) but brings 4% false
   positives, and the chain filter doesn't remove them because the false
   ones sit in chains too.
3. AND THE BRIDGES. The 38 CAD edges that escape ALL have a single mesh
   edge: they're corner contacts, where two faces touch almost at a point
   (dihedral 1-3 degrees, same texture: geometrically invisible). But
   topologically they're obvious - they're BRIDGES of the facet graph - and
   are found with Tarjan. Without cutting them, twenty CAD faces ended up
   in a single section. Only a bridge separating two parts of at least four
   facets is cut: a sliver attached on one side only isn't a face, and
   detaching it would leave loose mesh around.
   ⚠️ A face's "size" must be measured as the median distance between ALL
   its vertices, not as the length of its sides: a large planar face with a
   finely segmented outline (a plate's shaped edge) has sides as short as a
   facet's and would look as fine as one.
This mostly matters for TANGENT edges (plane-fillet), where the dihedral is
nearly zero and no angular criterion sees them: those are exactly the ones
a loose-tolerance Phase A used to erase, flattening the fillet into the
plane and taking away from Phase C the surface to recognize. Measured
damage, on test4 with Phase A alone at -a 0.005 and then -b -c: 171 faces
without the barrier, 138 with it. Disabled with --no-cad-edges.
Known limit: the signal is texture DENSITY, so it doesn't see a tangent
edge between two curved faces tessellated at the same rate (for those you
need geometry, i.e. the segmenter).

CURVED FILLETS, COUNTERSINKS AND OTHER CONICAL GEOMETRIES

A chamfer along a STRAIGHT edge is a plane, along an ARC it's a cone: both
are recognized. But along an arbitrary edge - a plate's shaped border, the
point where two fillets meet - it's neither one nor the other: it's a
blending surface that the CAD writes as a B-spline, and it's also
recognizable from the tessellation, which there is ten times denser. A
quadric fit can't help but split it into dozens of ten-degree osculating
little cylinders, each with a different radius and a jagged outline: those
are the "faces with too many sides" you see in the finished part. Those
fragments are recognized, grouped by contact and replaced by ONE free-form
surface.
The fit doesn't lower precision, it raises it: a B-spline has as many
degrees of freedom as it needs, and on these patches it reaches a tenth of
the tolerance the little cylinders were barely passing with. The check is
twofold: deviation on the mesh's VERTICES, and deviation on probe points,
i.e. the facet centroids projected onto the primitives of the fragments
being replaced (these are points that lie on the real surface, and are the
only way to look into the gaps BETWEEN vertices, where a spline with too
many poles wobbles). The sparsest pole grid that passes both tests wins; if
none passes, the patch stays tessellated.

HOW CLOSE IT GETS TO THE STARTING CAD (test4, measured)
The original part has 51 faces: 25 planes, 17 cylinders, 6 B-splines, 3
cones. With 'refit.py test4.stl -a -b -c -a 0.005' you get 54 faces: 28
planes, 18 cylinders, 4 B-splines, 4 cones, and the correspondence is
almost always one to one. What still doesn't match are the EDGES: 1,028
against 129, because boundaries that don't come from an exact intersection
stay the mesh's polyline (in the CAD they are 79 lines, 34 circles, 11
B-splines and 5 ellipses). Recognizing that a chain of segments is ONE CAD
curve is the next step.

The cone, moreover, is no longer a second-class surface. Its exact
intersections with what's around it - the counterbore that flares out
(coaxial cylinder), the fillet that closes it off (coaxial torus), the wall
that cuts it (plane: circle, ellipse or hyperbola depending on the tilt) -
are all there now. Without them, a countersink's boundary stayed the mesh's
polyline and, worse, the "circle through the vertices" fallback invented
one that doesn't actually lie on the surface: the pcurves crossed, the face
came out SelfIntersectingWire and the whole countersink went back to being
fully tessellated.

HOW COARSE THE MESH IS (and why it matters)
---------------------------------------------
Tolerances aren't an absolute number: they're scaled to what the mesh
itself can actually tell you.
  - the ACCEPTANCE threshold of a facet (--tol) is the wider of the
    absolute value and HALF the facet's sag, i.e. how far its chord lifts
    off the surface found; it never exceeds the growth tolerance. A facet
    cutting an R1.5 fillet in 22-degree steps is already two hundredths
    below the real surface: demanding its vertices land within a micron
    would be asking the mesh for a precision it doesn't have, and the
    result would be leaving every coarsely-tessellated fillet in strips.
    Above that limit, though, the facet does NOT lie on the surface and
    the region is rejected: that's the boundary between "reconstructing"
    and "deforming".
  - the seed's neighborhood crosses nearly-smooth edges (20 degrees). If
    nothing is found at that threshold, or a primitive is found resting on
    fewer than six facets, the neighborhood reopens up to 32 degrees:
    between two facets of the SAME fillet the angle is just the
    tessellation step, and on a coarse mesh it exceeds 20 degrees. The
    ceiling stays under the 45 degrees of chamfers, which is the case the
    tight threshold guards against.
  - in Phase A two twin strips (the two triangles of a warped
    quadrilateral, coplanar to within a tenth of a micron but separated by
    a few thousandths of a degree) get merged: without this, every
    tessellated fillet comes out with double the faces.
  - once the search is done, two RECOVERY passes put back what the seed
    order left behind: every region tries again to grow onto the facets
    still unclaimed, and clusters of loose facets fully surrounded by
    already-rebuilt surfaces are handed to whichever neighboring region
    explains them best. These are the slivers you used to see next to a
    hole, chamfer or fillet.

SPEED (-j N)
------------
The primitive search from the seeds, which is the biggest slice of Phases B
and C, runs on multiple PROCESSES: -j N uses up to N of them (default: the
machine's cores minus two, capped at 22). Not threads, because both
OpenCascade and the numpy code hold the GIL and threads gain nothing. The
number of processes is then scaled to the work at hand: on a small part,
spinning up twenty costs more than it saves. Only Phase A and the
in-solid replacements stay on a single core: they work on a shared
OpenCascade structure and must be checked one at a time.
With more processes the seeds are tried against a snapshot of the facets
already claimed, so the set of regions found can differ slightly from
-j 1; every region is still checked and rejected with the same criteria.
For a 100% reproducible result, use -j 1.

ONE CAD FACE, ONE FACE OF OURS
-------------------------------
Growth starts from different seeds, and a long wall often ends up as two or
three regions with the exact same surface: in the file they become two or
three faces separated by edges that don't exist in the CAD, with the
facets neither one claimed left in between - the "slivers". That's why
neighboring regions lying on the SAME surface are merged and refitted.
⚠️ THREE THINGS THAT LOOK LIKE DETAILS AND AREN'T.
 1. The merge runs TWICE: after the seeds, and again at the end. Two pieces
    separated by a few loose facets don't even touch, and pass the first
    round unscathed: they only become neighbors after the loose facets are
    recovered.
 2. Comparing the two primitives is a FILTER, not the verdict. Requiring
    axes within 0.05 degrees makes no sense for two pieces of the same
    wall: a cylinder fitted on ten facets spanning twenty degrees has a
    poorly conditioned fit (radius and center trade off against each
    other) and an axis uncertain by a few tenths of a degree. It's widened
    to one degree and 2% of the radius, and the union's refit, measured
    with the tol_grow yardstick - already the definition of "close" that
    the regions grew by - has the final say. A real step between two bores
    (two hundredths) blows through it sixfold and stays two faces.
 3. The union must not PINCH: if the merged boundary has an inner ring
    touching the outer one at a vertex, BRepCheck accepts it in memory but
    after the STEP write the face becomes "IntersectingWires" and the part
    is no longer valid. Those are left separate.
Measured on test4 against the real CAD file: the R=3 cylinder came out
split into two faces of 10.5 and 2.4 mm2 with three planar slivers in
between, where the CAD has ONE face of 13.697 mm2. Now it's one face of
13.722 and the slivers are gone (55 faces against the CAD's 51, it used to
be 59).

THE VOLUME CEILING ALSO COUNTS THE NEIGHBORS
-----------------------------------------------
Every replacement goes through a volume check. The ceiling used to be
computed on the region's area alone, but replacing a region doesn't just
move its own face: the polygonal boundaries touching it become curves and
the neighboring faces get rebuilt along with them. A spherical cap R0.5 of
half a mm2 wedged between three faces of 1, 5 and 12 mm2 overshot by a
hair (+0.0555 against 0.0513) and stayed tessellated - while two IDENTICAL
caps, elsewhere, passed. That's exactly the "identical geometry, different
results" you see when you open the file.
⚠️ On test8 the check was rejecting 66 conversions and NONE of them by more
than four times the ceiling: it wasn't catching real errors anymore, only
good conversions. The new ceiling is a TRUE bound, not an estimate - if no
face moves more than dev, the enclosed volume can't change by more than
the touched area times dev - and it never tightens the previous one.

AND THEN THE ARCS (last step, --no-arcs to turn it off)
----------------------------------------------------------
Once the phases are done, the FINISHED B-Rep is examined and every chain
of straight segments is asked: are you sitting on a circle? This is needed
because the engine only rebuilds with the exact curve the boundaries whose
intersection it can compute, and leaves everything else as it was - on
test8, 4,862 boundaries reused against 328 rebuilt. So a cylinder and the
plane it emerges from were still touching along a six-segment polyline.
Measured on test8: of the 249 replaceable chains (interior vertices of
degree 2, always the same two faces on either side) 110 sit on a circle to
within a MICRON, and the radii are the ones from the drawing - R0.5 fifty-
three times, R0.3 eighteen, R0.8 eight. It's not an approximation: it's the
CAD's real edge that the tessellation had split up, and putting it back
brings the part closer to the original.
Every arc is checked on BOTH faces touching it before being accepted, and
at the end the whole part is rechecked: if it doesn't hold up, all of them
are given up at once.
  part     edges             arcs    volume
  test8    5,962 -> 5,687     51     +0.0003%   (circles from 277 to 328)
  test4   1,039 ->   819     16     +0.0033%
  test1      666 ->   548     19     +0.0038%
  test5    1,106 ->   927     10     -0.0015%
  test2    1,203 -> 1,106      6     -0.0074%
No face moves, fidelity doesn't change (test8: mean deviation from the mesh
0.00082 mm before and after, maximum 0.044).
⚠️ TWO TRAPS, both of which cost hours. (1) The SVD's frame isn't
necessarily right-handed: gp_Ax2(P, N, X) measures the angle from X toward
N x X, so using vt[1] as the Y axis, half the arcs come out mirrored. (2)
MakeEdge(circle, V1, V2) takes the arc going from V1 to V2 at INCREASING
parameter: swapping the two vertices without also flipping the circle
takes the COMPLEMENTARY arc - 300 degrees instead of 60 - which in (u,v)
exits on the other side of the surface and produces "UnorientableShape".

THE TWO TOLERANCES, WHICH DON'T TALK TO EACH OTHER
------------------------------------------------------
There are two, independent, and they govern different things. Confusing
them is the easiest way to get a worse result while believing you've
improved it.

1. PHASE A'S TOLERANCE - the number after -a, or --lin-tol.
   How far the vertices of two faces can be from the common plane for them
   to be MERGED INTO ONE. It only concerns geometry that's already planar:
   it changes no fit. Default of the initial pass: 2e-6 x diagonal.

2. PHASES B AND C's TOLERANCE - --tol.
   How far a mesh vertex can be from the surface for the primitive to be
   ACCEPTED. It's the only tolerance of the two curved phases: all the
   others derive from it.
       region growth                  10 x tol   (never beyond 1e-3 x diagonal)
       intersection curves             4 x tol
       "same surface as the neighbor" 50 x tol
       edge-tolerance ceiling         20 x tol   (--max-edge-tol)
   Default: 1e-5 x diagonal, minimum 2e-4. On a part with a 148 mm
   diagonal that's 1.5e-3.
   ⚠️ WIDENING IT DOESN'T RECOGNIZE MORE: IT RECOGNIZES LESS.
   It's the opposite of what intuition says, and it's measured. test9, a
   coarse mesh (26,384 triangles over 112 mm), always
   'refit.py test9.stl -a -b -c -a 0.05', changing only --tol:

     --tol     C regions  cylinders  spheres  faces  >8 sides  mean dev    max
     0.02          97        41       15   2,444    112      0.00121 mm  0.496
     0.0015 (def) 390       304       39   1,281    139      0.00131 mm  0.150
     0.0005       413       346       20   1,594     95      0.00104 mm  0.150
     0.0002       389       342        7   1,881     94      0.00084 mm  0.053

   (dev = distance of the output's points from the mesh, over 20,000 samples.)
   At 0.02 the cylinders collapse from 304 to 41: with a wide threshold the
   curved band passes the PLANE test, gets taken as flat and never becomes
   a curved region. And the few that remain overrun past the CAD edge,
   pull in facets that don't belong to them and then the fit doesn't
   close: the whole region is lost. Peak recognition sits AROUND 5e-4, and
   below that you only lose what was fitted on noise (spheres from 39 down
   to 7, while the maximum deviation drops from 15 to 5 hundredths of a
   millimeter). You widen --tol only to accept a deliberate deformation,
   never to "catch more stuff".

WHEN A CURVE ISN'T RECOGNIZED, IT'S OFTEN THE MESH
------------------------------------------------------
On test9, 25% of the area still comes out tessellated. Looking at where:
619 of the 887 remaining faces are under half a square millimeter and add
up to 0.6% of the area (they're the engraved lettering, and there's
nothing to recognize there). The real area sits in a few large walls with
very little curvature, and there the problem is the mesh's step: a 520 mm2
band covering 11 degrees of wall is made of SIX chords and has seven
distinct vertices in cross-section. The best cylinder passes 8.3e-02 mm
from the vertices, the best cross-section circle (R 287) at 6.3e-02: the
cross-section isn't an arc, and with seven points it's indistinguishable
from anything else. Catching it would require raising --tol to 0.02-0.05,
i.e. accepting moving the wall by five to eight hundredths. This isn't
done by default: a tessellated wall in the right spot beats a smooth wall
that's been moved.

USAGE
-----
    python refit.py part.stl                  # A + B + C + final A re-merge
    python refit.py part.stl -a               # Phase A only
    python refit.py part.stl -a 0.01          # Phase A only, 10-micron merge
    python refit.py part.stl -b               # A + holes
    python refit.py part.stl -c               # A + fillets/machining
    python refit.py part.stl -b -c --report   # everything, with a .txt report
    python refit.py part.stl -b -c -a 0.01    # everything, final merge at 10 um
    python refit.py part.stl -b -c -i 2       # two B/C cycles before saving
    python refit.py part.stl -j 22            # up to 22 worker processes
    python refit.py part.stl -j 1             # everything sequential, reproducible
    python refit.py part.step -o out.step     # STEP input
    python refit.py --check                   # check the environment

DEPENDENCIES
------------
    pip install cadquery-ocp numpy       (OCP, recommended on Windows)
    or conda install -c conda-forge pythonocc-core numpy
"""

from __future__ import annotations

import argparse
import importlib
import math
import os
import sys
import time
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# =============================================================================
# 0. LOGGER
# =============================================================================


class _C:
    RESET = "\033[0m"
    DIM = "\033[2m"
    BOLD = "\033[1m"
    GREY = "\033[90m"
    RED = "\033[91m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    BLUE = "\033[94m"
    MAGENTA = "\033[95m"
    CYAN = "\033[96m"


def _enable_vt_windows() -> None:
    if os.name != "nt":
        return
    try:
        import ctypes

        k = ctypes.windll.kernel32
        for handle_id in (-11, -12):
            h = k.GetStdHandle(handle_id)
            mode = ctypes.c_uint32()
            if k.GetConsoleMode(h, ctypes.byref(mode)):
                k.SetConsoleMode(h, mode.value | 0x0004)
    except Exception:
        pass


_LEVELS = {"DEBUG": 10, "INFO": 20, "OK": 25, "WARN": 30, "ERROR": 40}
_STYLE = {
    "DEBUG": (_C.GREY, "···"),
    "INFO": (_C.CYAN, " i "),
    "OK": (_C.GREEN, " ✓ "),
    "WARN": (_C.YELLOW, " ! "),
    "ERROR": (_C.RED, " ✗ "),
}


def _out(text: str) -> None:
    """print() that doesn't die on consoles without UTF-8 (cp1252, pipes)."""
    try:
        print(text, flush=True)
    except UnicodeEncodeError:
        enc = getattr(sys.stdout, "encoding", None) or "ascii"
        print(text.encode(enc, "replace").decode(enc, "replace"), flush=True)


class Log:
    level = 20
    no_color = False
    _sink: List[str] = []

    @classmethod
    def _emit(cls, lvl: str, msg: str) -> None:
        cls._sink.append(f"[{lvl:<5}] {msg}")
        if _LEVELS[lvl] < cls.level:
            return
        col, tag = _STYLE[lvl]
        if cls.no_color:
            _out(f"[{tag.strip()}] {msg}")
        else:
            _out(f"{col}{tag}{_C.RESET} {msg}")

    @classmethod
    def debug(cls, m):
        cls._emit("DEBUG", m)

    @classmethod
    def info(cls, m):
        cls._emit("INFO", m)

    @classmethod
    def ok(cls, m):
        cls._emit("OK", m)

    @classmethod
    def warn(cls, m):
        cls._emit("WARN", m)

    @classmethod
    def error(cls, m):
        cls._emit("ERROR", m)

    @classmethod
    def banner(cls, title: str) -> None:
        line = "─" * max(10, 74 - len(title))
        if cls.no_color:
            _out(f"\n== {title} {line}")
        else:
            _out(f"\n{_C.BOLD}{_C.MAGENTA}══ {title} {line}{_C.RESET}")
        cls._sink.append(f"\n=== {title} ===")

    @classmethod
    def dump(cls, path: str) -> None:
        with open(path, "w", encoding="utf-8") as f:
            f.write("\n".join(cls._sink))


# =============================================================================
# 1. IMPORT OpenCascade (cadquery's OCP or pythonocc's OCC.Core)
# =============================================================================


def _load_occ():
    for ns in ("OCP", "OCC.Core"):
        try:
            importlib.import_module(ns + ".gp")
            return ns
        except ImportError:
            continue
    return None


_NS = _load_occ()
if _NS is None:
    print("\n[X] OpenCascade not found.\n    pip install cadquery-ocp   (or conda install -c conda-forge pythonocc-core)\n", file=sys.stderr)
    sys.exit(2)


def _m(name: str):
    return importlib.import_module(f"{_NS}.{name}")


def _st(cls, name: str):
    """
    STATIC method of an OCC class, whatever the binding:
      OCP        : BRep_Tool.Pnt_s
      pythonocc  : BRep_Tool.Pnt   (or brep_tool.Pnt / BRep_Tool_Pnt)
    """
    for nm in (name + "_s", name):
        f = getattr(cls, nm, None)
        if callable(f):
            return f
    mod = sys.modules.get(cls.__module__)
    for nm in (f"{cls.__name__}_{name}",):
        f = getattr(mod, nm, None) if mod else None
        if callable(f):
            return f
    low = getattr(mod, cls.__name__.lower(), None) if mod else None
    f = getattr(low, name, None) if low is not None else None
    if callable(f):
        return f
    raise AttributeError(f"{cls.__name__}.{name} not available in {_NS}")


_gp = _m("gp")
_Geom = _m("Geom")
_Geom2d = _m("Geom2d")
_TopAbs = _m("TopAbs")
_TopoDS = _m("TopoDS")
_TopExp = _m("TopExp")
_TopTools = _m("TopTools")
_BRep = _m("BRep")
_BRepAdaptor = _m("BRepAdaptor")
_BRepBuilderAPI = _m("BRepBuilderAPI")
_BRepGProp = _m("BRepGProp")
_GProp = _m("GProp")
_GeomAbs = _m("GeomAbs")
_ShapeUpgrade = _m("ShapeUpgrade")
_ShapeFix = _m("ShapeFix")
_BRepCheck = _m("BRepCheck")
_STEPControl = _m("STEPControl")
_IFSelect = _m("IFSelect")
_Interface = _m("Interface")
_BRepTools = _m("BRepTools")
_ShapeAnalysis = _m("ShapeAnalysis")
_BRepClass3d = _m("BRepClass3d")
_TColgp = _m("TColgp")
_TColStd = _m("TColStd")
_Geom2dAPI = _m("Geom2dAPI")
_TopLoc = _m("TopLoc")

gp_Pnt, gp_Dir, gp_Vec, gp_Ax2, gp_Ax3, gp_Pnt2d = (_gp.gp_Pnt, _gp.gp_Dir, _gp.gp_Vec, _gp.gp_Ax2, _gp.gp_Ax3, _gp.gp_Pnt2d)
TopAbs_FACE = _TopAbs.TopAbs_FACE
TopAbs_WIRE = _TopAbs.TopAbs_WIRE
TopAbs_EDGE = _TopAbs.TopAbs_EDGE
TopAbs_VERTEX = _TopAbs.TopAbs_VERTEX
TopAbs_SOLID = _TopAbs.TopAbs_SOLID
TopAbs_SHELL = _TopAbs.TopAbs_SHELL
TopAbs_REVERSED = _TopAbs.TopAbs_REVERSED
TopAbs_FORWARD = _TopAbs.TopAbs_FORWARD
TopAbs_OUT = _TopAbs.TopAbs_OUT
TopAbs_IN = _TopAbs.TopAbs_IN
TopExp_Explorer = _TopExp.TopExp_Explorer
TopTools_IndexedMapOfShape = _TopTools.TopTools_IndexedMapOfShape
TopTools_IndexedDataMapOfShapeListOfShape = _TopTools.TopTools_IndexedDataMapOfShapeListOfShape
BRep_Tool = _BRep.BRep_Tool
BRep_Builder = _BRep.BRep_Builder
BRepAdaptor_Surface = _BRepAdaptor.BRepAdaptor_Surface
BRepAdaptor_Curve = _BRepAdaptor.BRepAdaptor_Curve
GeomAbs_Plane = _GeomAbs.GeomAbs_Plane
GeomAbs_Cylinder = _GeomAbs.GeomAbs_Cylinder
GeomAbs_Cone = _GeomAbs.GeomAbs_Cone
GeomAbs_Sphere = _GeomAbs.GeomAbs_Sphere
GeomAbs_Torus = _GeomAbs.GeomAbs_Torus
GeomAbs_Line = _GeomAbs.GeomAbs_Line
GeomAbs_Circle = _GeomAbs.GeomAbs_Circle
GeomAbs_Ellipse = _GeomAbs.GeomAbs_Ellipse
GProp_GProps = _GProp.GProp_GProps
BRepCheck_Analyzer = _BRepCheck.BRepCheck_Analyzer
TopoDS_Compound = _TopoDS.TopoDS_Compound
TopoDS_Shell = _TopoDS.TopoDS_Shell
TopoDS_Face = _TopoDS.TopoDS_Face
TopoDS_Wire = _TopoDS.TopoDS_Wire
TopoDS_Edge = _TopoDS.TopoDS_Edge
TopoDS_Vertex = _TopoDS.TopoDS_Vertex
ShapeUpgrade_UnifySameDomain = _ShapeUpgrade.ShapeUpgrade_UnifySameDomain
ShapeFix_Shape = _ShapeFix.ShapeFix_Shape
ShapeFix_Solid = _ShapeFix.ShapeFix_Solid
STEPControl_Reader = _STEPControl.STEPControl_Reader
STEPControl_Writer = _STEPControl.STEPControl_Writer
STEPControl_AsIs = _STEPControl.STEPControl_AsIs
IFSelect_RetDone = _IFSelect.IFSelect_RetDone
BRepBuilderAPI_MakeEdge = _BRepBuilderAPI.BRepBuilderAPI_MakeEdge
BRepBuilderAPI_MakeWire = _BRepBuilderAPI.BRepBuilderAPI_MakeWire
BRepBuilderAPI_MakeFace = _BRepBuilderAPI.BRepBuilderAPI_MakeFace
BRepBuilderAPI_MakeSolid = _BRepBuilderAPI.BRepBuilderAPI_MakeSolid
BRepTools_ReShape = _BRepTools.BRepTools_ReShape
BRepClass3d_SolidClassifier = _BRepClass3d.BRepClass3d_SolidClassifier
ShapeAnalysis_ShapeTolerance = _ShapeAnalysis.ShapeAnalysis_ShapeTolerance
Geom_Plane = _Geom.Geom_Plane
Geom_CylindricalSurface = _Geom.Geom_CylindricalSurface
Geom_ConicalSurface = _Geom.Geom_ConicalSurface
Geom_SphericalSurface = _Geom.Geom_SphericalSurface
Geom_ToroidalSurface = _Geom.Geom_ToroidalSurface
Geom_Line = _Geom.Geom_Line
Geom_Circle = _Geom.Geom_Circle
Geom_Ellipse = _Geom.Geom_Ellipse
Geom_Hyperbola = _Geom.Geom_Hyperbola
Geom2d_BSplineCurve = _Geom2d.Geom2d_BSplineCurve
Geom_BSplineSurface = _Geom.Geom_BSplineSurface
TopLoc_Location = _TopLoc.TopLoc_Location


# --- cast TopoDS_Shape -> subtype -------------------------------------------
def _cast(kind: str):
    holder = getattr(_TopoDS, "topods", None)
    if holder is not None and hasattr(holder, kind):
        return getattr(holder, kind)
    fn = getattr(_TopoDS, "topods_" + kind, None)
    if fn is not None:
        return fn
    cls = getattr(_TopoDS, "TopoDS", None)
    if cls is not None:
        for nm in (kind + "_s", kind):
            f = getattr(cls, nm, None)
            if callable(f):
                return f
    raise AttributeError(f"cast TopoDS -> {kind} not available in {_NS}")


td_Face, td_Edge, td_Vertex = _cast("Face"), _cast("Edge"), _cast("Vertex")
td_Shell, td_Wire, td_Solid = _cast("Shell"), _cast("Wire"), _cast("Solid")

# --- static methods used everywhere -----------------------------------------
bt_Pnt = _st(BRep_Tool, "Pnt")
bt_Tolerance = _st(BRep_Tool, "Tolerance")
bt_Range = _st(BRep_Tool, "Range")
bt_Curve = _st(BRep_Tool, "Curve")
bt_Surface = _st(BRep_Tool, "Surface")
bt_Degenerated = _st(BRep_Tool, "Degenerated")
te_FirstVertex = _st(_TopExp.TopExp, "FirstVertex")
te_LastVertex = _st(_TopExp.TopExp, "LastVertex")
te_MapShapes = _st(_TopExp.TopExp, "MapShapes")
te_MapAncestors = _st(_TopExp.TopExp, "MapShapesAndAncestors")
gp_Surface = _st(_BRepGProp.BRepGProp, "SurfaceProperties")
gp_Volume = _st(_BRepGProp.BRepGProp, "VolumeProperties")
gp_Linear = _st(_BRepGProp.BRepGProp, "LinearProperties")
brt_OuterWire = _st(_BRepTools.BRepTools, "OuterWire")


def _iface_set(key: str, val: str) -> None:
    IS = _Interface.Interface_Static
    for fn in ("SetCVal_s", "SetCVal"):
        if hasattr(IS, fn):
            try:
                getattr(IS, fn)(key, val)
                return
            except Exception:
                pass


# =============================================================================
# 2. TOPOLOGY UTILITIES
# =============================================================================


def _size(m) -> int:
    for attr in ("Extent", "Size"):
        fn = getattr(m, attr, None)
        if callable(fn):
            return int(fn())
    return int(len(m))


def _iter_list(lst):
    """
    Iterates a TopTools_ListOfShape. ⚠️ With OCP, native iteration (for s in
    lst) is 50 times faster than the ListIterator called from Python: on a
    part with 10,000 edges that's the difference between 3 s and 0.1 s per
    index.
    """
    try:
        return list(lst)
    except TypeError:
        pass
    it_cls = getattr(_TopTools, "TopTools_ListIteratorOfListOfShape", None)
    out = []
    if it_cls is not None:
        it = it_cls(lst)
        while it.More():
            out.append(it.Value())
            it.Next()
    return out


def explore(shape, kind) -> List:
    """Unique sub-shapes (no duplicates, orientation ignored)."""
    m = TopTools_IndexedMapOfShape()
    te_MapShapes(shape, kind, m)
    return [m.FindKey(i) for i in range(1, _size(m) + 1)]


def count_sub(shape, kind) -> int:
    m = TopTools_IndexedMapOfShape()
    te_MapShapes(shape, kind, m)
    return _size(m)


def shape_stats(shape) -> Dict[str, int]:
    return {
        "solids": count_sub(shape, TopAbs_SOLID),
        "shells": count_sub(shape, TopAbs_SHELL),
        "faces": count_sub(shape, TopAbs_FACE),
        "edges": count_sub(shape, TopAbs_EDGE),
        "verts": count_sub(shape, TopAbs_VERTEX),
    }


def edge_face_map(shape):
    m = TopTools_IndexedDataMapOfShapeListOfShape()
    te_MapAncestors(shape, TopAbs_EDGE, TopAbs_FACE, m)
    return m


def count_free_edges(shape) -> int:
    """Edges with a single face: 0 = closed shell."""
    m = edge_face_map(shape)
    return sum(1 for i in range(1, _size(m) + 1) if _size(m.FindFromIndex(i)) == 1)


_SFT = _ShapeFix.ShapeFix_ShapeTolerance


def set_tolerance(sub, tol: float) -> None:
    """
    SETS the tolerance of a vertex/edge (even lower). ⚠️
    BRep_Builder.UpdateVertex/UpdateEdge(tol) can only RAISE it: a rollback
    done with those never brings anything back down.
    """
    try:
        _SFT().SetTolerance(sub, float(max(tol, 1e-9)), sub.ShapeType())
    except Exception:
        pass


def face_edges_map(shape):
    """edge(map index) -> [faces] built by walking the faces.
    ⚠️ ef.FindFromIndex(k) copies a C++ list on every call (0.3 ms): on
    10,000 edges x 30 rebuilds that's 80 s. Built this way it's almost
    free."""
    emap = TopTools_IndexedMapOfShape()
    te_MapShapes(shape, TopAbs_EDGE, emap)
    fmap = TopTools_IndexedMapOfShape()
    te_MapShapes(shape, TopAbs_FACE, fmap)
    e_faces = [[] for _ in range(_size(emap))]
    for i in range(1, _size(fmap) + 1):
        f = fmap.FindKey(i)
        ex = TopExp_Explorer(f, TopAbs_EDGE)
        seen = set()
        while ex.More():
            k = emap.FindIndex(ex.Current())
            if k > 0 and k not in seen:
                seen.add(k)
                e_faces[k - 1].append(i - 1)
            ex.Next()
    return emap, fmap, e_faces


BRepGProp_Face = _BRepGProp.BRepGProp_Face


def vpos(v) -> np.ndarray:
    p = bt_Pnt(td_Vertex(v))
    return np.array([p.X(), p.Y(), p.Z()], dtype=float)


def face_vertices(face) -> np.ndarray:
    return np.array([vpos(v) for v in explore(face, TopAbs_VERTEX)], dtype=float)


def face_area(face) -> float:
    g = GProp_GProps()
    gp_Surface(face, g)
    return abs(float(g.Mass()))


def shape_volume(shape) -> float:
    g = GProp_GProps()
    gp_Volume(shape, g)
    return float(g.Mass())


def edge_length(edge) -> float:
    g = GProp_GProps()
    gp_Linear(edge, g)
    return float(g.Mass())


def face_surface_type(face) -> int:
    return BRepAdaptor_Surface(td_Face(face), True).GetType()


def face_plane_normal(face) -> Optional[np.ndarray]:
    """OUTWARD normal of a planar face (None if not planar)."""
    f = td_Face(face)
    ad = BRepAdaptor_Surface(f, True)
    if ad.GetType() != GeomAbs_Plane:
        return None
    d = ad.Plane().Axis().Direction()
    n = np.array([d.X(), d.Y(), d.Z()], dtype=float)
    if f.Orientation() == TopAbs_REVERSED:
        n = -n
    nn = np.linalg.norm(n)
    return n / nn if nn > 1e-12 else None


def face_plane_point(face) -> np.ndarray:
    ad = BRepAdaptor_Surface(td_Face(face), True)
    p = ad.Plane().Location()
    return np.array([p.X(), p.Y(), p.Z()], dtype=float)


def max_tolerance(shape) -> float:
    try:
        return float(ShapeAnalysis_ShapeTolerance().Tolerance(shape, 1))
    except Exception:
        return float("nan")


_BC_NAMES: Dict[int, str] = {}
for _n in dir(_BRepCheck):
    if _n.startswith("BRepCheck_") and _n != "BRepCheck_Analyzer":
        try:
            _BC_NAMES[int(getattr(_BRepCheck, _n))] = _n[10:]
        except Exception:
            pass


def check_detail(shape, limit: int = 6) -> List[str]:
    """BRepCheck errors in plain text ([] = valid)."""
    try:
        an = BRepCheck_Analyzer(shape)
        if an.IsValid():
            return []
    except Exception as e:
        return [f"BRepCheck: {e}"]
    out: List[str] = []
    for kind, lab in ((TopAbs_FACE, "face"), (TopAbs_WIRE, "wire"), (TopAbs_EDGE, "edge"), (TopAbs_VERTEX, "vertex")):
        for s in explore(shape, kind):
            try:
                res = an.Result(s)
            except Exception:
                continue
            if res is None:
                continue
            try:
                sts = list(res.Status())
            except Exception:
                sts = []
            for st in sts:
                nm = str(st).split("_")[-1]
                if nm and nm != "NoError":
                    tag = f"{lab}:{nm}"
                    if tag not in out:
                        out.append(tag)
                    if len(out) >= limit:
                        return out
    return out


def is_valid(shape) -> bool:
    try:
        return bool(BRepCheck_Analyzer(shape).IsValid())
    except Exception:
        return False


def make_solid_from_faces(shape):
    """Shell (closed, hopefully) + solid, with the orientation checked."""
    b = BRep_Builder()
    shell = TopoDS_Shell()
    b.MakeShell(shell)
    for f in explore(shape, TopAbs_FACE):
        b.Add(shell, f)
    ms = BRepBuilderAPI_MakeSolid()
    ms.Add(shell)
    sol = ms.Solid()
    try:
        cl = BRepClass3d_SolidClassifier(sol)
        cl.PerformInfinitePoint(1e-6)
        if cl.State() == TopAbs_IN:
            Log.warn("Mesh normals point inward: flipping the shell.")
            sol = td_Solid(sol.Reversed())
    except Exception:
        pass
    return sol


def ensure_solid(shape):
    """If the shape isn't a solid (compound/shell), tries to make one."""
    if count_sub(shape, TopAbs_SOLID) >= 1:
        return shape
    return make_solid_from_faces(shape)


# =============================================================================
# 3. I/O
# =============================================================================


def read_stl(path: str):
    """Binary or ASCII STL -> solid of triangular faces with SHARED edges."""
    t0 = time.perf_counter()
    RWStl = _m("RWStl").RWStl
    tri = _st(RWStl, "ReadFile")(path)
    if tri is None or tri.NbTriangles() == 0:
        Log.error(f"Empty or unreadable STL: {path}")
        sys.exit(3)
    mk = _BRepBuilderAPI.BRepBuilderAPI_MakeShapeOnMesh(tri)
    mk.Build()
    sh = mk.Shape()
    sol = make_solid_from_faces(sh)
    st = shape_stats(sol)
    Log.ok(f"STL read in {time.perf_counter() - t0:.2f}s -> {os.path.basename(path)}  ({tri.NbTriangles():,} triangles · {tri.NbNodes():,} vertices)")
    fe = count_free_edges(sol)
    if fe:
        Log.warn(f"Mesh not closed: {fe:,} edges with a single face. They stay open in the output too (no geometry is invented).")
    return sol, st


def read_step(path: str):
    t0 = time.perf_counter()
    r = STEPControl_Reader()
    if r.ReadFile(path) != IFSelect_RetDone:
        Log.error(f"STEP read failed: {path}")
        sys.exit(3)
    r.TransferRoots()
    shape = r.OneShape()
    if shape.IsNull():
        Log.error("The STEP file contains no transferable geometry.")
        sys.exit(3)
    shape = ensure_solid(shape)
    st = shape_stats(shape)
    Log.ok(f"STEP read in {time.perf_counter() - t0:.2f}s -> {os.path.basename(path)}  ({st['faces']:,} faces · {st['solids']} solids)")
    return shape, st


def read_input(path: str):
    if not os.path.isfile(path):
        Log.error(f"File not found: {path}")
        sys.exit(2)
    if os.path.splitext(path)[1].lower() == ".stl":
        return read_stl(path)
    return read_step(path)


def write_step(shape, path: str, schema: str = "AP214IS") -> None:
    _iface_set("write.step.schema", schema)
    _iface_set("write.step.unit", "MM")
    _iface_set("write.precision.mode", "0")
    w = STEPControl_Writer()
    w.Transfer(shape, STEPControl_AsIs)
    if w.Write(path) != IFSelect_RetDone:
        Log.error(f"STEP write failed: {path}")
        return
    Log.ok(f"Saved: {path}  ({os.path.getsize(path) / 1024:,.0f} KB)")


# =============================================================================
# 4. PHASE A — MERGING COPLANAR FACES
# =============================================================================
#
# ⚠️ WHY UnifySameDomain WITH AN ANGULAR TOLERANCE ISN'T ENOUGH.
# The "normals within X degrees" criterion has no right calibration:
#   - at 0.05 degrees (the old default), on the test part the 38 mm-long
#     fillet strips were merged with the corner spheres' triangles
#     (0.01-0.05 degrees of difference) and the merged face carried away a
#     tolerance of 0.03 mm: deformed geometry, silently;
#   - at 0.005 degrees, on the second part (0.001 mm2 triangles, 0.01-0.03
#     degrees of float32 noise from the STL) almost nothing merges anymore.
# The right criterion is a DISTANCE: two facets lie on the same plane if
# ALL the vertices of the group are within lin_tol of the group's mean
# plane. Small, noisy facets merge (the deviation is microscopic), long
# crooked strips don't. The groups are computed here in numpy and handed to
# UnifySameDomain, locking (KeepShape) the edges between different groups:
# that way OCC builds the merged faces, but only merges what we tell it to.


def face_arrays(shape):
    """
    Data for the planar clustering of ANY shape with planar faces (even
    after phases B/C: curved faces remain their own groups).
    Returns V, lists of vertex indices per face, normals (None if curved),
    areas, edge map, faces per edge.
    """
    emap, fmap, e_faces = face_edges_map(shape)
    vmap = TopTools_IndexedMapOfShape()
    te_MapShapes(shape, TopAbs_VERTEX, vmap)
    V = np.array([vpos(vmap.FindKey(j)) for j in range(1, _size(vmap) + 1)], dtype=float).reshape(-1, 3)
    fv, nrm, areas = [], [], []
    for i in range(1, _size(fmap) + 1):
        f = td_Face(fmap.FindKey(i))
        m = TopTools_IndexedMapOfShape()
        te_MapShapes(f, TopAbs_VERTEX, m)
        fv.append([vmap.FindIndex(m.FindKey(k)) - 1 for k in range(1, _size(m) + 1)])
        nrm.append(face_plane_normal(f))
        areas.append(face_area(f))
    return V, fv, nrm, np.array(areas), emap, fmap, e_faces


def face_scale(V, idx, cap: int = 64) -> float:
    """Typical size of a face: median distance between its vertices."""
    if len(idx) < 2:
        return 1e-12
    P = V[idx]
    if len(P) > cap:
        P = P[np.linspace(0, len(P) - 1, cap).astype(int)]
    d = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)
    d = d[d > 0]
    return float(np.median(d)) if d.size else 1e-12


def texture_barrier(V, fv, nrm, areas, e_faces, size_cut: float = 1.5, ang_min_deg: float = 0.2):
    """
    Edges that used to be CAD edges and need to be protected from merging.

    ⚠️ THE MESH REMEMBERS THEM. A tessellator works one CAD face at a time:
    inside a face the texture is uniform, across a CAD edge it changes
    abruptly, because the two sides were tessellated by two independent
    passes (and with two different curvatures: the plane in huge facets,
    the fillet in thin strips). It's visible to the eye on the part and
    it's measurable: on test4 the facet-size jump is 0.00 (median) inside a
    planar face and 1.48 (p90, i.e. almost threefold) across different
    faces.
    Taken alone, the texture jump is sometimes wrong (inside a real planar
    face, Delaunay triangulation mixes slivers and huge triangles), but
    adding one condition is enough: inside a real planar face the facets
    are coplanar to within 0.04 degrees (p99 measured), so a texture jump
    that comes together with a half-degree angle isn't tessellation, it's
    an edge. With "texture > 1.5 AND angle > 0.2 degrees" on test4 the
    false cuts inside a planar face are ZERO and the true cuts 573.
    This matters for TANGENT edges (plane-fillet, where the dihedral is
    nearly zero and no angular criterion sees them): those are exactly the
    ones a loose-tolerance Phase A used to erase, flattening the fillet
    into the plane and taking away from Phase C the surface to recognize.
    """
    nE = len(e_faces)
    bar = np.zeros(nE, dtype=bool)
    if nE == 0:
        return bar
    siz = np.array([face_scale(V, idx) for idx in fv])
    N = np.array([n if n is not None else np.zeros(3) for n in nrm])
    cmin = math.cos(math.radians(ang_min_deg))
    for k, fs in enumerate(e_faces):
        if len(fs) != 2:
            continue
        a, b = fs
        if np.linalg.norm(N[a]) < 0.5 or np.linalg.norm(N[b]) < 0.5:
            continue
        if abs(math.log2(max(siz[a], 1e-12) / max(siz[b], 1e-12))) <= size_cut:
            continue
        if abs(float(N[a] @ N[b])) < cmin:
            bar[k] = True
    return bar


# ⚠️ ang_max 6 degrees: with -a 0.05 on an r=4 fillet, distance alone would
# merge facets up to 18 degrees apart from each other and the fillet would
# turn into a sharp-edged polygon. Truly coplanar faces differ by fractions
# of a degree, so the angular ceiling doesn't remove anything real.
def planar_clusters(V, fv, nrm, areas, e_faces, lin_tol: float, ang_max_deg: float = 6.0, narrow_frac: float = 0.5, barrier=None):
    """
    Planar-group label for every face (union-find).
    Two criteria, both mandatory:
      1. local : sin(angle between the normals) x max extent <= 2 lin_tol,
                 i.e. the REAL deviation the angle produces on the large
                 face. A tiny, noisy facet passes, a 38 mm-long strip
                 tilted 0.03 degrees doesn't;
      2. global: all vertices of the merged group within lin_tol of the
                 group's mean plane.
    With lin_tol = 0.01 you also merge the faces of a surface bulging by 10
    microns: that's the user's choice, and the resulting tolerance is
    recomputed and reported.
    """
    nF = len(fv)
    planar = [n is not None for n in nrm]
    N = np.array([n if n is not None else np.zeros(3) for n in nrm])
    cent = np.array([V[idx].mean(axis=0) if idx else np.zeros(3) for idx in fv])
    Lmax = np.array([float(np.linalg.norm(V[idx].max(axis=0) - V[idx].min(axis=0))) if idx else 0.0 for idx in fv])
    # "narrow" face: area much smaller than the square of its extent.
    # Its normal is barely meaningful, so the local criterion doesn't apply.
    narrow = np.array([a < 0.05 * L * L for a, L in zip(areas, np.maximum(Lmax, 1e-9))])
    narrow_tol = narrow_frac * lin_tol
    pairs = [(fs[0], fs[1]) for k, fs in enumerate(e_faces) if len(fs) == 2 and planar[fs[0]] and planar[fs[1]] and (barrier is None or not barrier[k])]
    parent = np.arange(nF)
    csum_n = N * areas[:, None]
    csum_c = cent * areas[:, None]
    carea = areas.copy()
    cverts = [set(idx) for idx in fv]

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    if pairs:
        pa = np.array(pairs)
        cosd = np.abs(np.einsum("ij,ij->i", N[pa[:, 0]], N[pa[:, 1]]))
        ang = np.degrees(np.arccos(np.clip(cosd, -1.0, 1.0)))
        sin_a = np.sqrt(np.maximum(0.0, 1.0 - cosd**2))
        dev_loc = sin_a * np.maximum(Lmax[pa[:, 0]], Lmax[pa[:, 1]])
        # pairs admitted ONLY because one of the two faces is a sliver:
        # they must be checked against the least-squares plane, not the mean.
        slim = dev_loc > 2.0 * lin_tol
        for k in np.argsort(ang, kind="stable"):
            if ang[k] > ang_max_deg:
                break
            if dev_loc[k] > 2.0 * lin_tol and not (narrow[pa[k, 0]] or narrow[pa[k, 1]]):
                continue
            a_, c_ = find(pa[k, 0]), find(pa[k, 1])
            if a_ == c_:
                continue
            ns = csum_n[a_] + csum_n[c_]
            nn = np.linalg.norm(ns)
            if nn < 1e-12:
                continue
            nrm_ = ns / nn
            cen = (csum_c[a_] + csum_c[c_]) / (carea[a_] + carea[c_])
            vs = np.fromiter(cverts[a_] | cverts[c_], dtype=np.int64)
            Q = V[vs]
            # ⚠️ the AVERAGE normal isn't the best plane: a long, narrow
            # sliver (0.3 x 47 mm) has a poorly determined normal, and one
            # degree is enough to throw away a merge that would be within a
            # micron by least squares. For those slivers we use the
            # least-squares plane of ALL the union's vertices, but with a
            # much tighter tolerance: absorbing a sliver must not cost the
            # rest of the part any precision.
            if float(np.abs((Q - cen) @ nrm_).max()) > lin_tol or slim[k]:
                # ⚠️ TWO TWIN SLIVERS TOO. The least-squares fallback exists to
                # absorb a sliver inside a large face, hence the area ratio.
                # But the tessellation of a fillet produces PAIRS of equal
                # strips (the quadrilateral is warped and the two triangles
                # stay separated by five thousandths of a degree): they are
                # coplanar to within a tenth of a micron, and the 1:1 area
                # ratio kept them apart, doubling the fillet's faces and
                # making the region irregular. If both are slivers the
                # fallback still applies, with its tight tolerance.
                if min(carea[a_], carea[c_]) > 0.2 * max(carea[a_], carea[c_]) and not (narrow[pa[k, 0]] and narrow[pa[k, 1]]):
                    continue
                cen2 = Q.mean(axis=0)
                try:
                    _, _, Vt = np.linalg.svd(Q - cen2, full_matrices=False)
                except np.linalg.LinAlgError:
                    continue
                n2 = Vt[2]
                if float(np.abs((Q - cen2) @ n2).max()) > narrow_tol:
                    continue
                if float(np.abs(n2 @ nrm_)) < math.cos(math.radians(ang_max_deg)):
                    continue
            if len(cverts[a_]) < len(cverts[c_]):
                a_, c_ = c_, a_
            parent[c_] = a_
            csum_n[a_] = ns
            csum_c[a_] = csum_c[a_] + csum_c[c_]
            carea[a_] += carea[c_]
            cverts[a_] |= cverts[c_]
            cverts[c_] = set()
    return np.array([find(i) for i in range(nF)])


def _unify(shape, unify_edges: bool, unify_faces: bool, lin: float, ang_deg: float, keep=None):
    u = ShapeUpgrade_UnifySameDomain(shape, unify_edges, unify_faces, True)
    u.SetLinearTolerance(lin)
    u.SetAngularTolerance(math.radians(ang_deg))
    if hasattr(u, "SetSafeInputMode"):
        u.SetSafeInputMode(True)
    if keep:
        for s in keep:
            u.KeepShape(s)
    # ⚠️ With wide tolerances (e.g. -a 0.05) UnifySameDomain can fail
    # ("Courbes non jointives") merging nearly-collinear edges: in that case
    # the shape is kept as it was, which is still valid.
    try:
        u.Build()
        out = u.Shape()
    except Exception as e:
        Log.warn(f"UnifySameDomain failed ({e}): step skipped, shape unchanged.")
        return shape
    return shape if (out is None or out.IsNull()) else out


def copy_shape(shape):
    """DEEP copy: new faces, new edges, new vertices."""
    cp = _BRepBuilderAPI.BRepBuilderAPI_Copy(shape)
    cp.Perform(shape)
    return cp.Shape()


def refit_face_planes(shape, lin_tol: float) -> int:
    """
    ⚠️ UnifySameDomain gives the merged face the FIRST triangle's plane: on
    a tiny triangle the float32 normal is off by 0.03 degrees and at 4 mm
    away the vertices come out 2e-3 mm off. Here every planar face with
    more than 3 vertices gets the least-squares plane of ITS OWN vertices
    (in place, same edges, same orientation).
    """
    b = BRep_Builder()
    n = 0
    for f in explore(shape, TopAbs_FACE):
        f = td_Face(f)
        ad = BRepAdaptor_Surface(f, True)
        if ad.GetType() != GeomAbs_Plane:
            continue
        V = face_vertices(f)
        if len(V) < 4:
            continue
        c = V.mean(axis=0)
        _, _, Vt = np.linalg.svd(V - c, full_matrices=False)
        nrm = Vt[-1]
        d0 = ad.Plane().Axis().Direction()
        old = np.array([d0.X(), d0.Y(), d0.Z()])
        if float(nrm @ old) < 0:
            nrm = -nrm
        dev_new = float(np.abs((V - c) @ nrm).max())
        p0 = ad.Plane().Location()
        dev_old = float(np.abs((V - np.array([p0.X(), p0.Y(), p0.Z()])) @ old).max())
        if dev_new >= dev_old:
            continue
        loc = TopLoc_Location()
        surf = bt_Surface(f, loc)
        if not loc.IsIdentity():
            continue
        b.UpdateFace(f, Geom_Plane(_mk_pnt(c), _mk_dir(nrm)), loc, 1e-7)
        n += 1
    return n


def removable_vertices(shape, lin_tol: float):
    """
    Vertices that sit on a STRAIGHT edge between just two faces and are
    aligned (within lin_tol) with the chain: the ones UnifySameDomain can
    remove. The criterion is the distance from the chord of the whole
    merged stretch, not the angle between consecutive segments (which on
    short segments is just float32 noise and on long segments deforms).
    Returns (vertices to LOCK, number of removable ones).
    """
    ef, fmap, e_faces = face_edges_map(shape)
    ve = TopTools_IndexedMapOfShape()
    te_MapShapes(shape, TopAbs_VERTEX, ve)
    nE, nV = _size(ef), _size(ve)
    edges = [td_Edge(ef.FindKey(k)) for k in range(1, nE + 1)]
    fkey = [tuple(sorted(fs)) for fs in e_faces]
    is_line = [BRepAdaptor_Curve(e).GetType() == GeomAbs_Line for e in edges]
    vp = np.array([vpos(ve.FindKey(j)) for j in range(1, nV + 1)]).reshape(-1, 3)
    e_verts = [(ve.FindIndex(te_FirstVertex(e)) - 1, ve.FindIndex(te_LastVertex(e)) - 1) for e in edges]
    v_edges = [set() for _ in range(nV)]
    for k, (a_, b_) in enumerate(e_verts):
        for j in (a_, b_):
            if j >= 0:
                v_edges[j].add(k)
    v_edges = [sorted(s) for s in v_edges]
    cand = set()
    for j in range(nV):
        es = v_edges[j]
        if len(es) == 2 and is_line[es[0]] and is_line[es[1]] and fkey[es[0]] == fkey[es[1]] and len(fkey[es[0]]) == 2:
            cand.add(j)

    def other(e, j):
        a_, c_ = e_verts[e]
        return c_ if a_ == j else a_

    removable = set()
    seen = set()
    for j in sorted(cand):
        if j in seen:
            continue
        seen.add(j)
        chain = [j]
        for k_dir, e0 in enumerate(v_edges[j]):
            cur, e = j, e0
            while True:
                nxt = other(e, cur)
                if k_dir == 0:
                    chain.append(nxt)
                else:
                    chain.insert(0, nxt)
                if nxt not in cand or nxt in seen:
                    break
                seen.add(nxt)
                es = [x for x in v_edges[nxt] if x != e]
                if not es:
                    break
                cur, e = nxt, es[0]
        if len(chain) < 3:
            continue
        start, k = 0, 1
        while k < len(chain) - 1:
            P0, P1 = vp[chain[start]], vp[chain[k + 1]]
            d = P1 - P0
            L = float(np.linalg.norm(d))
            ok = L > 1e-12
            if ok:
                d /= L
                q = vp[chain[start + 1 : k + 1]] - P0
                ok = bool(np.linalg.norm(q - np.outer(q @ d, d), axis=1).max() <= lin_tol)
            if ok:
                for m in chain[start + 1 : k + 1]:
                    removable.add(m)
                k += 1
            else:
                start = k
                k += 1
    keep = [ve.FindKey(j + 1) for j in range(nV) if j not in removable]
    return keep, len(removable)


def recompute_tolerances(shape) -> float:
    """
    Tolerances RECOMPUTED from the real geometry (vertex-plane, vertex-curve,
    curve-plane) instead of the ones UnifySameDomain inflated. Set, not just
    raised.
    """
    emap, fmap, e_faces = face_edges_map(shape)
    vmap = TopTools_IndexedMapOfShape()
    te_MapShapes(shape, TopAbs_VERTEX, vmap)
    nV, nE, nF = _size(vmap), _size(emap), _size(fmap)
    planes = []
    for i in range(1, nF + 1):
        ad = BRepAdaptor_Surface(td_Face(fmap.FindKey(i)), True)
        if ad.GetType() == GeomAbs_Plane:
            p0 = ad.Plane().Location()
            d0 = ad.Plane().Axis().Direction()
            planes.append((np.array([p0.X(), p0.Y(), p0.Z()]), np.array([d0.X(), d0.Y(), d0.Z()])))
        else:
            planes.append(None)
    v_faces = [set() for _ in range(nV)]
    v_edges = [set() for _ in range(nE and nV)]
    e_verts = []
    for k in range(1, nE + 1):
        e = td_Edge(emap.FindKey(k))
        a_ = vmap.FindIndex(te_FirstVertex(e)) - 1
        b_ = vmap.FindIndex(te_LastVertex(e)) - 1
        e_verts.append((a_, b_))
        for j in (a_, b_):
            if j >= 0:
                v_edges[j].add(k - 1)
                for i in e_faces[k - 1]:
                    v_faces[j].add(i)
    curves = []
    for k in range(1, nE + 1):
        e = td_Edge(emap.FindKey(k))
        try:
            c = bt_Curve(e, 0.0, 0.0)
            t0, t1 = bt_Range(e)
            # ⚠️ a DEGENERATE edge (a sphere's pole point, a seam that closes
            # up) has no 3D curve: bt_Curve returns None, it doesn't raise.
            curves.append(None if c is None else (c, float(t0), float(t1)))
        except Exception:
            curves.append(None)
    worst = 0.0
    vtol = [1e-7] * nV
    for j in range(nV):
        v = td_Vertex(vmap.FindKey(j + 1))
        P = vpos(v)
        d = 0.0
        for i in v_faces[j]:
            pl = planes[i]
            if pl:
                d = max(d, abs(float((P - pl[0]) @ pl[1])))
        for k in v_edges[j]:
            cv = curves[k]
            if cv is None:
                continue
            c, t0, t1 = cv
            ends = []
            for t in (t0, t1):
                q = c.Value(t)
                ends.append(float(np.linalg.norm(P - np.array([q.X(), q.Y(), q.Z()]))))
            d = max(d, min(ends))
        tol = max(1e-7, 1.2 * d + 1e-9)
        # ⚠️ near CURVED faces the tolerance is needed by the pcurves: never lower it
        if any(planes[i] is None for i in v_faces[j]):
            tol = max(tol, float(bt_Tolerance(v)))
        set_tolerance(v, tol)
        vtol[j] = tol
        worst = max(worst, tol)
    for k in range(nE):
        e = td_Edge(emap.FindKey(k + 1))
        tol = 1e-7
        for j in e_verts[k]:
            if j >= 0:
                tol = max(tol, vtol[j])
        cv = curves[k]
        if cv is not None:
            c, t0, t1 = cv
            ts = np.linspace(t0, t1, 5)
            Q = np.array([[c.Value(t).X(), c.Value(t).Y(), c.Value(t).Z()] for t in ts])
            for i in e_faces[k]:
                pl = planes[i]
                if pl:
                    tol = max(tol, 1.2 * float(np.abs((Q - pl[0]) @ pl[1]).max()) + 1e-9)
        if any(planes[i] is None for i in e_faces[k]):
            tol = max(tol, float(bt_Tolerance(e)))
        set_tolerance(e, tol)
        worst = max(worst, tol)
    return worst


def phase_a(shape, lin_tol: Optional[float] = None, ang_tol_deg: float = 0.005, validate: bool = False, title: str = "PHASE A — merging coplanar faces", use_barrier: bool = True):
    Log.banner(title)
    before = shape_stats(shape)
    Log.info(f"Input    : {before['faces']:,} faces · {before['edges']:,} edges · {before['verts']:,} vertices · {before['solids']} solid")
    free0 = count_free_edges(shape)
    base_ok = is_valid(shape)
    originale = shape
    t0 = time.perf_counter()
    stato: Dict[str, object] = {}

    def _attempt(block_pts):
        """
        ⚠️ WORKING ON A COPY. refit_face_planes, ShapeFix and
        recompute_tolerances modify planes, edges and vertices IN PLACE,
        and they're the same objects as the starting shape's: without a
        copy, after a failed attempt the original is ruined too (measured:
        the input shape, after the round trip, also came out invalid) and
        there's no way to step back.
        block_pts = points inside the planar groups NOT to merge.
        """
        work = copy_shape(originale) if base_ok else originale
        V, fv, nrm, areas, emap, fmap, e_faces = face_arrays(work)
        diag = float(np.linalg.norm(V.max(axis=0) - V.min(axis=0))) if len(V) else 1.0
        lt = lin_tol if lin_tol is not None else max(1e-5, 2e-6 * diag)
        bar = texture_barrier(V, fv, nrm, areas, e_faces) if use_barrier else None
        lab = planar_clusters(V, fv, nrm, areas, e_faces, lt, barrier=bar)
        cen = np.array([V[idx].mean(axis=0) if len(idx) else np.zeros(3) for idx in fv])
        bloccati = set()
        if block_pts is not None and len(block_pts) and len(cen):
            for q in block_pts:
                bloccati.add(int(lab[int(np.argmin(np.linalg.norm(cen - q, axis=1)))]))
        keep = []
        for k, fs in enumerate(e_faces):
            if len(fs) != 2 or lab[fs[0]] != lab[fs[1]] or int(lab[fs[0]]) in bloccati:
                keep.append(emap.FindKey(k + 1))
        if block_pts is None:
            Log.info(
                f"Planar groups (vertices within {lt:.1e} mm of the plane): "
                f"{len(np.unique(lab)):,} · edges locked {len(keep):,}"
                f"{'' if bar is None else f' · of which {int(bar.sum()):,} are CAD edges (texture)'}"
                f"   [{time.perf_counter() - t0:.2f}s]"
            )
        m = _unify(work, False, True, lt, 6.0, keep=keep)
        nfix = refit_face_planes(m, lt)
        keep_v, nrem = removable_vertices(m, lt)
        if block_pts is None:
            Log.info(f"Re-fitted planes {nfix:,} · removable collinear vertices {nrem:,}")
        m = _unify(m, True, False, lt, 30.0, keep=keep_v)
        # ⚠️ UnifySameDomain leaves merged faces with the wire marked
        # "UnorientableShape". ShapeFix_Shape fixes it in place.
        try:
            sf = ShapeFix_Shape(m)
            sf.SetPrecision(1e-7)
            sf.SetMaxTolerance(max(lt, 1e-6))
            sf.Perform()
            fixed = sf.Shape()
            # ⚠️ BUT IT CAN ALSO OPEN THE SHELL: a 0.8-micron edge in the
            # contour removes it from ONE face and leaves it in the other.
            # An "UnorientableShape" wire is a nuisance, an open shell is a
            # defect: if it opens the shell, the uncorrected shape is kept.
            if fixed is not None and not fixed.IsNull():
                if count_free_edges(fixed) <= count_free_edges(m):
                    m = fixed
                else:
                    Log.debug("ShapeFix was opening the shell: discarded")
        except Exception as e:
            Log.debug(f"ShapeFix after Unify skipped: {e}")
        stato["lt"] = lt
        return m

    merged = _attempt(None)
    # ⚠️ PHASE A MUST NOT MAKE THE SOLID WORSE. UnifySameDomain and ShapeFix
    # are robust but not infallible: on a free-form face with a contour of
    # hundreds of segments, the merge can occasionally leave a wire that
    # BRepCheck rejects. Merging faces is cosmetic, an invalid solid is a
    # defect. But throwing away the WHOLE merge for a single face is
    # disproportionate: we retry, blocking only the planar groups the bad
    # faces came from, and only give up if that isn't enough.
    if base_ok and not is_valid(merged):
        guasti = []
        for _ in range(1):
            nuovi = []
            for f in explore(merged, TopAbs_FACE):
                f = td_Face(f)
                if check_detail(f):
                    W = face_vertices(f)
                    if len(W):
                        nuovi.append(W.mean(axis=0))
            if not nuovi:
                break
            guasti.extend(nuovi)
            Log.debug(f"Phase A: {len(nuovi)} invalid merged faces, retrying leaving their groups alone ({len(guasti)} total)")
            merged = _attempt(guasti)
            if is_valid(merged):
                break
        if not is_valid(merged):
            Log.warn("Phase A: the merge made the solid invalid, kept the starting shape")
            after0 = shape_stats(originale)
            Log.ok(f"Output   : {after0['faces']:,} faces (no merge)   [{time.perf_counter() - t0:.2f}s]")
            return originale, before, after0
    worst = recompute_tolerances(merged)
    merged = ensure_solid(merged)
    dt = time.perf_counter() - t0

    after = shape_stats(merged)
    free1 = count_free_edges(merged)
    red_f = 100.0 * (1 - after["faces"] / max(1, before["faces"]))
    Log.ok(f"Output   : {after['faces']:,} faces · {after['edges']:,} edges · {after['verts']:,} vertices   [{dt:.2f}s]   face reduction -{red_f:.1f}%")
    Log.info(f"Free edges: {free1:,} (input {free0:,}) · recomputed max tolerance {worst:.1e} mm")
    if free1 > free0:
        Log.warn("Phase A opened up the shell somewhere (mesh non-manifold there).")
    if validate:
        ok = is_valid(merged)
        (Log.ok if ok else Log.warn)(f"BRepCheck after Phase A: {'OK' if ok else 'NOT valid'}")
    return merged, before, after


# =============================================================================
# 5. PRIMITIVES: algebraic fit + nonlinear refinement (from the previous engine)
# =============================================================================


def taubin_circle(x: np.ndarray, y: np.ndarray, iters: int = 40, eps: float = 1e-12) -> Tuple[float, float, float]:
    """Least-squares circle fit (Taubin). Returns (cx, cy, r)."""
    n = x.size
    if n < 3:
        raise ValueError("need at least 3 points")
    mx, my = x.mean(), y.mean()
    u, v = x - mx, y - my
    z = u * u + v * v

    Mz = z.mean()
    Mxy = (u * v).mean()
    Mxx = (u * u).mean()
    Myy = (v * v).mean()
    Mxz = (u * z).mean()
    Myz = (v * z).mean()
    Mzz = (z * z).mean()

    Cov_xy = Mxx * Myy - Mxy * Mxy
    Var_z = Mzz - Mz * Mz

    A3 = 4.0 * Mz
    A2 = -3.0 * Mz * Mz - Mzz
    A1 = Var_z * Mz + 4.0 * Cov_xy * Mz - Mxz * Mxz - Myz * Myz
    A0 = Mxz * (Mxz * Myy - Myz * Mxy) + Myz * (Myz * Mxx - Mxz * Mxy) - Var_z * Cov_xy
    A22, A33 = 2.0 * A2, 3.0 * A3

    xn, yn = 0.0, 1e30
    for _ in range(iters):
        yo = yn
        yn = A0 + xn * (A1 + xn * (A2 + xn * A3))
        if abs(yn) > abs(yo):
            xn = 0.0
            break
        dy = A1 + xn * (A22 + xn * A33)
        if abs(dy) < eps:
            break
        xo, xn = xn, xn - yn / dy
        if xn == 0.0 or abs((xn - xo) / xn) < eps:
            break
        if xn < 0.0:
            xn = 0.0
            break

    det = xn * xn - xn * Mz + Cov_xy
    if abs(det) < eps:  # degenerate -> Kasa
        A = np.column_stack([u, v, np.ones(n)])
        sol, *_ = np.linalg.lstsq(A, z, rcond=None)
        cx, cy = sol[0] / 2.0, sol[1] / 2.0
        r = math.sqrt(max(0.0, sol[2] + cx * cx + cy * cy))
        return cx + mx, cy + my, r

    cx = (Mxz * (Myy - xn) - Myz * Mxy) / det / 2.0
    cy = (Myz * (Mxx - xn) - Mxz * Mxy) / det / 2.0
    r = math.sqrt(max(0.0, cx * cx + cy * cy + Mz))
    return cx + mx, cy + my, r


def ortho_frame(axis: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    tmp = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(axis, tmp)
    u /= np.linalg.norm(u)
    v = np.cross(axis, u)
    return u, v


# --- 5.2  geometric primitives ------------------------------------------------

PLANE, AXIAL, SPHERE, TORUS = "plane", "axial", "sphere", "torus"
FREE = "free"  # free-form surface (B-spline) — see 5.4bis


@dataclass
class Prim:
    """
    kind == PLANE  : center (point), axis (normal)
    kind == AXIAL  : center (point on the axis at t=0), axis, r0, slope
                     radius(t) = r0 + slope*t   ->  slope==0 cylinder, else cone
    kind == SPHERE : center, r0
    """

    kind: str
    center: np.ndarray
    axis: Optional[np.ndarray] = None
    r0: float = 0.0  # TORUS: MAJOR radius (tube's axis)
    slope: float = 0.0
    r1: float = 0.0  # TORUS: MINOR radius (of the tube)
    rms: float = 1e30
    free: Optional["FreeForm"] = None  # FREE: the free-form surface

    # --- geometry -------------------------------------------------------------
    def _tr(self, P: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Cylindrical coordinates: (axial t, radial rho, radial unit vector)."""
        d = P - self.center
        t = d @ self.axis
        rad = d - np.outer(t, self.axis)
        rho = np.linalg.norm(rad, axis=1)
        uh = rad / np.maximum(rho, 1e-12)[:, None]
        return t, rho, uh

    def dist(self, P: np.ndarray) -> np.ndarray:
        """(Signed) distance of the points from the surface."""
        if self.kind == FREE:
            return self.free.dist(P)
        if self.kind == PLANE:
            return (P - self.center) @ self.axis
        if self.kind == SPHERE:
            return np.linalg.norm(P - self.center, axis=1) - self.r0
        if self.kind == TORUS:
            t, rho, _ = self._tr(P)
            return np.hypot(rho - self.r0, t) - self.r1
        t, rho, _ = self._tr(P)
        # for the cone, the true distance is the radial deviation * cos(semi-angle)
        return (rho - (self.r0 + self.slope * t)) / math.hypot(1.0, self.slope)

    def normal_at(self, P: np.ndarray) -> np.ndarray:
        """(Unoriented) normal of the primitive at the given points."""
        if self.kind == FREE:
            return self.free.normal_at(P)
        if self.kind == PLANE:
            return np.tile(self.axis, (len(P), 1))
        if self.kind == SPHERE:
            d = P - self.center
            return d / np.maximum(np.linalg.norm(d, axis=1), 1e-12)[:, None]
        if self.kind == TORUS:
            t, rho, uh = self._tr(P)
            n = (rho - self.r0)[:, None] * uh + t[:, None] * self.axis
            return n / np.maximum(np.linalg.norm(n, axis=1), 1e-12)[:, None]
        _, _, uh = self._tr(P)
        n = uh - self.slope * self.axis
        return n / np.maximum(np.linalg.norm(n, axis=1), 1e-12)[:, None]

    def label(self) -> str:
        if self.kind == FREE:
            return "FREE"
        if self.kind == PLANE:
            return "PLANE"
        if self.kind == SPHERE:
            return "SPHERE"
        if self.kind == TORUS:
            return "TORUS"
        return "CYL" if abs(self.slope) < 1e-3 else "CONE"


# --- 5.3  primitive fitting ---------------------------------------------------


def fit_plane(P: np.ndarray) -> Optional[Prim]:
    if len(P) < 3:
        return None
    c = P.mean(axis=0)
    _, _, Vt = np.linalg.svd(P - c, full_matrices=False)
    n = Vt[-1]
    p = Prim(PLANE, c, n / np.linalg.norm(n))
    p.rms = float(np.sqrt(np.mean(p.dist(P) ** 2)))
    return p


def fit_sphere(P: np.ndarray) -> Optional[Prim]:
    if len(P) < 4:
        return None
    A = np.column_stack([2.0 * P, np.ones(len(P))])
    b = (P**2).sum(axis=1)
    try:
        sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    except np.linalg.LinAlgError:
        return None
    c = sol[:3]
    rr = sol[3] + c @ c
    if not np.isfinite(rr) or rr <= 1e-12:
        return None
    p = Prim(SPHERE, c, None, math.sqrt(rr))
    p.rms = float(np.sqrt(np.mean(p.dist(P) ** 2)))
    return p


def _cone_algebraic(px, py, t, max_slope: float):
    """
    Cone with a known axis, in closed form.

    ⚠️ A SINGLE CIRCLE WON'T DO. Projecting a cone onto the plane
    orthogonal to the axis, the points sit on circles of DIFFERENT radius
    (one per height): fitting just one shifts the center, the radii come
    out wrong and the taper estimated by regression is off by an order of
    magnitude (-0.07 instead of -1). Here we use the cone's equation
        (x-cx)^2 + (y-cy)^2 = (r0 + s t)^2
    which, expanded, is LINEAR in the unknowns (cx, cy, A, B, C) with
    A = cx^2+cy^2-r0^2, B = 2 r0 s, C = s^2. Taper and radius, sign
    included, are recovered from C and B.
    """
    if len(t) < 6 or float(t.max() - t.min()) < 1e-9:
        return None
    M = np.column_stack([2.0 * px, 2.0 * py, -np.ones_like(px), t, t * t])
    rhs = px * px + py * py
    try:
        sol, *_ = np.linalg.lstsq(M, rhs, rcond=None)
    except np.linalg.LinAlgError:
        return None
    cx, cy, A, B, C = (float(x) for x in sol)
    if not all(np.isfinite(v) for v in (cx, cy, A, B, C)) or C < 0.0:
        return None
    s = math.sqrt(C)
    if s < 1e-9:
        return None
    if B < 0:
        s = -s
    r0 = 0.5 * B / s
    if not np.isfinite(r0) or r0 <= 1e-9 or abs(s) > max_slope:
        return None
    return cx, cy, r0, s


def _cone_rings(px, py, t, max_slope: float):
    """
    Cone with a known axis, derived from RINGS.

    ⚠️ A band tessellated with a SINGLE ROW of facets has vertices at only
    two heights: the cone's algebraic form becomes singular (with two
    values of t, the t^2 column is a linear combination of t and the
    constant) and the taper comes out random. But those two rings are
    exactly what's needed: one circle per ring gives center and radius, and
    two radii at two heights give the exact taper.
    """
    n = len(t)
    if n < 6:
        return None
    span = float(t.max() - t.min())
    if span < 1e-9:
        return None
    order = np.argsort(t)
    groups, cur = [], [int(order[0])]
    for i in order[1:]:
        i = int(i)
        if t[i] - t[cur[-1]] > 0.05 * span:
            groups.append(cur)
            cur = [i]
        else:
            cur.append(i)
    groups.append(cur)
    groups = [g for g in groups if len(g) >= 3]
    # ⚠️ needed for bands with one or two rows of facets: with many rings the
    # normal fit already has all the data it needs, and here we'd pay for a
    # circle fit per ring on every call.
    if not 2 <= len(groups) <= 4:
        return None
    cs = []
    for g in groups:
        try:
            cx, cy, _ = taubin_circle(px[g], py[g])
        except Exception:
            continue
        if np.isfinite(cx) and np.isfinite(cy):
            cs.append((cx, cy, len(g)))
    if len(cs) < 2:
        return None
    wsum = sum(k for _, _, k in cs)
    cx = sum(a * k for a, _, k in cs) / wsum
    cy = sum(b * k for _, b, k in cs) / wsum
    tt = np.array([float(t[g].mean()) for g in groups])
    rr = np.array([float(np.hypot(px[g] - cx, py[g] - cy).mean()) for g in groups])
    if len(tt) == 2:
        slope = float((rr[1] - rr[0]) / (tt[1] - tt[0]))
        r0 = float(rr[0] - slope * tt[0])
    else:
        slope, r0 = (float(x) for x in np.polyfit(tt, rr, 1))
    if not (np.isfinite(r0) and np.isfinite(slope)) or abs(slope) > max_slope:
        return None
    return cx, cy, r0, slope


def _axial_on_axis(P: np.ndarray, axis: np.ndarray, force_cyl: bool, max_slope: float) -> Optional[Prim]:
    """Least-squares radius and taper around a given axis."""
    O = P.mean(axis=0)
    u, v = ortho_frame(axis)
    Q = P - O
    px, py, t = Q @ u, Q @ v, Q @ axis
    best = None
    variants = []
    try:
        cx, cy, _ = taubin_circle(px, py)
        if np.isfinite(cx) and np.isfinite(cy):
            rho = np.hypot(px - cx, py - cy)
            if (not force_cyl) and t.max() - t.min() > 1e-9 and len(t) > 3:
                slope, r0 = np.polyfit(t, rho, 1)
            else:
                slope, r0 = 0.0, float(rho.mean())
            variants.append((cx, cy, float(r0), float(slope)))
    except Exception:
        pass
    # ⚠️ shortcut: if the first attempt already hugs the points there's
    # nothing to gain by trying the other two forms, which cost a linear
    # system and a couple of circle fits per call.
    good = float(np.ptp(np.hypot(px, py)) + (t.max() - t.min())) * 1e-5 + 1e-12
    if variants and not force_cyl:
        cx0, cy0, r00, s0 = variants[0]
        q0 = Prim(AXIAL, P.mean(axis=0) + cx0 * u + cy0 * v, axis, float(r00), float(s0))
        if float(np.abs(q0.dist(P)).max()) > good:
            for fn in (_cone_algebraic, _cone_rings):
                alt = fn(px, py, t, max_slope)
                if alt is not None:
                    variants.append(alt)
    elif not force_cyl:
        for fn in (_cone_algebraic, _cone_rings):
            alt = fn(px, py, t, max_slope)
            if alt is not None:
                variants.append(alt)
    tm = 0.5 * float(t.min() + t.max())
    for cx, cy, r0, slope in variants:
        if not np.isfinite(r0) or abs(slope) > max_slope:
            continue
        # origin at the data's center: this way the reference radius is the
        # band's real one, not the (maybe negative) one at the cone's apex
        cen = O + cx * u + cy * v + tm * axis
        rm = float(r0 + slope * tm)
        if not np.isfinite(rm) or rm <= 1e-6:
            continue
        q = Prim(AXIAL, cen, axis, rm, float(slope))
        q.rms = float(np.sqrt(np.mean(q.dist(P) ** 2)))
        if best is None or q.rms < best.rms:
            best = q
    return best


def fit_axial(P: np.ndarray, N: np.ndarray, W: np.ndarray, max_slope: float = 3.0, ratio_1d: float = 0.05) -> Optional[Prim]:
    """
    Cylinder or cone.

    The normals satisfy  n·a = s  constant  (s = 0 cylinder, s = sin(alpha)
    cone): i.e. the points {n_i} in normal-space lie on a PLANE. The fit is
    therefore the smallest eigenvector of the CENTERED normals' covariance.
    Centering is mandatory: without it, partial arcs and cones get the axis
    wrong.

    ⚠️ DEGENERATE CASE. A SMALL cylinder neighborhood has nearly COLLINEAR
    normals (they sit on a very short arc). Infinitely many planes contain
    a line, so the smallest eigenvector is arbitrary and the axis comes out
    random. In that case the degeneracy is resolved by picking the
    cylinder interpretation (s = 0): the axis must be orthogonal to both
    the mean normal and the spread direction, so  a = mean_n x w.

    ⚠️ BUT NOT ALWAYS. A CURVED CHAMFER (a 45-degree cone around a rounded
    edge) also has its normals on a short arc, yet its axis is very well
    determined: the normals' component along the axis is CONSTANT, and the
    smallest eigenvector finds it. If the cylinder is forced, the axis
    comes out orthogonal to the true one and the fit is off by hundredths;
    then a sphere through the two boundary circles "explains" the band
    better than the fake cylinder, and a cap turns up in place of the
    chamfer. So in the degenerate case BOTH interpretations are tried and
    the one that hugs the points more closely is kept.
    """
    if len(P) < 6 or len(N) < 3:
        return None
    w = W / max(W.sum(), 1e-300)
    nb = (N * w[:, None]).sum(axis=0)
    D = N - nb
    M = (D * w[:, None]).T @ D
    evals, evecs = np.linalg.eigh(M)

    if evals[2] <= 1e-16:
        return None  # all normals identical

    cands = []
    if evals[1] <= ratio_1d * evals[2]:
        # --- 1D spread: degenerate ---
        spread = evecs[:, 2]
        a = np.cross(nb, spread)
        na = np.linalg.norm(a)
        if na > 1e-9:
            cands.append((a / na, True))  # cylinder interpretation
        n0 = np.linalg.norm(evecs[:, 0])
        if n0 > 1e-9:
            cands.append((evecs[:, 0] / n0, False))  # cone interpretation
    else:
        # --- 2D spread: the plane fit in normal-space is valid ---
        if evals[0] / evals[1] > 0.15:
            return None  # 3D spread -> sphere, not axial
        n0 = np.linalg.norm(evecs[:, 0])
        if n0 < 1e-9:
            return None
        cands.append((evecs[:, 0] / n0, False))

    best = None
    for axis, force_cyl in cands:
        q = _axial_on_axis(P, axis, force_cyl, max_slope)
        if q is not None and (best is None or q.rms < best.rms):
            best = q
    return best


def fit_axial_rev(P: np.ndarray, Nrep: np.ndarray, Wp: np.ndarray) -> Optional[Prim]:
    """
    Cone/cylinder with the axis taken from revolution_axis instead of the
    normals' covariance.

    ⚠️ CURVED CHAMFER, A SINGLE ROW OF FACETS. A 45-degree chamfer around a
    rounded edge is a CONE, but tessellated it has vertices only on the two
    boundary circles and normals on a short arc: the normals' covariance is
    degenerate (1D spread), fit_axial falls back to the cylinder and gets
    the axis badly wrong. Then the same surface also gets "explained" very
    well by a SPHERE through the two circles, and the part fills up with
    spherical caps in place of the chamfers. The positions do carry the
    information, though: n . (a x (p - c)) = 0 is linear in (a, a x c) and
    gives the axis with no degenerate cases. From there, two least-squares
    fits are enough for radius and taper. The final verdict is still given
    by the normals' deviation, which on the cone is a third of the
    sphere's.
    """
    ra = revolution_axis(P, Nrep, Wp)
    if ra is None:
        return None
    a, c = ra
    na = float(np.linalg.norm(a))
    if na < 1e-9:
        return None
    a = a / na
    t = (P - c) @ a
    rho = np.linalg.norm((P - c) - np.outer(t, a), axis=1)
    if t.max() - t.min() < 1e-9 or len(t) < 4:
        return None
    try:
        slope, r0 = np.polyfit(t, rho, 1)
    except Exception:
        return None
    if not (np.isfinite(r0) and np.isfinite(slope)) or abs(slope) > 3.0:
        return None
    # ⚠️ the axis's reference point can fall beyond the cone's apex, and
    # there the "radius at t=0" is NEGATIVE: it's not a wrong fit, it's just
    # an inconvenient origin. The origin is moved to the middle of the data.
    tm = 0.5 * float(t.min() + t.max())
    c = c + tm * a
    r0 = float(r0 + slope * tm)
    if r0 <= 1e-6:
        return None
    q = Prim(AXIAL, c, a, float(r0), float(slope))
    q.rms = float(np.sqrt(np.mean(q.dist(P) ** 2)))
    return q


def revolution_axis(P: np.ndarray, N: np.ndarray, W: np.ndarray):
    """
    Axis of a surface of revolution, in closed form.

    For ANY surface of revolution (cylinder, cone, sphere, torus) the
    normal at a point has no tangential component:  n . (a x (p - c)) = 0.
    Written as  a.(p x n) = (a x c).n , it's LINEAR and homogeneous in
    (a, a x c): the axis comes out of the SVD, with no seeding and no
    degenerate cases. It's much more robust than the normals' covariance
    on thin strips.
    """
    if len(P) < 6:
        return None
    A = np.column_stack([np.cross(P, N), -N]) * np.sqrt(np.maximum(W, 1e-12))[:, None]
    try:
        _, _, Vt = np.linalg.svd(A, full_matrices=False)
    except np.linalg.LinAlgError:
        return None
    m, g = Vt[-1][:3], Vt[-1][3:]
    nm = float(np.linalg.norm(m))
    if nm < 1e-9:
        return None
    m = m / nm
    return m, -np.cross(m, g / nm)


def fit_torus(P: np.ndarray, N: np.ndarray, W: np.ndarray) -> Optional[Prim]:
    """
    Torus: the FILLET ALONG A CURVED EDGE.

    ⚠️ Without this primitive a real mechanical part won't convert. Every
    chamfer or fillet around a hole, a hub or a rounded edge is a TORUS,
    not a cylinder: on its triangles the best cylinder is off by a tenth of
    a millimeter and the region gets expelled every pass. What stays
    tessellated are exactly the round edges, which are the ones you see.
    In the meridian plane (rho, z) the torus's profile is a CIRCLE: axis
    from the SVD, then a Taubin circle on the profile.
    """
    q = revolution_axis(P, N, W)
    if q is None:
        return None
    a, c = q
    d = P - c
    z = d @ a
    rho = np.linalg.norm(d - np.outer(z, a), axis=1)
    try:
        R, z0, r = taubin_circle(rho, z)
    except Exception:
        return None
    if not (np.isfinite(R) and np.isfinite(r) and np.isfinite(z0)):
        return None
    if r <= 1e-9 or R <= 1e-9:
        return None
    if R < 0.15 * r or R > 60.0 * r:
        return None  # degenerates into a sphere or a cylinder
    p = Prim(TORUS, c + z0 * a, a, float(R), 0.0, float(r))
    p.rms = float(np.sqrt(np.mean(p.dist(P) ** 2)))
    return p


# --- 5.4bis  the FREE-FORM surface --------------------------------------------
# ⚠️ WHY IT'S NEEDED. A chamfer or fillet running along a CURVED edge is
# neither a cone nor a torus: it's a blending surface that the CAD writes
# as a B-spline. The primitive fit splits it into dozens of 10-degree
# osculating little cylinders, each with a jagged outline: those are
# exactly the "faces with too many sides" you see in the finished part.
# Here the patch is taken WHOLE and approximated with a least-squares
# tensor B-spline. It's not a compromise on precision: unlike a quadric, a
# B-spline has as many degrees of freedom as it needs, and on these
# patches it reaches a hundredth of the tolerance the little cylinders were
# barely passing with. One face instead of forty, and closer to the mesh.
#
# The surface is a HEIGHT FIELD over a base plane:
#     S(u, v) = C + u*X + v*Y + h(u, v)*Z
# so (u, v) are Cartesian coordinates on the base plane and the
# parametrization matches the one SurfParam uses for quadrics. It only
# holds for patches that don't fold back onto the base plane: the check is
# the normals' tilt (past ~70 degrees the patch is rejected and stays
# tessellated).

_BS_DEG = 3


def _bs_knots(t0: float, t1: float, n: int, deg: int = _BS_DEG) -> np.ndarray:
    """Uniform clamped knots for n poles of degree deg (full vector)."""
    inner = np.linspace(t0, t1, n - deg + 1)
    return np.concatenate([np.full(deg, t0), inner, np.full(deg, t1)])


def _bs_basis(t: np.ndarray, kv: np.ndarray, n: int, deg: int = _BS_DEG) -> np.ndarray:
    """(len(t), n) matrix of B-spline basis functions (Cox-de Boor recurrence)."""
    t = np.clip(np.asarray(t, float), kv[0], kv[-1])
    N = np.zeros((len(t), n + deg))
    idx = np.clip(np.searchsorted(kv, t, side="right") - 1, deg, n - 1)
    N[np.arange(len(t)), idx] = 1.0
    for d in range(1, deg + 1):
        Nn = np.zeros_like(N)
        for i in range(n + deg - d):
            d1 = kv[i + d] - kv[i]
            d2 = kv[i + d + 1] - kv[i + 1]
            a = (t - kv[i]) / d1 * N[:, i] if d1 > 0 else 0.0
            b = (kv[i + d + 1] - t) / d2 * N[:, i + 1] if d2 > 0 else 0.0
            Nn[:, i] = a + b
        N = Nn
    return N[:, :n]


def _bs_greville(kv: np.ndarray, n: int, deg: int = _BS_DEG) -> np.ndarray:
    """Greville abscissae: with the poles placed there the B-spline reproduces t."""
    return np.array([kv[i + 1 : i + deg + 1].mean() for i in range(n)])


def _d2_matrix(n: int) -> np.ndarray:
    if n < 3:
        return np.zeros((0, n))
    D = np.zeros((n - 2, n))
    for i in range(n - 2):
        D[i, i], D[i, i + 1], D[i, i + 2] = 1.0, -2.0, 1.0
    return D


@dataclass
class FreeForm:
    """B-spline height field over a plane: S(u,v) = C + uX + vY + h(u,v)Z."""

    C: np.ndarray
    X: np.ndarray
    Y: np.ndarray
    Z: np.ndarray
    poles: np.ndarray  # (nu, nv): pole heights
    ku: np.ndarray
    kv: np.ndarray

    def height(self, u, v) -> np.ndarray:
        u = np.atleast_1d(np.asarray(u, float))
        v = np.atleast_1d(np.asarray(v, float))
        Bu = _bs_basis(u, self.ku, self.poles.shape[0])
        Bv = _bs_basis(v, self.kv, self.poles.shape[1])
        return np.einsum("ni,ij,nj->n", Bu, self.poles, Bv)

    def uv(self, P: np.ndarray):
        d = np.atleast_2d(P) - self.C
        return d @ self.X, d @ self.Y

    def point(self, u, v) -> np.ndarray:
        u = np.atleast_1d(np.asarray(u, float))
        v = np.atleast_1d(np.asarray(v, float))
        return self.C + u[:, None] * self.X + v[:, None] * self.Y + self.height(u, v)[:, None] * self.Z

    def _grad(self, u, v):
        h = 1e-4 * max(float(self.ku[-1] - self.ku[0]), 1e-9)
        return ((self.height(u + h, v) - self.height(u - h, v)) / (2.0 * h), (self.height(u, v + h) - self.height(u, v - h)) / (2.0 * h))

    def dist(self, P: np.ndarray) -> np.ndarray:
        """
        (Signed) PERPENDICULAR distance from the surface.
        ⚠️ The deviation along Z would only be fine for a flat patch: on a
        strip tilted 70 degrees it overestimates threefold, and a
        comparison between two surfaces measured that way means nothing.
        It's corrected with the cosine of the local slope, which is exact
        to first order - and we need microns on slopes that change plane.
        """
        d = np.atleast_2d(P) - self.C
        u, v = d @ self.X, d @ self.Y
        dz = (d @ self.Z) - self.height(u, v)
        gu, gv = self._grad(u, v)
        return dz / np.sqrt(1.0 + gu * gu + gv * gv)

    def normal_at(self, P: np.ndarray) -> np.ndarray:
        u, v = self.uv(P)
        gu, gv = self._grad(u, v)
        n = -gu[:, None] * self.X - gv[:, None] * self.Y + np.ones((len(np.atleast_1d(u)), 1)) * self.Z
        return n / np.maximum(np.linalg.norm(n, axis=1), 1e-12)[:, None]

    def geom(self):
        """Geom_BSplineSurface with the SAME parametrization (poles at the Greville points)."""
        nu, nv = self.poles.shape
        gu = _bs_greville(self.ku, nu)
        gv = _bs_greville(self.kv, nv)
        arr = _TColgp.TColgp_Array2OfPnt(1, nu, 1, nv)
        for i in range(nu):
            for j in range(nv):
                q = self.C + gu[i] * self.X + gv[j] * self.Y + self.poles[i, j] * self.Z
                arr.SetValue(i + 1, j + 1, gp_Pnt(float(q[0]), float(q[1]), float(q[2])))

        def _kn(kk):
            vals = np.unique(np.round(kk, 12))
            ka = _TColStd.TColStd_Array1OfReal(1, len(vals))
            ma = _TColStd.TColStd_Array1OfInteger(1, len(vals))
            for i, t in enumerate(vals):
                ka.SetValue(i + 1, float(t))
                ma.SetValue(i + 1, int((np.abs(kk - t) < 1e-12).sum()))
            return ka, ma

        kau, mau = _kn(self.ku)
        kav, mav = _kn(self.kv)
        return Geom_BSplineSurface(arr, kau, kav, mau, mav, _BS_DEG, _BS_DEG, False, False)


def _free_tame(ff: "FreeForm", x, y, z, ngrid: int = 25) -> bool:
    """The surface doesn't rear up BETWEEN the data points: a grid is
    sampled and only the nodes with nearby data are checked (outside the
    patch the spline is free to go wherever, nobody uses it there)."""
    ex, ey = float(np.ptp(x)), float(np.ptp(y))
    hz = float(np.ptp(z))
    gx = np.linspace(x.min(), x.max(), ngrid)
    gy = np.linspace(y.min(), y.max(), ngrid)
    GX, GY = np.meshgrid(gx, gy, indexing="ij")
    d2 = ((GX.ravel()[:, None] - x[None, :]) ** 2 + (GY.ravel()[:, None] - y[None, :]) ** 2).min(axis=1)
    rad = (1.5 * max(ex / ngrid, ey / ngrid)) ** 2
    m = d2 <= rad
    if not m.any():
        return True
    H = ff.height(GX.ravel()[m], GY.ravel()[m])
    lim = 0.5 * hz + 1e-6
    return bool(H.max() <= z.max() + lim and H.min() >= z.min() - lim)


def fit_free(P: np.ndarray, Nrep: np.ndarray, tol: float, check: Optional[np.ndarray] = None, check_tol: Optional[float] = None, max_tilt_deg: float = 70.0) -> Optional[Prim]:
    """
    Least-squares B-spline over a point cloud with normals.

    tol       : maximum deviation allowed on the mesh's VERTICES.
    check     : points INSIDE the facets (centroids and side midpoints).
    check_tol : how far those are allowed to be.

    ⚠️ THE VERTICES TELL YOU NOTHING. A B-spline with enough poles passes
    through every vertex and does whatever it wants between one vertex and
    the next: on a test4 patch, 19x17 poles give 9e-05 on the vertices and
    7.6e-02 (seventy-five times the tolerance) in between. The two knobs
    are the grid's DENSITY and the REGULARIZATION, and they don't pull the
    same way: more poles bring the vertices closer and make it wobble, more
    regularization smooths and pushes it away. Here both are tried and the
    combination that passes the test even IN THE EMPTY SPACES is kept, the
    sparsest one possible. If none passes, None is returned: better
    tessellated than invented.
    """
    P = np.unique(np.round(np.asarray(P, float), 9), axis=0)
    if len(P) < 30:
        return None
    pl = fit_plane(P)
    if pl is None:
        return None
    Z = pl.axis / np.linalg.norm(pl.axis)
    Nrep = np.atleast_2d(np.asarray(Nrep, float))
    if float(np.sum(Nrep @ Z)) < 0.0:
        Z = -Z
    ct = Nrep @ Z
    if float(np.quantile(ct, 0.02)) < math.cos(math.radians(max_tilt_deg)):
        return None  # the patch folds back on itself
    X, Y = ortho_frame(Z)
    C = P.mean(axis=0)
    d = P - C
    x, y, z = d @ X, d @ Y, d @ Z
    ex, ey = float(np.ptp(x)), float(np.ptp(y))
    if min(ex, ey) < 1e-9:
        return None
    x0, x1 = x.min() - 0.06 * ex, x.max() + 0.06 * ex
    y0, y1 = y.min() - 0.06 * ey, y.max() + 0.06 * ey
    nmax2 = max(16, len(P) // 2)
    # ⚠️ THE GRID MUST BE SCALED TO THE PATCH. A SQUARE grid on a 1.4 x 4.3
    # mm strip puts the poles 0.16 mm apart crosswise, where the data sits
    # 0.3 mm apart: under-determined, and the spline wobbles. With square
    # cells - same degrees of freedom, well distributed - the same patch
    # fits to 2.6e-04 with no wave.
    rat = math.sqrt(max(ex, 1e-9) / max(ey, 1e-9))
    best = None
    for n in (5, 7, 9, 11, 14, 18):
        nu = max(4, int(round(n * rat)))
        nv = max(4, int(round(n / rat)))
        if nu * nv > nmax2:
            break
        try:
            ku = _bs_knots(x0, x1, nu)
            kv = _bs_knots(y0, y1, nv)
            Bu = _bs_basis(x, ku, nu)
            Bv = _bs_basis(y, kv, nv)
            A = (Bu[:, :, None] * Bv[:, None, :]).reshape(len(P), nu * nv)
            Du, Dv = _d2_matrix(nu), _d2_matrix(nv)
            Rg = np.vstack([np.kron(Du, np.eye(nv)), np.kron(np.eye(nu), Dv)])
            sc = float(np.linalg.norm(A)) / max(float(np.linalg.norm(Rg)), 1e-12)
            AtA = A.T @ A
            Atz = A.T @ z
            RtR = Rg.T @ Rg
            I = np.eye(nu * nv)
        except Exception:
            continue
        for lam in (1e-4, 1e-3, 1e-5, 1e-6, 1e-8):
            try:
                c = np.linalg.solve(AtA + (lam * sc**2) * RtR + 1e-10 * I, Atz).reshape(nu, nv)
            except Exception:
                continue
            ff = FreeForm(C, X, Y, Z, c, ku, kv)
            if float(np.abs(ff.dist(P)).max()) > tol:
                continue
            if not _free_tame(ff, x, y, z):
                continue
            pr = Prim(FREE, C, Z, 0.0, 0.0, 0.0, float(np.abs(ff.dist(P)).max()))
            pr.free = ff
            if check is None or not len(check):
                return pr
            cr = float(np.abs(ff.dist(np.asarray(check, float))).max())
            if best is None or cr < best[0]:
                best = (cr, pr)
            if check_tol is not None and cr <= check_tol:
                return pr
    if best is not None and check_tol is not None and best[0] <= check_tol:
        return best[1]
    return None


def normal_deviation(prim: Prim, P: np.ndarray, Nrep: np.ndarray) -> float:
    """
    Mean angular deviation (degrees) between the faces' normals and the
    normal the primitive predicts at their centroids.

    ⚠️ THIS IS THE MAIN DISCRIMINANT, more than the points' distance.
    Real example: a hole wall tessellated in strips as tall as the whole
    part has vertices only on the two edges. Those points lie EXACTLY on
    both a cylinder and a sphere (r = sqrt(50^2+250^2)), and the
    positional residual can't decide: floating-point noise wins. The
    normals, though, are horizontal, while the sphere would want them
    tilted 79 degrees.

    ⚠️ It's evaluated on the VERTICES, not the face's centroid. In the
    case above, the strip's centroid falls on the sphere's equator, where
    the sphere's normal is also horizontal: the comparison there
    distinguishes nothing. On the vertices (at the two edges) the sphere
    is off by 79 degrees and gets rejected.
    """
    if len(P) == 0:
        return 0.0
    pn = prim.normal_at(P)
    d = np.abs(np.einsum("ij,ij->i", pn, Nrep))
    return float(np.degrees(np.arccos(np.clip(d, -1.0, 1.0))).mean())


def rank_primitives(P, Nrep, N, W, allow_sphere=True, allow_cone=True, flat_slope: float = 1e-3, max_ndev: float = 25.0, allow_torus: bool = True) -> List[Prim]:
    """
    All plausible primitives, SORTED by residual (progressive penalty to
    prefer the simpler models).

    ⚠️ Returns a LIST, not just the winner. On a small neighborhood the
    residual is nearly identical for different primitives: a hole wall can
    look like a sphere. If the first choice grows badly, the caller must
    be able to fall back to the second instead of wasting the faces.
    """
    cands = []
    pl = fit_plane(P)
    if pl:
        cands.append((pl.rms * 1.00, pl))
    ax = fit_axial(P, N, W)
    span_ = float(np.linalg.norm(P.max(axis=0) - P.min(axis=0))) or 1.0
    if ax is None or ax.rms > 1e-4 * span_:
        # thin band: the axis from the normals is degenerate, try the
        # revolution one instead (see fit_axial_rev)
        alt = fit_axial_rev(P, Nrep, np.ones(len(P)))
        if alt is not None and (ax is None or alt.rms < ax.rms):
            ax = alt
    if ax:
        if abs(ax.slope) < flat_slope:
            ax.slope = 0.0
            cands.append((ax.rms * 1.06, ax))
        elif allow_cone:
            cands.append((ax.rms * 1.15, ax))
    if allow_sphere:
        sp = fit_sphere(P)
        if sp:
            cands.append((sp.rms * 1.10, sp))
    # ⚠️ the torus must be fitted on points, not faces: it needs the
    #    per-vertex repeated normals (Nrep), not one per face.
    tr = fit_torus(P, Nrep, np.ones(len(P))) if allow_torus else None
    if tr:
        cands.append((tr.rms * 1.25, tr))  # penalty: only wins if it must
    keep = []
    for sc, pr in cands:
        nd = normal_deviation(pr, P, Nrep)
        if nd <= max_ndev:
            keep.append((sc, pr))
    keep.sort(key=lambda c: c[0])
    return [c[1] for c in keep]


_ALIVE: List[object] = []  # OCC builders to keep alive (SWIG/pybind: references)


def _keep(obj):
    _ALIVE.append(obj)
    return obj


_EDGE_ERR = (
    "ok",
    "point projection failed",
    "parameter out of range",
    "different points on a closed curve",
    "infinite parameter",
    "inconsistent point and parameter",
    "line through coincident points",
)


def _mk_vertex(p: np.ndarray, tol: float):
    """
    Vertex with an EXPLICIT tolerance.
    BRepBuilderAPI_MakeVertex uses Precision::Confusion() = 1e-7 mm. The
    recomputed nodes sit on the curves within the mesh's noise (~1e-5 mm on
    a 500 mm part): with 1e-7, every MakeEdge fails with "inconsistent
    point and parameter". The tolerance has to be sized to the real
    deviation.
    """
    v = TopoDS_Vertex()
    BRep_Builder().MakeVertex(v, _mk_pnt(p), float(max(tol, 1e-7)))
    return v


def _mk_dir(v: np.ndarray):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if n < 1e-12:
        raise ValueError("null direction")
    v = v / n
    return gp_Dir(float(v[0]), float(v[1]), float(v[2]))


def _mk_pnt(p: np.ndarray):
    p = np.asarray(p, dtype=float)
    return gp_Pnt(float(p[0]), float(p[1]), float(p[2]))


# --- 5C.1  nonlinear primitive refinement (Levenberg-Marquardt) ---------------
#
# ⚠️ THIS IS WHERE THE PREVIOUS DISASTER USED TO COME FROM.
# The algebraic fit (normals' covariance + Taubin) is only a SEED: on a
# 90-degree arc it gets the radius wrong by 0.1 mm and the axis by 1e-4
# rad. With that error the cylinder is no longer tangent to the plane and
# every boolean degenerates. The true minimization of the orthogonal
# distance brings the residual down to the mesh's noise level (1e-5 mm),
# i.e. the exact nominal value.

LM_MAX_PTS = 900
LM_MAX_SLOPE = 3.0  # = the fitters' max_slope: atan(3) = 71.6 degrees


def lm_refine(prim: "Prim", P: np.ndarray, W: Optional[np.ndarray] = None, iters: int = 80) -> "Prim":
    """Minimizes the orthogonal points-to-surface distance. Returns a new Prim."""
    if prim is None or len(P) < 4:
        return prim
    if prim.kind == PLANE:
        return prim  # SVD is already optimal for the plane
    P_all = P
    if len(P) > LM_MAX_PTS:  # on huge regions the sample is enough
        sel = np.linspace(0, len(P) - 1, LM_MAX_PTS).astype(int)
        P = P[sel]
        if W is not None:
            W = W[sel]

    w = np.ones(len(P)) if W is None else np.sqrt(np.maximum(W, 1e-12))
    w = w / w.mean()

    if prim.kind == SPHERE:
        c0 = prim.center.astype(float).copy()
        x = np.array([0.0, 0.0, 0.0, float(prim.r0)])

        def build(v):
            return Prim(SPHERE, c0 + v[:3], None, float(v[3]))

    elif prim.kind == TORUS:
        a0 = prim.axis / np.linalg.norm(prim.axis)
        u0, v0 = ortho_frame(a0)
        c0 = prim.center.astype(float).copy()
        x = np.array([0.0, 0.0, 0.0, 0.0, 0.0, float(prim.r0), float(prim.r1)])

        def build(v):
            a = a0 + v[0] * u0 + v[1] * v0
            a = a / np.linalg.norm(a)
            c = c0 + v[2] * u0 + v[3] * v0 + v[4] * a0
            return Prim(TORUS, c, a, float(v[5]), 0.0, float(v[6]))

    else:
        a0 = prim.axis / np.linalg.norm(prim.axis)
        u0, v0 = ortho_frame(a0)
        c0 = prim.center.astype(float).copy()
        is_cyl = abs(prim.slope) < 1e-9
        x = np.array([0.0, 0.0, 0.0, 0.0, float(prim.r0)]) if is_cyl else np.array([0.0, 0.0, 0.0, 0.0, float(prim.r0), float(prim.slope)])

        def build(v):
            a = a0 + v[0] * u0 + v[1] * v0
            a = a / np.linalg.norm(a)
            c = c0 + v[2] * u0 + v[3] * v0
            return Prim(AXIAL, c, a, float(v[4]), 0.0 if is_cyl else float(v[5]))

    def resid(v):
        try:
            return w * build(v).dist(P)
        except Exception:
            return np.full(len(P), 1e12)

    r = resid(x)
    cost = float(r @ r)
    lam = 1e-3
    m = len(x)
    converged = False
    for _ in range(iters):
        if converged:
            break
        J = np.empty((len(P), m))
        for k in range(m):
            h = 1e-7 * max(1.0, abs(x[k]))
            xp = x.copy()
            xp[k] += h
            J[:, k] = (resid(xp) - r) / h
        A = J.T @ J
        g = J.T @ r
        for _try in range(30):
            try:
                dx = np.linalg.solve(A + lam * np.diag(np.maximum(np.diag(A), 1e-12)), -g)
            except np.linalg.LinAlgError:
                lam *= 10.0
                continue
            xn = x + dx
            rn = resid(xn)
            cn = float(rn @ rn)
            if cn < cost:
                converged = (cost - cn) < 1e-9 * max(cost, 1e-300)
                x, r, cost = xn, rn, cn
                lam = max(lam * 0.3, 1e-12)
                break
            lam *= 10.0
            if lam > 1e10:
                break
        if lam > 1e10:
            break
        if np.linalg.norm(dx) < 1e-14 * (1.0 + np.linalg.norm(x)):
            break

    out = build(x)
    # ⚠️ LM CAN RUN AWAY TOWARD THE DEGENERATE CONE, and the script used to
    # die. A cone's residual is  (rho - r0 - s*t) / hypot(1, s):  as s -> oo
    # the cone becomes a PLANE and the cost keeps dropping, so the minimum
    # isn't where you'd expect. On test9, LM reached s = 77.8 (semi-angle
    # 89.26 degrees) with a NEGATIVE reference radius (-0.81); OCC refuses
    # to build a Geom_ConicalSurface with radius < 0 and raised
    # Standard_ConstructionError, which nothing caught: end of the run and
    # of all the work already done. The fitters search for the cone within
    # |slope| <= max_slope with a positive radius: the refinement has no
    # right to leave that window. If it does, the refinement gets thrown
    # away, not the starting primitive.
    if out is not None and out.kind == AXIAL and out.slope != 0.0:
        if not (np.isfinite(out.r0) and np.isfinite(out.slope) and out.r0 > 1e-9 and abs(out.slope) <= LM_MAX_SLOPE):
            out = prim
    out.rms = float(np.sqrt(np.mean(out.dist(P_all) ** 2)))
    return out


def _refit(region, verts, norms, areas, kind) -> Optional[Prim]:
    idx = list(region)
    P = np.vstack([verts[i] for i in idx])
    N = np.array([norms[i] for i in idx])
    W = np.array([max(areas[i], 1e-9) for i in idx])
    ok = np.einsum("ij,ij->i", N, N) > 0.5
    if kind == PLANE:
        return fit_plane(P)
    if kind == SPHERE:
        return fit_sphere(P)
    if ok.sum() < 3:
        return None
    if kind == TORUS:
        Wp = np.concatenate([np.full(len(verts[i]), max(areas[i], 1e-12)) for i in idx])
        Nr = np.vstack([np.tile(norms[i], (len(verts[i]), 1)) for i in idx])
        return fit_torus(P, Nr, Wp)
    return fit_axial(P, N[ok], W[ok])


# --- 5.7  full segmentation ----------------------------------------------------


def refit_exact(faces_idx, verts, norms, areas, kind) -> Optional["Prim"]:
    """_refit() followed by the LM refinement."""
    p = _refit(set(faces_idx), verts, norms, areas, kind)
    P = np.vstack([verts[i] for i in faces_idx])
    Wt = np.concatenate([np.full(len(verts[i]), max(areas[i], 1e-12) / max(len(verts[i]), 1)) for i in faces_idx])
    q = lm_refine(p, P, Wt) if p is not None else None
    # ⚠️ second path for cylinders and cones: if the normals-based fit is
    # missing or poor (thin band, normals on a short arc), retry with the
    # revolution axis, which on those bands is the only one that holds up.
    if kind == AXIAL:
        span = float(np.linalg.norm(P.max(axis=0) - P.min(axis=0))) or 1.0
        dq = float(np.abs(q.dist(P)).max()) if q is not None else math.inf
        # ⚠️ deliberately wide threshold: the second path costs an extra SVD
        # and LM per fit, and is only needed when the first one really went
        # wrong (on thin bands it's off by a hundredth, not a micron). At
        # 1e-5 every part paid double.
        if dq > 1e-4 * span:
            Nr = np.vstack([np.tile(norms[i], (len(verts[i]), 1)) for i in faces_idx])
            r = fit_axial_rev(P, Nr, Wt)
            if r is not None:
                r = lm_refine(r, P, Wt)
                if r is not None and float(np.abs(r.dist(P)).max()) < dq:
                    q = r
    if q is None:
        return None
    # ⚠️ TAPER BELOW THE NOISE: a reamed hole often comes out as a cone with
    # a taper of a few hundredths of a micron (ø7.9999-8.0000). It becomes
    # a CONICAL_SURFACE instead of a CYLINDRICAL_SURFACE, and the face
    # misbehaves with its neighbors. If the radius variation along the
    # region doesn't exceed the fit's residual, it's noise: zeroed out.
    if q is not None and q.kind == AXIAL and q.slope != 0.0:
        t = (P - q.center) @ q.axis
        span = float(t.max() - t.min())
        if abs(q.slope) * span <= max(2.0 * float(q.rms), 1e-4):
            tm = 0.5 * float(t.max() + t.min())
            q.center = q.center + tm * q.axis
            q.r0 = float(q.r0 + q.slope * tm)
            q.slope = 0.0
    return q


# --- 5C.1b  GLOBAL primitive extraction (Schnabel-style RANSAC) ---------------
#
# ⚠️ WHY IT'S NEEDED, AND WHY GROWING FROM SEEDS ISN'T ENOUGH.
# Seed growth decides the primitive's type by looking at two rings of
# triangles and then carries that choice forward. On a real part it
# produces dozens of tiny regions (3-5-facet cylinders with 3% coverage):
# these are the "overly dense zones" you see in the converted model.
# Global RANSAC reasons the other way around: it proposes a primitive from
# a handful of nearby faces and then counts HOW MANY faces of the whole
# part support it. The primitive with the biggest support wins, which by
# construction is the large, true one, not the fragment. It's the idea
# behind CGAL Shape Detection (Schnabel's Efficient RANSAC), here in numpy
# to avoid adding dependencies.


def refit_best(faces_idx, verts, norms, areas, current_kind: Optional[str] = None) -> Optional["Prim"]:
    """
    Reopens the question of the primitive's TYPE, not just its parameters.

    ⚠️ THIS IS WHY ROUND FILLETS STAYED TESSELLATED.
    The type is chosen at the SEED, i.e. on two rings of triangles: there,
    plane, cylinder, sphere and torus are all equivalent and the simplest
    one almost always wins. Then the region grows into a 300-facet
    toroidal band... and keeps getting refitted AS A SPHERE, because
    _refit preserves the kind. The residual stays a hundred times above
    the gate, the purge crumbles it every pass and those facets never
    become a surface.
    Here, once the region has grown, ALL types are tried again and the
    best one is kept (with a complexity penalty: at equal residual the
    simpler one wins).
    """
    idx = list(faces_idx)
    if len(idx) < 3:
        return None
    P = np.vstack([verts[i] for i in idx])
    if P.size == 0:
        return None
    Nrep = np.vstack([np.tile(norms[i], (len(verts[i]), 1)) for i in idx])
    Wp = np.concatenate([np.full(len(verts[i]), max(areas[i], 1e-12)) for i in idx])
    N = np.array([norms[i] for i in idx])
    W = np.array([max(areas[i], 1e-9) for i in idx])
    ok = np.einsum("ij,ij->i", N, N) > 0.5

    cands = []
    pl = fit_plane(P)
    if pl is not None:
        cands.append((1.00, pl))
    if ok.sum() >= 3:
        ax = fit_axial(P, N[ok], W[ok])
        if ax is not None:
            cands.append((1.06 if abs(ax.slope) < 1e-3 else 1.15, ax))
    sp = fit_sphere(P)
    if sp is not None:
        cands.append((1.10, sp))
    tr = fit_torus(P, Nrep, Wp)
    if tr is not None:
        cands.append((1.25, tr))
    if not cands:
        return None

    best, bs = None, math.inf
    for pen, p in cands:
        p2 = lm_refine(p, P, Wp)
        if p2 is None:
            continue
        d = p2.dist(P)
        rms = float(np.sqrt(np.mean(d**2)))
        # small bonus for whatever was already the chosen type: avoids oscillation
        s = rms * pen * (0.95 if p2.kind == current_kind else 1.0)
        if s < bs:
            best, bs = p2, s
    return best


# --- 5C.2b  merging co-surface regions -----------------------------------------
#
# ⚠️ WITHOUT THIS THE MODEL FALLS APART. Growth starts from different seeds
# and a long fillet often ends up as TWO regions with the exact same
# primitive. They'd become two overlapping cylindrical faces separated by
# a nonexistent edge. Here, neighboring regions lying on the SAME surface
# are merged and refitted.


def prims_same(pa: "Prim", pb: "Prim", tol_len: float, cos_ang: float) -> bool:
    if pa is None or pb is None or pa.kind != pb.kind:
        return False
    if pa.kind == FREE:
        return False  # two free-form shapes never merge on sight
    if pa.kind == PLANE:
        if float(pa.axis @ pb.axis) < cos_ang:
            return False
        return abs(float((pb.center - pa.center) @ pa.axis)) <= tol_len
    if pa.kind == SPHERE:
        return float(np.linalg.norm(pa.center - pb.center)) <= tol_len and abs(pa.r0 - pb.r0) <= tol_len
    if abs(float(pa.axis @ pb.axis)) < cos_ang:
        return False
    d = pb.center - pa.center
    if float(np.linalg.norm(d - float(d @ pa.axis) * pa.axis)) > tol_len:
        return False
    if pa.kind == TORUS:
        return abs(pa.r0 - pb.r0) <= tol_len and abs(pa.r1 - pb.r1) <= tol_len and abs(float(d @ pa.axis)) <= tol_len
    if abs(pa.slope - abs(pb.slope) * (1.0 if float(pa.axis @ pb.axis) > 0 else -1.0)) > 1e-4:
        return False
    t = float(d @ pa.axis)
    return abs((pa.r0 + pa.slope * t) - pb.r0) <= tol_len


# =============================================================================
# 6. ANALYTIC CURVES and surface-surface intersections
# =============================================================================
class _Curve:
    period = None

    def param(self, P: np.ndarray) -> np.ndarray:  # pragma: no cover
        raise NotImplementedError

    def point(self, t: np.ndarray) -> np.ndarray:  # pragma: no cover
        raise NotImplementedError

    def dist(self, P: np.ndarray) -> np.ndarray:
        return np.linalg.norm(P - self.point(self.param(P)), axis=1)

    def to_geom(self):  # pragma: no cover
        raise NotImplementedError


@dataclass
class CLine(_Curve):
    p: np.ndarray
    d: np.ndarray

    def param(self, P):
        return np.atleast_1d((np.atleast_2d(P) - self.p) @ self.d)

    def point(self, t):
        t = np.atleast_1d(t)
        return self.p + np.outer(t, self.d)

    def to_geom(self):
        return Geom_Line(_mk_pnt(self.p), _mk_dir(self.d))

    def label(self):
        return "line"


@dataclass
class CCircle(_Curve):
    c: np.ndarray
    n: np.ndarray
    u: np.ndarray
    v: np.ndarray
    r: float
    period = 2.0 * math.pi

    def param(self, P):
        d = np.atleast_2d(P) - self.c
        return np.arctan2(d @ self.v, d @ self.u)

    def point(self, t):
        t = np.atleast_1d(t)
        return self.c + self.r * (np.cos(t)[:, None] * self.u + np.sin(t)[:, None] * self.v)

    def to_geom(self):
        ax2 = gp_Ax2(_mk_pnt(self.c), _mk_dir(self.n), _mk_dir(self.u))
        return Geom_Circle(ax2, float(self.r))

    def label(self):
        return f"circle r={self.r:.4f}"


@dataclass
class CEllipse(_Curve):
    c: np.ndarray
    n: np.ndarray
    u: np.ndarray  # MAJOR axis
    v: np.ndarray  # minor axis
    a: float
    b: float
    period = 2.0 * math.pi

    def param(self, P):
        d = np.atleast_2d(P) - self.c
        return np.arctan2((d @ self.v) / max(self.b, 1e-12), (d @ self.u) / max(self.a, 1e-12))

    def point(self, t):
        t = np.atleast_1d(t)
        return self.c + (self.a * np.cos(t))[:, None] * self.u + (self.b * np.sin(t))[:, None] * self.v

    def to_geom(self):
        ax2 = gp_Ax2(_mk_pnt(self.c), _mk_dir(self.n), _mk_dir(self.u))
        return Geom_Ellipse(ax2, float(self.a), float(self.b))

    def label(self):
        return f"ellipse {self.a:.4f}x{self.b:.4f}"


@dataclass
class CHyperbola(_Curve):
    """
    Hyperbola branch, parametrized like Geom_Hyperbola:
        P(t) = c + a*cosh(t)*u + b*sinh(t)*v
    `u` points toward the VERTEX of the represented branch (the other
    branch is the same curve with u flipped in sign). Not periodic.
    """

    c: np.ndarray
    n: np.ndarray
    u: np.ndarray  # transverse axis, toward the vertex
    v: np.ndarray  # conjugate axis
    a: float
    b: float
    period = None

    def param(self, P):
        P = np.atleast_2d(P)
        d = P - self.c
        t = np.arcsinh((d @ self.v) / max(self.b, 1e-12))
        # ⚠️ the closest point is NOT the one with the same ordinate: on a
        # steep branch the difference is tens of microns, i.e. more than
        # the tolerance the curve is judged by. Three Newton steps on
        # (P - S(t)) . S'(t) = 0 are enough and cost nothing.
        for _ in range(3):
            ch, sh = np.cosh(t), np.sinh(t)
            S = self.c + (self.a * ch)[:, None] * self.u + (self.b * sh)[:, None] * self.v
            T = (self.a * sh)[:, None] * self.u + (self.b * ch)[:, None] * self.v
            T2 = (self.a * ch)[:, None] * self.u + (self.b * sh)[:, None] * self.v
            W = P - S
            f = np.einsum("ij,ij->i", W, T)
            fp = np.einsum("ij,ij->i", W, T2) - np.einsum("ij,ij->i", T, T)
            step = np.where(np.abs(fp) > 1e-30, f / fp, 0.0)
            t = t - np.clip(step, -1.0, 1.0)
        return t

    def point(self, t):
        t = np.atleast_1d(np.asarray(t, float))
        return self.c + (self.a * np.cosh(t))[:, None] * self.u + (self.b * np.sinh(t))[:, None] * self.v

    def to_geom(self):
        ax2 = gp_Ax2(_mk_pnt(self.c), _mk_dir(self.n), _mk_dir(self.u))
        return Geom_Hyperbola(ax2, float(self.a), float(self.b))

    def label(self):
        return f"hyperbola {self.a:.4f}x{self.b:.4f}"


def _hyperbola(c, n, u, a: float, b: float) -> Optional[CHyperbola]:
    n = np.asarray(n, float)
    nn = np.linalg.norm(n)
    if nn < 1e-12 or not (a > 0.0 and b > 0.0):
        return None
    n = n / nn
    u = np.asarray(u, float)
    u = u - float(u @ n) * n
    nu = np.linalg.norm(u)
    if nu < 1e-12:
        return None
    u = u / nu
    return CHyperbola(np.asarray(c, float), n, u, np.cross(n, u), float(a), float(b))


def _circle(c, n, r) -> CCircle:
    n = n / np.linalg.norm(n)
    u, v = ortho_frame(n)
    return CCircle(np.asarray(c, float), n, u, np.cross(n, u), float(r))


def _ellipse(c, n, u, a: float, b: float) -> Optional[CEllipse]:
    """
    ⚠️ Orthonormalization is MANDATORY before building the ellipse.
    Geom_Ellipse requires the major axis >= minor and derives Y = Z x X: if
    the minor unit vector passed here were the opposite one, the Python
    parameter and OCC's would have opposite signs and the trim would take
    the wrong arc.
    """
    n = np.asarray(n, float)
    nn = np.linalg.norm(n)
    u = np.asarray(u, float) - 0.0
    if nn < 1e-12:
        return None
    n = n / nn
    u = u - float(u @ n) * n
    nu = np.linalg.norm(u)
    if nu < 1e-12:
        return None
    u = u / nu
    if b > a:
        a, b = b, a
        u = np.cross(n, u)
    return CEllipse(np.asarray(c, float), n, u, np.cross(n, u), float(a), float(b))


# --- 5C.6a  the cone: planar sections and coaxial circles ---------------------
# ⚠️ THE CONE USED TO BE THE BLACK HOLE OF INTERSECTIONS. Until now the only
# section computed was cone x plane ORTHOGONAL to the axis (a circle). But
# a countersink or a conical chamfer almost always borders something else:
# the counterbore it flares from (coaxial cylinder), the fillet that
# closes it off (coaxial torus), a slot's wall cutting across it (plane
# PARALLEL to the axis -> hyperbola), a tilted face (oblique plane ->
# ellipse). Without these curves the boundary stayed the mesh's polyline,
# and when the "circle through the vertices" fallback invented a circle
# outside the cone, the face came out SelfIntersectingWire and the
# countersink went back to being tessellated.


def _cone_frame(cone: "Prim"):
    """(apex, axis oriented toward the opening, semi-angle) of the cone."""
    a = cone.axis / np.linalg.norm(cone.axis)
    sg = 1.0 if cone.slope > 0 else -1.0
    a = sg * a
    s = abs(cone.slope)
    apex = cone.center - sg * (cone.r0 / s) * a
    return apex, a, math.atan(s)


def _cone_plane_curves(cone: "Prim", plane: "Prim", tol: float) -> List["_Curve"]:
    """
    Section of a cone with an arbitrary plane: circle, ellipse or
    hyperbola. In the plane, with the origin at the apex's projection and
    x along the axis's projection, the conic is
        (sin^2(psi) - cos^2(alpha)) x^2 - cos^2(alpha) y^2
        + 2 D cos(psi) sin(psi) x + D^2 (cos^2(psi) - cos^2(alpha)) = 0
    with psi = angle between the plane's normal and the axis, D =
    apex-to-plane distance. The sign of (sin^2(psi) - cos^2(alpha))
    decides: <0 ellipse, >0 hyperbola, ~0 parabola (skipped: a
    measure-zero case, stays the polyline).
    """
    apex, a, alpha = _cone_frame(cone)
    n = plane.axis / np.linalg.norm(plane.axis)
    D = float((plane.center - apex) @ n)
    if n @ a < 0.0:
        n, D = -n, -D
    if abs(D) <= max(10.0 * tol, 1e-9):
        return []  # plane through the apex: degenerate
    cps = float(n @ a)
    W = n - cps * a
    sps = float(np.linalg.norm(W))
    ca, sa = math.cos(alpha), math.sin(alpha)
    out: List[_Curve] = []
    P0 = apex + D * n  # apex's foot on the plane
    if sps < 1e-9:  # ORTHOGONAL plane -> circle
        r = abs(D) * sa / ca
        if r > 10.0 * tol:
            out.append(_circle(P0, n, r))
        return out
    E1 = W / sps
    A2 = sps * sps - ca * ca
    if abs(A2) < 1e-6:
        return []  # parabola
    x0 = -D * cps * sps / A2
    C = P0 + x0 * E1
    if A2 < 0.0:  # ELLIPSE
        ax = abs(D) * ca * sa / abs(A2)
        ay = abs(D) * sa / math.sqrt(-A2)
        if min(ax, ay) > 10.0 * tol:
            el = _ellipse(C, n, E1, ax, ay)
            if el is not None:
                out.append(el)
        return out
    # HYPERBOLA: two branches, one per nappe. Both are emitted and the
    # residual on the vertices decides (the wrong branch sits on the other
    # nappe, far away).
    ax = abs(D) * ca * sa / A2
    ay = abs(D) * sa / math.sqrt(A2)
    if min(ax, ay) <= 10.0 * tol:
        return out
    for sgn in (+1.0, -1.0):
        hy = _hyperbola(C, n, sgn * E1, ax, ay)
        if hy is None:
            continue
        # keeps only the branch on the cone's REAL nappe (the one with radii > 0)
        vtx = hy.point(np.array([0.0]))[0]
        if float((vtx - apex) @ a) > 0.0:
            out.append(hy)
    return out


def _coax_t(pa: "Prim", pb: "Prim", tol: float):
    """If pb is coaxial with pa: (t of pb's center on pa's axis, axes' sign)."""
    a = pa.axis / np.linalg.norm(pa.axis)
    if pb.kind == SPHERE:
        d = pb.center - pa.center
        t = float(d @ a)
        if float(np.linalg.norm(d - t * a)) > max(tol, 1e-6):
            return None
        return t, 1.0
    b = pb.axis / np.linalg.norm(pb.axis)
    cs = float(a @ b)
    if abs(cs) < 1.0 - 1e-6:
        return None
    d = pb.center - pa.center
    t = float(d @ a)
    if float(np.linalg.norm(d - t * a)) > max(tol, 1e-6):
        return None
    return t, (1.0 if cs > 0 else -1.0)


def _cone_coax_circles(cone: "Prim", other: "Prim", tol: float) -> List["_Curve"]:
    """
    Cone x COAXIAL surface of revolution -> circles, at the points where the
    two radius(t) profiles meet. Covers countersink-hole, chamfer-chamfer,
    chamfer-fillet and chamfer-cap, which are everyday cases on a turned or
    drilled part.
    """
    co = _coax_t(cone, other, tol)
    if co is None:
        return []
    tc, sg = co
    a = cone.axis / np.linalg.norm(cone.axis)
    s, r0 = cone.slope, cone.r0  # radius(t) = r0 + s*t
    roots: List[float] = []
    if other.kind == AXIAL:
        s2 = other.slope * sg
        b0 = other.r0 - s2 * tc  # radius(t) = b0 + s2*t
        den = s - s2
        if abs(den) < 1e-9:
            return []  # parallel profiles
        roots.append((b0 - r0) / den)
    elif other.kind in (SPHERE, TORUS):
        # (K + s t)^2 + (t - tc)^2 = rr^2   with K = r0 - Rmajor, rr = radius
        K = r0 if other.kind == SPHERE else r0 - other.r0
        rr = other.r0 if other.kind == SPHERE else other.r1
        qa = s * s + 1.0
        qb = 2.0 * (K * s - tc)
        qc = K * K + tc * tc - rr * rr
        disc = qb * qb - 4.0 * qa * qc
        if disc < 0.0:
            if disc < -(4.0 * qa * max(tol, 1e-9) * rr):
                return []
            disc = 0.0  # tangency
        sq = math.sqrt(disc)
        roots.append((-qb + sq) / (2.0 * qa))
        if sq > 1e-12:
            roots.append((-qb - sq) / (2.0 * qa))
    else:
        return []
    out: List[_Curve] = []
    for t in roots:
        r = r0 + s * t
        if r > 10.0 * tol:
            out.append(_circle(cone.center + t * a, a, r))
    return out


# --- 5C.6  EXACT intersection between two primitives ---------------------------


def surf_surf_curves(pa: "Prim", pb: "Prim", tol: float) -> List[_Curve]:
    """
    Analytic intersection curves between two primitives.
    Returns ALL possible solutions (0, 1 or 2): the comparison against the
    mesh's vertices will pick the right one.
    """
    if pa is None or pb is None:
        return []
    if pa.kind == FREE or pb.kind == FREE:
        return []  # free-form: no exact curve
    ka, kb = pa.kind, pb.kind
    if (ka, kb) in ((AXIAL, PLANE), (SPHERE, PLANE), (SPHERE, AXIAL), (TORUS, PLANE), (TORUS, AXIAL), (TORUS, SPHERE)):
        pa, pb = pb, pa
        ka, kb = pa.kind, pb.kind
    out: List[_Curve] = []
    try:
        # ---------- plane x plane -> line ----------
        if ka == PLANE and kb == PLANE:
            d = np.cross(pa.axis, pb.axis)
            nd = np.linalg.norm(d)
            if nd < 1e-9:
                return []
            d /= nd
            A = np.vstack([pa.axis, pb.axis, d])
            rhs = np.array([pa.axis @ pa.center, pb.axis @ pb.center, 0.0])
            p = np.linalg.solve(A, rhs)
            out.append(CLine(p, d))

        # ---------- plane x sphere -> circle ----------
        elif ka == PLANE and kb == SPHERE:
            h = float((pb.center - pa.center) @ pa.axis)
            if abs(h) > pb.r0 + tol:
                return []
            rr = math.sqrt(max(pb.r0**2 - h**2, 0.0))
            if rr < 10 * tol:
                return []
            out.append(_circle(pb.center - h * pa.axis, pa.axis, rr))

        # ---------- plane x cylinder/cone ----------
        elif ka == PLANE and kb == AXIAL:
            # ⚠️ No exclusive branch: ALL plausible interpretations are
            # generated and the residual on the vertices will choose. A
            # fitted axis is parallel to the plane to within 1e-6 rad,
            # never within 1e-9: a rigid test would pick the ellipse with a
            # 4e6 mm semi-axis instead of the tangency line.
            n, a = pa.axis, pb.axis
            cph = float(a @ n)
            if abs(cph) > 1.0 - 1e-6:  # ORTHOGONAL plane
                tstar = float((pa.center - pb.center) @ n) / cph
                r = pb.r0 + pb.slope * tstar
                if r > 10 * tol:
                    out.append(_circle(pb.center + tstar * a, a, r))
            if abs(pb.slope) >= 1e-9:  # CONE: generic conic
                out += _cone_plane_curves(pb, pa, tol)
            if abs(pb.slope) < 1e-9:
                if abs(cph) < 1e-2:  # ~PARALLEL plane
                    h = float((pb.center - pa.center) @ n)
                    foot = pb.center - h * n
                    w = np.cross(n, a)
                    nw = np.linalg.norm(w)
                    if nw > 1e-12 and abs(h) <= pb.r0 * 1.02 + tol:
                        w /= nw
                        s = math.sqrt(max(pb.r0**2 - h**2, 0.0))
                        out.append(CLine(foot, a))  # TANGENCY
                        if s > max(tol, 1e-9):
                            out.append(CLine(foot + s * w, a))
                            out.append(CLine(foot - s * w, a))
                if 1e-9 < abs(cph) < 1.0 - 1e-9:  # OBLIQUE plane
                    tstar = float((pa.center - pb.center) @ n) / cph
                    c = pb.center + tstar * a
                    umaj = a - cph * n
                    nu = np.linalg.norm(umaj)
                    if nu > 1e-12:
                        umaj /= nu
                        el = _ellipse(c, n, umaj, pb.r0 / abs(cph), pb.r0)
                        if el is not None:
                            out.append(el)

        # ---------- plane x torus ----------
        elif ka == PLANE and kb == TORUS:
            n, a = pa.axis, pb.axis
            cph = float(a @ n)
            if abs(cph) > 1.0 - 1e-6:  # plane ORTHOGONAL to the axis
                t = float((pa.center - pb.center) @ n) / cph
                if abs(t) <= pb.r1 + tol:
                    s = math.sqrt(max(pb.r1**2 - t**2, 0.0))
                    ctr = pb.center + t * a
                    for rr in (pb.r0 + s, pb.r0 - s):
                        if rr > 10 * tol:
                            out.append(_circle(ctr, a, rr))
            elif abs(cph) < 1e-2:  # MERIDIAN plane
                h = float((pb.center - pa.center) @ n)
                if abs(h) <= max(tol, 1e-6):
                    w = np.cross(n, a)
                    nw = np.linalg.norm(w)
                    if nw > 1e-12:
                        w /= nw
                        for sg in (+1.0, -1.0):
                            out.append(_circle(pb.center + sg * pb.r0 * w, n, pb.r1))

        # ---------- cylinder x torus ----------
        elif ka == AXIAL and kb == TORUS:
            if abs(pa.slope) > 1e-9:
                return _cone_coax_circles(pa, pb, tol)
            d = pb.center - pa.center
            t0 = float(d @ pa.axis)
            off = d - t0 * pa.axis
            # ⚠️ A STRAIGHT FILLET CONTINUING INTO A CORNER FILLET: ORTHOGONAL
            # axes, same tube radius, cylinder axis at distance R from the
            # torus center. The two surfaces are tangent along the torus's
            # meridian: a circle of radius r in the plane orthogonal to the
            # cylinder's axis, passing through the torus's center. Without
            # this case the boundary stayed the mesh's staircase.
            if abs(float(pa.axis @ pb.axis)) < 1e-3:
                rt = pb.r1
                if abs(pa.r0 - rt) <= max(tol, 1e-3 * pa.r0) and abs(float(np.linalg.norm(off)) - pb.r0) <= max(10.0 * tol, 1e-3 * pb.r0):
                    cj = pa.center + t0 * pa.axis  # point of the cylinder's axis closest to C
                    out.append(_circle(cj, pa.axis, pa.r0))
                return out
            if float(np.linalg.norm(off)) > max(tol, 1e-6):
                return []
            if abs(float(pa.axis @ pb.axis)) < 1.0 - 1e-6:
                return []
            dr = pa.r0 - pb.r0
            if abs(dr) > pb.r1 + tol:
                return []
            s = math.sqrt(max(pb.r1**2 - dr**2, 0.0))
            for sg in (0.0,) if s <= max(tol, 1e-9) else (+1.0, -1.0):
                out.append(_circle(pa.center + (t0 + sg * s) * pa.axis, pa.axis, pa.r0))

        # ---------- torus x torus (coaxial) -> circles ----------
        elif ka == TORUS and kb == TORUS:
            if abs(float(pa.axis @ pb.axis)) < 1.0 - 1e-6:
                return []
            d = pb.center - pa.center
            t0 = float(d @ pa.axis)
            if float(np.linalg.norm(d - t0 * pa.axis)) > max(tol, 1e-6):
                return []
            return []  # quartic: handled by the point-based fit

        # ---------- cylinder x sphere (coaxial) -> circles ----------
        elif ka == AXIAL and kb == SPHERE:
            if abs(pa.slope) > 1e-9:
                return _cone_coax_circles(pa, pb, tol)
            d = pb.center - pa.center
            t0 = float(d @ pa.axis)
            off = float(np.linalg.norm(d - t0 * pa.axis))
            if off > max(tol, 1e-6):
                return []
            # ⚠️ sqrt(rs^2 - rc^2) is unusable as a tangency test: with rs
            # and rc equal to within 1e-6 mm, the radicand is -8e-6 and the
            # tangent branch never fires. The comparison must be on the RADII.
            dr = pb.r0 - pa.r0
            if abs(dr) <= max(tol, 1e-3 * pa.r0):
                out.append(_circle(pa.center + t0 * pa.axis, pa.axis, pa.r0))
            if dr > 0.0:
                h = math.sqrt(max(pb.r0**2 - pa.r0**2, 0.0))
                if h > max(tol, 1e-9):
                    out.append(_circle(pa.center + (t0 + h) * pa.axis, pa.axis, pa.r0))
                    out.append(_circle(pa.center + (t0 - h) * pa.axis, pa.axis, pa.r0))

        # ---------- sphere x sphere -> circle ----------
        elif ka == SPHERE and kb == SPHERE:
            d = pb.center - pa.center
            dd = float(np.linalg.norm(d))
            if dd < 1e-9:
                return []
            x = (dd**2 + pa.r0**2 - pb.r0**2) / (2.0 * dd)
            rr2 = pa.r0**2 - x**2
            if rr2 < -(tol**2):
                return []
            out.append(_circle(pa.center + x * (d / dd), d / dd, math.sqrt(max(rr2, 0.0))))

        # ---------- cylinder x cylinder (parallel axes) -> lines ----------
        elif ka == AXIAL and kb == AXIAL:
            if abs(pa.slope) > 1e-9 or abs(pb.slope) > 1e-9:
                # at least one is a cone: coaxial -> circles, otherwise
                # quartic, handled by the point-based fit
                cn, ot = (pa, pb) if abs(pa.slope) > 1e-9 else (pb, pa)
                return _cone_coax_circles(cn, ot, tol)
            if abs(float(pa.axis @ pb.axis)) < 1.0 - 1e-9:
                # ⚠️ SHARP EDGE BETWEEN TWO FILLETS. Where two fillets of the
                # same radius meet with no corner sphere, the two cylinders
                # have INTERSECTING axes: the intersection isn't a quartic
                # but TWO planar ELLIPSES (Steinmetz). Without this case,
                # that edge stays polygonal and the face doesn't close.
                if abs(pa.r0 - pb.r0) > max(tol, 1e-3 * pa.r0):
                    return []
                a1 = pa.axis
                a2 = pb.axis * (1.0 if float(pa.axis @ pb.axis) > 0 else -1.0)
                w0 = pa.center - pb.center
                b_ = float(a1 @ a2)
                den = 1.0 - b_ * b_
                if abs(den) < 1e-12:
                    return []
                dd_ = float(a1 @ w0)
                ee_ = float(a2 @ w0)
                s_ = (b_ * ee_ - dd_) / den
                t_ = (ee_ - b_ * dd_) / den
                Q1 = pa.center + s_ * a1
                Q2 = pb.center + t_ * a2
                if float(np.linalg.norm(Q1 - Q2)) > max(10.0 * tol, 1e-6 * pa.r0):
                    return []
                Q = 0.5 * (Q1 + Q2)
                w = np.cross(a1, a2)
                nw = np.linalg.norm(w)
                if nw < 1e-12:
                    return []
                w /= nw
                th = math.acos(max(-1.0, min(1.0, b_)))
                r = 0.5 * (pa.r0 + pb.r0)
                for sgn in (-1.0, +1.0):
                    nrm = a1 + sgn * a2
                    maj = a1 - sgn * a2
                    if np.linalg.norm(nrm) < 1e-9 or np.linalg.norm(maj) < 1e-9:
                        continue
                    nrm = nrm / np.linalg.norm(nrm)
                    maj = maj / np.linalg.norm(maj)
                    half = math.cos(th / 2.0) if sgn > 0 else math.sin(th / 2.0)
                    if abs(half) < 1e-9:
                        continue
                    el = _ellipse(Q, nrm, maj, r / abs(half), r)
                    if el is not None:
                        out.append(el)
                return out
            a = pa.axis
            d = pb.center - pa.center
            perp = d - float(d @ a) * a
            dd = float(np.linalg.norm(perp))
            if dd < 1e-9:
                return []
            ex = perp / dd
            ey = np.cross(a, ex)
            x = (dd**2 + pa.r0**2 - pb.r0**2) / (2.0 * dd)
            y2 = pa.r0**2 - x**2
            if y2 < -(tol**2):
                return []
            y = math.sqrt(max(y2, 0.0))
            base = pa.center + x * ex
            if y <= max(tol, 1e-9):
                out.append(CLine(base, a))
            else:
                out.append(CLine(base + y * ey, a))
                out.append(CLine(base - y * ey, a))
    except Exception as e:
        Log.debug(f"intersection {ka}x{kb} failed: {e}")
        return []
    return out


# --- 5C.7  fallback fit (when the analytic intersection isn't available) ------


def fit_curve(P: np.ndarray, tol: float, arc_tol: Optional[float] = None) -> Optional[_Curve]:
    """Line or circle through the vertices (tol on the vertices, arc_tol on the arc between them)."""
    if arc_tol is None:
        arc_tol = tol
    if len(P) < 2:
        return None
    c = P.mean(axis=0)
    Q = P - c
    _, S, Vt = np.linalg.svd(Q, full_matrices=False)
    line = CLine(c, Vt[0])
    if float(np.abs(line.dist(P)).max()) <= tol:
        return line
    if len(P) < 4:
        return None
    n = Vt[-1]
    if float(np.abs(Q @ n).max()) > tol:
        return None
    u, v = ortho_frame(n)
    try:
        cx, cy, r = taubin_circle(Q @ u, Q @ v)
    except Exception:
        return None
    if not (np.isfinite(cx) and np.isfinite(cy) and np.isfinite(r)) or r <= 1e-9:
        return None
    circ = CCircle(c + cx * u + cy * v, n, u, np.cross(n, u), float(r))
    if float(np.abs(circ.dist(P)).max()) <= tol and arc_deviation(circ, P) <= arc_tol:
        return circ
    return None


def arc_deviation(cv: _Curve, P: np.ndarray, dens: int = 6) -> float:
    """
    Deviation of the ARC TRAVELED from the vertices' polyline, not just the
    vertices.

    ⚠️ THIS IS THE CHECK THAT WAS MISSING, AND IT'S WHY CURVES USED TO
    STICK OUT. Judging a curve only AT THE VERTICES is like judging a
    bridge by looking at the piers: an ellipse with an 889 mm semi-axis on
    a 49 mm part can pass exactly through every vertex of a short chain and
    then, BETWEEN one vertex and the next, drift millimeters outside the
    part. Here the curve is sampled densely along the whole arc and
    compared against the vertices' polyline: if it strays from it, it's
    rejected.
    """
    P = np.atleast_2d(P)
    if len(P) < 2:
        return 0.0
    try:
        t = np.asarray(cv.param(P), dtype=float)
        if cv.period:
            t = np.unwrap(t)
        ts = [np.linspace(t[k], t[k + 1], dens, endpoint=False) for k in range(len(t) - 1)]
        ts.append(np.array([t[-1]]))
        S_ = cv.point(np.concatenate(ts))
    except Exception:
        return math.inf
    A = P[:-1]
    B = P[1:]
    AB = B - A
    L2 = np.maximum(np.einsum("ij,ij->i", AB, AB), 1e-30)
    worst = 0.0
    for s in S_:
        u = np.clip(np.einsum("ij,ij->i", s - A, AB) / L2, 0.0, 1.0)
        d = float(np.linalg.norm(A + u[:, None] * AB - s, axis=1).min())
        if d > worst:
            worst = d
    return worst


def choose_curve(cands: List[_Curve], P: np.ndarray, tol: float, scale: float = 0.0, arc_tol: Optional[float] = None) -> Optional[_Curve]:
    """
    tol     : max deviation of the VERTICES from the curve.
    arc_tol : max deviation of the ARC from the vertices' polyline (must
              allow for the mesh chords' sag, otherwise no circle fits).
    """
    if arc_tol is None:
        arc_tol = tol
    best, bs = None, math.inf
    for cv in cands:
        if scale > 0:
            big = getattr(cv, "a", None)
            if big is None:
                big = getattr(cv, "r", None)
            if big is not None and float(big) > 20.0 * scale:
                continue  # huge primitive: a fit artifact
        try:
            s = float(np.abs(cv.dist(P)).max())
        except Exception:
            continue
        if s > tol:
            continue
        s2 = arc_deviation(cv, P)
        if s2 > arc_tol:
            continue
        sc = max(s, s2)
        if sc < bs:
            best, bs = cv, sc
    return best


# --- 5C.8  nodes: exact position as the intersection of the incident curves ---


# =============================================================================
# 7. TOPOLOGY INDEX of the current shape
# =============================================================================


class Topo:
    """
    Faces, edges, vertices of the shape with integer indices stable as long
    as the shape doesn't change. Rebuilt (costs ~0.1 s on 5,000 faces)
    after every accepted replacement.
    """

    __slots__ = ("shape", "faces", "fmap", "nF", "verts", "norms", "areas", "cents", "nverts", "edges", "emap", "nE", "e_faces", "f_edges", "e_verts", "vmap", "nV", "vpos", "adj", "planar")

    def __init__(self, shape, prev: Optional["Topo"] = None):
        self.shape = shape
        fmap = TopTools_IndexedMapOfShape()
        te_MapShapes(shape, TopAbs_FACE, fmap)
        self.fmap = fmap
        self.nF = _size(fmap)
        self.faces = [td_Face(fmap.FindKey(i)) for i in range(1, self.nF + 1)]

        vmap = TopTools_IndexedMapOfShape()
        te_MapShapes(shape, TopAbs_VERTEX, vmap)
        self.vmap = vmap
        self.nV = _size(vmap)
        self.vpos = np.array([vpos(vmap.FindKey(j)) for j in range(1, self.nV + 1)], dtype=float).reshape(-1, 3)

        self.verts, self.cents, self.norms, self.areas, self.nverts, self.planar = {}, {}, {}, {}, {}, {}
        for i, f in enumerate(self.faces):
            # ⚠️ faces untouched by the last replacement are the same
            # (IsSame): their data is copied from the previous index
            if prev is not None:
                k = prev.fmap.FindIndex(f) - 1
                if k >= 0 and prev.faces[k].Orientation() == f.Orientation():
                    self.verts[i] = prev.verts[k]
                    self.nverts[i] = prev.nverts[k]
                    self.cents[i] = prev.cents[k]
                    self.areas[i] = prev.areas[k]
                    self.planar[i] = prev.planar[k]
                    self.norms[i] = prev.norms[k]
                    continue
            m = TopTools_IndexedMapOfShape()
            te_MapShapes(f, TopAbs_VERTEX, m)
            idx = [vmap.FindIndex(m.FindKey(k)) - 1 for k in range(1, _size(m) + 1)]
            idx = [j for j in idx if j >= 0]
            V = self.vpos[idx] if idx else np.zeros((0, 3))
            self.verts[i] = V
            self.nverts[i] = len(idx)
            self.cents[i] = V.mean(axis=0) if len(V) else np.zeros(3)
            self.areas[i] = face_area(f)
            n = face_plane_normal(f)
            self.planar[i] = n is not None
            self.norms[i] = n if n is not None else np.zeros(3)

        ef, _, e_faces = face_edges_map(shape)
        self.emap = ef
        self.nE = _size(ef)
        # ⚠️ the map returns every edge with the orientation of the FIRST
        # use encountered: here they're normalized to FORWARD, so "fwd"
        # always means "from FirstVertex to LastVertex" and Reversed() does
        # what it says.
        self.edges = [td_Edge(ef.FindKey(k).Oriented(TopAbs_FORWARD)) for k in range(1, self.nE + 1)]
        self.e_faces = [sorted(fs) for fs in e_faces]
        self.f_edges = [[] for _ in range(self.nF)]
        for k, fs in enumerate(self.e_faces):
            for i in fs:
                self.f_edges[i].append(k)
        self.e_verts = []
        for e in self.edges:
            a = vmap.FindIndex(te_FirstVertex(e)) - 1
            b = vmap.FindIndex(te_LastVertex(e)) - 1
            self.e_verts.append((a, b))
        self.adj = [set() for _ in range(self.nF)]
        for fs in self.e_faces:
            for a in range(len(fs)):
                for b in range(a + 1, len(fs)):
                    self.adj[fs[a]].add(fs[b])
                    self.adj[fs[b]].add(fs[a])

    def face_index(self, face) -> int:
        return self.fmap.FindIndex(face) - 1

    def edge_index(self, edge) -> int:
        return self.emap.FindIndex(edge) - 1

    def vertex_index(self, v) -> int:
        return self.vmap.FindIndex(v) - 1


def plane_prim_of_face(topo: Topo, i: int) -> Optional["Prim"]:
    if not topo.planar[i]:
        return None
    return Prim(PLANE, face_plane_point(topo.faces[i]), topo.norms[i].copy())


# =============================================================================
# 8. SEGMENTATION: regions of facets lying on ONE curved primitive
# =============================================================================
#
# Model-guided seed growth (as in the previous engine), but with a
# "refine -> regrow" cycle and a final purge at the TIGHT tolerance: a
# region only enters play if EVERY one of its vertices lies on the surface
# within the mesh's noise. Faces that don't hold up stay free, and stay
# free: nothing is forced.


@dataclass
class Region:
    prim: "Prim"
    faces: List[int]
    rms: float = 0.0
    max_res: float = 0.0
    sag: float = 0.0  # max facet-to-surface deviation (chord sag)
    repeated: bool = False  # the same primitive appears elsewhere in the part
    closed_u: bool = False  # closed 360 degrees around the axis
    concave: bool = False  # material on the OUTSIDE of the surface (a hole)
    t_lo: float = 0.0
    t_hi: float = 0.0
    coverage: float = 0.0  # fraction of the full angle covered
    status: str = ""  # conversion outcome
    note: str = ""

    def label(self) -> str:
        p = self.prim
        if p.kind == FREE:
            ff = p.free
            return f"{'FREE':<6} {f'{ff.poles.shape[0]}x{ff.poles.shape[1]} poles':<22} {'CONCAVE' if self.concave else 'CONVEX':<8} {'--':>5} {len(self.faces):>5} faces"
        if p.kind == SPHERE:
            dim = f"ø{2 * p.r0:.4f}"
        elif p.kind == TORUS:
            dim = f"R{p.r0:.4f} r{p.r1:.4f}"
        elif abs(p.slope) < 1e-9:
            dim = f"ø{2 * p.r0:.4f}"
        else:
            dim = f"ø{2 * (p.r0 + p.slope * self.t_lo):.4f}-{2 * (p.r0 + p.slope * self.t_hi):.4f}"
        kind = "CONCAVE" if self.concave else "CONVEX"
        cl = "360°" if self.closed_u else f"{self.coverage * 360:.0f}°"
        return f"{p.label():<6} {dim:<22} {kind:<8} {cl:>5} {len(self.faces):>5} faces"


def _ring(seed: int, adj, taken, rings: int = 2, cap: int = 80, norms=None, cos_smooth: float = None) -> List[int]:
    """Complete adjacency rings around the seed, only through smooth edges."""
    out, seen, frontier = [seed], {seed}, [seed]
    for _ in range(max(1, rings)):
        nxt = []
        for f in frontier:
            for nb in adj[f]:
                if nb in seen or taken[nb]:
                    continue
                if norms is not None and cos_smooth is not None:
                    na, nb_ = norms[f], norms[nb]
                    if np.dot(na, na) > 0.5 and np.dot(nb_, nb_) > 0.5:
                        if abs(float(na @ nb_)) < cos_smooth:
                            continue
                seen.add(nb)
                out.append(nb)
                nxt.append(nb)
        frontier = nxt
        if not frontier or len(out) >= cap:
            break
    return out


def grow_region(seed_faces, prim, adj, taken, verts, norms, areas, tol: float, cos_ang: float, max_faces: int = 200000, refit: bool = True):
    region = set(seed_faces)
    frontier = list(seed_faces)
    next_refit = max(12, 2 * len(seed_faces)) if refit else 10**9
    while frontier:
        f = frontier.pop()
        for nb in adj[f]:
            if nb in region or taken[nb]:
                continue
            V = verts[nb]
            if V.size == 0:
                continue
            if float(np.abs(prim.dist(V)).max()) > tol:
                continue
            nrm = norms[nb]
            if np.dot(nrm, nrm) < 0.5:
                continue
            pn = prim.normal_at(V.mean(axis=0)[None, :])[0]
            if abs(float(pn @ nrm)) < cos_ang:
                continue
            region.add(nb)
            frontier.append(nb)
            if len(region) >= max_faces:
                frontier = []
                break
        if len(region) >= next_refit:
            newp = _refit(region, verts, norms, areas, prim.kind)
            if newp is not None:
                prim = newp
            next_refit = int(next_refit * 2)
    return sorted(region), prim


def largest_component(faces_idx, adj) -> List[int]:
    fs = set(faces_idx)
    best: List[int] = []
    while fs:
        s = fs.pop()
        comp, stack = {s}, [s]
        while stack:
            f = stack.pop()
            for nb in adj[f]:
                if nb in fs:
                    fs.discard(nb)
                    comp.add(nb)
                    stack.append(nb)
        if len(comp) > len(best):
            best = sorted(comp)
    return best


def sag_points(topo: Topo, faces: List[int]) -> np.ndarray:
    """INTERIOR points of the facets: centroids and side midpoints."""
    pts = []
    for i in faces:
        V = topo.verts[i]
        if len(V) < 3:
            continue
        pts.append(V.mean(axis=0))
        pts.append(0.5 * (V[0] + V[1]))
        pts.append(0.5 * (V[1] + V[2]))
        pts.append(0.5 * (V[0] + V[2]))
        if len(V) > 3:
            pts.append(0.5 * (V[-1] + V[0]))
    return np.array(pts) if pts else np.zeros((0, 3))


def region_sag(prim: "Prim", topo: Topo, faces: List[int]) -> float:
    """Sag: max distance between the facets' INTERIOR points and the surface."""
    pts = []
    for i in faces:
        V = topo.verts[i]
        if len(V) < 3:
            continue
        pts.append(V.mean(axis=0))
        pts.append(0.5 * (V[0] + V[1]))
        pts.append(0.5 * (V[1] + V[2]))
        pts.append(0.5 * (V[0] + V[2]))
        if len(V) > 3:
            pts.append(0.5 * (V[-1] + V[0]))
    if not pts:
        return 0.0
    return float(np.abs(prim.dist(np.array(pts))).max())


def describe_region(prim: "Prim", topo: Topo, faces: List[int]) -> Region:
    P = np.vstack([topo.verts[i] for i in faces])
    d = prim.dist(P)
    R = Region(prim=prim, faces=list(faces), rms=float(np.sqrt(np.mean(d**2))), max_res=float(np.abs(d).max()))
    R.sag = region_sag(prim, topo, faces)
    C = np.array([topo.cents[i] for i in faces])
    N = np.array([topo.norms[i] for i in faces])
    W = np.array([max(topo.areas[i], 1e-12) for i in faces])
    nat = prim.normal_at(C)
    R.concave = float(np.sum(W * np.einsum("ij,ij->i", nat, N))) < 0.0
    if prim.kind == FREE:
        R.coverage = 0.0
        R.closed_u = False
        return R
    if prim.kind in (AXIAL, TORUS):
        dd = P - prim.center
        t = dd @ prim.axis
        R.t_lo, R.t_hi = float(t.min()), float(t.max())
        u, v = ortho_frame(prim.axis)
        ang = np.sort(np.unique(np.round(np.arctan2(dd @ v, dd @ u), 6)))
        if len(ang) >= 3:
            gaps = np.diff(np.concatenate([ang, [ang[0] + 2 * math.pi]]))
            med = float(np.median(gaps))
            gmax = float(gaps.max())
            R.coverage = float(min(1.0, (2 * math.pi - gmax) / (2 * math.pi)))
            R.closed_u = gmax <= max(4.0 * med, math.radians(2.0)) and R.coverage > 0.9
        else:
            R.coverage = 0.0
    elif prim.kind == SPHERE:
        dd = C - prim.center
        dd = dd / np.maximum(np.linalg.norm(dd, axis=1), 1e-12)[:, None]
        m = (dd * W[:, None]).sum(axis=0) / W.sum()
        R.coverage = float(1.0 - np.linalg.norm(m))  # 0 = small cap, 1 = full sphere
        R.closed_u = False
    return R


_COS_QUASI = math.cos(math.radians(1.0))


def prims_equal(pa: "Prim", pb: "Prim", tol_len: float) -> bool:
    return prims_same(pa, pb, tol_len, math.cos(math.radians(0.05)))


_KIND_RANK = {AXIAL: 0, SPHERE: 2, TORUS: 3}


def _rank(p: "Prim") -> int:
    if p.kind == AXIAL:
        return 0 if abs(p.slope) < 1e-9 else 1
    return _KIND_RANK.get(p.kind, 9)


def prim_radius(p: "Prim") -> float:
    if p.kind in (PLANE, FREE):
        return 0.0
    if p.kind == TORUS:
        return float(p.r0 + p.r1)
    return float(abs(p.r0))


class Segmenter:
    """
    Regions of facets sitting EXACTLY on cylinders/cones/spheres/tori.
    only_cyl: only searches for cylinders (Phase B).
    """

    # ⚠️ seed_smooth = 20 degrees: the seed's neighborhood only crosses
    # nearly-smooth edges. At 50 degrees, the 45-degree CHAMFER facets next
    # to a cylinder strip got pulled in, the seed's fit was garbage, and
    # the part's rounded ends never got recognized.
    #
    # ⚠️ BUT 20 DEGREES IS A TESSELLATION THRESHOLD, NOT A GEOMETRIC ONE.
    # Between two facets of the SAME fillet, the angle is just the step it
    # was tessellated at: fine on good meshes (5-15 degrees), but a fillet
    # split into 8 strips over 180 degrees has 22.5-degree steps and stays
    # ENTIRELY outside: no neighbor, no candidate, the surface isn't even
    # visible. That's why the threshold is a ladder: try it tight first,
    # and only if the seed produces NO candidate at all does it reopen up
    # to seed_smooth_wide. The ceiling stays under the 45 degrees of
    # chamfers, which is the case the tight threshold protects us from.
    def __init__(
        self,
        topo,
        tol_fit: float,
        tol_grow: float,
        diag: float,
        min_faces: int = 4,
        grow_angle: float = 35.0,
        seed_smooth: float = 20.0,
        seed_smooth_wide: float = 32.0,
        allow_sphere: bool = True,
        allow_cone: bool = True,
        allow_torus: bool = True,
        allow_free: bool = True,
        blend_angle: float = 20.0,
        blend_min: int = 3,
        only_cyl: bool = False,
        threads: int = 1,
    ):
        self.topo = topo
        self.threads = max(1, int(threads))
        # exact constructor parameters: needed to rebuild the same
        # Segmenter inside the worker processes
        self._kw = dict(
            tol_fit=tol_fit,
            tol_grow=tol_grow,
            diag=diag,
            min_faces=min_faces,
            grow_angle=grow_angle,
            seed_smooth=seed_smooth,
            seed_smooth_wide=seed_smooth_wide,
            allow_sphere=allow_sphere,
            allow_cone=allow_cone,
            allow_torus=allow_torus,
            allow_free=allow_free,
            blend_angle=blend_angle,
            blend_min=blend_min,
            only_cyl=only_cyl,
            threads=1,
        )
        self.tol_fit, self.tol_grow, self.diag = tol_fit, tol_grow, diag
        self.min_faces = min_faces
        self.cos_ang = math.cos(math.radians(grow_angle))
        self.cos_seed = math.cos(math.radians(seed_smooth))
        self.cos_seed_wide = math.cos(math.radians(max(seed_smooth, seed_smooth_wide)))
        self.allow_sphere, self.allow_cone, self.allow_torus = allow_sphere, allow_cone, allow_torus
        self.allow_free = allow_free and not only_cyl
        self.blend_angle, self.blend_min = blend_angle, int(blend_min)
        self.only_cyl = only_cyl
        self.max_radius = 1.5 * diag
        self.taken = [False] * topo.nF

    # --- criteria ------------------------------------------------------------------
    def face_tol(self, p, i) -> float:
        """
        How far a facet can be from the surface, GIVEN HOW COARSE IT IS.

        ⚠️ A FACET KNOWS NO MORE THAN ITS CHORD LETS IT. A facet cutting an
        R1.5 fillet in 22-degree steps is already two or three hundredths
        below the true surface: demanding that its VERTICES land within a
        micron of that surface is asking the mesh for a precision it
        doesn't have. On test7 that's exactly what was happening to the
        compartments' vertical edges: the cone with the draft angle
        (r 1.585 at the bottom, 1.495 at the top) left a residual of a
        hundredth - a THIRD of the facets' own sag - and the region kept
        getting emptied out, leaving seven flat strips in place of one
        conical face.
        So the threshold is the wider of the absolute one and HALF the
        facet's sag, and it never exceeds the growth tolerance, which
        stays the definition of "close" for the rest of the code.
        """
        V = self.topo.verts[i]
        if V.shape[0] < 3:
            return self.tol_fit
        try:
            mids = np.vstack(((V + np.roll(V, -1, axis=0)) * 0.5, self.topo.cents[i][None, :]))
            sag = float(np.abs(p.dist(mids)).max())
        except Exception:
            return self.tol_fit
        if not np.isfinite(sag):
            return self.tol_fit
        return max(self.tol_fit, min(0.5 * sag, self.tol_grow))

    def within(self, p, i) -> bool:
        V = self.topo.verts[i]
        return bool(V.size) and float(np.abs(p.dist(V)).max()) <= self.face_tol(p, i)

    def _excess(self, p, reg):
        """residual of each facet MEASURED IN ITS OWN threshold (1.0 = at the limit)."""
        topo = self.topo
        out = []
        for i in reg:
            r_ = float(np.abs(p.dist(topo.verts[i])).max())
            out.append(r_ / max(self.face_tol(p, i), 1e-12))
        return np.array(out)

    def kind_ok(self, p) -> bool:
        if p is None or p.kind == PLANE:
            return False
        if self.only_cyl:
            return p.kind == AXIAL and abs(p.slope) < 1e-9
        return True

    def is_flat(self, reg) -> bool:
        """
        Region indistinguishable from a plane: not a curved feature.
        The curvature signal (deviation from the plane) has to clearly
        exceed the tolerance, otherwise it's noise that a huge cylinder
        "explains" by chance.
        """
        P = np.vstack([self.topo.verts[i] for i in reg])
        pl = fit_plane(P)
        return pl is not None and float(np.abs(pl.dist(P)).max()) <= 3.0 * self.tol_fit

    # --- type choice -----------------------------------------------------------------
    def choose_prim(self, faces: List[int], max_ndev: float = 15.0) -> Optional["Prim"]:
        """
        Reopens the question of the primitive's TYPE on a grown region,
        trying THE SIMPLEST FIRST and stopping at the first one that holds
        up (cylinder < cone < sphere < torus).

        ⚠️ A hole wall tessellated in strips has vertices ONLY on the two
        boundary circles: those points also lie exactly on a sphere. The
        residual can't tell them apart, the NORMALS can: the sphere would
        want them tilted by tens of degrees. That's why every candidate
        also has to pass the normals check.
        """
        topo = self.topo
        idx = list(faces)
        if len(idx) < 3:
            return None
        P = np.vstack([topo.verts[i] for i in idx])
        Nrep = np.vstack([np.tile(topo.norms[i], (len(topo.verts[i]), 1)) for i in idx])
        fallback = [None]

        def judge(p):
            if p is None or prim_radius(p) > self.max_radius:
                return False
            try:
                mx = float(np.abs(p.dist(P)).max())
                nd = normal_deviation(p, P, Nrep)
            except Exception:
                return False
            if not (np.isfinite(mx) and np.isfinite(nd)) or nd > max_ndev:
                return False
            if mx <= self.tol_fit:
                return True
            if fallback[0] is None or mx < fallback[0][0]:
                fallback[0] = (mx, p)
            return False

        ax = refit_exact(idx, topo.verts, topo.norms, topo.areas, AXIAL)
        if ax is not None:
            cyl = ax
            if abs(ax.slope) > 1e-9:
                cyl = lm_refine(Prim(AXIAL, ax.center.copy(), ax.axis.copy(), ax.r0, 0.0), P)
            if judge(cyl):
                return cyl
            if self.allow_cone and not self.only_cyl and abs(ax.slope) > 1e-9 and judge(ax):
                return ax
        if not self.only_cyl:
            if self.allow_sphere:
                sp = refit_exact(idx, topo.verts, topo.norms, topo.areas, SPHERE)
                if judge(sp):
                    return sp
            if self.allow_torus:
                tr = refit_exact(idx, topo.verts, topo.norms, topo.areas, TORUS)
                if judge(tr):
                    return tr
        return fallback[0][1] if fallback[0] else None

    # --- seeds ---------------------------------------------------------------------
    def seed_candidates(self, seed: List[int]):
        """
        Candidate primitives from the seed's neighborhood, ROBUST to
        unrelated neighbors.

        ⚠️ The neighborhood crosses smooth edges, so on a fillet (tangent
        by definition to its neighbors) it also picks up facets from OTHER
        surfaces, maybe even the big planar face next to it: the fit on
        that mixture is garbage and the seed dies. Here the large faces are
        removed, a fit is made, the facets within the growth tolerance are
        kept and it's refitted on the "inliers" alone.
        """
        topo = self.topo
        med = float(np.median([topo.areas[i] for i in seed]))
        core = [i for i in seed if topo.nverts[i] <= 8 and topo.areas[i] <= 12.0 * med]
        if len(core) < 4:
            return []
        vs = [topo.verts[i] for i in core]
        P = np.vstack(vs)
        Nrep = np.vstack([np.tile(topo.norms[i], (len(v), 1)) for i, v in zip(core, vs)])
        N = np.array([topo.norms[i] for i in core])
        W = np.array([max(topo.areas[i], 1e-9) for i in core])
        pl = fit_plane(P)
        if pl is not None and float(np.abs(pl.dist(P)).max()) <= 3.0 * self.tol_fit:
            return []
        out = []
        raw = rank_primitives(
            P, Nrep, N, W, self.allow_sphere and not self.only_cyl, self.allow_cone and not self.only_cyl, max_ndev=60.0, allow_torus=self.allow_torus and not self.only_cyl
        )
        raw = [c for c in raw if c.kind != PLANE]
        # ⚠️ algebraic fits on 4-10 facets can be off by 40% on the radius:
        # without an LM refinement, the true sphere never gets within
        # tolerance and the seed dies. Here every candidate (plus a sphere
        # and a cylinder tried regardless) gets refined on the inliers
        # alone, three times.
        extra = []
        if not self.only_cyl and self.allow_sphere and not any(c.kind == SPHERE for c in raw):
            sp0 = fit_sphere(P)
            if sp0 is not None:
                extra.append(sp0)
        if not any(c.kind == AXIAL for c in raw):
            ax0 = fit_axial(P, N, W)
            if ax0 is not None:
                if self.only_cyl:
                    ax0.slope = 0.0
                extra.append(ax0)
        for c in raw + extra:
            if self.only_cyl:
                if c.kind != AXIAL:
                    continue
                c.slope = 0.0
            try:
                c = lm_refine(c, P, iters=12)
            except Exception:
                continue
            inl = core
            for _ in range(3):
                inl2 = [i for i in core if float(np.abs(c.dist(topo.verts[i])).max()) <= self.tol_grow]
                if len(inl2) < 3:
                    c = None
                    break
                if set(inl2) == set(inl) and _ > 0:
                    break
                inl = inl2
                c2 = refit_exact(inl, topo.verts, topo.norms, topo.areas, c.kind)
                if c2 is None:
                    break
                if self.only_cyl:
                    c2.slope = 0.0
                c = c2
            if c is not None and prim_radius(c) <= self.max_radius:
                out.append((c, inl))
        out.sort(key=lambda t: (-len(t[1]), _rank(t[0])))
        return out[:3]

    # --- consolidation -------------------------------------------------------------
    def settle(self, reg, prim):
        """refine -> regrow -> tight purge, until the region is stable."""
        topo = self.topo
        reg = list(reg)
        p = prim
        last_n = 0
        if self.is_flat(reg):
            return None, "flat"
        for _ in range(6):
            if len(reg) > 1.25 * last_n or last_n == 0:
                p2 = self.choose_prim(reg)
                last_n = len(reg)
            else:
                p2 = refit_exact(reg, topo.verts, topo.norms, topo.areas, p.kind)
                if p2 is not None and prim_radius(p2) > self.max_radius:
                    p2 = None
            if not self.kind_ok(p2) and len(reg) >= 6:
                # ⚠️ region contaminated by facets from neighboring surfaces
                # (typical with the growth tolerance on long strips): keep
                # the half the SEED's fit explains better and try again
                res = self._excess(p, reg)
                cut = max(1.0, float(np.median(res)))
                trimmed = largest_component([i for i, r_ in zip(reg, res) if r_ <= cut], topo.adj)
                if len(trimmed) >= 4 and len(trimmed) < len(reg):
                    reg = trimmed
                    p2 = self.choose_prim(reg)
            if not self.kind_ok(p2):
                return None, "no primitive holds up"
            p = p2
            # ⚠️ CONTAMINATED REGION: growth from the wrong seed also
            # picked up facets from neighboring surfaces and the fit is a
            # compromise. The better half is kept, refitted, and restarted
            # from there.
            res = self._excess(p, reg)
            if res.max() > 1.0 and len(reg) >= 6:
                cut = max(1.0, float(np.median(res)))
                trimmed = largest_component([i for i, r_ in zip(reg, res) if r_ <= cut], topo.adj)
                if len(trimmed) >= 3:
                    p3 = refit_exact(trimmed, topo.verts, topo.norms, topo.areas, p.kind)
                    if p3 is not None and prim_radius(p3) <= self.max_radius:
                        p = p3
                        reg = trimmed
            reg2, _ = grow_region(reg, p, topo.adj, self.taken, topo.verts, topo.norms, topo.areas, self.tol_grow, self.cos_ang, refit=False)
            keep = largest_component([i for i in reg2 if self.within(p, i)], topo.adj)
            if len(keep) < self.min_faces:
                return None, f"too small ({len(keep)} faces within tolerance)"
            if set(keep) == set(reg):
                break
            reg = keep
        p = refit_exact(reg, topo.verts, topo.norms, topo.areas, p.kind)
        if not self.kind_ok(p) or prim_radius(p) > self.max_radius:
            return None, "final fit invalid"
        if abs(p.slope) < 1e-9:
            p.slope = 0.0
        if not all(self.within(p, i) for i in reg):
            keep = largest_component([i for i in reg if self.within(p, i)], topo.adj)
            if len(keep) < self.min_faces:
                return None, "too small after the final fit"
            reg = keep
        if self.is_flat(reg):
            return None, "flat"
        return reg, p

    def try_seed(self, s: int):
        """(prim, faces) from seed face s, or (None, reason)."""
        topo = self.topo
        # ⚠️ ADAPTIVE NEIGHBORHOOD. On a very fine mesh (0.03 mm facets on a
        # 5 mm radius) two rings of neighbors are flat within tolerance:
        # the curvature isn't visible and the seed dies ("no candidate").
        # The neighborhood widens until the curvature emerges or it runs out.
        # ⚠️ ONE CANDIDATE ISN'T ENOUGH: IT NEEDS SUPPORT. A primitive fitted
        # on three strips spanning twenty degrees of arc has a radius wrong
        # by a few hundredths, and growth stops right away: the fillet comes
        # out split into two or three faces instead of whole (this is
        # exactly what was happening to test7's vertical edges). So, as
        # long as the best support stays under MIN_SUPPORT facets, the
        # search keeps widening - first the rings, then the smoothing
        # threshold - and in the end the candidate with the most support is
        # tried first.
        MIN_SUPPORT = 6
        cands = []
        for cs in (self.cos_seed, self.cos_seed_wide):
            last_n = 0
            for rings, cap in ((1, 40), (2, 80), (4, 250)):
                seed = _ring(s, topo.adj, self.taken, rings, cap=cap, norms=topo.norms, cos_smooth=cs)
                seed = [i for i in seed if topo.planar[i] and topo.verts[i].size]
                if len(seed) < 4 or len(seed) == last_n:
                    if len(seed) == last_n:
                        break
                    continue
                last_n = len(seed)
                cands += self.seed_candidates(seed)
                if max((len(inl) for _, inl in cands), default=0) >= MIN_SUPPORT:
                    break
            if max((len(inl) for _, inl in cands), default=0) >= MIN_SUPPORT or cs == self.cos_seed_wide:
                break
        if not cands:
            return None, "no candidate on the seed"
        cands.sort(key=lambda t: (-len(t[1]), _rank(t[0])))
        why = ""
        for prim, inl in cands:
            reg, p2 = grow_region(inl, prim, topo.adj, self.taken, topo.verts, topo.norms, topo.areas, self.tol_grow, self.cos_ang)
            if len(reg) < self.min_faces:
                why = "insufficient growth"
                continue
            keep, pk = self.settle(reg, p2)
            if keep is None:
                why = pk
                continue
            return (pk, keep), ""
        return None, why

    # --- main loop -----------------------------------------------------------------
    # --- seed loop, sequential ------------------------------------------------------
    def _seed_loop(self, order, tries) -> Tuple[list, int]:
        topo = self.topo
        found: List[Tuple["Prim", List[int]]] = []
        tried = 0
        for s in order:
            if self.taken[s] or not topo.planar[s] or tries[s] >= 2:
                continue
            tried += 1
            res, why = self.try_seed(s)
            if res is None:
                tries[s] += 2
                continue
            pk, keep = res
            for i in keep:
                self.taken[i] = True
            found.append((pk, keep))
        return found, tried

    # --- seed loop, on multiple processes --------------------------------------------
    def _seed_loop_parallel(self, order, tries, pool) -> Tuple[list, int]:
        """
        Same search, but the seeds are split across the processes: each one
        runs its own greedy pass on its own slice (marking the facets it
        takes, so it doesn't pass over the same region twenty times), and
        then the main process accepts the regions ONE AT A TIME, largest
        first.

        ⚠️ Seed growth reads which facets are already taken: in parallel,
        each process works on a SNAPSHOT of that list, so two processes can
        claim the same facets. The check is here: a region that touches
        facets already assigned is thrown away and its seed goes back in
        the queue for the next round, with the updated list. In the end,
        the remaining seeds are done sequentially, so the result doesn't
        depend on the number of processes any more than it depends on seed
        order.
        """
        import pickle
        import tempfile

        topo = self.topo
        found: List[Tuple["Prim", List[int]]] = []
        tried = 0
        fd, path = tempfile.mkstemp(suffix=".refit.pkl")
        with os.fdopen(fd, "wb") as fh:
            pickle.dump((TopoLite(topo), self._kw), fh, protocol=pickle.HIGHEST_PROTOCOL)
        key = os.path.basename(path)
        pending = [s for s in order if topo.planar[s]]
        try:
            for _round in range(4):
                pending = [s for s in pending if not self.taken[s] and tries[s] < 2]
                if len(pending) < 24:
                    break
                taken_b = bytes(1 if t else 0 for t in self.taken)
                # interleaved slices: each process sees seeds scattered
                # across the whole part, so two processes rarely chase the
                # same region, and the load stays balanced
                n = min(_POOL_N or self.threads, max(1, len(pending) // 8))
                tasks = [(key, path, taken_b, pending[i::n]) for i in range(n)]
                try:
                    res = []
                    for part in pool.map(_worker_seeds, tasks):
                        res.extend(part)
                except Exception as ex:
                    Log.warn(f"worker processes unavailable ({type(ex).__name__}: {ex}): continuing on a single core")
                    break
                # larger regions first: in case of overlap, the one the
                # sequential loop would have found first wins
                res.sort(key=lambda t: (t[1] is None, -sum(topo.areas[i] for i in t[2]) if t[1] is not None else 0.0))
                for s, prim, keep in res:
                    tried += 1
                    if prim is None:
                        tries[s] += 2
                        continue
                    if any(self.taken[i] for i in keep):
                        continue  # conflict: retried later
                    for i in keep:
                        self.taken[i] = True
                    found.append((prim, keep))
        finally:
            try:
                os.unlink(path)
            except OSError:
                pass
        f2, t2 = self._seed_loop([s for s in pending if not self.taken[s]], tries)
        return found + f2, tried + t2

    @staticmethod
    def _pinch_count(topo, rf) -> int:
        """
        Vertices where the region's boundary PINCHES: more than two
        boundary edges at the same point. It's an inner ring touching the
        outer one.
        """
        rset = set(rf)
        cnt: Dict[int, int] = defaultdict(int)
        for i in rf:
            for k in topo.f_edges[i]:
                fs = topo.e_faces[k]
                if len(fs) == 2 and sum(1 for f in fs if f in rset) == 1:
                    for j in topo.e_verts[k]:
                        cnt[j] += 1
        return sum(1 for c in cnt.values() if c > 2)

    def _worth_merging(self, pa: "Prim", pb: "Prim", tol_len: float) -> bool:
        """
        Is it worth TRYING to merge these two? The final word is given by
        the union's refit: this is only the filter that avoids trying them
        all.
        ⚠️ prims_equal ALONE IS TOO STRICT FOR TWO PIECES OF THE SAME WALL.
        It asks for axes within 0.05 degrees, but a cylinder fitted on ten
        facets spanning 20 degrees has its axis uncertain by a few tenths
        of a degree: the two pieces of test4's R=3 cylinder come out 0.14
        and 0.28 degrees apart from each other and wouldn't merge. Here it
        widens to one degree and 2% of the radius - enough for the fit's
        uncertainty, too little to merge two different walls (in free-form
        fillet zones the neighboring "cylinders" sit at 2-13 degrees).
        """
        if prims_equal(pa, pb, tol_len):
            return True
        if pa is None or pb is None or pa.kind != pb.kind:
            return False
        if pa.kind in (PLANE, FREE):
            return False
        r = max(prim_radius(pa), prim_radius(pb), 1e-9)
        lim = max(tol_len, 0.02 * r)
        if pa.kind == SPHERE:
            return float(np.linalg.norm(pa.center - pb.center)) <= lim and abs(pa.r0 - pb.r0) <= lim
        if pa.axis is None or pb.axis is None:
            return False
        cos_ax = abs(float(pa.axis @ pb.axis))
        if cos_ax < _COS_QUASI:
            return False
        if pa.kind == TORUS:
            if abs(pa.r0 - pb.r0) > lim or abs(pa.r1 - pb.r1) > lim:
                return False
        else:
            sgn = 1.0 if float(pa.axis @ pb.axis) > 0 else -1.0
            if abs(pa.slope - sgn * pb.slope) > 0.02:
                return False
            if abs(prim_radius(pa) - prim_radius(pb)) > lim:
                return False
        d = pb.center - pa.center
        return float(np.linalg.norm(d - float(d @ pa.axis) * pa.axis)) <= lim

    def _merge_cosurface(self, found, tag: str = ""):
        """
        NEIGHBORING regions lying on the SAME surface: made into one.
        ⚠️ WITHOUT THIS THE MODEL FALLS APART: growth starts from different
        seeds and a long wall often ends up as two regions with the exact
        same primitive, which in the file become two faces separated by an
        edge that doesn't exist in the CAD.
        The union is refitted and accepted only if the single surface
        explains ALL the vertices within tol_fit: if it doesn't, they truly
        were two.
        """
        topo = self.topo
        tol_len = max(50.0 * self.tol_fit, 1e-5 * self.diag)
        n0 = len(found)
        changed = True
        while changed and len(found) > 1:
            changed = False
            owner = {}
            for k, (_, rf) in enumerate(found):
                for i in rf:
                    owner[i] = k
            for k in range(len(found)):
                pa, ra = found[k]
                if pa is None:
                    continue
                nbrs = {owner[j] for i in ra for j in topo.adj[i] if j in owner and owner[j] != k}
                for q in sorted(nbrs):
                    pb, rb = found[q]
                    if pb is None or not self._worth_merging(pa, pb, tol_len):
                        continue
                    union = sorted(set(ra) | set(rb))
                    p = refit_exact(union, topo.verts, topo.norms, topo.areas, pa.kind)
                    if p is None:
                        continue
                    if abs(p.slope) < 1e-9:
                        p.slope = 0.0
                    p = self._snap_cylinder(p, union)
                    # ⚠️ THE YARDSTICK IS tol_grow, NOT tol_fit. tol_fit is
                    # the threshold for saying "this facet is on this
                    # surface": the distance between two pieces of the SAME
                    # wall fitted separately is a different matter. A
                    # twenty-degree arc has a poorly conditioned fit -
                    # radius and center trade off against each other - and
                    # its vertices sit a micron from the cylinder that
                    # explains everything else: with tol_fit the merge
                    # almost never fired (zero times on test4) and the R=3
                    # cylinder stayed split into three faces with three
                    # slivers in between, where the CAD has ONE face.
                    # tol_grow is already the definition of "close" the
                    # regions grew by, and it's tight enough: a real step
                    # between two bores (2 hundredths) blows through it
                    # sixfold and stays two faces.
                    P = np.vstack([topo.verts[i] for i in union])
                    if float(np.abs(p.dist(P)).max()) > self.tol_grow:
                        continue
                    # ⚠️ AND IT MUST NOT PINCH. Merging four little sphere
                    # fragments of 0.015 mm2 used to produce a face whose
                    # inner ring touches the outer one at a vertex:
                    # BRepCheck accepts it in memory, but after the STEP
                    # write-and-reread round trip it becomes
                    # "IntersectingWires" and the part is no longer valid
                    # (measured on test8). If the union pinches where the
                    # pieces didn't, they're left separate.
                    if self._pinch_count(topo, union) > max(
                        self._pinch_count(topo, ra), self._pinch_count(topo, rb)
                    ):
                        continue
                    found[k] = (p, union)
                    found[q] = (None, [])
                    changed = True
                    break
                if changed:
                    break
            found = [(p, rf) for p, rf in found if p is not None]
        if tag and len(found) < n0:
            Log.debug(f"Co-surface merge ({tag}): {n0} -> {len(found)} regions")
        return found

    def run(self) -> List[Region]:
        topo = self.topo
        nF = topo.nF
        order = sorted(range(nF), key=lambda i: -topo.areas[i])
        tries = [0] * nF
        t_seed = time.perf_counter()
        pool = _seed_pool(self.threads, want=nF // 200) if nF >= 400 else None
        if pool is not None:
            found, tried = self._seed_loop_parallel(order, tries, pool)
        else:
            found, tried = self._seed_loop(order, tries)
        Log.debug(f"Seeds: {tried:,} tried in {time.perf_counter() - t_seed:.2f}s ({'1 process' if pool is None else f'{_POOL_N} processes'})")

        # ⚠️ STRAIGHTEN BEFORE MERGING. Two seeds on the same hole wall come
        # out as two CONES with semi-angle 0.03 and 0.35 degrees: it's the
        # same surface, but prims_equal sees them as different and doesn't
        # merge them, and the file ends up with two faces where the CAD has
        # one (measured on test4: a 13.7 mm2 cylinder split into two cones
        # of 10.4 and 2.4). Snapping to a cylinder BEFORE the comparison
        # makes them merge on their own.
        found = [(self._snap_cylinder(pp, rf), rf) for pp, rf in found]
        found = self._merge_cosurface(found, "after the seeds")

        found = self.relabel(found)
        regions = [describe_region(p, topo, rf) for p, rf in found]
        # ⚠️ facets each covering more than 30 degrees aren't a tessellation
        # of that surface: it's a random fit on 4 faces
        regions = [R for R in regions if not (R.prim.kind in (AXIAL, TORUS) and R.coverage > 0 and R.coverage * 360.0 / max(len(R.faces), 1) > 30.0)]

        # ⚠️ fragments: a 5-facet "cylinder" covering 7 degrees, or a torus
        # with a major radius under 0.3 mm, isn't a real feature: it's a
        # random fit on noise, and converted it produces a curved chip
        # between smooth faces (the step artifacts). They stay tessellated.
        def _fragment(R):
            # ⚠️ a narrow but LONG band (60 facets over 10 degrees) is a
            # real fillet: the facet count matters as much as the coverage
            if R.prim.kind in (AXIAL, TORUS) and not R.closed_u and R.coverage * 360.0 < 15.0 and len(R.faces) < 12:
                return True
            if R.prim.kind == TORUS and (R.prim.r0 < 0.3 * R.prim.r1 or (len(R.faces) < 6 and R.coverage * 360.0 < 30.0)):
                return True
            if R.prim.kind == SPHERE and R.coverage < 0.005 and len(R.faces) < 10:
                return True
            return False

        regions = [R for R in regions if not _fragment(R)]
        # ⚠️ RECOVERY MUST BE DONE ON THE FINAL MAP. Regions thrown out by
        # the filters (fragments, random fits) keep their facets occupied
        # up to here: recovering before the filters would mean finding
        # them still claimed by a region that then disappears, leaving
        # them unclaimed.
        found = self._absorb_free([(R.prim, list(R.faces)) for R in regions])
        found = self._wedged(found)
        found = [(self._snap_cylinder(pp, rf), rf) for pp, rf in found]
        # ⚠️ AND IT MERGES AGAIN. The first pass looks at the regions AS
        # THEY COME OUT OF THE SEEDS: two pieces of the same wall separated
        # by a few unclaimed facets don't even touch, and pass unscathed.
        # Only after the unclaimed facets are recovered (_absorb_free) and
        # the second straightening to a cylinder do they become neighbors
        # and comparable. Measured on test4: the R=3 cylinder came out
        # split into two faces of 10.5 and 2.4 mm2 with three slivers in
        # between, and in the CAD it's ONE single face.
        found = self._merge_cosurface(found, "after recovery")
        regions = [describe_region(pp, topo, rf) for pp, rf in found]
        regions = self._blend_patches(regions)
        regions.sort(key=lambda R: -sum(topo.areas[i] for i in R.faces))
        Log.info(f"Seeds tried {tried:,} · curved regions found {len(regions)} · facets involved {sum(len(R.faces) for R in regions):,} / {nF:,}")
        return regions

    def _face_scale(self, i) -> float:
        """Face extent: median distance BETWEEN ALL its vertices.
        ⚠️ NOT the side lengths: a large planar face with a finely
        segmented outline (a plate's shaped edge) would have sides as
        short as a facet's and would look just as fine."""
        V = self.topo.verts[i]
        if V.shape[0] < 2:
            return 1e-12
        if V.shape[0] > 64:
            V = V[np.linspace(0, V.shape[0] - 1, 64).astype(int)]
        d = np.linalg.norm(V[:, None, :] - V[None, :, :], axis=2)
        d = d[d > 0]
        return float(np.median(d)) if d.size else 1e-12

    @staticmethod
    def _bridges(nF: int, adj):
        """
        BRIDGE edges of the facet graph (Tarjan, iterative).
        ⚠️ THESE ARE NEEDED BECAUSE THE TEXTURE SIGNAL IS ALREADY NEARLY
        PERFECT. Measured on test4 against the real CAD file: the rule
        "dihedral > 30 or texture > 1.5" finds 1,285 of the CAD's 1,323
        edges with 18 false positives out of 7,200. The ones that escape
        all have a SINGLE mesh edge: they're corner contacts, where two CAD
        faces touch almost at a point. Geometrically invisible (dihedral
        1-3 degrees, same texture), but topologically obvious: they're
        bridges, and without cutting them, twenty CAD faces ended up in a
        single section. Cutting them, the sections spanning more than one
        original face drop from 7 to 3.
        """
        disc = [-1] * nF
        low = [0] * nF
        out = set()
        timer = 0
        for s0 in range(nF):
            if disc[s0] >= 0:
                continue
            disc[s0] = low[s0] = timer
            timer += 1
            stack = [(s0, -1, iter(adj[s0]))]
            while stack:
                u, par, it = stack[-1]
                avanti = False
                for v in it:
                    if v == par:
                        continue
                    if disc[v] < 0:
                        disc[v] = low[v] = timer
                        timer += 1
                        stack.append((v, u, iter(adj[v])))
                        avanti = True
                        break
                    low[u] = min(low[u], disc[v])
                if not avanti:
                    stack.pop()
                    if stack:
                        q = stack[-1][0]
                        low[q] = min(low[q], low[u])
                        if low[u] > disc[q]:
                            out.add((min(q, u), max(q, u)))
        return out

    def sections(self, ang_cut: float = 30.0, size_cut: float = 1.5):
        """
        The original CAD's FACES, cut out from the mesh.

        ⚠️ THE TESSELLATOR WORKS ONE FACE AT A TIME, and the mesh remembers
        it: inside a face the texture is uniform, across a CAD edge it
        changes abruptly. Cuts are made where there's a sharp edge (large
        dihedral) or where the texture changes, and what's left are the
        starting part's faces. Measured on test4 against the real CAD
        file: the sections cut out EXACTLY the countersink's cone (26.99
        mm2), the cylinders (105.52, 60.31) and - what matters here - the
        free-form surfaces (5.21 and 5.21) that the quadric fit was
        covering with five spheres and a torus.
        """
        topo = self.topo
        nF = topo.nF
        N = np.array([topo.norms[i] for i in range(nF)])
        siz = np.array([self._face_scale(i) for i in range(nF)])
        cang = math.cos(math.radians(ang_cut))
        adj = [[] for _ in range(nF)]
        for i in range(nF):
            for j in topo.adj[i]:
                if float(N[i] @ N[j]) < cang:
                    continue
                if abs(math.log2(max(siz[i], 1e-12) / max(siz[j], 1e-12))) > size_cut:
                    continue
                adj[i].append(j)
        br = self._bridges(nF, adj)
        if br:
            adj = [[j for j in adj[i] if (min(i, j), max(i, j)) not in br] for i in range(nF)]
        lab = np.full(nF, -1)
        nc = 0
        for s0 in range(nF):
            if lab[s0] >= 0:
                continue
            lab[s0] = nc
            stack = [s0]
            while stack:
                x = stack.pop()
                for y in adj[x]:
                    if lab[y] < 0:
                        lab[y] = nc
                        stack.append(y)
            nc += 1
        return lab, nc

    def _snap_cylinder(self, p, rf):
        """
        ⚠️ THE CONE THAT'S REALLY A CYLINDER. The axial fit has one more
        degree of freedom than the cylinder and always uses it: on a hole
        wall it comes out as a cone with a 0.02-degree semi-angle, i.e. a
        cylinder with a half-micron taper over the whole face.
        Geometrically it's the same surface, but in the STEP file it's a
        CONICAL_SURFACE where the CAD had a CYLINDRICAL_SURFACE, and anyone
        who opens the model notices. The cylinder is tried: if it explains
        the facets as well as the cone, it wins. (on test4: six cones out
        of nine turn back into cylinders, and the three that remain are
        the part's three real ones, all at 45 degrees)
        """
        if p is None or p.kind != AXIAL or abs(p.slope) < 1e-12:
            return p
        try:
            if abs(p.slope) > 0.05:
                return p  # real taper (>2.9 degrees)
            P = np.vstack([self.topo.verts[i] for i in rf])
            t = (P - p.center) @ p.axis
            rad = P - p.center - np.outer(t, p.axis)
            c = Prim(AXIAL, p.center.copy(), p.axis.copy(), float(np.linalg.norm(rad, axis=1).mean()), 0.0)
            q = lm_refine(c, P) or c
            if abs(q.slope) > 1e-12:
                q = Prim(AXIAL, q.center, q.axis, q.r0, 0.0)
            if all(self.within(q, i) for i in rf):
                return q
            # ⚠️ OR: THE CYLINDER DOESN'T EXPLAIN IT WORSE THAN THE CONE.
            # After a merge, the union's refit comes out conical almost
            # every time (the axial fit has one more degree of freedom and
            # uses it): a 0.04-degree semi-angle is a cylinder, but its
            # vertices sit just past face_tol and the snap wouldn't fire,
            # leaving a CONICAL_SURFACE in the file where the CAD has a
            # CYLINDRICAL_SURFACE.
            d_cono = float(np.abs(p.dist(P)).max())
            d_cil = float(np.abs(q.dist(P)).max())
            if d_cil <= max(1.05 * d_cono, self.tol_fit):
                return q
        except Exception:
            pass
        return p

    def _blend_patches(self, regions):
        """
        ⚠️ THE CHAMFER RUNNING ALONG A CURVED EDGE. Where the edge to be
        chamfered is neither straight nor a circular arc, the blending
        surface is neither a cone nor a torus: it's a free-form surface,
        and the CAD writes it as a B-spline (also recognizable from the
        tessellation, which there is ten times denser). A quadric fit
        can't help but split it into dozens of ten-degree OSCULATING
        little cylinders, each with a different radius and a jagged
        outline: those are the "faces with too many sides" you see in the
        finished part.
        Here those fragments are recognized (a curved band covering a few
        degrees), grouped by contact, stitched together with the loose
        islands they pull in, and replaced by ONE free-form surface. The
        fit must stay within the SAME tolerance as the fragments: if it
        doesn't fit, nothing is touched.
        """
        if not self.allow_free or len(regions) < self.blend_min:
            return regions
        topo = self.topo
        owner = {}
        for k, R in enumerate(regions):
            for i in R.faces:
                owner[i] = k
        lab, nc = self.sections()
        n_in = np.bincount(lab, minlength=nc)
        drop, extra = set(), []
        for c in range(nc):
            m = np.where(lab == c)[0]
            if len(m) < 4 * self.min_faces:
                continue
            # regions that lie ALMOST ENTIRELY inside this section: only
            # those can be dissolved. A straddling region is left alone,
            # otherwise a hole would open up in it.
            cnt, tot = {}, {}
            for i in m:
                q = owner.get(i)
                if q is not None:
                    cnt[q] = cnt.get(q, 0) + 1
            for q in cnt:
                tot[q] = len(regions[q].faces)
            dentro = [q for q in cnt if cnt[q] >= 0.8 * tot[q]]
            if len(dentro) < self.blend_min:
                continue
            fuori = {i for i in m if owner.get(i) is not None and owner[i] not in dentro}
            faces = [int(i) for i in m if i not in fuori]
            if len(faces) < 4 * self.min_faces:
                continue
            faces = self._fill_islands(faces, owner)
            P = np.vstack([topo.verts[i] for i in faces])
            N = np.vstack([np.tile(topo.norms[i], (len(topo.verts[i]), 1)) for i in faces])
            # ⚠️ THE REAL TEST IS AGAINST THE MESH, NOT AGAINST THE
            # FRAGMENTS. Probe points taken on the fragments' primitives
            # seemed like a good idea, but on a 1 mm-radius patch covered
            # by ø6 spheres, the fragments are the ones getting it wrong,
            # and the B-spline was rejected while being CLOSER to the part
            # than they were. The sag - distance of the facets' INTERIOR
            # points from the surface - is measured the same way for
            # everyone and tells you how far it strays from the mesh
            # between one vertex and the next: the new surface only needs
            # to not be worse than the ones it replaces.
            sag_old = max(region_sag(regions[q].prim, topo, list(regions[q].faces)) for q in dentro)
            pr = fit_free(P, N, self.tol_fit, check=sag_points(topo, faces), check_tol=2.0 * sag_old + 2.0 * self.tol_fit)
            if pr is None:
                continue
            # ⚠️ THE REAL TEST IS AGAINST THE MESH, NOT AGAINST THE
            # FRAGMENTS. Probe points taken on the fragments' primitives
            # seemed like a good idea, but on a 1 mm-radius patch covered
            # by ø6 spheres the fragments are the ones getting it wrong:
            # the probe said 2e-03 while the B-spline was CLOSER to the
            # part than they were. The sag - distance of the facets'
            # INTERIOR points from the surface - is measured the same way
            # for everyone, and tells you how far it strays from the mesh
            # between one vertex and the next: the new surface must not be
            # worse than the ones it replaces.
            drop |= set(dentro)
            extra.append(describe_region(pr, topo, faces))
        if not extra:
            return regions
        Log.debug(f"Free-form patches: {len(extra)} sections (in place of {len(drop)} fragments)")
        return [R for k, R in enumerate(regions) if k not in drop] + extra

    def _fill_islands(self, faces, owner, max_frac: float = 0.25):
        """UNCLAIMED facets fully surrounded by the patch: they belong to
        it. Without this, the new face's contour would have interior rings
        of four triangles and OpenCascade rejects it (UnorientableShape)."""
        topo = self.topo
        gs = set(faces)
        seen, comps = set(), []
        for s0 in range(topo.nF):
            if s0 in gs or s0 in seen:
                continue
            stack = [s0]
            seen.add(s0)
            c = []
            while stack:
                x = stack.pop()
                c.append(x)
                for y in topo.adj[x]:
                    if y not in gs and y not in seen:
                        seen.add(y)
                        stack.append(y)
            comps.append(c)
        if len(comps) <= 1:
            return sorted(gs)
        comps.sort(key=len, reverse=True)
        cap = max(4, int(max_frac * len(gs)))
        for c in comps[1:]:
            if len(c) > cap:
                continue
            if any(i in owner for i in c):
                continue  # the island is a real region: leave it
            gs |= set(c)
        return sorted(gs)

    def _belongs(self, p, i, base) -> bool:
        """
        Facet i lies on surface p, or is at least an island inside base.

        ⚠️ VERTEX NOISE ISN'T THE SAME AS CHORD SAG. face_tol widens the
        threshold when the facet is COARSE, but a long, thin strip has
        almost no sag and stays at the bare micron: on test6 the ø4
        fillet's strips have vertices 1.3 microns outside the cylinder
        (tol_fit allows 0.63) and dozens stayed outside, one next to the
        other, like flat slivers in the middle of a cylindrical face. That
        1.3 micron is the noise the tessellator wrote the vertices with,
        not a different surface.
        For these, and ONLY during recovery, topology counts: if the facet
        is surrounded by the region (at least two neighbors inside), faces
        the same way within 5 degrees, and is within the GROWTH tolerance,
        it's theirs. Three conditions together: close, parallel and closed
        in.
        """
        V = self.topo.verts[i]
        if V.size == 0:
            return False
        d = float(np.abs(p.dist(V)).max())
        if d <= self.face_tol(p, i):
            return True
        if d > self.tol_grow or base is None:
            return False
        if sum(1 for j in self.topo.adj[i] if j in base) < 2:
            return False
        nd = normal_deviation(p, V, np.tile(self.topo.norms[i], (len(V), 1)))
        return bool(np.isfinite(nd) and nd <= 5.0)

    def _wedged(self, found):
        """
        Islands of unclaimed facets fully surrounded by rebuilt surfaces.

        ⚠️ ON A TANGENT JOINT THE MESH SITS ON NEITHER SIDE. Where a ø4
        fillet meets a ø8 cylinder it's tangent to, the tessellator puts a
        strip with one edge on the fillet and the other on the cylinder: it
        doesn't sit EXACTLY on either one (on test6 it's 1.4 microns off
        from both, and the tolerance allows 0.63). Neither region claims
        it, and it stays a long, thin flat sliver between two analytic
        faces - the defect you see on test6 and test2 near a hole, chamfer
        and fillet.
        Here, GROUPS of facets left unclaimed that don't touch any other
        unclaimed facet (so they're an island surrounded by already-rebuilt
        surfaces), small ones, are picked up and handed to whichever
        neighboring region explains them best - as long as it explains them
        within the growth tolerance and with normals within 5 degrees.
        Choosing "one or the other" doesn't matter beyond a micron; leaving
        them unclaimed does.
        """
        topo = self.topo
        owner = {}
        for k, (_, rf) in enumerate(found):
            for i in rf:
                owner[i] = k
        loose = [i for i in range(topo.nF) if i not in owner and topo.verts[i].size]
        if not loose:
            return found
        lset = set(loose)
        par = {i: i for i in loose}

        def find(a):
            while par[a] != a:
                par[a] = par[par[a]]
                a = par[a]
            return a

        for i in loose:
            for j in topo.adj[i]:
                if j in lset:
                    a, b = find(i), find(j)
                    if a != b:
                        par[a] = b
        comps = defaultdict(list)
        for i in loose:
            comps[find(i)].append(i)
        add = defaultdict(list)
        for g in comps.values():
            if len(g) > 12:
                continue
            nbrs = {owner[j] for i in g for j in topo.adj[i] if j in owner}
            if not nbrs:
                continue
            a_g = sum(topo.areas[i] for i in g)
            best, best_d = None, None
            for k in nbrs:
                pk, rf = found[k]
                if pk is None or a_g > 0.25 * sum(topo.areas[i] for i in rf):
                    continue
                dmax = 0.0
                ok = True
                for i in g:
                    V = topo.verts[i]
                    d = float(np.abs(pk.dist(V)).max())
                    nd = normal_deviation(pk, V, np.tile(topo.norms[i], (len(V), 1)))
                    if d > self.tol_grow or not np.isfinite(nd) or nd > 5.0:
                        ok = False
                        break
                    dmax = max(dmax, d)
                if ok and (best_d is None or dmax < best_d):
                    best, best_d = k, dmax
            if best is not None:
                add[best] += g
        if not add:
            return found
        n_add = 0
        for k, g in add.items():
            p0, rf = found[k]
            union = sorted(set(rf) | set(g))
            p2 = refit_exact(union, topo.verts, topo.norms, topo.areas, p0.kind)
            if p2 is None or prim_radius(p2) > self.max_radius:
                continue
            if any(not self.within(p2, i) for i in rf):
                continue
            if abs(p2.slope) < 1e-9:
                p2.slope = 0.0
            found[k] = (p2, union)
            n_add += len(g)
        if n_add:
            Log.debug(f"Wedged slivers recovered: +{n_add}")
        return found

    def _absorb_free(self, found):
        """
        Facets left unclaimed that lie on a region's surface.

        ⚠️ WHOEVER GROWS FIRST WINS, AND THAT ISN'T NECESSARILY THE BEST
        ONE. Seed order decides who claims what, and with multiple
        processes the seeds start from a SNAPSHOT of the facets already
        claimed: two seeds on the SAME surface grow in parallel, the
        smaller one gets rejected at commit and its facets stay orphaned.
        On test2 the ø1 fillet's band came out with 92 facets instead of
        117, and the 25 orphans turned into flat strips around the holes -
        the "hole + chamfer + fillet" defect. Here, with the map frozen,
        every region tries again to grow onto the facets STILL unclaimed:
        only additions, never removals, and the refitted primitive must
        keep explaining all the starting facets.
        """
        topo = self.topo
        taken = [False] * topo.nF
        for _, rf in found:
            for i in rf:
                taken[i] = True
        n_add = 0
        for k in range(len(found)):
            p, rf = found[k]
            if p is None or not rf:
                continue
            reg2, _ = grow_region(rf, p, topo.adj, taken, topo.verts, topo.norms, topo.areas, self.tol_grow, self.cos_ang, refit=False)
            base = set(rf)
            add = [i for i in reg2 if i not in base and not taken[i] and self._belongs(p, i, base)]
            if not add:
                continue
            union = sorted(base | set(add))
            p2 = refit_exact(union, topo.verts, topo.norms, topo.areas, p.kind)
            if p2 is None or prim_radius(p2) > self.max_radius:
                continue
            # ⚠️ the STARTING facets must stay explained by the refitted
            # primitive: recovery can only add, never make things worse.
            if any(not self.within(p2, i) for i in rf):
                continue
            keep = set(largest_component([i for i in union if self._belongs(p2, i, base)], topo.adj))
            if len(keep) <= len(rf) or not base <= keep:
                continue
            if abs(p2.slope) < 1e-9:
                p2.slope = 0.0
            n_add += len(keep) - len(rf)
            found[k] = (p2, sorted(keep))
            for i in keep:
                taken[i] = True
        if n_add:
            Log.debug(f"Unclaimed facets recovered: +{n_add}")
        return found

    def relabel(self, found):
        """
        Competitive reassignment of boundary facets.

        ⚠️ Two TANGENT surfaces (a corner sphere and a fillet) are
        indistinguishable near the joint: the sphere's row of thin little
        triangles is within tolerance of the cylinder too, and the first
        region to grow there claims it. The boundary between the two
        regions shifts by one row and the joint curve is no longer the
        exact circle. Here every boundary facet goes to the region whose
        primitive explains it BEST (lowest residual), then it's refitted.
        """
        topo = self.topo
        if len(found) < 2:
            return found
        prims = [p for p, _ in found]
        sets = [set(rf) for _, rf in found]
        owner = {}
        for k, s in enumerate(sets):
            for i in s:
                owner[i] = k
        moved_total = 0
        for _ in range(4):
            moved = 0
            for i in list(owner):
                a = owner[i]
                nb_regions = {owner[j] for j in topo.adj[i] if j in owner and owner[j] != a}
                if not nb_regions:
                    continue
                V = topo.verts[i]
                ra = float(np.abs(prims[a].dist(V)).max())
                best, rb = a, ra
                for b in nb_regions:
                    r_ = float(np.abs(prims[b].dist(V)).max())
                    if r_ <= self.tol_fit and r_ < 0.5 * rb:
                        best, rb = b, r_
                if best != a:
                    sets[a].discard(i)
                    sets[best].add(i)
                    owner[i] = best
                    moved += 1
            if not moved:
                break
            moved_total += moved
            for k in range(len(found)):
                if len(sets[k]) >= 3:
                    p2 = refit_exact(sorted(sets[k]), topo.verts, topo.norms, topo.areas, prims[k].kind)
                    if p2 is not None:
                        if abs(p2.slope) < 1e-9:
                            p2.slope = 0.0
                        prims[k] = p2
        if moved_total:
            Log.debug(f"Reassigned {moved_total} boundary facets between tangent regions")
        out = []
        for k in range(len(found)):
            rf = largest_component(sorted(sets[k]), topo.adj)
            if len(rf) >= self.min_faces:
                out.append((prims[k], rf))
        return out


# =============================================================================
# 8b-bis. PARALLEL SEEDS
#
# ⚠️ WHY PROCESSES AND NOT THREADS. The calls into OpenCascade and this
# script's numpy code hold the GIL (measured: eight volume computations on
# four threads cost as much as in sequence), so threads gain nothing. The
# only pieces that are truly parallelizable are the ones working purely on
# numpy arrays: the primitive search from the seeds, which is the biggest
# slice of Phases B and C. The in-solid replacements stay sequential: they
# modify a shared OCC structure and must be checked one at a time.
# =============================================================================


class TopoLite:
    """Only the topology fields the segmentation needs: numpy and lists,
    so it can be shipped to another process."""

    __slots__ = ("nF", "verts", "norms", "areas", "nverts", "planar", "adj")

    def __init__(self, topo):
        n = topo.nF
        self.nF = n
        self.verts = [topo.verts[i] for i in range(n)]
        self.norms = [topo.norms[i] for i in range(n)]
        self.areas = [topo.areas[i] for i in range(n)]
        self.nverts = [topo.nverts[i] for i in range(n)]
        self.planar = [topo.planar[i] for i in range(n)]
        self.adj = [sorted(topo.adj[i]) for i in range(n)]

    def __getstate__(self):
        return {k: getattr(self, k) for k in self.__slots__}

    def __setstate__(self, st):
        for k, v in st.items():
            setattr(self, k, v)


_POOL = None
_POOL_N = 0
_W_KEY = None
_W_SEG = None


def _seed_pool(threads: int, want: Optional[int] = None):
    """
    Process pool, created once and reused by every phase.

    ⚠️ HOW MANY PROCESSES. Not "every one the machine has": spinning up an
    interpreter that imports OpenCascade costs nearly a second, and on a
    part with half a second of seed work, twenty-two processes take longer
    to start than they save (measured: test3.stl 4.2 s with four
    processes, 6.1 s with twenty-two). So the count is scaled to the work
    at hand, with -j as the ceiling.
    """
    global _POOL, _POOL_N
    if threads is None or threads < 2:
        return None
    if _POOL is not None:
        return _POOL
    k = int(threads) if want is None else int(max(2, min(int(threads), int(want))))
    try:
        from concurrent.futures import ProcessPoolExecutor

        _POOL = ProcessPoolExecutor(max_workers=k)
        _POOL_N = k
    except Exception as ex:
        Log.warn(f"no multiprocessing available ({type(ex).__name__}: {ex}): using a single core")
        _POOL, _POOL_N = None, 0
    return _POOL


def _worker_boot():
    return os.getpid()


def warm_pool(threads: int, n_faces: int = 0) -> None:
    """
    Starts the worker processes BEFORE they're needed.

    ⚠️ Every worker process is a fresh interpreter that has to import
    OpenCascade: nearly a second. Starting them when needed means paying
    that second all at once, and on a small part it eats the gain. Started
    here, they boot up while Phase A works on a single core.
    """
    pool = _seed_pool(threads, want=n_faces // 5000)
    if pool is None:
        return
    try:
        for _ in range(_POOL_N):
            pool.submit(_worker_boot)
    except Exception:
        pass
    Log.debug(f"Worker processes started: {_POOL_N}")


def close_pool() -> None:
    global _POOL, _POOL_N
    if _POOL is not None:
        try:
            _POOL.shutdown(wait=False, cancel_futures=True)
        except Exception:
            pass
    _POOL, _POOL_N = None, 0


def _worker_seeds(task):
    """Tries a list of seeds against a snapshot of the facets already claimed."""
    import pickle

    global _W_KEY, _W_SEG
    key, path, taken_b, seeds = task
    if _W_KEY != key or _W_SEG is None:
        with open(path, "rb") as fh:
            lite, kw = pickle.load(fh)
        _W_SEG = Segmenter(lite, **kw)
        _W_KEY = key
    seg = _W_SEG
    seg.taken = [b != 0 for b in taken_b]
    out = []
    for s in seeds:
        if seg.taken[s]:
            continue
        try:
            res, why = seg.try_seed(s)
        except Exception as ex:
            out.append((s, None, [f"{type(ex).__name__}: {ex}"]))
            continue
        if res is None:
            out.append((s, None, []))
            continue
        prim, keep = res
        # ⚠️ the process marks what it took: without this it would redo the
        # same region starting from every facet that's part of it
        for i in keep:
            seg.taken[i] = True
        out.append((s, prim, keep))
    return out


def segment_curved(
    topo: Topo,
    tol_fit: float,
    tol_grow: float,
    diag: float,
    min_faces: int = 4,
    grow_angle: float = 35.0,
    seed_smooth: float = 20.0,
    seed_smooth_wide: float = 32.0,
    allow_sphere: bool = True,
    allow_cone: bool = True,
    allow_torus: bool = True,
    allow_free: bool = True,
    only_cyl: bool = False,
    threads: int = 1,
) -> List[Region]:
    return Segmenter(
        topo,
        tol_fit,
        tol_grow,
        diag,
        min_faces=min_faces,
        grow_angle=grow_angle,
        seed_smooth=seed_smooth,
        seed_smooth_wide=seed_smooth_wide,
        allow_sphere=allow_sphere,
        allow_cone=allow_cone,
        allow_torus=allow_torus,
        allow_free=allow_free,
        only_cyl=only_cyl,
        threads=threads,
    ).run()


# =============================================================================
# 9. ANALYTIC SURFACES: parametrization, pcurves, edges
# =============================================================================
#
# ⚠️ We compute the pcurves (curves in the surface's (u,v) space) OURSELVES,
# not OCC's projector. The projector reports u in [0, 2pi): a boundary
# crossing the seam (u = 0) gets split, the wire doesn't close in parameter
# space and the face comes out invalid or, worse, as the region's
# COMPLEMENT (an 18 mm2 cap turning into a whole sphere). By sampling the
# 3D curve, projecting with our own parametrization and UNWRAPPING u and v
# along the wire, the seam no longer exists: every wire is continuous by
# construction.


class SurfParam:
    """Explicit parametrization (identical to Geom_*Surface's)."""

    def __init__(self, prim: "Prim", ref: Optional[np.ndarray] = None, pole_axis: Optional[np.ndarray] = None):
        # ⚠️ cone with a negative semi-angle: STEP writes it as a generic
        # SURFACE_OF_REVOLUTION. The axis (and the slope's sign) is flipped
        # on a COPY: same geometry, positive semi-angle.
        if prim.kind == AXIAL and prim.slope < 0:
            prim = Prim(AXIAL, prim.center.copy(), -prim.axis, prim.r0, -prim.slope)
        self.prim = prim
        k = prim.kind
        C = np.asarray(prim.center, float)
        if k == FREE:
            ff = prim.free
            self.C, self.X, self.Y, self.Z = ff.C, ff.X, ff.Y, ff.Z
            self.periodic_u = self.periodic_v = False
            return
        if k == PLANE:
            Z = prim.axis / np.linalg.norm(prim.axis)
            X, Y = ortho_frame(Z)
        elif k == SPHERE:
            if pole_axis is None:
                pole_axis = np.array([0.0, 0.0, 1.0])
            Z = np.asarray(pole_axis, float)
            Z = Z / np.linalg.norm(Z)
            if ref is not None:
                X = np.asarray(ref, float) - float(np.asarray(ref, float) @ Z) * Z
            else:
                X = ortho_frame(Z)[0]
            if np.linalg.norm(X) < 1e-9:
                X = ortho_frame(Z)[0]
            X = X / np.linalg.norm(X)
            Y = np.cross(Z, X)
        else:
            Z = prim.axis / np.linalg.norm(prim.axis)
            if ref is not None:
                X = np.asarray(ref, float) - float(np.asarray(ref, float) @ Z) * Z
            else:
                X = ortho_frame(Z)[0]
            if np.linalg.norm(X) < 1e-9:
                X = ortho_frame(Z)[0]
            X = X / np.linalg.norm(X)
            Y = np.cross(Z, X)
        self.C, self.X, self.Y, self.Z = C, X, Y, Z
        self.periodic_u = k != PLANE
        self.periodic_v = k == TORUS
        if k == AXIAL:
            self.alpha = math.atan(prim.slope)

    # --- 3D -> (u, v) --------------------------------------------------------
    def uv(self, P: np.ndarray):
        P = np.atleast_2d(P)
        d = P - self.C
        k = self.prim.kind
        if k == FREE or k == PLANE:
            return d @ self.X, d @ self.Y
        t = d @ self.Z
        x, y = d @ self.X, d @ self.Y
        if k == SPHERE:
            r = np.maximum(np.linalg.norm(d, axis=1), 1e-300)
            return np.arctan2(y, x), np.arcsin(np.clip(t / r, -1.0, 1.0))
        u = np.arctan2(y, x)
        if k == TORUS:
            rho = np.hypot(x, y)
            return u, np.arctan2(t, rho - self.prim.r0)
        if abs(self.prim.slope) < 1e-12:
            return u, t
        return u, t / math.cos(self.alpha)

    # --- (u, v) -> 3D --------------------------------------------------------
    def point(self, u, v):
        u = np.atleast_1d(np.asarray(u, float))
        v = np.atleast_1d(np.asarray(v, float))
        k = self.prim.kind
        if k == FREE:
            return self.prim.free.point(u, v)
        cu, su = np.cos(u)[:, None], np.sin(u)[:, None]
        radial = cu * self.X + su * self.Y
        if k == PLANE:
            return self.C + u[:, None] * self.X + v[:, None] * self.Y
        if k == SPHERE:
            r = self.prim.r0
            return self.C + r * np.cos(v)[:, None] * radial + r * np.sin(v)[:, None] * self.Z
        if k == TORUS:
            R, r = self.prim.r0, self.prim.r1
            return self.C + (R + r * np.cos(v))[:, None] * radial + r * np.sin(v)[:, None] * self.Z
        if abs(self.prim.slope) < 1e-12:
            return self.C + self.prim.r0 * radial + v[:, None] * self.Z
        sa, ca = math.sin(self.alpha), math.cos(self.alpha)
        return self.C + (self.prim.r0 + v * sa)[:, None] * radial + (v * ca)[:, None] * self.Z

    def geom(self):
        k = self.prim.kind
        if k == FREE:
            return self.prim.free.geom()
        ax3 = gp_Ax3(_mk_pnt(self.C), _mk_dir(self.Z), _mk_dir(self.X))
        if k == PLANE:
            return Geom_Plane(ax3)
        if k == SPHERE:
            return Geom_SphericalSurface(ax3, float(self.prim.r0))
        if k == TORUS:
            return Geom_ToroidalSurface(ax3, float(self.prim.r0), float(self.prim.r1))
        if abs(self.prim.slope) < 1e-12:
            return Geom_CylindricalSurface(ax3, float(self.prim.r0))
        return Geom_ConicalSurface(ax3, float(self.alpha), float(self.prim.r0))

    def iso_u_curve(self, u: float):
        """3D curve at constant u (generatrix or meridian): (Geom_Curve, param==v)."""
        k = self.prim.kind
        if k == FREE:
            return None
        radial = math.cos(u) * self.X + math.sin(u) * self.Y
        if k == AXIAL:
            if abs(self.prim.slope) < 1e-12:
                return Geom_Line(_mk_pnt(self.C + self.prim.r0 * radial), _mk_dir(self.Z))
            sa, ca = math.sin(self.alpha), math.cos(self.alpha)
            return Geom_Line(_mk_pnt(self.C + self.prim.r0 * radial), _mk_dir(sa * radial + ca * self.Z))
        if k == TORUS:
            c = self.C + self.prim.r0 * radial
            ax2 = gp_Ax2(_mk_pnt(c), _mk_dir(np.cross(radial, self.Z)), _mk_dir(radial))
            return Geom_Circle(ax2, float(self.prim.r1))
        if k == SPHERE:
            ax2 = gp_Ax2(_mk_pnt(self.C), _mk_dir(np.cross(radial, self.Z)), _mk_dir(radial))
            return Geom_Circle(ax2, float(self.prim.r0))
        return None


class CGeom(_Curve):
    """
    GENERIC intersection curve (a GeomAPI_IntSS B-spline) wrapped with the
    same interface as the analytic curves: needed for the boundary between
    two fillets at different angles, or between a fillet and a round
    chamfer, where the intersection is neither a line nor a circle.
    Without this the boundary stayed the mesh's staircase between two
    smooth faces.
    """

    period = None

    def __init__(self, curve):
        self.c = curve
        self._proj = _m("GeomAPI").GeomAPI_ProjectPointOnCurve

    def param(self, P):
        out = []
        for p in np.atleast_2d(P):
            pr = self._proj(_mk_pnt(p), self.c)
            out.append(float(pr.LowerDistanceParameter()) if pr.NbPoints() > 0 else float("nan"))
        return np.array(out)

    def point(self, t):
        t = np.atleast_1d(t)
        return np.array([[self.c.Value(float(x)).X(), self.c.Value(float(x)).Y(), self.c.Value(float(x)).Z()] for x in t])

    def dist(self, P):
        out = []
        for p in np.atleast_2d(P):
            pr = self._proj(_mk_pnt(p), self.c)
            out.append(float(pr.LowerDistance()) if pr.NbPoints() > 0 else float("inf"))
        return np.array(out)

    def to_geom(self):
        return self.c

    def label(self):
        return "intersection"


def intersection_curves(surf_a, surf_b, tol: float) -> List[_Curve]:
    """Intersection curves between two Geom_Surface (GeomAPI_IntSS)."""
    out: List[_Curve] = []
    try:
        ss = _m("GeomAPI").GeomAPI_IntSS(surf_a, surf_b, float(tol))
        if ss.IsDone():
            for i in range(1, ss.NbLines() + 1):
                c = ss.Line(i)
                if c is not None and not c.IsNull() if hasattr(c, "IsNull") else c is not None:
                    out.append(CGeom(_keep(c)))
    except Exception as e:
        Log.debug(f"IntSS: {e}")
    return out


def _unwrap_to(vals: np.ndarray, start: float, period: float) -> np.ndarray:
    """Unwraps a periodic sequence so the first value is close to start."""
    out = np.unwrap(vals, period=period)
    k = round((start - out[0]) / period)
    return out + k * period


def curve_points(curve, t0: float, t1: float, n: int) -> np.ndarray:
    ts = np.linspace(t0, t1, n)
    return np.array([[curve.Value(t).X(), curve.Value(t).Y(), curve.Value(t).Z()] for t in ts]), ts


def make_pcurve(sp: SurfParam, curve, t0: float, t1: float, traverse_fwd: bool, u_prev: Optional[float], v_prev: Optional[float], n: int = 25):
    """
    Pcurve of the edge (3D curve `curve` between t0 and t1) on surface sp.
    Returns (Geom2d_Curve, (u_end, v_end) in the direction of travel, deviation).
    u_prev/v_prev = end of the previous edge in the wire: the pcurve is
    unwrapped to start from there.
    """
    P, ts = curve_points(curve, t0, t1, n)
    u, v = sp.uv(P)
    if not traverse_fwd:
        u, v = u[::-1], v[::-1]
    if sp.periodic_u:
        u = _unwrap_to(u, u_prev if u_prev is not None else float(u[0]), 2 * math.pi)
    if sp.periodic_v:
        v = _unwrap_to(v, v_prev if v_prev is not None else float(v[0]), 2 * math.pi)
    u_end, v_end = float(u[-1]), float(v[-1])
    if not traverse_fwd:
        u, v = u[::-1], v[::-1]
    # a line in (u,v) and linear in t?  -> degree-1 B-spline (exact)
    lin_u = u[0] + (u[-1] - u[0]) * (ts - t0) / max(t1 - t0, 1e-300)
    lin_v = v[0] + (v[-1] - v[0]) * (ts - t0) / max(t1 - t0, 1e-300)
    if float(np.abs(lin_u - u).max()) < 1e-9 and float(np.abs(lin_v - v).max()) < 1e-9:
        poles = _TColgp.TColgp_Array1OfPnt2d(1, 2)
        poles.SetValue(1, gp_Pnt2d(float(u[0]), float(v[0])))
        poles.SetValue(2, gp_Pnt2d(float(u[-1]), float(v[-1])))
        knots = _TColStd.TColStd_Array1OfReal(1, 2)
        knots.SetValue(1, float(t0))
        knots.SetValue(2, float(t1))
        mults = _TColStd.TColStd_Array1OfInteger(1, 2)
        mults.SetValue(1, 2)
        mults.SetValue(2, 2)
        c2d = Geom2d_BSplineCurve(poles, knots, mults, 1)
    else:
        pts = _TColgp.TColgp_HArray1OfPnt2d(1, n)
        prm = _TColStd.TColStd_HArray1OfReal(1, n)
        for i in range(n):
            pts.SetValue(i + 1, gp_Pnt2d(float(u[i]), float(v[i])))
            prm.SetValue(i + 1, float(ts[i]))
        it = _Geom2dAPI.Geom2dAPI_Interpolate(pts, prm, False, 1e-12)
        it.Perform()
        if not it.IsDone():
            return None, (u_end, v_end), math.inf
        c2d = it.Curve()
    # "same parameter" deviation: |C3d(t) - S(pcurve(t))| over a dense sample
    tt = np.linspace(t0, t1, 2 * n + 1)
    Q = np.array([[curve.Value(t).X(), curve.Value(t).Y(), curve.Value(t).Z()] for t in tt])
    uu = np.array([c2d.Value(t).X() for t in tt])
    vv = np.array([c2d.Value(t).Y() for t in tt])
    dev = float(np.linalg.norm(Q - sp.point(uu, vv), axis=1).max())
    return c2d, (u_end, v_end), dev


def edge_curve(edge):
    """(Geom_Curve, t0, t1) of the edge, with the location applied."""
    e = td_Edge(edge)
    loc = TopLoc_Location()
    c = bt_Curve(e, loc, 0.0, 0.0)
    t0, t1 = bt_Range(e)
    if c is None:
        return None, t0, t1
    if not loc.IsIdentity():
        c = c.Transformed(loc.Transformation())
    return c, float(t0), float(t1)


def edge_endpoints(edge) -> Tuple[np.ndarray, np.ndarray]:
    e = td_Edge(edge)
    return vpos(te_FirstVertex(e)), vpos(te_LastVertex(e))


def edge_is_analytic(edge) -> bool:
    try:
        return BRepAdaptor_Curve(td_Edge(edge)).GetType() in (GeomAbs_Line, GeomAbs_Circle, GeomAbs_Ellipse)
    except Exception:
        return False


def make_edge_on_curve(cv: "_Curve", V1, V2, P1: np.ndarray, P2: np.ndarray, Pmid: Optional[np.ndarray], closed: bool, Pnext: Optional[np.ndarray] = None):
    """
    Analytic edge between two EXISTING vertices (shared with the
    neighbors). For a circle/ellipse the right arc is the one passing
    through the chain's midpoint vertex (Pmid). Returns (edge, fwd): fwd =
    True if the FORWARD edge follows the chain's direction of travel
    (from V1). For a CLOSED chain the direction is read from the second
    vertex (Pnext).
    """
    g = _keep(cv.to_geom())
    if closed:
        t0 = float(cv.param(P1[None, :])[0])
        me = _keep(BRepBuilderAPI_MakeEdge(g, V1, V1, t0, t0 + cv.period))
        fwd = True
        if Pnext is not None:
            dt = (float(cv.param(Pnext[None, :])[0]) - t0) % cv.period
            fwd = dt < 0.5 * cv.period
        return (me.Edge() if me.IsDone() else None), fwd
    t1 = float(cv.param(P1[None, :])[0])
    t2 = float(cv.param(P2[None, :])[0])
    if cv.period:
        while t2 <= t1:
            t2 += cv.period
        if Pmid is not None:
            tm = float(cv.param(Pmid[None, :])[0])
            while tm < t1:
                tm += cv.period
            if tm > t2:  # the short arc doesn't pass through Pmid
                t1, t2 = t2 - cv.period, t1
                V1, V2 = V2, V1
                fwd_is_v1 = False
            else:
                fwd_is_v1 = True
        else:
            if t2 - t1 > math.pi:  # prefer the short arc
                t1, t2 = t2 - cv.period, t1
                V1, V2 = V2, V1
                fwd_is_v1 = False
            else:
                fwd_is_v1 = True
    else:
        if t2 < t1:
            t1, t2 = t2, t1
            V1, V2 = V2, V1
            fwd_is_v1 = False
        else:
            fwd_is_v1 = True
    if t2 - t1 < 1e-12:
        return None, fwd_is_v1
    me = _keep(BRepBuilderAPI_MakeEdge(g, V1, V2, t1, t2))
    if not me.IsDone():
        return None, fwd_is_v1
    return me.Edge(), fwd_is_v1


def composed_edge_orientations(face, edge) -> List:
    """
    Edge orientations in the face, COMPOSED with the wire's and the face's
    (TopExp_Explorer only returns the one stored in the wire). In a
    consistent shell, every edge is FORWARD in one face and REVERSED in the
    other, in composed terms.
    """
    out = []
    f_rev = face.Orientation() == TopAbs_REVERSED
    it_w = _TopoDS.TopoDS_Iterator(face, False, True)
    while it_w.More():
        w = it_w.Value()
        w_rev = w.Orientation() == TopAbs_REVERSED
        it_e = _TopoDS.TopoDS_Iterator(w, False, True)
        while it_e.More():
            ed = it_e.Value()
            if ed.IsSame(edge):
                e_rev = ed.Orientation() == TopAbs_REVERSED
                out.append(TopAbs_REVERSED if (f_rev ^ w_rev ^ e_rev) else TopAbs_FORWARD)
            it_e.Next()
        it_w.Next()
    return out


def grow_vertex_tolerance(V, dist: float) -> None:
    cur = float(bt_Tolerance(td_Vertex(V)))
    need = 1.2 * dist + 1e-7
    if need > cur:
        BRep_Builder().UpdateVertex(td_Vertex(V), need)


# =============================================================================
# 10. LOCAL REPLACEMENT ENGINE (one region at a time, with rollback)
# =============================================================================
#
# For each region:
#   1. boundary -> chains (consecutive edges with the same neighbor);
#   2. each chain: if the neighbor is analytic and the vertices lie on the
#      exact intersection curve -> ONE new analytic edge between the
#      EXISTING vertices; otherwise the polygonal edges are kept as they
#      are (shared with the tessellated neighbor, pcurve added in place);
#   3. new face = surface + wire (+ seam if closed 360 degrees);
#   4. checks on the face (validity, area);
#   5. BRepTools_ReShape: region's faces -> new face, chains -> new edges in
#      the neighbors too;
#   6. checks on the solid (free edges unchanged, valid neighboring faces,
#      consistent orientation, volume within the chords' slack);
#   7. if a check fails: the shape stays as it was before.


@dataclass
class Chain:
    edges: List[int]  # edge indices (Topo), in traversal order
    fwd: List[bool]  # True = traversed from FirstVertex to LastVertex
    verts: List[int]  # vertices in order (len = edges+1; closed: first==last)
    nb: int  # neighboring face (Topo index), -1 = free boundary
    closed: bool = False


@dataclass
class Loop:
    chains: List[Chain]


class Engine:
    def __init__(self, shape, tol_fit: float, diag: float, max_edge_tol: float, allow_polyline: bool = True, verbose: bool = False):
        self.shape = shape
        self.tol_fit = tol_fit
        self.tol_curve = 4.0 * tol_fit
        self.diag = diag
        self.max_edge_tol = max_edge_tol
        self.allow_polyline = allow_polyline
        self.verbose = verbose
        self.registry = TopTools_IndexedMapOfShape()  # analytic curved faces
        self.analytic: Dict[int, "Prim"] = {}  # registry id -> Prim
        self.topo = Topo(shape)
        self.free0 = count_free_edges(shape)
        self.vol0 = shape_volume(shape)
        self.n_ok = 0
        self.n_fail = 0
        self.log: List[str] = []

    # --- analytic primitives of the neighbors ---------------------------------
    def prim_of_face(self, i: int) -> Optional["Prim"]:
        f = self.topo.faces[i]
        k = self.registry.FindIndex(f)
        if k > 0 and k in self.analytic:
            return self.analytic[k]
        return plane_prim_of_face(self.topo, i)

    def _register(self, face, prim):
        k = self.registry.Add(face)
        self.analytic[k] = prim

    # --- boundary chains ---------------------------------------------------------
    def region_loops(self, rf: List[int]) -> Tuple[Optional[List[Loop]], str]:
        topo = self.topo
        rset = set(rf)
        bnd: Dict[int, int] = {}  # edge -> region face
        for i in rf:
            for k in topo.f_edges[i]:
                fs = topo.e_faces[k]
                inside = [f for f in fs if f in rset]
                if len(fs) < 2:
                    return None, "free boundary of the mesh"
                if len(inside) == 1:
                    bnd[k] = i
        if not bnd:
            return None, "no boundary"
        # direction of travel: the face's interior to the left relative to
        # the outward normal (independent of OCC's orientation conventions)
        dir_edge: Dict[int, Tuple[int, int]] = {}
        out_v: Dict[int, List[int]] = defaultdict(list)
        for k, i in bnd.items():
            a, b = topo.e_verts[k]
            A, B = topo.vpos[a], topo.vpos[b]
            n = topo.norms[i]
            inward = np.cross(n, B - A)
            if float((topo.cents[i] - 0.5 * (A + B)) @ inward) < 0:
                a, b = b, a
            dir_edge[k] = (a, b)
            out_v[a].append(k)

        # traversal by rotating around the vertex (handles degree-4 nodes)
        def next_edge(k_in: int) -> Optional[int]:
            a, b = dir_edge[k_in]
            cands = out_v.get(b, [])
            if len(cands) == 1:
                return cands[0]
            if not cands:
                return None
            f = bnd[k_in]
            seen = set()
            for _ in range(64):
                es_at_b = [e for e in topo.f_edges[f] if b in topo.e_verts[e] and e != k_in]
                es_at_b = [e for e in es_at_b if e not in seen]
                if not es_at_b:
                    return None
                e = es_at_b[0]
                seen.add(e)
                if e in bnd:
                    return e if bnd[e] == f and dir_edge[e][0] == b else None
                nxt = [g for g in topo.e_faces[e] if g in rset and g != f]
                if not nxt:
                    return None
                f, k_in = nxt[0], e
            return None

        used = set()
        loops: List[List[int]] = []
        for k0 in sorted(bnd):
            if k0 in used:
                continue
            seq, k = [], k0
            while k is not None and k not in used:
                used.add(k)
                seq.append(k)
                k = next_edge(k)
            if k != k0 and (not seq or dir_edge[seq[-1]][1] != dir_edge[seq[0]][0]):
                return None, "boundary can't be closed"
            loops.append(seq)
        # chains: split where the neighbor changes
        out: List[Loop] = []
        for seq in loops:
            nbs = []
            for k in seq:
                fs = [f for f in topo.e_faces[k] if f not in rset]
                nbs.append(fs[0] if fs else -1)
            n = len(seq)
            starts = [j for j in range(n) if nbs[j] != nbs[j - 1]]
            chains: List[Chain] = []
            if not starts:
                ch = Chain([], [], [dir_edge[seq[0]][0]], nbs[0], closed=True)
                for k in seq:
                    ch.edges.append(k)
                    ch.fwd.append(dir_edge[k][0] == topo.e_verts[k][0])
                    ch.verts.append(dir_edge[k][1])
                chains.append(ch)
            else:
                for si, s in enumerate(starts):
                    e_ = starts[(si + 1) % len(starts)]
                    idxs = list(range(s, e_)) if e_ > s else list(range(s, n)) + list(range(0, e_))
                    ch = Chain([], [], [dir_edge[seq[idxs[0]]][0]], nbs[s])
                    for j in idxs:
                        k = seq[j]
                        ch.edges.append(k)
                        ch.fwd.append(dir_edge[k][0] == topo.e_verts[k][0])
                        ch.verts.append(dir_edge[k][1])
                    chains.append(ch)
            out.append(Loop(chains))
        return out, ""

    # --- converting a region -----------------------------------------------------
    def _safe(self, *a, **k) -> Tuple[bool, str]:
        """
        _convert() with the safety net.
        ⚠️ NOTHING DIES HERE. The engine works on a COPY (self.shape only
        changes after validation), so a region that blows up halfway
        through is simply a discarded region: the part stays as it was and
        the run continues. Without this, a single sick primitive used to
        throw away all the work already done - on test9 with --tol 0.0005
        it was a cone degenerated by the refinement (see LM_MAX_SLOPE) and
        the script died after twenty minutes of Phase C.
        The only things the engine touches IN PLACE are the tolerances of
        shared vertices and edges, and those have a backup: they're
        restored here.
        """
        try:
            return self._convert(*a, **k)
        except KeyboardInterrupt:
            raise
        except Exception as ex:
            try:
                self._restore_all(getattr(self, "_vtol_backup", None) or {})
            except Exception:
                pass
            self._vtol_backup = {}
            Log.debug(f"    !! exception during conversion: {type(ex).__name__}: {ex}")
            return False, f"exception {type(ex).__name__}"

    def convert(self, R: Region, strict_hole: bool = False) -> Tuple[bool, str]:
        mt0 = max_tolerance(self.shape) if self.verbose else 0.0
        ok, why = self._safe(R, strict_hole)
        # ⚠️ SECOND STRATEGY: if the planar neighbor rebuilt with the exact
        # arc comes out invalid (or the orientation doesn't work out, which
        # is the same case seen from the other attempt), retry leaving the
        # boundaries with the planes POLYGONAL: the analytic face still
        # converts, the boundary stays the mesh's.
        retry = "neighbor" in why or "orientation" in why or "volume" in why
        if not ok and not strict_hole and retry:
            ok2, why2 = self._safe(R, strict_hole, analytic_planes=False)
            if ok2:
                ok, why = ok2, why2
                R.note += " · boundaries with the planes left polygonal"
            elif "neighbor" in why2 or "orientation" in why2 or "volume" in why2:
                # ⚠️ THIRD STRATEGY: no new edge, anywhere. Needed when the
                # neighbor is a CLOSED analytic face (an already-converted
                # hole, with its own seam): rebuilding its boundary breaks
                # it ("wire:NotConnected"). This way the region still
                # becomes an exact surface and the contour stays IDENTICAL
                # to the mesh: no neighbor is touched, no new sliver.
                ok3, why3 = self._safe(R, strict_hole, analytic_planes=False, analytic_curved=False)
                if ok3:
                    ok, why = ok3, why3
                    R.note += " · contour left polygonal"
        # ⚠️ FOURTH STRATEGY: no curves "stretched through the vertices".
        # When the exact intersection doesn't exist, the code falls back
        # to a circle (or a line) passing through the chain's vertices.
        # That circle sits within a tenth of the tolerance of the surface,
        # but does NOT lie on it: projected in (u,v) it pokes a couple of
        # tenths of a millimeter outside the contour and the wire comes
        # out self-intersecting. It costs little to retry leaving that
        # chain polygonal: the analytic face still converts and it's much
        # better than a tessellated region.
        if not ok and not strict_hole and ("SelfIntersecting" in why or "Intersecting" in why):
            ok4, why4 = self._safe(R, strict_hole, fitted_curves=False)
            if ok4:
                ok, why = ok4, why4
                R.note += " · approximated boundaries left polygonal"
        if self.verbose:
            mt1 = max_tolerance(self.shape)
            if mt1 > max(mt0 * 1.5, self.max_edge_tol):
                Log.debug(f"    !! max tolerance rose from {mt0:.1e} to {mt1:.1e} ({'accepted' if ok else 'discarded'})")
        if ok:
            self.n_ok += 1
            R.status = "OK"
        else:
            self.n_fail += 1
            R.status = "rejected: " + why
        if self.verbose or not ok:
            Log.debug(f"  {R.label()} -> {R.status}")
        return ok, why

    def _convert(self, R: Region, strict_hole: bool, analytic_planes: bool = True, analytic_curved: bool = True, fitted_curves: bool = True) -> Tuple[bool, str]:
        topo = self.topo
        prim = R.prim
        rf = [topo.face_index(f) for f in R.fobjs]
        if any(i < 0 for i in rf):
            return False, "region's faces no longer present"
        rset = set(rf)
        self._rset = rset
        self._rf = rf
        loops, why = self.region_loops(rf)
        # ⚠️ BOUNDARY CAN'T BE CLOSED: almost always a pinch vertex, where
        # two stretches of the same region's boundary touch (sliver
        # facets). The facets passing through those vertices are removed,
        # the largest component is kept and it's retried: the bulk of the
        # region converts, the slivers stay tessellated.
        for _ in range(3):
            if loops is not None:
                break
            pinch = self._pinch_vertices(rf)
            if not pinch:
                break
            keep_ = [i for i in rf if not (set(topo.e_verts[k][0] for k in topo.f_edges[i]) | set(topo.e_verts[k][1] for k in topo.f_edges[i])) & pinch]
            keep_ = largest_component(keep_, topo.adj)
            if len(keep_) < max(4, len(rf) // 2):
                break
            rf = keep_
            rset = set(rf)
            self._rset = rset
            self._rf = rf
            R.faces = [i for i in rf]
            R.fobjs = [topo.faces[i] for i in rf]
            loops, why = self.region_loops(rf)
        if loops is None:
            return False, why
        if strict_hole:
            if not (prim.kind == AXIAL and abs(prim.slope) < 1e-12 and R.closed_u and R.concave):
                return False, "not a through/blind cylindrical hole"
            if len(loops) != 2 or any(len(L.chains) != 1 or not L.chains[0].closed for L in loops):
                desc = "; ".join(
                    f"ring {k}: {len(L.chains)} chains, neighbors "
                    + ",".join(("plane" if (ch.nb >= 0 and topo.planar[ch.nb] and topo.nverts[ch.nb] > 4) else f"facet{topo.nverts[ch.nb] if ch.nb >= 0 else ''}") for ch in L.chains[:6])
                    for k, L in enumerate(loops)
                )
                return False, f"the hole's boundary isn't made of two simple rings ({desc})"
            for L in loops:
                nb = L.chains[0].nb
                if nb < 0 or not topo.planar[nb]:
                    return False, "the hole doesn't end on planar faces"
                if abs(abs(float(topo.norms[nb] @ prim.axis)) - 1.0) > 1e-4:
                    return False, "opening face not orthogonal to the axis"
        if R.closed_u and len(loops) != 2:
            return False, f"region closed 360 degrees with {len(loops)} rings (expected 2)"
        if not R.closed_u and prim.kind in (AXIAL, TORUS) and R.coverage > 0.97:
            return False, "ambiguous angular coverage"

        # --- surface: frame chosen away from the seam ----------------------------
        C = np.array([topo.cents[i] for i in rf])
        W = np.array([max(topo.areas[i], 1e-12) for i in rf])
        ref = None
        pole_axis = None
        if prim.kind == SPHERE:
            d = C - prim.center
            d = d / np.maximum(np.linalg.norm(d, axis=1), 1e-12)[:, None]
            m = (d * W[:, None]).sum(axis=0) / W.sum()
            if float(np.linalg.norm(m)) < 0.2:
                return False, "sphere: region too wide (more than a hemisphere)"
            ref = m / np.linalg.norm(m)
            # polar axis: orthogonal to the mean direction, as far as
            # possible from the region's points (poles outside the face)
            u_, v_ = ortho_frame(ref)
            phi = np.linspace(0.0, math.pi, 91)
            A = np.cos(phi)[:, None] * u_ + np.sin(phi)[:, None] * v_
            worst = np.abs(A @ d.T).max(axis=1)
            kbest = int(np.argmin(worst))
            if worst[kbest] > math.cos(math.radians(20.0)):
                return False, "sphere: pole too close to the region"
            pole_axis = A[kbest]
        elif prim.kind in (AXIAL, TORUS):
            d = C - prim.center
            d = d - np.outer(d @ prim.axis, prim.axis)
            nn = np.linalg.norm(d, axis=1)
            good = nn > 1e-9
            if good.any():
                m = ((d[good] / nn[good, None]) * W[good, None]).sum(axis=0) / W[good].sum()
                if float(np.linalg.norm(m)) > 1e-6:
                    ref = m / np.linalg.norm(m)
        sp = SurfParam(prim, ref, pole_axis)
        surf = _keep(sp.geom())
        if R.closed_u:
            why = self._align_closed_chains(loops, sp)
            if why:
                return False, why

        # --- chains -> edges -----------------------------------------------------
        new_edges: Dict[int, Tuple[object, bool]] = {}  # id(chain) -> (edge, fwd)
        replaced: List[Tuple[Chain, object, bool]] = []
        n_analytic = n_poly = n_reused = 0
        vtol_backup: Dict[int, float] = {}
        self._vtol_backup = vtol_backup  # same dict: _safe() needs it

        def vertex_obj(j):
            return td_Vertex(topo.vmap.FindKey(j + 1))

        def bump_vertex(j, dist):
            V = vertex_obj(j)
            if j not in vtol_backup:
                vtol_backup[j] = float(bt_Tolerance(V))
            grow_vertex_tolerance(V, dist)

        med_area = float(np.median([topo.areas[i] for i in rf]))
        len_tot = len_debris = 0.0
        for L in loops:
            for ch in L.chains:
                Lc = float(sum(edge_length(topo.edges[k]) for k in ch.edges))
                len_tot += Lc
                if ch.nb >= 0 and topo.nverts[ch.nb] <= 4 and topo.areas[ch.nb] < 3.0 * med_area:
                    len_debris += Lc
        # ⚠️ FRAGMENT-REGION AMID TESSELLATION. A primitive fitted on four
        # or five facets, with its boundary resting almost entirely on
        # other loose facets, isn't a real feature of the part: it's the
        # noise of a fillet zone. Converting it produces exactly the
        # "lots of irregularly shaped faces" (random radii, ø1.4589,
        # ø2.0204...) that clutter the model. Better to leave the mesh
        # there: it gets cleaned up by hand.
        if not strict_hole and len(rf) < 12 and not R.repeated and len_debris > 0.6 * max(len_tot, 1e-9):
            return False, (f"small region isolated among tessellated facets ({len(rf)} faces, {100.0 * len_debris / max(len_tot, 1e-9):.0f}% of the boundary): left as mesh")
        for L in loops:
            for ch in L.chains:
                if ch.nb < 0:
                    return False, "free boundary"
                if len(ch.edges) == 1 and not ch.closed:
                    # ⚠️ CHAMFER NEXT TO A FILLET. A single edge is usually
                    # already good and gets reused as is. But when the
                    # neighbor is a wide band tessellated with a single
                    # chord (the end of a chamfer, the 10 mm-tall quad of a
                    # vertical fillet), that chord CUTS through the new
                    # surface: half a hundredth of a millimeter below the
                    # sphere. The face then gets rejected for tolerance,
                    # and the tessellation stays in its place. Here we
                    # check how far the edge really is from the surface:
                    # if it's more than the fit tolerance, we go on and try
                    # to rebuild it with the exact curve (sphere x plane =
                    # circle), which lies both on the new face and on the
                    # neighbor's plane.
                    _dv = self._edge_dev(ch.edges[0], prim)
                    if self.verbose:
                        Log.debug(
                            f"    spigolo singolo L={edge_length(topo.edges[ch.edges[0]]):.4f} "
                            f"scarto dalla superficie {_dv:.2e} (tol_fit {self.tol_fit:.2e}) "
                            f"vicino f{ch.nb} area {topo.areas[ch.nb]:.3f}"
                        )
                    if _dv <= self.tol_fit:
                        n_reused += 1
                        continue
                if len(ch.edges) == 1 and ch.closed and edge_is_analytic(topo.edges[ch.edges[0]]):
                    n_reused += 1
                    continue
                nbp = self.prim_of_face(ch.nb)
                if nbp is not None and nbp.kind == FREE:
                    nbp = None  # no exact curves against a free-form shape
                if not analytic_curved and nbp is not None and nbp.kind != PLANE:
                    nbp = None
                # ⚠️ a planar neighbor as small as a single facet isn't a
                # real plane of the part: it's a scrap of curved surface
                # flattened by Phase A. Replacing its boundary with an
                # exact arc makes it self-intersect: the polyline is kept.
                if nbp is not None and nbp.kind == PLANE and (not analytic_planes or topo.areas[ch.nb] < 20.0 * med_area):
                    nbp = None
                cv = None
                P = topo.vpos[ch.verts]
                if nbp is not None:
                    # ⚠️ HOW FAR THE ARC CAN STICK OUT. The yardstick is the
                    # mesh's sag in the region: the boundary's polyline is a
                    # chord of the true curve, so it deviates the way the
                    # facets deviate from the surface. With a wider
                    # yardstick even the wrong curve passes: on a 0.65 mm
                    # boundary, a radius-0.36 arc sticking out by 0.2 mm
                    # (almost a semicircle) was being accepted, the pcurves
                    # crossed and the face came out "SelfIntersectingWire".
                    cands = surf_surf_curves(prim, nbp, self.tol_fit)
                    cv = choose_curve(cands, P, self.tol_curve, scale=self.diag, arc_tol=max(self.tol_curve, 2.0 * R.sag))
                    # ⚠️ TANGENT SURFACES (fillet-plane, sphere-fillet): the
                    # exact intersection curve is poorly conditioned, a
                    # vertex 1e-6 from both surfaces can be 3e-3 from the
                    # curve. So the curve passing through the mesh's
                    # vertices is taken instead: it's consistent with the
                    # input, and the edge's tolerance measures its true
                    # deviation.
                    if cv is None and fitted_curves and len(ch.verts) >= 4:
                        Pq = P[:-1] if ch.closed else P
                        cv = fit_curve(Pq, self.tol_curve, arc_tol=max(self.tol_curve, 2.0 * R.sag))
                    # ⚠️ CURVED neighbor already converted and no analytic
                    # curve: generic intersection of the two surfaces (B-spline)
                    if cv is None and nbp.kind != PLANE and not ch.closed and len(ch.verts) >= 3:
                        # ⚠️ the mesh's boundary between two regions is a
                        # STAIRCASE of facets: interior vertices sit up to
                        # a facet away from the true curve. Only the
                        # ENDPOINTS matter (which stay as vertices):
                        # interior vertices disappear along with the curve.
                        try:
                            sn = bt_Surface(topo.faces[ch.nb], TopLoc_Location())
                            gen = intersection_curves(surf, sn, self.tol_fit)
                            seg_len = float(np.median(np.linalg.norm(np.diff(P, axis=0), axis=1)))
                            loose = max(self.tol_curve, 2.0 * seg_len, 2.0 * R.sag)
                            cv = choose_curve(gen, P, loose, arc_tol=loose)
                            if self.verbose:
                                Log.debug(
                                    f"    IntSS: {len(gen)} curves, deviations "
                                    f"{[f'{float(np.abs(g.dist(P)).max()):.1e}/{arc_deviation(g, P):.1e}' for g in gen]} "
                                    f"tol {loose:.1e} -> {None if cv is None else 'ok'}"
                                )
                            if cv is not None:
                                d_ends = float(np.abs(cv.dist(P[[0, -1]])).max())
                                if d_ends > self.max_edge_tol:
                                    cv = None
                        except Exception as ex:
                            Log.debug(f"generic intersection: {ex}")
                if self.verbose:
                    _c = surf_surf_curves(prim, nbp, self.tol_fit) if nbp is not None else []
                    _d = [f"{c.label()}:{float(np.abs(c.dist(P)).max()):.1e}/{arc_deviation(c, P):.1e}" for c in _c]
                    _f = fit_curve(P[:-1] if ch.closed else P, self.tol_curve, arc_tol=max(self.tol_curve, 2.0 * R.sag))
                    Log.debug(
                        f"    chain {len(ch.edges)} edges, neighbor {'-' if nbp is None else nbp.label()}"
                        f"{'' if nbp is None or nbp.kind != PLANE else f' area {topo.areas[ch.nb]:.2f}'}: "
                        f"candidates {_d} fit {None if _f is None else _f.label()} "
                        f"-> {None if cv is None else cv.label()}  (tol {self.tol_curve:.1e}, arc {max(self.tol_curve, 4.0 * R.sag):.1e})"
                    )
                if cv is None:
                    if strict_hole:
                        return False, "the hole's boundary isn't a circle"
                    if not self.allow_polyline:
                        return False, "polygonal boundary not convertible"
                    if len(ch.edges) == 1:
                        n_reused += 1  # the mesh's edge stays
                    else:
                        n_poly += 1
                    continue
                if ch.closed and cv.period is None:
                    n_poly += 1
                    continue
                j1, j2 = ch.verts[0], ch.verts[-1]
                V1, V2 = vertex_obj(j1), vertex_obj(j2)
                for j in (j1, j2):
                    bump_vertex(j, float(np.abs(cv.dist(topo.vpos[j][None, :]))[0]))
                Pmid = topo.vpos[ch.verts[len(ch.verts) // 2]] if len(ch.verts) > 2 else None
                Pnext = topo.vpos[ch.verts[1]] if len(ch.verts) > 1 else None
                e, fwd_is_v1 = make_edge_on_curve(cv, V1, V2, topo.vpos[j1], topo.vpos[j2], Pmid, ch.closed, Pnext)
                if e is None:
                    self._restore_vertices(vtol_backup)
                    return False, f"MakeEdge failed on {cv.label()}"
                e = td_Edge(e)
                # ⚠️ if the neighbor is an already-converted CURVED analytic
                # face, the new edge needs the pcurve on ITS surface too,
                # and needs it RIGHT AWAY: BRepCheck on the rebuilt
                # neighboring face, without the pcurve, answers
                # "UnorientableShape".
                # ⚠️ even a PLANAR neighbor needs the pcurve when the new
                # edge isn't a line/circle/ellipse: OCC's plane projection
                # can't compute a hyperbola's on its own, and the
                # neighboring face comes out "UnorientableShape".
                if nbp is not None and (nbp.kind != PLANE or isinstance(cv, CHyperbola)):
                    why = self._pcurve_on_neighbor(e, ch.nb, nbp)
                    if why:
                        self._restore_vertices(vtol_backup)
                        return False, why
                new_edges[id(ch)] = (e, fwd_is_v1)
                replaced.append((ch, e, fwd_is_v1))
                n_analytic += 1

        if self.verbose:
            Log.debug(f"    BOUNDARY {R.label()[:34]}: {len_tot:.2f} mm, of which {len_debris:.2f} mm ({100.0 * len_debris / max(len_tot, 1e-9):.0f}%) against tessellated facets")
        # --- seam for closed regions --------------------------------------------------
        seam = None
        if not R.closed_u:
            self._seam_path = None
        if R.closed_u:
            seam, why = self._make_seam(R, sp, loops, new_edges, bump_vertex)
            if seam is None:
                self._restore_vertices(vtol_backup)
                return False, why

        # --- wire and face -----------------------------------------------------------
        bb = BRep_Builder()
        F = TopoDS_Face()
        bb.MakeFace(F, surf, 1e-7)
        # wire direction: counterclockwise relative to the surface's
        # NATURAL normal. The mesh's traversal is counterclockwise relative
        # to the OUTWARD normal: if the surface is concave (a hole) it must
        # be flipped.
        flip = R.concave
        self._etol_backup = {}
        wires = self._build_wires(F, sp, loops, new_edges, seam, flip)
        if isinstance(wires, str):
            self._restore_all(vtol_backup)
            return False, wires
        for w, _, _ in wires:
            bb.Add(F, w)
        # ⚠️ THE CEILING MUST BE MEASURED AGAINST THE MESH, NOT THE PART. A
        # boundary left as it was deviates from the new surface by exactly
        # the mesh's chord sag: on a coarse tessellation (an R1.5 fillet
        # split into four 22-degree facets) that's three hundredths, and an
        # absolute ceiling rejects all of them. But those three hundredths
        # are already in the input: rejecting the region gains no
        # precision, it only loses the exact surface and the facets
        # remain. So the ceiling is the wider of the absolute one and the
        # region's sag, which is exactly the deviation the mesh carries
        # with it. ...but with a limit: on a botched fit (a 143 mm
        # "cylinder" stretched over four huge facets) the sag is
        # millimeter-sized, and without a hard ceiling the region would
        # pass through, deforming the part.
        tol_cap = max(self.max_edge_tol, min(2.0 * R.sag, 5.0 * self.max_edge_tol))
        worst_edge_tol = max(t for _, _, t in wires) if wires else 0.0
        if worst_edge_tol > tol_cap:
            self._restore_all(vtol_backup)
            return False, (f"polygonal edge tolerance {worst_edge_tol:.1e} > ceiling {tol_cap:.1e}")

        det = check_detail(F)
        if det:
            if self.verbose:
                self._dump_face(F)
            self._restore_all(vtol_backup)
            return False, "invalid face: " + ", ".join(det)
        mt = max_tolerance(F)
        # the face can carry along the tolerance the edges already had:
        # we're not the ones who put it there
        if mt > max(tol_cap, max((t for _, t, _ in wires), default=0.0)):
            self._restore_all(vtol_backup)
            return False, f"tolerance {mt:.1e} on the new face beyond the ceiling {tol_cap:.1e}"
        a_mesh = float(sum(topo.areas[i] for i in rf))
        a_new = face_area(F)
        perim = sum(edge_length(topo.edges[k]) for L in loops for ch in L.chains for k in ch.edges)
        if abs(a_new - a_mesh) > 0.03 * a_mesh + 2.0 * R.sag * perim + 1e-9:
            self._restore_all(vtol_backup)
            return False, f"inconsistent area: {a_new:.4f} against the mesh's {a_mesh:.4f} mm2"

        # --- replacement in the solid --------------------------------------------
        # neighbors' wires BEFORE the replacement: a neighbor that ends up
        # with more of them has a spurious ring (old edges left + new arc = notch)
        self._nb_wires = {}
        for ch, _, _ in replaced:
            if ch.nb not in self._nb_wires:
                self._nb_wires[ch.nb] = (topo.faces[ch.nb], count_sub(topo.faces[ch.nb], TopAbs_WIRE))
        for attempt in range(2):
            Fo = F if attempt == 0 else td_Face(F.Reversed())
            rs = BRepTools_ReShape()
            # keys always FORWARD: ReShape composes the new orientation with
            # the RELATIVE one between the key and its occurrence in the solid
            rs.Replace(td_Face(topo.faces[rf[0]].Oriented(TopAbs_FORWARD)), Fo)
            for i in rf[1:]:
                rs.Remove(td_Face(topo.faces[i].Oriented(TopAbs_FORWARD)))
            for ch, e, fwd_is_v1 in replaced:
                first = topo.edges[ch.edges[0]]
                # e FORWARD goes from V1 to V2 if fwd_is_v1; the chain
                # traverses first in direction ch.fwd[0]. The new edge's
                # orientation relative to the old one: the same if both go
                # the same way.
                same = fwd_is_v1 == ch.fwd[0]
                rs.Replace(first, e if same else td_Edge(e.Reversed()))
                for k in ch.edges[1:]:
                    rs.Remove(topo.edges[k])
            try:
                new_shape = rs.Apply(self.shape)
            except Exception as ex:
                self._restore_all(vtol_backup)
                return False, f"ReShape: {type(ex).__name__}: {ex}"
            ok, why = self._validate(new_shape, Fo, rf, replaced, R, a_mesh, sp)
            if self.verbose:
                Log.debug(f"    attempt {attempt + 1} (Fo {Fo.Orientation()}): {'ok' if ok else why}")
            if ok:
                break
            # ⚠️ when NO boundary was replaced (all chains polygonal) the
            # orientation check has no new edges to look at and says
            # nothing: the error only shows up on the VOLUME. Even in that
            # case the flipped face has to be tried.
            if attempt == 1 or ("orient" not in why and "volume" not in why):
                self._restore_all(vtol_backup)
                return False, why
        # --- accepted --------------------------------------------------------------
        # ⚠️ and the backups are thrown away: the inflated tolerances are
        # now final, and their keys are edges and vertices of the OLD
        # topology. If they stayed there, _safe()'s emergency rollback on
        # the next region would bring them back and break this face.
        self._etol_backup = {}
        self._vtol_backup = {}
        self._carry_registry(rs, new_shape)
        self._register(Fo, prim)
        self.shape = new_shape
        self.topo = Topo(new_shape, prev=self.topo)
        R.note = f"{n_analytic} analytic edges · {n_reused} reused · {n_poly} polygonal"
        return True, ""

    def _seam_path_search(self, targets_a: set, targets_b: set, sp: SurfParam, r_eff: float):
        """
        A path of edges INTERNAL to the region (both faces inside the
        region) from a vertex of the first ring to one of the second, as
        "vertical" as possible (cost = u deviation x radius). Returns
        (A, B, [(edge index, fwd), ...]) or None.
        """
        import heapq

        topo = self.topo
        rset = self._rset
        internal = [k for i in self._rf for k in topo.f_edges[i] if len(topo.e_faces[k]) == 2 and all(f in rset for f in topo.e_faces[k])]
        internal = sorted(set(internal))
        if not internal:
            return None
        adjv: Dict[int, List[Tuple[int, int]]] = defaultdict(list)
        for k in internal:
            a, b = topo.e_verts[k]
            adjv[a].append((b, k))
            adjv[b].append((a, k))
        best = None
        for A in targets_a:
            if A not in adjv:
                continue
            uA = float(sp.uv(topo.vpos[A][None, :])[0][0])
            dist = {A: 0.0}
            prev: Dict[int, Tuple[int, int]] = {}
            heap = [(0.0, A)]
            found = None
            while heap:
                d, v = heapq.heappop(heap)
                if d > dist.get(v, math.inf):
                    continue
                if v in targets_b and v != A:
                    found = v
                    break
                for w, k in adjv[v]:
                    uw = float(sp.uv(topo.vpos[w][None, :])[0][0])
                    du = abs(float(np.angle(np.exp(1j * (uw - uA)))))
                    nd = d + du * r_eff + 1e-6
                    if nd < dist.get(w, math.inf):
                        dist[w] = nd
                        prev[w] = (v, k)
                        heapq.heappush(heap, (nd, w))
            if found is None:
                continue
            if best is None or dist[found] < best[0]:
                edges = []
                v = found
                while v != A:
                    pv, k = prev[v]
                    edges.append((k, topo.e_verts[k][0] == pv))
                    v = pv
                edges.reverse()
                best = (dist[found], A, found, edges)
        if best is None:
            return None
        _, A, B, edges = best
        # the maximum u deviation along the path must stay small
        uA = float(sp.uv(topo.vpos[A][None, :])[0][0])
        for k, _ in edges:
            for j in topo.e_verts[k]:
                uj = float(sp.uv(topo.vpos[j][None, :])[0][0])
                # the path's edges are REAL mesh edges, they lie on the
                # surface: a deviation in u isn't an error, the path just
                # must not go all the way around (a quarter turn at most)
                if abs(float(np.angle(np.exp(1j * (uj - uA))))) > math.pi / 4:
                    return None
        return A, B, edges

    def _edge_dev(self, k: int, prim: "Prim") -> float:
        """How far the edge (its 3D curve, not just the endpoints) strays
        from the region's surface."""
        try:
            curve, t0, t1 = edge_curve(self.topo.edges[k])
            if curve is None:
                return 0.0
            P, _ = curve_points(curve, t0, t1, 7)
            return float(np.abs(prim.dist(P)).max())
        except Exception:
            return 0.0

    def _pinch_vertices(self, rf: List[int]) -> set:
        """Vertices touched by more than two of the region's boundary edges."""
        topo = self.topo
        rset = set(rf)
        cnt: Dict[int, int] = defaultdict(int)
        for i in rf:
            for k in topo.f_edges[i]:
                fs = topo.e_faces[k]
                if len(fs) == 2 and sum(1 for f in fs if f in rset) == 1:
                    for j in topo.e_verts[k]:
                        cnt[j] += 1
        return {j for j, c in cnt.items() if c > 2}

    def _align_closed_chains(self, loops, sp: SurfParam) -> str:
        """
        Region closed 360 degrees: the seam will run from vertex A of the
        first ring to vertex B of the second, at the SAME u. A closed chain
        that will become a circle has just ONE vertex (the starting one):
        here it's chosen so the two rings line up, BEFORE the edges are
        created.
        """
        topo = self.topo

        def rotate(ch: Chain, j: int) -> bool:
            if not ch.closed or j not in ch.verts[:-1]:
                return False
            s = ch.verts.index(j)
            ch.edges = ch.edges[s:] + ch.edges[:s]
            ch.fwd = ch.fwd[s:] + ch.fwd[:s]
            core = ch.verts[:-1]
            core = core[s:] + core[:s]
            ch.verts = core + [core[0]]
            return True

        def candidates(L: Loop) -> List[int]:
            out = []
            for ch in L.chains:
                out.extend(ch.verts[:-1] if ch.closed else [ch.verts[0], ch.verts[-1]])
            return sorted(set(out))

        L1, L2 = loops
        # A: if the first ring is a closed chain, any of its vertices will
        # do: the one that aligns best with the second ring is picked
        c1, c2 = candidates(L1), candidates(L2)
        if not c1 or not c2:
            return "rings with no vertices"
        u1, _ = sp.uv(topo.vpos[c1])
        u2, _ = sp.uv(topo.vpos[c2])
        D = np.abs(np.angle(np.exp(1j * (u2[None, :] - u1[:, None]))))
        ia, ib = np.unravel_index(int(np.argmin(D)), D.shape)
        A, B = c1[ia], c2[ib]
        # ⚠️ HOLES WITH MULTIPLE ROWS OF FACETS: the two rings' vertices
        # aren't aligned (deviation up to 0.1 mm) and a straight seam
        # doesn't exist. Then the seam is a PATH of mesh edges internal to
        # the region, from A to a vertex of the second ring: real edges, on
        # the surface, with both pcurves computed the same way as any other
        # edge.
        self._seam_path = None
        r_eff = max(prim_radius(sp.prim), 1e-6)
        if D[ia, ib] * r_eff > self.max_edge_tol:
            path = self._seam_path_search(set(c1), set(c2), sp, r_eff)
            if path is None:
                return f"seam: no aligned vertex on the second ring (deviation {D[ia, ib] * r_eff:.1e} mm)"
            A, B, edges = path
            self._seam_path = (A, B, edges)
        for ch in L1.chains:
            if ch.closed:
                rotate(ch, A)
        for ch in L2.chains:
            if ch.closed:
                rotate(ch, B)
        return ""

    def _dump_face(self, F) -> None:
        """Diagnostics: pcurve of every edge of every wire (endpoints in u,v)."""
        cos_ = _st(BRep_Tool, "CurveOnSurface")
        for wi, w in enumerate(explore(F, TopAbs_WIRE)):
            Log.debug(f"    wire {wi} orient {w.Orientation()} closed3d={_st(BRep_Tool, 'IsClosed')(w) if hasattr(BRep_Tool, 'IsClosed_s') else '?'}")
            ex = TopExp_Explorer(w, TopAbs_EDGE)
            while ex.More():
                e = td_Edge(ex.Current())
                try:
                    c2d = cos_(e, F, 0.0, 0.0)
                    t0, t1 = bt_Range(e)
                    a, b = c2d.Value(t0), c2d.Value(t1)
                    Pa, Pb = edge_endpoints(e)
                    Log.debug(f"      edge {e.Orientation()} tol={bt_Tolerance(e):.1e} uv({a.X():+.4f},{a.Y():+.4f})->({b.X():+.4f},{b.Y():+.4f}) 3d {np.round(Pa, 3)}->{np.round(Pb, 3)}")
                except Exception as ex_:
                    Log.debug(f"      edge {e.Orientation()}: no pcurve ({ex_})")
                ex.Next()

    def _restore_all(self, vbackup: Dict[int, float]) -> None:
        """Rolls back the tolerances inflated in place (shared vertices and edges)."""
        self._restore_vertices(vbackup)
        for e, t in getattr(self, "_etol_backup", {}).values():
            set_tolerance(e, t)
        self._etol_backup = {}

    def _restore_vertices(self, backup: Dict[int, float]) -> None:
        for j, t in backup.items():
            set_tolerance(td_Vertex(self.topo.vmap.FindKey(j + 1)), t)

    # --- seam ----------------------------------------------------------------------
    def _make_seam(self, R: Region, sp: SurfParam, loops, new_edges, bump_vertex):
        """
        Region closed 360 degrees: needs a seam edge from vertex A of the
        first ring to vertex B of the second, at the SAME u.
        """
        topo = self.topo

        def surviving_vertices(L: Loop) -> List[int]:
            out = []
            for ch in L.chains:
                if id(ch) in new_edges:
                    out.append(ch.verts[0])
                else:
                    out.extend(ch.verts[:-1])
            return out

        L1, L2 = loops
        path = getattr(self, "_seam_path", None)
        if path is not None:
            A, B, edges = path
            seam_edges = [(topo.edges[k], fwd) for k, fwd in edges]
            return (seam_edges, A, B, True), ""
        v1 = surviving_vertices(L1)
        v2 = surviving_vertices(L2)
        if not v1 or not v2:
            return None, "rings with no vertices"
        # A: L1's first surviving vertex (forced if L1 is a closed circle);
        # B: the one in L2 with the closest u
        u1, _ = sp.uv(topo.vpos[v1])
        u2, _ = sp.uv(topo.vpos[v2])
        best = None
        for ia, ja in enumerate(v1):
            du = np.abs(np.angle(np.exp(1j * (u2 - u1[ia]))))
            ib = int(np.argmin(du))
            if best is None or du[ib] < best[0]:
                best = (float(du[ib]), ja, v2[ib])
            if id(L1.chains[0]) in new_edges and L1.chains[0].closed:
                break  # A is forced
        du, A, B = best
        r_eff = prim_radius(R.prim)
        if du * r_eff > self.max_edge_tol:
            return None, f"seam: no aligned vertex on the second ring (deviation {du * r_eff:.1e} mm)"
        uA = float(sp.uv(topo.vpos[A][None, :])[0][0])
        vA = float(sp.uv(topo.vpos[A][None, :])[1][0])
        vB = float(sp.uv(topo.vpos[B][None, :])[1][0])
        if sp.periodic_v:
            vB = vA + float(np.angle(np.exp(1j * (vB - vA))))
        curve = sp.iso_u_curve(uA)
        if curve is None:
            return None, "seam not available for this surface"
        curve = _keep(curve)
        # real distances of the vertices from the seam curve
        for j, t in ((A, vA), (B, vB)):
            q = curve.Value(t)
            bump_vertex(j, float(np.linalg.norm(topo.vpos[j] - np.array([q.X(), q.Y(), q.Z()]))))
        VA = td_Vertex(topo.vmap.FindKey(A + 1))
        VB = td_Vertex(topo.vmap.FindKey(B + 1))
        if vB > vA:
            me = _keep(BRepBuilderAPI_MakeEdge(curve, VA, VB, vA, vB))
            fwd_from_A = True
        else:
            me = _keep(BRepBuilderAPI_MakeEdge(curve, VB, VA, vB, vA))
            fwd_from_A = False
        if not me.IsDone():
            return None, "seam's MakeEdge failed"
        return ([(td_Edge(me.Edge()), fwd_from_A)], A, B, True), ""

    # --- wire ----------------------------------------------------------------------
    def _build_wires(self, F, sp: SurfParam, loops, new_edges, seam, flip: bool):
        """
        Builds the face's wires with the pcurves unwrapped along the traversal.
        Returns [(wire, max edge tolerance)] or an error string.
        """
        topo = self.topo
        bb = BRep_Builder()

        def chain_items(ch: Chain):
            """[(edge, fwd, is_new)] in the chain's direction of travel."""
            if id(ch) in new_edges:
                e, fwd_is_v1 = new_edges[id(ch)]
                return [(e, fwd_is_v1, True)]
            return [(topo.edges[k], f, False) for k, f in zip(ch.edges, ch.fwd)]

        def loop_items(L: Loop, start_vertex: Optional[int] = None):
            items = []
            for ch in L.chains:
                items.extend(chain_items(ch))
            if start_vertex is not None:
                # rotate so the first edge starts at start_vertex
                for s in range(len(items)):
                    e, f, _ = items[s]
                    v0 = te_FirstVertex(td_Edge(e)) if f else te_LastVertex(td_Edge(e))
                    if topo.vmap.FindIndex(v0) - 1 == start_vertex or v0.IsSame(topo.vmap.FindKey(start_vertex + 1)):
                        items = items[s:] + items[:s]
                        break
                else:
                    return None
            return items

        sequences = []
        if seam is None:
            for L in loops:
                it = loop_items(L)
                sequences.append(it)
        else:
            seam_edges, A, B, _ = seam
            L1, L2 = loops
            it1 = loop_items(L1, A)
            it2 = loop_items(L2, B)
            if it1 is None or it2 is None:
                return "seam: starting vertex not found on the ring"
            # [seam A->B, ring 2 from B to B, seam B->A, ring 1 from A to A]
            up = [(e, f, ("seam", k)) for k, (e, f) in enumerate(seam_edges)]
            down = [(e, not f, ("seam", k)) for k, (e, f) in reversed(list(enumerate(seam_edges)))]
            seq = up + it2 + down + it1
            sequences.append(seq)

        out = []
        for seq in sequences:
            if flip:
                seq = [(e, not f, tag) for e, f, tag in reversed(seq)]
            w = TopoDS_Wire()
            bb.MakeWire(w)
            u_prev = v_prev = None
            u_start = v_start = None
            worst = worst_own = 0.0
            seam_pc = {}
            for e, fwd, tag in seq:
                e = td_Edge(e)
                curve, t0, t1 = edge_curve(e)
                if curve is None:
                    return "edge without a 3D curve"
                c2d, (u_end, v_end), dev = make_pcurve(sp, curve, t0, t1, fwd, u_prev, v_prev)
                if c2d is None:
                    return "pcurve can't be computed"
                if u_start is None:
                    uu, vv = sp.uv(np.array([[curve.Value(t0 if fwd else t1).X(), curve.Value(t0 if fwd else t1).Y(), curve.Value(t0 if fwd else t1).Z()]]))
                    u_start, v_start = float(c2d.Value(t0 if fwd else t1).X()), float(c2d.Value(t0 if fwd else t1).Y())
                u_prev, v_prev = u_end, v_end
                t_old = float(bt_Tolerance(e))
                tol_e = max(t_old, 1.2 * dev + 1e-7)
                key = self.topo.emap.FindIndex(e)
                if key > 0 and key not in self._etol_backup:
                    self._etol_backup[key] = (e, t_old)  # shared edge: rollback
                if isinstance(tag, tuple) and tag[0] == "seam":
                    # every seam edge appears twice: one FORWARD and one
                    # REVERSED. UpdateEdge(E, C1, C2, F)'s first pcurve is
                    # the FORWARD use's, the second the REVERSED use's.
                    kk = tag[1]
                    seam_pc[(kk, bool(fwd))] = (c2d, tol_e)
                    if (kk, True) in seam_pc and (kk, False) in seam_pc:
                        c_f, t_f = seam_pc[(kk, True)]
                        c_r, t_r = seam_pc[(kk, False)]
                        bb.UpdateEdge(e, c_f, c_r, F, max(t_f, t_r))
                        worst = max(worst, t_f, t_r)
                        worst_own = max(worst_own, t_f, t_r)
                else:
                    bb.UpdateEdge(e, c2d, F, tol_e)
                    # ⚠️ what matters for the ceiling is the deviation WE
                    # ADD: if the edge already arrives with a high
                    # tolerance, the neighbor's conversion gave it that and
                    # it's already in the model. Rejecting this region too
                    # doesn't remove it, it only removes another exact
                    # surface.
                    worst_own = max(worst_own, 1.2 * dev + 1e-7)
                    if self.verbose and tol_e > self.max_edge_tol:
                        Log.debug(
                            f"      edge beyond the ceiling: tol {tol_e:.2e} "
                            f"(pcurve deviation {dev:.2e}, pre-existing tolerance "
                            f"{t_old:.2e}, length {edge_length(e):.4f}, "
                            f"{'analytic' if edge_is_analytic(e) else 'polyline'})"
                        )
                    worst = max(worst, tol_e)
                bb.Add(w, e if fwd else td_Edge(e.Reversed()))
            # closure in (u,v) space
            if u_start is not None and u_prev is not None:
                gap = math.hypot(u_prev - u_start, v_prev - v_start)
                # deviation in parameters: the vertices lie on the surface
                # within the noise, so a small deviation is normal and the
                # tolerance covers it; a deviation of ~2pi is an unwrapping error
                r_eff = max(prim_radius(sp.prim), 1e-3)
                if gap * r_eff > 4.0 * self.max_edge_tol:
                    return f"wire not closed in parameter space (deviation {gap:.1e}, du {u_prev - u_start:+.3f} dv {v_prev - v_start:+.3f})"
            out.append((w, worst, worst_own))
        return out

    # --- validation on the solid ---------------------------------------------------
    def _normal_ok(self, Fo, rf, sp: "SurfParam") -> bool:
        """
        The new face faces the same way as the facets it replaces.

        ⚠️ NEEDED WHEN THE VOLUME DOESN'T SPEAK. The classic orientation
        check looks at the NEW edges: if we haven't created any (boundary
        entirely from the mesh) it says nothing, and on an open shell
        there isn't even the volume to act as a safety net. Here the new
        face's normal (composed with its orientation) is evaluated at the
        point of the largest facet and compared against that facet's
        normal.
        """
        topo = self.topo
        try:
            i = max(rf, key=lambda j: topo.areas[j])
            n_mesh = topo.norms[i]
            if float(n_mesh @ n_mesh) < 0.5:
                return True
            u, v = sp.uv(topo.cents[i][None, :])
            bf = BRepGProp_Face(Fo)
            P, Vn = gp_Pnt(), gp_Vec()
            bf.Normal(float(u[0]), float(v[0]), P, Vn)
            n = np.array([Vn.X(), Vn.Y(), Vn.Z()])
            nn = float(np.linalg.norm(n))
            if nn < 1e-12:
                return True
            return float((n / nn) @ n_mesh) > 0.0
        except Exception as ex:
            Log.debug(f"    normal check failed: {type(ex).__name__}: {ex}")
            return True

    def _validate(self, new_shape, Fo, rf, replaced, R: Region, a_mesh: float, sp: "SurfParam" = None):
        fe = count_free_edges(new_shape)
        if fe != self.free0:
            return False, f"free edges {fe} (was {self.free0})"
        emap = edge_face_map(new_shape)
        # consistent orientation: every new edge must be traversed in
        # opposite directions by the two faces sharing it
        check_edges = [e for _, e, _ in replaced]
        for e in check_edges:
            if not emap.Contains(e):
                return False, "new edge absent from the solid"
            faces = list(_iter_list(emap.FindFromKey(e)))
            ors = []
            for f in faces:
                ors.extend(composed_edge_orientations(f, e))
            if len(ors) == 2 and ors[0] == ors[1]:
                if self.verbose:
                    Log.debug(
                        f"    new edge: faces {[str(f.Orientation()) for f in faces]} "
                        f"edge orientations {[str(o) for o in ors]} "
                        f"Fo={Fo.Orientation()} in faces={[f.IsSame(Fo) for f in faces]}"
                    )
                return False, "new face's orientation inconsistent with its neighbors"
        # rebuilt neighboring faces valid and without spurious rings
        touched = set()
        for ch, e, _ in replaced:
            for f in _iter_list(emap.FindFromKey(e)):
                if not f.IsSame(Fo):
                    touched.add(self.registry.Add(f))
                    rec = self._nb_wires.get(ch.nb)
                    if rec is not None and count_sub(f, TopAbs_WIRE) > rec[1]:
                        return False, "neighboring face with a spurious ring (notch)"
        for k in touched:
            f = self.registry.FindKey(k)
            det = check_detail(f)
            if det:
                if self.verbose:
                    Log.debug(f"    invalid neighbor: type {face_surface_type(f)} area {face_area(f):.4f} -> {det}")
                    self._dump_face(td_Face(f))
                return False, "invalid neighboring face: " + ", ".join(det)
        # the neighbors of reused polygonal chains also share edges with an
        # updated tolerance: the new face itself has already been validated
        # direction of the new face (also works without new edges)
        if sp is not None and not self._normal_ok(Fo, rf, sp):
            return False, "new face's orientation inconsistent with its neighbors"
        # ⚠️ OPEN SHELL: no volume check. With free edges, OpenCascade's
        # "volume" is the flux of a shell that doesn't close, and it
        # changes by tens of mm3 even when replacing a tenth-sized facet:
        # it used to reject almost everything on a non-closed mesh. The
        # other checks remain (free edges unchanged, valid faces,
        # consistent area, normal direction).
        if self.free0 > 0:
            return True, ""
        # volume
        v1 = shape_volume(new_shape)
        v0 = self.vol0 if self.vol0 else shape_volume(self.shape)
        # ⚠️ THE CEILING MUST ALSO COUNT THE NEIGHBORS. Replacing a region
        # doesn't just move its own face: the polygonal boundaries touching
        # it become curves, and the neighboring faces get rebuilt along
        # with them. With the ceiling computed on the region's area alone,
        # an R0.5 spherical cap of half a mm2 wedged between three faces of
        # 1, 5 and 12 mm2 overshot by a hair (+0.0555 against 0.0513) and
        # stayed tessellated - while two identical caps, elsewhere, passed.
        # On test8 the volume check was rejecting 66 conversions and NONE
        # of them by more than four times the ceiling: it wasn't catching
        # real errors anymore, only good conversions. The new ceiling is a
        # TRUE bound, not an estimate: if no face moves more than dev, the
        # enclosed volume can't change by more than the touched area times
        # dev. And it never tightens the previous one.
        dev = max(R.sag, self.tol_fit)
        a_touch = face_area(Fo) + sum(face_area(self.registry.FindKey(k)) for k in touched)
        bound = max(10.0 * a_mesh, a_touch) * dev + 1e-9 * abs(v0) + 2e-3
        if abs(v1 - v0) > bound:
            if self.verbose:
                Log.debug(
                    f"    volume: new face area {face_area(Fo):.4f} (mesh {a_mesh:.4f}); "
                    f"touched neighbors: " + ", ".join(f"{face_surface_type(self.registry.FindKey(k))!s:.12} {face_area(self.registry.FindKey(k)):.4f}" for k in touched)
                )
                for ch, e, _ in replaced:
                    Log.debug(
                        f"      chain -> neighbor area before {self.topo.areas[ch.nb]:.4f} "
                        f"nverts {self.topo.nverts[ch.nb]} edges {len(ch.edges)} "
                        f"closed={ch.closed} new edge length {edge_length(e):.4f}"
                    )
            return False, f"volume changed by {v1 - v0:+.4f} mm3 (limit {bound:.4f})"
        self.vol0 = v1
        self.free0 = fe
        return True, ""

    def _pcurve_on_neighbor(self, e, nb: int, prim: "Prim") -> str:
        """Pcurve of the new edge on the curved analytic neighboring face (in place)."""
        f = self.topo.faces[nb]
        try:
            spn = SurfParam(prim, self._frame_ref_of_face(f, prim), self._pole_axis_of_face(f))
            # ⚠️ THE PLANE HAS ITS OWN FRAME. For cylinders and cones the X
            # axis was enough, but a Geom_Plane also has an ORIGIN, and the
            # fitted primitive's has nothing to do with the already-built
            # face's: the pcurve would end up translated by millimeters and
            # the neighbor's wire would no longer close.
            if prim.kind == PLANE:
                ad = BRepAdaptor_Surface(f, True)
                if ad.GetType() != GeomAbs_Plane:
                    return "the neighbor isn't a plane"
                pos = ad.Plane().Position()
                o, dx, dy, dz = pos.Location(), pos.XDirection(), pos.YDirection(), pos.Direction()
                spn.C = np.array([o.X(), o.Y(), o.Z()])
                spn.X = np.array([dx.X(), dx.Y(), dx.Z()])
                spn.Y = np.array([dy.X(), dy.Y(), dy.Z()])
                spn.Z = np.array([dz.X(), dz.Y(), dz.Z()])
            curve, t0, t1 = edge_curve(e)
            c2d, _, dev = make_pcurve(spn, curve, t0, t1, True, None, None)
            if c2d is None:
                return "pcurve on the neighbor can't be computed"
            tol = max(float(bt_Tolerance(e)), 1.2 * dev + 1e-7)
            if tol > self.max_edge_tol:
                return f"pcurve on the neighbor: deviation {dev:.1e} beyond the ceiling"
            BRep_Builder().UpdateEdge(e, c2d, f, tol)
        except Exception as ex:
            return f"pcurve on the neighbor: {type(ex).__name__}: {ex}"
        return ""

    def _pole_axis_of_face(self, f):
        ad = BRepAdaptor_Surface(f, True)
        try:
            if ad.GetType() == GeomAbs_Sphere:
                d = ad.Sphere().Position().Direction()
                return np.array([d.X(), d.Y(), d.Z()])
        except Exception:
            pass
        return None

    def _frame_ref_of_face(self, f, prim):
        """Frame of the already-built surface: X from the face's Geom_Surface."""
        ad = BRepAdaptor_Surface(f, True)
        t = ad.GetType()
        try:
            if t == GeomAbs_Cylinder:
                pos = ad.Cylinder().Position()
            elif t == GeomAbs_Cone:
                pos = ad.Cone().Position()
            elif t == GeomAbs_Sphere:
                pos = ad.Sphere().Position()
            elif t == GeomAbs_Torus:
                pos = ad.Torus().Position()
            else:
                return None
            d = pos.XDirection()
            x = np.array([d.X(), d.Y(), d.Z()])
            # the primitive's axis might be opposite to the built surface's:
            # SurfParam derives Y = Z x X, so X and a consistent axis are
            # enough. The prim's axis is aligned to the face's.
            dz = pos.Direction()
            z = np.array([dz.X(), dz.Y(), dz.Z()])
            if prim.kind != SPHERE and prim.axis is not None and float(prim.axis @ z) < 0:
                prim.axis = -prim.axis
                if prim.kind == AXIAL:
                    prim.slope = -prim.slope
            return x
        except Exception:
            return None

    def _carry_registry(self, rs, new_shape) -> None:
        """Analytic faces rebuilt by ReShape keep their primitive."""
        upd = {}
        for k, prim in list(self.analytic.items()):
            old = self.registry.FindKey(k)
            try:
                if rs.IsRecorded(old):
                    new = rs.Value(old)
                    if new is not None and not new.IsNull() and not new.IsSame(old):
                        upd[self.registry.Add(new)] = prim
            except Exception:
                pass
        self.analytic.update(upd)


# =============================================================================
# 11. PHASES B and C
# =============================================================================


@dataclass
class PhaseResult:
    shape: object
    regions: List[Region] = field(default_factory=list)
    n_ok: int = 0
    n_fail: int = 0
    faces_before: int = 0
    faces_after: int = 0
    free_before: int = 0
    free_after: int = 0
    vol_before: float = 0.0
    vol_after: float = 0.0
    max_tol: float = 0.0
    seconds: float = 0.0


def worst_tolerance_entity(shape) -> str:
    worst = (0.0, "")
    for v in explore(shape, TopAbs_VERTEX):
        t = float(bt_Tolerance(td_Vertex(v)))
        if t > worst[0]:
            worst = (t, f"vertex {np.round(vpos(v), 3)}")
    for e in explore(shape, TopAbs_EDGE):
        t = float(bt_Tolerance(td_Edge(e)))
        if t > worst[0]:
            worst = (t, f"edge from {np.round(vpos(te_FirstVertex(td_Edge(e))), 3)} ({str(BRepAdaptor_Curve(td_Edge(e)).GetType()).split('_')[-1]})")
    for f in explore(shape, TopAbs_FACE):
        t = float(bt_Tolerance(td_Face(f)))
        if t > worst[0]:
            worst = (t, f"face {str(face_surface_type(f)).split('_')[-1]} area {face_area(f):.3f}")
    return f"{worst[0]:.2e} mm on {worst[1]}"


def _print_regions(regions: List[Region], title: str) -> None:
    if not regions:
        return
    Log.info(title)
    hdr = f"   {'#':>3}  {'shape':<6} {'dimensions':<22} {'type':<8} {'ang.':>5} {'faces':>11}  {'RMS':>8} {'max':>8}  outcome"
    _out(hdr)
    for k, R in enumerate(regions):
        _out(f"   {k:>3}  {R.label()}  {R.rms:>8.1e} {R.max_res:>8.1e}  {R.status}" + (f"  [{R.note}]" if R.note and R.status == "OK" else ""))


def _prim_signature(p: "Prim"):
    """Primitive's identity regardless of where it sits."""
    if p.kind == AXIAL:
        return (AXIAL, round(float(p.r0), 3), round(float(p.slope), 3))
    if p.kind == TORUS:
        return (TORUS, round(float(p.r0), 3), round(float(p.r1), 3))
    return (p.kind, round(prim_radius(p), 3))


# ---------------------------------------------------------------------------
# 10b. POLYLINES -> ARCS (the CAD curves that had stayed polygonal)
# ---------------------------------------------------------------------------
#
# When a region is converted, the engine only rebuilds with the exact curve
# the boundaries for which it CAN compute the intersection between the two
# surfaces; everything else is left as it was ("reused"). On a real part
# those are the majority: on test8, 4,862 boundaries reused against 328
# rebuilt. So two already-fine analytic faces - a cylinder and the plane it
# emerges from - still touch along a six-segment polyline.
# Here the FINISHED B-Rep is examined and the question is asked: does this
# chain of segments sit on a circle? Measured on test8: 110 of 249 chains sit
# within a micron, and the radii are the ones from the drawing - R0.5 fifty-
# three times, R0.3 eighteen, R0.8 eight. It's not an approximation being
# allowed: it's the CAD's real edge, which the tessellation had split up.


def _is_line_edge(e) -> bool:
    try:
        return BRepAdaptor_Curve(td_Edge(e)).GetType() == GeomAbs_Line
    except Exception:
        return False


class _ArcSurf:
    """Minimal SurfParam over any Geom_Surface, for make_pcurve."""

    def __init__(self, surf):
        self.s = surf
        self.periodic_u = bool(surf.IsUPeriodic())
        self.periodic_v = bool(surf.IsVPeriodic())
        self._pr = _m("GeomAPI").GeomAPI_ProjectPointOnSurf()

    def uv(self, P):
        P = np.atleast_2d(P)
        u = np.empty(len(P))
        v = np.empty(len(P))
        for i, q in enumerate(P):
            self._pr.Init(_mk_pnt(q), self.s)
            if self._pr.NbPoints() < 1:
                u[i] = v[i] = np.nan
            else:
                u[i], v[i] = self._pr.LowerDistanceParameters()
        return u, v

    def point(self, u, v):
        u = np.atleast_1d(np.asarray(u, float))
        v = np.atleast_1d(np.asarray(v, float))
        out = np.empty((len(u), 3))
        for i in range(len(u)):
            q = self.s.Value(float(u[i]), float(v[i]))
            out[i] = (q.X(), q.Y(), q.Z())
        return out


def _arc_fit_circle(P: np.ndarray):
    """(center, radius, normal, X, Y, max deviation) of the circle through points P."""
    c0 = P.mean(axis=0)
    Q = P - c0
    _, _, vt = np.linalg.svd(Q, full_matrices=False)
    nrm = vt[2]
    fuori_piano = float(np.abs(Q @ nrm).max())
    # (!) RIGHT-HANDED FRAME, mandatory. The SVD doesn't guarantee vt[0] x
    # vt[1] equals vt[2]: it comes out left-handed about half the time.
    # gp_Ax2(P, N, X), instead, measures the angle from X toward N x X, so
    # if vt[1] is used as the Y axis to compute the angles, half the time
    # the parameters passed to MakeEdge are the mirrored arc's and the arc
    # gets rejected.
    X = vt[0]
    Y = np.cross(nrm, X)
    Y = Y / max(float(np.linalg.norm(Y)), 1e-300)
    xy = np.c_[Q @ X, Q @ Y]
    A = np.c_[2 * xy, np.ones(len(xy))]
    b = (xy ** 2).sum(axis=1)
    try:
        sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    except np.linalg.LinAlgError:
        return None
    cen = sol[:2]
    r2 = sol[2] + cen @ cen
    if not np.isfinite(r2) or r2 <= 0:
        return None
    rad = math.sqrt(r2)
    dev = float(np.abs(np.linalg.norm(xy - cen, axis=1) - rad).max())
    return c0 + cen[0] * X + cen[1] * Y, rad, nrm, X, Y, max(dev, fuori_piano)


def _arc_chains(topo: "Topo"):
    """
    Maximal chains of consecutive STRAIGHT edges replaceable by a single
    edge: interior vertices must have degree 2 (no other face touches
    there) and the two side faces must always be the same.
    Returns [(edges, vertices in order, [face, face]), ...].
    """
    nE = len(topo.edges)
    lin = [_is_line_edge(topo.edges[k]) for k in range(nE)]
    ef = {}
    for i in range(topo.nF):
        for k in topo.f_edges[i]:
            ef.setdefault(k, []).append(i)
    ve = {}
    for k in range(nE):
        a, b = topo.e_verts[k]
        ve.setdefault(a, []).append(k)
        ve.setdefault(b, []).append(k)
    usato = set()
    out = []
    for v0 in list(ve):
        if len(ve[v0]) == 2:
            continue                       # interior vertex: not an endpoint
        for e0 in ve[v0]:
            if e0 in usato or not lin[e0]:
                continue
            facce = frozenset(ef.get(e0, ()))
            if len(facce) != 2:
                continue
            cat = [e0]
            usato.add(e0)
            a, b = topo.e_verts[e0]
            v = b if a == v0 else a
            seq = [v0, v]
            while len(ve[v]) == 2:
                nxt = [k for k in ve[v] if k != cat[-1]]
                if not nxt or nxt[0] in usato or not lin[nxt[0]]:
                    break
                if frozenset(ef.get(nxt[0], ())) != facce:
                    break
                k = nxt[0]
                cat.append(k)
                usato.add(k)
                a, b = topo.e_verts[k]
                v = b if a == v else a
                seq.append(v)
            if len(cat) >= 3:
                out.append((cat, seq, sorted(facce)))
    return out


def _arc_wire_anchor(F_new, E_new):
    """(u_prev, v_prev, traversed forward) for E_new in the face's wire."""
    cos_ = _st(BRep_Tool, "CurveOnSurface")
    for w in explore(F_new, TopAbs_WIRE):
        ex = _BRepTools.BRepTools_WireExplorer(_TopoDS.TopoDS.Wire_s(w), F_new)
        seq = []
        while ex.More():
            seq.append(td_Edge(ex.Current()))
            ex.Next()
        for i, e in enumerate(seq):
            if not e.IsSame(E_new):
                continue
            fwd = e.Orientation() == TopAbs_FORWARD
            if len(seq) == 1:
                return None, None, fwd
            pre = seq[i - 1]
            try:
                c2d = cos_(pre, F_new, 0.0, 0.0)
                a, b = bt_Range(pre)
                q = c2d.Value(b if pre.Orientation() == TopAbs_FORWARD else a)
                return float(q.X()), float(q.Y()), fwd
            except Exception:
                return None, None, fwd
    return None, None, None


def snap_arcs(shape, tol: float, title: str = "Arcs"):
    """
    Replaces with an exact circular arc every chain of segments that truly
    sits on a circle. Every replacement is checked on BOTH FACES touching
    it before being accepted, and at the end the whole part is rechecked:
    if something doesn't hold up, all of them are given up at once.
    """
    t0 = time.perf_counter()
    topo = Topo(shape)
    cands = _arc_chains(topo)
    if not cands:
        Log.info(f"{title}: no replaceable polyline   "
                 f"[{time.perf_counter() - t0:.2f}s]")
        return shape, 0
    bb = BRep_Builder()
    rs = BRepTools_ReShape()
    fatti = 0
    for cat, seq, facce in cands:
        P = topo.vpos[np.array(seq)]
        dr = P[-1] - P[0]
        L = float(np.linalg.norm(dr))
        if L < 1e-12:
            continue                       # chain closed on itself
        d0 = P - P[0]
        dr = dr / L
        dev_retta = float(np.linalg.norm(d0 - np.outer(d0 @ dr, dr), axis=1).max())
        fit = _arc_fit_circle(P)
        if fit is None:
            continue
        cen, rad, nrm, X, Y, dev = fit
        # (!) the circle has to explain the chain MUCH better than the
        # line: an almost-straight polyline fits any old circle and it
        # doesn't mean anything.
        if not (dev <= tol and dev < 0.2 * dev_retta and rad < 1e4):
            continue
        ang = np.unwrap(np.arctan2((P - cen) @ Y, (P - cen) @ X))
        d_ang = np.diff(ang)
        if not (np.all(d_ang > 0) or np.all(d_ang < 0)):
            continue                       # doesn't always turn the same way
        if abs(ang[-1] - ang[0]) >= 2 * math.pi - 1e-9:
            continue
        # (!) THE DIRECTION. The new edge takes the place of the chain's
        # FIRST one and inherits its orientation flag in the two wires: its
        # natural direction has to match. And the arc must be taken
        # between the RIGHT two parameters: with the vertices swapped and
        # the circle left as is, you get the COMPLEMENTARY arc - 300
        # degrees instead of 60 - which in (u,v) pokes out on the other
        # side of the surface and produces "UnorientableShape".
        va, _vb = topo.e_verts[cat[0]]
        testa = va == seq[0]
        i1, i2 = (seq[0], seq[-1]) if testa else (seq[-1], seq[0])
        p1, p2 = (ang[0], ang[-1]) if testa else (ang[-1], ang[0])
        if p2 < p1:
            nrm, p1, p2 = -nrm, -p1, -p2
        ax2 = gp_Ax2(_mk_pnt(cen), _mk_dir(nrm), _mk_dir(X))
        circ = _keep(Geom_Circle(ax2, float(rad)))
        V1 = td_Vertex(topo.vmap.FindKey(i1 + 1))
        V2 = td_Vertex(topo.vmap.FindKey(i2 + 1))
        me = _keep(BRepBuilderAPI_MakeEdge(circ, V1, V2, float(p1), float(p2)))
        if not me.IsDone():
            continue
        E = td_Edge(me.Edge())
        rl = BRepTools_ReShape()
        rl.Replace(td_Edge(topo.edges[cat[0]].Oriented(TopAbs_FORWARD)),
                   td_Edge(E.Oriented(TopAbs_FORWARD)))
        for k in cat[1:]:
            rl.Remove(td_Edge(topo.edges[k].Oriented(TopAbs_FORWARD)))
        buono = True
        for fi in facce:
            F_new = td_Face(rl.Apply(topo.faces[fi]))
            up, vp, fwd = _arc_wire_anchor(F_new, E)
            if fwd is None:
                buono = False
                break
            try:
                sp = _ArcSurf(_st(BRep_Tool, "Surface")(topo.faces[fi]))
                c2d, _, devp = make_pcurve(sp, circ, float(p1), float(p2), fwd, up, vp)
            except Exception:
                buono = False
                break
            if c2d is None or not np.isfinite(devp) or devp > tol:
                buono = False
                break
            bb.UpdateEdge(E, c2d, topo.faces[fi], float(tol))
        if buono:
            for fi in facce:
                if not is_valid(td_Face(rl.Apply(topo.faces[fi]))):
                    buono = False
                    break
        if not buono:
            continue
        rs.Replace(td_Edge(topo.edges[cat[0]].Oriented(TopAbs_FORWARD)),
                   td_Edge(E.Oriented(TopAbs_FORWARD)))
        for k in cat[1:]:
            rs.Remove(td_Edge(topo.edges[k].Oriented(TopAbs_FORWARD)))
        fatti += 1
    if not fatti:
        Log.info(f"{title}: no polyline to promote "
                 f"({len(cands)} chains examined)   [{time.perf_counter() - t0:.2f}s]")
        return shape, 0
    nuovo = rs.Apply(shape)
    # (!) and now the usual check: if the whole part doesn't hold up,
    # EVERYTHING is given up. A prettier edge isn't worth a broken solid.
    fe0, fe1 = count_free_edges(shape), count_free_edges(nuovo)
    v0, v1 = shape_volume(shape), shape_volume(nuovo)
    male = ""
    if fe1 > fe0:
        male = f"free edges {fe0} -> {fe1}"
    elif not is_valid(nuovo):
        male = "BRepCheck: " + ", ".join(check_detail(nuovo, 3))
    elif abs(v1 - v0) > max(1e-4 * abs(v0), 1e-9):
        # (!) the volume DOES change, and rightly so: the polyline was
        # cutting inside the arc, the true arc lies outside by a sag. On
        # test8 that's 0.05 mm3 out of 2,754, i.e. 2 parts in a hundred
        # thousand - a twentieth of what Phase C allows itself. The ceiling
        # is only there to catch disasters.
        male = f"volume {v0:.4f} -> {v1:.4f}"
    if male:
        Log.warn(f"{title}: {fatti} arcs rejected all together ({male})")
        return shape, 0
    n0, n1 = len(topo.edges), len(Topo(nuovo).edges)
    Log.ok(f"{title}: {fatti:,} polylines promoted to exact arcs · "
           f"edges {n0:,} -> {n1:,} · volume {100 * (v1 - v0) / max(abs(v0), 1e-12):+.4f}%"
           f"   [{time.perf_counter() - t0:.2f}s]")
    return nuovo, fatti


def run_phase(
    shape,
    which: str,
    tol: Optional[float],
    max_edge_tol: Optional[float],
    min_faces: int,
    allow_sphere: bool,
    allow_cone: bool,
    allow_torus: bool,
    allow_free: bool = True,
    validate: bool = False,
    threads: int = 1,
) -> PhaseResult:
    """which = 'B' (holes only) or 'C' (all curved features)."""
    is_b = which == "B"
    Log.banner("PHASE B — circular holes" if is_b else "PHASE C — fillets, chamfers, countersinks, spheres, bosses")
    t_all = time.perf_counter()
    res = PhaseResult(shape=shape)
    res.faces_before = count_sub(shape, TopAbs_FACE)
    res.free_before = count_free_edges(shape)
    res.vol_before = shape_volume(shape)

    topo = Topo(shape)
    diag = float(np.linalg.norm(topo.vpos.max(axis=0) - topo.vpos.min(axis=0)))
    tol_fit = tol if (tol and tol > 0) else max(2e-4, 1e-5 * diag)
    tol_grow = min(10.0 * tol_fit, 1e-3 * diag)
    if max_edge_tol is None or max_edge_tol <= 0:
        max_edge_tol = max(20.0 * tol_fit, 2e-4 * diag)
    Log.info(f"Faces {topo.nF:,} · diagonal {diag:.1f} mm · vertex-surface tolerance {tol_fit:.1e} mm · edge tolerance ceiling {max_edge_tol:.1e} mm")

    t0 = time.perf_counter()
    regions = segment_curved(
        topo, tol_fit, tol_grow, diag, min_faces=min_faces, allow_sphere=allow_sphere, allow_cone=allow_cone, allow_torus=allow_torus, allow_free=allow_free, only_cyl=is_b, threads=threads
    )
    for R in regions:
        R.fobjs = [topo.faces[i] for i in R.faces]
    # ⚠️ REPETITION = FEATURE. A weird radius that shows up ONCE among four
    # facets is noise from a fillet zone; the same radius showing up twelve
    # times across the part is a repeated feature (a mill pass repeated,
    # a series of grooves) and should be rebuilt even if it sits amid loose
    # facets.
    sig = defaultdict(int)
    for R in regions:
        sig[_prim_signature(R.prim)] += 1
    for R in regions:
        R.repeated = sig[_prim_signature(R.prim)] >= 3
    if is_b:
        regions = [R for R in regions if R.closed_u and R.concave]
        Log.info(f"Candidate holes (cylinders closed 360 degrees, concave): {len(regions)}")
    Log.debug(f"Segmentation in {time.perf_counter() - t0:.2f}s")
    if not regions:
        Log.warn("No region to convert.")
        res.shape = shape
        res.faces_after = res.faces_before
        res.free_after = res.free_before
        res.vol_after = res.vol_before
        res.max_tol = max_tolerance(shape)
        res.seconds = time.perf_counter() - t_all
        return res

    eng = Engine(shape, tol_fit, diag, max_edge_tol, allow_polyline=not is_b, verbose=(Log.level <= 10))
    # order: large regions first (more analytic boundary for the following ones)
    order = sorted(range(len(regions)), key=lambda k: -sum(topo.areas[i] for i in regions[k].faces))
    pending = list(order)
    # ⚠️ two passes: a region rejected only because a neighbor was still
    # tessellated can succeed once the neighbor has been converted.
    for round_ in range(2):
        again = []
        for k in pending:
            R = regions[k]
            ok, why = eng.convert(R, strict_hole=is_b)
            # ⚠️ "amid tessellated facets" depends on WHEN you look: if the
            # neighbors get converted in the meantime, the boundary is no
            # longer a field of loose facets and the region should be
            # retried.
            if not ok and round_ == 0 and ("polygonal" in why or "tolerance" in why or "seam" in why or "tessellated" in why):
                again.append(k)
        if not again:
            break
        pending = again
        Log.info(f"Second pass on {len(pending)} regions rejected for polygonal boundaries")

    res.shape = eng.shape
    res.regions = regions
    res.n_ok = sum(1 for R in regions if R.status == "OK")
    res.n_fail = len(regions) - res.n_ok
    res.faces_after = count_sub(res.shape, TopAbs_FACE)
    res.free_after = count_free_edges(res.shape)
    res.vol_after = shape_volume(res.shape)
    res.max_tol = max_tolerance(res.shape)
    res.seconds = time.perf_counter() - t_all
    if res.max_tol > max_edge_tol:
        Log.warn("Maximum tolerance beyond the ceiling: " + worst_tolerance_entity(res.shape))
    _print_regions(regions, "Regions:")
    Log.ok(
        f"Converted {res.n_ok}/{len(regions)} regions · faces {res.faces_before:,} -> "
        f"{res.faces_after:,} · free edges {res.free_before} -> {res.free_after} · "
        f"volume {res.vol_before:.3f} -> {res.vol_after:.3f} mm3 "
        f"({100 * (res.vol_after - res.vol_before) / max(abs(res.vol_before), 1e-9):+.4f}%) · "
        f"max tolerance {res.max_tol:.1e} mm   [{res.seconds:.1f}s]"
    )
    if res.free_after != res.free_before:
        Log.error("The number of free edges changed: this should NOT happen, please report it.")
    if validate:
        ok = is_valid(res.shape)
        (Log.ok if ok else Log.warn)(f"BRepCheck: {'OK' if ok else 'NOT valid'}")
    return res


# =============================================================================
# 12. REPORT
# =============================================================================


def write_report(path: str, src: str, before, after_a, results: List[PhaseResult]) -> None:
    L = [
        "=" * 78,
        "refit.py — REPORT",
        "=" * 78,
        f"Source   : {src}",
        f"Date     : {time.strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "--- PHASE A ---------------------------------------------------------",
        f"  faces    : {before['faces']:,}  ->  {after_a['faces']:,}",
        f"  edges    : {before['edges']:,}  ->  {after_a['edges']:,}",
        f"  vertices : {before['verts']:,}  ->  {after_a['verts']:,}",
    ]
    for name, res in results:
        L += [
            "",
            f"--- PHASE {name} " + "-" * (65 - len(name)),
            f"  regions converted  : {res.n_ok} / {len(res.regions)}",
            f"  faces              : {res.faces_before:,} -> {res.faces_after:,}",
            f"  free edges         : {res.free_before} -> {res.free_after}",
            f"  volume             : {res.vol_before:.4f} -> {res.vol_after:.4f} mm3",
            f"  max tolerance      : {res.max_tol:.2e} mm",
            f"  time               : {res.seconds:.1f} s",
            "",
        ]
        for k, R in enumerate(res.regions):
            L.append(f"  {k:>3}  {R.label()}  rms {R.rms:.1e}  max {R.max_res:.1e}  {R.status}" + (f"  [{R.note}]" if R.note and R.status == "OK" else ""))
    L.append("")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
    Log.ok(f"Report: {path}")


# =============================================================================
# 13. CLI
# =============================================================================


def default_threads() -> int:
    """Two cores left free for the system, never more than 22: beyond that
    nothing more is gained because there are only a few dozen seeds to try."""
    n = os.cpu_count() or 2
    return max(1, min(22, n - 2))


class _PhaseArg(argparse.Action):
    """
    ⚠️ PHASES RUN IN THE ORDER THEY'RE WRITTEN. argparse only gives the
    final values, not the order: this action records every -a/-b/-c into a
    list as it appears. This way '-c' alone means ONLY Phase C (useful for
    figuring out whether a defect comes from an earlier phase), and
    '-a -b -c -a 0.005' means exactly what's written.
    """

    def __call__(self, parser, ns, values, option_string=None):
        seq = list(getattr(ns, "sequence", None) or [])
        # con nargs=0 argparse passa una lista vuota, non None
        val = None if isinstance(values, (list, tuple)) else values
        seq.append((option_string.lstrip("-")[0].upper(), val))
        ns.sequence = seq


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="refit.py",
        description=(
            "Soft reconversion of a mesh (STL, or STEP born from a mesh) into an analytic B-Rep.\n"
            "\n"
            "  Phase A (-a)  merging of coplanar faces and collinear edges.\n"
            "  Phase B (-b)  through or blind circular holes: a cylinder closed 360\n"
            "                degrees, concave, opening onto two planes orthogonal to\n"
            "                the axis with two circular boundaries. Only those: it's\n"
            "                the cautious phase.\n"
            "  Phase C (-c)  fillets, chamfers, countersinks, counterbores, corner\n"
            "                spheres, bosses: cylinders, cones, spheres and tori (even\n"
            "                partial), with boundaries left polygonal where no exact\n"
            "                curve exists.\n"
            "\n"
            "Every region is replaced on its own and checked immediately (solid still\n"
            "closed, faces valid, consistent orientation, area and volume consistent\n"
            "with the mesh): if a check fails it's rolled back and that piece stays\n"
            "tessellated. The output is always consistent with the input, at worst\n"
            "it's less clean. Without -a/-b/-c all three phases run."
        ),
        epilog=(
            "examples\n"
            "  python refit.py part.stl                  all phases\n"
            "  python refit.py part.stl -a               only merging coplanar faces\n"
            "  python refit.py part.stl -a 0.01          same, merging up to 10 microns\n"
            "  python refit.py part.stl -b -c --report   holes and features, with a .txt report\n"
            "  python refit.py part.stl -b -c -a 0.01    final merge at 10 microns\n"
            "  python refit.py part.stl -b -c -i 2       two B/C cycles before saving\n"
            "  python refit.py part.step -o out.step     STEP input\n"
            "  python refit.py part.stl -j 1             everything sequential (reproducible)\n"
            "\n"
            "report notes (--report), one line per region\n"
            "  [N analytic . N reused . N polygonal]  boundaries rebuilt with the exact\n"
            "      curve; boundaries already good taken from the mesh; boundaries left as a polyline.\n"
            "  boundaries with the planes left polygonal  the face is exact, the\n"
            "      boundaries with the neighboring planes stay the mesh's (second attempt).\n"
            "  contour left polygonal                no new edge at all: the contour\n"
            "      stays identical to the mesh (third attempt, to avoid breaking an\n"
            "      already-closed analytic neighbor).\n"
            "  rejected: ...                         why the region stays tessellated.\n"
            "      'region too small amid tessellated facets' means a primitive fitted\n"
            "      on a fillet zone's noise: better to leave the mesh and finish by hand.\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("input", nargs="?", help=".stl or .step/.stp file")
    p.add_argument("-o", "--output", metavar="FILE", help="output STEP file (default: <name>.step, or <name>_phaseA.step with Phase A alone)")
    p.add_argument("--check", action="store_true", help="check the environment (OpenCascade, numpy, required functions) and exit")

    g = p.add_argument_group("phases (none = A B C A)")
    g.add_argument(
        "-a",
        nargs="?",
        const=-1.0,
        type=float,
        default=None,
        dest="ph_a",
        metavar="TOL",
        action=_PhaseArg,
        help="Phase A: merge coplanar faces. Optional value = max distance of "
        "the vertices from the common plane, in mm (e.g. -a 0.01 merges "
        "faces with deviations up to 10 microns). Without a value: 2e-6 x "
        "diagonal. Can be repeated: 'refit.py x.stl -a -b -c -a 0.005' "
        "does a tight A, B, C and a wide final re-merge",
    )
    g.add_argument("-b", nargs=0, dest="ph_b", action=_PhaseArg, help="Phase B: only through/blind circular holes")
    g.add_argument("-c", nargs=0, dest="ph_c", action=_PhaseArg, help="Phase C: fillets, chamfers, countersinks, spheres, tori, free-form shapes")
    g.add_argument("-i", "--iterations", type=int, default=1, metavar="N", help="repeat the whole phase sequence N times before saving (default 1)")
    g.add_argument("--keep-a", action="store_true", help="also save the intermediate <name>_phaseA.step after the first Phase A")

    g = p.add_argument_group("Phase A")
    g.add_argument(
        "--lin-tol",
        type=float,
        default=None,
        help="max distance of the vertices from the common plane to merge two "
        "faces (mm): this is the INITIAL Phase A's tolerance, the one -a "
        "only changes for the final pass. Default: 2e-6 x diagonal, "
        "minimum 1e-5",
    )
    g.add_argument("--no-cad-edges", action="store_true", help="[A] don't protect the CAD edges recognized by the mesh's texture change (see the note in the source)")
    g.add_argument("--ang-tol", type=float, default=0.005, help="angular tolerance for the classic merge of non-triangular STEP (degrees)")

    g = p.add_argument_group("Phases B and C")
    g.add_argument(
        "--tol",
        type=float,
        default=None,
        help="max distance between a mesh vertex and the surface for the "
        "primitive to be accepted (mm). It's the ONLY tolerance of phases "
        "B and C: all the others derive from it - region growth 10x (and "
        "never beyond 1e-3 x diagonal), intersection curves 4x, "
        "'is it the same surface as the neighbor?' 50x, edge-tolerance "
        "ceiling 20x (see --max-edge-tol). SMALLER = the part stays "
        "closer to the mesh and regions don't overrun past the CAD's "
        "edges, but what the mesh doesn't support stays tessellated. "
        "LARGER = more faces recognized, but the region widens past the "
        "true edge and the part deforms. Does NOT touch Phase A, which "
        "has --lin-tol and the value after -a. Default: 1e-5 x diagonal, "
        "minimum 2e-4 (on a part with a 150 mm diagonal: 1.5e-3)",
    )
    g.add_argument(
        "--max-edge-tol",
        type=float,
        default=None,
        help="edge tolerance ceiling: a polygonal boundary left on an analytic "
        "face that deviates more than this makes the region get rejected "
        "(mm). Default: 20 x --tol, and in any case at least 2e-4 x "
        "diagonal",
    )
    g.add_argument(
        "--no-arcs",
        action="store_true",
        help="don't promote to an exact arc the polylines that sit on a circle "
        "(it's the last step, after all phases: on test8 it removes 288 "
        "edges without moving a single face)",
    )
    g.add_argument(
        "--arc-tol",
        type=float,
        default=None,
        help="how far a polyline's vertices can be from the circle for it to "
        "become an arc (mm). Default: --tol, and failing that 1e-5 x diagonal",
    )
    g.add_argument(
        "--min-faces",
        type=int,
        default=4,
        help="minimum facets for a group to become a region (default 4). "
        "Besides this, Phase C rejects regions under 12 facets that have "
        "more than 60%% of their boundary resting on other tessellated "
        "facets: those are primitives fitted on noise",
    )
    g.add_argument("--no-sphere", action="store_true", help="[C] don't recognize spheres")
    g.add_argument("--no-cone", action="store_true", help="[C] don't recognize cones")
    g.add_argument("--no-torus", action="store_true", help="[C] don't recognize tori")
    g.add_argument("--no-free", action="store_true", help="[C] don't rebuild free-form fillet patches (B-spline): they stay tessellated or split into little cylinders")

    g = p.add_argument_group("system / output")
    g.add_argument(
        "-j",
        "--threads",
        type=int,
        default=default_threads(),
        metavar="N",
        help="ceiling on worker processes for the primitive search in "
        f"Phases B and C (default {default_threads()} on this machine, 1 = "
        "everything sequential). How many are actually used depends on "
        "the work at hand. Phase A and the in-solid replacements stay on "
        "a single core. With more processes the regions found can change "
        "slightly: with -j 1 the result is reproducible",
    )
    g.add_argument("--validate", action="store_true", help="full BRepCheck at the end of every phase (slow on large shapes)")
    g.add_argument("--report", action="store_true", help="write <name>_report.txt with the outcome of every region")
    g.add_argument("--log", metavar="FILE", help="write the log to a file")
    g.add_argument("-v", "--verbose", action="store_true", help="per-region detail: curves tried, attempts, rejection reasons")
    g.add_argument("--quiet", action="store_true", help="only warnings and errors")
    g.add_argument("--no-color", action="store_true", help="no ANSI colors")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    _enable_vt_windows()
    try:
        sys.stdout.reconfigure(errors="replace")
    except Exception:
        pass
    Log.no_color = args.no_color or not sys.stdout.isatty()
    Log.level = 10 if args.verbose else (30 if args.quiet else 20)

    if args.check:
        Log.banner("Environment check")
        Log.ok(f"OpenCascade  : {_NS}")
        try:
            import OCP  # noqa

            Log.ok(f"cadquery-ocp : {importlib.metadata.version('cadquery-ocp')}")
        except Exception:
            pass
        Log.ok(f"numpy        : {np.__version__}")
        Log.ok(f"python       : {sys.version.split()[0]} ({sys.platform})")
        for name, fn in (
            ("BRep_Tool.Pnt", bt_Pnt),
            ("TopExp.FirstVertex", te_FirstVertex),
            ("BRepTools_ReShape", BRepTools_ReShape),
            ("MakeShapeOnMesh", getattr(_BRepBuilderAPI, "BRepBuilderAPI_MakeShapeOnMesh", None)),
            ("UnifySameDomain.KeepShape", getattr(ShapeUpgrade_UnifySameDomain, "KeepShape", None)),
        ):
            (Log.ok if fn is not None else Log.error)(f"{name:<28}: {'ok' if fn is not None else 'MISSING'}")
        return 0

    if not args.input:
        Log.error("Missing input file (.stl or .step). Use -h for help.")
        return 2

    # phase sequence, in the order written. No flags: A B C A.
    seq = list(getattr(args, "sequence", None) or [])
    if not seq:
        seq = [("A", None), ("B", None), ("C", None), ("A", None)]
    stem, _ = os.path.splitext(args.input)
    only_a = all(w == "A" for w, _ in seq)
    out_path = args.output or (f"{stem}_phaseA.step" if only_a else f"{stem}.step")
    t_start = time.perf_counter()

    shape, before = read_input(args.input)
    # the Phase B/C worker processes are started right away: they take
    # about a second to spin up, and they do it while the first phase runs
    if (not only_a) and args.threads > 1 and before["faces"] >= 2000:
        warm_pool(args.threads, before["faces"])

    Log.info("Sequence: " + " ".join(w + ("" if v is None or v <= 0 else f"({v:g})") for w, v in seq))
    results: List[Tuple[str, PhaseResult]] = []
    after_a = before
    n_a = 0
    for it in range(max(1, args.iterations)):
        tag = f" (cycle {it + 1})" if args.iterations > 1 else ""
        for w, val in seq:
            if w == "A":
                n_a += 1
                lt = val if (val is not None and val > 0) else args.lin_tol
                titolo = "PHASE A — merging coplanar faces" if n_a == 1 else "PHASE A — re-merging faces"
                shape, bef_a, after_a = phase_a(shape, lt, args.ang_tol, validate=args.validate, title=titolo + tag, use_barrier=not args.no_cad_edges)
                # ⚠️ 'before' stays the FILE's count, not the first Phase A's:
                # with Phase A not first in sequence (e.g. -b -c -a) otherwise the
                # report would say "994 input faces" instead of 5,682.
                if n_a == 1:
                    if args.keep_a and not only_a:
                        write_step(shape, f"{stem}_phaseA.step")
            elif w == "B":
                rb = run_phase(shape, "B", args.tol, args.max_edge_tol, args.min_faces, True, True, True, validate=args.validate, threads=args.threads)
                shape = rb.shape
                results.append(("B" + tag, rb))
            else:
                rc = run_phase(
                    shape,
                    "C",
                    args.tol,
                    args.max_edge_tol,
                    args.min_faces,
                    not args.no_sphere,
                    not args.no_cone,
                    not args.no_torus,
                    allow_free=not args.no_free,
                    validate=args.validate,
                    threads=args.threads,
                )
                shape = rc.shape
                results.append(("C" + tag, rc))
    # (!) LAST STEP: polylines that are actually CAD arcs turn back into arcs.
    # This runs at the end because it needs the finished B-Rep: an edge becomes
    # a replaceable chain only once the faces on both its sides are in place.
    if not (only_a or args.no_arcs):
        Log.banner("Arcs · the polylines that used to be CAD curves")
        _V = Topo(shape).vpos
        _diag = float(np.linalg.norm(_V.max(axis=0) - _V.min(axis=0))) if len(_V) else 1.0
        _tolA = args.arc_tol if args.arc_tol and args.arc_tol > 0 else (
            args.tol if args.tol and args.tol > 0 else max(2e-4, 1e-5 * _diag))
        shape, _ = snap_arcs(shape, _tolA)
    Log.banner("Saving")
    # ⚠️ THE LAST WORD BEFORE WRITING. Every conversion validates itself,
    # but the checks are local: if something slipped through, this is where it shows.
    if not is_valid(shape):
        Log.warn("the finished part doesn't pass BRepCheck: " + ", ".join(check_detail(shape, 4)))
    write_step(shape, out_path)
    _ri, _ = read_step(out_path)
    if is_valid(shape) and not is_valid(_ri):
        Log.warn("valid in memory but not after the STEP write: " + ", ".join(check_detail(_ri, 4)))
    if args.report:
        write_report(f"{stem}_report.txt", args.input, before, after_a, results)
    if args.log:
        Log.dump(args.log)

    Log.banner("Done")
    st = shape_stats(shape)
    Log.ok(
        f"{before['faces']:,} input faces -> {st['faces']:,} output faces · "
        f"free edges {count_free_edges(shape)} · max tolerance {max_tolerance(shape):.1e} mm "
        f"· {time.perf_counter() - t_start:.1f}s"
    )
    for name, r in results:
        Log.ok(f"Phase {name}: {r.n_ok}/{len(r.regions)} regions converted")
    close_pool()
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\nInterrupted.")
    finally:
        close_pool()
