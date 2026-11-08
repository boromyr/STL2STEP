# How Refit works — technical notes

The full technical reference (phases, tolerances, mesh repair, output details). For the short version go back to the [README](../README.md).

"Soft" reconversion of a mesh (STL, or a STEP born from a mesh) into an analytic B-Rep using OpenCascade.

The principle is a single one: **every change to the geometry is local**, is verified right away
(solid still closed, faces valid, consistent orientation, area and volume consistent with the mesh)
and is discarded if it doesn't pass the check. Whatever can't be converted stays tessellated as it
was: the output file is always consistent with the input one, at worst it's just less "clean". There
is no global Sewing and nothing is rebuilt globally.

## Getting started

### Requirements

- Python 3.9+
- [OpenCascade bindings](https://github.com/CadQuery/OCP) (`cadquery-ocp`), `numpy`, `scipy`
- `pymeshlab` (optional but recommended: enables mesh repair, see below)

```bash
pip install cadquery-ocp numpy scipy pymeshlab
```

Works with OCP 7.x and 8.x. `cadquery-ocp` 8.0.1 is built against `vtk==9.6.2`: with a newer vtk
(e.g. after a blanket `pip` upgrade) `import OCP` fails with "DLL load failed" — pin `vtk==9.6.2` if
that happens. `scipy` is needed by the mesh check and by `tools/qa.py`; without `pymeshlab` the mesh
check still runs (on `numpy` alone) and defects are reported but left unrepaired.

### First run

```bash
python refit.py part.stl
```

This runs the default phase sequence (`A C A(0.005) C A(0.005)`, see below) and writes `part.step`
next to the input. A few of the bundled test meshes are a good first try:

```bash
python refit.py test/test0.stl
```

Then check the result against the input mesh:

```bash
python tools/qa.py test/test0.stl test/test0.step
```

and, for the test files that ship with their original CAD (`test0`, `test4`, `test10`, `test12`,
`test13`), compare face-by-face against the ground truth:

```bash
python tools/cadcmp.py test/test0_original.step test/test0.step -v
```

See the docstring at the top of [`refit.py`](../refit.py) for the full guide (every flag, how to read
the report, the final arc-snapping step).

### Viewer

```bash
pip install pywebview
python viewer/server.py              # or double-click "Refit Viewer.pyw" in the repo root
python viewer/server.py part.stl     # with part.stl already open
```

opens a desktop window (pywebview: a normal Windows window drawn by the system's WebView2, the way
Tauri apps are; three.js comes from a CDN) where you pick a mesh, set the phases and tolerances, run `refit.py` with its log streamed live, and look at
the result: input mesh, converted STEP and, when there is one, the `<name>_original.step` CAD, side
by side. Faces are coloured by surface type (the triangles still left from the mesh in yellow), the
B-Rep edges can be coloured by curve type, hovering a face shows its type, radius or angles and
area, and the panel compares face and edge counts with the original CAD. The command it runs is
shown in the page and is exactly what you'd type; the line under it warns when the output file
already exists and will be overwritten. `--browser` shows the same page in the default browser
instead (that's also the fallback without pywebview), `--serve` only runs the local server; in the
browser, files dropped on the view are copied into `viewer_uploads/`.

## Mesh check (STL input)

Before anything else the mesh is checked: triangles with a repeated vertex or zero area, duplicate
triangles, holes, non-manifold edges and vertices, inconsistent winding, self-intersections, separate
pieces (bodies, or internal cavities: a piece of negative volume inside the body becomes a void of
the solid, as in the CAD's `BREP_WITH_VOIDS`). The defects that would break the conversion are
repaired with [MeshLab](https://www.meshlab.net/) through `pymeshlab`, one filter per defect and never
by remeshing, so the vertices of a CAD export stay exact:

- zero-area triangles with a vertex on their side: flipped with the neighbor (no crack); strips of
  them (a zero-width sliver the exporter tessellated as a fan of flat triangles, which flipping
  can't take away) removed and the crack zipped, the triangles on both sides split at each
  other's vertices (`test14`: 441 free edges → 11);
- duplicate triangles: removed;
- non-manifold edges and vertices: vertices split (bodies that touch come apart, nothing moves);
- inconsistent winding: re-oriented coherently;
- holes up to `--close-holes N` boundary edges (default 10): closed; bigger ones stay open;
- self-intersections: only reported (repairing them would move the surface).

`--no-mesh-repair` only checks. Without `pymeshlab` the check still runs (numpy) and the defects
are reported but left in.

## Phases

- **Phase A** (`-a`): merging of coplanar faces and collinear edges.
- **Phase B** (`-b`): circular holes — cylindrical walls closed 360°, concave, with edges that are
  exact circles.
- **Phase C** (`-c`): fillets, chamfers, countersinks, spotfaces, corner spheres, bosses — every
  region of facets that sits on a cylinder / cone / sphere / torus gets replaced by the analytic
  surface. Where no quadric describes the surface, the region is rebuilt with a B-spline.
  A fillet along a curved edge that the tessellator cut into strips (each strip alone looks like a
  cylinder or a cone) is merged back into one torus when the torus explains the strips better
  than their own cylinders did. The walls of an extruded outline — a sketch with splines, embossed
  or engraved text — become one surface swept along the extrusion direction: a cylinder when the
  profile is a circular arc, otherwise a B-spline profile. A band that doesn't take one profile is
  cut where its curvature jumps (the end of a fillet, the corner between two curves), and a spline
  wall the seeds had already cut into a few osculating cylinders is taken whole again when one
  profile holds it (`test9` 1,062 → 954 faces, `test8` 479 → 468). Stretches tessellated so
  coarsely that the mesh doesn't say where the curve runs stay planar, and so does a face Phase A
  already made planar that meets a couple of slivers of the round along its rim: a merged polygon
  carrying twenty times the area of all the other strips together is the face, and the band is
  tried without it (`test14`'s 485 mm² side, a plane in the CAD, had become an 8×4 B-spline; any
  strip four times bigger took ten good bands off `test9`'s spline walls). Where a free-form band runs
  into a round it's tangent to, the facets along the contact line sit on both: one with more sides
  on its neighbour than on its own region, and on the neighbour's surface, goes over to it (`test4`:
  no more notches at the ends of the chamfer bands, edges 137 → 126). Every facet is a chord of
  the surface under it: a region whose facets sag below its surface further than the mesh's own
  chords ever do (measured on the regions that fit exactly) has the wrong surface - a sphere
  standing over a flat disc - and is rejected, or only the facets that don't belong are left out.
  Some exporters cap the angle between chords rather than their sag, and there a bigger radius sags
  deeper: the sag per mm of radius is measured too, on the exact regions of the largest radii
  (`test16`: 10° chords, its r 18 rounds sag 0.07 mm where its fillets sag 0.01; 1,581 → 1,170
  faces), and a region may also sag 1.5 times what the exact regions of a similar radius do, cones
  included (`test17`, from OpenSCAD, whose chords have a fixed length: 3,078 → 2,643 faces). A band swept along a direction checks its fit inside every big planar polygon it takes
  in, not just at its centre (`test16`'s 25 mm handle underside came out a wavy B-spline, 0.4 mm
  between its checks, and read back broken). A constant-radius fillet between two
  curved walls (a hole meeting another hole, whose width changes all along) is rebuilt from the
  two walls and the ball rolling between them, as one B-spline face whose edges follow the ball's
  contact lines; the straight and circular pieces of the same chain become the exact cylinders and
  tori the ball sweeps there, and each fillet ends where its ball leaves the walls' facets. A facet
  left over at a fillet's corner, wedged among faces already found, joins the fillet when its
  vertices are on it and its middle sags no deeper than a chord of the ball does; a fillet that then
  fails to convert is converted without it (`test13` reproduces all 57 faces of its CAD one to one, and has no others; `test8` 468 → 463
  faces).
  A thread is a profile swept by a screw motion, which no quadric and no extrusion describes: its
  flanks used to come out as hundreds of small cylinders, spheres and loose triangles. The screw is
  found from the facets' normals (a screw's velocity is tangent to every surface it sweeps onto
  itself, so one least-squares eigenvector gives axis, point and pitch), refined on how the
  vertices of every turn fall onto one profile, and each straight piece of that profile becomes
  one face: a B-spline swept along the helix for a flank, cylinders cut into pieces of at most
  300° for a crest or a root (`test15`: an internal thread of 4.4 turns, pitch 2.14 mm, two faces
  where there were 190 regions and 3,400 triangles, 1,431 → 795 faces; `test17`: the lid's 1/8"
  thread, 2,642 → 1,730 faces). A rounded crest or root is a circular arc of the profile, swept the
  same way; a flank takes the triangles wedged in its rim whose corners are all on its sweep (taking
  only those standing on a row left a sawtooth rim that touched itself at every vertex, and the face
  was invalid); a 300° piece of crest or root cylinder that is refused (a helical band went 354°
  round with its facets' width: "ambiguous angular coverage") gives way to its halves (`test14`'s M10 nut, pitch 0.907 mm, 42,000 triangles: on the hole cut
  out of the part, 5,515 → 3,735 faces, its flanks from none to almost all converted).
  The profile is then cut once for the whole thread, not strip by strip: in screw coordinates every
  vertex of every turn falls on one curve, the radius a periodic function of the position along the
  pitch, and the fewest straight pieces and arcs that hold it within `tol_fit` are its faces - the
  way the CAD drew it. Each stretch of the pitch is represented by its densest row of vertices (the
  turns pile up on it, strays from a run-out or another feature don't), arcs go up to a fifth of the
  pitch (a wider "arc" is a line and a rounding taken together), a line is preferred to an arc, a
  piece takes the tangent rows it shares with its neighbours, and two pieces with no row between
  them are joined by the straight piece between their ends (a flat crest has rows only on its edges).
  Cut strip by strip, `test14`'s G1/8 Whitworth thread lost half of every root to the flank's
  strip, where straight segments took it as chords; cut once, its six pieces come out as in its
  CAD (crest arcs r 0.1249, flat crest, flanks at 27.5°, root arc r 0.1245 - 0.137 of the pitch,
  Whitworth's). A piece of a turn is a band too (from 45°): the profile is proven by the whole
  thread, and this one runs into a pocket for eight turns, 0.24 to 0.83 of each left. The cut once
  is kept when it puts more facets in the thread's bands than the strips do (`test15`'s V thread:
  the same two flanks either way). On a closed 4-turn section of `test14`'s thread: 1,183 → 245
  faces, polyline edges 2,887 → 697, mean distance from the CAD 0.0049 → 0.0036 mm (maximum
  unchanged, 0.038), valid.
  A tube that is no quadric - a loft narrowing along its axis, its sections ellipses on tilted
  planes - used to come out as dozens of osculating cylinders and cones, each a strip of a few
  degrees. Strips that touch, with nearly parallel axes and similar radii, going round a common
  axis, are rebuilt as one B-spline in the tube's own coordinates (the angle around the axis and
  the height along it), every facet it keeps held within its own threshold: refitted while the
  facets it can't hold are let go, then grown over the ones around it it does hold (`test15`: the
  sharpener's guide arm, 140 strips and 5,600 facets, one face; 795 → 545 faces, deviation from
  the mesh 0.018 → 0.010 mm).
  OpenSCAD writes circles as regular polygons, every corner exactly on the circle: the flat sides
  are chords, and where another body cuts them its vertices lie on the chord, under the circle. A
  flat facet whose plane holds two generators of a cylinder or cone, with its two extreme vertices
  on it, none outside, at most 15° wide and none deeper than such a chord sags, belongs to it. A
  region whose vertices are exact and whose widest facets are several and all as wide proves its
  own sag, and may sag half as much again (tori too, around their axis and around the tube). An
  ASCII STL is exact only to its digits - OpenSCAD writes six, 21.4919 - and that is the precision
  the exact-vertex tests use there. A surface going (almost) all the way round whose boundary isn't
  two rings, or whose seam finds no vertices in line on both rims, is converted in pieces around
  its axis, and torus strips found only late are merged once more at the end (`test17`: the boss
  wall, its rounded rim, the bore and the bowl; 1,730 → 940 faces).
  Some exporters tessellate an edge shared by two faces twice, each face with its own steps, and
  make the mesh conform by dropping one side's points onto the other side's chords: those
  T-vertices lie on a straight line between two neighbours, inside the curve by the chord's sag.
  Where that line is the chord of an arc (the polyline bends the same way at both its ends) and the
  vertex stands off the arc by more than the tolerance, the fits and the vertex checks take the
  chord's far end in its place; the facet's own shape is still used to measure a free-form patch's
  sag (`test16`: the r 4 round along the handle, 55 mm strips all ending on such points, was 13
  planes and 2 cylinders of the wrong radius, now one face; 1,202 → 813 faces. `test15`: 545 → 390,
  `test17`: 940 → 931). A sphere's seam is put behind the face rather than through its middle, and
  an arc that would cross a neighbouring edge just past their shared vertex - a needle between two
  nearly tangent arcs, hidden in memory by the vertex's tolerance - is left a polyline: both came
  back invalid from the STEP file. On a curved face the pcurves are asked the same in (u, v): on
  `test18`'s window the r 1.2 round, once merged, ends in a needle between a cylinder and a patch it
  meets almost tangent, and the curve where they cross wandered over the needle's other edge 1.6 mm
  from their vertex - valid in memory and read back, but the mesher, given a wire crossing itself,
  left a corner of the round uncovered at the patch's far end, 0.50 mm off the part. Refused, the
  chain is cut and its pieces away from the needle still become curves. There a crossing inside
  the vertex's tolerance doesn't count: `test8`'s ellipses on r 0.5 rounds meet a line at both
  ends, nearly tangent, and cross 3e-4 from the vertex (tolerance 9.6e-3) - read back valid and
  meshed; refused, 7 of its 17 ellipses went.
  Where rounds of different radii meet (`test16`'s r 4 along the handle running into the r 2 along
  the scoop's rim), the CAD fills the corners between them with patches of its own: free-form,
  folding over 90°, bounded by a dozen pieces of round, their facets on none of their neighbours. A
  group of loose facets no single B-spline holds is filled the way a plate is, once everything
  else has taken its facets: flattened (onto its mean plane, or conformally when it folds), a
  coarse spline through its vertices, the rest a polyharmonic interpolation (r⁵), sampled densely
  and fitted by the sparsest B-spline keeping the vertices within the tolerance, the facets within
  what the rounds around them sag and the engine's cap on the sag, and the normals with the
  facets; a group that doesn't hold is cut in two across its length and each half tried again,
  and a piece the engine still refuses is cut again there. A free-form patch over a whole smooth
  section that the engine refuses gives way to the regions it had replaced (`test16`: 788 → 660
  faces; `test17`: 866 → 733; `test15`: 388 → 390 faces, deviation from the mesh 0.014 → 0.010 mm).
  Weak regions - not exact (a surface exported exactly sits on its vertices to a few hundredths of
  the tolerance), small or a few degrees wide, or free-form already - that meet across a smooth
  edge are merged when one such filled B-spline holds them all: the strips of osculating
  cylinders, cones and little B-splines a loft, a variable-radius round or a round rolling between
  two lofts break into become one face, and a junction where the curvature jumps stays as it was
  (a merged face the engine refuses gives way to its pieces). Where the vertices come in rows - a
  CAD tessellates a round along a gently curved edge with its cross-section cut finely and the rows
  millimetres apart, long facets between them - the interpolation doesn't follow the facets: across
  a 3.5 mm gap between two rows of `test18`'s r 1.2 window rounds it swung 0.16 mm off facets that
  sag 0.009, and two strips of the same round, each an osculating cylinder, never became one face.
  There, when no interpolation holds, the spline is fitted to the vertices and, with a twentieth of
  their weight in all, to points spread over the facets themselves, kept on the same terms
  (vertices within the tolerance, the facets' inner points within what they sag, normals with the
  facets); a patch that held already is the same. Only there: slices of one round - quadrics all,
  their section radii within a twenty-fifth of their median (a cylinder's radius, a cone's where
  its facets are, a torus' tube, a sphere's; `test18`'s within 2.2%, while within a tenth
  `test11`'s exact r 3 cylinder and a 4-facet r 3.35 piece beside it made a patch and three
  ellipses went) - whose vertices come in rows (the widest gap from a facet's
  inner point to the nearest vertex six times the vertices' spacing or more: 6 to 14 on those
  rounds, 1.5 for the median region). Everywhere, it held unions the interpolation had refused for
  good reason - a round against a flat, the curvature jumping across - and exact pieces went
  (`test8`: 452 → 384 faces, but 38 more polylines, 33 circles and 7 ellipses fewer; `test6` 4 mm³
  of volume). And the rows no farther apart than sixteen spacings: the gap between a round's rows
  grows as the square root of its spine's radius, and past that the spine is straight - a CAD's
  exact cylinders along straight or polygonal edges, their vertices on their end circles
  (`test4`'s r 3 and r 5; `test8`'s r 0.5 pieces, 19 to 95 spacings apart; `test16`'s r 1 ones, 100
  to 680), which went into patches otherwise (`test8`: 9 ellipses gone; `test16`: a sphere beside
  them refused, 801 → 837 faces). On a closed section of `test18`'s window (11,164
  triangles): 260 → 235 faces, valid, polylines 995 → 926 mm in all (between curved faces 300 →
  309, between a curved face and a plane 320 → 290), 0.050 and 0.048 mm from the mesh (0.019 and
  0.057 before, the worst now where a merged patch meets the section's cut); 9 short edges a
  quadric piece met as a circle are polylines along a patch now, 4 of them on the cut. A face
  going all the way round whose seam is a path of mesh edges, slanting, can cover a strip of its
  surface twice once its neighbours' edges are in: valid to BRepCheck, its tessellation inside the
  part; it is tessellated as it sits and, if it strays, built in pieces round the axis; a
  neighbour's new edge gets its
  pcurve on such a face by the face's own u at the edge's two vertices. A torus whose vertex rows
  round its tube have one gap far wider than the others (two rounds at the ends of a wall, read as
  one barrel through both) gets no allowance for a big radius's sag, and the facets it isn't the
  chord surface of are converted as the cylinder or cone they lie on (`test16`, Autodesk exporter:
  995 → 874 faces, deviation 0.39 → 0.058 mm; `test17`, TurboCAD: 1,195 → 1,131 faces, 0.54 →
  0.095 mm; `test18`: 1,307 → 1,212; `test9` 931 → 896, `test15` 390 → 381).
  A big polygon from Phase A has all its vertices on its rim, and lies on any surface through
  that rim: a facet of five vertices or more whose normal is 45° off the surface's at its vertices
  is no chord of it. A region whose surface misses the facets it keeps once those that aren't its
  chords are left out is fitted again as a cylinder or cone on them, the facets it doesn't hold let
  go, three times (`test17`'s flange wall, an r 50 cylinder two millimetres tall, came out of the
  seeds as the r 50 sphere through its two rings, with the foot's whole flat face in it). And every
  group of loose facets left at the end that is curved - small polygons whose vertices are off one
  plane - is filled too, held to the mesh's own chord sag, not only those of the smooth sections no
  single B-spline holds (`test17`'s eight windows blend into the bore's floor round in sections a
  plane owns: 200 loose triangles). `test17`: 1,131 → 657 faces; `test18`: 1,212 → 766; `test16`
  874 → 864; `test8` 457 → 452.
  Only the facets the patch would take can make a plane own a section: on `test17`'s slots the
  section of a shoulder (where the slot's r 1.5 round turns) held three big triangles of the r 21.5
  wall it is tangent to, the biggest 26% of the section - "a plane" - and the shoulder was left to
  tori merged late into a stepped spline, while its mirror twin, whose triangle was 25%, became
  one patch. The wall's cylinder straddles the section and the patch leaves its triangles out.
  A piece of a chord is still under the circle: a flat facet with two vertices on a cylinder and
  the others inside it, each ten tolerances deep at least and no deeper than a quarter of what the
  cylinder drops over its distance to the nearest vertex on it (L²/4r), is a chord of it - a round
  tangent to the cylinder drops far faster and keeps its rows. `test17`'s 47-sided r 21.5 wall
  meets the flange in a ring broken by the window blends' feet, 0.02-0.045 inside the circle; its
  windowed half was left to tilted "cylinders" r 16-31 across three strips and a B-spline. Where a
  cylinder region has three such chords - beside it, or past the facets the growth stops at - it
  grows through vertices as deep as its widest chord sags, and is fitted on the vertices of its
  chords that lie on it (with the chords' inner points the r 21.5 wall came out r 21.4989-21.4997
  from pass to pass, its corners a tolerance outside, and the region swung between 67 and 258
  facets). The mesh's chord-sag reference leaves a big sag out only if it also stands alone, half
  again the next one down: `test17`'s two big rounds, 0.050 and 0.063, sat just under three times
  the 90th percentile, and two small spheres more or less put the reference at 0.023 and the part
  at 989 faces. `test17`: 657 → 619 faces (deviation 0.096 → 0.067 mm); `test18`: 766 → 623 (0.028
  → 0.025); `test16` 864 → 863; `test15` 383 → 391; the CAD comparisons unchanged.
- **Filled patches on thin strips.** `test18`'s leg rounds come out of the seeds as strips of
  osculating cylinders (r 1.2 to 6, 20-80° each) whose facets run 1.6 mm along the leg and 0.07
  across, and no merge held: the filled patch's grid, its margin and its sampling step were sized
  by the median edge, 1.6 mm across a band 0.86 wide, and the samples off the band, where the
  interpolation runs wild, outweighed the data (3e-2 off at a tolerance of 7.7e-4). Past a 3:1
  spread the grid now follows the facets' spacing in each direction, and a strip that still doesn't
  hold on its mean plane, where it lies bent, is laid out conformally and tried again (2.8e-4).
  `test18`: 623 → 549 faces; `test17`: 619 → 589 (deviation 0.095 → 0.057 mm); `test16` 863 → 855;
  `test15` 391 → 383; `test9` 896 → 891; the rest and the CAD comparisons unchanged.
- **A coarse wall's notches, its holes.** At the flange, `test17`'s 47-gon wall ends in fans of
  long strips running down to the window blends' feet a few hundredths inside the circle, too
  shallow for a partial chord: the seeds covered them with three- and four-facet cylinders r 10 to
  27 and tori, or left them loose. On a cylinder that is a coarse polygon and goes half round at
  least, a neighbouring facet, loose or of a region under 16 facets, with no vertex outside, none
  deeper than the widest chord, two on it, its normal the wall's at its middle, and two sides or
  more on the wall already (a notch), goes over to it, without a refit. Taken one side at a time,
  the facets ran up the window blends' rims: 50 faces fewer, but some forty circles and ellipses of
  the wall's boundary became staircases of chords, and the boundary matters more than the count. A
  face going round can also have holes besides its two rings (the loops that go once round the axis
  are the rings, the others holes, their pcurves unwrapped into the rings' period). (Tried: a
  free-form patch's rim ceiling raised to the mesh's chord sag - 33 faces fewer, 9 circles fewer.)
  `test17`: 589 → 571 faces, the most circles and the fewest polyline edges between curved faces of
  any version; the rest unchanged.
  Nor a round's facet: a round leaving the wall tangentially has facets a few hundredths inside it
  with the wall's normal, and the top of a window arch lost two of its blend's to the wall, a knee
  of another one a third - notches in the arches' contours. A chord of an extruded wall runs along
  it (the one this is for, 37.6 mm), a round's facet doesn't (0.27-0.56 mm): a facet is taken only
  if it spans along the axis at least as far as the widest chord is wide. `test17` 489 → 480, the
  arches smooth again, the rest unchanged.
- **Latitude bands of a coarse round.** `test17`'s bore floor round (R8.5 r2, OpenSCAD's
  rotate_extrude in 47 sectors) came out as four spheres r 10.5-14, each exact through its two
  rings, and the strips between them: with 18 degrees between bands they don't meet tangentially,
  and two spheres were only merged when they were the same sphere. Two spheres with different
  centres are now tried as one torus, its axis the line through the centres, and kept only if the
  torus holds every vertex to a tenth of the tolerance (spheres r 29 and 40 osculating a free-form
  blend made tori within the tolerance and 0.09 off the part between the vertices). The torus then
  went round a window blend's island in a figure eight, and the file split the ring at the vertex
  into an outer wire and a hole touching it (InvalidImbricationOfWires): off the planes, axial
  surfaces and free-form patches a figure eight counts as rings touching, and the facets at the
  vertex go. `test17` 571 → 537, `test18` 549 → 517, `test16` 855 → 809 (deviation 0.053 → 0.024),
  `test15` 383 → 360 (0.012 → 0.008), `test6` 90 → 88; fewer polyline edges between curved faces
  everywhere, the CAD comparisons unchanged.
- **Slits and holes in free-form patches.** `test17`'s window bottoms came out as B-splines with
  thirteen holes in them (planar triangles of 0.01-0.1 mm2, and two B-splines of 0.3 mm2) and
  edges notched facet by facet, two patches meeting along a saw-tooth: loose facets their stage
  never handed them. Last of all, each free-form patch takes the loose facets with two sides on it,
  the loose islands beside it (bounded otherwise by big quadrics only: three rim facets between a
  window bottom and the wall, one side on each, stayed a grey triangle; bounded by any region, a
  bore-floor plane came out unmeshable) and the small patches (40 facets at most) it encloses,
  and the facets of the
  free-form patch next to it with two sides on it, if that one stays in one piece; only where its
  surface as fitted already holds them (every vertex within tol_fit, the normal within 20 degrees,
  the sag no deeper than the patch's own or the mesh's chords). Nothing is refitted; if the engine
  refuses the patch, the patch as it was is tried. `test17` 537 → 489 faces (holes in B-splines 13
  → 1, polyline edges between curved faces and planes 528 → 419), `test18` 517 → 489 (8 → 2, 460 →
  382), `test15` 360 → 354, `test16` 809 → 808, the rest unchanged. (Tried: the vertices as close
  as the patch's own worst one - `test15`'s thread patches 0.0089 off the mesh, 0.008 before, and no
  hole closed; the wall below the windows taking facets one side at a time where no other exact
  curved region touches them - 90 polyline edges more, 5 circles and 4 ellipses fewer; the strips
  under the windows are the circle's, finer than the 47-gon's sides, with one corner on the window
  bottom's rim 0.04 inside it.)
- **A window bottom in one piece.** `test17`'s window bottoms are one round wrapped 180 degrees
  round the slot's end, exported as corners and a middle whose seams carry vertex pairs 0.01 mm
  apart. A 180-degree bend has no mean plane, so the union of a bottom's halves is laid out
  conformally - and there each pair's offset comes out a few percent off its 3D one: the
  interpolation through the points, exact at them, swung 2.2 mm off the facets between them, and no
  grid held (3e-2 off at 42x23). When the conformal layout fails, it is smoothed (least squares of
  (u, v) on a quadratic in the points' coordinates) and tried again, if every triangle then turns
  the way it faces in 3D, with one more grid between 0.5 and 0.35 facets a span (everywhere, that
  grid made a strip under another window a face the readers couldn't mesh). A merge takes the loose
  islands its union encloses (40 facets at most): left out, the spline over them was held by
  nothing, and the late fill made each a patch of its own. And the late fill's patches - the wedge
  at a slot's corner, the row where the side round hands over - get the merges once more, after the
  walls have taken their fans back; two free-form patches only if the union needs no more than half
  again their poles (`test15`'s 6-facet corner made a 5x11 thread patch 18x48, and four faces more
  after it), a quadric's piece always (held to it, a torus r 23 osculating one of `test16`'s blends
  stayed, 0.032 off between its vertices). The small pieces the first merges left (40 facets at
  most, free-form or a quadric off its facets - an exact one stays) and what the walls cut get that
  pass too, each union with the loose facets in its notches: a merge is tried against the patch as
  it is at that moment, and a 23-facet strip along a slot edge that failed against the half-bottom
  it met first fits the whole bottom. A rejected merge gives way to its pieces, and a rejected
  piece to its own - except where the engine refused it for an edge left polygonal beyond the
  ceiling: `test17`'s bore-floor corners meet the floor's exact circle along facets 3.7 mm long,
  the circle a chord's sag (0.13) off the patch between their vertices, and each came back as the
  twenty tori, cylinders and loose facets it was made of. On the last pass such a patch (a merge or
  a section's) is cut in halves first and each filled again; if one goes in, the halves stay, and a
  half refused in turn stays mesh. (The halves were filled on facet indices that converting the
  first one renumbered: each is now held by its faces.)
- **A free-form patch held at its rim.** A B-spline patch is fitted to its facets' vertices and sag
  points, so along a long edge nothing holds it between two vertices: `test14`'s variable blend by
  the frame's front edge (CAD face 650) has a 17 mm edge whose polyline lies on the CAD to 0.004,
  while the spline bulged 0.07-0.15 off it, and a neighbour's arc r 13.35 ran outside the vertices'
  footprint, where the spline is only extrapolated - refused for an edge polygonal beyond the
  ceiling. Such a patch is now fitted once more with points along its edges every 0.5 mm as data of
  their own (weight 0.15, 0.3, then 0.6, the grid widened to cover them), kept only if the vertices
  stay within the tolerance, the interior no more than a quarter beyond the first fit and every edge
  under the ceiling, and the face is built again; if that fails too, the first fit's refusal stands.
  And an edge's tolerance on a free-form patch was measured to its vertical foot on the height
  field, not to the nearest point - on a steep flank 2-3 times the distance (`test17`'s slot ends,
  0.05-0.085 against a ceiling of 0.035): the pcurve now goes through the orthogonal foot (Newton,
  five steps). Between them, a 0.001 mm2 triangle beside a `test17` slot patch kept two edges, an
  arc and a chord 0.29 mm long a few microns apart: valid in memory, one edge once the STEP reader
  had merged them, and the part read back invalid. A face that would be left with two edges and is
  thinner than twice the tolerance keeps its polyline. `test17` 322 → 306 faces (deviation 0.067 →
  0.057), `test12` 3083 → 3046, `test14` 4240 → 4188 (0.274 → 0.260), the rest unchanged. Face 650
  itself is still mostly facets: the analytic edges of a 37-facet piece of it stand 0.3-0.5 off the
  patch, and fitted to them the spline folds (the face's (u, v) loop crosses itself) - the rollback
  leaves the mesh.
- **Spheres and tori as neighbours.** A new edge against a neighbour already built on a curved
  surface gets its pcurve there in a frame rebuilt from that face, and on the two closed surfaces
  the frame was wrong. A sphere is built with its X axis turned away from the reference it is
  given, so that its seam runs behind the face; rebuilt from the face's own X, the frame came out
  half a turn round the pole, every pcurve on a sphere neighbour sat on the far side of the
  sphere, π off in u, and the neighbour's wire no longer closed (`BadOrientationOfSubshape`,
  `Invalid`): all 219 refusals against a sphere in `test9`, `test12`, `test13` and `test14`. A
  torus goes round in v as well, and the new pcurve's v was left where the unwrapping put it, in
  (−π, π]: a tube face lying across v = π (`test8`'s r 0.5 corners, v 3.0 to 3.3) got its arc a
  whole turn away (36 refusals in `test8`, `test12` and `test14`). Those regions went in at a later
  try with every edge of their contour left a polyline, or not at all (`test14`: 94 and 3 refused,
  `test12`: 17, `test13`: 8, `test8`: 6, `test9`: 1). The sphere's frame is now read from the face
  as it is, like a plane's, and v is shifted into the face's own domain as u already was.
  And the patch refitted to its rim was judged only at the points it was fitted on, 0.5 mm apart
  and three on a short edge: on `test14`'s refused patches the boundary stood up to five times
  further off between them (0.005 at the samples, 0.027 in between), and `test17`'s slot patches
  passed to be refused right after; it is now judged every 0.1 mm, at least 17 points an edge.
  `test13`: still 57 of 57 CAD faces, no refusal against a sphere (16 before) and no region left
  with a polygonal contour (8), edge tolerance at most 0.0044 (0.0054); `test8`: no refusal
  against a torus (12), polygonal contours 8 → 2, the same 452 faces; `test9`: no refusal against a
  sphere (2), polyline edges 1,853 → 1,850, the same 886 faces.
- **Where u is cut.** A cylinder's or torus's u is unwrapped from the principal angle, so it is cut
  at −X, and X was the facets' mean direction: right while the face is an arc, but round a ring that
  is almost closed the mean is a short vector the facets' unevenness points anywhere. A face that
  covers more than half a turn now has its cut in the middle of the widest gap between its vertices.
  And a ring tessellated in uneven strips (`test12`'s 3.5 mm hole, 124 facets) has its widest gap
  between vertices under one strip - 11°, "covered" 0.969 - so it was taken for an open face and
  the cut fell under that strip: "wire not closed in parameter space, du 6.283". Whether the
  boundary goes round the axis was asked only above 0.97 of coverage; it is now asked from 0.9, and
  that hole goes in as a closed face. It is still refused there, by the plane whose contour runs
  along its rim, but no longer by its own parametrization.
- **A few facets on the surface to the float are a face.** A cylinder, cone, sphere or torus of
  fewer than 12 facets whose boundary lies mostly against loose facets is left as mesh: it is the
  noise of a fillet zone, fitted to a fraction of the tolerance (rms 1-2e-4 on the suite). A coarse
  export of a real face has few facets too: OpenCascade at 28° gave `test4`'s r 5.6 cylinder 7
  strips, 16 vertices on it to 1.3e-7, and it was left as mesh. A region whose vertices sit on its
  surface within a thousandth of the tolerance, three of them or more for every parameter of the
  primitive (three facets of anything fit "exactly"), now goes in (`test4` tessellated by
  OpenCascade at 0.05 mm / 28°: 117 → 111 faces; the 13-part suite unchanged).
- **Under a window the coarse wall is a fan.** Below each window bottom the 47-gon wall runs in a
  fan of thirty-odd slivers, from the polygon's sides at the flange (down to 0.047 inside) up to the
  round's foot on the true circle: twisted up to 17 degrees, one corner on the circle or two, never
  a chord, so the seeds covered them with cylinders r 5-25 and B-splines of three to seven facets.
  Once every stage has had its go, a coarse polygon cylinder takes the groups beside it that lie in
  its chords' band - no vertex outside, none deeper than its own deepest chord, one on it, the normal
  within 30 degrees, each facet as long along the axis as half the polygon's side (a round's rim
  facets at an arch are 0.3-0.6 mm; at half the side the pair of slivers where an arch's slope meets
  its knee came apart and stayed two grey triangles, at 0.35 it goes) - loose, or the part of a
  small region in the band if the rest holds together (taken only whole, a sliver beside a window
  stayed between the wall and a torus r 11 of four fan facets and four of the side round's), two
  facets or more. The
  wall then met the flange along the cutters' points, between the polygon and the circle: the exact
  circle stood off them and 105 segments stayed polylines. Last of all, on three segments or more,
  an exact curve may stand off the vertices as far as the edges it replaces stand off the faces
  (tried first, or on two segments, it took the place of the circles beside the arches and failed
  their pcurves).
  The two together: `test17` 480 → 329 faces (each window bottom one B-spline, the skirts under
  them and the slivers beside the arches gone, the bore floor's corners in a few pieces; polyline
  edges between curved faces 765 → 495, between curved faces and planes 405 → 266, between planes
  129 → 78), `test18` 489 → 316 (382 → 254 between curved faces and planes, 196 → 158 between
  planes; deviation 0.025 → 0.024), `test16` 808 → 801, `test9` 890 → 884, the rest and the CAD
  comparisons unchanged. (Tried: a patch taking the teeth of a quadric next to it, as it takes a
  free-form one's - the steps where a window's side round meets its bottom - 25 polyline edges
  fewer between curved faces, but the bore floor's last merges then left a plane bounded by two
  lines on top of each other, a face no reader meshes; a big patch taking a cylinder's teeth only,
  `test16` 801 → 833.)
- **Five windows, one bottom each.** `test17`'s five windows are copies of one mesh, and their
  bottoms still came out three ways. On two of them a strip one facet wide across the bottom fitted
  a sphere to 2.3e-5 (r 2.08 on one, 2.45 on the other, a torus on the rest): every vertex on one
  of its two long sides, and two samples across say nothing of the curvature across. Such a sphere
  or torus, no vertex inside its region, is no longer exact for the merges (a cylinder's band of
  full-length facets is held by its two circles and stays). On the window at -20 degrees, three
  rows of the bottom's middle had come out a tube patch (`_tube_patches`, three strips a quarter
  turn round the slot, 0.45 of the tolerance off its facets), and a tube never merged: the bottom
  stayed in three faces. A tube merges now, unless it holds its facets as a CAD's own tube does,
  within a fifth of the tolerance (`test9`'s pipes, 0.01-0.1 of it: merged with the pieces around
  them they were refused at their coarse edge, and the pieces given back last left 22 planes of mesh
  beside them). And a patch takes the teeth of an exact quadric next to it, as it takes a free-form
  one's (an arch leg's cylinder, on its vertices to 1e-7, grown a facet past its seam into the
  knee), now that the final Phase A can't leave a face two edges: two segments in line between the
  same two faces are one line only if each face keeps three edges, wherever every face around it is
  curved (beside a plane the lens still collapses, as the raw mesh's needles should: held there too,
  `test17` came out with a hundred faces more). Not a noisy quadric's: `test16` gave 50 facets and
  came out with 32 faces more. `test17` 329 → 323 (each of the five bottoms one B-spline; polyline
  edges between curved faces 495 → 461), `test18` 316 → 312 (between curved faces 1044 → 1007, curved
  and planes 254 → 240; deviation 0.024 → 0.023), `test15` 354 → 353, the rest and the CAD
  comparisons unchanged. (Tried: halving a refused merge only when it was made of four pieces or
  more - `test17` 329 → 338, a bore-floor corner of three came back as its pieces.) Still stepped:
  each window's knee where the arch leg's side wall meets it - Phase A merged the first triangles of
  the knee's round into that plane, within its tolerance, and the plane is a single facet to the
  segmenter; the right knee at 140 degrees is meshed coarser than its twins and keeps the exact torus
  and cylinder its facets lie on.

The phases run in the order they're written on the command line. With no flags the sequence is
`A C A(0.005) C A(0.005)`, the one that came out best on the test parts, every sequence run on all
of them:

- Phase B before C adds nothing (Phase C finds the holes itself) and costs time;
- re-merging the planar faces at 0.005 mm leaves the part exactly as close to the mesh as the tight
  default, with far fewer faces;
- a second Phase C after that re-merge finds what the first one couldn't, because the merged planes
  give the curved regions their true boundaries (`test2` 96 → 82 faces, `test6` 98 → 89); a third
  finds nothing;
- a wide first Phase A flattens fillets into planes before Phase C can see them.

```bash
python refit.py part.stl                     # the default: A C A(0.005) C A(0.005)
python refit.py part.stl -a -c -a 0.005      # one Phase C only, faster
python refit.py part.stl -a -b -c -a 0.005   # holes first, the previous recipe
```

### Speed and `-j`

`-j N` starts N worker processes at once (default: the machine's cores minus two, at most 22) and
every phase that follows uses them: the seeds of the curved phases (both rounds), the type review,
the exact consensus and the refits of the segmentation, the fits of the promoted merges and the
checks of the arcs. Phase A and the in-solid replacements stay on one core: they change a shared
OpenCascade structure and are checked one at a time. The result doesn't depend on `-j`: the workers
compute against a snapshot and the main process replays their answers in the sequential order,
taking one only where nothing it depended on has changed since, otherwise computing it again
itself (the same regions with `-j 1` and `-j 22` on `test15`, `test16` and `test17`). On a
24-core machine, against the previous version: `test16` 2,055 → 409 s, `test15` 1,414 → 900 s,
`test17` 546 → 228 s, each with the very same file, and `test14` 1,971 → 1,486 s; the 13-part
suite comes out identical, in about half the time. With `-v` the log lists where each segmentation spent its time
(`Segmentation steps: ...`).
The merges of weak regions into one B-spline try a group's neighbours in turn until one union
holds - a fit of 0.1-0.3 s each, 766 of them on the closed section of `test18`'s window, 618 that
didn't hold. Nothing a union reads changes until a merge, which ends the step, so all of a step's
unions are built first, fitted on the workers together and read in the same order: the first that
holds is the merge, the fits after it are dropped unread. The same merges as on one core (the
section's regions compared facet by facet), 175 → 65 s; each facet's polygon is walked once per
merge pass instead of once per union.
The mesh's topology is rebuilt at every phase and for every engine, and each build looks for the
vertices that sit on a chord of the mesh (T-vertices, see the curved phases): that search tried
every pair of a vertex's neighbours once per facet, with a numpy cross product each, so every vertex
six times over. It is now done once per vertex and kept, with the facet's own preference applied
after, with the very same numpy operations (redone in plain Python, the sums in another order moved
`test5`'s volume by a thousandth of a mm³: a vertex on a straight CAD edge sits at the float
threshold); the cache goes with every change of the mesh. The same vertices to the last bit
(`test16`, `test5`, `test4`, `test11`, `test0` and a section of `test14` compared facet by facet),
the topology built 4-6 times faster (`test16` 20.0 → 4.8 s, a 15,000-triangle section of `test14`
13.4 → 2.4 s).

## Output file

After the phases, the edges that are still polylines of the mesh get the curve they really sit
on: the exact intersection of the two faces, a circle, or a spline through its vertices (checked on
both faces). A chain that doesn't take one curve as a whole is cut at its corners, or halved when
it has none (a short one loses its last or first segment), and tried again in pieces, and the step
repeats while it finds something. A chain of only two segments gets the circle through its three
vertices when that circle lies on both faces (where a rolling-ball fillet hands over to the
cylinder its ball sweeps next, the two meet along the ball's own section). Between two tangent
faces the spline leaves out the vertices where the two faces' normals part by more than 1.5°:
they're off the line of contact, a row of the fillet's facets left with the wall, and the spline
through them waved (`test13`: its fillets' rails within 6 µm of the CAD's edges instead of 67).
Before that, a needle triangle left on a rim - two edges between the same two faces running out to
a tip and back - is collapsed, so the rim runs through as one chain (`test13`: edges 163 → 139, no
more kinks at the fillets' corners).

The STEP file is written without pcurves (the edges' 2D curves on the faces), as CAD systems do:
the reader rebuilds them, and the file is 3-7 times smaller. The output is read back right away
and checked (validity, free edges, every face triangulated against its own area); the faces that
come back worse get their pcurves written after all. `--pcurves` always writes them all.

STEP keeps no tolerances, so a face that was valid in memory only thanks to a vertex's tolerance
(two edges reaching their vertex at a narrow angle and crossing just past it, or a lens-shaped
sliver whose ring touches the face's rim) can read back invalid with or without pcurves. When that
happens the part as read back is mended — the vertex moved to the edges' crossing, the sliver
covered by the face around it — and written again, kept only if it reads back better and with the
same volume, also when the reader has cut a face in two. On a planar face the crossing is looked
for in the face's own (u, v), at any angle and between two straight edges too: two skew lines never
meet in space, but `test12`'s 5.3 and 13.8 mm edges at a right-angled corner touched in a T 2.1e-4
from their vertex once read back. A hole may be filled by two or three slivers together (an
ellipse, its chord and a 0.03 mm edge): they go as one. `test12`, whose file read back with one or
three invalid faces, now reads back valid.
The reader's own pcurves can also break a face the part had whole: `test12`'s 55×26 patch over a
blend (115 mm²) has an edge 6.7e-3 off its surface, its projected pcurve came back ending 0.003
off its neighbours' in (u, v), the reader closed both gaps with degenerate edges and no mesher
could triangulate the face — a hole in any viewer, 3.3 mm off the mesh. Since the part read back
invalid too and was mended on what the reader made of it, the broken face had become the reference
and the file passed as faithful. Before writing, every B-spline face is now written alone without
pcurves and read back (0.04 s each); those that come back unable to mesh carry their pcurves in
every file written (`test12`: 1 face, valid, 0.055 mm from the mesh, 12.4 → 16.5 MB). And a
region whose new face the mesher can't triangulate at all is refused in the engine, instead of
passing on the global volume. A gap in (u, v) that a vertex's tolerance covers in memory can also
read back as a pole: on the section of `test18`'s window a 41×36 patch met a circle 7e-3 off it
whose pcurve ended 0.001 in (u, v) from the next line's - 1.2 µm on the surface, under the 5e-3 the
file declares as its precision - and the reader closed it with a degenerate edge, a hole of 16 mm²
with pcurves or without. When the pcurves don't help, the faces to blame have their pcurves' ends
moved onto one point at each vertex, on a copy (a face left invalid gets its pcurves back), and the
copy is written if it reads back faithful. A
degenerate edge on a surface that has no pole there (a plane, a cylinder, a torus with R > r) is
dropped before writing: the reader made a wire of its own out of one, and a 0.33 mm² face of a
`test14` spindle torus came back as the whole torus, 36,342 mm².

## Repository layout

- [`refit.py`](../refit.py) — the main tool, run directly from the repo root.
- [`Refit Viewer.pyw`](../Refit%20Viewer.pyw) — double-click launcher for the viewer (no console).
- [`viewer/`](../viewer) — the local viewer (`python viewer/server.py`, see above).
- [`tools/`](../tools) — QA, analysis and debug scripts, all run from the repo root
  (e.g. `python tools/qa.py ...`); each imports `refit.py` by locating it relative to its own path,
  so the repo can be moved around freely as long as this layout is kept:
  - [`tools/qa.py`](../tools/qa.py) — checks an output against its input mesh: exact two-sided deviation
    (output → mesh and mesh → output, so holes in the output show up too), validity, free edges,
    face and edge types, tolerances. `python tools/qa.py part.stl part.step`
  - [`tools/verify.py`](../tools/verify.py) — measures the deviation between an output and the starting
    mesh (one-sided, sampled).
  - [`tools/archi.py`](../tools/archi.py) — counts how many of the mesh's polylines are actually
    circular arcs.
  - [`tools/tassellate.py`](../tools/tassellate.py) — estimates how much curvature is left tessellated
    in the output.
  - [`tools/cadcmp.py`](../tools/cadcmp.py) — compares an output face by face with the CAD file the mesh
    came from: faces reproduced one to one, split or of the wrong type, radii and angles, mesh left
    over, polylines where the CAD has a curve. `python tools/cadcmp.py test/test10_original.step test/test10.step -v`
  - other scripts (`align.py`, `cadfit.py`, `cadinfo.py`, `extra.py`, `icp.py`, `leak.py`,
    `mancanti.py`, `rbrep.py`, `segnali.py`, `strips.py`, `vstrips.py`) — analysis and debug tools
    used during development.
- [`test/`](../test) — `test0.stl` … `test18.stl`, the test meshes (`test0` the simplest). `test0`, `test4`, `test10`,
  `test12`, `test13` and `test14` come with the CAD file they were exported from
  (`testN_original.step`): the reference for `tools/cadcmp.py`. `test12` is the hardest one: 1,328
  CAD faces (169 free-form), six internal cavities. `test13` has rolling-ball fillets between two
  crossing holes. `test14` is the biggest (189k triangles, 965 CAD faces, most of its 229 B-spline
  faces extruded spline walls); its `_original.step` is the part alone, a surface model left in the
  assembly's frame (`cadcmp.py` aligns it).

## Branches

- `main` — documentation and comments in English.
- `ita` — documentation and comments in Italian (original development language).

## License

MIT — see [LICENSE](../LICENSE).
