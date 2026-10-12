# STL2STEP

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
python refit.py test0.stl
```

Then check the result against the input mesh:

```bash
python tools/qa.py test0.stl test0.step
```

and, for the test files that ship with their original CAD (`test0`, `test4`, `test10`, `test12`,
`test13`), compare face-by-face against the ground truth:

```bash
python tools/cadcmp.py test0_original.step test0.step -v
```

See the docstring at the top of [`refit.py`](refit.py) for the full guide (every flag, how to read
the report, the final arc-snapping step).

### Viewer

```bash
pip install pywebview
python viewer/server.py              # or double-click viewer/Refit Viewer.pyw
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
  coarsely that the mesh doesn't say where the curve runs stay planar. Where a free-form band runs
  into a round it's tangent to, the facets along the contact line sit on both: one with more sides
  on its neighbour than on its own region, and on the neighbour's surface, goes over to it (`test4`:
  no more notches at the ends of the chamfer bands, edges 137 → 126). Every facet is a chord of
  the surface under it: a region whose facets sag below its surface further than the mesh's own
  chords ever do (measured on the regions that fit exactly) has the wrong surface - a sphere
  standing over a flat disc - and is rejected, or only the facets that don't belong are left out. A constant-radius fillet between two
  curved walls (a hole meeting another hole, whose width changes all along) is rebuilt from the
  two walls and the ball rolling between them, as one B-spline face whose edges follow the ball's
  contact lines; the straight and circular pieces of the same chain become the exact cylinders and
  tori the ball sweeps there, and each fillet ends where its ball leaves the walls' facets. A facet
  left over at a fillet's corner, wedged among faces already found, joins the fillet when its
  vertices are on it and its middle sags no deeper than a chord of the ball does; a fillet that then
  fails to convert is converted without it (`test13` reproduces all 57 faces of its CAD one to one, and has no others; `test8` 468 → 463
  faces).

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
same volume (`test12` now reads back valid), also when the reader has cut a face in two. A
degenerate edge on a surface that has no pole there (a plane, a cylinder, a torus with R > r) is
dropped before writing: the reader made a wire of its own out of one, and a 0.33 mm² face of a
`test14` spindle torus came back as the whole torus, 36,342 mm².

## Repository layout

- [`refit.py`](refit.py) — the main tool, run directly from the repo root.
- [`viewer/`](viewer) — the local viewer (`python viewer/server.py`, see above).
- [`tools/`](tools) — QA, analysis and debug scripts, all run from the repo root
  (e.g. `python tools/qa.py ...`); each imports `refit.py` by locating it relative to its own path,
  so the repo can be moved around freely as long as this layout is kept:
  - [`tools/qa.py`](tools/qa.py) — checks an output against its input mesh: exact two-sided deviation
    (output → mesh and mesh → output, so holes in the output show up too), validity, free edges,
    face and edge types, tolerances. `python tools/qa.py part.stl part.step`
  - [`tools/verify.py`](tools/verify.py) — measures the deviation between an output and the starting
    mesh (one-sided, sampled).
  - [`tools/archi.py`](tools/archi.py) — counts how many of the mesh's polylines are actually
    circular arcs.
  - [`tools/tassellate.py`](tools/tassellate.py) — estimates how much curvature is left tessellated
    in the output.
  - [`tools/cadcmp.py`](tools/cadcmp.py) — compares an output face by face with the CAD file the mesh
    came from: faces reproduced one to one, split or of the wrong type, radii and angles, mesh left
    over, polylines where the CAD has a curve. `python tools/cadcmp.py test10_original.step test10.step -v`
  - other scripts (`align.py`, `cadfit.py`, `cadinfo.py`, `extra.py`, `icp.py`, `leak.py`,
    `mancanti.py`, `rbrep.py`, `segnali.py`, `strips.py`, `vstrips.py`) — analysis and debug tools
    used during development.
- `test0.stl` … `test14.stl` — test meshes (`test0` the simplest). `test0`, `test4`, `test10`,
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

MIT — see [LICENSE](LICENSE).
