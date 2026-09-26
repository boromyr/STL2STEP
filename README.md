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
python refit.py part.stl -a -b -c -a 0.005
```

This runs the default phase sequence (`A B C A`, see below) with tolerance `0.005` and writes
`part.step` next to the input. A few of the bundled test meshes are a good first try:

```bash
python refit.py test0.stl -a -b -c -a 0.005
```

Then check the result against the input mesh:

```bash
python tools/qa.py test0.stl test0.step
```

and, for the test files that ship with their original CAD (`test0`, `test4`, `test10`, `test12`),
compare face-by-face against the ground truth:

```bash
python tools/cadcmp.py test0_original.step test0.step -v
```

See the docstring at the top of [`refit.py`](refit.py) for the full guide (every flag, how to read
the report, the final arc-snapping step).

## Mesh check (STL input)

Before anything else the mesh is checked: triangles with a repeated vertex or zero area, duplicate
triangles, holes, non-manifold edges and vertices, inconsistent winding, self-intersections, separate
pieces (bodies, or internal cavities: a piece of negative volume inside the body becomes a void of
the solid, as in the CAD's `BREP_WITH_VOIDS`). The defects that would break the conversion are
repaired with [MeshLab](https://www.meshlab.net/) through `pymeshlab`, one filter per defect and never
by remeshing, so the vertices of a CAD export stay exact:

- zero-area triangles with a vertex on their side: flipped with the neighbor (no crack);
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
  profile is a circular arc, otherwise a B-spline profile. Stretches tessellated so coarsely that
  the mesh doesn't say where the curve runs stay planar.

The phases run in the order they're written on the command line. With no flags the sequence is
`A B C A` (the last one re-merges the planar faces split by the replacements).

```bash
python refit.py part.stl -a -b -c -a 0.005
```

## Repository layout

- [`refit.py`](refit.py) — the main tool, run directly from the repo root.
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
- `test0.stl` … `test13.stl` — test meshes (`test0` the simplest). `test0`, `test4`, `test10` and
  `test12` come with the CAD file they were exported from (`testN_original.step`): the reference
  for `tools/cadcmp.py`. `test12` is the hardest one: 1,328 CAD faces (169 free-form), six internal
  cavities.

## Branches

- `main` — documentation and comments in English.
- `ita` — documentation and comments in Italian (original development language).

## License

MIT — see [LICENSE](LICENSE).
