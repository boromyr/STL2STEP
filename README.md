# STL2STEP

"Soft" reconversion of a mesh (STL, or a STEP born from a mesh) into an analytic B-Rep using OpenCascade.

The principle is a single one: **every change to the geometry is local**, is verified right away
(solid still closed, faces valid, consistent orientation, area and volume consistent with the mesh)
and is discarded if it doesn't pass the check. Whatever can't be converted stays tessellated as it
was: the output file is always consistent with the input one, at worst it's just less "clean". There
is no global Sewing and nothing is rebuilt globally.

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

The phases run in the order they're written on the command line. With no flags the sequence is
`A B C A` (the last one re-merges the planar faces split by the replacements).

```
python refit.py part.stl -a -b -c -a 0.005
```

See the docstring at the top of `refit.py` for the full guide (tolerances, how to read the report,
the final arc-snapping step).

## Requirements

```
pip install cadquery-ocp numpy scipy pymeshlab
```

Works with OCP 7.x and 8.x. `cadquery-ocp` 8.0.1 is built against `vtk==9.6.2`: with a newer vtk
(e.g. after a blanket `pip` upgrade) `import OCP` fails with "DLL load failed". `scipy` is needed by the mesh
check and by `qa.py`; `pymeshlab` (optional) repairs the mesh defects.

## Files

- `refit.py` — the main tool.
- `qa.py` — checks an output against its input mesh: exact two-sided deviation (output → mesh and
  mesh → output, so holes in the output show up too), validity, free edges, face and edge types,
  tolerances. `python qa.py part.stl part.step`
- `verify.py` — measures the deviation between an output and the starting mesh (one-sided, sampled).
- `archi.py` — counts how many of the mesh's polylines are actually circular arcs.
- `tassellate.py` — estimates how much curvature is left tessellated in the output.
- `cadcmp.py` — compares an output face by face with the CAD file the mesh came from: faces
  reproduced one to one, split or of the wrong type, radii and angles, mesh left over, polylines
  where the CAD has a curve. `python cadcmp.py test10_original.step test10.step -v`
- `test0.stl` … `test13.stl` — test meshes (`test0` the simplest). `test0`, `test4`, `test10` and
  `test12` come with the CAD file they were exported from (`testN_original.step`): the reference
  for `cadcmp.py`. `test12` is the hardest one: 1,328 CAD faces (169 free-form), six internal
  cavities.
- other scripts (`align.py`, `cadfit.py`, `cadinfo.py`, `extra.py`, `icp.py`, `leak.py`,
  `mancanti.py`, `rbrep.py`, `segnali.py`, `strips.py`, `vstrips.py`) — analysis and debug tools
  used during development.

## Branches

- `main` — documentation and comments in English.
- `ita` — documentation and comments in Italian (original development language).

## License

MIT — see [LICENSE](LICENSE).
