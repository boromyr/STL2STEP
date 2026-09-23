# STL2STEP

"Soft" reconversion of a mesh (STL, or a STEP born from a mesh) into an analytic B-Rep using OpenCascade.

The principle is a single one: **every change to the geometry is local**, is verified right away
(solid still closed, faces valid, consistent orientation, area and volume consistent with the mesh)
and is discarded if it doesn't pass the check. Whatever can't be converted stays tessellated as it
was: the output file is always consistent with the input one, at worst it's just less "clean". There
is no global Sewing and nothing is rebuilt globally.

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
pip install cadquery-ocp numpy scipy
```

Works with OCP 7.x and 8.x. `cadquery-ocp` 8.0.1 is built against `vtk==9.6.2`: with a newer vtk
(e.g. after a blanket `pip` upgrade) `import OCP` fails with "DLL load failed". `scipy` is only
needed by `qa.py`.

## Files

- `refit.py` — the main tool.
- `qa.py` — checks an output against its input mesh: exact two-sided deviation (output → mesh and
  mesh → output, so holes in the output show up too), validity, free edges, face and edge types,
  tolerances. `python qa.py part.stl part.step`
- `verify.py` — measures the deviation between an output and the starting mesh (one-sided, sampled).
- `archi.py` — counts how many of the mesh's polylines are actually circular arcs.
- `tassellate.py` — estimates how much curvature is left tessellated in the output.
- `test0.stl` … `test9.stl` — test meshes, ordered by increasing facet count (`test0` the
  simplest). `test4_original.step` is the CAD file `test4.stl` was exported from: the reference
  for comparing an output face by face (`extra.py test4_original.step test4.step`).
- other scripts (`align.py`, `cadfit.py`, `cadinfo.py`, `extra.py`, `icp.py`, `leak.py`,
  `mancanti.py`, `rbrep.py`, `segnali.py`, `strips.py`, `vstrips.py`) — analysis and debug tools
  used during development.

## Branches

- `main` — documentation and comments in English.
- `ita` — documentation and comments in Italian (original development language).

## License

MIT — see [LICENSE](LICENSE).
