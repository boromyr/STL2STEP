# Refit — turn a mesh into a real CAD model

**Got an STL and need a STEP you can actually edit?** Refit rebuilds your mesh as clean, analytic
CAD geometry: flat faces become planes, holes become true cylinders, fillets and chamfers become
real fillets and chamfers. Free-form surfaces become smooth B-splines instead of thousands of
tiny triangles.

| Before: the mesh (2,064 triangles) | After: the CAD model (72 faces) |
| :---: | :---: |
| ![Input mesh of test0.stl](docs/before.webp) | ![STEP output of test0.stl](docs/after.webp) |

## Why you'll like it

- **Safe by design.** Every change is small, local, and checked right away: the solid must stay
  closed, valid and faithful to your mesh. If a change doesn't pass, it's thrown away. You never
  get a broken file, at worst a part stays a little less "clean".
- **Faithful to your part.** The output stays within a tight tolerance of the original mesh. Run
  the built-in checks to see exactly how close.
- **Smart about curves.** Holes, fillets, chamfers, countersinks, spheres, threads, embossed text:
  Refit recognises them and gives each one its exact surface.
- **Fixes broken meshes for you.** Holes, duplicate triangles, flipped normals and similar defects
  are detected and repaired before conversion, without moving your vertices.
- **See it before you trust it.** A built-in viewer shows the mesh and the result side by side,
  colours faces by type, and runs the conversion for you.
- **Fast.** It uses all your CPU cores, and the result is identical whether it runs on one core or
  twenty.

## Quick start

**1. Install** (Python 3.9 or newer):

```bash
pip install cadquery-ocp numpy scipy pymeshlab
```

`pymeshlab` is optional: it enables automatic mesh repair. For the viewer also run
`pip install pywebview`.

**2. Convert** your mesh:

```bash
python refit.py part.stl
```

You get `part.step` next to your STL. Want to try it first? Use one of the bundled examples:

```bash
python refit.py test/test0.stl
```

**3. Check the result** against the original mesh:

```bash
python tools/qa.py test/test0.stl test/test0.step
```

## The viewer

Prefer clicking to typing? Double-click **`Refit Viewer.pyw`** in this folder (or run
`python viewer/server.py`). Pick a mesh, press *Convert to STEP*, and watch the result appear next
to your input. Hover any face to see its type, radius and area.

## Good to know

- Refit never rebuilds your model from scratch. Whatever it can't convert with confidence simply
  stays as it was, so the output always matches the input.
- The conversion runs in phases (flat faces, holes, curved surfaces). The default sequence is a
  good choice for most parts; every option is explained at the top of [`refit.py`](refit.py).
- If `import OCP` fails with "DLL load failed", pin `vtk==9.6.2` (`pip install vtk==9.6.2`).

## What's in the box

| Folder / file | What it is |
| --- | --- |
| [`refit.py`](refit.py) | The converter. |
| [`Refit Viewer.pyw`](Refit%20Viewer.pyw), [`viewer/`](viewer) | The desktop viewer. |
| [`test/`](test) | Example meshes, some with the original CAD file for comparison. |
| [`tools/`](tools) | Quality checks and analysis scripts (`qa.py`, `cadcmp.py`, ...). |
| [`docs/technical.md`](docs/technical.md) | How it works, in depth: phases, tolerances, mesh repair, output details. |

## License

MIT — see [LICENSE](LICENSE).
