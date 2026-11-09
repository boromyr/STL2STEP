# Refit

**Reconstruction of analytic CAD geometry from triangle meshes.**

Refit converts an STL mesh (or a STEP file that originated from a mesh) into a STEP model made of
true CAD surfaces. Flat regions become planes, holes become cylinders, fillets and chamfers become
fillets and chamfers, and free-form areas become B-spline surfaces, in place of thousands of
triangles.

| Input mesh (2,064 triangles) | Output B-Rep (72 faces) |
| :---: | :---: |
| ![Input mesh of test0.stl](docs/before.webp) | ![STEP output of test0.stl](docs/after.webp) |

## Features

- **Reliable output.** Every modification is local and verified immediately: the solid must remain
  closed and valid, with consistent orientation and with area and volume matching the mesh. A
  modification that fails the check is discarded, so the result is always a consistent model.
- **Geometric fidelity.** The output stays within a tight tolerance of the source mesh, and
  bundled quality tools report the measured deviation.
- **Feature recognition.** Holes, fillets, chamfers, countersinks, spheres, bosses, threads and
  embossed or engraved text are identified and given their exact surface.
- **Automatic mesh repair.** Holes, duplicate and degenerate triangles, non-manifold edges and
  inconsistent orientation are detected and repaired before conversion, without moving vertices.
- **Integrated viewer.** A desktop application displays the input mesh, the converted model and
  the original CAD side by side, with faces colored by surface type, and lets you fix by hand the
  faces that came out in pieces.
- **Parallel execution.** All CPU cores are used, and the result is identical regardless of the
  number of processes.

## Installation

Requires Python 3.9 or later. Missing dependencies are checked at startup (`refit.py`, the viewer or
`python deps.py`) and installed automatically with pip; set `REFIT_NO_AUTOINSTALL=1` to disable that.
To install them by hand:

```bash
pip install cadquery-ocp numpy scipy pymeshlab
```

`pymeshlab` is optional and enables automatic mesh repair. The viewer additionally requires
`pywebview`:

```bash
pip install pywebview
```

If `import OCP` fails with "DLL load failed", install the matching version of VTK
(`pip install vtk==9.6.2`).

## Usage

Convert a mesh:

```bash
python refit.py part.stl
```

The result is written to `part.step`, next to the input file. The `test/` folder contains example
meshes:

```bash
python refit.py test/test0.stl
```

Verify the result against the source mesh:

```bash
python tools/qa.py test/test0.stl test/test0.step
```

Phases, tolerances and all other options are documented in the header of [`refit.py`](refit.py).

## Viewer

Launch the viewer with `Refit Viewer.pyw` (double-click) or:

```bash
python viewer/server.py
```

The viewer allows selecting a mesh, running the conversion and inspecting the result. The type,
radius and area of each face are shown on hover.

### Manual face correction

Faces that came out in pieces can be fixed by hand: press **Fix faces** (E) on the STEP, click one
of the pieces (the selection grows to the pieces like it; Shift+click or Shift+drag paints single
faces in or out, and only the selection is coloured), choose what they are — Auto, Plane,
Cylinder, Cone, Torus, Sphere, Tube or Free-form — and press **Build**. The choices are saved next to the mesh as `<name>.hints.json`, and the conversion
(`--hints`) turns each of them into one face of that type (two, for a Free-form or Tube round that
runs all the way round a hole).

## Repository contents

| Path | Description |
| --- | --- |
| [`refit.py`](refit.py) | Converter. |
| [`deps.py`](deps.py) | Startup dependency check and automatic install. |
| [`Refit Viewer.pyw`](Refit%20Viewer.pyw), [`viewer/`](viewer) | Desktop viewer. |
| [`test/`](test) | Example meshes, some with the original CAD file for comparison. |
| [`tools/`](tools) | Quality checks and analysis scripts (`qa.py`, `cadcmp.py`, ...). |
| [`docs/technical.md`](docs/technical.md) | Technical reference: phases, tolerances, mesh repair, output details. |

## License

MIT, see [LICENSE](LICENSE).
