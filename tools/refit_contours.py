"""
refit_contours.py - refit.py helped by contours.py: convert, find the faces of the mesh the conversion cut into pieces,
convert again with those faces as hints.

    python tools/refit_contours.py MESH.stl [-o OUT.step] [--from FIRST.step] [refit.py options...]

1. refit.py MESH -o OUT ...                     the plain conversion
2. contours.py MESH --step OUT --hints H        the faces contours.py reads from the mesh (whole) that the conversion
                                                cut (see Contours.hints): written beside OUT, <OUT>.contours.json
3. refit.py MESH -o OUT --hints H ...           again with them (a --hints file given is kept: its groups first)

Nothing to hint -> the first conversion stays. The first one is kept as <OUT>_pass1.step.
--from FIRST.step: a conversion of MESH already made (the viewer's "Refine with contours") is the first: steps 2-3
only, half the time. --from OUT where OUT is this script's own second pass: its <OUT>_pass1.step is the first.
Measured 20/11/2026 against the true CAD (whole faces / CAD edges caught / false contours / faces):
  test12  943 -> 969 / 89.6 -> 89.8% / 5128 -> 3139 mm / 2811 -> 2276, largest tolerance the same (0.12 mm)
  test12_foro (the hexagonal hole) 54 -> 59 / 64.6 -> 64.6% / 1498 -> 995 mm / 503 -> 278
  test0, test4, test13: no hint (the same conversion); test10: 3 hints, the same faces
No true CAD (edges of the output on the mesh's contours / faces / largest tolerance):
  test18  65.9 -> 66.7% / 297 -> 289 (cylinders 104 -> 98, tori 21 -> 19) / 0.025 -> 0.026 mm
  test15  81.9 -> 83.1% / 177 -> 168 (cylinders 34, tori 10 -> 9) / 0.061 mm the same
  test9   87.5 -> 99.5% / 883 -> 526 (planes 531 -> 186, cylinders 122 -> 159, cones 37 -> 53) / 0.070 -> 0.102 mm:
          its coarse cylinders and chamfer cones (r 250-400 mm, three strips of two triangles, vertices within 1e-7)
          one face each - its lettered walls one face each, not a fan of flats (contours.grow_axial) - its coarse
          blends one B-spline each (hinted in few pieces when those are B-splines or planes); their rims' single
          chords made arcs at the end (snap_arcs), the largest left a chord to a B-spline. (Its first conversion
          1041 -> 883 faces once fit_axial tried the cylinder across a coarse wall's skewed facets too)
          on the viewer's options (-a -c -a 0.005): 1041 -> 563 faces / 0.15 mm, no single-triangle flat left
Time: the conversion twice (test12 690 + 810 s, test18 310 + 470 s).
"""
import sys, os, json, shutil, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
REFIT = os.path.join(os.path.dirname(HERE), "refit.py")


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0].startswith("-"):
        print(__doc__); return 2
    mesh = argv[0]
    rest = argv[1:]
    out = None
    for k in ("-o", "--output"):
        if k in rest:
            i = rest.index(k); out = rest[i + 1]; del rest[i:i + 2]
    first = None
    if "--from" in rest:
        i = rest.index("--from"); first = rest[i + 1]; del rest[i:i + 2]
    user_hints = None
    if "--hints" in rest:
        i = rest.index("--hints"); user_hints = rest[i + 1]; del rest[i:i + 2]
    out = out or os.path.splitext(mesh)[0] + ".step"
    base = os.path.splitext(out)[0]
    run = lambda extra: subprocess.call([sys.executable, REFIT, mesh, "-o", out] + rest + extra)
    p1, hf = base + "_pass1.step", base + ".contours.json"

    def stamp():
        """out's size and mtime: written into the hints file once the second pass made it"""
        st = os.stat(out)
        return [st.st_size, st.st_mtime_ns]

    # ⚠️ REFINING A REFINED RESULT STARTS FROM ITS FIRST PASS: the viewer's "Refine" again (--from OUT, -o OUT) would
    # copy the second pass over <OUT>_pass1.step - the first conversion lost, 5-12 minutes to make again
    if first and os.path.abspath(first) == os.path.abspath(out) and os.path.isfile(p1) and os.path.isfile(out):
        try:
            with open(hf, encoding="utf-8") as fh:
                if json.load(fh).get("pass2") == stamp():
                    first = p1
        except (OSError, ValueError):
            pass

    if first:
        if not os.path.isfile(first):
            print(f"[x] --from {first}: no such file"); return 1
        print(f"== refit_contours: 1/3 the conversion already made: {first}", flush=True)
        if os.path.abspath(first) != os.path.abspath(out):
            shutil.copyfile(first, out)
    else:
        print(f"== refit_contours: 1/3 conversion -> {out}", flush=True)
        if run(["--hints", user_hints] if user_hints else []) != 0:
            return 1
    print("== refit_contours: 2/3 the faces it cut (contours.py)", flush=True)
    sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
    try:
        import contours as C
        c = C.Contours.from_file(mesh)
        hs = c.hints(C.R.read_input(out)[0])
    except Exception as e:      # (the first conversion is there: a failure here costs only the second)
        print(f"   contours.py failed ({type(e).__name__}: {e}): the conversion stays")
        return 0
    if not hs:
        print("   none: the conversion stays")
        return 0
    groups = (json.load(open(user_hints, encoding="utf-8")).get("groups", []) if user_hints else []) + hs
    doc = {"refit_hints": 1, "groups": groups}
    with open(hf, "w", encoding="utf-8") as fh:
        json.dump(doc, fh)
    kinds = {}
    for h in hs:
        kinds[h["type"]] = kinds.get(h["type"], 0) + 1
    print(f"   {len(hs)} faces ({', '.join(f'{k} {n}' for k, n in kinds.items())}) -> {hf}", flush=True)
    if os.path.abspath(first or "") != os.path.abspath(p1):
        shutil.copyfile(out, p1)
    print(f"== refit_contours: 3/3 conversion with them -> {out} (the first kept as {p1})", flush=True)
    if run(["--hints", hf]) != 0:
        return 1
    doc["pass2"] = stamp()
    with open(hf, "w", encoding="utf-8") as fh:
        json.dump(doc, fh)
    return 0


if __name__ == "__main__":
    sys.exit(main())
