# -*- coding: utf-8 -*-
"""
viewer/server.py — local viewer for refit.py
=============================================================================
A small web app (Python standard library + three.js in the page) that:

  * lists the .stl / .step files in test/ and next to refit.py (or opens any other one),
  * shows the input mesh and the converted STEP side by side (and the
    <name>_original.step CAD, when there is one), with the faces coloured by
    surface type, the B-Rep edges, the triangles still left from the mesh,
  * launches refit.py with the chosen options and streams its log live.

    python viewer/server.py              a desktop window (pywebview)
    python viewer/server.py part.stl     the same, with part.stl already open
    python viewer/server.py --browser    the same page in the default browser
    python viewer/server.py --serve      only the server (http://127.0.0.1:8765)

The window is a normal Windows window (pip install pywebview): the page is
drawn by the system's WebView2, the same engine Tauri apps use, but there is
no browser around it. Without pywebview it falls back to the browser.

The conversion runs in its own process (the same `python refit.py ...` you'd
type), so the viewer never changes what refit.py does: the command line it
runs is shown in the page and can be copied as is.

Only listens on 127.0.0.1. Nothing is sent anywhere: three.js is the only
thing the page loads from the network (jsdelivr CDN).
"""
import argparse
import json
import os
import re
import shlex
import signal
import struct
import subprocess
import sys
import threading
import time
import uuid
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, unquote, urlparse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import refit  # noqa: E402,F401  (runs the dependency check)

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
REFIT = os.path.join(ROOT, "refit.py")
TESTS = os.path.join(ROOT, "test")
UPLOADS = os.path.join(ROOT, "viewer_uploads")
EXTS = (".stl", ".step", ".stp")

sys.path.insert(0, ROOT)


def _python_exe():
    """The interpreter for refit.py: python.exe even when the viewer runs under
    pythonw.exe (double-clicked .pyw), so the conversion still has a stdout."""
    exe = sys.executable
    if os.path.basename(exe).lower() == "pythonw.exe":
        alt = os.path.join(os.path.dirname(exe), "python.exe")
        if os.path.isfile(alt):
            return alt
    return exe

# =============================================================================
# OpenCascade (through refit.py's own wrappers: works with OCP and pythonocc)
# =============================================================================
_R = None
_R_lock = threading.Lock()
_OCC_lock = threading.Lock()  # OpenCascade holds the GIL anyway: one at a time


def refit_module():
    """refit.py imported lazily (it pulls in OpenCascade: a couple of seconds)."""
    global _R
    with _R_lock:
        if _R is None:
            import refit  # noqa: E402

            _R = refit
    return _R


# =============================================================================
# STL: read directly with numpy (no OpenCascade, no repair: what's in the file)
# =============================================================================
def read_stl_triangles(path):
    with open(path, "rb") as f:
        data = f.read()
    if len(data) >= 84:
        n = struct.unpack_from("<I", data, 80)[0]
        if 84 + 50 * n == len(data):
            rec = np.frombuffer(data, dtype=np.dtype([("n", "<f4", 3), ("v", "<f4", 9), ("a", "<u2")]), count=n, offset=84)
            return np.ascontiguousarray(rec["v"]).reshape(-1, 3, 3).astype(np.float32)
    txt = data.decode("latin-1")
    nums = re.findall(r"vertex\s+(\S+)\s+(\S+)\s+(\S+)", txt)
    V = np.array(nums, dtype=np.float32).reshape(-1, 3, 3)
    return V


def stl_payload(path):
    tris = read_stl_triangles(path)
    pos = tris.reshape(-1, 3)
    lo, hi = pos.min(0), pos.max(0)
    a = tris[:, 1] - tris[:, 0]
    b = tris[:, 2] - tris[:, 0]
    area = float(0.5 * np.linalg.norm(np.cross(a, b), axis=1).sum())
    # unique vertices / edges / boundary edges: same counts as a mesh tool
    q = np.round(pos.astype(np.float64), 6)
    _, inv = np.unique(q, axis=0, return_inverse=True)
    inv = inv.reshape(-1, 3)
    E = np.sort(np.concatenate([inv[:, [0, 1]], inv[:, [1, 2]], inv[:, [2, 0]]]), axis=1)
    _, cnt = np.unique(E, axis=0, return_counts=True)
    header = {
        "kind": "stl",
        "name": os.path.basename(path),
        "bbox": [lo.tolist(), hi.tolist()],
        "stats": {
            "triangles": int(len(tris)),
            "vertices": int(inv.max() + 1) if len(inv) else 0,
            "edges": int(len(cnt)),
            "boundaryEdges": int((cnt == 1).sum()),
            "nonManifoldEdges": int((cnt > 2).sum()),
            "area": area,
        },
    }
    return header, {"pos": pos.astype(np.float32)}


# =============================================================================
# STEP: faces (triangulated, one block of vertices per face) + B-Rep edges
# =============================================================================
SURF_NAMES = {
    "GeomAbs_Plane": "plane",
    "GeomAbs_Cylinder": "cylinder",
    "GeomAbs_Cone": "cone",
    "GeomAbs_Sphere": "sphere",
    "GeomAbs_Torus": "torus",
    "GeomAbs_BezierSurface": "bspline",
    "GeomAbs_BSplineSurface": "bspline",
    "GeomAbs_SurfaceOfRevolution": "revolution",
    "GeomAbs_SurfaceOfExtrusion": "extrusion",
    "GeomAbs_OffsetSurface": "offset",
}
CURVE_NAMES = {
    "GeomAbs_Line": "line",
    "GeomAbs_Circle": "circle",
    "GeomAbs_Ellipse": "ellipse",
    "GeomAbs_Hyperbola": "hyperbola",
    "GeomAbs_Parabola": "parabola",
    "GeomAbs_BezierCurve": "bspline",
    "GeomAbs_BSplineCurve": "bspline",
    "GeomAbs_OffsetCurve": "offset",
}


def _enum_name(v):
    s = getattr(v, "name", None) or str(v)
    return s.split(".")[-1]


def _surface_info(R, face):
    s = R.BRepAdaptor_Surface(face)
    kind = SURF_NAMES.get(_enum_name(s.GetType()), "other")
    p = {}
    try:
        if kind == "cylinder":
            p["radius"] = s.Cylinder().Radius()
        elif kind == "cone":
            c = s.Cone()
            p["halfAngle"] = abs(np.degrees(c.SemiAngle()))
            p["refRadius"] = c.RefRadius()
        elif kind == "sphere":
            p["radius"] = s.Sphere().Radius()
        elif kind == "torus":
            t = s.Torus()
            p["majorRadius"] = t.MajorRadius()
            p["minorRadius"] = t.MinorRadius()
        elif kind == "bspline":
            p["degree"] = [s.UDegree(), s.VDegree()]
            p["poles"] = [s.NbUPoles(), s.NbVPoles()]
    except Exception:
        pass
    return kind, p


def _transform(P, loc):
    if loc.IsIdentity():
        return P
    tr = loc.Transformation()
    M = np.array([[tr.Value(i, j) for j in range(1, 5)] for i in range(1, 4)])
    return P @ M[:, :3].T + M[:, 3]


def step_payload(path):
    R = refit_module()
    GCPnts = R._m("GCPnts")
    with _OCC_lock:
        r = R.STEPControl_Reader()
        if r.ReadFile(path) != R.IFSelect_RetDone:
            raise RuntimeError(f"STEP read failed: {path}")
        r.TransferRoots()
        shape = r.OneShape()
        if shape.IsNull():
            raise RuntimeError("the STEP file contains no geometry")

        # bounding box from the vertices -> deflection proportional to the part
        V = np.array([[q.X(), q.Y(), q.Z()] for q in (R.bt_Pnt(R.td_Vertex(v)) for v in R.explore(shape, R.TopAbs_VERTEX))])
        if not len(V):
            raise RuntimeError("the STEP file has no vertices")
        lo, hi = V.min(0), V.max(0)
        diag = float(np.linalg.norm(hi - lo)) or 1.0
        lin = max(diag * 4e-4, 1e-4)
        # angular deflection: what shows when zoomed in is the angle between two
        # chords, not their sag - at 0.3 rad a fillet's smooth B-spline rim
        # read as a polyline (test13); 0.1 on faces, 0.05 (~3 deg) on edges
        R._BRepMesh.BRepMesh_IncrementalMesh(shape, lin, False, 0.1, True)

        fmap = R.TopTools_IndexedMapOfShape()
        R.te_MapShapes(shape, R.TopAbs_FACE, fmap)
        nf = R._size(fmap)
        pos_l, idx_l, faces = [], [], []
        tess = []  # (points, triangles) per face, as drawn: the hand fixes pick on them
        nv = nt = 0
        for i in range(1, nf + 1):
            f = R.td_Face(fmap.FindKey(i))
            kind, params = _surface_info(R, f)
            n_edges = len(R.explore(f, R.TopAbs_EDGE))
            n_wires = len(R.explore(f, R.TopAbs_WIRE))
            loc = R.TopLoc_Location()
            tri = R.bt_Triangulation(f, loc)
            P = np.zeros((0, 3))
            T = np.zeros((0, 3), dtype=np.int64)
            if tri is not None and tri.NbTriangles():
                P, T = R._tri_arrays(tri)
                P = _transform(P, loc)
                if f.Orientation() == R.TopAbs_REVERSED:
                    T = T[:, [0, 2, 1]]
            area = 0.0
            if len(T):
                A, B, C = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
                area = float(0.5 * np.linalg.norm(np.cross(B - A, C - A), axis=1).sum())
            # a planar triangle with a single wire is a facet left from the mesh
            facet = kind == "plane" and n_edges == 3 and n_wires == 1
            faces.append({
                "type": kind, "p": {k: (round(v, 6) if isinstance(v, float) else v) for k, v in params.items()},
                "area": area, "edges": n_edges, "wires": n_wires, "facet": facet,
                "v0": nv, "nv": len(P), "t0": nt, "nt": len(T),
            })
            pos_l.append(P.astype(np.float32))
            idx_l.append((T + nv).astype(np.uint32))
            tess.append((np.asarray(P, dtype=float), np.asarray(T, dtype=np.int64)))
            nv += len(P)
            nt += len(T)

        # edges: one polyline each, flattened into segments
        emap = R.TopTools_IndexedMapOfShape()
        R.te_MapShapes(shape, R.TopAbs_EDGE, emap)
        anc = R.edge_face_map(shape)
        ne = R._size(emap)
        seg_l, edges = [], []
        ns = 0
        etypes = {}
        free = 0
        for i in range(1, ne + 1):
            e = R.td_Edge(emap.FindKey(i))
            if R.bt_Degenerated(e):
                continue
            try:
                nfa = len(R._iter_list(anc.FindFromKey(e)))
            except Exception:
                nfa = 2
            if nfa == 1:
                free += 1
            try:
                c = R.BRepAdaptor_Curve(e)
                kind = CURVE_NAMES.get(_enum_name(c.GetType()), "other")
                if kind == "line":
                    a, b = c.Value(c.FirstParameter()), c.Value(c.LastParameter())
                    pts = np.array([[a.X(), a.Y(), a.Z()], [b.X(), b.Y(), b.Z()]])
                else:
                    d = GCPnts.GCPnts_TangentialDeflection(c, 0.05, lin, 2)
                    pts = np.array([[q.X(), q.Y(), q.Z()] for q in (d.Value(k) for k in range(1, d.NbPoints() + 1))])
            except Exception:
                continue
            if len(pts) < 2:
                continue
            etypes[kind] = etypes.get(kind, 0) + 1
            seg = np.stack([pts[:-1], pts[1:]], 1).reshape(-1, 3)
            edges.append([ns, len(pts) - 1, kind, 1 if nfa == 1 else 0])
            seg_l.append(seg.astype(np.float32))
            ns += len(pts) - 1

    stypes = {}
    for f in faces:
        k = "facet" if f["facet"] else f["type"]
        stypes[k] = stypes.get(k, 0) + 1
    header = {
        "kind": "step",
        "name": os.path.basename(path),
        "bbox": [lo.tolist(), hi.tolist()],
        "faces": faces,
        "edges": edges,
        "stats": {
            "faces": nf,
            "edges": len(edges),
            "vertices": len(V),
            "freeEdges": free,
            "faceTypes": stypes,
            "edgeTypes": etypes,
            "area": float(sum(f["area"] for f in faces)),
            "triangles": nt,
        },
    }
    bufs = {
        "pos": np.concatenate(pos_l) if pos_l else np.zeros((0, 3), np.float32),
        "idx": np.concatenate(idx_l) if idx_l else np.zeros((0, 3), np.uint32),
        "seg": np.concatenate(seg_l) if seg_l else np.zeros((0, 3), np.float32),
    }
    with _info_lock:
        if len(_STEP_INFO) > 6:
            _STEP_INFO.pop(next(iter(_STEP_INFO)))
        _STEP_INFO[_file_key(path)] = {
            "shape": shape, "fmap": fmap, "tess": tess, "faces": faces, "diag": diag,
        }
    return header, bufs


def pack(header, bufs):
    """[u32 header length][header JSON][pad to 4][buffers, each 4-aligned]."""
    blobs, off = [], 0
    header["buffers"] = {}
    for name, arr in bufs.items():
        arr = np.ascontiguousarray(arr)
        b = arr.tobytes()
        header["buffers"][name] = {"offset": off, "length": int(arr.size), "dtype": str(arr.dtype)}
        blobs.append(b)
        off += len(b)
    hj = json.dumps(header, separators=(",", ":")).encode("utf-8")
    hj += b" " * ((-(4 + len(hj))) % 4)
    return struct.pack("<I", len(hj)) + hj + b"".join(blobs)


_cache = {}
_cache_lock = threading.Lock()


def mesh_bytes(path):
    st = os.stat(path)
    key = (os.path.normcase(os.path.abspath(path)), st.st_mtime_ns, st.st_size)
    with _cache_lock:
        if key in _cache:
            return _cache[key]
    t0 = time.perf_counter()
    if path.lower().endswith(".stl"):
        header, bufs = stl_payload(path)
    else:
        header, bufs = step_payload(path)
    header["loadTime"] = time.perf_counter() - t0
    header["path"] = os.path.abspath(path)
    header["mtime"] = st.st_mtime
    data = pack(header, bufs)
    with _cache_lock:
        if len(_cache) > 12:
            _cache.pop(next(iter(_cache)))
        _cache[key] = data
    return data


# =============================================================================
# Hand fixes: pick faces of the converted STEP, say what they are
# =============================================================================
# Some rounds come out in pieces (strips of cylinders and tori, little
# B-splines, triangles left between them). A click on one picks the piece and
# grows across smooth edges into the pieces like it; the user says what the
# whole is (auto, plane, cylinder, cone, sphere, torus, tube, free-form) and
# the input's triangles under those faces go to <name>.hints.json, one point
# per triangle. refit.py --hints then converts each group as ONE face of that
# type (see "HAND FIXES" in refit.py). The file is the mesh's, so groups
# picked on one conversion still hold on the next.
_STEP_INFO = {}
_STL_INFO = {}
_info_lock = threading.Lock()
HINT_TYPES = ("auto", "plane", "cylinder", "cone", "sphere", "torus", "tube", "freeform")


def _file_key(path):
    st = os.stat(path)
    return (os.path.normcase(os.path.abspath(path)), st.st_mtime_ns, st.st_size)


def step_info(path):
    """The STEP as the page shows it: shape, faces (same order), tessellation."""
    key = _file_key(path)
    with _info_lock:
        info = _STEP_INFO.get(key)
    if info is None:
        step_payload(path)
        with _info_lock:
            info = _STEP_INFO.get(key)
    if info is None:
        raise RuntimeError(f"could not read {path}")
    return info


def stl_info(path):
    """The input's triangles (float64), their centroids and a KD-tree on them."""
    from scipy.spatial import cKDTree

    key = _file_key(path)
    with _info_lock:
        info = _STL_INFO.get(key)
    if info is None:
        tris = read_stl_triangles(path).astype(np.float64)
        cen = tris.mean(axis=1)
        lo, hi = tris.reshape(-1, 3).min(0), tris.reshape(-1, 3).max(0)
        info = {"tris": tris, "cen": cen, "tree": cKDTree(cen), "diag": float(np.linalg.norm(hi - lo)) or 1.0, "t2f": {}}
        with _info_lock:
            if len(_STL_INFO) > 4:
                _STL_INFO.pop(next(iter(_STL_INFO)))
            _STL_INFO[key] = info
    return info


def _close_pick(si, pick):
    """
    The pick (a mask over the input's triangles) with the few triangles it
    encloses - a triangle of another face left inside a round would make it
    a ring, cut in two (refit.py: _hint_holes, the same rule): a connected
    piece of the triangles outside the pick, not the largest, at most 5% of
    the pick's triangles and 2% of its area.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    if "edges" not in si:
        V, inv = np.unique(np.round(si["tris"].reshape(-1, 3), 9), axis=0, return_inverse=True)
        T = inv.reshape(-1, 3)
        E = np.sort(np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]]), axis=1)
        u, ei = np.unique(E, axis=0, return_inverse=True)
        ei = ei.ravel()
        tid = np.tile(np.arange(len(T)), 3)
        M = coo_matrix((np.ones(len(tid)), (tid, ei)), shape=(len(T), len(u))).tocsr()
        t = si["tris"]
        area = 0.5 * np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1)
        si["edges"] = ((M @ M.T).tocsr(), area)
    A, area = si["edges"]
    out = ~pick
    if not out.any() or not pick.any():
        return pick
    idx = np.nonzero(out)[0]
    lab = connected_components(A[idx][:, idx], directed=False)[1]
    size = np.bincount(lab)
    a_c = np.bincount(lab, weights=area[idx])
    ok = (size <= 0.05 * pick.sum()) & (a_c <= 0.02 * area[pick].sum()) & (np.arange(len(size)) != np.argmax(size))
    if not ok.any():
        return pick
    res = pick.copy()
    res[idx[ok[lab]]] = True
    return res


def _tri_normals(P, T):
    A, B, C = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
    n = np.cross(B - A, C - A)
    a = np.linalg.norm(n, axis=1)
    return (A + B + C) / 3.0, n / np.maximum(a, 1e-300)[:, None], 0.5 * a


def face_graph(info):
    """
    Neighbours of every face across its edges, with the angle between the two
    faces there (their tessellations' normals next to the edge's middle: 0 =
    tangent), and every face's curvature (1/radius: from its surface for the
    quadrics, from its own tessellation for the rest, None for a lone
    triangle).
    """
    if "graph" in info:
        return info["graph"]
    from scipy.spatial import cKDTree

    R = refit_module()
    tess, faces, fmap, shape = info["tess"], info["faces"], info["fmap"], info["shape"]
    nf = len(faces)
    C_l, N_l, O_l = [], [], []
    kappa = [None] * nf
    for i, (P, T) in enumerate(tess):
        if not len(T):
            continue
        c, n, a = _tri_normals(P, T)
        ok = a > 1e-18
        C_l.append(c[ok]); N_l.append(n[ok]); O_l.append(np.full(int(ok.sum()), i))
        f, p = faces[i], faces[i]["p"]
        kind = f["type"]
        if kind == "plane":
            kappa[i] = None if f["facet"] else 0.0
        elif kind in ("cylinder", "sphere") and p.get("radius"):
            kappa[i] = 1.0 / p["radius"]
        elif kind == "torus" and p.get("minorRadius"):
            kappa[i] = 1.0 / p["minorRadius"]
        elif ok.sum() >= 2:
            # the turn between neighbouring triangles over the distance of their centres
            E = {}
            turns = []
            for t, tri in enumerate(T):
                if not ok[t]:
                    continue
                for a_, b_ in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
                    k_ = (min(a_, b_), max(a_, b_))
                    if k_ in E:
                        s = E.pop(k_)
                        d = float(np.linalg.norm(c[t] - c[s]))
                        if d > 1e-12:
                            turns.append(float(np.arccos(np.clip(n[t] @ n[s], -1.0, 1.0))) / d)
                    else:
                        E[k_] = t
            kappa[i] = float(np.median(turns)) if turns else 0.0
    TC = np.concatenate(C_l) if C_l else np.zeros((0, 3))
    TN = np.concatenate(N_l) if N_l else np.zeros((0, 3))
    TO = np.concatenate(O_l) if O_l else np.zeros(0, int)
    tree = cKDTree(TC) if len(TC) else None
    cent = [tess[i][0].mean(axis=0) if len(tess[i][0]) else np.zeros(3) for i in range(nf)]
    nbr = [dict() for _ in range(nf)]
    with _OCC_lock:
        em = R.edge_face_map(shape)
        for k in range(1, R._size(em) + 1):
            e = R.td_Edge(em.FindKey(k))
            if R.bt_Degenerated(e):
                continue
            fl = sorted({fmap.FindIndex(g) - 1 for g in R._iter_list(em.FindFromIndex(k))})
            fl = [x for x in fl if x >= 0]
            if len(fl) != 2 or tree is None:
                continue
            a, b = fl
            try:
                c = R.BRepAdaptor_Curve(e)
                q = c.Value(0.5 * (c.FirstParameter() + c.LastParameter()))
                m = np.array([q.X(), q.Y(), q.Z()])
            except Exception:
                continue
            _, idx = tree.query(m, k=min(64, len(TC)))
            idx = np.atleast_1d(idx)
            na = next((TN[j] for j in idx if TO[j] == a), None)
            nb = next((TN[j] for j in idx if TO[j] == b), None)
            ang = 90.0 if na is None or nb is None else float(np.degrees(np.arccos(np.clip(na @ nb, -1.0, 1.0))))
            # (several edges between the same two faces: the sharpest counts)
            ang = max(ang, nbr[a].get(b, 0.0))
            nbr[a][b] = nbr[b][a] = ang
    # a lone triangle: the turn to its smooth neighbours over the distance
    for i in range(nf):
        if kappa[i] is None:
            ks = [np.radians(ang) / max(float(np.linalg.norm(cent[i] - cent[j])), 1e-12) for j, ang in nbr[i].items() if ang < 30.0]
            kappa[i] = ("facet", float(np.median(ks)) if ks else 0.0)
    info["graph"] = {"nbr": nbr, "kappa": kappa, "cent": cent}
    return info["graph"]


def grow_selection(info, seed, max_angle=30.0):
    """
    The face clicked and the pieces like it: across edges where the faces are
    (nearly) tangent, into faces of a similar curvature (within 2.5x - a round's
    strips; not the wall or the bore it rounds), the triangles left between
    them, never into a face far bigger than the pieces so far (a face of its
    own, not a piece). From a plane: the pieces of that plane only.
    Returns (faces, label).
    """
    G = face_graph(info)
    faces, nbr, kap = info["faces"], G["nbr"], G["kappa"]
    flat = 3.0 / max(info["diag"], 1e-9)

    def k_of(i):
        k = kap[i]
        return k[1] if isinstance(k, tuple) else k

    def is_facet(i):
        return isinstance(kap[i], tuple)

    ks = k_of(seed)
    planar = ks is None or ks < flat
    sel, stack = {seed}, [seed]
    amax = faces[seed]["area"]
    while stack and len(sel) < 4000:
        f = stack.pop()
        for g, ang in nbr[f].items():
            if g in sel:
                continue
            if ang > (2.0 if planar else max_angle):
                continue
            kg = k_of(g)
            if planar:
                if not is_facet(g) and kg is not None and kg >= flat:
                    continue
            elif not is_facet(g):
                if kg is None or not (ks / 2.5 <= kg <= ks * 2.5):
                    continue
                if faces[g]["area"] > 20.0 * amax:
                    continue
            sel.add(g)
            amax = max(amax, faces[g]["area"])
            stack.append(g)
    n = len(sel)
    label = (f"One round, {n} faces" if not planar else f"One plane, {n} faces") if n > 1 else "1 face"
    return sorted(sel), label


def tri_faces(stl_path, step_path):
    """For every input triangle the output face it lies on (the nearest of the
    STEP's tessellation triangles to its centre)."""
    from scipy.spatial import cKDTree

    si = stl_info(stl_path)
    key = _file_key(step_path)
    if key in si["t2f"]:
        return si["t2f"][key]
    info = step_info(step_path)
    C_l, O_l, A_l, B_l, D_l = [], [], [], [], []
    for i, (P, T) in enumerate(info["tess"]):
        if not len(T):
            continue
        A_l.append(P[T[:, 0]]); B_l.append(P[T[:, 1]]); D_l.append(P[T[:, 2]])
        C_l.append(P[T].mean(axis=1)); O_l.append(np.full(len(T), i))
    if not C_l:
        raise RuntimeError("the STEP has no tessellation")
    TC, TO = np.concatenate(C_l), np.concatenate(O_l)
    TA, TB, TD = np.concatenate(A_l), np.concatenate(B_l), np.concatenate(D_l)
    k = min(12, len(TC))
    _, idx = cKDTree(TC).query(si["cen"], k=k)
    idx = idx.reshape(len(si["cen"]), -1)
    # exact distance to each candidate triangle, the nearest wins
    best = np.full(len(si["cen"]), np.inf)
    out = np.full(len(si["cen"]), -1)
    for c in range(idx.shape[1]):
        j = idx[:, c]
        d = _pt_tri_dist(si["cen"], TA[j], TB[j], TD[j])
        w = d < best
        best[w] = d[w]
        out[w] = TO[j][w]
    si["t2f"] = {key: out}
    return out


def _pt_tri_dist(P, A, B, C):
    """Distance from points P to triangles (A, B, C), row by row."""
    AB, AC, AP = B - A, C - A, P - A
    n = np.cross(AB, AC)
    nn = np.maximum(np.einsum("ij,ij->i", n, n), 1e-300)
    # barycentric of the projection
    d00, d01, d11 = np.einsum("ij,ij->i", AB, AB), np.einsum("ij,ij->i", AB, AC), np.einsum("ij,ij->i", AC, AC)
    d20, d21 = np.einsum("ij,ij->i", AP, AB), np.einsum("ij,ij->i", AP, AC)
    den = np.maximum(d00 * d11 - d01 * d01, 1e-300)
    v = (d11 * d20 - d01 * d21) / den
    w = (d00 * d21 - d01 * d20) / den
    inside = (v >= 0) & (w >= 0) & (v + w <= 1)
    d_plane = np.abs(np.einsum("ij,ij->i", AP, n)) / np.sqrt(nn)

    def seg(P_, X, Y):
        XY = Y - X
        t = np.clip(np.einsum("ij,ij->i", P_ - X, XY) / np.maximum(np.einsum("ij,ij->i", XY, XY), 1e-300), 0.0, 1.0)
        return np.linalg.norm(P_ - (X + t[:, None] * XY), axis=1)

    d_edge = np.minimum(np.minimum(seg(P, A, B), seg(P, B, C)), seg(P, C, A))
    return np.where(inside, d_plane, d_edge)


def hints_path(inp):
    return os.path.splitext(inp)[0] + ".hints.json"


def read_hints(inp):
    p = hints_path(inp)
    if not os.path.isfile(p):
        return {"refit_hints": 1, "input": os.path.basename(inp), "groups": []}
    with open(p, "r", encoding="utf-8") as f:
        d = json.load(f)
    d.setdefault("groups", [])
    return d


def write_hints(inp, d):
    p = hints_path(inp)
    if not d.get("groups"):
        if os.path.isfile(p):
            os.remove(p)
        return None
    tmp = p + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(d, f, separators=(",", ":"))
    os.replace(tmp, p)
    return p


def hints_state(inp, out=None):
    """The groups of the input's hints file, each with the faces of the current
    output it lies on (the faces holding most of its triangles)."""
    d = read_hints(inp)
    t2f = None
    if out and os.path.isfile(out) and out.lower().endswith((".step", ".stp")) and inp.lower().endswith(".stl"):
        try:
            t2f = tri_faces(inp, out)
            si = stl_info(inp)
        except Exception:
            t2f = None
    groups = []
    for g in d["groups"]:
        e = {k: g.get(k) for k in ("type", "label", "result")}
        e["n"] = len(g.get("points", []))
        if t2f is not None and e["n"]:
            dist, ti = si["tree"].query(np.asarray(g["points"], dtype=float))
            fids = t2f[ti[dist < 1e-6 * si["diag"] + 1e-6]]
            # a face counts if most of its triangles are the group's
            cnt = np.bincount(fids[fids >= 0], minlength=1)
            tot = np.bincount(t2f[t2f >= 0], minlength=len(cnt))[: len(cnt)]
            e["faces"] = [int(i) for i in np.nonzero(cnt >= 0.5 * np.maximum(tot, 1))[0]]
        groups.append(e)
    return {"path": hints_path(inp) if d["groups"] else None, "groups": groups}


def _boundary_loops(tris):
    """How many boundary loops the largest connected piece of a set of
    triangles (n, 3, 3) has (refit keeps only that piece of a pick): 1 a
    patch, 2 or more a ring or a patch with holes."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    V, inv = np.unique(np.round(tris.reshape(-1, 3), 6), axis=0, return_inverse=True)
    T = inv.reshape(-1, 3)
    E = np.sort(np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]]), axis=1)
    u, ei, c = np.unique(E, axis=0, return_inverse=True, return_counts=True)
    ei = ei.ravel()
    # (triangles joined through shared edges: the largest piece)
    tid = np.tile(np.arange(len(T)), 3)
    M = coo_matrix((np.ones(len(tid)), (tid, ei)), shape=(len(T), len(u))).tocsr()
    lab_t = connected_components(M @ M.T, directed=False)[1]
    big = np.argmax(np.bincount(lab_t))
    keep = lab_t[tid] == big
    e_in, n_in = np.unique(ei[keep], return_counts=True)
    b = u[e_in[n_in == 1]]
    if not len(b):
        return 0
    g = coo_matrix((np.ones(len(b)), (b[:, 0], b[:, 1])), shape=(len(V), len(V)))
    lab = connected_components(g, directed=False)[1]
    return len(np.unique(lab[np.unique(b)]))


def _fit_preview(tris, kind, diag):
    """What the type the user picked makes of those triangles - before any
    Build: the quadric's size and how far the triangles' corners are from it."""
    R = refit_module()
    names = {"plane": "plane", "cylinder": "cylinder", "cone": "cone", "sphere": "sphere", "torus": "torus", "tube": "tube", "freeform": "free-form face", "auto": "face"}

    def sheets(k):
        # (a closed ring doesn't lay flat: refit cuts it in two - see fit_hint)
        many = {"tube": "tubes", "freeform": "free-form faces"}[k]
        return f"two {many} (a ring, or a patch with a hole, is cut in two)" if _boundary_loops(tris) >= 2 else f"one {names[k]}"

    if kind in ("tube", "freeform"):
        return {"kind": kind, "text": sheets(kind)}
    verts = [t for t in tris]
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    a = np.linalg.norm(n, axis=1)
    norms = n / np.maximum(a, 1e-300)[:, None]
    areas = 0.5 * a
    tol = R.auto_tol(diag) if hasattr(R, "auto_tol") else max(2e-4, 1e-5 * diag)
    with _OCC_lock:
        p = R.fit_hint_quadric(kind, verts, norms, areas, list(range(len(tris))), tol)
    if p is None:
        return {"kind": "freeform" if kind == "auto" else kind, "text": sheets("freeform") if kind == "auto" else f"no {names[kind]} fits"}
    P = tris.reshape(-1, 3)
    dev = float(np.abs(p.dist(P)).max())
    lab = p.label().lower()
    lab = {"cyl": "cylinder", "free": "free-form face"}.get(lab, lab)
    size = ""
    if p.kind == R.TORUS:
        size = f"R {p.r0:.4g} r {p.r1:.4g}"
    elif p.kind == R.SPHERE or (p.kind == R.AXIAL and abs(p.slope) < 1e-12):
        size = f"r {p.r0:.4g}"
    elif p.kind == R.AXIAL:
        size = f"half-angle {np.degrees(np.arctan(abs(p.slope))):.3g}°"
    # (refit refuses a hand fix standing further off than this - HINT_MAX_DEV
    # edge tolerance ceilings, Phase C's defaults)
    cap = R.HINT_MAX_DEV * max(20.0 * tol, 2e-4 * diag)
    return {"kind": lab, "text": f"one {lab}" + (f" · {size}" if size else "") + f" · within {dev:.2g} mm", "dev": dev, "cap": cap, "far": dev > cap}


def hint_add(inp, out, faces, kind):
    if kind not in HINT_TYPES:
        raise ValueError(f"unknown type {kind}")
    t2f = tri_faces(inp, out)
    si = stl_info(inp)
    pick = np.isin(t2f, np.asarray(faces, dtype=int))
    if not pick.any():
        raise ValueError("no input triangle under those faces")
    pick = _close_pick(si, pick)
    pts = np.round(si["cen"][pick], 6)
    t = si["tris"][pick]
    # (each point's triangle area: refit takes a merged facet only when the
    # pick covers half of it - see _hint_owner)
    areas = np.round(0.5 * np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1), 9)
    pv = _fit_preview(t, kind, si["diag"])
    d = read_hints(inp)
    # (a triangle belongs to one group: the newest)
    newset = set(map(tuple, pts.tolist()))
    for g in d["groups"]:
        keep = [tuple(q) not in newset for q in g.get("points", [])]
        g["points"] = [q for q, k in zip(g.get("points", []), keep) if k]
        if len(g.get("areas") or []) == len(keep):
            g["areas"] = [a for a, k in zip(g["areas"], keep) if k]
        else:
            g.pop("areas", None)
    d["groups"] = [g for g in d["groups"] if g.get("points")]
    label = f"{len(faces)} face{'s' if len(faces) != 1 else ''} → {pv['text']}"
    d["groups"].append({"type": kind, "label": label, "faces_picked": len(faces), "points": pts.tolist(), "areas": areas.tolist()})
    d["input"] = os.path.basename(inp)
    write_hints(inp, d)
    st = hints_state(inp, out)
    st["preview"] = pv
    st["label"] = label
    return st


# =============================================================================
# Jobs: refit.py in its own process, the log kept line by line
# =============================================================================
class Job:
    def __init__(self, cmd, inp, out):
        self.id = uuid.uuid4().hex[:10]
        self.cmd = cmd
        self.input = inp
        self.output = out
        self.lines = []
        self.rc = None
        self.t0 = time.time()
        self.t1 = None
        self.cancelled = False
        env = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONUNBUFFERED="1")
        kw = {}
        if os.name == "nt":
            kw["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP | getattr(subprocess, "CREATE_NO_WINDOW", 0)
        else:
            kw["start_new_session"] = True
        self.proc = subprocess.Popen(cmd, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, env=env, **kw)
        threading.Thread(target=self._pump, daemon=True).start()

    def _pump(self):
        for raw in iter(self.proc.stdout.readline, b""):
            self.lines.append(raw.decode("utf-8", "replace").rstrip("\r\n"))
        self.proc.wait()
        self.rc = self.proc.returncode
        self.t1 = time.time()

    def cancel(self):
        if self.rc is not None:
            return
        self.cancelled = True
        # the Phase B/C worker processes are children: kill the whole tree
        if os.name == "nt":
            subprocess.run(["taskkill", "/F", "/T", "/PID", str(self.proc.pid)], capture_output=True)
        else:
            try:
                os.killpg(self.proc.pid, signal.SIGTERM)
            except Exception:
                self.proc.terminate()

    def state(self, since=0):
        out_ok = self.rc == 0 and os.path.isfile(self.output)
        rep = os.path.splitext(self.output)[0] + "_report.txt"
        return {
            "id": self.id,
            "cmd": self.cmd,
            "lines": self.lines[since:],
            "next": len(self.lines),
            "running": self.rc is None,
            "rc": self.rc,
            "cancelled": self.cancelled,
            "elapsed": (self.t1 or time.time()) - self.t0,
            "output": self.output if out_ok else None,
            "report": rep if (out_ok and os.path.isfile(rep)) else None,
        }


JOBS = {}


def source_mesh(step_path):
    """The mesh a STEP was converted from: <name>.stl beside it (any case of
    the extension; <name>_refit / _phaseA.step count as <name>), None for a
    <name>_original.step - the CAD, not a conversion."""
    d, b = os.path.split(os.path.splitext(step_path)[0])
    if b.endswith("_original"):
        return None
    for suf in ("_refit", "_phaseA"):
        if b.endswith(suf):
            b = b[: -len(suf)]
    try:
        names = os.listdir(d or ".")
    except OSError:
        return None
    for n in names:
        s, e = os.path.splitext(n)
        if e.lower() == ".stl" and os.path.normcase(s) == os.path.normcase(b) and os.path.isfile(os.path.join(d, n)):
            return os.path.join(d, n)
    return None


def default_output(inp, seq_args):
    stem = os.path.splitext(inp)[0]
    phases = [a for a in seq_args if a in ("-a", "-b", "-c")]
    only_a = bool(phases) and all(a == "-a" for a in phases)
    out = f"{stem}_phaseA.step" if only_a else f"{stem}.step"
    if os.path.normcase(os.path.abspath(out)) == os.path.normcase(os.path.abspath(inp)):
        out = f"{stem}_refit.step"
    return out


def resolve(p):
    if not p:
        raise ValueError("missing path")
    p = unquote(p)
    if not os.path.isabs(p):
        p = os.path.join(ROOT, p)
    return os.path.abspath(p)


def list_files():
    out = []
    for d in (ROOT, TESTS, UPLOADS):
        if not os.path.isdir(d):
            continue
        for n in sorted(os.listdir(d), key=lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]):
            fp = os.path.join(d, n)
            if os.path.isfile(fp) and n.lower().endswith(EXTS):
                st = os.stat(fp)
                out.append({
                    "name": n if d == ROOT else os.path.relpath(fp, ROOT).replace("\\", "/"),
                    "path": fp,
                    "size": st.st_size,
                    "mtime": st.st_mtime,
                })
    return out


# =============================================================================
# HTTP
# =============================================================================
class Handler(BaseHTTPRequestHandler):
    server_version = "refit-viewer/1.0"

    def log_message(self, fmt, *a):
        pass

    def _send(self, code, body, ctype="application/json; charset=utf-8"):
        if isinstance(body, (dict, list)):
            body = json.dumps(body).encode("utf-8")
        elif isinstance(body, str):
            body = body.encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _json_body(self):
        n = int(self.headers.get("Content-Length") or 0)
        return json.loads(self.rfile.read(n) or b"{}")

    def do_GET(self):
        u = urlparse(self.path)
        q = {k: v[0] for k, v in parse_qs(u.query).items()}
        try:
            if u.path in ("/", "/index.html"):
                with open(os.path.join(HERE, "index.html"), "rb") as f:
                    return self._send(200, f.read(), "text/html; charset=utf-8")
            if u.path == "/api/info":
                return self._send(200, {"root": ROOT, "python": _python_exe(), "refit": REFIT, "cpus": os.cpu_count()})
            if u.path == "/api/exists":
                p = resolve(q.get("path"))
                return self._send(200, {"exists": os.path.isfile(p) and p.lower().endswith(EXTS)})
            if u.path == "/api/files":
                return self._send(200, list_files())
            if u.path == "/api/related":
                # the output refit.py would write for this input (if it's already
                # there) and the <name>_original.step CAD, if there is one
                p = resolve(q.get("path"))
                args = shlex.split(q.get("args", ""), posix=True)
                out = resolve(q["output"]) if q.get("output") else default_output(p, args)
                stem = os.path.splitext(p)[0]
                for suf in ("_refit", "_phaseA", "_original"):
                    if stem.endswith(suf):
                        stem = stem[: -len(suf)]
                orig = stem + "_original.step"
                return self._send(200, {
                    "output": out if (os.path.isfile(out) and out != p) else None,
                    "original": orig if os.path.isfile(orig) else None,
                    "mesh": source_mesh(p) if p.lower().endswith((".step", ".stp")) else None,
                })
            if u.path == "/api/default-output":
                inp = resolve(q.get("input"))
                args = shlex.split(q.get("args", ""), posix=True)
                out = resolve(q["output"]) if q.get("output") else default_output(inp, args)
                return self._send(200, {"output": out, "exists": os.path.isfile(out)})
            if u.path == "/api/mesh":
                p = resolve(q.get("path"))
                if not p.lower().endswith(EXTS) or not os.path.isfile(p):
                    return self._send(404, {"error": f"not found: {p}"})
                return self._send(200, mesh_bytes(p), "application/octet-stream")
            if u.path == "/api/text":
                p = resolve(q.get("path"))
                if not p.lower().endswith(".txt") or not os.path.isfile(p):
                    return self._send(404, {"error": "not found"})
                with open(p, "r", encoding="utf-8", errors="replace") as f:
                    return self._send(200, f.read(), "text/plain; charset=utf-8")
            if u.path == "/api/job":
                j = JOBS.get(q.get("id"))
                if not j:
                    return self._send(404, {"error": "no such job"})
                return self._send(200, j.state(int(q.get("since", 0))))
            if u.path == "/api/hints":
                inp = resolve(q.get("input"))
                out = resolve(q["output"]) if q.get("output") else None
                return self._send(200, hints_state(inp, out))
            return self._send(404, {"error": "not found"})
        except Exception as ex:  # noqa: BLE001
            return self._send(500, {"error": f"{type(ex).__name__}: {ex}"})

    def do_POST(self):
        u = urlparse(self.path)
        q = {k: v[0] for k, v in parse_qs(u.query).items()}
        try:
            if u.path == "/api/upload":
                name = os.path.basename(q.get("name", "upload.stl"))
                if not name.lower().endswith(EXTS):
                    return self._send(400, {"error": "only .stl / .step / .stp"})
                os.makedirs(UPLOADS, exist_ok=True)
                n = int(self.headers.get("Content-Length") or 0)
                dst = os.path.join(UPLOADS, name)
                with open(dst, "wb") as f:
                    left = n
                    while left > 0:
                        chunk = self.rfile.read(min(left, 1 << 20))
                        if not chunk:
                            break
                        f.write(chunk)
                        left -= len(chunk)
                return self._send(200, {"path": dst})
            if u.path == "/api/run":
                body = self._json_body()
                inp = resolve(body.get("input"))
                if not os.path.isfile(inp):
                    return self._send(400, {"error": f"input not found: {inp}"})
                args = shlex.split(body.get("args", ""), posix=True)
                out = body.get("output") or default_output(inp, args)
                out = resolve(out)
                if not out.lower().endswith((".step", ".stp")):
                    return self._send(400, {"error": "the output must be a .step / .stp file"})
                if any(j.rc is None for j in JOBS.values()):
                    return self._send(409, {"error": "a conversion is already running"})
                # (the hand fixes picked in the page, if asked and there are any)
                if body.get("hints") and os.path.isfile(hints_path(inp)) and "--hints" not in args:
                    args += ["--hints", hints_path(inp)]
                cmd = [_python_exe(), "-u", REFIT, inp, *args, "-o", out, "--no-color"]
                j = Job(cmd, inp, out)
                JOBS[j.id] = j
                return self._send(200, {"id": j.id, "cmd": cmd, "output": out})
            if u.path == "/api/cancel":
                j = JOBS.get(self._json_body().get("id"))
                if j:
                    j.cancel()
                return self._send(200, {"ok": True})
            if u.path == "/api/select":
                # the face clicked and the pieces like it (see grow_selection)
                body = self._json_body()
                out = resolve(body.get("output"))
                faces, label = grow_selection(step_info(out), int(body["face"]))
                return self._send(200, {"faces": faces, "label": label})
            if u.path in ("/api/hints/add", "/api/hints/remove", "/api/hints/clear"):
                body = self._json_body()
                inp = resolve(body.get("input"))
                out = resolve(body["output"]) if body.get("output") else None
                if not inp.lower().endswith(".stl"):
                    return self._send(400, {"error": "hand fixes need an .stl input"})
                if u.path == "/api/hints/add":
                    if not out or not os.path.isfile(out):
                        return self._send(400, {"error": "convert it first: the faces are picked on the STEP"})
                    return self._send(200, hint_add(inp, out, [int(x) for x in body.get("faces", [])], str(body.get("type", "auto"))))
                d = read_hints(inp)
                if u.path == "/api/hints/remove":
                    k = int(body.get("index", -1))
                    if 0 <= k < len(d["groups"]):
                        d["groups"].pop(k)
                else:
                    d["groups"] = []
                write_hints(inp, d)
                return self._send(200, hints_state(inp, out))
            return self._send(404, {"error": "not found"})
        except Exception as ex:  # noqa: BLE001
            return self._send(500, {"error": f"{type(ex).__name__}: {ex}"})


# =============================================================================
# Desktop window (pywebview)
# =============================================================================
class NativeApi:
    """What the page can ask the window for (window.pywebview.api.*)."""

    def __init__(self):
        # ⚠️ private: pywebview exposes every PUBLIC attribute to the page and
        # walks into it - through the window it ends up in the .NET objects
        self._window = None

    def open_file(self):
        import webview

        kind = getattr(getattr(webview, "FileDialog", None), "OPEN", None) or webview.OPEN_DIALOG
        # ⚠️ only letters, digits and spaces before the parenthesis: with "Mesh / STEP"
        # pywebview raised ValueError and the dialog never opened
        res = self._window.create_file_dialog(kind, directory=TESTS if os.path.isdir(TESTS) else ROOT, allow_multiple=False,
                                             file_types=("Mesh and STEP (*.stl;*.step;*.stp)", "All files (*.*)"))
        if not res:
            return None
        return res[0] if isinstance(res, (list, tuple)) else res

    def save_file(self, suggested):
        import webview

        kind = getattr(getattr(webview, "FileDialog", None), "SAVE", None) or webview.SAVE_DIALOG
        d, n = os.path.split(suggested or "")
        res = self._window.create_file_dialog(kind, directory=d or ROOT, save_filename=n or "out.step",
                                             file_types=("STEP (*.step;*.stp)",))
        if not res:
            return None
        return res[0] if isinstance(res, (list, tuple)) else res

    def reveal(self, path):
        if os.name == "nt" and path and os.path.exists(path):
            subprocess.Popen(["explorer", "/select,", os.path.normpath(path)])

    def set_title(self, title):
        self._window.set_title(title)


def _dark_title_bar(window):
    """Windows 10/11: the title bar follows the page's dark theme."""
    if os.name != "nt":
        return
    try:
        import ctypes

        hwnd = window.native.Handle.ToInt32()
        on = ctypes.c_int(1)
        for attr in (20, 19):  # DWMWA_USE_IMMERSIVE_DARK_MODE (and its pre-20H1 number)
            if ctypes.windll.dwmapi.DwmSetWindowAttribute(hwnd, attr, ctypes.byref(on), ctypes.sizeof(on)) == 0:
                break
    except Exception:
        pass


def run_window(url):
    import webview
    from webview.menu import Menu, MenuAction, MenuSeparator

    api = NativeApi()
    # The window as the user laid it out: 2547 x 1495 at (6, 15), in the
    # display's own pixels (2560 x 1600 at 175%). pywebview takes logical
    # pixels and multiplies them by the scale, so they go in divided by it
    # (+0.5: its int() must land on the very pixel).
    W_PX, H_PX, X_PX, Y_PX = 2547, 1495, 6, 15
    w, h, x, y = 1480, 920, None, None
    try:
        sc = webview.screens[0]  # logical pixels
        try:
            import ctypes

            k = ctypes.windll.shcore.GetScaleFactorForDevice(0) / 100.0
        except Exception:
            k = 1.0
        if X_PX + W_PX <= (sc.x + sc.width) * k + 1 and Y_PX + H_PX <= (sc.y + sc.height) * k + 1:
            w, h, x, y = ((v + 0.5) / k for v in (W_PX, H_PX, X_PX, Y_PX))
        else:
            # a smaller display: 85% of it, never past it, centred a bit
            # above the middle (the taskbar takes the bottom)
            w, h = min(w, int(sc.width * 0.85)), min(h, int(sc.height * 0.85))
            x, y = sc.x + (sc.width - w) // 2, sc.y + max(0, (sc.height - h) // 2 - 24)
    except Exception:
        pass
    win = webview.create_window(
        "Refit", url, js_api=api, width=w, height=h, x=x, y=y, min_size=(min(900, w), min(600, h)),
        background_color="#121418", text_select=True,
    )
    api._window = win
    js = lambda code: (lambda: win.evaluate_js(code))  # noqa: E731
    menu = [
        Menu("File", [
            MenuAction("Open…", js("window.refitUI && refitUI.openFile()")),
            MenuAction("Show in folder", js("window.refitUI && refitUI.reveal()")),
            MenuSeparator(),
            MenuAction("Exit", win.destroy),
        ]),
        Menu("View", [
            MenuAction("Fit (F)", js("refitUI.view(null)")),
            MenuAction("Isometric (I)", js("refitUI.view('iso')")),
            MenuAction("Top", js("refitUI.view('top')")),
            MenuAction("Front", js("refitUI.view('front')")),
            MenuAction("Right", js("refitUI.view('right')")),
            MenuAction("Free / Z-up orbit (O)", js("refitUI.orbit()")),
            Menu("Section", [
                MenuAction("X", js("refitUI.section('x')")),
                MenuAction("Y", js("refitUI.section('y')")),
                MenuAction("Z", js("refitUI.section('z')")),
                MenuAction("Flip", js("refitUI.flipSection()")),
                MenuAction("Off", js("refitUI.section(null)")),
            ]),
            MenuSeparator(),
            MenuAction("Input (1)", js("refitUI.tab('input')")),
            MenuAction("STEP (2)", js("refitUI.tab('output')")),
            MenuAction("Original CAD (3)", js("refitUI.tab('cad')")),
            MenuSeparator(),
            MenuAction("Log (L)", js("refitUI.log()")),
        ]),
        Menu("Help", [
            MenuAction("README", lambda: os.startfile(os.path.join(ROOT, "README.md")) if os.name == "nt" else None),
            MenuAction("Open the refit folder", lambda: subprocess.Popen(["explorer", ROOT]) if os.name == "nt" else None),
        ]),
    ]
    win.events.shown += lambda: _dark_title_bar(win)
    # private_mode=False + a storage folder: the page remembers the last file
    # and the options between sessions, as it does in the browser
    store = os.path.join(os.environ.get("LOCALAPPDATA") or os.path.expanduser("~"), "refit-viewer")
    webview.start(menu=menu, private_mode=False, storage_path=store)


class Server(ThreadingHTTPServer):
    # ⚠️ on Windows SO_REUSEADDR lets a SECOND server bind a port already in
    # use, and requests then land on either one at random: a viewer left open
    # answered for the new one. Without it the bind fails and we pick another port.
    allow_reuse_address = os.name != "nt"
    daemon_threads = True


def main():
    ap = argparse.ArgumentParser(description="Local viewer for refit.py")
    ap.add_argument("file", nargs="?", help=".stl / .step to open right away")
    ap.add_argument("--port", type=int, default=8765)
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--browser", action="store_true", help="open the page in the default browser instead of a window")
    g.add_argument("--serve", "--no-browser", action="store_true", dest="serve", help="only run the server")
    a = ap.parse_args()
    try:
        srv = Server(("127.0.0.1", a.port), Handler)
    except OSError:
        # port taken (another viewer open): any free port will do
        srv = Server(("127.0.0.1", 0), Handler)
    url = f"http://127.0.0.1:{srv.server_address[1]}/"
    if a.file:
        from urllib.parse import quote

        url += "?open=" + quote(os.path.abspath(a.file))
    # OpenCascade loads in the background while the page opens
    threading.Thread(target=refit_module, daemon=True).start()

    window = not (a.browser or a.serve)
    if window:
        try:
            import webview  # noqa: F401
        except ImportError:
            print("pywebview is not installed (pip install pywebview): opening the browser instead", flush=True)
            window = False
            a.browser = True
    try:
        if window:
            threading.Thread(target=srv.serve_forever, daemon=True).start()
            run_window(url)
        else:
            print(f"refit viewer on {url}  (Ctrl+C to stop)", flush=True)
            if a.browser:
                threading.Timer(0.6, lambda: webbrowser.open(url)).start()
            srv.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        for j in JOBS.values():
            j.cancel()
        srv.server_close()


if __name__ == "__main__":
    main()
