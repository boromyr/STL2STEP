# -*- coding: utf-8 -*-
"""Renders a B-Rep the way a CAD viewer would: shaded faces + black edges."""
import sys, numpy as np
sys.path.insert(0, "D:/Users/PC/Desktop/refit")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
import refit as R
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.BRep import BRep_Tool
from OCP.TopLoc import TopLoc_Location
from OCP.TopAbs import TopAbs_FACE, TopAbs_EDGE
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.GCPnts import GCPnts_QuasiUniformAbscissa

def tri_and_edges(shape, defl):
    BRepMesh_IncrementalMesh(shape, defl, False, 0.2, True)
    tris = []
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        f = R.td_Face(ex.Current()); ex.Next()
        loc = TopLoc_Location()
        t = R.bt_Triangulation(f, loc)
        if t is None: continue
        tr = loc.Transformation()
        pts = []
        for i in range(1, t.NbNodes() + 1):
            p = t.Node(i).Transformed(tr)
            pts.append((p.X(), p.Y(), p.Z()))
        pts = np.array(pts)
        for i in range(1, t.NbTriangles() + 1):
            a, b, c = t.Triangle(i).Get()
            tris.append(pts[[a-1, b-1, c-1]])
    segs = []
    ex = TopExp_Explorer(shape, TopAbs_EDGE)
    seen = set()
    while ex.More():
        e = R.td_Edge(ex.Current()); ex.Next()
        h = e.TShape().This()
        if h in seen: continue
        seen.add(h)
        try:
            ad = BRepAdaptor_Curve(e)
            n = 24
            ab = GCPnts_QuasiUniformAbscissa(ad, n)
            P = [ad.Value(ab.Parameter(i)) for i in range(1, ab.NbPoints()+1)]
            segs.append(np.array([[p.X(), p.Y(), p.Z()] for p in P]))
        except Exception:
            pass
    return np.array(tris), segs

def cam(eye, target, up=np.array([0.,0.,1.])):
    f = target - eye; f /= np.linalg.norm(f)
    r = np.cross(f, up); r /= np.linalg.norm(r)
    u = np.cross(r, f)
    return np.stack([r, u, -f])

def main(path, out, eye_dir, target=None, zoom=None):
    sh, _ = R.read_input(path)
    topo = R.Topo(sh)
    lo = topo.vpos.min(axis=0); hi = topo.vpos.max(axis=0)
    ctr = target if target is not None else (lo + hi) / 2
    diag = float(np.linalg.norm(hi - lo))
    eye = ctr + eye_dir / np.linalg.norm(eye_dir) * diag * 2
    M = cam(eye, ctr)
    tris, segs = tri_and_edges(sh, diag * 2e-4)
    T = (tris - eye) @ M.T
    nr = np.cross(T[:,1]-T[:,0], T[:,2]-T[:,0]); L = np.linalg.norm(nr, axis=1)
    nr /= np.maximum(L, 1e-12)[:, None]
    light = np.array([0.35, 0.45, 0.82]); light /= np.linalg.norm(light)
    shade = 0.30 + 0.62 * np.clip(np.abs(nr @ light), 0, 1)
    z = T[:, :, 2].mean(axis=1)
    order = np.argsort(z)
    fig, ax = plt.subplots(figsize=(13, 10))
    ax.add_collection(PolyCollection(T[order][:, :, :2],
                                     facecolors=np.stack([shade]*3, axis=1)[order],
                                     edgecolors="none"))
    X = T[:, :, 0]; Y = T[:, :, 1]
    x0, x1 = (zoom[0], zoom[1]) if zoom else (X.min(), X.max())
    y0, y1 = (zoom[2], zoom[3]) if zoom else (Y.min(), Y.max())
    # --- z-buffer to remove hidden edges ---------------------------------------
    W = 1400; H = max(8, int(W * (y1 - y0) / max(x1 - x0, 1e-9)))
    zb = np.full((H, W), -1e30)
    sx = (W - 1) / (x1 - x0); sy = (H - 1) / (y1 - y0)
    def px(P):
        return (P[..., 0] - x0) * sx, (P[..., 1] - y0) * sy
    tu, tv = px(T)
    tz = T[..., 2]
    for k in range(len(T)):
        u = tu[k]; v = tv[k]
        iu0 = max(0, int(np.floor(u.min()))); iu1 = min(W - 1, int(np.ceil(u.max())))
        iv0 = max(0, int(np.floor(v.min()))); iv1 = min(H - 1, int(np.ceil(v.max())))
        if iu1 < iu0 or iv1 < iv0: continue
        gu, gv = np.meshgrid(np.arange(iu0, iu1 + 1), np.arange(iv0, iv1 + 1))
        d = ((v[1]-v[2])*(u[0]-u[2]) + (u[2]-u[1])*(v[0]-v[2]))
        if abs(d) < 1e-12: continue
        l0 = ((v[1]-v[2])*(gu-u[2]) + (u[2]-u[1])*(gv-v[2])) / d
        l1 = ((v[2]-v[0])*(gu-u[2]) + (u[0]-u[2])*(gv-v[2])) / d
        l2 = 1 - l0 - l1
        m = (l0 >= -0.002) & (l1 >= -0.002) & (l2 >= -0.002)
        if not m.any(): continue
        zz = l0*tz[k,0] + l1*tz[k,1] + l2*tz[k,2]
        sub = zb[iv0:iv1+1, iu0:iu1+1]
        np.maximum(sub, np.where(m, zz, -1e30), out=sub)
    eps = (T[..., 2].max() - T[..., 2].min()) * 2e-3 + 1e-9
    S = []
    for s in segs:
        Q = (s - eye) @ M.T
        u, v = px(Q)
        iu = np.clip(np.round(u).astype(int), 0, W - 1)
        iv = np.clip(np.round(v).astype(int), 0, H - 1)
        vis = Q[:, 2] > zb[iv, iu] - eps
        run = []
        for i in range(len(Q)):
            if vis[i]:
                run.append(Q[i, :2])
            else:
                if len(run) > 1: S.append(np.array(run))
                run = []
        if len(run) > 1: S.append(np.array(run))
    ax.add_collection(LineCollection(S, colors="k", linewidths=0.6))
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_facecolor("#202020"); fig.patch.set_facecolor("#202020")
    plt.tight_layout(); plt.savefig(out, dpi=100)
    print("saved", out, len(tris), "triangles", len(segs), "edges")

if __name__ == "__main__":
    tgt = None; zm = None
    if len(sys.argv) > 3:
        cx, cy, cz, hw = [float(x) for x in sys.argv[3].split(",")]
        tgt = np.array([cx, cy, cz])
        zm = (-hw, hw, -hw * 0.77, hw * 0.77)
    d = np.array([-0.55, -0.6, 0.58])
    if len(sys.argv) > 4:
        d = np.array([float(x) for x in sys.argv[4].split(",")])
    main(sys.argv[1], sys.argv[2], d, target=tgt, zoom=zm)
