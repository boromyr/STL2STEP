import sys, math
sys.path.insert(0, "D:/Users/PC/Desktop/refit")
import numpy as np
from refit import *
import refit as _rf
_st=_rf._st
Log.level=30

def tri_dist(P, A, B, C):
    """vectorized point-to-triangle distance (P: n x 3, A/B/C: m x 3) -> n x m"""
    out=np.empty((len(P),len(A)))
    AB=B-A; AC=C-A; N=np.cross(AB,AC)
    nn=np.linalg.norm(N,axis=1); nn[nn<1e-20]=1e-20; Nu=N/nn[:,None]
    for i,p in enumerate(P):
        w=p-A
        d_plane=np.einsum('ij,ij->i',w,Nu)
        proj=p-d_plane[:,None]*Nu
        # barycentric coordinates
        v0=AB; v1=AC; v2=proj-A
        d00=np.einsum('ij,ij->i',v0,v0); d01=np.einsum('ij,ij->i',v0,v1)
        d11=np.einsum('ij,ij->i',v1,v1); d20=np.einsum('ij,ij->i',v2,v0); d21=np.einsum('ij,ij->i',v2,v1)
        den=d00*d11-d01*d01; den[np.abs(den)<1e-20]=1e-20
        v=(d11*d20-d01*d21)/den; w2=(d00*d21-d01*d20)/den; u=1-v-w2
        inside=(u>=-1e-9)&(v>=-1e-9)&(w2>=-1e-9)
        d=np.abs(d_plane)
        # outside: distance to the three segments
        if not inside.all():
            idx=~inside
            best=np.full(idx.sum(), np.inf)
            for X,Y in ((A,B),(B,C),(C,A)):
                e=(Y-X)[idx]; f=p-X[idx]
                t=np.clip(np.einsum('ij,ij->i',f,e)/np.maximum(np.einsum('ij,ij->i',e,e),1e-20),0,1)
                best=np.minimum(best, np.linalg.norm(f-t[:,None]*e,axis=1))
            d=d.copy(); d[idx]=best
        out[i]=d
    return out

def main():
    mesh,_=read_stl(sys.argv[1])
    tm=Topo(mesh)
    A=[];B=[];C=[]
    for i in range(tm.nF):
        V=tm.verts[i]
        for k in range(1,len(V)-1):
            A.append(V[0]);B.append(V[k]);C.append(V[k+1])
    A=np.array(A);B=np.array(B);C=np.array(C)
    out,_=read_step(sys.argv[2])
    to=Topo(out)
    # sample the vertices + centers of the output faces, and interior points via triangulation
    from OCP.BRepMesh import BRepMesh_IncrementalMesh
    BRepMesh_IncrementalMesh(out, 0.05, False, 0.5, True)
    pts=[]
    for i,f in enumerate(to.faces):
        loc=TopLoc_Location()
        tri=_st(BRep_Tool,"Triangulation")(f, loc)
        if tri is None: continue
        tr=loc.Transformation()
        n=tri.NbNodes()
        for k in range(1, n+1):
            q=tri.Node(k).Transformed(tr)
            pts.append([q.X(),q.Y(),q.Z()])
    P=np.array(pts)
    if len(sys.argv)>3:
        n=int(sys.argv[3])
        if len(P)>n:
            P=P[np.random.default_rng(0).choice(len(P),n,replace=False)]
    print("points sampled on the output:", len(P), " mesh triangles:", len(A))
    d=np.empty(len(P))
    step=200
    for s in range(0,len(P),step):
        d[s:s+step]=tri_dist(P[s:s+step],A,B,C).min(axis=1)
    print(f"output->mesh distance:  mean {d.mean():.5f}  p99 {np.percentile(d,99):.5f}  max {d.max():.5f} mm")
    print(f"points beyond 0.05 mm: {(d>0.05).sum()}   beyond 0.1 mm: {(d>0.1).sum()}")

if __name__=="__main__":
    main()
