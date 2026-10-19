"""Independent P1 triangular finite elements; no FVM stencil/face flux reuse."""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from scipy.interpolate import RegularGridInterpolator
from electrostatics import Device, EPS0, Q, edges


class FiniteElement:
    def __init__(self, device=Device(), h_nm=0.5):
        self.device = device
        self.x,self.z = edges(device,h_nm)
        xx,zz = np.meshgrid(self.x,self.z)
        points = np.column_stack((xx.ravel(),zz.ravel()))
        nx,nz = len(self.x),len(self.z)
        rows,cols,data = [],[],[]
        rhs = np.zeros(nx*nz)
        for j in range(nz-1):
            for i in range(nx-1):
                a=j*nx+i; b=a+1; c=a+nx; d=c+1
                triangles = ((a,b,d),(a,d,c)) if (i+j)%2 == 0 else ((a,b,c),(b,d,c))
                for ids in triangles:
                    ids=np.array(ids); xy=points[ids]
                    area=abs(np.linalg.det(np.column_stack((xy[1]-xy[0],xy[2]-xy[0]))))/2
                    coefficients=np.linalg.inv(np.column_stack((np.ones(3),xy)))
                    gradients=coefficients[1:,:]
                    silicon=abs(np.mean(xy[:,1])) < device.silicon_nm/2
                    k=device.silicon_k if silicon else device.oxide_k
                    stiffness=k*area*(gradients.T@gradients)
                    rows.extend(np.repeat(ids,3)); cols.extend(np.tile(ids,3)); data.extend(stiffness.ravel())
                    rho=-Q*device.acceptors_cm3*1e6 if silicon else 0
                    rhs[ids] += rho/EPS0*1e-18*area/3
        self.A=coo_matrix((data,(rows,cols)),shape=(nx*nz,nx*nz)).tocsc()
        self.rhs=rhs
        self.points=points
        # The interface endpoints belong to the ideal source/drain contacts.
        si=np.abs(points[:,1]) <= device.silicon_nm/2+1e-10
        self.groups=[np.flatnonzero((points[:,0]==0)&si),
                     np.flatnonzero((points[:,0]==device.length_nm)&si),
                     np.flatnonzero(points[:,1]==self.z[0]),
                     np.flatnonzero(points[:,1]==self.z[-1])]
        self.fixed=np.concatenate(self.groups)
        self.free=np.setdiff1d(np.arange(nx*nz),self.fixed)

    def solve(self, vg=0.5,vd=0.05):
        d=self.device
        voltages=[d.contact_V,d.contact_V+vd,vg-d.flatband_V,vg-d.flatband_V]
        u=np.zeros(len(self.points))
        for ids,v in zip(self.groups,voltages): u[ids]=v
        rhs=self.rhs[self.free]-self.A[self.free][:,self.fixed]@u[self.fixed]
        u[self.free]=spsolve(self.A[self.free][:,self.free],rhs)
        return u.reshape(len(self.z),len(self.x))

    def interpolate(self, u, x, z):
        xx,zz=np.meshgrid(x,z)
        return RegularGridInterpolator((self.z,self.x),u)(np.column_stack((zz.ravel(),xx.ravel()))).reshape(xx.shape)

    def barrier(self, u):
        line=RegularGridInterpolator((self.z,self.x),u)(np.column_stack((np.zeros(len(self.x)),self.x)))
        i=int(np.argmin(line)); value=line[i]; xmin=self.x[i]
        if 0<i<len(line)-1:
            c=np.polyfit(self.x[i-1:i+2]-self.x[i],line[i-1:i+2],2)
            shift=-c[1]/(2*c[0]) if c[0]>0 else np.inf
            if self.x[i-1]-self.x[i] <= shift <= self.x[i+1]-self.x[i]:
                value=np.polyval(c,shift); xmin+=shift
        return {"barrier_eV":float(self.device.contact_V-value),"x_barrier_nm":float(xmin)}
