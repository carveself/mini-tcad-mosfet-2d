"""Conservative cell-centred finite volumes for a depleted double-gate MOSFET.

Coordinates are in nm; potentials in V; terminal charges are C/m of width.
Only fixed ionized acceptors are included. No transport/current is inferred.
"""
from dataclasses import dataclass
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from scipy.optimize import brentq

EPS0 = 8.8541878128e-12
Q = 1.602176634e-19
TERMINALS = ("source", "drain", "bottom_gate", "top_gate")


@dataclass(frozen=True)
class Device:
    length_nm: float = 30.0
    silicon_nm: float = 10.0
    oxide_nm: float = 2.0
    oxide_k: float = 3.9
    silicon_k: float = 11.7
    acceptors_cm3: float = 1e17
    contact_V: float = 0.9
    flatband_V: float = 0.0

    def __post_init__(self):
        values = list(self.__dict__.values())
        if not np.all(np.isfinite(values)):
            raise ValueError("All device parameters must be finite")
        if min(self.length_nm, self.silicon_nm, self.oxide_nm,
               self.oxide_k, self.silicon_k) <= 0 or self.acceptors_cm3 < 0:
            raise ValueError("Lengths/permittivities must be positive; doping nonnegative")

    @property
    def scale_length_nm(self):
        """Thin-oxide parabolic estimate, NOT an exact eigenmode solution."""
        return np.sqrt(self.silicon_k * self.silicon_nm * self.oxide_nm /
                       (2 * self.oxide_k))

    @property
    def eigen_length_nm(self):
        """Fundamental symmetric layered Laplace mode, including oxide geometry."""
        a, t = self.silicon_nm/2, self.oxide_nm
        f = lambda k: np.tan(k*a)*np.tan(k*t)-self.oxide_k/self.silicon_k
        k = brentq(f,1e-12,(np.pi/(2*max(a,t)))*(1-1e-9),xtol=1e-14)
        return 1/k


def edges(device, h_nm):
    if not np.isfinite(h_nm) or h_nm <= 0:
        raise ValueError("Mesh spacing must be positive and finite")
    nx = max(8, int(np.ceil(device.length_nm / h_nm)))
    ns = max(4, 2 * int(np.ceil(device.silicon_nm / h_nm / 2)))
    no = max(2, int(np.ceil(device.oxide_nm / h_nm)))
    if nx * (ns + 2 * no) > 500_000:
        raise ValueError("Requested grid exceeds the teaching solver's cell limit")
    t = device.silicon_nm / 2
    z = np.concatenate((np.linspace(-t-device.oxide_nm, -t, no+1)[:-1],
                        np.linspace(-t, t, ns+1)[:-1],
                        np.linspace(t, t+device.oxide_nm, no+1)))
    return np.linspace(0, device.length_nm, nx+1), z


def vertical_exact(device, z_nm, gate_V):
    """Long-channel symmetric layered Poisson solution, body charge -q NA."""
    z = np.asarray(z_nm)
    a = device.silicon_nm / 2
    r = Q * device.acceptors_cm3 * 1e6 / EPS0 * 1e-18
    surface = gate_V-device.flatband_V-r*a*device.oxide_nm/device.oxide_k
    return np.where(np.abs(z) <= a,
                    surface + r*(z*z-a*a)/(2*device.silicon_k),
                    gate_V-device.flatband_V-r*a*(a+device.oxide_nm-np.abs(z))/device.oxide_k)


class FiniteVolume:
    def __init__(self, device=Device(), h_nm=0.5):
        self.device = device
        self.xe, self.ze = edges(device, h_nm)
        self.x = (self.xe[:-1]+self.xe[1:])/2
        self.z = (self.ze[:-1]+self.ze[1:])/2
        self.dx, self.dz = np.diff(self.xe), np.diff(self.ze)
        self.nx, self.nz = len(self.x), len(self.z)
        self.si = np.abs(self.z) < device.silicon_nm/2
        self.k = np.where(self.si, device.silicon_k, device.oxide_k)
        n = self.nx*self.nz
        rows, cols, data = [], [], []
        diag = np.zeros(n)
        self.boundaries = []
        self.B = np.zeros((n, 4))

        def edge(p, q, conductance):
            diag[p] += conductance; diag[q] += conductance
            rows.extend((p,q)); cols.extend((q,p)); data.extend((-conductance,-conductance))

        def boundary(p, g, terminal, x, z):
            diag[p] += g; self.B[p, terminal] += g
            self.boundaries.append((p,g,terminal,x,z))

        for j in range(self.nz):
            for i in range(self.nx):
                p = j*self.nx+i
                if i+1 < self.nx:
                    edge(p,p+1,self.k[j]*self.dz[j]/((self.dx[i]+self.dx[i+1])/2))
                if j+1 < self.nz:
                    # Series dielectric resistance: weighted harmonic face coefficient.
                    g = self.dx[i]/(self.dz[j]/(2*self.k[j])+self.dz[j+1]/(2*self.k[j+1]))
                    edge(p,p+self.nx,g)
                if self.si[j] and i == 0:
                    boundary(p,self.k[j]*self.dz[j]/(self.dx[i]/2),0,0,self.z[j])
                if self.si[j] and i == self.nx-1:
                    boundary(p,self.k[j]*self.dz[j]/(self.dx[i]/2),1,device.length_nm,self.z[j])
                if j == 0:
                    boundary(p,self.k[j]*self.dx[i]/(self.dz[j]/2),2,self.x[i],self.ze[0])
                if j == self.nz-1:
                    boundary(p,self.k[j]*self.dx[i]/(self.dz[j]/2),3,self.x[i],self.ze[-1])
        rows.extend(range(n)); cols.extend(range(n)); data.extend(diag)
        self.A = coo_matrix((data,(rows,cols)),shape=(n,n)).tocsc()
        self.lu = splu(self.A)
        self.volume_nm2 = self.dz[:,None]*self.dx[None,:]
        self.rho = np.broadcast_to(np.where(self.si[:,None],
                       -Q*device.acceptors_cm3*1e6,0), (self.nz,self.nx)).copy()
        self.rhs_charge = self.source_rhs(self.rho)
        self.offset = self.lu.solve(self.rhs_charge)
        self.response = self.lu.solve(self.B)

    def source_rhs(self, rho):
        return (np.broadcast_to(rho,(self.nz,self.nx))*self.volume_nm2/EPS0*1e-18).ravel()

    def voltages(self, vg, vd, vs=0.0, bottom_vg=None):
        d = self.device
        v = np.array([d.contact_V+vs,d.contact_V+vd,
                      (vg if bottom_vg is None else bottom_vg)-d.flatband_V,
                      vg-d.flatband_V])
        if not np.all(np.isfinite(v)):
            raise ValueError("Voltages must be finite")
        return v

    def solve(self, vg=0.5, vd=0.05, vs=0.0, bottom_vg=None):
        return (self.offset+self.response@self.voltages(vg,vd,vs,bottom_vg)).reshape(self.nz,self.nx)

    def solve_custom(self, rho, boundary_value):
        rhs = self.source_rhs(rho)
        for p,g,t,x,z in self.boundaries:
            rhs[p] += g*boundary_value(t,x,z)
        return self.lu.solve(rhs).reshape(self.nz,self.nx)

    def centerline(self, potential, voltage):
        # Even number of silicon rows: interpolate symmetrically to z=0.
        middle = np.flatnonzero(self.si)
        j = middle[len(middle)//2]
        line = (potential[j-1]+potential[j])/2
        return np.r_[0,self.x,self.device.length_nm], np.r_[voltage[0],line,voltage[1]]

    def barrier(self, vg, vd, criterion_eV=0.2):
        if not np.isfinite(criterion_eV) or criterion_eV <= 0:
            raise ValueError("Barrier criterion must be positive")
        v = self.voltages(vg,vd)
        x, psi = self.centerline(self.solve(vg,vd),v)
        i = int(np.argmin(psi))
        xmin, value = x[i], psi[i]
        if 0 < i < len(x)-1:
            # Local minimum improves barrier extraction without a grid-node jump.
            c = np.polyfit(x[i-1:i+2]-x[i],psi[i-1:i+2],2)
            shift = -c[1]/(2*c[0]) if c[0] > 0 else np.inf
            if x[i-1]-x[i] <= shift <= x[i+1]-x[i]:
                xmin, value = x[i]+shift, np.polyval(c,shift)
        return {"barrier_eV": float(v[0]-value), "x_barrier_nm": float(xmin),
                "minimum_potential_V":float(value)}

    def threshold(self, vd, criterion_eV=0.2, bracket=(-1.0,1.3)):
        f = lambda vg: self.barrier(vg,vd,criterion_eV)["barrier_eV"]-criterion_eV
        if f(bracket[0])*f(bracket[1]) > 0:
            raise ValueError("Barrier threshold is not bracketed; change the gate interval")
        return float(brentq(f,*bracket,xtol=1e-11))

    def dibl(self, low=0.05, high=0.7, criterion_eV=0.2):
        if not high > low:
            raise ValueError("High drain voltage must exceed low drain voltage")
        a, b = self.threshold(low,criterion_eV),self.threshold(high,criterion_eV)
        return {"threshold_low_V":a,"threshold_high_V":b,
                "dibl_mV_V":float(1000*(a-b)/(high-low))}

    def terminal_charges(self, potential, voltages):
        q = np.zeros(4)
        for p,g,t,_,_ in self.boundaries:
            q[t] += EPS0*g*(voltages[t]-potential.ravel()[p])
        return q

    def capacitance(self):
        c = np.zeros((4,4))
        for p,g,t,_,_ in self.boundaries:
            c[t] += EPS0*g*(np.eye(4)[t]-self.response[p])
        return c

    def residual(self, potential, voltages):
        rhs = self.rhs_charge+self.B@voltages
        return float(np.max(np.abs(self.A@potential.ravel()-rhs))/max(np.max(np.abs(rhs)),1e-30))


def manufactured(device, x, z):
    """Exact layered solution with continuous normal dielectric displacement."""
    a, tox = device.silicon_nm/2,device.oxide_nm
    k = np.where(np.abs(z) <= a,device.silicon_k,device.oxide_k)
    resistance = np.where(z < -a,(z+a+tox)/device.oxide_k,
                  np.where(z <= a,tox/device.oxide_k+(z+a)/device.silicon_k,
                  tox/device.oxide_k+2*a/device.silicon_k+(z-a)/device.oxide_k))
    total = 2*tox/device.oxide_k+2*a/device.silicon_k
    t = resistance/total
    wave = 2*np.pi/device.length_nm
    psi = np.cos(wave*x)*t*(1-t)
    rho_normalized = np.cos(wave*x)*(k*wave**2*t*(1-t)+2/(k*total**2))
    return psi, rho_normalized*EPS0/1e-18
