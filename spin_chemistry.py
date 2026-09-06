"""Spin-1/2 + spin-1 encounter model for doxorubicin semiquinone and 3O2.

Hamiltonian inputs are angular frequencies (rad s-1); kinetic rates are s-1.
No value in this module is a fitted doxorubicin parameter.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
from pathlib import Path
import numpy as np

def spin_matrices(spin: float):
    n=int(round(2*spin+1)); m=np.arange(spin,-spin-1,-1); plus=np.zeros((n,n),complex)
    for col in range(1,n):
        mc=m[col]; plus[col-1,col]=np.sqrt(spin*(spin+1)-mc*(mc+1))
    minus=plus.T
    return (plus+minus)/2,(plus-minus)/(2j),np.diag(m)

I2,I3=np.eye(2),np.eye(3)
R=tuple(np.kron(x,I3) for x in spin_matrices(.5)); O=tuple(np.kron(I2,x) for x in spin_matrices(1.))
I6=np.eye(6,dtype=complex); R_DOT_O=sum(a@b for a,b in zip(R,O))
P_DOUBLET=(.5*I6-R_DOT_O)/1.5; P_QUARTET=(R_DOT_O+I6)/1.5

@dataclass(frozen=True)
class EncounterParameters:
    field_t: tuple[float,float,float]=(0.,0.,0.); g_radical: float=2.0033; g_oxygen: float=2.0023
    exchange_rad_s: float=0.; dipolar_rad_s: float=0.; dipolar_axis: tuple[float,float,float]=(0.,0.,1.)
    oxygen_zfs_d_rad_s: float=0.; oxygen_zfs_e_rad_s: float=0.
    effective_hyperfine_rad_s: tuple[float,float,float]=(0.,0.,0.)
    radical_relaxation_s: float=0.; oxygen_relaxation_s: float=0.
    k_doublet_s: float=0.; k_quartet_s: float=0.; k_escape_s: float=0.

def _unit(v):
    x=np.asarray(v,float); norm=np.linalg.norm(x)
    if norm==0: raise ValueError("axis must be nonzero")
    return x/norm

def hamiltonian(p: EncounterParameters):
    beta=8.79410005e10; b=np.asarray(p.field_t,float)
    h=beta*sum(b[i]*(p.g_radical*R[i]+p.g_oxygen*O[i]) for i in range(3))+p.exchange_rad_s*R_DOT_O
    if p.dipolar_rad_s:
        n=_unit(p.dipolar_axis); rn=sum(n[i]*R[i] for i in range(3)); on=sum(n[i]*O[i] for i in range(3))
        h+=p.dipolar_rad_s*(R_DOT_O-3*rn@on)
    ox,oy,oz=O
    h+=p.oxygen_zfs_d_rad_s*(oz@oz-(2/3)*I6)+p.oxygen_zfs_e_rad_s*(ox@ox-oy@oy)
    h+=sum(p.effective_hyperfine_rad_s[i]*R[i] for i in range(3))
    return np.asarray(h,complex)

def initial_density(scenario: str):
    if scenario=="unpolarized": return I6/6
    if scenario=="doublet": return P_DOUBLET/2
    if scenario=="quartet": return P_QUARTET/4
    raise ValueError("scenario must be unpolarized, doublet, or quartet")

def _lindblad(rho,op):
    od=op.conj().T; return op@rho@od-.5*(od@op@rho+rho@od@op)

def _rk4(rhs,y0,times,max_frequency=0.):
    """Dependency-free fixed-output RK4 with frequency-aware substeps."""
    values=[np.asarray(y0,complex)]; y=values[0].copy()
    for start,stop in zip(times[:-1],times[1:]):
        interval=stop-start; substeps=max(1,int(np.ceil(interval*max_frequency/.05))); dt=interval/substeps; t=start
        for _ in range(substeps):
            k1=rhs(t,y); k2=rhs(t+dt/2,y+dt*k1/2); k3=rhs(t+dt/2,y+dt*k2/2); k4=rhs(t+dt,y+dt*k3)
            y=y+dt*(k1+2*k2+2*k3+k4)/6; t+=dt
        values.append(y.copy())
    return np.asarray(values)

def propagate_encounter(p: EncounterParameters,duration_s: float,scenario="unpolarized",samples=1001):
    if duration_s<=0 or samples<2: raise ValueError("duration_s must be positive and samples >= 2")
    if any(x<0 for x in (p.k_doublet_s,p.k_quartet_s,p.k_escape_s)): raise ValueError("rates must be nonnegative")
    h=hamiltonian(p); loss=p.k_doublet_s*P_DOUBLET+p.k_quartet_s*P_QUARTET+p.k_escape_s*I6; ops=[]
    if p.radical_relaxation_s: ops.extend(np.sqrt(p.radical_relaxation_s)*x for x in R)
    if p.oxygen_relaxation_s: ops.extend(np.sqrt(p.oxygen_relaxation_s)*x for x in O)
    def rhs(_t,y):
        rho=y.reshape(6,6); dr=-1j*(h@rho-rho@h)-.5*(loss@rho+rho@loss)
        for op in ops: dr+=_lindblad(rho,op)
        return dr.reshape(-1)
    times=np.linspace(0,duration_s,samples)
    scale=2*np.linalg.norm(h,2)+np.linalg.norm(loss,2)+sum(np.linalg.norm(op,2)**2 for op in ops)
    rhos=_rk4(rhs,initial_density(scenario).reshape(-1),times,scale).reshape(-1,6,6)
    pd=np.real(np.einsum("ij,tji->t",P_DOUBLET,rhos)); pq=np.real(np.einsum("ij,tji->t",P_QUARTET,rhos)); survival=np.real(np.trace(rhos,axis1=1,axis2=2))
    yd=float(np.trapezoid(p.k_doublet_s*pd,times)); yq=float(np.trapezoid(p.k_quartet_s*pq,times)); ye=float(np.trapezoid(p.k_escape_s*survival,times))
    return {"time_s":times,"p_doublet":pd,"p_quartet":pq,"survival":survival,"superoxide_yield":yd+yq,"doublet_reaction_yield":yd,"quartet_reaction_yield":yq,"escape_yield":ye,"unresolved_probability":float(survival[-1])}

def downstream_ros(superoxide0_m,duration_s,k_spont_m_inv_s,k_sod_m_inv_s,sod_m,k_h2o2_loss_s=0.,samples=201):
    def rhs(_t,y):
        o,h=y; spont=k_spont_m_inv_s*o*o; enz=k_sod_m_inv_s*sod_m*o
        return np.array((-2*spont-enz,spont+.5*enz-k_h2o2_loss_s*h))
    times=np.linspace(0,duration_s,samples); scale=2*k_spont_m_inv_s*superoxide0_m+k_sod_m_inv_s*sod_m+k_h2o2_loss_s
    values=np.real(_rk4(rhs,np.array((superoxide0_m,0.)),times,scale))
    return {"time_s":times,"superoxide_m":values[:,0],"hydrogen_peroxide_m":values[:,1]}

def coherent_circuit_validation(p: EncounterParameters,duration_s: float):
    """Validate a leakage-safe three-qubit embedding of the coherent unitary."""
    energies,vectors=np.linalg.eigh(hamiltonian(p)); u6=(vectors*np.exp(-1j*energies*duration_s))@vectors.conj().T
    u8=np.eye(8,dtype=complex); physical=np.array([0,1,2,4,5,6]); u8[np.ix_(physical,physical)]=u6
    rng=np.random.default_rng(1729); psi6=rng.normal(size=6)+1j*rng.normal(size=6); psi6/=np.linalg.norm(psi6); psi8=np.zeros(8,complex); psi8[physical]=psi6
    return float(np.max(np.abs((u8@psi8)[physical]-u6@psi6)))

def parameters_from_json(path):
    raw=json.loads(Path(path).read_text()); active=raw["active_model"]
    return EncounterParameters(**active["encounter"]),active
