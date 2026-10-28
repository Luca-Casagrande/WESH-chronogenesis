#!/usr/bin/env python3
"""Checks for the bounded Section-2 revision, using explicit finite matrices.

No continuum-gravity closure, stochastic cosmology or BBN validation is claimed.
The smooth-generator argument is analytic; its proof is in the manuscript, Section 2.4.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

TOL=3e-11
RESULTS={}

def adj(a):
    return a.conj().T

def comm(a,b):
    return a@b-b@a

def check(name,residual):
    val=float(np.linalg.norm(residual))
    RESULTS[name]=val
    if not np.isfinite(val) or val>TOL:
        raise AssertionError(f'{name}: {val}')

def nonzero(name,a):
    v=float(np.linalg.norm(a))
    RESULTS[name]=v
    if not np.isfinite(v) or v<1e-5:
        raise AssertionError(f'Unexpected zero in {name}: {v}')

def constraint_checks():
    sx=np.array([[0,1],[1,0]],complex)
    sz=np.diag([1.,-1.]).astype(complex)
    i=np.eye(2)
    c1=np.kron(sx,np.eye(4))+np.kron(i,np.kron(sz,i))
    c2=np.kron(sz,np.eye(4))+np.kron(i,np.kron(i,sz))
    check('additive_nonabelian_commutator_identity',comm(c1,c2)-np.kron(comm(sx,sz),np.eye(4)))
    nonzero('additive_nonabelian_defect_norm',comm(c1,c2))
    # Independent local constraint blocks, exactly the class used in the text.
    local=np.kron(-sz,i)+np.kron(i,sz)
    lc1=np.kron(local,np.eye(4)); lc2=np.kron(np.eye(4),local)
    check('commuting_independent_block_constraints',comm(lc1,lc2))

def exchange_checks(d):
    qdim=2
    p=np.arange(-(d//2),d//2+1,dtype=float)
    nodes=p/(d//2)
    f=np.exp(2j*np.pi*np.outer(np.arange(d),np.arange(d))/d)/np.sqrt(d)
    t=adj(f)@np.diag(nodes)@f
    j=t@t
    iq=np.eye(qdim)
    jr=np.kron(j,iq)
    pr=np.kron(np.diag(p),iq)
    hr=-pr
    v=np.zeros((d*qdim*d,d*qdim),complex)
    for n in range(d):
        for q in range(qdim):
            v[(n*qdim+q)*d+n,n*qdim+q]=1.
    lift=lambda a:v@a@adj(v)
    jl,hl,pl=map(lift,(jr,hr,pr))
    check(f'd{d}_sector_exchange_identity',comm(jl,hl)-lift(np.kron(comm(np.diag(p),j),iq)))
    check(f'd{d}_opposite_sector_commutators',comm(jl,hl)+comm(jl,pl))
    nonzero(f'd{d}_gm_noncommutation_norm',comm(jl,hl))
    charges=lift(np.kron(np.eye(d),np.diag([-1.,1.])))
    check(f'd{d}_spectator_charge_conservation',comm(jl,charges))
    # A diagonal block at a fixed charge is the complete nondegenerate commutant.
    a=np.arange(d*d).reshape(d,d).astype(complex)
    elem=(p[None,:]-p[:,None])*a
    check(f'd{d}_nondegenerate_commutant_components',comm(a,np.diag(p))-elem)
    # Demonstrate that noncommutation is not a nonzero mean drift in every state.
    rho=np.eye(d*qdim)/(d*qdim)
    ladj_h=jr@hr@jr-(jr@jr@hr+hr@jr@jr)/2
    nonzero(f'd{d}_gm_adjoint_action_norm',ladj_h)
    check(f'd{d}_maximally_mixed_mean_drift_zero',np.trace(rho@ladj_h))

def profile_checks():
    edge0=np.array([0.4,2.,3.6])
    mean=edge0.mean()
    kernel=edge0/mean
    old_gamma=0.7; nu=1.3
    new_gamma=old_gamma*mean
    check('profile_mean_one',kernel.mean()-1.)
    check('profile_rescaling_preserves_edge_coefficients',old_gamma*edge0-new_gamma*kernel)
    ratio=new_gamma/nu
    check('edge_specific_block_ratio',new_gamma*kernel/nu-ratio*kernel)
    check('single_edge_ratio',new_gamma/nu-ratio)
    RESULTS['ratio_not_fixed_to_unity']=float(ratio)

def polynomial_checks():
    # Exact symbolic checks of the counterexamples; not substitutes for
    # the general smooth-generator proof in the text.
    import sympy as sp
    u,r=sp.symbols('u r',positive=True)
    check('homogeneous_nongenerator_example',float(sp.simplify((r*u)**2-r**2*u**2)))
    if sp.diff(u**2,u).subs(u,0)!=0: raise AssertionError('Derivative condition')
    defect=sp.expand((r*u+(r*u)**2)-r*(u+u**2))
    if sp.simplify(defect-r*(r-1)*u**2)!=0: raise AssertionError('Homogeneity defect')
    RESULTS['u_plus_u2_homogeneity_defect']=str(defect)
    # Nonlinear primitive f(u)=u+u^2 has gain 1+u+v in coordinate u.
    v=sp.symbols('v')
    if sp.cancel(((u+u*u)-(v+v*v))/(u-v)) != u+v+1:
        raise AssertionError('Coordinate-dependent gain')
    RESULTS['nonlinear_gain']='1 + u + v'

if __name__=='__main__':
    constraint_checks()
    for d in (3,5,7): exchange_checks(d)
    profile_checks()
    polynomial_checks()
    RESULTS['status']='PASS'
    RESULTS['scope']='Finite algebra and symbolic identities; no continuum or empirical inference.'
    Path(__file__).with_name('revision_checks.json').write_text(json.dumps(RESULTS,indent=2)+'\n')
    print(json.dumps(RESULTS,indent=2))
