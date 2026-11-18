#!/usr/bin/env python3
"""Exact checks for the Section 5 regional identity and the aligned endpoint.

Covers the finite-array volume identity, the q=1 law, the uniqueness of the
aligned active stationary sector on the two-site block, and the arithmetic of
the amplitude bound once the documented endpoint is inserted. The volume
assignment, full collection and the product initial law are hypotheses of the
effective realization, not outputs of these algebraic checks.
"""
from pathlib import Path
import importlib.util
import json
import sympy as sp

here=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('base_checks',here/'verify_section5.py')
base=importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
passed=[]
def check(a,b,name):
    diff=a-b
    terms=list(diff) if isinstance(diff,sp.MatrixBase) else [diff]
    if not all(sp.simplify(x)==0 for x in terms):
        raise AssertionError(name)
    passed.append(name)

q,r=sp.symbols('q r_bi', nonnegative=True)
n=sp.symbols('n_b',integer=True,positive=True)
LP,dt,nu=sp.symbols('L_P Delta_tau nu',positive=True)
wb=sp.symbols('w_bi',real=True)
K=sp.symbols('K',integer=True,nonnegative=True)
sig=r*wb**2*q*(1-q)/2
V=2*n*LP**3*dt
rate=2*r*q*wb**2*K/LP
meanrate=sp.simplify(rate.subs(K,n*(1-q)/2))
Qblk=2*sig/LP
check(meanrate,n*Qblk,'Ensemble weighted rate = n_b Q_blk')
conditionalvar=(8*sp.pi*LP**2)**2*dt*rate/V**2
check(conditionalvar,(8*sp.pi)**2*(r*q*wb**2*K/n)/V,
      'Conditional regional variance at fixed parities')
ensemblevar=conditionalvar.subs(K,n*(1-q)/2)
check(ensemblevar,(8*sp.pi)**2*sig/V,'Exact finite-array four-volume identity')
check(sig.subs(q,1),0,'q=1 extinction of the response')
check(sig.subs(q,0),0,'q=0 extinction of the response')
check(base.weights.subs(base.q,1),sp.Matrix([sp.Rational(1,2),sp.Rational(1,2),0,0]),
      'q=1 law has two equiprobable even Bell states')
rho1=base.rho.subs(base.q,1)
check(rho1,sp.diag(sp.Rational(1,2),0,0,sp.Rational(1,2)),
      'q=1 barycenter is a mixed separable state')
check(sp.trace(rho1*rho1),sp.Rational(1,2),'q=1 barycentric purity')
check(base.Q.subs(base.q,1)*base.weights.subs(base.q,1),sp.zeros(4,1),
      'q=1 full-law stationarity')
check((base.tau*sum(base.weights[k]*sum(base.state_rates[k]) for k in range(4))).subs(base.q,1),
      base.tau*base.nu/2,'q=1 mean chronogenetic activity')
# Sum exact binomial parity weights for several finite sizes; no large-N limit.
p=(1-q)/2
for nb in range(1,7):
    weighted=sum(sp.binomial(nb,k)*p**k*(1-p)**(nb-k)*conditionalvar.subs({n:nb,K:k})
                 for k in range(nb+1))
    check(weighted,((8*sp.pi)**2*sig/V).subs(n,nb),f'Finite parity mixture: {nb} block(s)')
# The discarded pathwise equality is not an identity: one generic configuration suffices.
example=(rate-meanrate).subs({n:2,K:0,q:sp.Rational(1,2),r:1,wb:1,LP:1})
assert example!=0
passed.append('Mean relation must not be asserted for each parity realization')

# --- Uniqueness of the aligned active stationary sector on C^2 x C^2 ---------
# Operators of the block: J_x=|0><0| x 1, J_y=1 x |0><0|, J_xy=J_x-J_y, H=0.
Jx,Jy,Jxy,I4=base.Jx,base.Jy,base.L,base.I4
check(Jxy,sp.diag(0,1,-1,0),'Mismatch operator is diag(0,1,-1,0)')
e=[I4[:,k] for k in range(4)]
# Vectors with sharp mismatch are eigenvectors of J_xy: |01>, |10>, or even.
for k in (1,2):
    proj=e[k]*e[k].T
    for J in (Jx,Jy,Jxy):
        m=sp.trace(J*proj)
        check(sp.trace(J*J*proj)-m**2,0,f'Basis state {k} is sharp in every channel')
# Even subspace: psi = cos(th)|00> + sin(th) e^{i phi}|11>, p = cos(th)^2.
th,phi=sp.symbols('theta phi',real=True)
nu_=base.nu
pp=sp.cos(th)**2
psi=sp.cos(th)*e[0]+sp.sin(th)*sp.exp(sp.I*phi)*e[3]
pi_=psi*psi.H
def zero(expr):
    return sp.simplify(sp.expand_complex(sp.expand_trig(sp.simplify(expr))))
def tcheck(a,b,name):
    diff=a-b
    terms=list(diff) if isinstance(diff,sp.MatrixBase) else [diff]
    if not all(zero(x)==0 for x in terms):
        raise AssertionError(name)
    passed.append(name)
def var(J,state):
    m=sp.trace(J*state)
    return sp.trace(J*J*state)-m**2
g=pp*(1-pp)
tcheck(var(Jxy,pi_),0,'Even subspace has sharp mismatch')
tcheck(var(Jx,pi_),g,'Local variance x on the even subspace')
tcheck(var(Jy,pi_),g,'Local variance y on the even subspace')
for J,name in ((Jx,'x'),(Jy,'y')):
    d=J-pp*I4
    post=d*psi                       # unnormalized centered update, norm^2 = g
    tcheck((post.H*post)[0],g,f'Centered update {name} has norm squared p(1-p)')
    tcheck((post.H*Jx*post)[0],g*(1-pp),f'Local update {name} sends p to 1-p')
    tcheck((post.H*Jxy*Jxy*post)[0],0,f'Local update {name} stays in the even subspace')
# No-event flow: d psi/dv = -(1/2) sum_alpha c_alpha (Delta^2 - V) psi, local channels at rate nu.
flow=sp.zeros(4,1)
for J in (Jx,Jy):
    d=J-pp*I4
    flow+=-nu_*(d*d-g*I4)*psi/2
dp=(flow.H*Jx*psi+psi.H*Jx*flow)[0]
tcheck(dp,-2*nu_*pp*(1-pp)*(1-2*pp),'No-event flow of p')
tcheck((1-2*pp)*dp,-2*nu_*pp*(1-pp)*(1-2*pp)**2,'No-event flow of g=p(1-p) is nonpositive')
half=sp.Rational(1,2)
tcheck(flow.subs(th,sp.pi/4),sp.zeros(4,1),'No-event flow vanishes at p=1/2')
aligned=psi.subs(th,sp.pi/4)
proj_aligned=aligned*aligned.H
for J,name in ((Jx,'x'),(Jy,'y')):
    post=2*(J-I4/2)*aligned
    tcheck(post*post.H,proj_aligned.subs(phi,phi+sp.pi),f'Local event {name} shifts phi by pi')
sym_bary=(proj_aligned+proj_aligned.subs(phi,phi+sp.pi))/2
tcheck(sym_bary,sp.diag(half,0,0,half),'Aligned barycenter for every phase')
tcheck(sym_bary,base.rho.subs(base.q,1),'Aligned barycenter equals the q=1 barycenter')
# Unit negativity and maximally mixed marginals: R_ent = barycenter - 1/4, gate 1/2.
pt=sp.Matrix(4,4,lambda i,j:proj_aligned[2*(j//2)+i%2,2*(i//2)+j%2])
tcheck(pt*pt,I4/4,'Aligned partial transpose squares to 1/4: unit negativity')
R=sp.diag(half,0,0,half)-I4/4
check(sum(abs(x) for x in R.eigenvals(multiple=True))/2,half,'Aligned gate equals 1/2')
check(base.tau*2*nu_*half*(1-half),base.tau*nu_/2,'Aligned activity')
# Product states without mismatch are basis states: no time production.
for k in range(4):
    proj=e[k]*e[k].T
    check(sum(var(J,proj) for J in (Jx,Jy,Jxy)),0,f'Computational state {k} has zero activity')

# --- Arithmetic of the amplitude bound with the documented endpoint ----------
import numpy as np
spec9=importlib.util.spec_from_file_location('bbn_source',here.parent/'section09/verify_wesh_section9.py')
bbn=importlib.util.module_from_spec(spec9);spec9.loader.exec_module(bbn)
endpoint=float(bbn.upper_boundary(50)['upper'])
budget=endpoint**2
numbers={'sigma_95':endpoint,'sigma_95_squared':budget,
         'coupling_bound_at_q_half':8*budget,
         'q_times_one_minus_q_bound_at_unit_coupling':2*budget,
         'misaligned_occupation_bound_near_alignment':budget}
assert abs(budget-1.24756976365e-8)<1e-17
assert abs(8*budget-1.0e-7)<5e-9          # manuscript: 1.0e-7
assert abs(2*budget-2.5e-8)<5e-10         # manuscript: 2.5e-8
assert abs(budget-1.25e-8)<5e-11          # manuscript: 1.25e-8
# Exact root of (1/2) q (1-q) = budget near q=1 at unit coupling.
qroot=(1+np.sqrt(1-8*budget))/2
assert abs((1-qroot)/2-budget)<1e-14
numbers['one_minus_q_exact_root_at_unit_coupling']=float(1-qroot)
passed.append('Amplitude-bound arithmetic with the documented endpoint')
report={'all_checks_passed':True,'additional_exact_checks':len(passed),
        'checks':passed,'base_exact_checks':len(base.checks),
        'amplitude_bound_numbers':numbers,
        'scope':'Finite algebra and declared deterministic-volume, full-collection matching. The amplitude-bound numbers reuse the documented endpoint of Section 9; no new data are introduced.'}
(here/'VERIFICHE_AGGIUNTE.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
