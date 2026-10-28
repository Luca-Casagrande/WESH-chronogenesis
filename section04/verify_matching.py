"""Algebraic consistency checks for the Section 4 matching identities.

The finite arrays below check changes of integration variable and covariance
bookkeeping. They are not simulations of WESH, cosmology, or observations.
"""
from pathlib import Path
import json
import numpy as np

WORK = Path(__file__).resolve().parent
rng = np.random.default_rng(20261014)
res = {}

def relative(name, left, right, tol=5e-13):
    a, b = np.asarray(left), np.asarray(right)
    err = float(np.max(np.abs(a-b)) / max(1.,float(np.max(np.abs(a))),float(np.max(np.abs(b)))))
    assert err < tol, (name,err)
    res[name] = err

# Arbitrary positive piecewise-constant intensity and matching data.
n, c, z = 11, 5, 7
ell = rng.uniform(.4, 2., n)
ds = rng.uniform(.01,.3,n)
lam = rng.uniform(.02, 3., (n,c))
tauT = .73
G, hbar = .37, 1.31
LP4=(G*hbar)**2
jac = rng.uniform(.2,2.,n) # d tau / ds: independent effective matching
p=lam/lam.sum(axis=1)[:,None]
dA=ell[:,None]*lam*ds[:,None]
dt=tauT*dA.sum(axis=1)
dtau=jac*ds
beta=tauT*ell*lam.sum(axis=1)/jac
lam_tau=ell[:,None]*lam/jac[:,None]
relative('sum_channel_proportions',p.sum(axis=1),np.ones(n))
relative('time_changed_compensators',dA,p*dt[:,None]/tauT)
relative('proper_time_compensators',dA,lam_tau*dtau[:,None])
relative('jacobian_identity',beta,tauT*lam_tau.sum(axis=1))

# Normalized regularized kernels on a finite set of volume elements.
dVz=rng.uniform(.1, 1.1,z)
f=rng.uniform(.01,1.,(z,n,c))
f/=np.einsum('z,znc->nc',dVz,f)[None,:,:]
relative('kernel_normalization',np.einsum('z,znc->nc',dVz,f),np.ones((n,c)))
w=rng.uniform(-1.8,1.8,(n,c))
selected=np.array([True,False,True,True,False,True,False])
V=float(dVz[selected].sum())
F=np.einsum('z,znc->nc',dVz[selected],f[selected])
W2=np.sum(p*w*w*F*F,axis=1)
pref=(8*np.pi*G*hbar)**2/(V*V)
var_s=pref*np.sum(w*w*F*F*dA)
var_t=pref/tauT*np.sum(W2*dt)
var_tau=pref/tauT*np.sum(beta*W2*dtau)
relative('regional_identity_s_vs_t',var_s,var_t)
relative('regional_identity_s_vs_tau',var_s,var_tau)
K=(8*np.pi*G*hbar)**2*np.einsum('znc,ync,nc->zy',f,f,w*w*dA)
v=np.where(selected,dVz,0)
relative('spatial_fubini_covariance_check',v@K@v/(V*V),var_s)
# Conditional target inequality and its saturation case.
wstar=float(np.max(np.abs(w)))
upper=pref*wstar*wstar*float(dA.sum())
assert var_s <= upper*(1+1e-13)
sigma2=var_s*V/(8*np.pi)**2
lower=(sigma2/wstar**2)*V/LP4
assert float(dA.sum()) >= lower*(1-1e-13)
res['target_activity_bound_ratio']=float(lower/dA.sum())
var_sat=pref*wstar*wstar*float(dA.sum())
sigma2_sat=var_sat*V/(8*np.pi)**2
relative('target_bound_saturation',(sigma2_sat/wstar**2)*V/LP4,dA.sum())

# Exact local/global calibration with bilocal activity retained.
N=4
local=rng.uniform(.2,2.,n)
bi=rng.uniform(.1,1.,(n,3))
gamma=rng.uniform(.1,3.,n)
dtau0=tauT*ell*local*ds
dA_bi=(ell[:,None]*bi*ds[:,None]).sum(axis=1)
dt0=tauT*ell*(N*local+bi.sum(axis=1))*ds
relative('global_time_includes_bilocal',dt0,N*dtau0+tauT*dA_bi)
relative('local_proper_rate',ell*local*ds/dtau0,np.ones(n)/tauT)
relative('common_rate_cancellation',bi/(tauT*local[:,None]),(gamma[:,None]*bi)/(tauT*(gamma*local)[:,None]))
res['bilocal_increment_strictly_positive']=bool(np.all(tauT*dA_bi>0))

# Constant local quadratic density gives the stated volume scaling.
Ncells, duration, kappa, weight = 17, 4.2, .81, -.43
LP=np.sqrt(G*hbar)
V4=Ncells*LP**3*duration
Q=Ncells*(kappa/LP)*duration*weight**2
relative('local_extensivity_Q',Q,(kappa*weight**2)*V4/LP**4)
relative('local_volume_variance',(8*np.pi*G*hbar)**2*Q/V4**2,(8*np.pi)**2*kappa*weight**2/V4)

# A signed kernel can have unit integral but |F_R|>1.
res['signed_unit_integral_does_not_bound_subregion']=bool(sum([2.,-1.])==1. and abs(2.)>1)
# Anticipating selection of a terminal Poisson realization is not centered.
a=.8
res['anticipating_window_mean_counterexample']=float(a*np.exp(-a))


(WORK/'matching_checks.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
