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
# Channels 0-2 are vertex channels, 3-4 edge channels. Vertex events carry
# no response weight; the source is fed by the edge innovations alone.
vertex=np.array([True,True,True,False,False])
edge=~vertex
w=rng.uniform(-1.8,1.8,(n,c))
w[:,vertex]=0.
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
relative('source_reduces_to_edge_sum',var_s,pref*np.sum((w*w*F*F*dA)[:,edge]))
# Conditional target inequality for the bilocal activity, and saturation.
# The bound presupposes |F_R|<=1; these normalized positive kernels satisfy it.
assert float(np.max(np.abs(F))) <= 1+1e-13
wstar=float(np.max(np.abs(w)))
A_bi=float(dA[:,edge].sum()); A_tot=float(dA.sum())
upper=pref*wstar*wstar*A_bi
assert var_s <= upper*(1+1e-13)
sigma2=var_s*V/(8*np.pi)**2
lower=(sigma2/wstar**2)*V/LP4
assert A_tot >= A_bi >= lower*(1-1e-13)
res['target_bilocal_activity_bound_ratio']=float(lower/A_bi)
res['target_total_activity_bound_ratio']=float(lower/A_tot)
var_sat=pref*wstar*wstar*A_bi
sigma2_sat=var_sat*V/(8*np.pi)**2
relative('target_bound_saturation',(sigma2_sat/wstar**2)*V/LP4,A_bi)

# Vertex neutrality and unentangled laws: with every edge intensity zero the
# regional variance is a sum of nonnegative vertex terms. It vanishes for
# zero vertex weights and is positive for any nonzero one on an active vertex.
dA_unent=dA.copy(); dA_unent[:,edge]=0.
w_general=rng.uniform(.2,1.8,(n,c))*rng.choice([-1.,1.],(n,c))
vertex_terms=pref*(w_general*w_general*F*F*dA_unent).sum(axis=0)
assert np.all(vertex_terms[vertex]>0) and np.all(vertex_terms[edge]==0)
res['unentangled_variance_general_vertex_weights']=float(vertex_terms.sum())
res['unentangled_variance_neutral_vertices']=float(pref*np.sum(w*w*F*F*dA_unent))
assert res['unentangled_variance_neutral_vertices']==0.
for k in np.flatnonzero(vertex):
    single=np.zeros((n,c)); single[:,k]=.7
    assert pref*np.sum(single*single*F*F*dA_unent)>0
res['one_active_vertex_weight_suffices']=True

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

# Constant quadratic density of the edge activity gives the volume scaling.
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
