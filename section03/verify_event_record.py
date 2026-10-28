"""Independent checks of the Section 3 reviewer claims.
Uses exact finite-state algebra and deterministic matrix computations, not
Monte Carlo trajectories or an observational fit.
"""
from pathlib import Path
import json
import numpy as np
import sympy as sp

OUT = Path(__file__).resolve().parent
p = sp.symbols('p', real=True, positive=True)
cx, cy = sp.symbols('c_x c_y', real=True, positive=True)
I2 = sp.eye(2)
J = sp.diag(1,0)
Jx, Jy = sp.kronecker_product(J,I2), sp.kronecker_product(I2,J)
I4 = sp.eye(4)
H = sp.Matrix([[1,1],[1,-1]]) / sp.sqrt(2)
ket = sp.Matrix([sp.sqrt(p),0,0,sp.sqrt(1-p)])
# Real transpose suffices for these real states, with 0<p<1 understood.

def simp(x):
    return sp.factor(sp.simplify(x))

def rate_generator(psi):
    rho = psi*psi.T
    py = simp(sp.trace(Jy*rho))
    f = cy*py*(1-py)
    drift = sp.zeros(4)
    jump_sum=0
    pieces=[]
    for A,c in [(Jx,cx),(Jy,cy)]:
        j=simp(sp.trace(A*rho))
        delta=A-j*I4
        var=simp(sp.trace(rho*delta**2))
        dB=-c*(delta**2*rho+rho*delta**2)/2+c*var*rho
        drift += dB
        post=sp.simplify(delta*rho*delta/var)
        post_p=simp(sp.trace(Jy*post))
        jump_sum += c*var*(cy*post_p*(1-post_p)-f)
        pieces.append({'j':str(j),'var':str(var),'post_p_y':str(post_p)})
    direct=simp(cy*(1-2*py)*sp.trace(Jy*drift)+jump_sum)
    return direct,pieces

A0,pieces0=rate_generator(ket)
AH,piecesH=rate_generator(sp.kronecker_product(H,I2)*ket)
v=p*(1-p); z=1-2*p
assert sp.simplify(A0+cy*(cx+cy)*v*z*z)==0
assert sp.simplify(AH+cy*cy*v*z*z)==0
assert sp.simplify(AH-A0-cx*cy*v*z*z)==0

# Exactly the same statement for the associated count at small auxiliary s.
# E N_y(s)=c_y v s+0.5 A lambda_y(0) s^2+o(s^2).
countdiff=simp((AH-A0)/2)

# Bell sector checks, throughout the original nonlinear coefficient closure.
I=np.eye(2);Z=np.diag([1.,-1.]); pauliP=(I+Z)/2
Ax=np.kron(pauliP,I); Ay=np.kron(I,pauliP); L=Ax-Ay
basis=np.eye(4)
kets=[(basis[0]+basis[3])/np.sqrt(2), (basis[0]-basis[3])/np.sqrt(2),
      (basis[1]+basis[2])/np.sqrt(2), (basis[1]-basis[2])/np.sqrt(2)]
states=[np.outer(k,k.conj()) for k in kets]
def ptr(rho,keep):
    t=rho.reshape(2,2,2,2)
    return np.trace(t,axis1=1,axis2=3) if keep==0 else np.trace(t,axis1=0,axis2=2)
def gate(rho):
    diff=rho-np.kron(ptr(rho,0),ptr(rho,1))
    return np.abs(np.linalg.eigvalsh(diff)).sum()/2

def Var(r,A):
    return float(np.real(np.trace(A@A@r)-np.trace(A@r)**2))

def D(A,r):
    return A@r@A-(A@A@r+r@A@A)/2

maxres=0
sigma_results=[]
for q in [0.,.1,.25,.6,.99,1.]:
    probs=np.array([(1+q)/4]*2+[(1-q)/4]*2)
    rho=sum(w*r for w,r in zip(probs,states))
    maxres=max(maxres,abs(gate(rho)-q/2),np.linalg.norm(rho-(np.eye(4)+q*np.kron(Z,Z))/4))
    coeff=[1.,1.,2.3*gate(rho)]
    drifts=[]; rates=[]
    for r in states:
        drift=np.zeros((4,4));vr=[]
        for c,A in zip(coeff,[Ax,Ay,L]):
            d=A-np.trace(A@r)*np.eye(4);vv=Var(r,A)
            drift+=-c*(d@d@r+r@d@d)/2+c*vv*r
            vr.append(c*vv)
        drifts.append(float(np.linalg.norm(drift)))
        rates.append(vr)
    maxres=max(maxres,max(drifts),np.linalg.norm(sum(c*D(A,rho) for c,A in zip(coeff,[Ax,Ay,L]))))
    meanbi=float(np.dot(probs,np.array(rates)[:,2]))
    maxres=max(maxres,abs(meanbi-2.3*q*(1-q)/4))
    # Local calibration d tau = (c/4) ds, tau_T=1 and block volume=2.
    wloc=.3;wbi=.4;r_bi=2.3
    per_sector=[wloc**2, wloc**2+r_bi*q*wbi**2]
    sigma_eff=(1+q)/2*per_sector[0]+(1-q)/2*per_sector[1]
    expected=wloc**2+.5*r_bi*wbi**2*q*(1-q)
    maxres=max(maxres,abs(sigma_eff-expected))
    sigma_results.append({'q':q,'gate_global':float(gate(rho)), 'conditional_gate_even':float(gate((states[0]+states[1])/2)),
         'conditional_gate_odd':float(gate((states[2]+states[3])/2)), 'mean_bilocal_aux_rate':meanbi,
         'amplitude_even_fixed_global_law':per_sector[0], 'amplitude_odd_fixed_global_law':per_sector[1], 'amplitude_law_average':sigma_eff})
assert maxres<1e-12

# Explicit demonstration that conditioning the covariance != changing the gate.
# q=.6 has global gate .3 but an independently reclosed odd law has gate .5.
q=sp.symbols('q', nonnegative=True)
r,wl,wb=sp.symbols('r w_l w_b', real=True)
old=wl**2+r*wb**2*q*(1-q)/2
new_mix=wl**2+r*wb**2*(1-q)/2
assert sp.simplify(new_mix-old-r*wb**2*(1-q)**2/2)==0

# Branch preservation is weaker than commutation with the constraints.
C=np.diag([0.,1.]);P=np.diag([1.,0.])
K0=np.diag([1.,0.]);K1=np.array([[0.,1.],[0.,0.]])
assert np.linalg.norm(K0.T@K0+K1.T@K1-np.eye(2))<1e-14
assert np.linalg.norm((np.eye(2)-P)@K1@P)<1e-14
assert np.linalg.norm(C@K1-K1@C)>0.5

out={'exact':{'rate_derivative_original':str(A0),'rate_derivative_hadamard':str(AH),
 'rate_derivative_difference':str(simp(AH-A0)), 'count_difference_s_squared_coefficient':str(countdiff),
 'update_data_original':pieces0,'update_data_hadamard':piecesH,
 'difference_reclosed_vs_original_mixture_variance':str(sp.factor(new_mix-old))},
 'bell_samples':sigma_results,'max_bell_residual':maxres,
 'constraint_commutant_sufficient_not_necessary':True,
 'remarks': ['The total conditional count-rate derivative includes both local channels.',
 'These are algebraic checks with no Monte Carlo trajectories.',
 'All comparisons are at fixed auxiliary exposure; no operational-signalling conclusion is asserted.']}
(OUT/'review_checks.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out['exact'],indent=2));print('Maximum Bell residual:',maxres)
