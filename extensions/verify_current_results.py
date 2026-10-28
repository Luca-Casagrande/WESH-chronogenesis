#!/usr/bin/env python3
"""Regression checks of results added after the first reproducibility bundle.

Targets: physical algebra and product embeddings (Appendix A); the sector-weight
identities (B); finite flow bounds (C); activity/purity; Bell-law relaxation (F);
and cross-epoch inherited radiation. Exact finite identities and selected
numerical examples supplement the manuscript's proofs. They do not establish
the infinite-time martingale theorems, a continuum limit, or new empirical fits.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag, expm

ROOT = Path(__file__).resolve().parent.parent
RNG = np.random.default_rng(20261028)
CHECKS: list[dict] = []
DIAGNOSTICS: dict = {}


def check(name: str, residual: float, tolerance: float = 2e-10) -> None:
    value = float(abs(residual))
    if not np.isfinite(value) or value > tolerance:
        raise AssertionError(f'{name}: {value} > {tolerance}')
    CHECKS.append({'name': name, 'passed': True, 'residual': value,
                   'tolerance': tolerance})


def exact(name: str, expr) -> None:
    terms = list(expr) if isinstance(expr, sp.MatrixBase) else [expr]
    if any(sp.simplify(sp.expand_complex(t)) != 0 for t in terms):
        raise AssertionError(name)
    CHECKS.append({'name': name, 'passed': True, 'method': 'exact symbolic'})


def tn(a):
    return float(np.linalg.svd(a, compute_uv=False).sum())


def adj(a):
    return a.conj().T


def comm(a, b):
    return a @ b - b @ a


def rand_h(d):
    a = RNG.normal(size=(d,d)) + 1j*RNG.normal(size=(d,d))
    return (a + adj(a))/2


def density(d):
    a = RNG.normal(size=(d,d)) + 1j*RNG.normal(size=(d,d))
    a = a @ adj(a)
    return a/np.trace(a)


def pure(d):
    p = RNG.normal(size=d) + 1j*RNG.normal(size=d)
    p /= np.linalg.norm(p)
    return np.outer(p,p.conj())


def diss(j, r):
    return j @ r @ j - (j @ j @ r + r @ j @ j)/2


def gate(r):
    t = r.reshape(2,2,2,2)
    x = np.trace(t,axis1=1,axis2=3)
    y = np.trace(t,axis1=0,axis2=2)
    return .5*tn(r-np.kron(x,y))


I2 = np.eye(2, dtype=complex)
Z = np.diag([1.,-1.]).astype(complex)
X = np.array([[0.,1.],[1.,0.]], dtype=complex)
JX = np.kron((I2+Z)/2,I2)
JY = np.kron(I2,(I2+Z)/2)
JS = [JX,JY,JX-JY]
E = np.eye(4,dtype=complex)
BELL = [(E[:,0]+E[:,3])/np.sqrt(2), (E[:,0]-E[:,3])/np.sqrt(2),
        (E[:,1]+E[:,2])/np.sqrt(2), (E[:,1]-E[:,2])/np.sqrt(2)]
PIS = [np.outer(v,v.conj()) for v in BELL]


def physical_algebra():
    om = -sp.Rational(1,2)+sp.I*sp.sqrt(3)/2
    inds = [-1,0,1]
    f = sp.Matrix(3,3,lambda n,m: om**((inds[n]*inds[m])%3))/sp.sqrt(3)
    t = sp.I/sp.sqrt(3)*sp.Matrix([[0,1,-1],[-1,0,1],[1,-1,0]])
    j = sp.Matrix([[2,-1,-1],[-1,2,-1],[-1,-1,2]])/3
    w = sp.Matrix.hstack(sp.Matrix([-1,2,-1])/sp.sqrt(6),
                         sp.ones(3,1)/sp.sqrt(3))
    exact('A: declared Fourier convention', f.H*sp.diag(-1,0,1)*f-t)
    exact('A: square of temporal matrix', t*t-j)
    exact('A: even-sector isometry', w.H*w-sp.eye(2))
    exact('A: two-level restriction', j*w-w*sp.diag(1,0))
    rev = sp.Matrix([[0,0,1],[0,1,0],[1,0,0]])
    exact('A: even image', rev*w-w)
    exact('A: orientation reversal',rev*t*rev+t)
    p = np.diag([-1.,0.,1.]); q = np.diag([-1.,1.])
    v = np.zeros((18,6),complex)
    for n in range(3):
        for k in range(2):
            v[(2*n+k)*3+n,2*n+k]=1
    c = np.kron(np.kron(-p,I2),np.eye(3))+np.kron(np.eye(6),p)
    proj = v@adj(v)
    b1,b2 = rand_h(6),rand_h(6)
    lift = lambda b:v@b@adj(v)
    check('A: complete algebra multiplication',tn(lift(b1@b2)-lift(b1)@lift(b2)))
    check('A: reduced inverse',tn(adj(v)@lift(b1)@v-b1))
    check('A: zero-extended unit',tn(lift(np.eye(6))-proj))
    check('A: Dirac commutation',tn(comm(c,lift(b1))))
    jn=np.array(j,complex);wn=np.array(w,complex)
    charge=np.kron(np.kron(np.eye(3),q),np.eye(3))
    check('A: neutral event subalgebra',tn(comm(charge,lift(np.kron(jn,I2)))))
    check('A: bare compression',tn(adj(v)@np.kron(np.eye(6),jn)@v-(2/3)*np.eye(6)))
    enc=v@np.kron(wn,np.array([[1.],[0.]]))
    enc2=np.kron(enc,enc)
    for i in range(3):
        r=density(4);tens=r.reshape(2,2,2,2)
        rx=np.trace(tens,axis1=1,axis2=3);ry=np.trace(tens,axis1=0,axis2=2)
        rp=enc2@r@adj(enc2)
        ax=enc@rx@adj(enc);ay=enc@ry@adj(enc)
        # Compute the full encoded pair norm, not a reduced-space shortcut.
        check(f'A: mixed-state gate isometry {i}',
              .5*tn(rp-np.kron(ax,ay))-gate(r),2e-9)


def sector_weights():
    p=np.diag([1.,1.,0.,0.]).astype(complex)
    h=block_diag(rand_h(2),rand_h(2))
    js=[block_diag(rand_h(2),rand_h(2)) for _ in range(2)]
    cs=[.7,.3];ident=np.eye(4)
    for i in range(15):
        pi=pure(4);w=float(np.trace(p@pi).real)
        drift=-1j*comm(h,pi); jump_mean=0.;jump_square=0.;jump_f=0.;jump_w2=0.
        for j,c in zip(js,cs):
            delta=j-np.trace(j@pi).real*ident
            variance=float(np.trace(delta@delta@pi).real)
            update=delta@pi@delta/variance
            dw=float(np.trace(p@update).real)-w
            rate=c*variance
            drift+=c*(-(delta@delta@pi+pi@delta@delta)/2+variance*pi)
            jump_mean+=rate*dw;jump_square+=rate*dw**2
            jump_w2+=rate*((w+dw)**2-w**2)
            jump_f+=rate*((w+dw)*(1-w-dw)-w*(1-w))
        dw_drift=float(np.trace(p@drift).real)
        check(f'B: compensated sector drift {i}',dw_drift+jump_mean)
        check(f'B: quadratic-variation density {i}',
              2*w*dw_drift+jump_w2-jump_square)
        check(f'B: uncertainty budget generator {i}',
              (1-2*w)*dw_drift+jump_f+jump_square)
    for vector in [E[:,0],E[:,3]]:
        pi=np.outer(vector,vector.conj());w=float(np.trace(p@pi).real)
        for j in js:
            d=j-np.trace(j@pi).real*ident;var=float(np.trace(d@d@pi).real)
            if var>1e-14:
                check('B: sharp-sector event invariance',np.trace(p@d@pi@d/var).real-w)


def purity_activity():
    maximum=0.
    for count in [1,3,7]:
        for _ in range(10):
            pis=[pure(4) for _ in range(count)]
            probs=RNG.dirichlet(np.ones(count))
            rho=sum(p*pi for p,pi in zip(probs,pis))
            h=.1*rand_h(4);cs=[.6,.6,.9*gate(rho)];tau=.7
            dr=-1j*comm(h,rho)+sum(c*diss(j,rho) for j,c in zip(JS,cs))
            activity=0.;spread=0.
            for j,c in zip(JS,cs):
                for p,pi in zip(probs,pis):
                    v=float((np.trace(pi@j@j)-np.trace(pi@j)**2).real)
                    activity+=tau*c*p*v
                    spread+=tau/2*c*p*np.linalg.norm(comm(j,pi-rho))**2
            loss=-tau*np.trace(rho@dr).real
            residual=abs(activity-loss-spread);maximum=max(maximum,residual)
            check('Activity-purity identity',residual)
            check('Nonnegative purity contribution',max(0.,-loss))
    q=.4;probs=np.array([(1+q)/4]*2+[(1-q)/4]*2)
    rho=sum(p*pi for p,pi in zip(probs,PIS));nu=.7;gam=.8;tau=.9
    cs=[nu,nu,gam*q/2]
    check('Stationary Bell purity loss',tn(sum(c*diss(j,rho) for j,c in zip(JS,cs))))
    activity=tau*sum(p*c*(np.trace(pi@j@j)-np.trace(pi@j)**2).real
                      for p,pi in zip(probs,PIS) for j,c in zip(JS,cs))
    check('Stationary Bell continuing activity',activity-tau*(nu/2+gam*q*(1-q)/4))
    DIAGNOSTICS['activity_purity_max_residual']=maximum


def superop(js,cs,h):
    d=len(h);ident=np.eye(d)
    a=-1j*(np.kron(ident,h)-np.kron(h.T,ident))
    for j,c in zip(js,cs):
        jj=j@j
        a+=c*(np.kron(j.conj(),j)-.5*np.kron(ident,jj)-.5*np.kron(jj.T,ident))
    return a


def flow_control():
    nu=.4;gam=.25;h=.07*(np.kron(X,I2)+np.kron(I2,X));V=.8
    mstar=2*np.linalg.norm(h,2)+2*sum(c*np.linalg.norm(j,2)**2
                                         for c,j in zip([nu,nu,gam],JS))
    kappa=3*gam*np.linalg.norm(JS[2],2)**2
    def generator(r):return superop(JS,[nu,nu,gam*gate(r)],h)
    def rhs(t,x):
        r=x.reshape(4,4,order='F');return generator(r)@x
    rho0=density(4)
    sol=solve_ivp(rhs,(0,V),rho0.reshape(-1,order='F'),rtol=1e-11,atol=1e-13,
                  method='DOP853',dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    reference=lambda t:sol.sol(t).reshape(4,4,order='F')
    errors=[]
    for n in [4,8,16,32]:
        step=V/n;r=rho0.copy();err=0.
        for i in range(n):
            a=generator(r)
            mid=(expm(step*a/2)@r.reshape(-1,order='F')).reshape(4,4,order='F')
            err=max(err,tn(mid-reference((i+.5)*step)))
            r=(expm(step*a)@r.reshape(-1,order='F')).reshape(4,4,order='F')
            err=max(err,tn(r-reference((i+1)*step)))
        bound=.5*mstar*step*np.expm1(kappa*V)
        check(f'C: frozen approximation bound n={n}',max(0.,err-bound),2e-9)
        check(f'C: trace n={n}',np.trace(r)-1)
        check(f'C: positivity n={n}',max(0.,-np.linalg.eigvalsh((r+adj(r))/2).min()))
        errors.append({'steps':n,'error':err,'upper_bound':bound})
    if errors[-1]['error']>=errors[0]['error']:raise AssertionError('refinement convergence')
    eta0=.8*rho0+.2*np.eye(4)/4
    sol2=solve_ivp(rhs,(0,V),eta0.reshape(-1,order='F'),rtol=1e-11,atol=1e-13,
                   method='DOP853',dense_output=True)
    if not sol2.success:raise RuntimeError(sol2.message)
    for t in np.linspace(0,V,21):
        distance=tn(reference(t)-sol2.sol(t).reshape(4,4,order='F'))
        check('C: gate-controlled initial-data bound',
              max(0.,distance-np.exp(kappa*t)*tn(rho0-eta0)),2e-9)
    # Measurable path with one discontinuity; exact cell averages of c_bi.
    cut=.317;c1=.08;c2=.23
    A1=superop(JS,[nu,nu,c1],h);A2=superop(JS,[nu,nu,c2],h)
    U=expm((V-cut)*A2)@expm(cut*A1)
    joint=RNG.normal(size=(8,8))+1j*RNG.normal(size=(8,8));joint/=tn(joint)
    def extend(a,x):
        tensor=x.reshape(4,2,4,2);out=np.empty_like(tensor)
        for i in range(2):
            for j in range(2):
                out[:,i,:,j]=(a@tensor[:,i,:,j].reshape(-1,order='F')).reshape(4,4,order='F')
        return out.reshape(8,8)
    product_checks=[]
    for n in [4,8,16,32]:
        step=V/n;Up=np.eye(16,dtype=complex);eps=0.
        for i in range(n):
            left=i*step;right=(i+1)*step
            d1=max(0.,min(right,cut)-left);d2=step-d1
            avg=(d1*c1+d2*c2)/step
            Up=expm(step*superop(JS,[nu,nu,avg],h))@Up
            eps+=2*np.linalg.norm(JS[2],2)**2*(d1*abs(avg-c1)+d2*abs(avg-c2))
        error=tn(extend(Up-U,joint))
        check(f'C: averaged-product finite spectator n={n}',max(0.,error-eps),2e-9)
        product_checks.append({'steps':n,'spectator_test_error':error,'epsilon':eps})
    DIAGNOSTICS['flow']={'M_star':mstar,'kappa_star':kappa,'frozen':errors,
                         'averaged_products':product_checks}


def bell_relaxation():
    nu=.7;gam=.9;records=[]
    for q in [0.,.3,.8,1.]:
        even=(1+q)/2;odd=(1-q)/2
        p0=np.array([.8*even,.2*even,.15*odd,.85*odd])
        target=np.array([even/2]*2+[odd/2]*2)
        def rhs(v,y):
            p=y[:4];c=.5*np.abs(p-.25).sum()
            k1=nu/2;k2=nu/2+gam*c
            return [k1*(p[1]-p[0]),k1*(p[0]-p[1]),
                    k2*(p[3]-p[2]),k2*(p[2]-p[3]),c]
        grid=np.linspace(0,8,81)
        sol=solve_ivp(rhs,(0,8),np.r_[p0,0.],t_eval=grid,rtol=2e-11,atol=2e-13,
                      method='DOP853')
        if not sol.success:raise RuntimeError(sol.message)
        p=sol.y[:4];integral=sol.y[4]
        tv=.5*np.abs(p-target[:,None]).sum(axis=0)
        bound=np.exp(-nu*grid)*.5*np.abs(p0-target).sum()
        check(f'F: parity-weight conservation q={q}',np.max(abs(p[:2].sum(axis=0)-even)))
        check(f'F: normalization q={q}',np.max(abs(p.sum(axis=0)-1)))
        check(f'F: total-variation relaxation q={q}',max(0.,np.max(tv-bound)),2e-9)
        check(f'F: even imbalance q={q}',np.max(abs(p[0]-p[1]-(p0[0]-p0[1])*np.exp(-nu*grid))),2e-9)
        check(f'F: odd imbalance q={q}',np.max(abs(p[2]-p[3]-(p0[2]-p0[3])*np.exp(-nu*grid-2*gam*integral))),2e-9)
        records.append({'q':q,'initial_gate':float(.5*np.abs(p0-.25).sum()),
                        'final_gate':float(.5*np.abs(p[:,-1]-.25).sum()),
                        'stationary_gate':q/2,'final_TV':float(tv[-1]),'upper_TV':float(bound[-1])})
    DIAGNOSTICS['bell_relaxation']=records


def inherited_memory():
    a,astar,rhor,kappa,u = sp.symbols('a a_star rho_r kappa u_perp',positive=True)
    U=kappa*rhor*u*(astar/a)**4
    rho1=sp.Rational(7,43)*rhor*(astar/a)**4
    exact('Cross-epoch fixed amplitude',U/(kappa*rho1)-sp.Rational(43,7)*u)
    exact('Cross-epoch derivative',sp.diff(U/(kappa*rho1),a))
    X2,T,F=sp.symbols('X2 T F',positive=True)
    variance=sp.Rational(43,7)**2*T**2*X2*F
    covariance=variance*sp.ones(3)
    if covariance.rank()!=1:raise AssertionError('rank-one cross-epoch covariance')
    CHECKS.append({'name':'Rank-one inherited covariance','passed':True,'method':'exact symbolic'})
    # Ensure the documented observational inputs and active numerical module agree.
    import importlib.util
    path=ROOT/'section09/verify_wesh_section9.py'
    spec=importlib.util.spec_from_file_location('bbn_source',path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    inputs=json.loads((ROOT/'inputs/BBN_inputs.json').read_text(encoding='utf-8'))
    check('Documented BBN mean',mod.N_HAT-inputs['data']['N_nu_mean'])
    check('Documented BBN width',mod.S_N-inputs['data']['N_nu_error_used_in_gaussian_compression'])
    check('Documented BBN alpha',mod.ALPHA-inputs['inference']['alpha'])
    endpoint=mod.upper_boundary(50)['upper']
    DIAGNOSTICS['BBN_endpoint']=endpoint
    DIAGNOSTICS['continued_source_coefficient']=8*np.pi*endpoint/3


def main():
    for f in [physical_algebra,sector_weights,purity_activity,flow_control,
              bell_relaxation,inherited_memory]:
        f()
    report={'all_passed':True,'seed':20261028,'number_of_checks':len(CHECKS),
            'checks':CHECKS,'diagnostics':DIAGNOSTICS,
            'scope':'Exact finite identities and selected regression cases for the current manuscript; analytic proofs remain in the manuscript.'}
    (Path(__file__).with_name('current_results_checks.json')).write_text(
        json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'all_passed':True,'checks':len(CHECKS),'diagnostics':DIAGNOSTICS},indent=2))


if __name__=='__main__':main()
