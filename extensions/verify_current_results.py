#!/usr/bin/env python3
"""Regression checks of results added after the first reproducibility bundle.

Targets: physical algebra and product embeddings (Appendix A); the sector-weight
identities (B); coupling and flow bounds (C); activity/purity; Bell-law relaxation (F);
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


def negativity(pair, dx, dy):
    """Negativity of a pair state, capped at one."""
    pt = pair.reshape(dx,dy,dx,dy).transpose(2,1,0,3).reshape(dx*dy,dx*dy)
    return min(1., tn(pt)-1.)


def upsilon(pair, dx=2, dy=2):
    """Negativity-weighted pair correlation of one conditional pair state."""
    t = pair.reshape(dx,dy,dx,dy)
    x = np.trace(t,axis1=1,axis2=3)
    y = np.trace(t,axis1=0,axis2=2)
    return negativity(pair,dx,dy)*(pair-np.kron(x,y))


def gate(probs, pures, dx=2, dy=2):
    """Entanglement gate of a law of pure pair states."""
    return .5*tn(sum(w*upsilon(r,dx,dy) for w,r in zip(probs,pures)))


def bell_gate(p):
    """Gate of a Bell-supported law with weights p: unit negativity and
    maximally mixed marginals reduce it to (1/2) sum_k |p_k - 1/4|."""
    return .5*np.abs(np.asarray(p)-.25).sum()


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
    dreg=enc.shape[0]
    for i in range(3):
        r=density(4);tens=r.reshape(2,2,2,2)
        rx=np.trace(tens,axis1=1,axis2=3);ry=np.trace(tens,axis1=0,axis2=2)
        rp=enc2@r@adj(enc2)
        ax=enc@rx@adj(enc);ay=enc@ry@adj(enc)
        # First line of the embedding identity, on the full encoded pair space.
        check(f'A: mixed-state correlation isometry {i}',
              tn(rp-np.kron(ax,ay))-tn(r-np.kron(rx,ry)),2e-9)
        # Second line: the partial transpose keeps its nonzero singular values.
        check(f'A: negativity under embedding {i}',
              negativity(rp,dreg,dreg)-negativity(r,2,2),2e-9)
    probs=RNG.dirichlet(np.ones(3));pures=[pure(4) for _ in range(3)]
    embedded=[enc2@r@adj(enc2) for r in pures]
    check('A: gate of a law equals gate of its embedded image',
          gate(probs,embedded,dreg,dreg)-gate(probs,pures),2e-9)


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
            h=.1*rand_h(4);cs=[.6,.6,.9*gate(probs,pis)];tau=.7
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


def bell_rate_matrix(nu, c_bi):
    """Unit-lapse transition generator on (Phi+,Phi-,Psi+,Psi-)."""
    k1=nu/2;k2=nu/2+c_bi
    return np.array([[-k1,k1,0,0],[k1,-k1,0,0],[0,0,-k2,k2],[0,0,k2,-k2]])


def flow_control():
    nu=.4;gam=.25;V=.8
    norms=[np.linalg.norm(j,2)**2 for j in JS]
    lam=nu*norms[0]+nu*norms[1]+gam*norms[2]
    dxy=2
    # --- Constants of the appendix (hbar=1 in these dimensionless checks) ---
    def constants(hnorm):
        return (2*hnorm+7*lam, 2*hnorm+29*lam, 7*gam*(dxy+1.5)*norms[2])
    # --- Per-state estimates behind the coupling inequality ------------------
    worst={'V':0.,'G':0.,'N':0.,'jump':0.,'speed':0.}
    for trial in range(400):
        d=int(RNG.integers(2,6));j=rand_h(d);a=np.linalg.norm(j,2)**2
        p1=pure(d);vec=np.linalg.eigh(p1)[1][:,-1]
        w0=vec+RNG.choice([1e-3,1e-2,1e-1,1.])*(RNG.normal(size=d)+1j*RNG.normal(size=d))
        p2=np.outer(w0,w0.conj())/np.vdot(w0,w0).real
        dist=tn(p1-p2);ident=np.eye(d)
        def parts(pi):
            delta=j-np.trace(j@pi).real*ident
            var=float(np.trace(pi@delta@delta).real)
            return var,delta@pi@delta,-(delta@delta@pi+pi@delta@delta)/2+var*pi
        v1,g1,n1=parts(p1);v2,g2,n2=parts(p2)
        worst['V']=max(worst['V'],abs(v1-v2)-3*a*dist)
        worst['G']=max(worst['G'],tn(g1-g2)-8*a*dist)
        worst['N']=max(worst['N'],tn(n1-n2)-12*a*dist)
        worst['speed']=max(worst['speed'],tn(n1)-5*a)
        if min(v1,v2)>1e-12:
            worst['jump']=max(worst['jump'],min(v1,v2)*tn(g1/v1-g2/v2)-11*a*dist)
    check('C: variance estimate |V-V~| <= 3 a d',max(0.,worst['V']))
    check('C: update estimate <= 8 a d',max(0.,worst['G']))
    check('C: no-event estimate <= 12 a d',max(0.,worst['N']))
    check('C: no-event norm <= 5 a',max(0.,worst['speed']))
    check('C: common-jump estimate <= 11 a d',max(0.,worst['jump']))
    # --- Generator bound with a nonzero Hamiltonian --------------------------
    h=.07*(np.kron(X,I2)+np.kron(I2,X));hnorm=np.linalg.norm(h,2)
    mstar,lstar,kappa=constants(hnorm)
    worst_gen=0.
    for trial in range(20):
        probs=RNG.dirichlet(np.ones(3));pures=[pure(4) for _ in range(3)]
        a=superop(JS,[nu,nu,gam*gate(probs,pures)],h)
        for _ in range(5):
            x=RNG.normal(size=(4,4))+1j*RNG.normal(size=(4,4))
            y=(a@x.reshape(-1,order='F')).reshape(4,4,order='F')
            worst_gen=max(worst_gen,tn(y)/tn(x)-mstar)
    check('C: generator trace-norm bound',max(0.,worst_gen))
    # --- Self-consistent flow on Bell-supported laws (zero Hamiltonian) ------
    # There the no-event flow vanishes and the law is a four-state chain;
    # distinct Bell projectors are at trace distance 2, so W_1 = 2 TV.
    # The constants are convenient, unoptimized upper bounds; a short exposure
    # keeps the exponential factor of the stated bounds moderate.
    mstar,lstar,kappa=constants(0.);VC=.1
    def rhs(v,p):return bell_rate_matrix(nu,gam*bell_gate(p))@p
    def w1(p,q):return float(np.abs(p-q).sum())
    p0=np.array([.45,.15,.07,.33])
    sol=solve_ivp(rhs,(0,VC),p0,rtol=1e-12,atol=1e-14,method='DOP853',dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    errors=[]
    for n in [4,8,16,32]:
        step=VC/n;p=p0.copy();err=0.
        for i in range(n):
            q_gen=bell_rate_matrix(nu,gam*bell_gate(p))
            err=max(err,w1(expm(step*q_gen/2)@p,sol.sol((i+.5)*step)))
            p=expm(step*q_gen)@p
            err=max(err,w1(p,sol.sol((i+1)*step)))
        bound=kappa/(lstar+kappa)*mstar*step*np.expm1((lstar+kappa)*VC)
        check(f'C: frozen approximation bound n={n}',max(0.,err-bound),2e-9)
        check(f'C: frozen law normalization n={n}',p.sum()-1)
        check(f'C: frozen law positivity n={n}',max(0.,-p.min()))
        rho=sum(w*pi for w,pi in zip(p,PIS))
        check(f'C: frozen barycenter positivity n={n}',max(0.,-np.linalg.eigvalsh(rho).min()))
        errors.append({'steps':n,'W1_error':err,'upper_bound':float(bound)})
    if errors[-1]['W1_error']>=errors[0]['W1_error']:raise AssertionError('refinement convergence')
    # Mesh convergence at first order, as the bound predicts.
    check('C: first-order mesh convergence',
          max(0.,errors[-1]['W1_error']/errors[0]['W1_error']-.2))
    q0=np.array([.30,.28,.22,.20])
    sol2=solve_ivp(rhs,(0,VC),q0,rtol=1e-12,atol=1e-14,method='DOP853',dense_output=True)
    if not sol2.success:raise RuntimeError(sol2.message)
    stability=[]
    for tgrid in np.linspace(0,VC,21):
        distance=w1(sol.sol(tgrid),sol2.sol(tgrid))
        check('C: initial-data stability in W_1',
              max(0.,distance-np.exp((lstar+kappa)*tgrid)*w1(p0,q0)),2e-9)
        bary=tn(sum((a-b)*pi for a,b,pi in zip(sol.sol(tgrid),sol2.sol(tgrid),PIS)))
        check('C: barycenter distance bounded by W_1',max(0.,bary-distance))
        stability.append(distance)
    # Coefficient Lipschitz bound on the same pair of laws.
    for tgrid in np.linspace(0,VC,5):
        pa,pb=sol.sol(tgrid),sol2.sol(tgrid)
        check('C: coefficient Lipschitz bound',
              max(0.,gam*abs(bell_gate(pa)-bell_gate(pb))-gam*(dxy+1.5)*w1(pa,pb)))
    # kappa=0: coefficients do not depend on the law and the barycentric flow contracts.
    lin=expm(V*superop(JS,[nu,nu,0.],h))
    r1,r2=density(4),density(4)
    out=lambda r:(lin@r.reshape(-1,order='F')).reshape(4,4,order='F')
    check('C: contraction without bilocal coupling',max(0.,tn(out(r1)-out(r2))-tn(r1-r2)))
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
    DIAGNOSTICS['flow']={'M_star':mstar,'L_star':lstar,'kappa_star':kappa,
                         'coupling_exposure':VC,
                         'frozen':errors,'stability_W1_max':max(stability),
                         'averaged_products':product_checks}


def bell_relaxation():
    nu=.7;gam=.9;records=[]
    for _ in range(10):
        weights=RNG.dirichlet(np.ones(4))
        check('F: Bell-supported gate formula',gate(weights,PIS)-bell_gate(weights))
    for q in [0.,.3,.8,1.]:
        even=(1+q)/2;odd=(1-q)/2
        p0=np.array([.8*even,.2*even,.15*odd,.85*odd])
        target=np.array([even/2]*2+[odd/2]*2)
        def rhs(v,y):
            p=y[:4];c=bell_gate(p)
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
        records.append({'q':q,'initial_gate':float(bell_gate(p0)),
                        'final_gate':float(bell_gate(p[:,-1])),
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
