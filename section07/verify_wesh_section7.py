#!/usr/bin/env python3
"""Checks of the WESH homogeneous trace-and-conservation closure.

This verifies identities in the declared effective model; it does not validate
its physical matching or supply observational evidence. The geometry check
starts from the flat FRW metric. Numerical tests compare an independently
implemented trace equation with the coupled fluid equations, using arbitrary
dimensionless smooth test sources (not cosmological fits).

Python >=3.10; dependencies: sympy, numpy, scipy.
Run: python verify_wesh_section7.py --output section7_checks.json
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from time import perf_counter
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).with_name('section7_checks.json'))
    args = parser.parse_args()
    start = perf_counter()
    checks: dict[str, bool] = {}
    diagnostics: dict[str, float] = {}

    def zero(name: str, expr: sp.Expr | sp.MatrixBase) -> None:
        if isinstance(expr, sp.MatrixBase):
            checks[name] = all(sp.simplify(v) == 0 for v in expr)
        else:
            checks[name] = sp.simplify(expr) == 0

    def small(name: str, value: float, bound: float = 2e-8) -> None:
        diagnostics[name] = float(value)
        checks[name] = bool(np.isfinite(value) and value < bound)

    t, xx, yy, zz = sp.symbols('t x y z', real=True)
    coords = (t, xx, yy, zz)
    a = sp.Function('a')(t)
    g = sp.diag(-1, a*a, a*a, a*a)
    gi = g.inv()
    gamma = [[[sp.simplify(sp.Rational(1,2)*sum(
        gi[k,l]*(sp.diff(g[l,j],coords[i])+sp.diff(g[l,i],coords[j])
                 -sp.diff(g[i,j],coords[l])) for l in range(4)))
        for j in range(4)] for i in range(4)] for k in range(4)]
    Ric = sp.Matrix(4,4,lambda i,j: sp.simplify(sum(
        sp.diff(gamma[k][i][j],coords[k])-sp.diff(gamma[k][i][k],coords[j])
        +sum(gamma[k][k][l]*gamma[l][i][j]-gamma[k][j][l]*gamma[l][i][k]
             for l in range(4)) for k in range(4))))
    R = sp.simplify(sum(gi[i,j]*Ric[i,j] for i in range(4) for j in range(4)))
    Ht = sp.diff(a,t)/a
    Ein = sp.simplify(Ric-sp.Rational(1,2)*g*R)
    zero('metric_R', R-6*(sp.diff(Ht,t)+2*Ht**2))
    zero('metric_G00', Ein[0,0]-3*Ht**2)
    zero('metric_spatial_Einstein', Ein[1,1]/a**2+2*sp.diff(Ht,t)+3*Ht**2)
    zero('metric_mixed_Einstein', sp.Matrix([Ein[0,i] for i in range(1,4)]))

    U, V, zeta = (sp.Function(name)(t) for name in ('U','V','zeta'))
    S = sp.diag(U,V*a*a,V*a*a,V*a*a)
    def divergence(T: sp.MatrixBase) -> sp.MatrixBase:
        return sp.Matrix([sp.simplify(sum(gi[i,i]*(sp.diff(T[i,j],coords[i])
            -sum(gamma[k][i][i]*T[k,j]+gamma[k][i][j]*T[i,k] for k in range(4)))
            for i in range(4))) for j in range(4)])
    divS = divergence(S)
    zero('fluid_trace', sum(gi[i,i]*S[i,i] for i in range(4))+U-3*V)
    zero('fluid_divergence_time', divS[0]+sp.diff(U,t)+3*Ht*(U+V))
    zero('fluid_divergence_space', divS[1:4,0])
    zero('trace_to_filter', (-divS[0]).subs(V,(U-4*zeta)/3)
         -(sp.diff(U,t)+4*Ht*(U-zeta)))
    D = U-zeta
    C = sp.simplify(S.subs(V,(U-4*zeta)/3)+zeta*g)
    zero('compensating_density', C[0,0]-D)
    zero('compensating_pressure', C[1,1]/a**2-D/3)
    zero('compensating_trace', sum(gi[i,i]*C[i,i] for i in range(4)))
    zero('compensating_divergence', divergence(C)[0]+sp.diff(D,t)+4*Ht*D)
    zero('D_driven_balance', (sp.diff(D,t)+4*Ht*D+sp.diff(zeta,t))
         -(sp.diff(U,t)+4*Ht*(U-zeta)))
    zero('U_integrating_factor', sp.diff(a**4*U,t)/a**4
         -(sp.diff(U,t)+4*Ht*U))
    zero('D_integrating_factor', sp.diff(a**4*D,t)/a**4
         -(sp.diff(D,t)+4*Ht*D))

    h,hd,rho,pp,uu,vs,ze,lam,kappa,mu,nu,dz = sp.symbols(
        'h hd rho pp uu vs ze lam kappa mu nu dz', real=True)
    ud = -4*h*(uu-ze)
    rhod = -3*h*(rho+pp)
    hdot = -kappa*(rho+pp)/2-sp.Rational(2,3)*(uu-ze)
    zero('Friedmann_constraint_propagation', 6*h*hdot-kappa*rhod-ud)
    trace = 6*(hdot+2*h*h)-4*lam-kappa*(rho-3*pp)-4*ze
    zero('Einstein_trace_on_constraint', trace.subs(h*h,(kappa*rho+lam+uu)/3))
    pressure_eq = 2*hdot+3*h*h-lam+kappa*pp+(uu-4*ze)/3
    zero('spatial_Einstein_on_constraint', pressure_eq.subs(h*h,(kappa*rho+lam+uu)/3))

    # Ahmed--Sorkin variables here use kappa=1; their coefficients are distinct
    # from WESH's primitive local rate nu.
    rhoAS,pAS,lambdaAS,DAS = sp.symbols('rhoAS pAS lambdaAS DAS',real=True)
    E0 = 3*h*h-lambdaAS-rhoAS
    Ea = 2*hd+3*h*h-lambdaAS+pAS
    mixed = nu*E0+mu*Ea
    Rtrace = 6*(hd+2*h*h)-4*lambdaAS-rhoAS+3*pAS
    zero('AS_trace_combination', mixed.subs(mu,3*nu)-nu*Rtrace)
    h_AS = sp.solve(mixed,hd)[0]
    dAS = 6*h*h_AS-dz+3*h*(rhoAS+pAS)
    zero('AS_general_residual', sp.factor(dAS+3*h*(1+nu/mu)*E0+dz))
    zero('AS_trace_residual', sp.factor((dAS+4*h*E0+dz).subs(mu,3*nu)))
    zero('AS_dictionary_residual', E0.subs({lambdaAS:lam+ze,rhoAS:kappa*rho})
         .subs(h*h,(kappa*rho+lam+uu)/3)-(uu-ze))
    rhoD,rhoR = sp.symbols('rhoD rhoR',real=True)
    zero('radiation_cancels_trace', (rhoD+rhoR)-3*(rhoR/3)-rhoD)
    zero('conformal_trace_dictionary', (a*sp.diff(a*sp.diff(a,t),t))/a**3
         -(sp.diff(Ht,t)+2*Ht**2))
    zero('trace_equation_reconstructs_U',
         6*h*((4*(lam+ze)+kappa*(rho-3*pp))/6-2*h*h)
         -kappa*rhod+4*h*(3*h*h-kappa*rho-lam)-4*h*ze)

    # First-order and exact reference-fraction algebra.
    ep,Hb,dH,rhob,drho,u1 = sp.symbols('ep Hb dH rhob drho u1',real=True)
    Fres = 3*(Hb+ep*dH)**2-kappa*(rhob+ep*drho)-lam-ep*u1
    zero('linearized_Friedmann', sp.diff(Fres,ep).subs(ep,0)-(6*Hb*dH-kappa*drho-u1))
    ub,xb,w,n,n0,alpha,U0,X = sp.symbols('ub xb w n n0 alpha U0 X', real=True)
    frac_rhs = 4*xb-4*ub+3*(1+w)*ub
    zero('general_fraction_equation',frac_rhs-(4*xb-(1-3*w)*ub))
    for name,w0,alpha0 in [('radiation',sp.Rational(1,3),0),('dust',0,1),('deSitter',-1,4)]:
        zero(name+'_fraction',frac_rhs.subs(w,w0)-(4*xb-alpha0*ub))
        free=U0*sp.exp(-alpha0*(n-n0))
        zero(name+'_free_mode',sp.diff(free,n)+alpha0*free)
        sol=U0+4*X*(n-n0) if alpha0==0 else (
            U0*sp.exp(-alpha0*(n-n0))+4*X/alpha0*(1-sp.exp(-alpha0*(n-n0))))
        zero(name+'_constant_source_solution',sp.diff(sol,n)+alpha0*sol-4*X)
    zero('radiation_variance_prefactor',sp.Rational(1,4)*4**2-4)
    vF = -ze-dz/(3*h)
    zero('instantaneous_conservation',dz+3*h*(ze+vF))
    zero('instantaneous_trace',-ze+3*vF+4*ze+dz/h)
    zero('coincidence_requires_constant_source',dz+4*h*ze-4*h*ze-dz)

    # Independent numerical integration, units kappa=1.
    Lambda0=0.18; dust0=0.34; rad0=0.80; initU=0.01
    initH=np.sqrt((Lambda0+dust0+rad0+initU)/3)
    def source(t: float) -> float:
        return 0.015*np.sin(2.7*t)+0.01*np.cos(0.4*t)
    def rhs(t: float,y: np.ndarray) -> np.ndarray:
        aa,HH,UU,rd,rr,ii=y; zz=source(t);rho=rd+rr;p=rr/3
        return np.array([aa*HH,-(rho+p)/2-2*(UU-zz)/3,-4*HH*(UU-zz),
                         -3*HH*rd,-4*HH*rr,4*aa**4*HH*zz])
    tgrid=np.linspace(0,5,301)
    coupled=solve_ivp(rhs,(0,5),[1,initH,initU,dust0,rad0,0],
                      t_eval=tgrid,rtol=2e-11,atol=2e-13,method='DOP853')
    if not coupled.success: raise RuntimeError(coupled.message)
    aa,hh,us,rd,rr,integral=coupled.y
    small('numeric_Friedmann_residual',np.max(np.abs(3*hh**2-(rd+rr)-Lambda0-us)))
    small('numeric_retarded_integral',np.max(np.abs(us-(initU+integral)/aa**4)))
    small('numeric_dust_normalization',np.max(np.abs(rd*aa**3-dust0)))
    small('numeric_radiation_normalization',np.max(np.abs(rr*aa**4-rad0)))
    def trace_rhs(t: float,y: np.ndarray) -> np.ndarray:
        aa,HH=y
        return np.array([aa*HH,(4*(Lambda0+source(t))+dust0/aa**3)/6-2*HH**2])
    tr=solve_ivp(trace_rhs,(0,5),[1,initH],t_eval=tgrid,
                 rtol=2e-11,atol=2e-13,method='DOP853')
    if not tr.success: raise RuntimeError(tr.message)
    small('numeric_independent_trace_a',np.max(np.abs(tr.y[0]-aa)))
    small('numeric_independent_trace_H',np.max(np.abs(tr.y[1]-hh)))
    uuFromTrace=3*tr.y[1]**2-Lambda0-dust0/tr.y[0]**3-rad0/tr.y[0]**4
    small('numeric_reconstructed_response',np.max(np.abs(uuFromTrace-us)))
    checks['numeric_expanding_domain']=bool(np.min(hh)>0)

    result={'checks':checks,'number_of_checks':len(checks),
            'all_passed':all(checks.values()),'numeric_residuals':diagnostics,
            'seconds':perf_counter()-start,
            'scope':'Algebra and numerical consistency of the stated effective closure; no empirical validation.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if not result['all_passed']: raise SystemExit(1)

if __name__=='__main__':
    main()
