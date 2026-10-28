#!/usr/bin/env python3
"""Check the real-mode normalization and component identities in WESH Appendix E.

Run with Python 3.10+ and SymPy. This script differentiates explicit functions;
it does not establish the global kernel classification, which is proved in text,
nor any stochastic, nonlinear or uniform-in-wavenumber stability result.
The existing section06/verify_wesh_section6.py provides the original independent
metric-variation and five-kernel-basis checks. Neither script supplies physical
matching assumptions.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sympy as s


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, default=Path(__file__).with_name('appendixE_checks.json'))
    args = ap.parse_args()
    checks: dict[str, bool] = {}

    def check(name: str, expr: s.Expr) -> None:
        checks[name] = s.simplify(s.expand_trig(expr)) == 0

    H, k = s.symbols('H k', positive=True)
    eta = s.symbols('eta', real=True, nonzero=True)
    z = k*eta
    j = (3/z**3-1/z)*s.sin(z)-3*s.cos(z)/z**2
    y = (1/z-3/z**3)*s.cos(z)-3*s.sin(z)/z**2
    cj, cy = eta**2*j, eta**2*y
    ode = lambda f: s.diff(f,eta,2)-2*s.diff(f,eta)/eta+(k*k-4/eta**2)*f
    check('chi_j_mode_equation', ode(cj))
    check('chi_y_mode_equation', ode(cy))
    check('real_basis_Wronskian', cj*s.diff(cy,eta)-s.diff(cj,eta)*cy-eta**2/k)
    oldF = -eta*s.exp(-s.I*k*eta)*(1-3*s.I/(k*eta)-3/(k*k*eta**2))
    check('relation_to_existing_complex_basis', s.expand_complex(oldF)+k*(cy+s.I*cj))

    energy = lambda f: H**2*((3-k*k*eta**2)*f+3*eta*s.diff(f,eta))
    flux = lambda f: H**2*eta**2*(s.diff(f,eta)+f/eta)  # S_0i /(i k_i)
    tf = lambda f: -H**2*eta**2*f  # coefficient of k_i k_j - delta_ij k^2/3
    check('chi_y_exact_energy', energy(cy)+H**2*k*eta**3*s.cos(z))
    check('chi_j_exact_energy', energy(cj)-H**2*k*eta**3*s.sin(z))
    check('chi_y_exact_flux_coefficient', flux(cy)+H**2*eta**2/k*(s.cos(z)+z*s.sin(z)))
    check('chi_j_exact_flux_coefficient', flux(cj)-H**2*eta**2/k*(s.sin(z)-z*s.cos(z)))
    check('chi_y_series', s.series(cy,eta,0,5).removeO()+3/(k**3*eta)+eta/(2*k)+k*eta**3/8)
    check('chi_j_series', s.series(cj,eta,0,6).removeO()-k*k*eta**4/15)
    for name,f,ep,ec,qp,qc,apow,ac in [
        ('y',cy,3,-H**2*k,2,-H**2/k,1,3*H**2/k**3),
        ('j',cj,4,H**2*k*k,5,H**2*k*k/3,6,-H**2*k*k/15)]:
        check(f'{name}_energy_leading',s.limit(energy(f)/eta**ep,eta,0,dir='-')-ec)
        check(f'{name}_flux_leading',s.limit(flux(f)/eta**qp,eta,0,dir='-')-qc)
        check(f'{name}_TF_leading',s.limit(tf(f)/eta**apow,eta,0,dir='-')-ac)
    check('energy_cancels_eta_inverse', s.limit(eta*energy(cy),eta,0))
    check('energy_cancels_eta_linear', s.limit(energy(cy)/eta,eta,0))

    # Compute the tensor from the conformal flat metric and its Christoffels,
    # then compare all components with the displayed orthonormal formulas.
    xx,yy,zz=s.symbols('x y z', real=True)
    k1,k2,k3=s.symbols('k1 k2 k3',real=True)
    ks=(k1,k2,k3); ksq=sum(q*q for q in ks)
    co=(eta,xx,yy,zz)
    a=-1/(H*eta)
    g=a*a*s.diag(-1,1,1,1); gi=g.inv()
    ga=[[[s.simplify(s.Rational(1,2)*sum(gi[m,l]*(s.diff(g[l,jj],co[ii])+s.diff(g[l,ii],co[jj])-s.diff(g[ii,jj],co[l])) for l in range(4)))
         for jj in range(4)] for ii in range(4)] for m in range(4)]
    f=s.Function('f')(eta)
    phase=s.exp(s.I*(k1*xx+k2*yy+k3*zz))
    ff=f*phase
    hes=s.Matrix(4,4,lambda i,jj:s.diff(ff,co[i],co[jj])-sum(ga[m][i][jj]*s.diff(ff,co[m]) for m in range(4)))
    box=sum(gi[i,jj]*hes[i,jj] for i in range(4) for jj in range(4))
    PP=hes-g*(box+3*H*H*ff)
    PPhat=PP/(a*a*phase)
    on_shell={s.diff(f,eta,2):2*s.diff(f,eta)/eta-(ksq-4/eta**2)*f}
    U=H*H*((3-ksq*eta*eta)*f+3*eta*s.diff(f,eta))
    target=s.zeros(4)
    target[0,0]=U
    for i in range(3):
        target[0,i+1]=target[i+1,0]=s.I*ks[i]*H*H*eta*eta*(s.diff(f,eta)+f/eta)
        for jj in range(3):
            dij=int(i==jj)
            target[i+1,jj+1]=dij*U/3-H*H*eta*eta*f*(ks[i]*ks[jj]-dij*ksq/3)
    for i in range(4):
        for jj in range(4):
            check(f'orthonormal_component_{i}{jj}_from_metric', (PPhat[i,jj]-target[i,jj]).subs(on_shell))
    check('free_orthonormal_trace',-target[0,0]+sum(target[i,i] for i in range(1,4)))

    result={'scope':'Appendix E real-basis and orthonormal-stress algebra; all limits at fixed k>0',
            'sympy_version':s.__version__,'checks':checks,'passed':sum(checks.values()),'total':len(checks)}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if not all(checks.values()):
        raise SystemExit('Appendix E algebra verification failed.')

if __name__=='__main__':
    main()
