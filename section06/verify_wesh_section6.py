#!/usr/bin/env python3
"""Reproducible checks of WESH Section 6 on fixed four-dimensional de Sitter.

The connection variation is obtained from the perturbed metric, not supplied
as the claimed Ricci identity. Further checks cover the retarded homogeneous
filter, kernel, finite-k free modes, conformal optics and test dust.

Requirements: Python 3.10+, SymPy and NumPy. No network access or package
installation. Run: python verify_wesh_section6.py [--output checks.json]
The script exits nonzero if any check fails. It does not test the physical
matching, quantum-gravity emergence, stochastic nonlinear stability or data.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np
import sympy as sp


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("section6_checks.json"))
    args = parser.parse_args()
    checks: dict[str, bool] = {}
    start = perf_counter()

    def check(name: str, expression) -> None:
        value = bool(expression) if isinstance(expression, (bool, np.bool_)) else sp.simplify(expression) == 0
        checks[name] = bool(value)

    def check_matrix(name: str, matrix: sp.Matrix) -> None:
        checks[name] = all(sp.simplify(v) == 0 for v in matrix)

    t, x, y, z = sp.symbols("t x y z", real=True)
    H = sp.symbols("H", positive=True)
    coords = (t, x, y, z)
    a = sp.exp(H*t)
    g = sp.diag(-1, a*a, a*a, a*a)
    gi = sp.diag(-1, 1/(a*a), 1/(a*a), 1/(a*a))
    chi = sp.Function("chi")(*coords)
    gamma = [[[sp.S.Zero for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for i in range(1, 4):
        gamma[0][i][i] = H*a*a
        gamma[i][0][i] = gamma[i][i][0] = H

    def hessian(f):
        return sp.Matrix(4, 4, lambda i,j: sp.diff(f,coords[i],coords[j])
             - sum(gamma[k][i][j]*sp.diff(f,coords[k]) for k in range(4)))

    def box(f):
        hes = hessian(f)
        return sp.simplify(sum(gi[i,i]*hes[i,i] for i in range(4)))

    def scalar_tensor(f):
        return hessian(f)-g*(box(f)+3*H*H*f)

    def divergence(tensor):
        return sp.Matrix([sp.simplify(sum(gi[i,i]*(
            sp.diff(tensor[i,j],coords[i])
            -sum(gamma[k][i][i]*tensor[k,j]+gamma[k][i][j]*tensor[i,k]
                 for k in range(4))) for i in range(4))) for j in range(4)])

    P = scalar_tensor(chi)
    trace = sum(gi[i,i]*P[i,i] for i in range(4))
    check_matrix("P_divergence_four_components", divergence(P))
    check("P_trace_off_shell", trace+3*(box(chi)+4*H*H*chi))

    # Direct derivative of Gamma(g_epsilon) at epsilon=0.
    eps = sp.symbols("eps", real=True)
    ge = (1-eps*chi)*g
    gie = gi/(1-eps*chi)
    dgamma = [[[sp.simplify(sp.diff(sp.Rational(1,2)*sum(
        gie[k,l]*(sp.diff(ge[l,j],coords[i])+sp.diff(ge[l,i],coords[j])
                  -sp.diff(ge[i,j],coords[l])) for l in range(4)), eps).subs(eps,0))
        for j in range(4)] for i in range(4)] for k in range(4)]
    dRic = sp.zeros(4)
    for i in range(4):
        for j in range(4):
            dRic[i,j] = sp.simplify(
                sum(sp.diff(dgamma[k][i][j],coords[k])
                    -sp.diff(dgamma[k][i][k],coords[j]) for k in range(4))
                +sum(dgamma[k][i][j]*gamma[l][k][l]
                     +gamma[k][i][j]*dgamma[l][k][l]
                     -dgamma[l][i][k]*gamma[k][j][l]
                     -gamma[l][i][k]*dgamma[k][j][l]
                     for k in range(4) for l in range(4)))
    check_matrix("Ricci_variation_from_metric_all_components",
                 dRic-hessian(chi)-g*box(chi)/2)
    dR = sp.simplify(sum(gi[i,i]*dRic[i,i] for i in range(4))+12*H*H*chi)
    check("scalar_R_including_inverse_metric_variation", dR-3*box(chi)-12*H*H*chi)
    dEin = dRic+6*H*H*chi*g-g*dR/2-3*H*H*chi*g
    check_matrix("Einstein_plus_Lambda_h_all_components", dEin-P)
    source = sp.Rational(3,4)*(box(chi)+4*H*H*chi)
    C = hessian(chi)-g*box(chi)/4
    check_matrix("conserved_source_decomposition", P+source*g-C)
    check("compensating_tensor_trace", sum(gi[i,i]*C[i,i] for i in range(4)))
    check_matrix("compensating_tensor_divergence", divergence(C)
                 -sp.Matrix([sp.diff(source,c) for c in coords]))
    check_matrix("gauge_invariant_Ricci_deviation", dRic+3*H*H*chi*g-P-2*source*g)
    check("curvature_source_identity", dR-4*source)

    f = sp.Function("f")(t)
    forcing = sp.Function("zeta")(t)
    U = 3*H*(H*f-sp.diff(f,t))
    V = sp.diff(f,t,2)+2*H*sp.diff(f,t)-3*H*H*f
    pf = scalar_tensor(f)
    check("homogeneous_U_from_tensor", pf[0,0]-U)
    check("homogeneous_V_from_tensor", pf[1,1]/(a*a)-V)
    check("homogeneous_Bianchi_off_shell", sp.diff(U,t)+3*H*(U+V))
    field_rule = {sp.diff(f,t,2): -3*H*sp.diff(f,t)+4*H*H*f-sp.Rational(4,3)*forcing}
    check("trace_filter_on_shell", (sp.diff(U,t)+4*H*U-4*H*forcing).subs(field_rule))
    dh = (H*f-sp.diff(f,t))/2
    check("proper_expansion_Friedmann_relation", U-6*H*dh)
    check("proper_expansion_filter", (sp.diff(dh,t)+4*H*dh-sp.Rational(2,3)*forcing).subs(field_rule))
    Om = sp.sqrt(1-eps*f)
    Hp = sp.diff(Om*a,t)/(Om*Om*a)
    check("expansion_from_perturbed_proper_time", sp.diff(Hp,eps).subs(eps,0)-dh)
    t0, h0, u = sp.symbols("t0 h0 u", real=True)
    response = sp.exp(-4*H*(t-t0))*h0+sp.Rational(2,3)*sp.Integral(
        sp.exp(-4*H*(t-u))*sp.Function('zeta')(u),(u,t0,t))
    check("retarded_expansion_solution", sp.diff(response,t)+4*H*response-sp.Rational(2,3)*forcing)
    check("retarded_expansion_initial_value", response.subs(t,t0).doit()-h0)
    freq = sp.symbols("frequency", real=True)
    transfer = sp.Rational(2,3)/(4*H-sp.I*freq)
    check("expansion_spectral_transfer", transfer*sp.conjugate(transfer)
          -4/(9*(freq*freq+16*H*H)))

    basis = [a, a*x, a*y, a*z, H*a*(x*x+y*y+z*z)-sp.exp(-H*t)/H]
    for j,b in enumerate(basis):
        check_matrix(f"kernel_basis_{j}_P_zero", scalar_tensor(b))
        check_matrix(f"kernel_basis_{j}_Hessian_identity", hessian(b)+H*H*b*g)
        xi = sp.Matrix([gi[i,i]*sp.diff(b,coords[i])/(2*H*H) for i in range(4)])
        lie = sp.Matrix(4,4,lambda i,j: sum(xi[k]*sp.diff(g[i,j],coords[k])
               +g[k,j]*sp.diff(xi[k],coords[i])+g[i,k]*sp.diff(xi[k],coords[j])
               for k in range(4)))
        check_matrix(f"kernel_basis_{j}_Lie_derivative", lie+b*g)
    jet = sp.Matrix(5,5,lambda i,j: basis[j] if i==0 else sp.diff(basis[j],coords[i-1]))
    check("kernel_five_independent_jets", jet.subs({t:0,x:0,y:0,z:0}).det()!=0)
    decay = sp.exp(-4*H*t)
    check("decaying_mode_U", U.subs(f,decay).doit()-15*H*H*decay)
    check("decaying_mode_V", V.subs(f,decay).doit()-5*H*H*decay)
    check("decaying_mode_deltaH", dh.subs(f,decay).doit()-sp.Rational(5,2)*H*decay)

    eta = sp.symbols("eta", real=True, nonzero=True)
    k = sp.symbols("k", positive=True)
    F = -eta*sp.exp(-sp.I*k*eta)*(1-3*sp.I/(k*eta)-3/(k*k*eta*eta))
    Fb = sp.conjugate(F)
    mode_eq = lambda ff: sp.diff(ff,eta,2)-2/eta*sp.diff(ff,eta)+(k*k-4/(eta*eta))*ff
    check("Fourier_mode_F_equation", mode_eq(F))
    check("Fourier_mode_F_conjugate_equation", mode_eq(Fb))
    check("Fourier_mode_Wronskian", F*sp.diff(Fb,eta)-sp.diff(F,eta)*Fb-2*sp.I*k*eta*eta)
    aa = -1/(H*eta)
    fdot = -H*eta*sp.diff(F,eta)
    fU = 3*H*(H*F-fdot)-k*k*F/(aa*aa)
    f01 = sp.I*k/aa*(fdot-H*F)
    check("Fourier_exact_U", fU-H*H*k*k*eta**3*sp.exp(-sp.I*k*eta))
    check("Fourier_exact_flux", f01-H*H*k*eta**2*(sp.I-k*eta)*sp.exp(-sp.I*k*eta))
    expansion = 3/(k*k*eta)+eta/2+k*k*eta**3/8-sp.I*k**3*eta**4/15
    check("Fourier_future_series", sp.series(F,eta,0,5).removeO()-expansion)
    check("Fourier_U_leading_order", sp.limit(fU/eta**3,eta,0)-H*H*k*k)
    check("Fourier_flux_leading_order", sp.limit(f01/eta**2,eta,0)-sp.I*H*H*k)
    # Spatial TF component along k: S_TF11=-(2/3) k^2 chi/a^2.
    tf11 = -sp.Rational(2,3)*k*k*F/(aa*aa)
    check("Fourier_anisotropic_leading_order", sp.limit(tf11/eta,eta,0)+2*H*H)

    # Check that the component formulas follow from P, independently of the exact F.
    plane = f*sp.exp(sp.I*k*x)
    pt = scalar_tensor(plane)
    check("Fourier_U_component_from_P", (pt[0,0]/sp.exp(sp.I*k*x)
         -3*H*(H*f-sp.diff(f,t))+k*k*f/a**2))
    check("Fourier_flux_component_from_P", pt[0,1]/(a*sp.exp(sp.I*k*x))
         -sp.I*k/a*(sp.diff(f,t)-H*f))
    spatial_trace = sum(pt[i,i] for i in range(1,4))/a**2
    check("Fourier_TF_component_from_P", (pt[1,1]/a**2-spatial_trace/3)
          /sp.exp(sp.I*k*x)+sp.Rational(2,3)*k*k*f/a**2)

    em, ob, OmE, OmO, zz, DA = sp.symbols('chi_em chi_obs Om_em Om_obs zbar DAbar', real=True)
    subOm = {OmE: 1-eps*em/2, OmO: 1-eps*ob/2}
    zexpr = OmO/OmE*(1+zz)-1
    check("redshift_endpoint_variation", sp.diff(zexpr.subs(subOm),eps).subs(eps,0)/(1+zz)-(em-ob)/2)
    check("area_distance_endpoint_variation", sp.diff((OmE*DA).subs(subOm),eps).subs(eps,0)/DA+em/2)
    dl_ratio = OmO**2/OmE
    check("luminosity_distance_endpoint_variation", sp.diff(dl_ratio.subs(subOm),eps).subs(eps,0)-(em/2-ob))
    check("distance_duality_conformal_transport", (1+zexpr)**2*OmE*DA-dl_ratio*(1+zz)**2*DA)
    cp = sp.Function('chi')(eta)
    Phi, Psi = -cp/2, cp/2
    check("integrated_optical_potential", Phi+Psi)
    confH = sp.Function('mathcalH')(eta)
    dm, theta = sp.Function('delta_m')(eta), sp.Function('theta_m')(eta)
    deriv_dm2 = -(-confH*theta+k*k*Phi)+3*sp.diff(Psi,eta,2)
    source_dm = -k*k*Phi+3*sp.diff(Psi,eta,2)+3*confH*sp.diff(Psi,eta)
    check("dust_second_order_equation", deriv_dm2+confH*(-theta+3*sp.diff(Psi,eta))-source_dm)
    check("dust_conformal_forcing", source_dm-(k*k*cp/2+3*sp.diff(cp,eta,2)/2+3*confH*sp.diff(cp,eta)/2))

    # A finite-dimensional diagnostic of linear covariance pushforward.
    rng = np.random.default_rng(602026)
    min_eval = 1.0
    for _ in range(50):
        A = rng.normal(size=(8,8))
        K = A@A.T
        L = rng.normal(size=(5,8))
        KL = L@K@L.T
        min_eval = min(min_eval,float(np.linalg.eigvalsh((KL+KL.T)/2).min()))
    check("finite_covariance_positive_pushforward_50_cases", min_eval > -1e-10)

    result = {
        "scope": "fixed-background, linear effective scalar response; algebraic checks only",
        "sympy_version": sp.__version__, "numpy_version": np.__version__,
        "total_groups": len(checks), "passed_groups": sum(checks.values()),
        "checks": checks, "minimum_sample_covariance_eigenvalue": min_eval,
        "not_established_by_checks": ["microscopic gravitational matching", "spacetime emergence",
          "Gaussianity of the source", "uniform stochastic/nonlinear stability", "observational likelihood"],
        "seconds": round(perf_counter()-start,3),
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
    if not all(checks.values()):
        raise SystemExit('One or more verification groups failed.')


if __name__ == '__main__':
    main()
