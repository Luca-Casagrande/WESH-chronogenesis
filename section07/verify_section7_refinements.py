#!/usr/bin/env python3
"""Checks of structural marginality and the fixed-scale-factor comparison.

Scope: algebra within the specified homogeneous effective closure.
These checks neither validate the stochastic matching nor establish global
nonlinear stability. Run with Python 3.10+ and SymPy.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sympy as sp


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).with_name('section7_refinements_checks.json'))
    args = parser.parse_args()
    checks: dict[str, bool] = {}

    def zero(name: str, expr: sp.Expr) -> None:
        checks[name] = sp.simplify(expr) == 0

    n, n0 = sp.symbols('n n0', real=True)
    mu, nu, rho0 = sp.symbols('mu nu rho0', positive=True)
    D0, w = sp.symbols('D0 w', real=True)
    b = 1 + nu / mu
    Dfree = D0 * sp.exp(-3*b*(n-n0))
    radiation = rho0 * sp.exp(-4*(n-n0))
    frac = Dfree / radiation
    zero('constant_source_residual_solution', sp.diff(Dfree, n) + 3*b*Dfree)
    zero('radiation_relative_exponent', sp.diff(frac,n)-(4-3*b)*frac)
    zero('non_growth_condition_numerator', (3*b-4)*mu-(3*nu-mu))
    zero('trace_choice_density_power', (3*b).subs(mu,3*nu)-4)
    zero('trace_choice_fraction_constant', sp.diff(frac,n).subs(mu,3*nu))
    checks['interior_relative_decay'] = bool((4-3*b).subs(mu,nu) < 0)
    checks['beyond_trace_relative_growth'] = bool((4-3*b).subs(mu,4*nu) > 0)

    # At equal n and equal conserved material normalization, density is equal
    # even when the two proper-time histories differ.
    rho_a = rho0*sp.exp(-3*(1+w)*n)
    rho_b = rho0*sp.exp(-3*(1+w)*n)
    zero('matched_material_density_at_fixed_n',rho_a-rho_b)
    kappa, U, Href = sp.symbols('kappa U Href', positive=True)
    rho_ref = 3*Href**2/kappa
    zero('fixed_n_expansion_fraction', U/(6*Href**2) - (U/(kappa*rho_ref))/2)

    result = {'number_of_checks':len(checks), 'checks':checks,
              'all_passed':all(checks.values()),
              'scope':'Homogeneous residual and comparison identities; no global stability or empirical inference.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if not result['all_passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
