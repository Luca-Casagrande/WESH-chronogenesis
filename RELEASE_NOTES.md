# Version 2.0.0

Align the reproducibility supplement with `main_40_armonizzato.tex`.
The manuscript source is not included; its SHA-256 is recorded in
`source_manifest.json`.

- Replace six verification programs with the supplied revised scripts. The
  other fourteen verification programs and the documented BBN inputs are unchanged.
- Evaluate the entanglement gate as one half the trace norm of the law average
  of conditional pair correlations weighted by capped negativity. Product-state
  realizations have zero contribution; the gate is not a function of the
  barycenter alone.
- Use zero vertex response weights. Vertex events still calibrate time; the
  curvature response is supplied by edge innovations. Include finite checks
  of vertex neutrality and Proposition 4.5.
- Verify `sigma_eff^2 = (r_bi * w_bi^2 / 2) * q * (1-q)`, with no local floor,
  and the aligned active stationary sector identities of Proposition 5.3.
- Update the physical embedding checks and Appendix C examples for the
  law-dependent gate, finite generator/coupling constants, frozen-flow
  approximation and initial-law stability in `W_1`. The short coupling
  exposure `0.1` is intentional.
- Refresh the manuscript-to-code map, execution records, source provenance and
  checksums after running all twenty programs with the tested dependencies.

Verification: **20/20 programs passed** with Python 3.12.10 and
`requirements-tested.txt`, using `python run_all.py --jobs 1`.
The nominal 95% one-history endpoint remains
`sigma_95 = 0.000111694662524`. Its amplitude-bound consequences reproduce the
manuscript's rounded values `1.0e-7`, `2.5e-8` and `1.25e-8`.

These are consistency checks of the declared finite models and numerical
examples; the general analytic proofs remain in the manuscript.

# Version 1.0.1

Maintenance refinement of the same reproducibility supplement.

- Correct the Section 5 command in the verification script's docstring.
- Remove the redundant random self-subtraction check in Section 9; retain the
  exact algebraic counterexample as a comment. Section 9 now has 68 checks.
- Rerun all twenty verification programs and refresh the execution record,
  tested dependencies, source manifest and file checksums.
- Include Git line-ending and ignore rules in the published repository.

No model, input value or numerical tolerance is changed. The reproduced nominal
95% one-history endpoint remains `sigma_95 = 0.000111694662524`.

# Version 1.0.0

Reproducibility supplement for the WESH chronogenesis manuscript.

- Twenty programs covering finite-model algebra, event statistics, geometric
  identities, cosmological integrals and the compressed BBN benchmark.
- Includes the current activity–purity identity, physical algebra, gate embedding,
  sector-weight identities, flow-bound examples, Bell relaxation and cross-epoch
  radiation memory.
- Documented input data, reference outputs, dependency versions and checksums.

This version records computational checks, not new observational experiments.
The manuscript-to-code map states which results are checked numerically and
which remain analytic proofs. The optional publisher-PDF check is outside the
default suite.
