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
