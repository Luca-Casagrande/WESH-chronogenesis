# WESH chronogenesis
## Reproducibility supplement · Online Resource 1

**Luca Casagrande** · Independent researcher, Portogruaro, Italy  
ORCID: [0009-0007-6645-2316](https://orcid.org/0009-0007-6645-2316)

Computational companion to **A dissipative time-field completion of
Wheeler–DeWitt dynamics: chronogenesis and a conserved cosmological response**.

This package reproduces the finite event construction, geometric response,
causal-volume statistics and compressed primordial-abundance benchmark used in
the manuscript. It contains **20 executable verification programs**, documented
inputs and the outputs of a complete run. The checks include symbolic identities,
finite matrices, ODE integrations and seeded sampling of the declared statistical
reference models.

## Version 2.0.0

This version aligns six programs with `main_40_armonizzato.tex`; the other
fourteen programs and the BBN inputs are unchanged. The manuscript source is
identified by SHA-256 in `source_manifest.json` and is not redistributed.

- The entanglement gate averages conditional pair correlations, weighted by
  capped negativity, over the law and then takes one half the trace norm.
  It depends on the law, not just its barycenter.
- Vertex response weights are zero: vertex events calibrate time while edge
  innovations supply the curvature response (Proposition 4.5).
- The Bell realization has
  `sigma_eff^2 = (r_bi * w_bi^2 / 2) * q * (1-q)`, with no local floor.
  The aligned active stationary sector is checked through the finite identities
  of Proposition 5.3; its unique barycenter does not imply a unique phase law.
- The Appendix C examples use the updated finite constants and `W_1` estimates
  for frozen-flow approximation and initial-law stability. The short coupling
  exposure is intentionally `0.1`.

The BBN endpoint is unchanged. Its amplitude-bound consequences are approximately
`r_bi*w_bi^2 <= 1.0e-7` at `q=1/2`, `q*(1-q) <= 2.5e-8` at unit coupling,
and `(1-q)/2 <= 1.25e-8` near alignment at unit coupling.
See `RELEASE_NOTES.md` and `VERIFICATION_REPORT.md` for the supplied execution.

## Run the checks

Python 3.10 or newer is required; the supplied run uses Python 3.12.10. Create a
virtual environment, install the dependencies, then run the suite:

```sh
python -m pip install -r requirements.txt
python run_all.py
```

For the exact tested dependency versions, use `requirements-tested.txt` with
Python 3.12.10. All scientific checks run locally without network access or GPU/QPU
hardware. `run_all.py` works from any working directory. It verifies the explicit
script inventory, writes a log beside each program, regenerates the numerical
outputs and records execution status in `run_summary.json`. It exits nonzero on
missing, failed or timed-out checks. `tested_environment.json` records the actual
versions and platform. Re-running overwrites the supplied generated outputs.

On a freshly extracted copy, `python verify_integrity.py` checks the delivered
snapshot against `SHA256SUMS.txt`. Regenerating outputs changes execution times
and may change last-place floating-point digits, so verify integrity before
re-running the suite. See `VERIFICATION_REPORT.md` for the supplied run.

For example, to reproduce the BBN calculations alone:

```sh
python section09/verify_wesh_section9.py
```

## What the package contains

| Path | Content |
|---|---|
| `section02`–`section05` | Physical lifts, event law, time changes, matching and Bell realization |
| `section06`–`section09` | Conserved response, cosmology, inherited radiation and BBN calculations |
| `appendixE`, `appendixG`, `appendixEG` | Geometric modes, cosmological integrals and statistical checks |
| `extensions` | Current-manuscript algebra, flow bounds, activity/purity and relaxation checks |
| `inputs/BBN_inputs.json` | Observational values and benchmark assumptions |
| `MANUSCRIPT_MAP.md` | Program-to-equation correspondence and exact coverage |
| `run_summary.json`, `tested_environment.json` | Execution record and tested environment |
| `SHA256SUMS.txt`, `source_manifest.json` | File integrity and input-source provenance |

## Primordial-abundance benchmark

The documented source is T.-H. Yeh et al., *The LBT Yp Project. V. Cosmological
Implications of a New Determination of Primordial 4He*, **Astrophysical Journal
1008, 241 (2026)**, DOI [10.3847/1538-4357/ae879c](https://doi.org/10.3847/1538-4357/ae879c).
The adopted input is Table 2, **BBN+Yp+D**, together with the central abundance
scalings in Eqs. (8)–(9). The BBN+CMB row is not used.

For fifty e-folds of ideal radiation, zero initial response and thermal transfer
factor one, the one-history Gaussian-compression test gives

```text
sigma_95 = 1.1169466252e-4
B(50)    = 8383.37677099
```

This is a nominal 95% confidence endpoint in the stated sampling model. The
separate typical-history tolerance is approximately `1.82351e-5`; the suite
reports both quantities and the associated coverage and template diagnostics.

The source PDF is not redistributed. Numerical inputs are explicitly recorded
in this package. An optional source-text check can be run on a separately
obtained publisher PDF:

```sh
python -m pip install -r requirements-optional.txt
python section09/verify_wesh_section9_refinement.py --published-pdf PATH_TO_PDF
```

The default suite does not re-fetch or re-verify the publisher PDF.

## Scope and citation

These are computational checks of the specified models and calculations. The
analytic proofs are in the manuscript. Synthetic samples are labelled in the
programs; the package contains no new observational measurements, hardware
experiments or nuclear-network/CMB likelihood calculation.

Please cite the manuscript and the exact archived release used for reproduction.
`CITATION.cff` supplies authorship and package metadata.
The package version is `2.0.0`; its
creation does not imply publication or acceptance of the manuscript.
