# Verification record

Package: **WESH chronogenesis 2.0.0**

Completed (UTC): `2026-11-18T00:38:46.581971+00:00`

Result: **20 of 20 verification programs completed successfully.**

Command:

```sh
python run_all.py --jobs 1
```

The run used a new Python 3.12.10 virtual environment installed from
`requirements-tested.txt`: NumPy 2.4.6, SciPy 1.18.0,
SymPy 1.14.0 and mpmath 1.3.0.
Platform: `Windows-11-10.0.26200-SP0`. All twenty programs exited with status zero.

The six revised programs are aligned with `main_40_armonizzato.tex`, identified
by SHA-256 in `source_manifest.json`. Their code is the supplied update, with
the repository's LF line-ending convention. The other fourteen verification
programs and `inputs/BBN_inputs.json` are byte-for-byte unchanged from v1.0.1.
The program-level logs and generated JSON/CSV files were regenerated.
`run_summary.json` records this execution; `tested_environment.json` records
the environment, whose versions match the previous supplied run.

Section 5 reports **109 exact symbolic checks** in the base program and
**48 additional checks** in the additions program (which also imports the
base program). The extensions program reports **241 checks**; Section 9
retains **68 checks**. These counts describe different kinds of checks and
are not independent tests of general theorems.

| Program | Result | Seconds |
|---|---|---:|
| `appendixE/verify_wesh_appendixE.py` | PASS | 5.782 |
| `appendixEG/verify_review_refinements.py` | PASS | 0.860 |
| `appendixG/verify_wesh_appendixG.py` | PASS | 4.968 |
| `extensions/verify_current_results.py` | PASS | 1.969 |
| `section02/verify_lifts_and_primitives.py` | PASS | 0.844 |
| `section02/verify_section2.py` | PASS | 0.250 |
| `section03/verify_event_record.py` | PASS | 1.093 |
| `section03/verify_section3.py` | PASS | 1.750 |
| `section04/verify_matching.py` | PASS | 0.250 |
| `section04/verify_time_change.py` | PASS | 0.250 |
| `section05/verify_additions.py` | PASS | 2.469 |
| `section05/verify_section5.py` | PASS | 1.125 |
| `section06/verify_section6_refinements.py` | PASS | 0.719 |
| `section06/verify_wesh_section6.py` | PASS | 4.875 |
| `section07/verify_section7_refinements.py` | PASS | 0.734 |
| `section07/verify_wesh_section7.py` | PASS | 1.812 |
| `section08/verify_section8_refinements.py` | PASS | 2.969 |
| `section08/verify_wesh_section8.py` | PASS | 1.875 |
| `section09/verify_wesh_section9.py` | PASS | 1.828 |
| `section09/verify_wesh_section9_refinement.py` | PASS | 1.500 |

## Reproduced BBN benchmark

| Quantity | Value |
|---|---:|
| B(50) | 8383.376770994 |
| Nominal 95% one-history endpoint | 0.000111694662524 |
| Inherited species-equivalent rms at endpoint | 0.936378439244 |
| Inherited radiation-fraction rms at endpoint | 0.152433699412 |
| Linear speed-up rms at endpoint | 0.0762168497059 |

The Section 9 result file, including its benchmark and sampling diagnostics,
is unchanged from v1.0.1. The six-script update uses this endpoint in
`sigma_eff^2 = (r_bi*w_bi^2/2)*q*(1-q)`:

| Amplitude-bound consequence | Reproduced value | Manuscript rounding |
|---|---:|---:|
| Coupling bound at q=1/2 | 9.980558108996e-08 | 1.0e-7 |
| q(1-q) bound at unit coupling | 2.495139527249e-08 | 2.5e-8 |
| Misaligned occupation bound near alignment, unit coupling | 1.247569763624e-08 | 1.25e-8 |

The last row is the near-alignment approximation stated in the manuscript;
the additions program also checks the exact root.

## Alignment and scope

The updated checks use the law average of negativity-weighted conditional pair
correlations, zero vertex response weights, no state-independent local floor,
and the finite identities for Propositions 4.5 and 5.3. Appendix C exercises
the updated finite constants and `W_1` flow estimates; the coupling exposure
`0.1` is intentional. `MANUSCRIPT_MAP.md` lists the corresponding stable labels.

The result refers to the input values and assumptions in `inputs/BBN_inputs.json`.
The optional publisher-PDF provenance check was not run in this execution.
The programs check consistency of the declared finite models and selected
numerical examples. The general analytic proofs remain in the manuscript.
