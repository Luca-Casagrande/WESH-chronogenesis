# Verification record

Package: **WESH chronogenesis 1.0.1**  
Completed (UTC): `2026-10-29T11:39:24.533432+00:00`  
Result: **20 of 20 verification programs completed successfully.**

Command:

```sh
python run_all.py --jobs 1
```

The program-level logs and generated JSON/CSV files are included next to their
programs. `run_summary.json` is the machine-readable execution record.
Section 9 reports **68 checks** after removing the redundant random
self-subtraction check. Its algebraic counterexample remains a comment.

| Program | Result | Seconds |
|---|---|---:|
| `appendixE/verify_wesh_appendixE.py` | PASS | 2.469 |
| `appendixEG/verify_review_refinements.py` | PASS | 0.750 |
| `appendixG/verify_wesh_appendixG.py` | PASS | 2.875 |
| `extensions/verify_current_results.py` | PASS | 1.891 |
| `section02/verify_lifts_and_primitives.py` | PASS | 0.656 |
| `section02/verify_section2.py` | PASS | 0.281 |
| `section03/verify_event_record.py` | PASS | 0.985 |
| `section03/verify_section3.py` | PASS | 1.563 |
| `section04/verify_matching.py` | PASS | 0.266 |
| `section04/verify_time_change.py` | PASS | 0.203 |
| `section05/verify_additions.py` | PASS | 1.515 |
| `section05/verify_section5.py` | PASS | 1.110 |
| `section06/verify_section6_refinements.py` | PASS | 0.594 |
| `section06/verify_wesh_section6.py` | PASS | 3.984 |
| `section07/verify_section7_refinements.py` | PASS | 0.781 |
| `section07/verify_wesh_section7.py` | PASS | 2.235 |
| `section08/verify_section8_refinements.py` | PASS | 2.782 |
| `section08/verify_wesh_section8.py` | PASS | 2.000 |
| `section09/verify_wesh_section9.py` | PASS | 1.875 |
| `section09/verify_wesh_section9_refinement.py` | PASS | 1.609 |

## Reproduced BBN benchmark

| Quantity | Value |
|---|---:|
| B(50) | 8383.376770994 |
| Nominal 95% one-history endpoint | 0.000111694662524 |
| Inherited species-equivalent rms at endpoint | 0.936378439244 |
| Inherited radiation-fraction rms at endpoint | 0.152433699412 |
| Linear speed-up rms at endpoint | 0.0762168497059 |

The result refers to the input values and assumptions in `inputs/BBN_inputs.json`.
The optional publisher-PDF provenance check was not run in this execution.
The numerical programs exercise the named identities and finite examples listed
in `MANUSCRIPT_MAP.md`; they do not certify the full manuscript or constitute a
new observational test.
