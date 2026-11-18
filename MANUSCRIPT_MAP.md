# Manuscript-to-code map

The table identifies the calculations exercised by each program. Equation and
proposition labels refer to the aligned manuscript source; the labels remain
stable under renumbering. The input source hash is recorded in
`source_manifest.json`. The manuscript itself is not duplicated in this bundle.

| Program | Location | Checks | Manuscript labels |
|---|---|---|---|
| [`section02/verify_section2.py`](section02/verify_section2.py) | 2.3, 2.5 | Noether identity, spectral lift, compression, Kraus completeness | `eq:noether-dirichlet`, `eq:spectral-event-lift`, `eq:finite-physical-instrument` |
| [`section02/verify_lifts_and_primitives.py`](section02/verify_lifts_and_primitives.py) | 2.1, 2.4, 2.5 | Constraint algebra, primitive coordinates, energy exchange and kernel normalization | `eq:additive-constraint-algebra`, `eq:primitive-rescaling`, `eq:lift-sector-exchange` |
| [`section03/verify_section3.py`](section03/verify_section3.py) | 2.2; 3.1–3.6 | Law-averaged negativity-weighted gate, product-law neutrality, Lipschitz bounds, centered events, covariance, stationarity, spectators and instantaneous locality | `eq:trace-norm-wesh-gate`, `eq:section2-product-faithfulness`, `eq:gate-trace-lipschitz`, `eq:conditional-drift-identity`, `eq:marginal-defined-extension` |
| [`section03/verify_event_record.py`](section03/verify_event_record.py) | 3.6; 5.3 | Conditional event records at fixed auxiliary exposure, law-dependent gate and edge-only weighted activity | `eq:event-record-preparation-dependence`, `eq:trace-norm-wesh-gate`, `eq:stationary-effective-amplitude` |
| [`section04/verify_matching.py`](section04/verify_matching.py) | 4.1–4.5 | Edge-only weighted activity, regional variance, bilocal activity bound, vertex neutrality and rate cancellation | `eq:regional-weighted-activity`, `eq:minimum-activity-target`, `prop:vertex-neutrality-equivalence`, `eq:local-global-calibration-bilocal` |
| [`section04/verify_time_change.py`](section04/verify_time_change.py) | 4.1; D.1 | Generalized inverse with plateaus and terminal activity, finite calculus examples | `eq:chronogenetic-channel-proportions`, `eq:app-time-changed-compensators` |
| [`section05/verify_section5.py`](section05/verify_section5.py) | 5; F | Bell negativity and gate, generator, stationary law, neutral vertices and amplitude without a local floor | `eq:trace-norm-wesh-gate`, `eq:bell-variances-intensities`, `eq:stationary-effective-amplitude`, `eq:stationary-local-floor`, `eq:app-bell-transition-generator` |
| [`section05/verify_additions.py`](section05/verify_additions.py) | 5.3–5.4; 9.4; F | Conditional/ensemble regional variance, aligned active stationary sector identities and arithmetic of the amplitude bound | `eq:stationary-block-volume-law`, `eq:app-conditional-block-variance`, `prop:aligned-stationary-uniqueness`, `eq:bbn-response-weight-region` |
| [`section06/verify_wesh_section6.py`](section06/verify_wesh_section6.py) | 6; E | Metric-derived curvature variation, conserved response, kernel, modes and optics | `eq:conserved-conformal-identities`, `eq:expansion-power-transfer`, `eq:conformal-endpoint-observables` |
| [`section06/verify_section6_refinements.py`](section06/verify_section6_refinements.py) | 6.3–6.5 | Compensating fluid and Ma–Bertschinger potential dictionary | `eq:de-sitter-trace-filter`, `eq:test-matter-conservation` |
| [`section07/verify_wesh_section7.py`](section07/verify_wesh_section7.py) | 7 | FRW closure, Friedmann propagation, Ahmed–Sorkin correspondence and numerical ODE comparison | `eq:exact-frw-response-system`, `eq:frw-retarded-response`, `eq:ahmed-sorkin-dictionary` |
| [`section07/verify_section7_refinements.py`](section07/verify_section7_refinements.py) | 7.2–7.3 | Trace-case marginality and comparison at equal scale factor | `subsec:ahmed-sorkin-identity`, `eq:reference-fraction-response` |
| [`section08/verify_wesh_section8.py`](section08/verify_wesh_section8.py) | 8 | Cone volume, covariance, integrated variance, inherited projection and Gaussian-reference sampling | `eq:radiation-past-cone-volume`, `eq:radiation-unequal-time-covariance`, `eq:independent-inherited-variance` |
| [`section08/verify_section8_refinements.py`](section08/verify_section8_refinements.py) | 8–9 | Background terms, duration/thermal sensitivity and BBN normalization handoff | `eq:inherited-thermal-transfer`, `eq:bbn-inherited-dispersion-scale` |
| [`section09/verify_wesh_section9.py`](section09/verify_wesh_section9.py) | 9 | One-history inversion, benchmark, sensitivity and coverage diagnostics | `eq:bbn-confidence-boundaries`, `eq:bbn-numerical-benchmark`, `eq:bbn-response-weight-region` |
| [`section09/verify_wesh_section9_refinement.py`](section09/verify_wesh_section9_refinement.py) | 9.1, 9.5; G.4 | Abundance contrast, helium curvature, optional published-PDF provenance | `eq:app-bbn-abundance-contrast`, `eq:app-helium-template-expansion` |
| [`appendixE/verify_wesh_appendixE.py`](appendixE/verify_wesh_appendixE.py) | E | Real mode basis, metric components and fixed-k asymptotics | `eq:app-real-free-mode-basis`, `eq:app-free-mode-stress-components`, `eq:app-real-mode-stress-asymptotics` |
| [`appendixG/verify_wesh_appendixG.py`](appendixG/verify_wesh_appendixG.py) | G | Residual equation, truncated cone, thermal transport, confidence and likelihood-ratio diagnostics | `eq:app-general-mixed-equation-residual`, `eq:app-truncated-radiation-cone`, `eq:app-bbn-likelihood-ratio-coverage` |
| [`appendixEG/verify_review_refinements.py`](appendixEG/verify_review_refinements.py) | E–G | Ambient embedding, Gauss sign, window remainder and likelihood monotonicity | `prop:app-conformal-kernel`, `eq:app-window-forcing-covariance`, `eq:app-bbn-likelihood-ratio-coverage` |
| [`extensions/verify_current_results.py`](extensions/verify_current_results.py) | 2.5, 3.5, 8.5; A–C, F | Physical algebra and negativity-preserving embeddings, sector weights, finite generator/coupling constants, W_1 flow examples, activity/purity, Bell relaxation and inherited amplitude | `eq:dirac-reduced-map`, `eq:app-isometric-gate`, `eq:app-sector-weight-dynamics`, `eq:app-sector-quadratic-budget`, `eq:app-generator-flow-bounds`, `eq:app-finite-flow-constants`, `eq:app-frozen-flow-error`, `eq:app-prescribed-product-error`, `eq:app-gate-controlled-stability`, `eq:activity-purity-decomposition`, `prop:app-bell-relaxation`, `eq:inherited-epoch-memory` |

## Coverage of the appendices

A and F have explicit matrix and finite-distribution checks, including
negativity-preserving embeddings and the Bell amplitude without a local floor. B checks the local
sector-weight drift and quadratic-variation identities; bounded-martingale
convergence remains an analytic result. C exercises the stated bounds on finite
examples in `W_1`, including a short (`0.1`) self-consistent Bell-law flow,
a discontinuous prescribed path and a finite spectator. The stable label
`eq:stationary-local-floor` now denotes the amplitude range with lower bound
zero and upper bound `r_bi*w_bi^2/8`; it does not assert a positive local floor.
D checks explicit time-change examples and covariance bookkeeping; its general
optional-sampling and stochastic-Fubini proofs remain in the manuscript.
E and G combine symbolic calculations with numerical diagnostics.

The additional `extensions` program checks formulas already in the manuscript.
It introduces no new physical assumptions or observations. Passed numerical
cases are regression checks, not general proofs. The source-to-geometry matching,
continuum limit, nuclear network and CMB likelihood are not computed by this suite.
