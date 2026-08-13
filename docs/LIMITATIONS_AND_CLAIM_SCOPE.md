# Limitations and Claim Scope

This file is a guardrail for manuscript editing, peer review, and repository reuse. The scientific package is strongest when its claims remain aligned with the evidence actually bundled here.

## Claims supported as primary results

1. **Open finite-domain closure error is algebraic and tracks the interaction exponent.** Across the six deterministic `(alpha,rho)` cases, the fitted open-boundary exponent is consistent with `beta_open approximately alpha` over the resolved asymptotic range.
2. **The deterministic finite-domain error grows on the modulational-instability timescale in the early MI regime.** The error-growth rates are tied to the independently computed MI spectrum; the package reports the measured ratios and their spread rather than asserting exact equality.
3. **The early deterministic fidelity horizon grows approximately logarithmically with domain size.** The tested leading relation is `Z_f approximately (beta/g) ln N + C`, with the small-wave-number `(alpha,rho)=(0.8,0.15)` case retained as a pre-asymptotic stress case.
4. **The background-tail correction improves finite-domain fidelity.** For fixed nonzero modes it substantially accelerates static convergence and delays deterministic error growth.

## Claims supported only as secondary/pilot evidence

5. **Paired noisy-background distributional discrepancy is reduced by the correction in the bundled pilot.** This is based on `R=64` paired realizations for `(alpha,rho)=(0.8,0.30)` and should remain a secondary observation.

## Claims this package does NOT support

- It does not establish converged `Q_0.99`, `Q_0.999`, block-maxima, or rare/extreme-event statistics.
- It does not establish a universal fidelity law for every nonlocal nonlinear system.
- It does not prove the nonlinear horizon relation globally in late nonlinear dynamics; the derivation/prediction is an early MI-regime scaling supported by deterministic numerics.
- It does not show that the corrected `1+alpha` exponent is uniform as `q -> 0`.
- It does not provide direct experimental validation or a device-specific photonic parameter mapping.
- It does not establish statistical-horizon prediction from a production Monte Carlo campaign.

## Required wording discipline

Prefer:

- `consistent with`, `tracks`, `numerically supports`, `early-regime scaling`, `paired stochastic pilot`.

Avoid unless new evidence is added:

- `universal`, `proves` for simulation-only statements, `rogue-wave statistics`, `converged extreme-event statistics`, `statistical fidelity horizon` as a headline result.

## If new calculations are added

Do not overwrite the frozen files silently. Create a new versioned result set and update `manifests/claims_and_evidence.csv`, `data/headline_results.csv`, figure-source CSVs, and checksums together.
