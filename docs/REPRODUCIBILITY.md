# Reproducibility

## Fast independent smoke test

Run:

`python scripts/reproducibility_smoke_test.py`

The supplied reference release passes all five checks:

| Check | Reference result | Status |
|---|---:|---|
| Kernel-tail local exponent at alpha=0.8 | 0.79985912 | PASS |
| Exact-periodization maximum absolute error | 4.44e-16 | PASS |
| FFT open operator vs direct summation | 1.88e-15 in the independent smoke test | PASS |
| Static open beta at alpha=0.8 | 0.80024977 | PASS |
| Frozen dynamic benchmark at alpha=0.8, rho=0.30, N=128, Z=10 | within 1 percent relative tolerance | PASS |

Machine-readable smoke-test output is stored in `logs/reproducibility_smoke_test.csv`.

## Additional bundled validation

- Exact periodization is also tabulated in `data/validation/operator_validation.csv`.
- Open-operator FFT/direct checks are included in `data/validation/operator_validation.csv` and `data/representative_raw/open_operator_fft_vs_direct_N64.csv`.
- Time-step convergence data are in `data/validation/time_step_convergence.csv`.
- Dynamic step sensitivity is in `data/validation/dynamic_step_sensitivity.csv`.
- Dense infinite-symbol plotting verification is in `logs/dispersion_symbol_verification.csv`.

## Fast rebuild

`python scripts/run_all.py`

This rebuild uses the bundled numerical results and regenerates summaries, figures, and panel-source CSVs.

## Full deterministic recomputation

`python scripts/run_all.py --recompute`

## Paired stochastic pilot

`python scripts/run_all.py --stochastic`

The stochastic calculation is deliberately limited to `R=64` paired realizations. It is secondary evidence only.
