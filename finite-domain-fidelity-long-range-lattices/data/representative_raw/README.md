# Representative run-level data

These compact files make the processing chain auditable without distributing every intermediate field snapshot.

- `deterministic_a080_r030_N128_timeseries.csv` - full saved deterministic diagnostic time series for the discovery case and `N=128`, including open/corrected field errors.
- `reference_periodic_a080_r030_timeseries.csv` - corresponding exact-periodic reference mode amplitude and diagnostics.
- `open_operator_fft_vs_direct_N64.csv` - componentwise validation of the FFT open operator against explicit direct summation.
- `stochastic_paired_maxima_N256_Z25.csv` - one per-realization paired stochastic diagnostic at `N=256`, `Z=25` for the pilot.
- `stochastic_pilot_seed_and_config.csv` - exact pilot seed and configuration.

The full deterministic propagation-map data remain in `data/dynamic/propagation_map_a080_r030_N256.csv`.
