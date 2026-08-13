# Fidelity Horizons in Modulationally Unstable Long-Range Nonlinear Lattices

Public computational repository supporting the study:

**Fidelity Horizons in Modulationally Unstable Long-Range Nonlinear Lattices**  
Othman H. Y. Zalloum  
Department of Applied Mathematics and Physics, Palestine Polytechnic University, Hebron, Palestine  
ORCID: 0000-0001-7282-1332

## Purpose

This repository contains the Python source code, numerical data, figure-source CSV files, validation outputs, and reproduction scripts needed to inspect and reproduce the principal computational results of the associated manuscript.

The study considers a normalized long-range discrete nonlinear Schrodinger lattice with algebraically decaying coupling. The main deterministic results connect:

1. algebraic finite-domain closure error;
2. amplification of that error by modulational instability; and
3. the resulting finite propagation distance over which a finite simulation remains within a specified tolerance of the infinite-periodic reference.

The leading tested horizon scaling is

`Z_f ~ (beta/g) ln(N) + C`,

with the open finite-domain closure exponent tracking the interaction exponent in the studied power-law kernel.

## Release status

This repository is prepared as **v1.0.0** for public release.

The deterministic study is the primary evidence base. The bundled stochastic calculation is an explicitly limited paired pilot with `R=64`; it is included as a secondary distributional diagnostic and must not be interpreted as a converged rare-event or extreme-event study.

## Repository contents

- `src/` - core model and numerical operators.
- `scripts/` - theory, static, dynamic, validation, stochastic-pilot, analysis, figure, and reproducibility scripts.
- `data/theory/` - kernel-tail, dispersion, modulational-instability, and selected-mode data.
- `data/static/` - finite-domain closure errors, fits, and local exponents.
- `data/dynamic/` - deterministic error growth, fidelity horizons, fitted growth rates, and propagation maps.
- `data/validation/` - operator and time-step validation tables.
- `data/stochastic/` - paired stochastic-pilot data and summaries.
- `data/representative_raw/` - compact run-level outputs supporting independent inspection.
- `figure_data/` - combined source CSVs and panel-level source CSVs for the eight main figures.
- `figures/` - publication figures as vector PDF files. Raster PNG duplicates are intentionally omitted from the public release.
- `logs/` - compact machine-readable validation outputs retained for reproducibility.
- `manifests/` - result provenance, data inventory, and repository checksums.
- `docs/` - claim-scope, reproducibility, and release guidance.

## Quick verification

Create an environment from `requirements.txt` or `environment.yml`, then run:

`python scripts/reproducibility_smoke_test.py`

The smoke test independently checks:

- the algebraic kernel-tail exponent;
- exact periodization against the infinite Fourier symbol;
- the FFT open operator against explicit direct summation;
- the open-boundary static scaling in a representative case; and
- one frozen deterministic dynamic benchmark.

A successful run ends with:

`All reproducibility smoke tests PASS.`

## Fast data-to-figure rebuild

Run:

`python scripts/run_all.py`

This rebuilds summary tables, publication figures, and panel-level CSV files from the bundled numerical results without rerunning the full deterministic parameter campaign.

## Full deterministic recomputation

Run:

`python scripts/run_all.py --recompute`

This reruns theory/static calculations, deterministic dynamics, analysis, validation, summaries, and figures.

## Stochastic pilot

To rerun the bundled paired pilot, use:

`python scripts/run_all.py --stochastic`

The bundled stochastic calculation uses `R=64` paired realizations and is secondary evidence only. A larger paired ensemble and reference-domain campaign would be required before making converged rare-event claims.

## Headline reproducibility files

- `data/headline_results.csv` - compact numerical results quoted in the study.
- `data/table01_parameters.csv` - principal computational parameters.
- `data/table02_main_results.csv` - one-row-per-case deterministic summary.
- `manifests/results_manifest.csv` - result-to-evidence provenance.
- `manifests/claims_and_evidence.csv` - claim-scope and supporting files.
- `logs/reproducibility_smoke_test.csv` - independent smoke-test output.
- `docs/REPRODUCIBILITY.md` - concise validation summary.
- `docs/LIMITATIONS_AND_CLAIM_SCOPE.md` - explicit limits on interpretation.

## Reproducibility notes

The package includes the complete public Python implementation required for the reported calculations. The public release is script-based; the author's full development/authoring notebook is not required to reproduce the principal results and is not included in this repository release.

## Citation

If you use this repository, cite the associated manuscript and the software release. Citation metadata are provided in `CITATION.cff`.

A persistent archive DOI should be added to `CITATION.cff` after an accepted/final research release is archived.

## License

The MIT License in `LICENSE` applies to the software source code in `src/` and `scripts/`. See `LICENSE_SCOPE.md` for the current scope of licensing for data, figures, and manuscript-related materials.

## Contact

Othman H. Y. Zalloum  
Department of Applied Mathematics and Physics  
Palestine Polytechnic University  
Hebron, Palestine  
ORCID: 0000-0001-7282-1332
