# GitHub release checklist

Repository target: `ozalloum/finite-domain-fidelity-long-range-lattices`
Release target: `v1.0.0`

## Before making the repository public

- [ ] Create the GitHub repository with the exact intended name.
- [ ] Upload the contents of this package at the repository root.
- [ ] Confirm that no internal audit, submission, reviewer, or private notebook files were added accidentally.
- [ ] Confirm the author email/contact information to be displayed publicly, if any.
- [ ] Decide whether numerical data and figures should receive a separate open license and update `LICENSE_SCOPE.md` if applicable.
- [ ] Run `python scripts/reproducibility_smoke_test.py` from a clean environment.
- [ ] Run `python scripts/run_all.py` and confirm that all eight figures regenerate.
- [ ] Check that `manifests/repository_file_checksums.csv` matches the uploaded release contents.
- [ ] Make the repository public.
- [ ] Create GitHub release `v1.0.0` from the exact public commit.
- [ ] Only after the repository and release are live, insert the repository URL and release tag into the manuscript availability statements.
- [ ] If the paper is accepted, archive the accepted computational release with a persistent DOI and update `CITATION.cff` and the manuscript if journal timing permits.
