# Reproducibility status

The scripts are preserved from the final analysis provenance and sanitized to replace local absolute paths with environment variables. After author authorization, a controlled end-to-end local reproduction was completed using legally accessible data and an independent output directory.

The core CHARLS and ELSA results and raster figures/tables matched the accepted aggregate sources. The documented crosswalk provenance difference is described in `audit/CONTROLLED_END_TO_END_REPRODUCTION_SUMMARY.md`. Public release remains subject to a separate author authorization.

Known requirements:

- Python; NumPy; pandas; SciPy; statsmodels; patsy; pyreadstat; Matplotlib.
- Local access to protocol/mapping inputs referenced by the scripts.
- Adequate memory for repeated-wave cohort files.

Random seeds recorded in the analysis scripts are retained. Models and scientific definitions were not altered during path sanitization.
