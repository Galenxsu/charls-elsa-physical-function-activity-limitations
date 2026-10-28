# Sensitivity analyses added in v1.2.0

These scripts reproduce the sensitivity analyses reported in Supplementary
Tables S16-S18. They require authorized local copies of CHARLS and ELSA data;
no participant-level data are distributed in this repository.

Run the scripts from the repository root after defining the paths documented
in `config/paths.example.env`:

1. `python src/sensitivity/01_sensitivity_analyses.py`
2. `python src/sensitivity/02_model_bc_comparisons.py`
3. `python src/sensitivity/03_charls_bmi_cesd10.py`
4. `python src/sensitivity/verify_reference_results.py`

By default, newly generated aggregate outputs are written below
`outputs/sensitivity`. Published aggregate reference outputs are retained in
`reference_results/sensitivity` for comparison.

The optional V1 cohort path is used only to reproduce an historical CES-D10
missingness audit. If `CHARLS_V1_DEIDENTIFIED_COHORT` is not supplied, that
historical audit is skipped; the S18 models remain available.
