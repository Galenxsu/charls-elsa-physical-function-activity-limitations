# CHARLS–ELSA objective physical function study: analysis code

This repository contains the analysis code supporting **Changes in Objective Physical Function and Subsequent Activity Limitations: Longitudinal Evidence from CHARLS and ELSA**.

## Release status

The authors confirmed that the five included Python files were developed within this research project, contain no unlicensed collaborator or third-party code, and may be released under the MIT License. A controlled local reproduction and prepublication privacy and integrity audit were completed before the public release.

## Intended scope

The code covers variable construction, complete-case selection, person-period construction, prespecified regression models, sensitivity analyses, and submission figures/tables. CHARLS is the primary finding cohort. ELSA is an independently implemented conceptual replication with harmonized but non-equivalent outcomes.

The ELSA code included here corresponds exclusively to the Wave 6 smoking structural-skip correction accepted on 3 September 2026. Earlier ELSA implementations and results are superseded and are not included.

All cross-cohort outputs use the current CHARLS primary estimate after correction of the ADL structural-skip implementation: HR 1.181677, 95% CI 1.075609–1.298206, P=0.000503 (submission display: HR 1.182, 95% CI 1.076–1.298, P<0.001). Earlier CHARLS crosswalk values are classified as `SUPERSEDED_BY_ADL_STRUCTURAL_SKIP_CORRECTION`.

## Data are not included

CHARLS and ELSA individual-level data cannot be redistributed. This repository contains no original or derived participant-level data, identifiers, member lists, restricted questionnaires, or credentials. Eligible researchers must obtain the data directly from the data custodians and place authorized files in local directories described in `docs/DATA_ACCESS.md`.

## Repository layout

```text
src/charls/   cleared CHARLS scripts
src/elsa/     cleared ELSA script
src/figures/  cleared figure/table scripts
config/       path-variable template only
docs/         data access, execution order, provenance, and limitations
audit/        release-candidate privacy, license, and integrity reports
data/         intentionally empty and ignored
outputs/      generated locally and ignored
```

## Environment

Create the environment with either:

```bash
conda env create -f environment.yml
conda activate charls-elsa-reproduction
```

or install `requirements.txt` in a clean Python environment. The supplied environment files document the supported software environment. Reproduction requires independent authorized access to the restricted cohort data.

## Execution order

1. Obtain authorized CHARLS and/or ELSA data.
2. Define the path variables shown in `config/paths.example.env` outside the repository.
3. Run `src/charls/01_primary_analysis.py`.
4. Run `src/charls/02_exploratory_extensions.py`.
5. Run `src/elsa/01_conceptual_replication.py`.
6. Run `src/figures/build_submission_figures_tables.py` after the required aggregate results are available.

The restored scripts must intentionally fail when required environment variables or files are absent and must not silently substitute data.

## License boundary

The MIT License applies only to original code in this repository. It does not apply to CHARLS or ELSA data, questionnaires, codebooks, data dictionaries, proprietary software, or third-party code. See `audit/LICENSE_AND_THIRD_PARTY_AUDIT.md` for the completed review.

## Citation

Repository: https://github.com/Galenxsu/charls-elsa-physical-function-activity-limitations

Archived release: https://doi.org/10.5281/zenodo.22699531
