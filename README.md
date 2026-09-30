# CHARLS–ELSA objective physical function study: analysis code

This repository contains the analysis code supporting **Decline in objective physical function and subsequent activity limitation: longitudinal associations after accounting for earlier function in CHARLS and ELSA**.

## Release status

The authors confirmed that the five included Python files were developed within this research project, contain no unlicensed collaborator or third-party code, and may be released under the MIT License. A controlled local reproduction and prepublication privacy and integrity audit were completed before the public release.

## Intended scope

The code covers variable construction, complete-case selection, person-period construction, prespecified regression models, sensitivity analyses, and submission figures/tables. CHARLS is the primary finding cohort. ELSA is an independently implemented conceptual replication with harmonized but non-equivalent outcomes.

Version 1.1.0 aligns the public CHARLS primary workflow with the grip-focused submission analysis. The ELSA code corresponds to the Wave 6 smoking structural-skip correction used in the conceptual replication. Earlier implementations and results are superseded and are not included.

The current grip-focused CHARLS primary analysis includes 3,767 participants, 6,689 person-periods and 1,010 events. The hazard ratio per one-standard-deviation greater grip decline is 1.213772 (95% CI 1.119082–1.316474; P=2.9406873e-06), displayed in the submission as HR 1.214 (95% CI 1.119–1.316). The earlier biomarker-complete estimate of HR 1.181677 is retained only as a secondary physiological analysis and is not the primary submission estimate.

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
3. Define `CHARLS_ITEM_MAP_FILE` as the local path to the validated ADL/IADL item-level coding map, then run `src/charls/01_primary_analysis.py`.
4. Run `src/charls/02_exploratory_extensions.py`.
5. Run `src/elsa/01_conceptual_replication.py`.
6. Run `src/figures/build_submission_figures_tables.py` after the required aggregate results are available.

The restored scripts must intentionally fail when required environment variables or files are absent and must not silently substitute data.

## License boundary

The MIT License applies only to original code in this repository. It does not apply to CHARLS or ELSA data, questionnaires, codebooks, data dictionaries, proprietary software, or third-party code. See `audit/LICENSE_AND_THIRD_PARTY_AUDIT.md` for the completed review.

## Citation

Repository: https://github.com/Galenxsu/charls-elsa-physical-function-activity-limitations

Archived release v1.1.0: https://doi.org/10.5281/zenodo.23050625

The previous v1.0.0 archive remains available at https://doi.org/10.5281/zenodo.22699531.
