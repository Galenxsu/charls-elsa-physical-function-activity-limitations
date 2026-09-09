# Controlled end-to-end reproduction summary

Run UUID: `4e668d39-22b6-40d9-82bf-3230533f5d8a`

Execution date: 2026-09-06. Outputs were written outside the repository to an independent local directory. No frozen result, manuscript, or supplement was overwritten.

## Environment

- Python 3.12.14
- NumPy 2.5.2 for analyses
- pandas 3.0.5
- SciPy 1.18.0
- statsmodels 0.14.6
- patsy 1.0.2
- pyreadstat 1.3.6
- Matplotlib 3.11.1

The figure-only process used the accessible document-library path and NumPy 2.3.5; no model or statistical estimate was calculated in that process.

## CHARLS

- Primary analysis: 2,499 participants, 4,449 person-period records, 674 events.
- Grip-decline conditional HR: 1.1816773140; 95% CI 1.0756086532–1.2982056907; P=0.0005033143.
- Primary continuous-result CSV: byte-identical to the accepted frozen file.
- Primary table CSV: byte-identical.
- Exploratory-extension all-coefficient CSV: byte-identical.
- Complete 21-test FDR CSV: byte-identical.

## ELSA

- Primary conceptual replication: 3,991 participants, 7,233 person-period records, 733 events.
- Grip-decline HR: 1.1152778685; 95% CI 1.0174594885–1.2225004908; P=0.0198307488.
- Seven of seven models completed successfully.
- Twenty-six of twenty-seven common aggregate CSV/XLSX files were byte-identical to the accepted Wave 6 smoking-corrected package.
- The only non-identical file was `CHARLS_ELSA_EVIDENCE_CROSSWALK.csv`. The ELSA row was unchanged. The reproduced CHARLS row correctly used the later ADL structural-skip-corrected primary estimate (HR 1.1816773140). The older estimate retained in the historical ELSA crosswalk is explicitly classified as `SUPERSEDED_BY_ADL_STRUCTURAL_SKIP_CORRECTION`. This is a documented cross-package provenance update, not an ELSA-model discrepancy.

## Figures and tables

- Figure 1 PNG: byte-identical.
- Figure 2 PNG: byte-identical.
- Table 1, Table 2, and Table 3 source-data CSV files: byte-identical.
- PDF, SVG, and DOCX byte hashes differed because their generated metadata changed; dimensions/content source and raster PNG outputs were unchanged.

## Privacy and release boundary

No participant-level output, identifier, member list, restricted document, or local execution log was copied into the repository. Only this aggregate summary and aggregate QA status are retained.

Status: `CONTROLLED_END_TO_END_REPRODUCTION_COMPLETED_WITH_DOCUMENTED_CROSSWALK_PROVENANCE_DIFFERENCE`.
