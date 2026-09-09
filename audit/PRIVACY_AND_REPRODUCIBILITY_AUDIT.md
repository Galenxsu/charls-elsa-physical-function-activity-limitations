# Privacy and reproducibility audit

## Outcome

The publishable tree currently contains repository documentation and empty source/data/output structure only. Five code files are held in a separate non-release quarantine pending authorship confirmation. No CHARLS or ELSA participant-level data, identifiers, member lists, restricted questionnaires, credentials, or historical ELSA result files were copied.

## Automated checks

- Restricted data extensions in candidate: 0.
- Local Windows user paths in candidate: 0.
- Credential/private-key patterns: 0.
- Superseded ELSA result anchors: 0.
- Quarantined Python syntax compilation before removal: 5/5 files passed.
- Models executed during packaging: 0.
- Public upload or release: 0.

Before quarantine, the scanner flagged one line because the CHARLS primary script writes a local de-identified person-period file after explicitly excluding `pid` and `ID`. No generated file was included. If the code is later restored, all derived row-level outputs must remain outside version control.

## Reproducibility classification

`CODE_QUARANTINED_PENDING_AUTHORSHIP_CONFIRMATION`

The package has not been rerun against restricted data in a clean environment. Therefore numerical reproduction remains pending an authorized local run and cannot be marked verified by this packaging audit.

## Scientific result source

Only the accepted ELSA implementation following the Wave 6 smoking structural-skip correction is represented. Earlier ELSA results are excluded and classified as `SUPERSEDED_BY_CORRECTED_WAVE6_SMOKING_MAPPING`.
