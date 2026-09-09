# Authorized end-to-end reproduction preparation

No analysis was executed in this stage. The following gate must be completed locally by a researcher with valid CHARLS and ELSA access.

## Before execution

1. Clear code authorship and restore only approved scripts from quarantine.
2. Create a clean environment and capture the exact Python and package versions.
3. Set paths through local environment variables; do not edit scientific definitions or hard-code user paths.
4. Verify hashes of the approved scripts, protocols, mappings, and authorized input files.
5. Confirm that the ELSA implementation is exclusively the accepted Wave 6 smoking structural-skip version.
6. Confirm that output locations are outside the Git repository.

## Execution boundary

- Do not change exposure, outcome, risk-set, covariate, standardization, model, missing-data, multiplicity, or sensitivity-analysis definitions.
- Do not copy participant-level output, IDs, member lists, row-level logs, or restricted documentation into the repository.
- Stop on input-hash, model-definition, or result-anchor conflict.

## Permitted public outputs after successful reproduction

- aggregate sample/event checks;
- aggregate model-result concordance status without participant rows;
- software and package versions;
- script hashes and overall PASS/FAIL status;
- figure/table hashes where no restricted information is embedded.

## Required final record

Create a non-sensitive aggregate reproduction report stating whether the independently generated manuscript-level values matched the frozen results. The report must not contain local absolute paths or protected data excerpts.

Status: `PREPARED_NOT_EXECUTED`.
