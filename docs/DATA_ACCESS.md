# Data access and local placement

## CHARLS

Eligible researchers may apply for CHARLS data through the official CHARLS website. Data are not redistributed in this repository. After authorization, preserve the original filenames and place the files under a local directory referenced by `CHARLS_DATA_ROOT`.

## ELSA

ELSA data are available through the UK Data Service, Study Number 5050, subject to the applicable End User Licence and project conditions. Data are not redistributed here. After authorization, set `ELSA_DATA_ROOT` to the local directory containing the required Stata files.

## Prohibited repository content

- original or derived individual-level data;
- participant identifiers or hashed member lists;
- row-level analysis datasets or person-period records;
- restricted questionnaires, codebooks, or user guides;
- authentication tokens, passwords, API keys, download credentials, or local `.env` files.

No real example data are supplied. Any future example dataset must be completely synthetic and must not be created by perturbing participant records.

