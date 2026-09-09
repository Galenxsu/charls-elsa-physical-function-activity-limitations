# Provenance and result-version policy

## ELSA

The sole accepted ELSA source is the analysis completed after the Wave 6 smoking structural-skip correction on 3 September 2026. In that implementation, valid HESKA responses take precedence; HESKA=-1 is classified as not currently smoking only when the preceding HESMK response equals 2 (never smoked). Official derived smoking variables support concordance checks and do not overwrite a valid HESKA response.

All earlier ELSA result sets are classified as `SUPERSEDED_BY_CORRECTED_WAVE6_SMOKING_MAPPING` and are excluded from this candidate.

## CHARLS

The included scripts correspond to the ADL structural-skip-corrected primary and exploratory analysis implementations. Earlier pre-correction results are not release sources.

## Candidate modifications

Only filesystem configuration was changed in copied scripts: hard-coded local roots were replaced by environment variables or repository-relative output directories. Original source files were not modified.

