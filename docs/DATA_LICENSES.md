# Data licensing and attribution

## Collagen provenance master and reconstructed membership

The recovered master contains UniProt-derived sequence windows and annotations,
with project-derived scores and provenance classifications. UniProt Consortium
database material is licensed under [CC BY 4.0](https://www.uniprot.org/help/license).
Retain UniProt Consortium attribution and the original accession identifiers.
The project-derived collagen master annotations and reconstructed membership
manifests are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
as part of this competition data disclosure, attributed to Scott A. Jarmusch.
Changes from UniProt include extraction of peptide windows, scoring, ranking,
provenance classification and selection. No UniProt endorsement is implied.
Historical UniProt release and full-protein snapshot remain unavailable.

## Hemolytik2

The official [Hemolytik2 database deposit](https://zenodo.org/records/19699377)
by Anand Singh Rathore and Gajendra Raghava explicitly identifies its data license
as CC BY 4.0. The ZIP's published MD5 is 190ec14e209aef44282660f8a2139149 and was
verified. Its complete-data CSV matches this project's training CSV exactly after
CRLF/LF normalization; parsed CSV rows also match. See DATA_LICENSE_VERIFICATION.json.
This provides a dataset-specific permissive source for the bundled training data.
The GitHub README/license disagreement is retained as historical context, not
used to relicense third-party software. No Hemolytik2 executable code is used.
Attribution: Singh, A., Raj SA, K., Rathore, A.S., and Raghava, G.P.S.,
“Hemolytic 2: An Updated Database of Hemolytic Peptides and Proteins,”
[doi:10.1021/acs.chemrestox.5c00322](https://doi.org/10.1021/acs.chemrestox.5c00322).

## Competition requirement boundary

The [organizer requirements](https://github.com/szczurek-lab/amp-challenge-2027)
require a training-data/filter summary for benchmark entry; full eligibility
additionally requires full disclosure and permissive release of non-public data.
Public-source data retain source terms; the project's MIT code license does not
replace them. Source URLs and preprocessing are in TRAINING_DATA.md. The recovered
collagen master and frozen manifests provide the previously missing corpus input.
No assertion of independently verified historical row identity or retraining of
the language-model checkpoint is made. These scientific limits are disclosed.
