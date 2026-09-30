# Training and inference data disclosure

All files are inherited from the user-provided package unless identified as new
organizer/taxonomy snapshots. `ASSET_MANIFEST.json` records sizes and SHA-256
hashes. This inventory distinguishes supplied claims from independently inspected
artifacts. No missing split, labeling threshold or historical version is invented.

| Source / file | Role | Preprocessing and modification status | Provenance / license |
|---|---|---|---|
| Veltri / AMP Scanner: `veltri_positive.csv`, `veltri_negative.csv` | Public AMP training and held-out evaluation | Handoff reports peptide-level canonical examples and removal of conflicting labels; training/split code absent. These files lack complete chemical modification metadata, so canonical letters do not establish unmodified chemistry. | Supplied via HydrAMP starter kit lineage; original source terms apply. |
| HydrAMP / GRAMPA: `mic_data.csv` | Public E. coli MIC training/evaluation | Reported target log10(MIC/uM). Exact censoring, aggregation, length filters and split IDs absent. File columns are sequence/value; modification exclusions cannot be independently reconstructed. | HydrAMP/GRAMPA public collection; original database provenance and terms apply. |
| Hemolytik2: `Hemolytik2_complete_data.csv` | Public toxicity training/evaluation | Handoff reports retaining linear canonical free-terminus peptides, conservative quantitative/qualitative binary labeling, and removing 18 conflicts. Raw file includes modified examples; exact retained rows and thresholds absent. | Upstream LICENSE.txt is GPL-3.0; README also contains conflicting MIT/noncommercial language. See notices. |
| HydrAMP experimental: `hemolysis.csv` | Additional toxicity training examples, not generated-peptide validation | HC50 in uM with censoring flag; exact binarization/censoring handling missing. Modification metadata not sufficient to independently confirm exclusion. | HydrAMP starter kit experimental results; inherited snapshot. |
| JNP: `jnp_collagen_amp_ml_anchors.csv` | Upweighted adaptation only | 11 canonical-letter anchor/control sequences; reported MICs and qualitative FMCA labels. No graphic digitization. Assay chemistry is not encoded; TFK-18 length caveat is documented separately. | Jarmusch et al. 2026, DOI 10.1021/acs.jnatprod.5c01318, CC BY 4.0 article. |
| UniProt: `collagen_seed_candidates.csv` | Natural inference parent windows; a reported subset also trained the language model | APV/helix-based 20-aa windows. Provenance exclusion is applied before inference; windows use standard residues. Parent IDs/coordinates are present, full protein FASTAs and historical release/query absent. | UniProt records; [UniProt license information](https://www.uniprot.org/help/license). Inherited snapshot date unknown. |
| `unlabelled_positive.csv`, `unlabelled_negative.csv` | Disclosed historical assets | Handoff explicitly says not used by final peptide AMP classifier; negative file includes long proteins. Exact positive-file contribution to the language model is not independently recoverable. | HydrAMP starter kit lineage; source terms apply. |
| `antibacterial.fasta` | Novelty filtering / organizer validation only | Organizer file retained unchanged. Exact overlap excluded for library; Levenshtein >0.80 excluded for top. Not used to fit model weights here. | Organizer commit 5c8a5d8e2551c8cf572d3d3bfcfe7633b109d91e, BSD-3-Clause repository notice retained. |
| `nematode_collagen_reference.tsv` | Inference exclusion audit only | UniProt taxonomy 6231 plus collagen protein-name query, paginated, retrieved 2026-09-30. Accessions and genera exclude nematode seeds. | New UniProt annotation snapshot; no model fitting. |

The handoff reports 3,327 AMP examples (2,648 train/679 test), 4,319 MIC examples
(3,494/825), and 1,461 toxicity examples (1,172/289). These are historical summary
counts, not reconstructed preprocessing results. JNP examples are not independent
validation. The exact 5,068 + 5,068 language-model training membership and validation
split were not provided. Public-only holdout metrics and predictions remain
available, but their independence cannot be re-audited without the original splits.

Generation itself requires only committed adapted weights, seed table, taxonomy
snapshot and organizer reference. No training CSV is silently used as an output
cache. Inference reproducibility does not imply training reproducibility.

Primary source links:
- [HydrAMP starter kit](https://github.com/szczurek-lab/hydramp-starter-kit)
- [HydrAMP](https://github.com/szczurek-lab/hydramp)
- [GRAMPA](https://github.com/zswitten/Antimicrobial-Peptides)
- [Hemolytik2](https://github.com/raghavagps/Hemolytik2)
- [JNP article](https://doi.org/10.1021/acs.jnatprod.5c01318)

Remaining disclosure gaps are the exact source versions, model-fitting scripts,
split IDs, retained row lists, modification decisions for datasets without
chemical metadata, and conflicting Hemolytik2 redistribution statements. These
limitations must remain visible in any submission; a FASTA-validator pass does
not resolve them.
