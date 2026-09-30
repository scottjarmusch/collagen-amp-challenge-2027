# Training, adaptation, inference and validation data

The user supplied a final inference archive and subsequently the earlier
`collagen_amp_ml_working_model.zip`. Archived code is retained under
`scripts/archived_training/`; its original absolute paths are not executed.
`scripts/reproduce_public_training.py` reuses the archived preprocessing and
split helpers, freezes sequence-level manifests, then trains from those files.
See `TRAINING_RECONSTRUCTION.md` and `TRAINING_REPRODUCTION.json` for evidence.
`ASSET_MANIFEST.json` records hashes; `data/public_data_manifest.csv` preserves
the recovered source URLs. Its historical role labels are superseded by the
verified uses below. All known data sources are public; none is identified as
proprietary or non-public in the supplied archives.

## Veltri / AMP Scanner

`data/training/veltri_positive.csv` and `veltri_negative.csv`, obtained via the
HydrAMP starter kit, train and evaluate the public AMP classifier. Strings are
stripped/uppercased; only the 20 standard residue letters and lengths 8-50 are
retained. Each class is deduplicated; cross-class conflicts are removed (zero
here). The public data are shuffled with pandas seed 42. There are 3,327 examples,
2,648 train and 679 holdout. The adapted script removes exact JNP matches before
splitting and deduplicates; its holdout has 677 examples. It does not perform the
public script's preliminary shuffle. The scripts' group-based splitting is
described below. Chemical modification metadata are absent, so alphabet filtering
does not prove every experimental source peptide was unmodified.

## HydrAMP / GRAMPA E. coli MIC

`data/training/mic_data.csv` is used for public MIC training and validation.
The same canonical 8-50 filter applies; values are converted to numeric and
missing values dropped. Duplicate sequence values are aggregated by median.
The provided value is already log10(MIC/uM); the scripts do not log it again.
There are 4,319 examples, 3,494 train and 825 holdout. Domain adaptation removes
exact JNP matches before grouping. The resulting holdout is identical to the
public model's holdout. Chemical modification and upstream assay/censoring
metadata are not present in this reduced sequence/value file; exclusions beyond
canonical sequence/length cannot be established from it.

## Hemolytik2

`data/training/Hemolytik2_complete_data.csv` supplies public toxicity training
and validation labels. In addition to canonical 8-50 residues, the code requires:
missing `non_nat`; `lyn_cyc` equal to linear; N-terminus free or H; C-terminus free,
OH, carboxyl (COO-), or carboxylic acid (case-insensitive). Other terminal
modifications and annotated non-natural/cyclic examples are excluded by this
filter. It does not explicitly filter the `ldmix` stereochemistry field; therefore
it must not be described as independently certifying all-L chemistry.

Labels follow the archived code exactly. LC50/HC50/HD50/EC50/LD50/HL50 values
with = or < and <=64 uM are toxic; = or > and >=128 uM are low-toxicity. Parsed
hemolysis >=50% at <=128 uM is toxic unless the percentage has a < sign; <=10%
at >=64 uM is low-toxicity unless the percentage has a > sign. The parser reads
the concentration inequality but does not use that sign in the percentage-rule
decision; this limitation is preserved rather than silently repaired. Explicit
Non-hemolytic/Low hemolytic fields or non/poor hemolytic text are low-toxicity.
Unresolved records are dropped.

## HydrAMP experimental hemolysis

`data/training/hemolysis.csv` adds public toxicity labels. Canonical 8-50 residue
filtering applies. Censored HC50 >=128 uM is low-toxicity. Uncensored HC50 <=64 uM
is toxic; uncensored >=128 uM is low-toxicity; intermediate records are dropped.
The reduced CSV lacks complete modification metadata. Combined with Hemolytik2,
18 conflicting sequences are removed, then sequences are deduplicated. The final
set contains 1,461 examples: 1,172 train and 289 holdout.

## JNP collagen peptide/control anchors

`data/jnp_collagen_amp_ml_anchors.csv` contains 11 examples from Jarmusch et al.,
J. Nat. Prod. 2026, DOI 10.1021/acs.jnatprod.5c01318. These are adaptation data,
not independent validation. AMP label is 1 when the minimum parsed E. coli /
S. aureus MIC is <=16 uM. Toxic label is 1 for intermediate/high/very_high
cytotoxicity categories. MIC adaptation uses uncensored parsed E. coli MIC,
converted to log10; `>` values are omitted from MIC fitting. Each adaptation
example has sample weight 20. Qualitative FMCA constraints were not digitized
into exact viability values. The stored sequences use canonical letters but
experimental termini/chemistry are not encoded. The preserved TFK-18 record has
33 residues and requires clarification against the assayed truncated peptide.
The article is CC BY 4.0; retain attribution. See JNP_DOMAIN_ADAPTATION.md.

## Learned proposal model and UniProt biological source

`unlabelled_positive.csv` supplies canonical 18-30-residue public AMP sequences
for the masked LM: 5,068 according to the archived metrics and training script.
`unlabelled_negative.csv` is disclosed but not used by the final public classifier
or masked LM; it includes long proteins. The LM balances public positives with
5,068 randomly shuffled collagen windows from the original
`master_all_candidates_with_refined_provenance.csv`, excluding the refined
collagen-like/repeat and nematode/cuticle categories. That master table is not in
either supplied archive. The training code is available, but exact historical LM
membership cannot be regenerated from the smaller final inference seed table.

`collagen_seed_candidates.csv` contains 25,415 inherited 20-aa UniProt windows,
accessions, parent names, organisms, protein lengths and peptide coordinates.
It is an inference asset, not an independent validation set. APV/helix selection
was performed upstream. Before inference, this release excludes 811 prohibited
source windows, leaving 24,604. Full-length parent FASTAs and the original UniProt
release/query are not bundled. The raw original seed table is preserved for audit.
`nematode_collagen_reference.tsv` is an added UniProt annotation snapshot queried
on 2026-09-30 (taxonomy 6231 AND protein_name:collagen), used only for inference
exclusions. No model was retrained on this reference.

## Splits and unbiased comparisons

All reconstructed predictor splits use GroupShuffleSplit(test_size=0.2,
n_splits=1, random_state=42). A group is length//5 plus the lexicographically
smallest sequence 4-mer. This is a heuristic grouping, not exhaustive homology
clustering. Explicit train/holdout manifests now contain sequence IDs, group IDs
and reconstructed row indices. They were reconstructed from archived code and
seed, not recovered as previously exported split files. Train/holdout sequence
and group separation was checked. Added JNP groups do not overlap adapted
holdout groups. The shared untouched holdouts are also exported.

Melittin was present in the original public AMP and MIC training sets, but no JNP
anchor was in the public toxicity training or holdout set. Public-only scores
remain useful unadapted comparisons; melittin's AMP/MIC scores are not independent
validation. JNP-adapted scores are always in-sample adaptation behavior.

## Organizer novelty reference

`antibacterial.fasta` is the exact organizer Git blob at commit
5c8a5d8e2551c8cf572d3d3bfcfe7633b109d91e. It is used only for novelty filtering
and validation, not fitting. Library exact matches are excluded; top Levenshtein
ratios >0.80 are excluded. The original BSD-3-Clause notice is retained.

## Licensing and provenance limits

Primary sources: [HydrAMP starter kit](https://github.com/szczurek-lab/hydramp-starter-kit),
[HydrAMP](https://github.com/szczurek-lab/hydramp),
[GRAMPA](https://github.com/zswitten/Antimicrobial-Peptides),
[Hemolytik2](https://github.com/raghavagps/Hemolytik2),
[JNP](https://doi.org/10.1021/acs.jnatprod.5c01318), and
[UniProt licensing](https://www.uniprot.org/help/license).

HydrAMP code is MIT, but original database terms still apply. Hemolytik2's
LICENSE.txt says GPL-3.0 while its README contains conflicting MIT/noncommercial
statements; that conflict is explicitly disclosed, not treated as a blanket
MIT data license. Code and original model assets have the project's MIT license;
third-party data retain their own terms. See THIRD_PARTY_NOTICES.md. Historical
download dates/releases were not recorded in the archives; exact supplied file
hashes are committed. Remaining provenance limits concern the original LM
collagen corpus/full-protein snapshot and unresolved assay-chemistry annotations.
