# Masked-LM corpus membership reconstruction

The supplied `amp_challenge_provenance_refined.zip` contains the historical
`master_all_candidates_with_refined_provenance.csv`. Its original bytes are
bundled under `data/training/`, with Git newline conversion disabled for this
file to preserve its source hash. The frozen row indices are zero-based parsed
data-row indices, excluding the header; they are sensitive to source row order.

## Exact replay

The reconstruction executes the corpus-selection statements from the archived
`train_masked_lm.py`, redirecting only file paths. A separate indexed replay
asserts identical AMP, selected collagen and shuffled corpus sequence order.
Two independent executions match; rerunning checks rather than overwrites the
frozen manifests.

1. Read public `Sequence`, convert to strings and uppercase, retain canonical
   sequences of length 18–30. No deduplication: 5,068 rows remain.
2. Read all 75,244 master rows in their supplied order. Exclude only labels
   `collagen-like/repeat` and `nematode/cuticle collagen`: 54,042 rows remain.
   The historical script does not add another collagen length/alphabet filter.
3. Seed Python random with 42, shuffle eligible collagen rows, select the first
   `len(amp)` rows: 5,068. Concatenate AMP then collagen and shuffle again using
   the continuing Python RNG state. The resulting 10,136 sequences are unique
   and compatible with the archived alphabet/maximum length.
4. Replay `random_split` with a separate Torch generator seeded 42 and
   `nval=max(1000,int(.1*N))`: 9,123 training and 1,013 validation rows.
   There are zero exact shared sequences, but no homology grouping was used.

## Historical collagen composition

| Archived provenance label | Selected rows |
|---|---:|
| structural/canonical collagen | 2,309 |
| other probable collagen | 1,538 |
| collagen-associated | 1,099 |
| other collagen-source annotation | 122 |

This broader training universe must not be described as exclusively strict
structural/probable collagen. Later inference filters do not retroactively change
the checkpoint's training corpus. The reconstruction preserves the original
labels and selection rather than silently applying newer exclusions.

## Frozen evidence and verification boundary

`data/splits/masked_lm_corpus.csv` records corpus order, source row index,
sequence hash, parent accession, original provenance and split. The train and
validation CSVs preserve each Torch subset's index order. Full source/manifest
hashes, software versions and counts are in `MASKED_LM_MEMBERSHIP.json`.

The recorded counts match `archived_masked_lm_metrics.json` exactly. However,
that historical file and the checkpoint contain no corpus hash or original
row-index manifest. The new hashes freeze this deterministic reconstruction;
there is no historical hash against which to authenticate it. Matching counts
alone cannot prove that an earlier source file had identical rows/order.
The weights and reported masked validation accuracy have not been retrained or
independently reproduced in this task. Original full-protein FASTA and UniProt
retrieval release remain missing. Inference weights, code and FASTAs are unchanged.

## SHA-256 fingerprints

- Selected collagen row-index order: `97ad6fd5ab3d24c21bc5b60d67456b65a7ce43e591dc2247b4ab8106594cea5f`
- Selected collagen sequence order: `76231eb30815dc37b57bb69e2d6d57f5025eb5375bdedc7ff1c47762bbc2dd60`
- Combined shuffled corpus order: `f3d744b8389851a01f8d5d6477ab207f1e3fc86fe3f835b0ed76ff3ffeebf5a6`

Each sequence/index fingerprint hashes UTF-8 values joined with LF and a final LF.
