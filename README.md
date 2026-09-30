# Collagen-domain-adapted AMP generation

A masked peptide language model proposes collagen-derived analogues. Public-data
AMP, E. coli MIC and hemolysis predictors, adapted using a small JNP collagen
peptide set, rank the proposals. The released checkpoints and original scientific
selection thresholds are preserved. This is an inference release, not a
reconstruction of the missing original model-training pipeline.

**Status:** see [validation report](docs/VALIDATION_REPORT.md). A local validator
pass does not establish public-GitHub validation or competition acceptance.

## Run from a clean checkout

Install Git and uv, then run from the repository root:

```sh
uv sync
uv run generate
```

Python 3.13 is specified in `.python-version`; `uv.lock` pins dependencies.
The checkpoints require scikit-learn 1.8.0. Generation uses CPU PyTorch 2.10.0,
one Torch thread and deterministic algorithms. Allow several minutes and enough
memory for approximately 60,000 feature vectors. No GPU or network is required
after dependency installation. All model and inference data assets are local.

The command freshly generates `generate/library.fasta` (50,000 unique peptides)
and `generate/top.fasta` (100 ranked peptides, each a library member). Default seed
42 is fixed by this release. Non-default counts or seeds are rejected because
the released quota policy is specifically 80 structural plus 20 probable collagen
designs. It also writes annotation tables, excluded seeds and a generation audit.
Generated files are ignored by Git and are never read as inference inputs.

## Selection and provenance

Top eligibility requires AMP probability >=0.90, predicted E. coli MIC <=16 uM,
toxicity probability <=0.35, and charge increase <=+1 relative to the natural
parent. Existing physicochemical guardrails remain. Top selection enforces <80%
internal similarity, at most two peptides per parent, and the 80/20 source quota.
Failure to meet any quota raises an error instead of emitting a shortened top list.

The organizer's reference snapshot removes exact library overlaps and top
Levenshtein ratios >0.80. Exact-overlap replacements preserve unaffected library
slots. Top replacements follow the next eligible rank and retain all constraints.
No thresholds are relaxed. See [methods](docs/METHODS.md) for the explicitly
documented biological-source correction discovered during this handoff.

## Evidence and limitations

The handoff reports public-only AMP AUC 0.9651, MIC MAE 0.4407 log10(uM), and
toxicity AUC 0.9272. Adapted results are 0.9635, 0.4446 and 0.9298 respectively.
These historical metrics are retained, not independently re-estimated: public-only
weights, split identifiers and training scripts were not supplied. Historical
public-only prediction/ranking CSVs remain available as the comparison without
JNP adaptation. They do not establish an unbiased external benchmark for newly
generated peptides. JNP-adapted rankings are in-sample adaptation behavior.

Read [JNP adaptation](docs/JNP_DOMAIN_ADAPTATION.md),
[training disclosure](docs/TRAINING_DATA.md), and [third-party notices](docs/THIRD_PARTY_NOTICES.md).
Predictions do not establish antimicrobial efficacy or clinical safety.

## Organizer validation

After pushing to the actual public repository:

```sh
uv run python scripts/verify_submission.py https://github.com/OWNER/collagen-amp-challenge-2027
```

The unmodified validator and reference are pinned to organizer commit
`5c8a5d8e2551c8cf572d3d3bfcfe7633b109d91e`. The validator clones, syncs, generates,
checks novelty and generates again. Store the full log, date and checked commit
in `docs/VALIDATION_REPORT.md`. This repository's MIT license covers original
code and released original model assets; third-party data retain their own terms.
