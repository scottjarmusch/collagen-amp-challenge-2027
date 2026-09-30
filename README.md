# Collagen-domain-adapted AMP generation

A masked peptide language model proposes collagen-derived analogues. Public-data
AMP, E. coli MIC and hemolysis predictors, adapted using a small JNP collagen
peptide set, rank the proposals. The released checkpoints and original scientific
selection thresholds are preserved. The original public training pipeline has also been recovered, its deterministic
splits reconstructed and frozen, and all three public predictors retrained with
exact reproduction of the archived holdout metrics and predictions.

**Status:** the official validator passed against the public GitHub repository.
See the [validation report](docs/VALIDATION_REPORT.md) and
[requirements audit](docs/REQUIREMENTS_AUDIT.md) for evidence and remaining
provenance/rights caveats. Validator success is not competition acceptance.

## What this model is doing

Many AMP models learn to recognize or reproduce features of known antimicrobial
peptides. This project starts with public AMP/non-AMP, MIC and hemolysis data to
learn that general landscape. It then uses the experimentally characterized
collagen-derived peptides from the JNP study as a small **domain-adaptation set**.
These examples teach a more specific distinction: some collagen peptides combine
strong antimicrobial activity with comparatively low mammalian-cell toxicity.

The public data teach the model what a typical AMP looks like; the JNP examples
help it recognize a potentially safer collagen AMP. This is a design objective,
not a proven mechanism or a claim that all collagen peptides are safe.

A learned sequence generator proposes novel analogues near natural collagen
peptides. Adapted models rank them by predicted AMP activity, MIC and toxicity.
Selection limits charge gain to **no more than +1 relative to the natural parent**,
so higher Lys/Arg content cannot become the easy route to generic hypercationic
designs. [The detailed rationale](docs/JNP_DOMAIN_ADAPTATION.md) explains the
adaptation, public holdout checks and scientific limitations.

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
These metrics have now been reproduced using archived training code and explicit
reconstructed split manifests. Public-only weights and comparisons are included.
The toxicity holdout is identical before/after adaptation. The AMP holdouts differ;
a shared 425-peptide comparison is separately reported. Melittin was in the public
AMP and MIC training sets, so that anchor's public-only scores are not independent.
JNP-adapted rankings are in-sample adaptation behavior. See
[training reconstruction](docs/TRAINING_RECONSTRUCTION.md).

To verify frozen splits and retrain the public models without overwriting released
inference weights:

```sh
uv run python scripts/reproduce_public_training.py
uv run python scripts/interpret_toxicity.py
```

Read [JNP adaptation](docs/JNP_DOMAIN_ADAPTATION.md),
[training disclosure](docs/TRAINING_DATA.md), and [third-party notices](docs/THIRD_PARTY_NOTICES.md).
Predictions do not establish antimicrobial efficacy or clinical safety.

## Organizer validation

Validate the public repository:

```sh
uv run python scripts/verify_submission.py https://github.com/scottjarmusch/collagen-amp-challenge-2027
```

The unmodified validator and reference are pinned to organizer commit
`5c8a5d8e2551c8cf572d3d3bfcfe7633b109d91e`. The validator clones, syncs, generates,
checks novelty and generates again. Store the full log, date and checked commit
in `docs/VALIDATION_REPORT.md`. This repository's MIT license covers original
code and released original model assets; third-party data retain their own terms.
