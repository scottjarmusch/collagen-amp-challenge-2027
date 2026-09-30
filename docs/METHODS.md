# Methods and reproducibility boundary

## Preserved learned architecture

The supplied public AMP and MIC predictors are histogram gradient boosting
models; toxicity is a random forest. Each uses 432 features: 12 physicochemical
or composition summaries, 20 amino-acid frequencies and 400 dipeptide frequencies.
The feature implementation is retained verbatim. Charge is calculated at pH 7;
MIC is predicted as log10(MIC/uM), then exponentiated. Toxicity is a fitted binary
class probability, not an assay-calibrated human-cell toxicity measurement.

The masked BiLSTM has a 22-token embedding of width 48, bidirectional hidden
width 64 and a 128-to-96-to-20 output head. The handoff reports training on 5,068
AMP plus 5,068 collagen sequences and masked validation accuracy 0.1688. The
checkpoint contains weights, alphabet and maximum length, but not optimizer
state or the original collagen master table. The recovered language-model script
filters the unlabelled public-positive CSV to canonical 18-30-mers, samples an equal
number of collagen candidates after excluding two refined provenance categories,
shuffles with seed 42, and uses a seeded random 10%/minimum-1,000 validation split.
The recovered collagen master table now permits deterministic replay of the
archived masked-LM corpus selection: 5,068 public AMP plus 5,068 collagen rows,
matching the archived counts. Explicit membership and split manifests are frozen.
Historical metrics contain no corpus hash, so this is a reconstruction from the
supplied inputs, not independently authenticated original row membership.
See [masked-LM membership](MASKED_LM_MEMBERSHIP.md). Full-length UniProt source
FASTA/release metadata remain unavailable; no LM retraining is claimed.

Seeds are precomputed 20-residue windows with UniProt accessions, parent names,
organisms and positions. Full-length protein FASTAs and the historical retrieval
query/release were not included. Thus inference is reproducible from committed
windows; rebuilding those windows from the original full-length source is not
currently independently reproducible. Original source classification remains
an annotation, particularly for unreviewed probable collagens.

## Explicit compliance correction to source eligibility

The inherited table had 25,415 windows and incorrectly retained prohibited source
classes. Before generation, exclude nematode accessions and genera identified in
the committed UniProt `taxonomy_id:6231 AND protein_name:collagen` snapshot, plus
protein names containing collagen-like, repeat, cuticle, collagenase,
collagen-binding, hydroxylase, metalloproteinase, MMP, or an enzyme EC annotation.
Canonical collagen alpha-chain names ending in `chain-like` remain eligible.
The taxonomy snapshot was retrieved 2026-09-30 in 500-record pages. Genus matching
also covers obsolete accessions or changes in species synonyms. This is an
explicit conservative correction to meet the requested biological exclusions,
not a change to learned weights, ranking thresholds or quotas. The original
seed table remains intact for audit; excluded rows are written on each run.
The old training corpus is unavailable, so the learned checkpoint's historical
exposure to excluded proteins cannot be retrospectively certified.

## Proposal and ranking

Eligible seeds are sorted by inherited adapted rank, chunked in groups of 2,000,
and sampled three times each. Existing per-chunk seeds `8000 + chunk` are retained
as part of the fixed seed-42 schedule. Proposals change 3-6 positions using
temperature 0.85 and top-8 masked-residue probabilities. Duplicates and implausible
variants are removed. The existing plausibility filter constrains hydrophobic
fraction/runs, aromatic content, cysteine and absolute charge.

Fitness = logit(P(AMP)) - predicted log10(MIC) + logit(1-P(toxicity)), using
the original 1e-6 stabilizers. Stable sorting uses score, AMP probability,
MIC and sequence as the final tie-break. The first 50,000 form the initial
library. An exact organizer-reference match is replaced only in its own slot
by the next ranked unused nonmatching proposal beyond the initial library.
Exhaustion is an error; replacement does not fabricate or heuristically mutate
a peptide. The audit records original and replacement proposal ranks.

Top-100 ranking retains the supplied revised score:
fitness - 0.35*max(charge gain,0) - 0.10*(1-parent identity)*20/5.
Eligibility retains AMP >=0.90, MIC <=16 uM, toxicity <=0.35, charge gain <=1,
charge in [2.5,8.25], APV <=0.245, alpha >=0.95, hydrophobic moment [0.40,1.15],
aromatic fraction <=0.15, hydrophobic fraction [0.15,0.50], maximum nine-residue
hydrophobic fraction <=6/9, <=2 cysteines and no adjacent cysteines.
Greedy ranked selection imposes 80 structural/20 probable, <=2 per parent, and
internal RapidFuzz ratio <80. Each candidate must also satisfy the organizer's
exact `Levenshtein.ratio` test (reject strictly >0.80) against every reference.
That screen uses the same implementation as the unmodified official validator.
Skipping a failing candidate selects the next eligible ranked candidate without
relaxing any constraint. All 100 must belong to the library.

## Determinism

Python/NumPy/Torch are seeded, Torch uses one thread and deterministic algorithms,
output newlines are LF, and ties have an explicit sequence key. The generator
does not read cached outputs. Runtime dependency versions are locked. Byte identity
is checked on repeated clean-output runs on the tested host; equivalence across
different CPU architectures or future dependency updates is not claimed.
