# Validation report

**The unmodified official validator passed against the public GitHub repository.**

Repository: https://github.com/scottjarmusch/collagen-amp-challenge-2027

Verified commit: `2a1a67d934aa01fa1d71f6ed9c49bc4586ecd776`. Run completed 2026-09-30; evidence recorded 2026-09-30T15:51:32.153359+00:00.
Organizer commit: `5c8a5d8e2551c8cf572d3d3bfcfe7633b109d91e`.
Environment: Windows x86-64, Python 3.13.15, uv 0.12.21, CPU PyTorch 2.10.0,
scikit-learn 1.8.0; all dependency versions are in `uv.lock`.

## Command and result

From the local repository root, with uv/Git on PATH and PYTHONUTF8=1:

```sh
uv run python scripts/verify_submission.py https://github.com/scottjarmusch/collagen-amp-challenge-2027 --dir ../../work/github-validation
```

The validator cloned GitHub into a new directory, ran `uv sync`, generated both
FASTAs, checked sequence constraints/novelty, and generated them again to compare
bytes. The full result is below and in `official-github-validation.log`.
Ordinary dependency caches were allowed; no generated-file cache was used.
The organizer validator and reference are unchanged Git-blob copies.

The later release commit adds training reconstruction, public-only weights, split
manifests, reports and submission notes. It does not change any generator code,
runtime dependency or inference asset. `VALIDATED_RUNTIME.json` proves those
inputs match the tested GitHub commit. The recorded SHA is the actual tested
commit, not a claim that a later evidence commit was independently rerun.

## Sequence and reproducibility checks

- Exactly 50,000 unique canonical library sequences and 100 ranked top members.
- All peptides are 20 residues; the output specification is linear, unmodified,
  all-L standard residues with free termini. FASTA itself does not encode chemistry.
- Zero exact overlaps with the organizer's 39,448 references.
- Maximum top/reference Levenshtein ratio 0.702702702703 (limit <=0.80).
- Maximum internal top similarity 0.70 (required <0.80).
- Top: AMP >=0.90, E. coli MIC <=16 uM, toxicity <=0.35, charge gain <=+1.
- 80 structural / 20 probable collagen; 84 parents; at most two peptides per parent.
- 811 prohibited-source windows were excluded; 24,604 seeds remain.
- 56,704 learned plausible proposals; no novelty replacements were needed.
- Original run, absent-output regeneration, local clean clone and public GitHub
  clean clone produced identical FASTA bytes. The official validator also verified
  its own two consecutive generations.

```text
c64fbf4dd67400ef835996b904b91d8721b15ea49be842ac3cf43c5cdc8074dc  generate/library.fasta
8de201de7c19c719ab09c8b4ff9e8d201af06cb1d034031b0c5dc754825e157e  generate/top.fasta
```

## Training and documentation checks

Recovered archived code was used to reconstruct the original deterministic splits
(seed 42), then explicit manifests were frozen. These are reconstructed splits,
not falsely described as previously exported originals. Public AMP/MIC/toxicity
models retrained from the manifests reproduce every archived metric and holdout
prediction exactly. Adapted holdout metrics also reproduce exactly. The public
and adapted toxicity holdout is the same 289 sequences. AMP holdouts differ;
their common 425-example comparison is reported separately. Melittin's public
AMP/MIC scores are not independent because it was in their training sets.

Public/adapted weights, source scripts, source URLs, preprocessing, split manifests,
methods, JNP limitations, MIT project license and third-party notices are included.
The optional toxicity interpretation is an exact tree-path decomposition, not
mechanistic evidence. It does not alter inference or selection.

## Remaining review items

The machine validator does not certify dataset rights or biological provenance.
The recovered collagen master table now permits deterministic replay of the
archived masked-LM corpus selection: 5,068 public AMP plus 5,068 collagen rows,
matching the archived counts. Explicit membership and split manifests are frozen.
Historical metrics contain no corpus hash, so this is a reconstruction from the
supplied inputs, not independently authenticated original row membership.
See [masked-LM membership](MASKED_LM_MEMBERSHIP.md). Full-length UniProt source
FASTA/release metadata remain unavailable; no LM retraining is claimed. Hemolytik2 data match the official CC BY 4.0 deposit (see DATA_LICENSES.md). The TFK-18 archive
record contains an extended 33-residue sequence requiring assay-level clarification.
Those caveats are visible in TRAINING_DATA.md and JNP_DOMAIN_ADAPTATION.md.
Accordingly, technical validity is verified; blanket full-compliance/co-authorship
eligibility is not asserted. No Kaggle submission or rule acceptance was performed.

## Full official GitHub validation output

```text
Seed audit: 24604 eligible, 811 excluded
Proposed chunk 1/13
Proposed chunk 2/13
Proposed chunk 3/13
Proposed chunk 4/13
Proposed chunk 5/13
Proposed chunk 6/13
Proposed chunk 7/13
Proposed chunk 8/13
Proposed chunk 9/13
Proposed chunk 10/13
Proposed chunk 11/13
Proposed chunk 12/13
Proposed chunk 13/13
Generated 50000 library sequences and 100 top sequences.
Seed audit: 24604 eligible, 811 excluded
Proposed chunk 1/13
Proposed chunk 2/13
Proposed chunk 3/13
Proposed chunk 4/13
Proposed chunk 5/13
Proposed chunk 6/13
Proposed chunk 7/13
Proposed chunk 8/13
Proposed chunk 9/13
Proposed chunk 10/13
Proposed chunk 11/13
Proposed chunk 12/13
Proposed chunk 13/13
Generated 50000 library sequences and 100 top sequences.
[1] Cloning https://github.com/scottjarmusch/collagen-amp-challenge-2027 → ..\..\work\github-validation
[2] Installing dependencies
[3] Generating library
Running: uv run --no-sync generate
[4] Verifying full library
[5] Verifying top list
[6] Checking library overlap with data\antibacterial.fasta
[7] Checking top similarity with data\antibacterial.fasta
[8] Checking reproducibility
Running: uv run --no-sync generate

All checks passed. Submission is valid!
```
