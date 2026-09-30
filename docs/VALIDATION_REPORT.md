# Validation report

**Local clean-clone validation passed. Public GitHub validation has not run.
This project is not declared submission-ready.**

Run completed (UTC): 2026-09-30T15:32:21.854986+00:00

Validated source commit: `f172120261791c42dd47ab3d08e422c518f19ad9`.
The final evidence commit adds documentation/results only; the generator,
weights, dependencies and inference assets remain those of this validated commit.
Organizer commit: `5c8a5d8e2551c8cf572d3d3bfcfe7633b109d91e`.
Host: Windows x86-64; Python 3.13.15; uv 0.12.21; CPU inference.

## Exact official-validator invocation

From the handoff workspace (a local repository path, not a public GitHub URL):

```powershell
$env:UV_PYTHON_INSTALL_DIR=Join-Path (Get-Location) 'work/python'
$env:UV_CACHE_DIR=Join-Path (Get-Location) 'work/uv-cache'
$env:PATH=(Join-Path (Get-Location) 'work/tooling/bin')+';'+$env:PATH
$env:PYTHONUTF8='1'
$localRepo=(Resolve-Path outputs/collagen-amp-challenge-2027).Path
& .\work\tooling\bin\uv.exe run --project outputs/collagen-amp-challenge-2027 python outputs/collagen-amp-challenge-2027/scripts/verify_submission.py $localRepo --dir work/clean-validation2 --antibacterial-fasta outputs/collagen-amp-challenge-2027/data/antibacterial.fasta
```

The supplied validator is unmodified and was compared byte-for-byte with the
organizer Git blob. It performed a fresh Git clone, `uv sync`, two generation
runs, count/alphabet/length/uniqueness checks, top membership, full-library exact
overlap, top Levenshtein novelty, and byte reproducibility checks. The clean
clone contained neither generated outputs nor a preexisting virtual environment.
Dependency downloads may use uv's ordinary package cache; generated sequences
do not use a cache. All committed data/model hashes matched the manifest in
the clean clone.

## Additional checks

`uv run python scripts/audit_submission.py` passed with:

- 50,000 unique library sequences; 100 unique ranked top members.
- 39,448 organizer reference records; zero exact library overlaps.
- Maximum top/reference Levenshtein ratio: 0.702702702703 (allowed <=0.80).
- Maximum internal top ratio: 0.700000000000 (required <0.80).
- Maximum charge gain: 0.999996837732 (required <=1).
- 80 structural and 20 probable collagen designs; 84 parent proteins; at most two per parent.
- All activity, MIC and toxicity thresholds passed.
- Adapted anchor predictions matched historical values within floating-point tolerance (maximum absolute discrepancy 3.55e-15).
- 811 prohibited-source windows excluded, leaving 24,604 inference seeds.
- 56,704 unique plausible learned proposals; 0 library novelty replacements and 0 top novelty rejections were needed.

The entire output directory was moved to `work/run1-generated` and generation
was repeated from absent outputs. Both FASTAs matched the first run and the
official validator's separate clean-clone outputs exactly. SHA-256:

```text
c64fbf4dd67400ef835996b904b91d8721b15ea49be842ac3cf43c5cdc8074dc  generate/library.fasta
8de201de7c19c719ab09c8b4ff9e8d201af06cb1d034031b0c5dc754825e157e  generate/top.fasta
```

## Corrections and limits

No model weights or scientific score/eligibility thresholds were changed.
The inherited seed universe did contain prohibited sources; the new explicit
source gate corrected that, changing candidate membership and historical ranks.
Novelty enforcement, serial random-forest inference, deterministic sorting/LF
output, count/quota failure checks and dependency pins were added. An initial
run failed because the sandbox denied multiprocessing pipes; serial evaluation
fixed inference. Git clone also needed an approved sandbox exception for process
pipes. No validator checks were bypassed. An initial CRLF-based reference digest
was corrected to the organizer's exact Git-blob bytes before clean-clone validation.

The public repository has since been created at
https://github.com/scottjarmusch/collagen-amp-challenge-2027 using the signed-in
browser with explicit user authorization. Public-URL validation is pending.
No Kaggle submission has been made. See `PUBLISHING_HANDOFF.md`.

Scientific disclosure limits remain visible in `TRAINING_DATA.md`: missing
original training/split/preprocessing records and public-only checkpoints,
unverified original full-protein snapshot, JNP TFK-18 sequence/chemistry caveat,
and conflicting upstream Hemolytik2 license statements. Local FASTA validity
does not certify these scientific/provenance questions. Public holdout claims
are inherited metrics, not newly reproduced independent validation. No
before/after feature-attribution comparison was invented.

## Full official validation result

Also preserved verbatim in `official-local-validation.log`:

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
[1] Cloning C:\Users\salja\Documents\Codex\2026-09-30\files-mentioned-by-the-user-codex\outputs\collagen-amp-challenge-2027 â†’ work\clean-validation2
[2] Installing dependencies
[3] Generating library
Running: uv run --no-sync generate
[4] Verifying full library
[5] Verifying top list
[6] Checking library overlap with outputs\collagen-amp-challenge-2027\data\antibacterial.fasta
[7] Checking top similarity with outputs\collagen-amp-challenge-2027\data\antibacterial.fasta
[8] Checking reproducibility
Running: uv run --no-sync generate

All checks passed. Submission is valid!
```
