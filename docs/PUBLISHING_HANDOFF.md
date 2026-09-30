# Publication and submission status

Public repository: https://github.com/scottjarmusch/collagen-amp-challenge-2027

The repository was created in the account owner's signed-in browser on
2026-09-30 with explicit authorization. GitHub-URL validation is in progress;
see VALIDATION_REPORT.md for the latest evidence. No Kaggle submission has been
made. Local source and model assets are committed, with uv.lock and documentation.

After any inference change, validate a fresh clone:

```sh
uv run python scripts/verify_submission.py https://github.com/scottjarmusch/collagen-amp-challenge-2027 --dir submission
```

Use a previously nonexistent clone directory. On Windows, use a UTF-8 terminal
or set PYTHONUTF8=1 for the organizer's Unicode log output. All individual
model/data files are below 100 MB; Git LFS is not required.

The machine validator does not verify training disclosure or dataset rights.
Read TRAINING_DATA.md and THIRD_PARTY_NOTICES.md before asserting full scientific
or licensing compliance. Original preprocessing/split records, public-only model
weights, and the assayed TFK-18 sequence/chemistry are not reconstructed here.

After public validation passes, prepare the library/top FASTAs, abstract,
selection explanation, training summary and repository URL for the current
Kaggle upload interface. Do not submit automatically.
