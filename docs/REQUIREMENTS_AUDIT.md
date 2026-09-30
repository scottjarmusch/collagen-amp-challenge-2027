# Submission requirements audit â€” 2026-09-30

Sources: [organizer requirements](https://github.com/szczurek-lab/amp-challenge-2027),
[Kaggle submission requirements](https://www.kaggle.com/competitions/amp-challenge/overview/submission-requirements),
and [competition rules](https://www.kaggle.com/competitions/amp-challenge/rules).

| Requirement | Status / evidence |
|---|---|
| Public GitHub repository | Verified public under scottjarmusch/collagen-amp-challenge-2027 |
| Learned generative method, model weights and code | Present; original masked LM and adapted predictors unchanged |
| Abstract, selection/ranking and usage docs | SUBMISSION_ABSTRACT.md, METHODS.md, README.md |
| Permissive OSI project license | MIT; third-party datasets separately attributed |
| Defined Python, uv, committed lockfile | Python 3.13, pyproject.toml and uv.lock |
| Clean sync and default generate entry point | Passed official GitHub-clone validator |
| 50,000 unique library and ranked top 100 | Passed; top is a subset of library |
| Standard alphabet and length 8-50 | Passed; generated sequences are 20-mers |
| Linear/free termini/no chemical modifications | Explicit design specification; do not request amidation or other modifications during synthesis |
| Fixed seed / byte reproducibility | Passed seed-42 clean-output and clean-clone checks |
| Exact library novelty and organizer top Levenshtein | Passed against unchanged organizer reference |
| Charge, activity, toxicity, diversity and provenance quotas | Passed independent audit with unchanged thresholds |
| Public training reconstruction | Passed; frozen reconstructed manifests and exact metric/prediction reproduction |
| JNP interpretation | In-sample adaptation explicitly disclosed; same toxicity holdout confirmed |
| Training sources and filters | Disclosed in TRAINING_DATA.md; LM master table recovered and membership reconstructed; original full-protein release and historical corpus fingerprint remain unavailable |
| Data rights | UniProt CC BY 4.0 verified; Hemolytik2 matches official CC BY 4.0 deposit; project-derived collagen data released with attribution (DATA_LICENSES.md) |
| Assay-to-sequence mapping | TFK-18 extended-record caveat disclosed; not silently relabeled |
| Kaggle package | FASTAs and supporting notes prepared after GitHub validation |
| Kaggle account/team eligibility, one entry/model and rule acceptance | User must verify; not inferred from repository validation |
| Actual Kaggle submission | Authorized by the owner after validation; Kaggle signed in and rules accepted; writeup submission in progress |

The owner removed the independently uploaded candidate CSV in commit `f3523c2`.
The validated generated FASTAs remain unchanged.

The current Kaggle page displayed a close time of October 1, 2026 at 00:00
Europe/Copenhagen (September 30 at 22:00 UTC). Recheck the live page before upload.
Competition-specific rules allow one submission per model and require a generative
method. The foundational rules include a public-code-sharing provision referring
to Kaggle discussion/notebooks; include the repository link in the appropriate
competition channel when completing your submission. No message was posted here.

Technical validator success is established. The unresolved provenance/rights and
assay caveats above prevent a blanket assertion that every scientific or
co-authorship eligibility condition has been independently certified.
