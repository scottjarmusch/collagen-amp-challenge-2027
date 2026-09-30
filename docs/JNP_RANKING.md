# JNP peptides in the de novo design space

Historical handoff results against the original library. Source exclusions in
this release change the generated comparison population; the ranks below have
not been recalculated for the current library. Model predictions are checked
separately by `scripts/audit_submission.py`.

Important: adapted-model ranks are in-sample because JNP labels were used for domain adaptation. Public-only ranks are the non-circular comparison.

|Peptide|Adapted constrained rank|Public-only rank vs 50k|Public-only percentile|Pred MIC uM|Pred toxicity|
|---|---:|---:|---:|---:|---:|
|TRR-26|1|4275|91.5%|3.69|0.105|
|LRS-21|3|13227|73.5%|3.82|0.082|
|REI-26|6|20276|59.5%|4.93|0.105|
|GFD-30|32|26082|47.8%|10.93|0.106|
|LL-37|949|16367|67.3%|1.54|0.678|
|TFK-18|2126|23106|53.8%|3.74|0.727|
|melittin|4121|24817|50.4%|1.49|0.962|
|LEL-28|>7225|37955|24.1%|40.11|0.058|
|GPE-19|>7225|45004|10.0%|31.10|0.094|
|SPE-22|>7225|20851|58.3%|29.82|0.046|
|GEK-25|>7225|49226|1.5%|111.00|0.064|
