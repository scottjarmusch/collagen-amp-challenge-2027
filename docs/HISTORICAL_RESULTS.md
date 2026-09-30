# Historical result provenance

The following inherited files describe the handoff experiment, not the current
regenerated population: `strict_top100_summary.json`, `lm_generated_summary.json`,
`generator_analysis.json`, `generator_feature_shifts.csv`, and all `jnp_*rank*.csv`
files. They are retained for the public-only comparison and audit trail. Source
eligibility was corrected in this release, so historical counts, diversity,
novelty maxima and ranks must not be substituted for current validation.

`public_model_metrics.json`, `domain_adapted_public_holdout_metrics.json` and
`masked_lm_metrics.json` are supplied training/evaluation summaries; their original
splits were subsequently reconstructed and frozen from the recovered scripts,
and the public/adapted holdout metrics reproduced (TRAINING_RECONSTRUCTION.md).
`jnp_public_vs_domain_adapted.csv` contains the supplied
before/after anchor predictions. The current checkpoint's adapted predictions
are independently compared against those values in the local audit.

Current evidence is in `VALIDATION_REPORT.md`, `LOCAL_AUDIT.json`, the recorded
validator log, and freshly written `generate/generation_audit.json`.
