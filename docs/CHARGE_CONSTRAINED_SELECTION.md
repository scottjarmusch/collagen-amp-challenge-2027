# Charge-constrained final selection

The original ML fitness strongly rewarded increased cationicity. To prevent a generic AMP shortcut, final top-100 selection now uses a constrained objective rather than unconstrained scalar ranking.

Eligibility requires:
- AMP probability >= 0.90
- predicted E. coli MIC <= 16 uM
- predicted toxicity probability <= 0.35
- generated net charge increase relative to its natural collagen parent <= +1.0
- absolute generated charge <= +8.25
- existing physicochemical/transmembrane guardrails

Within the eligible set, ranking uses:

`revised_score = ml_selectivity_score - 0.35 * max(delta_charge, 0) - 0.10 * mutation_count/5`

The top 100 additionally retain the 80 structural / 20 probable collagen provenance target, maximum two peptides per parent UniProt protein, and <80% pairwise internal similarity.

This changes the median charge gain in the top 100 from approximately +2 in the unconstrained selection to approximately 0, while preserving high predicted AMP probability, <=16 uM predicted MIC, and low predicted toxicity.
