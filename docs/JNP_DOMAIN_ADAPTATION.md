# JNP domain adaptation

Public AMP, E. coli MIC and hemolysis models were trained first. The handoff
reports subsequently adding the small JNP peptide/control set with 20x sample
weight. These are domain examples, not a standalone training set and not an
independent validation cohort. The intended distinction is potent low-cytotoxicity
collagen peptides versus inactive collagen peptides and more cytotoxic controls.

The supplied before/after toxicity predictions are:

| Anchor | Public-only | JNP-adapted |
|---|---:|---:|
| REI-26 | 0.297 | 0.105 |
| TRR-26 | 0.267 | 0.105 |
| LRS-21 | 0.227 | 0.082 |
| GFD-30 | 0.262 | 0.106 |
| TFK-18 record | 0.277 | 0.727 |
| LL-37 | 0.137 | 0.678 |
| melittin | 0.921 | 0.962 |

Values come from `jnp_public_vs_domain_adapted.csv`, preserved from the handoff.
The low-toxicity active collagen anchors retain high adapted AMP probability;
GPE-19 and GEK-25 remain poor AMP predictions. Public-only ranks and predictions
are retained in the original comparison CSVs. Adapted ranks measure fit to
adaptation examples, not prospective predictive performance.

The handoff's public toxicity holdout AUC changes 0.92719 -> 0.92984; average
precision 0.79724 -> 0.79689; balanced accuracy 0.77268 -> 0.76334. Thus performance
is broadly retained in these reported summaries, not literally unchanged. Split
IDs and public-only checkpoints were not supplied, so this release cannot
independently certify an untouched holdout or reproduce the original comparison.
The supplied disclosure says JNP matches were removed from public holdouts before
adaptation; the filtering and exact example identities require the missing
training pipeline to audit.

Public hemolysis and JNP FMCA human-cell viability are different endpoints.
The adapted model combines those labels as a practical proxy; its probability
must not be read as a measured viability fraction. Qualitative JNP viability
constraints were encoded as labels without digitizing graphical values.

Sequence identity caveat: the handoff's `TFK-18` record contains 33 residues.
The published table displays an extended region for this named truncated peptide.
The stored sequence is preserved because changing it would change the adaptation
experiment and require retraining. Resolve the experimentally assayed sequence
and termini against the paper/supplement before claiming sequence-level biological
validation of this anchor. A canonical-letter sequence alone cannot prove free
termini or unmodified experimental chemistry.

The scientific rationale is that public data may associate cationic/helical
features with toxicity, while collagen domain examples teach an exception.
That does not identify a molecular mechanism. A valid before/after feature
attribution requires both original fitted toxicity models and the same feature
background. Only adapted weights were supplied, so comparative attribution was
not fabricated. Charge, Boman index, hydrophobic moment, APV, alpha propensity,
composition and dipeptides are all present in the retained feature space.

Source: Jarmusch et al., J. Nat. Prod. 2026, 89, 242-250,
[doi:10.1021/acs.jnatprod.5c01318](https://doi.org/10.1021/acs.jnatprod.5c01318).
The article is CC BY 4.0; attribution is retained.
