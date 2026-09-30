# Method abstract

We present a collagen-domain-adapted generative antimicrobial peptide model.
A masked bidirectional peptide language model proposes sequence analogues around
natural collagen-derived windows. Public-data predictors estimate antimicrobial
activity, E. coli minimal inhibitory concentration and hemolysis/toxicity.
A small experimentally characterized collagen peptide set provides upweighted
domain adaptation, preserving a learned distinction between active low-toxicity
collagen peptides, inactive collagen examples and more cytotoxic controls.
Selection combines the three learned predictions with an explicit maximum +1
charge gain relative to each natural parent, diversity constraints and parent
protein quotas. Prohibited collagen-like/repeat, nematode/cuticle and
collagen-associated enzyme/binder sources are excluded from inference.
The released generator produces 50,000 unique canonical peptides and a ranked
top 100, screened against the organizer's exact and Levenshtein novelty rules.
JNP-adapted rankings are in-sample adaptation results, not independent validation.
Public-only comparison summaries are disclosed alongside limitations in the
inherited training provenance. Generated-peptide efficacy and toxicity remain
computational predictions requiring experimental evaluation.
