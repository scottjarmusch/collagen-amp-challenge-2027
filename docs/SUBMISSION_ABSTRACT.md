# Method abstract

We present a collagen-domain-adapted generative antimicrobial peptide model.
The submitted peptides are model-generated analogues of natural collagen-derived
sequence windows, not unchanged fragments extracted from collagen and not
unconstrained de novo sequences. The generator starts from recorded natural
20-residue parent windows and proposes amino-acid substitutions; the submitted
top 100 each contain 3–6 substitutions relative to their recorded parent window.
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
