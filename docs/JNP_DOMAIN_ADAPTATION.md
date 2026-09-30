# JNP collagen-domain adaptation

## Why domain adaptation was used

The JNP study identified collagen-derived peptides with antimicrobial activity
comparable to classical AMP controls and substantially lower cytotoxicity toward
U-937 cells under the reported assay conditions. This is a small but biologically
informative set representing an AMP phenotype that may be underrepresented in
public databases. Eleven examples are not sufficient to establish a standalone,
generally reliable toxicity model; their role here is domain adaptation after
public-data learning.

The public data teach the model what a typical AMP looks like. The collagen-derived
JNP examples then help it learn what a potentially safer collagen AMP looks like.
This conceptual description concerns learned model behavior, not a mechanistic rule.

## Public toxicity model

The public model is learned before JNP adaptation, using Hemolytik2 and HydrAMP
experimental hemolysis data. It is a 700-tree random forest over physicochemical,
amino-acid-composition and dipeptide features. Public holdout ROC-AUC is approximately
0.927. This is principally a hemolysis-trained proxy, not a direct model of every
kind of mammalian-cell cytotoxicity.

The recovered training scripts and seed 42 were used to reconstruct and explicitly
freeze train/holdout manifests. Retraining from those manifests reproduced every
archived public metric and holdout prediction exactly. These manifests are
reconstructed originals, not previously exported row-ID files. See
[training reconstruction](TRAINING_RECONSTRUCTION.md) and
[data disclosure](TRAINING_DATA.md).

## JNP collagen peptide anchors

The anchors provide contrasting activity/toxicity examples:

- **REI-26, TRR-26, LRS-21 and GFD-30:** active collagen-derived peptides with
  comparatively low cytotoxicity.
- **GPE-19:** a collagen peptide prioritized by another AMP predictor but
  experimentally inactive in the reported assays.
- **GEK-25:** a charged collagen-derived negative control that was experimentally
  inactive, illustrating why charge alone is insufficient.
- **TFK-18:** an active collagen-derived comparator with greater cytotoxicity than
  the newly identified low-cytotoxicity collagen peptides. The preserved model
  record is an extended 33-residue sequence; see the limitations below.
- **LL-37 and melittin:** potent classical AMP controls with greater cytotoxicity
  under the JNP assay conditions.
- **LEL-28 and SPE-22:** additional collagen examples contributing activity and
  low-cytotoxicity labels; low predicted toxicity does not imply strong E. coli activity.

These contrasts discourage simplistic rules such as “collagen-derived = safe” or
“high charge = good AMP.” The models use sequence features; a collagen-origin label
is not itself a predictor feature.

## Effect of JNP adaptation

The archived adaptation code refits the model on public training examples plus
JNP examples with sample weight 20. The adapted toxicity forest has 900 trees
versus 700 in the public baseline; MIC fitting also uses 450 rather than 400
iterations. It is not a warm-start update with every other setting held fixed.

| Peptide | Public model | JNP-adapted model |
|---|---:|---:|
| REI-26 | ~0.297 | ~0.105 |
| TRR-26 | ~0.267 | ~0.105 |
| LRS-21 | ~0.227 | ~0.082 |
| GFD-30 | ~0.262 | ~0.106 |
| LEL-28 | ~0.168 | ~0.058 |
| SPE-22 | ~0.128 | ~0.046 |
| TFK-18 archive record | ~0.277 | ~0.727 |
| LL-37 | ~0.137 | ~0.678 |
| melittin | ~0.921 | ~0.962 |

These are **adaptation behaviors**, not independent test-set predictions. The JNP
examples contributed to fitting, so their adapted scores and ranks are not
prospective validation. The released checkpoint reproduces the supplied numerical
predictions to floating-point precision. GPE-19 and GEK-25 retain low adapted AMP
probability; lowering their toxicity scores does not make them good AMP designs.

### Public holdout stability

The same 289 public toxicity holdout peptides were kept outside both model fits.
Neither JNP sequences nor their archived grouping keys overlap this holdout.
The following values were reproduced after reconstructing the original splits:

| Metric | Before JNP adaptation | After JNP adaptation |
|---|---:|---:|
| ROC-AUC | 0.9272 | 0.9298 |
| Average precision | 0.7972 | 0.7969 |
| Balanced accuracy | 0.7727 | 0.7633 |

These results are consistent with changing behavior in collagen-relevant sequence
space without materially degrading general public holdout discrimination. The
metrics are not literally unchanged, and no significance or formal equivalence
test is claimed. The AMP holdouts differ, so their shared 425-example comparison
is reported separately in TRAINING_RECONSTRUCTION.md. Grouping by minimum 4-mer
and length bin is a heuristic, not exhaustive homology exclusion.

## Interpretation

A useful interpretation is that the public model associates some combinations
of cationicity, helicity, amphiphilicity and hydrophobicity with membrane toxicity.
JNP supplies counterexamples: some collagen-derived peptides retain antimicrobial
behavior without equivalent mammalian-cell toxicity. The adapted model learns
a finer distinction within AMP-like sequence space instead of globally lowering
toxicity predictions. The changes for LL-37, TFK-18 and melittin illustrate that
toxicity can increase after adaptation as well as decrease.

The **TRR-26 versus TFK-18** contrast is particularly informative. Both are
collagen-derived and antimicrobial, but the reported experiments indicate
different mammalian-cell toxicity. Their public-model scores are similar
(~0.267 and ~0.277); the adapted archive-record predictions separate strongly
(~0.105 and ~0.727). Thus collagen origin alone is not an adequate explanation.
The TFK-18 sequence-mapping caveat limits biological interpretation of this contrast.

### A domain-specific objective

A generic AMP classifier often asks, “Does this sequence resemble known
antimicrobial peptides?” This project asks a narrower question: “Does it resemble
a potent, relatively low-cytotoxicity AMP near collagen-derived sequence space?”
That is a domain-specific design objective, not a claim of universal superiority
over other AMP models or universal safety of collagen-derived peptides.

## Limitations

The adaptation set is small and in-sample. Public-only results are retained as
the unadapted comparison, but melittin was already in public AMP and MIC training;
its public-only scores for those endpoints are not independent validation.
No JNP anchor was in the public toxicity training or holdout set.

Public hemolysis and JNP FMCA U-937 viability are different endpoints. Combining
their labels is a practical modeling choice; output probability is not a measured
viability fraction or clinical safety estimate. Qualitative viability constraints
were not digitized into falsely precise graphical measurements.

The supplied TFK-18 model record contains 33 residues, matching an extended region
displayed in the published table for this named truncated peptide. The record and
weights are preserved; resolving the experimentally assayed sequence and termini
could require retraining and would change the original adaptation experiment.
Canonical sequence letters alone do not establish experimental chemistry.

The recovered collagen master table now permits deterministic replay of the
archived masked-LM corpus selection: 5,068 public AMP plus 5,068 collagen rows,
matching the archived counts. Explicit membership and split manifests are frozen.
Historical metrics contain no corpus hash, so this is a reconstruction from the
supplied inputs, not independently authenticated original row membership.
See [masked-LM membership](MASKED_LM_MEMBERSHIP.md). Full-length UniProt source
FASTA/release metadata remain unavailable; no LM retraining is claimed.
Experimental testing of generated peptides remains necessary.

## Role in generation

The learned masked sequence model proposes nearby collagen-derived analogues.
JNP-adapted AMP, MIC and toxicity predictions enter the design fitness. Selection
favors antimicrobial probability >=0.90, predicted E. coli MIC <=16 uM and toxicity
probability <=0.35, while retaining diversity and biological provenance.

Earlier unconstrained optimization tended to increase positive charge because
that could improve AMP scores. Current selection therefore permits **no more than
+1 net charge relative to the natural collagen parent**. This explicit constraint
limits an easy route toward generic hypercationic synthetic AMP designs. It works
alongside learned fitness; it does not replace the learned models with a heuristic.
The existing 80/20 provenance quota, at most two peptides per parent and internal
similarity <80% remain unchanged.

## Model interpretation and future analysis

Before/after exact tree-path probability decompositions are now available in
[TOXICITY_INTERPRETATION.md](TOXICITY_INTERPRETATION.md). The complete feature CSVs
include charge, Boman index, hydrophobic moment, APV, alpha-helical propensity,
hydrophobic/aromatic fractions, amino-acid composition and dipeptide motifs.
The decomposition reconstructs each forest prediction to <1e-12 residual.
It is path-dependent model attribution, not SHAP or a biological mechanism.
Further shared-background attribution or prospective experiments could assess
the robustness and biological relevance of these learned distinctions.

## Citation

Scott A. Jarmusch, Taj Muhammad, Ulf Göransson, Adam A. Strömstedt.
“α-Helical Peptides Encoded in Collagen Exhibit Antimicrobial Activity with Low Cytotoxicity.”
*Journal of Natural Products*, 2026, **89**, 242–250.
[DOI: 10.1021/acs.jnatprod.5c01318](https://doi.org/10.1021/acs.jnatprod.5c01318).
The article is licensed under **CC BY 4.0**.
