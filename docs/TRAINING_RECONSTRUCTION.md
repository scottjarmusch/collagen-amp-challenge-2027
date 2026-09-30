# Reconstructed training splits and public-model reproduction

The earlier working-model archive was supplied after initial inference validation.
Its four inference checkpoints are byte-identical to those in the final handoff.
It also contains the original public-only checkpoints and training scripts.
The scripts have been preserved under `scripts/archived_training/`; their hard-coded
historical paths are not used by the portable reproduction entry point.

The original row-level split IDs were **not separately preserved**. We reconstructed
the deterministic splits from the archived preprocessing, input CSVs and seed 42,
then froze them as explicit manifests in `data/splits/`. A repeated reconstruction
must match the committed bytes or the script raises an error. No new random split
is silently substituted. Sequence IDs are SHA-256 of canonical sequence strings;
row indices refer to the reconstructed preprocessed table, not original assay rows.

```sh
uv run python scripts/reproduce_public_training.py
```

The script reloads frozen training/holdout CSVs before fitting and writes newly
trained public weights only to ignored `training_reproduction/`. Released inference
weights are not overwritten. All reported public metrics reproduce exactly, and
each retrained model's holdout predictions equal its archived checkpoint's
predictions (maximum difference 0). See TRAINING_REPRODUCTION.json.

| Task | Train | Holdout | Reproduced metric |
|---|---:|---:|---:|
| AMP | 2648 | 679 | AUC 0.9651251805069854 |
| E. coli MIC | 3494 | 825 | MAE 0.440726004347393 log10(uM) |
| Toxicity | 1172 | 289 | AUC 0.9271869750613775 |

Adapted checkpoints were also evaluated on reconstructed adapted-public holdouts;
all archived metrics reproduced. Toxicity uses the same 289 held-out peptides
before/after adaptation: AUC 0.92719 -> 0.92984, AP 0.79724 -> 0.79689, balanced
accuracy 0.77268 -> 0.76334. Performance is broadly retained, not numerically
identical and not a formal statistical equivalence test.

The MIC holdout is also identical (825 peptides). AMP holdouts differ: 679 public
versus 677 adapted, with 425 common untouched examples. On that common set,
public/adapted AUC is 0.96186/0.96009; AP 0.96431/0.96406; balanced accuracy
0.90606/0.89664. These shared comparisons avoid implying that two different AMP
holdouts are the same experiment. No JNP sequence or grouping key overlaps the
adapted holdouts. The grouping heuristic does not prove all homologues are absent.

Melittin occurs in public AMP and MIC training, so it is not an independent
public-only anchor for those endpoints. Public toxicity training contains no
exact JNP anchor. Adapted JNP examples remain in-sample by construction.

Public AMP: HistGradientBoostingClassifier, 350 iterations, learning rate .05,
31 leaves, L2 1, seed 42. MIC: HistGradientBoostingRegressor, 400 iterations,
learning rate .04, 23 leaves, L2 2, seed 42. Toxicity: 700-tree random forest,
depth 14, minimum leaf 2, sqrt features, balanced_subsample, seed 42. Serial
execution changes scheduling only. The original adaptation refits models with
20x JNP weights, uses 450 MIC iterations and 900 toxicity trees; it is not a
warm-start continuation with every hyperparameter held fixed.

The recovered collagen master table now permits deterministic replay of the
archived masked-LM corpus selection: 5,068 public AMP plus 5,068 collagen rows,
matching the archived counts. Explicit membership and split manifests are frozen.
Historical metrics contain no corpus hash, so this is a reconstruction from the
supplied inputs, not independently authenticated original row membership.
See [masked-LM membership](MASKED_LM_MEMBERSHIP.md). Full-length UniProt source
FASTA/release metadata remain unavailable; no LM retraining is claimed.
