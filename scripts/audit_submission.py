"""Independent constraint audit using the unmodified organizer implementation."""
from pathlib import Path
import hashlib
import json
import joblib
import numpy as np
import pandas as pd
import Levenshtein
from rapidfuzz import fuzz
import verify_submission as official
from collagen_amp_challenge.features import charge, featurize

ROOT = Path(__file__).resolve().parents[1]

def main():
    out = ROOT / 'generate'
    library = official._verify_sequences(out / 'library.fasta')
    official._verify_top(out / 'top.fasta', library, 100)
    _, refs = official._read_fasta(ROOT / 'data/antibacterial.fasta')
    _, top = official._read_fasta(out / 'top.fasta')
    official._verify_no_overlap(library, set(refs))
    official._veritfy_max_simularity(set(top), set(refs))
    table = pd.read_csv(out / 'top_annotated.csv')
    assert list(table.sequence) == top
    assert (table.amp_probability >= .90).all()
    assert (table.pred_ecoli_mic_uM <= 16).all()
    assert (table.toxicity_probability <= .35).all()
    gains = [charge(s)-charge(p) for s,p in zip(table.sequence,table.parent_sequence)]
    assert max(gains) <= 1.0
    assert table.parent_entry.value_counts().max() <= 2
    assert table.source_class.value_counts().to_dict() == {'structural_collagen':80,'probable_collagen':20}
    internal = max(fuzz.ratio(a,b)/100 for i,a in enumerate(top) for b in top[i+1:])
    assert internal < .80
    excluded = set(pd.read_csv(out/'excluded_seeds.csv').Entry)
    all_rows = pd.read_csv(out/'library_annotated.csv')
    assert set(all_rows.sequence) == library
    assert not (set(all_rows.parent_entry) & excluded)
    seeds = pd.read_csv(ROOT/'data/collagen_seed_candidates.csv')
    seed_pairs = set(zip(seeds.Entry, seeds.Lowest_Peptide))
    assert all((e,p) in seed_pairs for e,p in zip(all_rows.parent_entry,all_rows.parent_sequence))
    anchors = pd.read_csv(ROOT/'docs/jnp_public_vs_domain_adapted.csv')
    X = np.vstack([featurize(s) for s in anchors.sequence])
    deviations = {}
    for name,col in [('amp_classifier','adapted_amp'),('mic_ecoli','adapted_mic'),('toxicity','adapted_toxicity')]:
        model=joblib.load(ROOT/f'checkpoint/{name}_jnp_adapted.joblib')
        if name == 'toxicity': model.set_params(n_jobs=1)
        predicted=10**model.predict(X) if name=='mic_ecoli' else model.predict_proba(X)[:,1]
        deviation=float(np.max(np.abs(predicted-anchors[col].to_numpy())))
        assert np.allclose(predicted,anchors[col],rtol=1e-8,atol=1e-10), (name,deviation)
        deviations[name]=deviation
    result={
        'library_count':len(library), 'top_count':len(top),
        'reference_records':len(refs), 'unique_references':len(set(refs)),
        'library_exact_overlap':0,
        'top_max_reference_levenshtein_ratio':max(Levenshtein.ratio(s,r) for s in top for r in set(refs)),
        'top_max_internal_ratio':internal,
        'max_charge_gain':max(gains),
        'source_counts':table.source_class.value_counts().to_dict(),
        'parent_count':int(table.parent_entry.nunique()),
        'excluded_seed_count':len(pd.read_csv(out/'excluded_seeds.csv')),
        'anchor_max_absolute_prediction_deviation':deviations,
        'sha256':{n:hashlib.sha256((out/n).read_bytes()).hexdigest() for n in ['library.fasta','top.fasta']},
    }
    (ROOT/'docs/LOCAL_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
