"""Replay archived corpus selection; freeze membership without retraining the LM."""
from pathlib import Path
import ast
import csv
import hashlib
import json
import random
from collections import Counter
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
MASTER = ROOT/'data/training/master_all_candidates_with_refined_provenance.csv'
PUBLIC = ROOT/'data/training/unlabelled_positive.csv'
SCRIPT = ROOT/'scripts/archived_training/train_masked_lm.py'
OUT = ROOT/'data/splits'

def sha(data):
    return hashlib.sha256(data).hexdigest()

def replay():
    # Execute the archived selection assignments verbatim, redirecting only paths.
    tree = ast.parse(SCRIPT.read_text())
    selected = []
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'amp' for t in node.targets):
            selected = tree.body[tree.body.index(node):]
            break
    selected = selected[:next(i for i,n in enumerate(selected) if isinstance(n, ast.ClassDef))]
    class LocalCSV:
        @staticmethod
        def read_csv(path):
            return pd.read_csv(MASTER if Path(path).name == MASTER.name else PUBLIC)
    ns = {'pd':LocalCSV, 'PUB':PUBLIC.parent, 'AA':'ACDEFGHIKLMNPQRSTVWY', 'random':random.Random(42)}
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(SCRIPT), 'exec'), ns)
    return ns['amp'], ns['coll'], ns['seqs']

def main():
    amp, collagen, seqs = replay()
    assert (amp, collagen, seqs) == replay(), 'Non-deterministic selection'
    public = pd.read_csv(PUBLIC)
    master = pd.read_csv(MASTER)
    aa = set('ACDEFGHIKLMNPQRSTVWY')
    records = []
    for i,s in enumerate(public.Sequence.astype(str).str.upper()):
        if 18 <= len(s) <= 30 and set(s) <= aa:
            records.append(dict(source='public_amp', source_row_index=i, sequence=s, parent_entry='', provenance_refined=''))
    eligible = [dict(source='collagen', source_row_index=int(i), sequence=str(r.Lowest_Peptide), parent_entry=str(r.Entry), provenance_refined=str(r.Provenance_Refined))
                for i,r in master.iterrows() if r.Provenance_Refined not in ['collagen-like/repeat','nematode/cuticle collagen']]
    eligible_count = len(eligible)
    rng = random.Random(42)
    rng.shuffle(eligible)
    picked = eligible[:len(records)]
    assert [r['sequence'] for r in records] == amp
    assert [r['sequence'] for r in picked] == collagen
    records += picked
    rng.shuffle(records)
    assert [r['sequence'] for r in records] == seqs
    assert all(0 < len(s) <= 30 and set(s) <= aa for s in seqs)
    nval = max(1000, int(.1*len(seqs)))
    train, val = torch.utils.data.random_split(range(len(seqs)), [len(seqs)-nval,nval], generator=torch.Generator().manual_seed(42))
    train_ids = set(train.indices)
    for i,r in enumerate(records):
        r.update(corpus_index=i, sequence_sha256=sha(r['sequence'].encode()), split='train' if i in train_ids else 'validation')
    OUT.mkdir(exist_ok=True)
    hashes = {}
    for name, rows in [('masked_lm_corpus.csv',records), ('masked_lm_train.csv',[records[i] for i in train.indices]), ('masked_lm_validation.csv',[records[i] for i in val.indices])]:
        import io
        f=io.StringIO(newline='')
        writer=csv.DictWriter(f,fieldnames=list(records[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
        data=f.getvalue().encode()
        path=OUT/name
        if path.exists(): assert path.read_bytes()==data, f'Frozen manifest differs: {name}'
        else: path.write_bytes(data)
        hashes[name]=sha(data)
    counts={'n_amp':len(amp),'n_collagen':len(collagen),'n_total':len(seqs)}
    archived=json.loads((ROOT/'docs/archived_masked_lm_metrics.json').read_text())
    assert all(counts[k]==archived[k] for k in counts), (counts,archived)
    trseq={records[i]['sequence'] for i in train.indices};valseq={records[i]['sequence'] for i in val.indices}
    report={
        'status':'Deterministic reconstruction from supplied source files and archived code; not an independently preserved historical membership fingerprint.',
        'seed':42, 'master_rows':len(master), 'eligible_collagen_rows':eligible_count,
        'counts':counts, 'archived_counts_match':True, 'train_rows':len(train), 'validation_rows':len(val),
        'collagen_provenance_counts':dict(Counter(r['provenance_refined'] for r in picked)),
        'unique_collagen_sequences':len(set(collagen)), 'unique_corpus_sequences':len(set(seqs)),
        'train_validation_shared_unique_sequences':len(trseq&valseq),
        'selected_collagen_row_indices_sha256':sha(('\n'.join(str(r['source_row_index']) for r in picked)+'\n').encode()),
        'ordered_collagen_sequences_sha256':sha(('\n'.join(collagen)+'\n').encode()),
        'ordered_corpus_sequences_sha256':sha(('\n'.join(seqs)+'\n').encode()),
        'source_sha256':{p.relative_to(ROOT).as_posix():sha(p.read_bytes()) for p in [MASTER,PUBLIC,SCRIPT]},
        'manifest_sha256':hashes,
        'runtime':{'pandas':pd.__version__,'torch':torch.__version__},
        'limits':['Historical metrics contain counts but no corpus hash or row manifest; counts alone cannot prove original row identity.', 'Historical random_split is row-based, without deduplication or homology grouping.', 'No LM retraining, weight replacement, generation or candidate-selection change performed.']
    }
    (ROOT/'docs/MASKED_LM_MEMBERSHIP.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
