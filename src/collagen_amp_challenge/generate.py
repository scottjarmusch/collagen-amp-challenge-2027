from pathlib import Path
import argparse, random, hashlib, json
import joblib, numpy as np, pandas as pd, torch
from torch import nn
from rapidfuzz import fuzz
from Levenshtein import ratio
from .features import featurize, feature_names, charge

AA='ACDEFGHIKLMNPQRSTVWY'; VOCAB={a:i+1 for i,a in enumerate(AA)}; PAD=0; MASK=21; MAXLEN=30
HYD=set('AILMFWVYC'); ROOT=Path(__file__).resolve().parents[2]

class MaskedBiLSTM(nn.Module):
    def __init__(self):
        super().__init__(); self.emb=nn.Embedding(22,48,padding_idx=PAD); self.lstm=nn.LSTM(48,64,batch_first=True,bidirectional=True); self.head=nn.Sequential(nn.Linear(128,96),nn.ReLU(),nn.Linear(96,20))
    def forward(self,x,pos):
        h,_=self.lstm(self.emb(x)); b=torch.arange(x.size(0),device=x.device); return self.head(h[b,pos])

def frac(s,aset): return sum(a in aset for a in s)/len(s)
def longest(s,aset):
    best=cur=0
    for a in s:
        if a in aset:cur+=1;best=max(best,cur)
        else:cur=0
    return best

def plausible(s):
    return frac(s,HYD)<=.55 and longest(s,HYD)<=4 and frac(s,set('FWY'))<=.15 and s.count('C')<=2 and 'CC' not in s and 1.5<=charge(s)<=9.5

def local_hyd(s,w=9): return max(sum(a in HYD for a in s[i:i+w])/w for i in range(len(s)-w+1))

def propose_batch(lm,parents,nmut_arr,rng,temp=.85,topk=8,batch_size=768):
    seqs=[list(s) for s in parents]; poslists=[]
    for s,nm in zip(seqs,nmut_arr): poslists.append(sorted(rng.sample(range(len(s)),int(nm))))
    np_rng=np.random.default_rng(rng.randrange(2**32))
    for step in range(max(nmut_arr)):
        active=[i for i,p in enumerate(poslists) if step<len(p)]
        for st in range(0,len(active),batch_size):
            ids=[];poss=[];idxs=[];olds=[]
            for i in active[st:st+batch_size]:
                pos=poslists[i][step]; cur=seqs[i]; arr=[VOCAB[a] for a in cur]+[PAD]*(MAXLEN-len(cur)); old=cur[pos]; arr[pos]=MASK
                ids.append(arr);poss.append(pos);idxs.append(i);olds.append(old)
            with torch.no_grad(): logits=lm(torch.tensor(ids),torch.tensor(poss)).cpu().numpy()/temp
            for row,i,old in zip(logits,idxs,olds):
                row[VOCAB[old]-1]=-1e9; ix=np.argpartition(row,-topk)[-topk:]; pr=np.exp(row[ix]-row[ix].max());pr/=pr.sum(); seqs[i][poslists[i][step]]=AA[np_rng.choice(ix,p=pr)]
    return [''.join(s) for s in seqs],poslists

def generate_chunk(lm,seeds,chunk,start,nseeds,reps=3):
    rng=random.Random(8000+chunk); ss=seeds.iloc[start:start+nseeds]; parent_rows=[]
    for _,r in ss.iterrows():
        for _ in range(reps): parent_rows.append((r,rng.randint(3,6)))
    parents=[x[0].Lowest_Peptide for x in parent_rows]; nm=np.array([x[1] for x in parent_rows]); variants,positions=propose_batch(lm,parents,nm,rng)
    rows=[];seen=set()
    for (r,nmut),s,pos in zip(parent_rows,variants,positions):
        if s==r.Lowest_Peptide or s in seen or not plausible(s):continue
        seen.add(s);rows.append({'sequence':s,'parent_sequence':r.Lowest_Peptide,'parent_entry':r.Entry,'parent_protein':r.Protein_Name,'parent_organism':r.Organism,'source_class':r.Source_Class_Strict2,'mutation_count':nmut,'mutation_positions':';'.join(str(p+1) for p in pos),'parent_identity':1-nmut/20,'parent_adapted_rank':int(r.adapted_rank_strict2)})
    return pd.DataFrame(rows)

def write_fasta(seqs,path,prefix):
    with open(path,'w',encoding='utf-8',newline='\n') as f:
        for i,s in enumerate(seqs,1):f.write(f'>{prefix}_{i:05d}\n{s}\n')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--n-sequences',type=int,default=50000);ap.add_argument('--top-k',type=int,default=100);ap.add_argument('--seed',type=int,default=42);args=ap.parse_args()
    if args.seed!=42: raise ValueError('Released checkpoint fixes seed=42 for reproducibility.')
    if args.n_sequences != 50000 or args.top_k != 100:
        raise ValueError('Released selection policy requires exactly 50000 library and 100 top sequences.')
    random.seed(args.seed);np.random.seed(args.seed);torch.manual_seed(args.seed)
    torch.set_num_threads(1);torch.use_deterministic_algorithms(True)
    ref_path=ROOT/'data'/'antibacterial.fasta'
    ref_bytes=ref_path.read_bytes()
    if hashlib.sha256(ref_bytes).hexdigest() != 'cbbeac64ba95746d87961e8ad9dd0849ae8058d15a300b2e7f6990730ca521e9':
        raise ValueError('Organizer reference differs from the pinned snapshot; revalidate before updating it.')
    refs=[];parts=[]
    for line in ref_bytes.decode().splitlines():
        line=line.strip()
        if line.startswith('>'):
            if parts:refs.append(''.join(parts));parts=[]
        elif line:parts.append(line.upper())
    if parts:refs.append(''.join(parts))
    ref_set=set(refs);refs=sorted(ref_set)
    out=ROOT/'generate';out.mkdir(exist_ok=True)
    seeds=pd.read_csv(ROOT/'data'/'collagen_seed_candidates.csv').sort_values('adapted_rank_strict2').reset_index(drop=True)
    nematodes=pd.read_csv(ROOT/'data'/'nematode_collagen_reference.tsv',sep='\t')
    nematode_genera=set(nematodes.Organism.str.split().str[0])
    forbidden_name=seeds.Protein_Name.str.contains(r'collagen[ -]like|repeat|cuticl|collagenase|collagen[ -]binding|hydroxylase|metalloproteinase|\bMMP\b|\(EC ',case=False,regex=True,na=False)
    forbidden_taxon=seeds.Entry.isin(nematodes.Entry)|seeds.Organism.str.split().str[0].isin(nematode_genera)
    excluded=seeds[forbidden_name|forbidden_taxon].copy()
    excluded.to_csv(out/'excluded_seeds.csv',index=False)
    seeds=seeds[~(forbidden_name|forbidden_taxon)].reset_index(drop=True)
    if not seeds.Source_Class_Strict2.isin(['structural_collagen','probable_collagen']).all():
        raise ValueError('Unexpected seed provenance class.')
    print(f'Seed audit: {len(seeds)} eligible, {len(excluded)} excluded',flush=True)
    ck=torch.load(ROOT/'checkpoint'/'masked_peptide_lm.pt',map_location='cpu',weights_only=True);lm=MaskedBiLSTM();lm.load_state_dict(ck['state_dict']);lm.eval()
    amp=joblib.load(ROOT/'checkpoint'/'amp_classifier_jnp_adapted.joblib');mic=joblib.load(ROOT/'checkpoint'/'mic_ecoli_jnp_adapted.joblib');tox=joblib.load(ROOT/'checkpoint'/'toxicity_jnp_adapted.joblib')
    # Serial tree aggregation avoids both nondeterministic floating-point sums
    # and multiprocessing pipe requirements on restricted Windows hosts.
    tox.set_params(n_jobs=1)
    frames=[]; starts=list(range(0,len(seeds),2000))
    for chunk,start in enumerate(starts):
        frames.append(generate_chunk(lm,seeds,chunk,start,min(2000,len(seeds)-start),3))
        print(f'Proposed chunk {chunk+1}/{len(starts)}',flush=True)
    df=pd.concat(frames,ignore_index=True).drop_duplicates('sequence')
    X=np.vstack([featurize(s) for s in df.sequence]);pa=amp.predict_proba(X)[:,1];lv=mic.predict(X);pt=tox.predict_proba(X)[:,1]
    df['amp_probability']=pa;df['pred_ecoli_log10_mic']=lv;df['pred_ecoli_mic_uM']=10**lv;df['toxicity_probability']=pt;df['ml_selectivity_score']=np.log((pa+1e-6)/(1-pa+1e-6))-lv+np.log((1-pt+1e-6)/(pt+1e-6))
    df=df.sort_values(['ml_selectivity_score','amp_probability','pred_ecoli_mic_uM','sequence'],ascending=[False,False,True,True],kind='stable').reset_index(drop=True)
    if len(df)<args.n_sequences: raise RuntimeError(f'Only {len(df)} unique plausible variants generated; requested {args.n_sequences}.')
    chosen=list(range(args.n_sequences));next_candidate=args.n_sequences;replaced=[]
    for slot,idx in enumerate(chosen):
        if df.iloc[idx].sequence not in ref_set:continue
        while next_candidate<len(df) and df.iloc[next_candidate].sequence in ref_set:next_candidate+=1
        if next_candidate>=len(df):raise RuntimeError('Insufficient reference-exact-free replacement candidates.')
        replaced.append({'slot':slot+1,'original_candidate_rank':idx+1,'replacement_candidate_rank':next_candidate+1})
        chosen[slot]=next_candidate;next_candidate+=1
    lib=df.iloc[chosen].copy().reset_index(drop=True);lib['library_rank']=np.arange(1,len(lib)+1)
    F=pd.DataFrame(np.vstack([featurize(s) for s in lib.sequence]),columns=feature_names())
    for c in ['charge','boman','hmoment','apv','alpha','aromatic_frac','hydrophobic_frac']:lib[c]=F[c].values
    lib['max9_hydrophobic_fraction']=[local_hyd(s) for s in lib.sequence]
    lib.to_csv(out/'library_annotated.csv',index=False);write_fasta(lib.sequence,out/'library.fasta','collagen_lm')
    # Final top-100 selection uses a constrained rather than unconstrained objective.
    # The global predictors strongly reward cationicity, so we prohibit the cheap
    # optimization route of adding >+1 net charge relative to the natural parent.
    # This keeps designs close to the JNP collagen phenotype while still requiring
    # predicted antibacterial activity and safety.
    lib['parent_charge']=[charge(s) for s in lib.parent_sequence]
    lib['delta_charge']=lib.charge-lib.parent_charge
    lib['revised_score']=lib.ml_selectivity_score - 0.35*np.maximum(lib.delta_charge,0) - 0.10*(1-lib.parent_identity)*20/5
    mask=((lib.amp_probability>=.90)&(lib.pred_ecoli_mic_uM<=16)&(lib.toxicity_probability<=.35)&
          (lib.delta_charge<=1.0)&(lib.charge<=8.25)&
          (lib.charge>=2.5)&(lib.apv<=.245)&(lib.alpha>=.95)&(lib.hmoment>=.40)&(lib.hmoment<=1.15)&
          (lib.aromatic_frac<=.15)&(lib.hydrophobic_frac>=.15)&(lib.hydrophobic_frac<=.50)&
          (lib.max9_hydrophobic_fraction<=6/9)&(~lib.sequence.str.contains('CC'))&(lib.sequence.str.count('C')<=2))
    el=lib[mask].sort_values(['revised_score','ml_selectivity_score','amp_probability','sequence'],ascending=[False,False,False,True],kind='stable');quota={'structural_collagen':80,'probable_collagen':20};cnt={k:0 for k in quota};parents={};top=[];novelty_rejected=0
    for _,r in el.iterrows():
        c=r.source_class;e=str(r.parent_entry);s=r.sequence
        if c not in quota or cnt[c]>=quota[c] or parents.get(e,0)>=2:continue
        # Same function and strict >0.80 rejection as the unmodified organizer validator.
        if any(ratio(s,ref)>0.80 for ref in refs):
            novelty_rejected+=1;continue
        if all(fuzz.ratio(s,x.sequence)<80 for x in top):
            top.append(r);cnt[c]+=1;parents[e]=parents.get(e,0)+1
            if len(top)==args.top_k:break
    if len(top)!=100 or cnt!=quota:
        raise RuntimeError(f'Selection cannot satisfy the unchanged top-100 constraints: {len(top)} selected, {cnt}.')
    T=pd.DataFrame(top);T['top100_rank']=np.arange(1,len(T)+1);T.to_csv(out/'top_annotated.csv',index=False);write_fasta(T.sequence,out/'top.fasta','collagen_lm_top')
    (out/'generation_audit.json').write_text(json.dumps({'seed':args.seed,'proposal_count':len(df),'library_exact_replacements':replaced,'top_novelty_rejections':novelty_rejected,'source_counts':cnt,'parent_count':len(parents)},indent=2)+'\n',encoding='utf-8')
    print(f'Generated {len(lib)} library sequences and {len(T)} top sequences.')
if __name__=='__main__':main()
