import os, re, json, math, joblib, warnings
from pathlib import Path
from collections import Counter
import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor, RandomForestClassifier
from sklearn.metrics import roc_auc_score, average_precision_score, accuracy_score, balanced_accuracy_score, mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupShuffleSplit

ROOT=Path('/mnt/data/collagen_amp_ml')
PUB=Path('/mnt/data/ML-AMP')
MODEL=ROOT/'models'; OUT=ROOT/'outputs'; MODEL.mkdir(exist_ok=True); OUT.mkdir(exist_ok=True)
import sys
sys.path.insert(0,str(ROOT/'src'))
from collagen_amp_ml.features import featurize, canonical

AA=set('ACDEFGHIKLMNPQRSTVWY')

def seq_group(s,k=4):
    # deterministic low-complexity family grouping: minimum lexical k-mer + length bin.
    # conservative enough to separate many close analog families while remaining cheap.
    km=[s[i:i+k] for i in range(len(s)-k+1)] or [s]
    return f'{len(s)//5}:{min(km)}'

def Xmat(seqs): return np.vstack([featurize(s) for s in seqs])

def clean_seq(s, minlen=8,maxlen=50):
    s=canonical(s)
    if s is None or not(minlen<=len(s)<=maxlen): return None
    return s

def group_split(df,label_col,test_size=.2):
    groups=np.array([seq_group(s) for s in df.sequence])
    y=df[label_col].values
    sp=GroupShuffleSplit(n_splits=1,test_size=test_size,random_state=42)
    return next(sp.split(np.zeros(len(df)),y,groups))

# ---------- AMP classifier: peptide-level Veltri benchmark ----------
parts=[]
for fn,label in [('veltri_positive.csv',1),('veltri_negative.csv',0)]:
    d=pd.read_csv(PUB/fn)
    d['sequence']=d['Sequence'].map(clean_seq)
    d=d.dropna(subset=['sequence'])[['sequence']].drop_duplicates()
    d['amp_label']=label; d['source']=fn
    parts.append(d)
amp=pd.concat(parts,ignore_index=True)
conf=amp.groupby('sequence').amp_label.nunique(); bad=set(conf[conf>1].index)
amp_bal=amp[~amp.sequence.isin(bad)].drop_duplicates('sequence').sample(frac=1,random_state=42).reset_index(drop=True)
tr,te=group_split(amp_bal,'amp_label')
X=Xmat(amp_bal.sequence); y=amp_bal.amp_label.values
amp_model=HistGradientBoostingClassifier(max_iter=350,learning_rate=.05,max_leaf_nodes=31,l2_regularization=1.0,random_state=42)
amp_model.fit(X[tr],y[tr]); pp=amp_model.predict_proba(X[te])[:,1]
amp_metrics={'n_total':len(amp_bal),'n_train':len(tr),'n_test':len(te),'class_counts':amp_bal.amp_label.value_counts().to_dict(),'roc_auc':roc_auc_score(y[te],pp),'average_precision':average_precision_score(y[te],pp),'balanced_accuracy':balanced_accuracy_score(y[te],pp>=.5),'conflicting_sequences_removed':len(bad)}
joblib.dump(amp_model,MODEL/'amp_classifier_public.joblib')

# ---------- MIC regression: HydrAMP value = log10(MIC/uM), E coli ----------
mic=pd.read_csv(PUB/'mic_data.csv')
mic['sequence']=mic.sequence.map(clean_seq); mic=pd.to_numeric(mic.value,errors='coerce').to_frame('value').join(mic.sequence)
mic=mic.dropna().groupby('sequence',as_index=False).value.median()
trm,tem=group_split(mic,'value')
Xm=Xmat(mic.sequence); ym=mic.value.values
mic_model=HistGradientBoostingRegressor(max_iter=400,learning_rate=.04,max_leaf_nodes=23,l2_regularization=2.0,random_state=42)
mic_model.fit(Xm[trm],ym[trm]); pm=mic_model.predict(Xm[tem])
mic_metrics={'n_total':len(mic),'n_train':len(trm),'n_test':len(tem),'mae_log10_mic':mean_absolute_error(ym[tem],pm),'rmse_log10_mic':mean_squared_error(ym[tem],pm)**0.5,'r2':r2_score(ym[tem],pm)}
joblib.dump(mic_model,MODEL/'mic_ecoli_public.joblib')

# ---------- Hemolysis / cytotoxicity classifier ----------
h=pd.read_csv(PUB/'Hemolytik2_complete_data.csv')
h['sequence']=h.seq.map(clean_seq)
# Restrict to challenge-like linear, natural-AA, mostly free termini to reduce chemistry confounding.
flt=h.sequence.notna() & h.non_nat.isna() & h.lyn_cyc.astype(str).str.lower().eq('linear')
flt &= h.nter.astype(str).str.lower().isin(['free','h'])
flt &= h.cter.astype(str).str.lower().isin(['free','oh','carboxyl (coo-)','carboxylic acid'])
h=h[flt].copy()
pat50=re.compile(r'(?i)(LC50|HC50|HD50|EC50|LD50|HL50)\s*(=|>|<)\s*([0-9.]+)(?:\s*±\s*[0-9.]+)?\s*[µμu]M')
patpct=re.compile(r'(?i)([<>~]?)\s*([0-9.]+)(?:\s*±\s*[0-9.]+)?(?:\s*-\s*([0-9.]+))?\s*%\s*(?:hemolysis|hemolytic)\s*(?:at|upto|up to)?\s*([<>]?)\s*([0-9.]+)\s*[µμu]M')
def tox_label(r):
    x=str(r.activity)
    m=pat50.search(x)
    if m:
        op=m.group(2); v=float(m.group(3))
        if op in ['=','<'] and v<=64: return 1
        if op in ['=','>'] and v>=128: return 0
    m=patpct.search(x)
    if m:
        pop=m.group(1); p1=float(m.group(2)); p2=float(m.group(3)) if m.group(3) else p1
        pct_lo=min(p1,p2); pct_hi=max(p1,p2); conc=float(m.group(5)); cop=m.group(4)
        # high hemolysis by/below 128 µM => toxic; very low hemolysis through >=64 µM => safer.
        if pop != '<' and pct_lo>=50 and conc<=128: return 1
        if pop != '>' and pct_hi<=10 and conc>=64: return 0
    if str(r.non_hem) in ['Non-hemolytic','Low hemolytic']: return 0
    if re.search(r'(?i)\bnon[- ]?hemolytic\b|\bpoor hemolytic\b',x): return 0
    return np.nan
h['toxic_label']=h.apply(tox_label,axis=1)
# add HydrAMP HC50 exact/censored
hh=pd.read_csv(PUB/'hemolysis.csv'); hh['sequence']=hh.sequence.map(clean_seq)
def hlab(r):
    if pd.isna(r.sequence): return np.nan
    if bool(r.hc50_censored) and r.hc50_uM>=128: return 0
    if (not bool(r.hc50_censored)) and r.hc50_uM<=64: return 1
    if (not bool(r.hc50_censored)) and r.hc50_uM>=128: return 0
    return np.nan
hh['toxic_label']=hh.apply(hlab,axis=1)
tox=pd.concat([h[['sequence','toxic_label']],hh[['sequence','toxic_label']]],ignore_index=True).dropna()
nuniq=tox.groupby('sequence').toxic_label.nunique(); badtox=set(nuniq[nuniq>1].index)
tox=tox[~tox.sequence.isin(badtox)].drop_duplicates('sequence').reset_index(drop=True)
trt,tet=group_split(tox,'toxic_label')
Xt=Xmat(tox.sequence); yt=tox.toxic_label.astype(int).values
tox_model=RandomForestClassifier(n_estimators=700,max_depth=14,min_samples_leaf=2,max_features='sqrt',class_weight='balanced_subsample',n_jobs=-1,random_state=42)
tox_model.fit(Xt[trt],yt[trt]); pt=tox_model.predict_proba(Xt[tet])[:,1]
tox_metrics={'n_total':len(tox),'class_counts':{str(k):int(v) for k,v in tox.toxic_label.value_counts().items()},'n_train':len(trt),'n_test':len(tet),'roc_auc':roc_auc_score(yt[tet],pt) if len(set(yt[tet]))>1 else None,'average_precision':average_precision_score(yt[tet],pt),'balanced_accuracy':balanced_accuracy_score(yt[tet],pt>=.5),'conflicting_sequences_removed':len(badtox)}
joblib.dump(tox_model,MODEL/'toxicity_public.joblib')

# ---------- JNP held-out evaluation BEFORE adaptation ----------
j=pd.read_csv(ROOT/'data'/'jnp_collagen_amp_ml_anchors.csv')
rows=[]
for _,r in j.iterrows():
    s=clean_seq(r.sequence)
    xx=featurize(s).reshape(1,-1)
    amp_p=float(amp_model.predict_proba(xx)[0,1])
    logmic=float(mic_model.predict(xx)[0]); mic_uM=10**logmic
    tox_p=float(tox_model.predict_proba(xx)[0,1])
    rows.append({'name':r['name'],'sequence':s,'source':r['source'],'observed_mic_ecoli_uM':r['mic_ecoli_uM'],'observed_cytotox_label':r['cytotox_label'],'public_amp_probability':amp_p,'public_pred_ecoli_mic_uM':mic_uM,'public_toxicity_probability':tox_p,'public_selectivity_proxy':(1-tox_p)/(mic_uM+1e-6)})
jpred=pd.DataFrame(rows); jpred.to_csv(OUT/'jnp_public_only_predictions.csv',index=False)

metrics={'amp_classifier':amp_metrics,'mic_ecoli':mic_metrics,'toxicity':tox_metrics}
(ROOT/'outputs'/'public_model_metrics.json').write_text(json.dumps(metrics,indent=2))
print(json.dumps(metrics,indent=2))
print('\nJNP held-out public-only predictions:')
print(jpred[['name','observed_mic_ecoli_uM','public_amp_probability','public_pred_ecoli_mic_uM','observed_cytotox_label','public_toxicity_probability']].to_string(index=False))
