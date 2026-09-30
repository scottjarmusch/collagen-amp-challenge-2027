from pathlib import Path
import sys, re, json, joblib
import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor, RandomForestClassifier
from sklearn.metrics import roc_auc_score, average_precision_score, balanced_accuracy_score, mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupShuffleSplit
ROOT=Path('/mnt/data/collagen_amp_ml'); PUB=Path('/mnt/data/ML-AMP'); MODEL=ROOT/'models'; OUT=ROOT/'outputs'
sys.path.insert(0,str(ROOT/'src'))
from collagen_amp_ml.features import featurize, canonical

def clean(s,minlen=8,maxlen=50):
 s=canonical(s); return s if s and minlen<=len(s)<=maxlen else None
def Xmat(seqs):return np.vstack([featurize(s) for s in seqs])
def group(s,k=4):
 km=[s[i:i+k] for i in range(len(s)-k+1)] or [s];return f'{len(s)//5}:{min(km)}'
def split(df,ycol):
 sp=GroupShuffleSplit(n_splits=1,test_size=.2,random_state=42)
 return next(sp.split(np.zeros(len(df)),df[ycol],np.array([group(s) for s in df.sequence])))
def pmic(x):
 s=str(x); cens=s.startswith('>'); s=s.lstrip('<>')
 try:return float(s),cens
 except:return np.nan,cens

j=pd.read_csv(ROOT/'data'/'jnp_collagen_amp_ml_anchors.csv'); j['sequence']=j.sequence.map(clean)
jseq=set(j.sequence)
# labels
ja=[]
for _,r in j.iterrows():
 e,ec=pmic(r.mic_ecoli_uM); s,sc=pmic(r.mic_saureus_uM)
 ja.append({'sequence':r.sequence,'amp_label':int(min(e,s)<=16),'ecoli_mic':e if not ec else np.nan,
            'toxic_label':int(r.cytotox_label in ['intermediate','high','very_high']),'name':r['name']})
ja=pd.DataFrame(ja)

# AMP public peptide benchmark
parts=[]
for fn,label in [('veltri_positive.csv',1),('veltri_negative.csv',0)]:
 d=pd.read_csv(PUB/fn);d['sequence']=d.Sequence.map(clean);d=d.dropna(subset=['sequence'])[['sequence']].drop_duplicates();d['amp_label']=label;parts.append(d)
a=pd.concat(parts,ignore_index=True);a=a[~a.sequence.isin(jseq)].drop_duplicates('sequence').reset_index(drop=True)
tr,te=split(a,'amp_label'); X=Xmat(a.sequence); y=a.amp_label.values
Xtr=np.vstack([X[tr],Xmat(ja.sequence)]); ytr=np.r_[y[tr],ja.amp_label.values]
w=np.r_[np.ones(len(tr)),np.full(len(ja),20.0)]
am=HistGradientBoostingClassifier(max_iter=350,learning_rate=.05,max_leaf_nodes=31,l2_regularization=1.0,random_state=42)
am.fit(Xtr,ytr,sample_weight=w); p=am.predict_proba(X[te])[:,1]
ampm={'public_holdout_auc':roc_auc_score(y[te],p),'public_holdout_ap':average_precision_score(y[te],p),'public_holdout_bal_acc':balanced_accuracy_score(y[te],p>=.5)}
joblib.dump(am,MODEL/'amp_classifier_jnp_adapted.joblib')

# MIC public + exact JNP E coli upweighted
m=pd.read_csv(PUB/'mic_data.csv');m['sequence']=m.sequence.map(clean);m=m.dropna(subset=['sequence']);m['value']=pd.to_numeric(m.value,errors='coerce');m=m.dropna();m=m[~m.sequence.isin(jseq)].groupby('sequence',as_index=False).value.median()
trm,tem=split(m,'value');Xm=Xmat(m.sequence);ym=m.value.values
jm=ja.dropna(subset=['ecoli_mic']).copy();Xjm=Xmat(jm.sequence);yjm=np.log10(jm.ecoli_mic.values)
Xtrm=np.vstack([Xm[trm],Xjm]);ytrm=np.r_[ym[trm],yjm];wm=np.r_[np.ones(len(trm)),np.full(len(jm),20.0)]
mm=HistGradientBoostingRegressor(max_iter=450,learning_rate=.04,max_leaf_nodes=23,l2_regularization=2.0,random_state=42)
mm.fit(Xtrm,ytrm,sample_weight=wm);pm=mm.predict(Xm[tem])
micm={'public_holdout_mae_log10':mean_absolute_error(ym[tem],pm),'public_holdout_rmse_log10':mean_squared_error(ym[tem],pm)**.5,'public_holdout_r2':r2_score(ym[tem],pm)}
joblib.dump(mm,MODEL/'mic_ecoli_jnp_adapted.joblib')

# Reuse parsed public toxicity dataset construction from full script by reading source-independent logic here
h=pd.read_csv(PUB/'Hemolytik2_complete_data.csv'); h['sequence']=h.seq.map(clean)
flt=h.sequence.notna() & h.non_nat.isna() & h.lyn_cyc.astype(str).str.lower().eq('linear')
flt &= h.nter.astype(str).str.lower().isin(['free','h']);flt &= h.cter.astype(str).str.lower().isin(['free','oh','carboxyl (coo-)','carboxylic acid']);h=h[flt].copy()
pat50=re.compile(r'(?i)(LC50|HC50|HD50|EC50|LD50|HL50)\s*(=|>|<)\s*([0-9.]+)(?:\s*±\s*[0-9.]+)?\s*[µμu]M')
patpct=re.compile(r'(?i)([<>~]?)\s*([0-9.]+)(?:\s*±\s*[0-9.]+)?(?:\s*-\s*([0-9.]+))?\s*%\s*(?:hemolysis|hemolytic)\s*(?:at|upto|up to)?\s*([<>]?)\s*([0-9.]+)\s*[µμu]M')
def tlab(r):
 x=str(r.activity);m=pat50.search(x)
 if m:
  op=m.group(2);v=float(m.group(3))
  if op in ['=','<'] and v<=64:return 1
  if op in ['=','>'] and v>=128:return 0
 m=patpct.search(x)
 if m:
  pop=m.group(1);p1=float(m.group(2));p2=float(m.group(3)) if m.group(3) else p1;lo=min(p1,p2);hi=max(p1,p2);conc=float(m.group(5))
  if pop!='<' and lo>=50 and conc<=128:return 1
  if pop!='>' and hi<=10 and conc>=64:return 0
 if str(r.non_hem) in ['Non-hemolytic','Low hemolytic'] or re.search(r'(?i)\bnon[- ]?hemolytic\b|\bpoor hemolytic\b',x):return 0
 return np.nan
h['toxic_label']=h.apply(tlab,axis=1)
hh=pd.read_csv(PUB/'hemolysis.csv');hh['sequence']=hh.sequence.map(clean)
def hlab(r):
 if pd.isna(r.sequence):return np.nan
 if bool(r.hc50_censored) and r.hc50_uM>=128:return 0
 if not bool(r.hc50_censored) and r.hc50_uM<=64:return 1
 if not bool(r.hc50_censored) and r.hc50_uM>=128:return 0
 return np.nan
hh['toxic_label']=hh.apply(hlab,axis=1)
t=pd.concat([h[['sequence','toxic_label']],hh[['sequence','toxic_label']]],ignore_index=True).dropna();nu=t.groupby('sequence').toxic_label.nunique();bad=set(nu[nu>1].index);t=t[~t.sequence.isin(bad)&~t.sequence.isin(jseq)].drop_duplicates('sequence').reset_index(drop=True)
trt,tet=split(t,'toxic_label');Xt=Xmat(t.sequence);yt=t.toxic_label.astype(int).values
Xj=Xmat(ja.sequence);yj=ja.toxic_label.values
Xtrt=np.vstack([Xt[trt],Xj]);ytrt=np.r_[yt[trt],yj];wt=np.r_[np.ones(len(trt)),np.full(len(ja),20.0)]
tm=RandomForestClassifier(n_estimators=900,max_depth=14,min_samples_leaf=2,max_features='sqrt',class_weight='balanced_subsample',n_jobs=-1,random_state=42)
tm.fit(Xtrt,ytrt,sample_weight=wt);pt=tm.predict_proba(Xt[tet])[:,1]
toxm={'public_holdout_auc':roc_auc_score(yt[tet],pt),'public_holdout_ap':average_precision_score(yt[tet],pt),'public_holdout_bal_acc':balanced_accuracy_score(yt[tet],pt>=.5)}
joblib.dump(tm,MODEL/'toxicity_jnp_adapted.joblib')

# Compare public-only vs adapted on JNP
pubA=joblib.load(MODEL/'amp_classifier_public.joblib');pubM=joblib.load(MODEL/'mic_ecoli_public.joblib');pubT=joblib.load(MODEL/'toxicity_public.joblib')
rows=[]
for _,r in ja.iterrows():
 xx=featurize(r.sequence).reshape(1,-1)
 rows.append({'name':r['name'],'sequence':r.sequence,'observed_active':r.amp_label,'observed_ecoli_mic_uM':r.ecoli_mic,'observed_toxic':r.toxic_label,
 'public_amp':pubA.predict_proba(xx)[0,1],'adapted_amp':am.predict_proba(xx)[0,1],
 'public_mic':10**pubM.predict(xx)[0],'adapted_mic':10**mm.predict(xx)[0],
 'public_toxicity':pubT.predict_proba(xx)[0,1],'adapted_toxicity':tm.predict_proba(xx)[0,1]})
comp=pd.DataFrame(rows);comp.to_csv(OUT/'jnp_public_vs_domain_adapted.csv',index=False)
metrics={'jnp_weight':20,'amp':ampm,'mic':micm,'toxicity':toxm};(OUT/'domain_adapted_public_holdout_metrics.json').write_text(json.dumps(metrics,indent=2))
print(json.dumps(metrics,indent=2));print(comp.round(3).to_string(index=False))
