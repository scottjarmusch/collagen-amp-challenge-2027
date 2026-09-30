from pathlib import Path
import sys, re, json, math, joblib
import numpy as np, pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.metrics import roc_auc_score, balanced_accuracy_score, mean_absolute_error

ROOT=Path('/mnt/data/collagen_amp_ml'); sys.path.insert(0,str(ROOT/'src'))
from collagen_amp_ml.features import featurize, charge, boman, hmoment, apv, alpha

MODEL=ROOT/'models'; OUT=ROOT/'outputs'
amp=joblib.load(MODEL/'amp_classifier_public.joblib')
mic=joblib.load(MODEL/'mic_ecoli_public.joblib')
tox=joblib.load(MODEL/'toxicity_public.joblib')
j=pd.read_csv(ROOT/'data'/'jnp_collagen_amp_ml_anchors.csv')

def parse_mic(x):
    s=str(x); cens=s.startswith('>'); s=s.lstrip('<>')
    try:return float(s),cens
    except:return np.nan,cens

def compact_features(s):
    x=featurize(s).reshape(1,-1)
    p_amp=amp.predict_proba(x)[0,1]
    log_mic=mic.predict(x)[0]
    p_tox=tox.predict_proba(x)[0,1]
    n=len(s); arom=sum(s.count(a) for a in 'FWY')/n; hyd=sum(s.count(a) for a in 'AILMFWVY')/n
    return [p_amp,log_mic,p_tox,charge(s),boman(s),hmoment(s),apv(s),alpha(s),arom,hyd]

X=np.asarray([compact_features(s) for s in j.sequence],float)
# Activity = active against at least one tested bacterium at challenge threshold <=16 uM.
yact=[]; ytox=[]; ylog=[]; okmic=[]
for _,r in j.iterrows():
    e,ec=parse_mic(r.mic_ecoli_uM); s,sc=parse_mic(r.mic_saureus_uM)
    yact.append(int(min(e,s)<=16))
    ytox.append(int(r.cytotox_label in ['intermediate','high','very_high']))
    ylog.append(np.log10(e) if np.isfinite(e) and not ec else np.nan)
yact=np.array(yact); ytox=np.array(ytox); ylog=np.array(ylog)

act_model=Pipeline([('z',StandardScaler()),('lr',LogisticRegression(C=.15,class_weight='balanced',solver='liblinear',random_state=42))])
tox_model=Pipeline([('z',StandardScaler()),('lr',LogisticRegression(C=.08,class_weight='balanced',solver='liblinear',random_state=42))])
mic_cal=Pipeline([('z',StandardScaler()),('ridge',Ridge(alpha=8.0))])

loo=LeaveOneOut()
act_cv=cross_val_predict(act_model,X,yact,cv=loo,method='predict_proba')[:,1]
tox_cv=cross_val_predict(tox_model,X,ytox,cv=loo,method='predict_proba')[:,1]
ok=np.isfinite(ylog)
mic_cv=cross_val_predict(mic_cal,X[ok],ylog[ok],cv=LeaveOneOut())
metrics={
 'n_jnp':len(j),
 'activity_positive':int(yact.sum()),
 'toxicity_positive':int(ytox.sum()),
 'loo_activity_auc':float(roc_auc_score(yact,act_cv)),
 'loo_activity_balanced_accuracy':float(balanced_accuracy_score(yact,act_cv>=.5)),
 'loo_toxicity_auc':float(roc_auc_score(ytox,tox_cv)),
 'loo_toxicity_balanced_accuracy':float(balanced_accuracy_score(ytox,tox_cv>=.5)),
 'loo_ecoli_mic_mae_log10':float(mean_absolute_error(ylog[ok],mic_cv)),
 'warning':'JNP adapter validation is leave-one-out on n=11 and is descriptive, not an independent estimate of generalization.'
}
act_model.fit(X,yact); tox_model.fit(X,ytox); mic_cal.fit(X[ok],ylog[ok])
joblib.dump(act_model,MODEL/'jnp_activity_domain_adapter.joblib')
joblib.dump(tox_model,MODEL/'jnp_toxicity_domain_adapter.joblib')
joblib.dump(mic_cal,MODEL/'jnp_mic_domain_adapter.joblib')

rows=[]
for i,r in j.iterrows():
    rows.append({
      'name':r['name'],'sequence':r.sequence,
      'observed_active':int(yact[i]),'observed_cytotox':r.cytotox_label,
      'public_amp_probability':X[i,0],'public_pred_ecoli_mic_uM':10**X[i,1],'public_toxicity_probability':X[i,2],
      'JNP_adapted_activity_probability':float(act_model.predict_proba(X[i:i+1])[0,1]),
      'JNP_adapted_pred_ecoli_mic_uM':float(10**mic_cal.predict(X[i:i+1])[0]),
      'JNP_adapted_toxicity_probability':float(tox_model.predict_proba(X[i:i+1])[0,1]),
      'LOO_activity_probability':act_cv[i], 'LOO_toxicity_probability':tox_cv[i]
    })
pd.DataFrame(rows).to_csv(OUT/'jnp_domain_adaptation_predictions.csv',index=False)
(OUT/'jnp_adapter_metrics.json').write_text(json.dumps(metrics,indent=2))
print(json.dumps(metrics,indent=2))
print(pd.DataFrame(rows)[['name','public_amp_probability','JNP_adapted_activity_probability','public_pred_ecoli_mic_uM','JNP_adapted_pred_ecoli_mic_uM','public_toxicity_probability','JNP_adapted_toxicity_probability']].round(3).to_string(index=False))
