
from pathlib import Path
import json, joblib, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, ExtraTreesRegressor, RandomForestClassifier
from sklearn.metrics import roc_auc_score, average_precision_score, mean_absolute_error
from sklearn.model_selection import GroupShuffleSplit
from .features import featurize, feature_names
from .data import load_binary, jnp_activity_table

def seq_cluster_key(s,k=4):
    # deterministic coarse sequence family key used only as a leakage-resistant fallback split.
    # Final release should replace this with MMseqs2/CD-HIT clustering if available.
    return min([s[i:i+k] for i in range(max(1,len(s)-k+1))])

def matrix(seqs):
    return np.vstack([featurize(s) for s in seqs])

def train_amp_classifier(pos_csv,neg_csv,outdir):
    pos=load_binary(pos_csv,1); neg=load_binary(neg_csv,0)
    df=pd.concat([pos,neg],ignore_index=True).drop_duplicates("sequence")
    X=matrix(df.sequence); y=df.amp_label.values
    groups=np.array([seq_cluster_key(s) for s in df.sequence])
    sp=GroupShuffleSplit(n_splits=1,test_size=.2,random_state=42)
    tr,te=next(sp.split(X,y,groups))
    model=HistGradientBoostingClassifier(max_iter=300,learning_rate=.05,max_leaf_nodes=31,l2_regularization=1.0,random_state=42)
    model.fit(X[tr],y[tr])
    p=model.predict_proba(X[te])[:,1]
    metrics={"roc_auc":float(roc_auc_score(y[te],p)),"average_precision":float(average_precision_score(y[te],p)),
             "n_train":int(len(tr)),"n_test":int(len(te))}
    Path(outdir).mkdir(parents=True,exist_ok=True)
    joblib.dump(model,Path(outdir)/"amp_classifier.joblib")
    (Path(outdir)/"amp_classifier_metrics.json").write_text(json.dumps(metrics,indent=2))
    return metrics

def fit_jnp_domain_adapter(jnp_csv,outdir):
    # Small-domain adapter, intentionally not a standalone global predictor.
    d=jnp_activity_table(jnp_csv)
    cyto_map={"low":0,"low_or_neutral":0,"intermediate":1,"high":2,"very_high":3}
    d["cyto_ord"]=d.cytotox_label.map(cyto_map)
    X=matrix(d.sequence)
    ok=np.isfinite(d.mic_bacterial_geomean_uM)
    # Tiny JNP set: use a deliberately high-variance local ensemble only as a domain
    # calibration/contrast layer. It is never the final global MIC model.
    activity=ExtraTreesRegressor(n_estimators=500,max_features=0.7,min_samples_leaf=1,random_state=42)
    activity.fit(X[ok],np.log2(d.loc[ok,"mic_bacterial_geomean_uM"].values))
    toxicity=RandomForestClassifier(n_estimators=300,max_depth=4,min_samples_leaf=1,class_weight="balanced",random_state=42)
    toxicity.fit(X,d.cyto_ord.values)
    Path(outdir).mkdir(parents=True,exist_ok=True)
    joblib.dump(activity,Path(outdir)/"jnp_activity_adapter.joblib")
    joblib.dump(toxicity,Path(outdir)/"jnp_toxicity_adapter.joblib")
    return {"n_jnp":len(d),"note":"Domain adapters are calibration/contrast models only; full public models required for final inference."}
