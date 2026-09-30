"""Reconstruct archived deterministic splits, freeze them, and train from CSVs.

The original scripts are retained unchanged. Only their preprocessing sections
are executed here; archived hard-coded paths and fitting/output code are not.
Released inference checkpoints are never overwritten.
"""
import os
os.environ['OMP_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
from pathlib import Path
import ast, hashlib, json
import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor, RandomForestClassifier
from sklearn.metrics import roc_auc_score, average_precision_score, balanced_accuracy_score, mean_absolute_error, mean_squared_error, r2_score
from collagen_amp_challenge.features import featurize, canonical

R=Path(__file__).resolve().parents[1]
SPLITS=R/'data/splits'
OUT=R/'training_reproduction'

def between(source,start,end):
    return source[source.index(start):source.index(end)]

def namespace(source):
    # Load exactly the archived imports and helper function definitions.
    nodes=[]
    for node in ast.parse(source).body:
        if isinstance(node,(ast.Import,ast.ImportFrom,ast.FunctionDef)):
            if isinstance(node,ast.ImportFrom) and node.module=='collagen_amp_ml.features':continue
            nodes.append(node)
    ns={'ROOT':R,'PUB':R/'data/training','featurize':featurize,'canonical':canonical}
    exec(compile(ast.Module(body=nodes,type_ignores=[]),'archived_helpers','exec'),ns)
    return ns

def reconstruct(adapted=False):
    name='train_domain_adapted.py' if adapted else 'train_full_public.py'
    source=(R/'scripts/archived_training'/name).read_text(encoding='utf-8')
    ns=namespace(source)
    if adapted:
        exec(between(source,"j=pd.read_csv",'# AMP public peptide benchmark'),ns)
        exec(between(source,'parts=[]','Xtr=np.vstack'),ns)
        exec(between(source,"m=pd.read_csv(PUB/'mic_data.csv')",'jm=ja.dropna'),ns)
        exec(between(source,"h=pd.read_csv(PUB/'Hemolytik2_complete_data.csv')",'Xj=Xmat'),ns)
        return {'amp':(ns['a'],ns['tr'],ns['te'],'amp_label'),
                'mic':(ns['m'],ns['trm'],ns['tem'],'value'),
                'toxicity':(ns['t'],ns['trt'],ns['tet'],'toxic_label')},ns
    exec(between(source,'parts=[]','X=Xmat(amp_bal.sequence)'),ns)
    exec(between(source,"mic=pd.read_csv(PUB/'mic_data.csv')",'Xm=Xmat'),ns)
    exec(between(source,"h=pd.read_csv(PUB/'Hemolytik2_complete_data.csv')",'Xt=Xmat'),ns)
    return {'amp':(ns['amp_bal'],ns['tr'],ns['te'],'amp_label'),
            'mic':(ns['mic'],ns['trm'],ns['tem'],'value'),
            'toxicity':(ns['tox'],ns['trt'],ns['tet'],'toxic_label')},ns

def group(s):
    return f'{len(s)//5}:{min([s[i:i+4] for i in range(len(s)-3)] or [s])}'

def freeze(task,df,tr,te,label,prefix=''):
    df=df.copy().reset_index(drop=True)
    df['reconstructed_row_index']=np.arange(len(df))
    df['sequence_id']=[hashlib.sha256(s.encode()).hexdigest() for s in df.sequence]
    df['split_group']=[group(s) for s in df.sequence]
    assert not(set(df.iloc[tr].sequence)&set(df.iloc[te].sequence))
    assert not(set(df.iloc[tr].split_group)&set(df.iloc[te].split_group))
    for role,ix in [('train',tr),('holdout',te)]:
        p=SPLITS/f'{prefix}{task}_{role}.csv'
        data=df.iloc[ix].to_csv(index=False,lineterminator='\n').encode()
        if p.exists():
            assert p.read_bytes()==data,f'Frozen split changed: {p}; investigate rather than overwrite.'
        else:p.write_bytes(data)

def metrics(model,frame,label,regression=False):
    X=np.vstack([featurize(s) for s in frame.sequence]);y=frame[label].to_numpy()
    if regression:
        p=model.predict(X)
        return {'mae_log10_mic':float(mean_absolute_error(y,p)),
                'rmse_log10_mic':float(mean_squared_error(y,p)**.5),'r2':float(r2_score(y,p))}
    p=model.predict_proba(X)[:,1]
    return {'roc_auc':float(roc_auc_score(y,p)),
            'average_precision':float(average_precision_score(y,p)),
            'balanced_accuracy':float(balanced_accuracy_score(y,p>=.5))}

def main():
    SPLITS.mkdir(exist_ok=True);OUT.mkdir(exist_ok=True)
    public,ns=reconstruct();adapted,ans=reconstruct(True)
    for prefix,tasks in [('',public),('adapted_public_',adapted)]:
        for task,args in tasks.items():freeze(task,*args,prefix=prefix)
    anchors=ans['ja'].copy()
    anchors['adaptation_weight']=20.0
    anchors.to_csv(SPLITS/'jnp_adaptation.csv',index=False,lineterminator='\n')
    estimators={
        'amp':HistGradientBoostingClassifier(max_iter=350,learning_rate=.05,max_leaf_nodes=31,l2_regularization=1.,random_state=42),
        'mic':HistGradientBoostingRegressor(max_iter=400,learning_rate=.04,max_leaf_nodes=23,l2_regularization=2.,random_state=42),
        'toxicity':RandomForestClassifier(n_estimators=700,max_depth=14,min_samples_leaf=2,max_features='sqrt',class_weight='balanced_subsample',n_jobs=1,random_state=42),
    }
    names={'amp':'amp_classifier','mic':'mic_ecoli','toxicity':'toxicity'}
    expected=json.loads((R/'docs/public_model_metrics.json').read_text())
    adapted_expected=json.loads((R/'docs/domain_adapted_public_holdout_metrics.json').read_text())
    results={'split_origin':'Reconstructed from archived code and seed 42; not separately preserved original row IDs.',
             'public_retraining':{},'adapted_archived_holdouts':{},'shared_untouched_holdouts':{},'jnp_public_overlap':{}}
    for task,estimator in estimators.items():
        label=public[task][3]
        train=pd.read_csv(SPLITS/f'{task}_train.csv',float_precision='round_trip')
        test=pd.read_csv(SPLITS/f'{task}_holdout.csv',float_precision='round_trip')
        X=np.vstack([featurize(s) for s in train.sequence]);y=train[label].to_numpy()
        print(f'Retraining public {task}: {len(train)} train / {len(test)} holdout',flush=True)
        estimator.fit(X,y)
        got=metrics(estimator,test,label,task=='mic')
        archive=expected[names[task]]
        deltas={k:got[k]-archive[k] for k in got}
        assert all(abs(d)<.005 for d in deltas.values()),(task,deltas)
        archived_model=joblib.load(R/f'checkpoint/{names[task]}_public.joblib')
        if task=='toxicity':archived_model.set_params(n_jobs=1)
        archived_metrics=metrics(archived_model,test,label,task=='mic')
        Xt=np.vstack([featurize(s) for s in test.sequence])
        p=estimator.predict(Xt) if task=='mic' else estimator.predict_proba(Xt)[:,1]
        q=archived_model.predict(Xt) if task=='mic' else archived_model.predict_proba(Xt)[:,1]
        results['public_retraining'][task]={'n_train':len(train),'n_holdout':len(test),'metrics':got,
            'difference_from_archived_metrics':deltas,'archived_model_on_reconstructed_holdout':archived_metrics,
            'max_prediction_difference_from_archived_model':float(np.max(np.abs(p-q)))}
        joblib.dump(estimator,OUT/f'{names[task]}_retrained.joblib')
        am=joblib.load(R/f'checkpoint/{names[task]}_jnp_adapted.joblib')
        if task=='toxicity':am.set_params(n_jobs=1)
        atest=pd.read_csv(SPLITS/f'adapted_public_{task}_holdout.csv',float_precision='round_trip')
        atrain=pd.read_csv(SPLITS/f'adapted_public_{task}_train.csv',float_precision='round_trip')
        adapted_got=metrics(am,atest,label,task=='mic')
        key_map={'roc_auc':'public_holdout_auc','average_precision':'public_holdout_ap','balanced_accuracy':'public_holdout_bal_acc',
                 'mae_log10_mic':'public_holdout_mae_log10','rmse_log10_mic':'public_holdout_rmse_log10','r2':'public_holdout_r2'}
        expected_task=adapted_expected['toxicity' if task=='toxicity' else task]
        assert all(abs(v-expected_task[key_map[k]])<1e-12 for k,v in adapted_got.items())
        anchor_groups=set(group(s) for s in anchors.sequence)
        assert not(set(atest.split_group)&anchor_groups)
        results['adapted_archived_holdouts'][task]={'n_holdout':len(atest),'metrics':adapted_got,'jnp_group_overlap':0}
        overlap=set(test.sequence)&set(atest.sequence)
        common=test[test.sequence.isin(overlap)].copy()
        assert not(overlap & (set(train.sequence)|set(atrain.sequence)|set(anchors.sequence)))
        common.to_csv(SPLITS/f'{task}_shared_untouched_holdout.csv',index=False,lineterminator='\n')
        results['shared_untouched_holdouts'][task]={'n_common':len(common),
            'same_original_holdout_sequences':set(test.sequence)==set(atest.sequence),
            'public':metrics(archived_model,common,label,task=='mic'),
            'adapted':metrics(am,common,label,task=='mic')}
        results['jnp_public_overlap'][task]={'public_train_names':anchors[anchors.sequence.isin(train.sequence)]['name'].tolist(),
            'public_holdout_names':anchors[anchors.sequence.isin(test.sequence)]['name'].tolist()}
    results['removed_conflicts']={'amp':len(ns['bad']),'toxicity':len(ns['badtox'])}
    results['split_files_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(SPLITS.glob('*.csv'))}
    (R/'docs/TRAINING_REPRODUCTION.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(results,indent=2))

if __name__=='__main__':main()
