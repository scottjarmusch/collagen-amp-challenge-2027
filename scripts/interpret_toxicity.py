"""Exact tree-path probability decomposition; not SHAP or causal attribution."""
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
from collagen_amp_challenge.features import featurize,feature_names

R=Path(__file__).resolve().parents[1]

def decompose(model,x):
    contributions=np.zeros(len(x));baseline=0.
    assert list(model.classes_)==[0,1]
    for estimator in model.estimators_:
        tree=estimator.tree_
        positive=tree.value[:,0,1]/tree.value[:,0,:].sum(axis=1)
        baseline+=positive[0]
        node=0
        while tree.children_left[node]!=tree.children_right[node]:
            feature=tree.feature[node]
            child=tree.children_left[node] if x[feature]<=tree.threshold[node] else tree.children_right[node]
            contributions[feature]+=positive[child]-positive[node]
            node=child
    return baseline/len(model.estimators_),contributions/len(model.estimators_)

def main():
    anchors=pd.read_csv(R/'data/jnp_collagen_amp_ml_anchors.csv')
    rows=[];summaries=[];global_rows=[]
    names=feature_names()
    for suffix,tag in [('public','public'),('jnp_adapted','adapted')]:
        model=joblib.load(R/f'checkpoint/toxicity_{suffix}.joblib');model.set_params(n_jobs=1)
        global_rows.extend({'feature':f,'model':tag,'impurity_importance':v} for f,v in zip(names,model.feature_importances_))
        for _,anchor in anchors.iterrows():
            x=featurize(anchor.sequence);base,values=decompose(model,x)
            prediction=float(model.predict_proba(x.reshape(1,-1))[0,1])
            assert abs(base+values.sum()-prediction)<1e-12
            summaries.append({'name':anchor['name'],'model':tag,'baseline':base,'prediction':prediction,'reconstruction_error':base+values.sum()-prediction})
            rows.extend({'name':anchor['name'],'model':tag,'feature':f,'feature_value':float(v),'path_contribution':float(c)} for f,v,c in zip(names,x,values))
    contributions=pd.DataFrame(rows)
    contributions.to_csv(R/'docs/toxicity_path_contributions.csv',index=False)
    summary=pd.DataFrame(summaries);summary.to_csv(R/'docs/toxicity_path_summary.csv',index=False)
    pd.DataFrame(global_rows).to_csv(R/'docs/toxicity_global_importance.csv',index=False)
    delta=contributions.pivot(index=['name','feature','feature_value'],columns='model',values='path_contribution').reset_index()
    delta['adapted_minus_public']=delta.adapted-delta.public
    delta.to_csv(R/'docs/toxicity_path_changes.csv',index=False)
    text=['# Toxicity model interpretation','',
        'This is an exact decision-tree-path decomposition of the two archived random forests.',
        'Each prediction equals the average root probability plus feature contributions along',
        'the selected paths. The residual is checked below 1e-12. This is path-dependent,',
        'not SHAP, a counterfactual intervention, or evidence of a biological mechanism.',
        'Correlated physicochemical/composition/motif features can trade attribution.',
        'The forests differ in training data and tree count (700 versus 900), so the',
        'comparison describes the fitted models, not an isolated controlled causal effect.',
        'JNP anchors were used for adaptation; this analysis is in-sample interpretation.','']
    wanted=['REI-26','TRR-26','LRS-21','GFD-30','TFK-18','LL-37','melittin']
    for name in wanted:
        sub=summary[summary.name==name].set_index('model')
        d=delta[delta.name==name].copy()
        largest=d.loc[d.adapted_minus_public.abs().nlargest(8).index]
        text += [f'## {name}','',f'Probability {sub.loc["public","prediction"]:.6f} -> {sub.loc["adapted","prediction"]:.6f}; '
                 f'baseline change {sub.loc["adapted","baseline"]-sub.loc["public","baseline"]:+.6f}.','',
                 '| Feature | Public path contribution | Adapted path contribution | Change |',
                 '|---|---:|---:|---:|']
        text += [f'| {r.feature} | {r.public:+.6f} | {r.adapted:+.6f} | {r.adapted_minus_public:+.6f} |' for r in largest.itertuples()]
        text += ['','Negative changes lower the fitted toxicity probability; positive changes raise it.',
                 'Only the largest changes are shown; the complete CSV includes all 432 features.','']
    text += ['TFK-18 refers to the preserved 33-residue archive record; see the sequence caveat in JNP_DOMAIN_ADAPTATION.md.']
    (R/'docs/TOXICITY_INTERPRETATION.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
    print('Exact path decompositions verified for all 11 anchors in both archived models.')

if __name__=='__main__':main()
