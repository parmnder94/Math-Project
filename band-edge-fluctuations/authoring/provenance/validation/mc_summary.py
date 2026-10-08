# Summarize Monte Carlo of X_N (mc_lor.py outputs) against the exact constants.
import numpy as np, glob, json, sys
files=sys.argv[1:-2];truth=json.load(open(sys.argv[-2]));out=sys.argv[-1]
xs=np.concatenate([np.load(f) for f in files]);S=len(xs);N=200
mean=xs.mean();var=xs.var(ddof=1);c3=np.mean((xs-mean)**3)
# bootstrap standard errors
rng=np.random.default_rng(0);B=400;bm=[];bv=[];bk=[]
for _ in range(B):
    y=xs[rng.integers(0,S,S)];mu=y.mean();bm.append(mu);bv.append(y.var(ddof=1));bk.append(np.mean((y-mu)**3))
pred_mean=2*N*truth["z_inf"]+truth["mean_correction"]
res={"N":N,"samples":int(S),
     "E[X_N]":{"mc":mean,"se":float(np.std(bm)),"theory_2Nz+m":pred_mean,"z_score":(mean-pred_mean)/float(np.std(bm))},
     "m_estimate":{"mc":mean-2*N*truth["z_inf"],"se":float(np.std(bm)),"theory":truth["mean_correction"]},
     "Var":{"mc":var,"se":float(np.std(bv)),"theory":truth["variance"],"z_score":(var-truth["variance"])/float(np.std(bv))},
     "third_cumulant":{"mc":c3,"se":float(np.std(bk)),"theory":truth["third_cumulant"],"z_score":(c3-truth["third_cumulant"])/float(np.std(bk))},
     "note":"finite-N values carry O(1/N) corrections in addition to the statistical errors"}
json.dump(res,open(out,"w"),indent=2);print(json.dumps(res,indent=2))
