"""Parameter search used to choose synapse gain, inhibition gain and adaptation (Methods 2.4).

  python simulation/sweep.py 0.2,3,3 0.3,3,8   # ws,inh,adapt triplets
"""
import sys, numpy as np, pandas as pd, brain
from brain import load, stabilise, Brain
W, n = load(); t = n.type.fillna('').values
G = lambda *ts: np.flatnonzero(np.isin(t, ts))
taste = {'sugar': G('LB3c','LB3d','LB3b'), 'bitter': G('LB1a','LB1b','LB1c','LB1d'), 'tastepeg': G('claw_tpGRN','dorsal_tpGRN')}
acc = G('MN9','GNG588','GNG232','GNG540','GNG550','GNG056'); avr = G('MDN','DNg28')
pns = np.flatnonzero(pd.Series(t).str.contains(r'_(?:adPN|lPN|vPN|ilPN|l2PN|lvPN|l2PNm|lPNm)', regex=True).values)
orns = ['ORN_DM1','ORN_DM2','ORN_VA2','ORN_V','ORN_DA2','ORN_DL5']
for ws, inh, ad in [tuple(map(float, a.split(','))) for a in sys.argv[1:]]:
    brain.W_SCALE, brain.INH_GAIN, brain.ADAPT_INC = ws, inh, ad
    b = Brain(stabilise(W, n), n, seed=0)
    out = []
    for k, idx in taste.items():
        r = np.zeros(b.n, np.float32); r[idx] = 120; b.reset(); c = b.run(r, 250); off = b.run(None, 100)
        out.append(f"{k}: acc {c[acc].sum()} avr {c[avr].sum()} act {int((c>0).sum())}/{int((off>0).sum())}")
    vecs, offs = [], []
    for o in orns:
        r = np.zeros(b.n, np.float32); r[t == o] = 120; b.reset(); c = b.run(r, 250); off = b.run(None, 100)
        vecs.append(c[pns].astype(float)); offs.append(int((off > 0).sum()))
    C = np.corrcoef(vecs); corr = C[~np.eye(len(orns), dtype=bool)].mean()
    print(f"ws={ws} inh={inh} adapt={ad} | {' | '.join(out)} | odor PNcorr {corr:.2f} off {offs}", flush=True)
