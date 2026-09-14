"""Probe: odour onset latency (ORN->PN->LH->DN) and left/right antenna asymmetry -> lateralised output."""
import json
import numpy as np, pandas as pd
from brain import Brain, load, stabilise
from chemistry import composition
from physics import emission_mol_s, co2_emission_mol_s, ppm_at, PEAK_FACTOR
from receptors import Olfaction

W, n = load(); t = n.type.fillna('').values; side = n.side.fillna('').values
ol = Olfaction(n)
b = Brain(stabilise(W, n), n, seed=5)
pn = np.flatnonzero(pd.Series(t).str.contains(r'_(?:adPN|lPN|vPN|ilPN|l2PN|lvPN|l2PNm|lPNm)', regex=True).values)
lh = np.flatnonzero(pd.Series(t).str.match(r'^(LHPV|LHAV|LHAD|LHCENT|LHPD|LHMB)').values)
dn = np.flatnonzero((n.superclass == 'descending_neuron').values)
steer = {s: np.flatnonzero(np.isin(t, ['DNa01', 'DNa02']) & (side == s)) for s in 'LR'}
dnL, dnR = dn[side[dn] == 'L'], dn[side[dn] == 'R']
orn_mask = pd.Series(t).str.startswith('ORN_').values
print('ORN sides', pd.Series(side[orn_mask]).value_counts().to_dict(), 'steer', {k: len(v) for k, v in steer.items()})

c = composition('damaged', 7)
q = emission_mol_s(c['volatiles_mg_per_kg'], c['water_L_per_kg'], c['taste']['pH'])
ppm = {k: v * PEAK_FACTOR for k, v in ppm_at(q, 1.0).items()}
rates = ol.orn_rates(ppm, 0)

def vec(scaleL, scaleR):
    v = np.zeros(b.n, np.float32)
    for o, r in rates.items():
        idx = ol.orn_index[o]
        v[idx[side[idx] == 'L']] += r * scaleL
        v[idx[side[idx] == 'R']] += r * scaleR
        v[idx[(side[idx] != 'L') & (side[idx] != 'R')]] += r * (scaleL + scaleR) / 2
    return v

# latency: 50 ms silence then odour on
out = {}
rec = []
b.reset(); b.run(None, 50, record=rec); b.run(vec(1, 1), 250, record=rec)
def first(idx, thr=3):
    s = set(idx)
    cnt = np.array([sum(1 for i in r if i in s) for r in rec])
    cs = np.cumsum(cnt[50:])
    hit = np.flatnonzero(cs >= thr)
    return int(hit[0]) if hit.size else None
lat = {'PN': first(pn), 'LH': first(lh), 'DN': first(dn, 5)}
print('latency ms after onset', lat); out['latency_ms'] = lat

res = []
for name, (sl, sr) in {'left only': (1, 0), 'left 2:1': (1, .5), 'both': (1, 1), 'right 2:1': (.5, 1), 'right only': (0, 1)}.items():
    acc = np.zeros(4)
    for r in range(6):
        b.reset(); cnt = b.run(vec(sl, sr), 300)
        acc += [cnt[steer['L']].sum(), cnt[steer['R']].sum(), cnt[dnL].sum(), cnt[dnR].sum()]
    acc /= 6
    idx = (acc[0] - acc[1]) / (acc[0] + acc[1] + 1e-9)
    idx_dn = (acc[2] - acc[3]) / (acc[2] + acc[3] + 1e-9)
    print(f'{name:10s} DNa01/02 L {acc[0]:.1f} R {acc[1]:.1f} (L-R index {idx:+.2f}) | all DNs L {acc[2]:.0f} R {acc[3]:.0f} ({idx_dn:+.2f})', flush=True)
    res.append(dict(stim=name, steerL=acc[0], steerR=acc[1], steer_index=idx, dnL=acc[2], dnR=acc[3], dn_index=idx_dn))
out['bilateral'] = res
json.dump(out, open('results/probe_bilateral.json', 'w'), indent=1)
