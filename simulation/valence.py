"""Calibrate a smell-valence readout against real fly behaviour, then test it on unseen odorants.

Knaden et al. 2012 measured attraction indices (AI) for 110 odorants. For each odorant the DoOR receptor
profile drives the connectome's ORNs; brain output neurons (descending, lateral horn, mushroom body output)
are recorded. A ridge readout maps brain output -> AI, scored by 10-fold cross-validation.
Controls: the same readout on ORN rates directly (no brain), and on a scrambled-wiring brain.
Writes results/valence_readout.npz and results/valence_cv.json.
"""

import json

import numpy as np
import pandas as pd

import brain
from brain import Brain, load, scramble, stabilise
from receptors import Olfaction, to_input_vector

STIM_MS = 300
REPEATS = 2


def output_sets(n):
    t = n.type.fillna('')
    return {
        'DN': np.flatnonzero((n.superclass == 'descending_neuron').values),
        'LHON': np.flatnonzero(t.str.match(r'^(LHPV|LHAV|LHAD|LHCENT|LHPD|LHMB)').values),
        'MBON': np.flatnonzero(t.str.startswith('MBON').values),
    }


def ridge_cv(X, y, lam, folds=10, seed=0):
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(y))
    pred = np.zeros(len(y))
    for f in range(folds):
        te = idx[f::folds]
        tr = np.setdiff1d(idx, te)
        mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-6
        Z = (X[tr] - mu) / sd
        w = np.linalg.solve(Z.T @ Z + lam * np.eye(Z.shape[1]), Z.T @ (y[tr] - y[tr].mean()))
        pred[te] = ((X[te] - mu) / sd) @ w + y[tr].mean()
    return float(np.corrcoef(pred, y)[0, 1]), pred


def fit_ridge(X, y, lam):
    mu, sd = X.mean(0), X.std(0) + 1e-6
    Z = (X - mu) / sd
    w = np.linalg.solve(Z.T @ Z + lam * np.eye(Z.shape[1]), Z.T @ (y - y.mean()))
    return dict(w=w, mu=mu, sd=sd, b=y.mean())


def main():
    W, n = load()
    ol = Olfaction(n)
    beh = json.load(open('research/behavior.json'))['knaden2012_valence']['all_110_odorants']
    rows = []
    for o in beh:
        prof = ol.receptor_profile(o['CAS'])
        if prof is not None and prof.sum() > 0:
            rows.append((o['odorant'], o['AI_median'], ol.orn_rates_from_profile(prof)))
    print(f'{len(rows)} of {len(beh)} Knaden odorants have DoOR profiles')
    y = np.array([r[1] for r in rows])
    orn_X = np.array([[r[2].get(o, 0.0) for o in ol.orn_types] for r in rows])
    outs = output_sets(n)
    sel = np.concatenate(list(outs.values()))
    results = {'n_odorants': len(rows)}
    feats = {}
    for label in ['fly', 'scrambled']:
        Ws = stabilise(W, n)
        if label == 'scrambled':
            Ws = scramble(Ws, seed=0)
        b = Brain(Ws, n, seed=1)
        X = np.zeros((len(rows), len(sel)), np.float32)
        for i, (name, ai, rates) in enumerate(rows):
            v = to_input_vector(b.n, n, orn_rates=rates)
            for r in range(REPEATS):
                b.reset()
                X[i] += b.run(v, STIM_MS)[sel] / REPEATS
            if i % 20 == 0:
                print(f'  {label} {i}/{len(rows)} {name}: output spikes {X[i].sum():.0f}', flush=True)
        feats[label] = X
        del b, Ws
    np.savez('results/valence_features.npz', fly=feats['fly'], scrambled=feats['scrambled'], orn=orn_X, y=y,
             names=np.array([r[0] for r in rows]), sel=sel, orn_types=np.array(ol.orn_types))
    for lam in [10, 100, 1000, 10000]:
        entry = {}
        for label, X in [('no_brain_ORN_rates', orn_X), ('fly_brain', np.log1p(feats['fly'])),
                         ('scrambled_brain', np.log1p(feats['scrambled']))]:
            keep = X.std(0) > 0
            if keep.sum() == 0:
                entry[label] = None
                continue
            r, _ = ridge_cv(X[:, keep], y, lam)
            entry[label] = round(r, 3)
        results[f'lambda_{lam}'] = entry
        print(lam, entry, flush=True)
    best = max([k for k in results if k.startswith('lambda_')],
               key=lambda k: results[k]['fly_brain'] if results[k]['fly_brain'] is not None else -1)
    lam = int(best.split('_')[1])
    Xf = np.log1p(feats['fly'])
    keep = Xf.std(0) > 0
    model = fit_ridge(Xf[:, keep], y, lam)
    np.savez('results/valence_readout.npz', keep=keep, sel=sel, lam=lam, **model)
    results['chosen_lambda'] = lam
    json.dump(results, open('results/valence_cv.json', 'w'), indent=2)


if __name__ == '__main__':
    main()
