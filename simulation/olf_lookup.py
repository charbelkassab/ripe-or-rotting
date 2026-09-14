"""Precompute the connectome's olfactory response to each food's odour at every plume dilution.

For food odour mixtures (from chemistry + emission physics) scaled by a unit-emission concentration c_unit
(s/m3, the plume variable), the whole CNS is run for 100 ms (the decision window) with spontaneous ORN
firing (8 Hz) and the internal-state gains. Stored: projection-neuron (PN) spikes in total and per glomerulus,
descending-neuron spikes, plus a no-odour baseline -> results/olf_lookup.npz
"""
import numpy as np, pandas as pd

from brain import Brain, load, stabilise
from chemistry import composition
from physics import R_GAS, co2_emission_mol_s, emission_mol_s
from receptors import Olfaction, to_input_vector
from states import STATES, central_drive

WINDOW_MS = 100
SPONT_ORN = 8.0
LEVELS = np.logspace(-1, 5, 13)
REPEATS = 4
BASE_REPEATS = 16
CONDITIONS = [  # (food key, scenario, day, temp_C, states)
    ('ripe', 'intact', 0, 25, ['fed', 'starved_24h']),
    ('rotting', 'damaged', 7, 25, ['fed', 'starved_24h']),
    ('ripe', 'intact', 0, 18, ['fed']), ('rotting', 'damaged', 7, 18, ['fed']),
    ('ripe', 'intact', 0, 32, ['fed']), ('rotting', 'damaged', 7, 32, ['fed']),
]


def main():
    W, n = load()
    ol = Olfaction(n)
    b = Brain(stabilise(W, n), n, seed=11)
    t = n.type.fillna('').values
    pn_glom = {}
    for i, ty in enumerate(t):
        if 'PN' in ty and '_' in ty and f'ORN_{ty.split("_")[0]}' in ol.orn_index:
            pn_glom.setdefault(ty.split('_')[0], []).append(i)
    gloms = sorted(pn_glom)
    pn_all = np.concatenate([np.array(v) for v in pn_glom.values()])
    dn = np.flatnonzero((n.superclass == 'descending_neuron').values)
    spont = {o: SPONT_ORN for o in ol.orn_types}

    def run(rates, state, reps):
        v = to_input_vector(b.n, n, orn_rates=rates) + central_drive(n, state)
        out = []
        for _ in range(reps):
            b.reset()
            c = b.run(v, WINDOW_MS)
            out.append([c[pn_all].sum(), c[dn].sum()] + [c[pn_glom[g]].sum() for g in gloms])
        return np.array(out, float)

    save = {'levels': LEVELS, 'gloms': np.array(gloms)}
    for state in ['fed', 'starved_24h']:
        rates = {o: r * STATES[state]['orn'].get(o, 1.0) for o, r in spont.items()}
        base = run(rates, state, BASE_REPEATS)
        save[f'baseline_{state}'] = base
        print(state, 'baseline PN', base[:, 0].mean().round(1), '+/-', base[:, 0].std().round(1), flush=True)
    for food, sc, day, T, states in CONDITIONS:
        comp = composition(sc, day, temp_C=T)
        q = emission_mol_s(comp['volatiles_mg_per_kg'], comp['water_L_per_kg'], comp['taste']['pH'], T)
        qco2 = co2_emission_mol_s(comp['extra']['fermentation_rate_g_per_kg_h'])
        save[f'q_{food}_{T}'] = np.array([q.get(k, 0.0) for k in sorted(q)])
        save[f'qnames_{food}_{T}'] = np.array(sorted(q))
        save[f'qco2_{food}_{T}'] = qco2
        for state in states:
            res = np.zeros((len(LEVELS), REPEATS, 2 + len(gloms)))
            for li, c1 in enumerate(LEVELS):
                ppm = {k: v * c1 * R_GAS * 1e6 for k, v in q.items()}
                rates = ol.orn_rates(ppm, qco2 * c1 * R_GAS * 1e6, gains=STATES[state]['orn'])
                for o in ol.orn_types:
                    rates[o] = rates.get(o, 0.0) + SPONT_ORN * STATES[state]['orn'].get(o, 1.0)
                res[li] = run(rates, state, REPEATS)
            save[f'resp_{food}_{T}_{state}'] = res
            print(food, T, state, 'PN by level', res[:, :, 0].mean(1).round(0).tolist(), flush=True)
    np.savez('results/olf_lookup.npz', **save)


if __name__ == '__main__':
    main()
