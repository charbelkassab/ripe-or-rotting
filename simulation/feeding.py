"""Taste -> feeding decision readout, and validation against known fly feeding behaviour.

Acceptance neurons: MN9 (proboscis extension motor neuron), Fdg (GNG588, feeding command),
G2N-1 (GNG232), sugar SEL PN/LN (GNG540, GNG550, GNG056; Yao & Scott 2022).
Aversion neurons: MDN (moonwalker, backing away), DNg28 (Bitter-SEL; Yao & Scott 2022).
Validation (results/feeding_validation.json):
  1. sucrose dose-response: acceptance should rise with concentration
  2. bitter added to sucrose should lower acceptance
  3. starvation should raise acceptance of the same sucrose
"""

import json

import numpy as np

from brain import Brain, load, stabilise
from receptors import GRN_TYPES, taste_rates, to_input_vector
from states import STATES, central_drive

STIM_MS = 400
REPEATS = 6
ACCEPT = ['MN9', 'GNG588', 'GNG232', 'GNG540', 'GNG550', 'GNG056']
AVERSE = ['MDN', 'DNg28']


def solution(sucrose_mM=0.0, caffeine_like_bitter=0.0):
    """A plain test solution in the same format chemistry.composition() produces."""
    return {'sucrose_mM': sucrose_mM, 'glucose_mM': 0.0, 'fructose_mM': 0.0, 'ethanol_mM': 0.0,
            'acetic_acid_mM': 0.0, 'pH': 6.5, 'K_mM': 0.0, 'Na_mM': 0.0, 'tannin_mM': caffeine_like_bitter,
            'CO2_mM': 0.0, 'yeast_percent_w_v': 0.0, 'amino_acids_mM': 0.0, 'glycerol_mM': 0.0,
            'osmolality_mOsm': sucrose_mM + 50}


class Taster:
    def __init__(self, W=None, n=None):
        if W is None:
            W, n = load()
        self.n = n
        t = n.type.fillna('').values
        self.acc = np.flatnonzero(np.isin(t, ACCEPT))
        self.avr = np.flatnonzero(np.isin(t, AVERSE))
        self.brain = Brain(stabilise(W, n), n, seed=2)

    def taste(self, taste_dict, state='fed', extra_input=None, repeats=REPEATS):
        g = STATES[state]['grn']
        rates = taste_rates(taste_dict, gains=g)
        v = to_input_vector(self.brain.n, self.n, grn_rates=rates) + central_drive(self.n, state)
        if extra_input is not None:
            v = v + extra_input
        acc = avr = 0.0
        for _ in range(repeats):
            self.brain.reset()
            c = self.brain.run(v, STIM_MS)
            acc += c[self.acc].sum() / repeats
            avr += c[self.avr].sum() / repeats
        return {'accept': acc, 'aversion': avr, 'score': acc - avr, 'grn_rates': rates}


def main():
    T = Taster()
    out = {}
    curve = {}
    for s in [0, 10, 30, 100, 300, 1000]:
        r = T.taste(solution(s))
        curve[s] = r
        print('sucrose', s, round(r['accept'], 1), round(r['aversion'], 1), flush=True)
    vals = [curve[s]['score'] for s in sorted(curve)]
    out['sucrose_dose_response'] = {str(s): {k: round(v, 2) for k, v in curve[s].items() if k != 'grn_rates'} for s in curve}
    out['test1_monotonic_rise'] = bool(all(b >= a - 1.0 for a, b in zip(vals, vals[1:])) and vals[-1] > vals[0] + 2)
    base = T.taste(solution(100))
    bitter = T.taste(solution(100, caffeine_like_bitter=2.0))
    out['sucrose100_plus_bitter'] = {'without': round(base['score'], 2), 'with': round(bitter['score'], 2)}
    out['test2_bitter_suppresses'] = bool(bitter['score'] < base['score'] - 1.0)
    starved = T.taste(solution(30), state='starved_24h')
    fed = T.taste(solution(30), state='fed')
    out['sucrose30_fed_vs_starved'] = {'fed': round(fed['score'], 2), 'starved': round(starved['score'], 2)}
    out['test3_starvation_increases'] = bool(starved['score'] > fed['score'] + 1.0)
    print(json.dumps({k: v for k, v in out.items() if k.startswith('test') or 'vs' in k or 'bitter' in k}, indent=1))
    json.dump(out, open('results/feeding_validation.json', 'w'), indent=2)


if __name__ == '__main__':
    main()
