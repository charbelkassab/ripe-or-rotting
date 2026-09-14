"""Calibrate the acceptance -> feeding probability logistic on proboscis extension data, then test it.

Calibration (Inagaki et al. 2014): fed flies extend the proboscis 50-60% of the time to 800 mM sucrose.
Held-out test (same paper): 1-day starved flies reach 50-60% at 300 mM. The model predicts the starved
half-response concentration from its own starved-state acceptance curve.
Slope: set so that a fed fly goes from ~10% to ~90% over the model's 100 -> 1000 mM acceptance range (ESTIMATE).
"""

import json

import numpy as np

from feeding import Taster, solution

CONC = [10, 30, 100, 300, 800, 1000, 2000]


def main():
    T = Taster()
    curves = {}
    for state in ['fed', 'starved_24h']:
        curves[state] = [T.taste(solution(c), state=state, repeats=8)['accept'] for c in CONC]
        print(state, [round(a, 1) for a in curves[state]], flush=True)
    fed = np.array(curves['fed'])
    a800 = float(np.interp(np.log10(800), np.log10(CONC), fed))
    a100 = float(np.interp(2, np.log10(CONC), fed))
    a1000 = float(np.interp(3, np.log10(CONC), fed))
    slope = max((a1000 - a100) / (2 * np.log(9)), 0.3)
    A50 = a800 - slope * np.log(0.55 / 0.45)  # 55% at 800 mM
    st = np.array(curves['starved_24h'])
    p_st = 1 / (1 + np.exp(-(st - A50) / slope))
    p_fed = 1 / (1 + np.exp(-(fed - A50) / slope))
    above = np.flatnonzero(p_st >= 0.55)
    if above.size == 0:
        pred = None
    elif above[0] == 0:
        pred = CONC[0]
    else:
        i = above[0]
        pred = float(10 ** np.interp(0.55, [p_st[i - 1], p_st[i]], np.log10([CONC[i - 1], CONC[i]])))
    out = dict(A50=A50, slope=slope, conc_mM=CONC, accept_fed=curves['fed'], accept_starved=curves['starved_24h'],
               p_fed=p_fed.tolist(), p_starved=p_st.tolist(), predicted_starved_55pct_mM=pred,
               observed_starved_1d_55pct_mM=300, observed_fed_55pct_mM=800,
               note='fed point is the calibration target; starved point is a held-out prediction')
    print(json.dumps({k: v for k, v in out.items() if not isinstance(v, list)}, indent=1))
    json.dump(out, open('results/feeding_calibration.json', 'w'), indent=2)


if __name__ == '__main__':
    main()
