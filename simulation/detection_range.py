"""How far away can the simulated brain detect each banana? -> results/detection_range.json

Detection threshold c* (s/m3) = plume concentration per unit emission at which the connectome olfactory lookup
reaches 50% detection in a 100 ms window. Distance where the plume centreline reaches c*:
  mean:      1/(pi U sigma^2)            sigma = a r + 0.005
  filaments: PEAK x mean, PEAK = p95/mean measured in atmosphere.py (~6)
Spread coefficient a = 0.1 (calibrated indoor plume) and 0.2 (more turbulent outdoor air, ESTIMATE).
"""
import json
import numpy as np
from navigate import Sensory

exp = json.load(open('results/experiments.json'))
peak = float(np.mean([c['p95_over_mean'] for c in exp['plume_calibration']]))
sens = Sensory()
lv = sens.levels
out = {'peak_factor_p95': peak, 'rows': []}
for food, temp, state in [('ripe', 25, 'fed'), ('ripe', 25, 'starved_24h'), ('rotting', 25, 'fed'), ('rotting', 25, 'starved_24h'),
                          ('ripe', 18, 'fed'), ('rotting', 18, 'fed'), ('ripe', 32, 'fed'), ('rotting', 32, 'fed')]:
    p, _, _ = sens.table(food, temp, state)
    i = int(np.argmax(p >= 0.5))
    c_star = float(10 ** np.interp(0.5, [p[i - 1], p[i]], [lv[i - 1], lv[i]])) if i > 0 else float(10 ** lv[0])
    for U in [0.3, 1.0]:
        for a in [0.1, 0.2]:
            for label, factor in [('mean', 1.0), ('filament_p95', peak)]:
                sigma = np.sqrt(factor / (np.pi * U * c_star))
                r = max((sigma - 0.005) / a, 0.0)
                out['rows'].append(dict(food=food, temp=temp, state=state, c_star=c_star, wind=U, spread=a,
                                        basis=label, range_m=float(r)))
    print(food, temp, state, 'c* = %.2f s/m3' % c_star,
          {f"U{row['wind']}_a{row['spread']}_{row['basis']}": round(row['range_m'], 1) for row in out['rows'][-8:]})
json.dump(out, open('results/detection_range.json', 'w'), indent=1)
