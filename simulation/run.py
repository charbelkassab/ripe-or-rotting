"""Full grid: food stage x internal state, for taste on contact and smell at distance.

  python run.py taste   -> results/grid_taste.json
  python run.py smell   -> results/grid_smell.json
"""

import json
import sys

import numpy as np

from brain import Brain, load, stabilise
from chemistry import composition
from feeding import Taster
from physics import PEAK_FACTOR, co2_emission_mol_s, emission_mol_s, ppm_at
from receptors import Olfaction, to_input_vector
from states import STATES, central_drive

FOODS = [('ripe (day 0)', 'intact', 0), ('overripe (day 5)', 'intact', 5), ('rotting day 2', 'damaged', 2),
         ('rotting day 4', 'damaged', 4), ('rotting day 7', 'damaged', 7), ('rotting day 14', 'damaged', 14)]
STATE_LIST = ['fed', 'starved_24h', 'protein_deprived', 'starved_and_protein_deprived']
DISTANCES = [0.01, 0.1, 1.0, 10.0]
YEAST_SURFACE_ENRICHMENT = [1.0, 20.0]  # bulk average vs flies feeding on surface colonies (ESTIMATE)

# Literature glomerular valence, used only as a labelled estimate (our own readout failed validation).
VALENCE = {  # +1 attractive, -1 aversive
    'DM1': 1, 'VA2': 1,                       # Semmelhack & Wang 2009 (vinegar attraction)
    'DM2': 1, 'DM4': 1,                       # Knaden 2012 PN-level attractive bias
    'D': -1, 'DA4l': -1, 'DA4m': -1, 'DC3': -1, 'DL1': -1, 'DL4': -1, 'DL5': -1,  # Knaden 2012 aversive bias
    'V': -1,                                  # CO2 avoidance at rest (Suh 2004)
    'DA2': -1,                                # geosmin (Stensmyr 2012)
    'DC4': -1,                                # acid avoidance via Ir64a (Ai 2010)
}


def taste_grid():
    T = Taster()
    out = []
    for label, sc, day in FOODS:
        comp = composition(sc, day)
        for enrich in YEAST_SURFACE_ENRICHMENT:
            if sc == 'intact' and enrich != 1.0:
                continue
            t = dict(comp['taste'])
            t['yeast_percent_w_v'] *= enrich
            for state in STATE_LIST:
                r = T.taste(t, state=state)
                row = dict(food=label, yeast_enrichment=enrich, state=state, accept=r['accept'],
                           aversion=r['aversion'], score=r['score'],
                           grn_rates={k: round(v, 1) for k, v in r['grn_rates'].items()})
                out.append(row)
                print(f"{label:17s} yeast x{enrich:<4g} {state:28s} accept {r['accept']:5.1f} aversion {r['aversion']:4.1f} "
                      f"| GRN {row['grn_rates']}", flush=True)
    json.dump(out, open('results/grid_taste.json', 'w'), indent=1)


def smell_grid():
    W, n = load()
    ol = Olfaction(n)
    b = Brain(stabilise(W, n), n, seed=3)
    t = n.type.fillna('').values
    pn_glom = {}
    for i, ty in enumerate(t):
        if 'PN' in ty and '_' in ty:
            g = ty.split('_')[0]
            if f'ORN_{g}' in ol.orn_index:
                pn_glom.setdefault(g, []).append(i)
    pn_glom = {g: np.array(v) for g, v in pn_glom.items()}
    out = []
    for label, sc, day in FOODS:
        comp = composition(sc, day)
        q = emission_mol_s(comp['volatiles_mg_per_kg'], comp['water_L_per_kg'], comp['taste']['pH'])
        qco2 = co2_emission_mol_s(comp['extra']['fermentation_rate_g_per_kg_h'])
        for r_m in DISTANCES:
            ppm = {k: v * PEAK_FACTOR for k, v in ppm_at(q, r_m).items()}
            co2 = ppm_at({'CO2': qco2}, r_m)['CO2'] * PEAK_FACTOR
            for state in STATE_LIST:
                rates = ol.orn_rates(ppm, co2, gains=STATES[state]['orn'])
                v = to_input_vector(b.n, n, orn_rates=rates) + central_drive(n, state)
                pn = {g: 0.0 for g in pn_glom}
                active = 0.0
                for _ in range(3):
                    b.reset()
                    c = b.run(v, 300)
                    for g, idx in pn_glom.items():
                        pn[g] += c[idx].sum() / 3
                    active += (c > 0).sum() / 3
                num = sum(VALENCE.get(g, 0) * a for g, a in pn.items())
                den = sum(abs(VALENCE.get(g, 0)) * a for g, a in pn.items()) + 1e-9
                row = dict(food=label, distance_m=r_m, state=state, orn_total_hz=sum(rates.values()),
                           brain_active=active, pn_total=sum(pn.values()), valence_estimate=num / den if den > 1 else 0.0,
                           top_orn={k[4:]: round(val) for k, val in sorted(rates.items(), key=lambda x: -x[1])[:5]},
                           ppm={k: float(f'{val:.3g}') for k, val in ppm.items()}, co2_excess_ppm=co2)
                out.append(row)
                print(f"{label:17s} {r_m:5.2f} m {state:28s} ORN {row['orn_total_hz']:6.0f} Hz  PN {row['pn_total']:6.0f} "
                      f"active {active:6.0f}  valence {row['valence_estimate']:+.2f}", flush=True)
    json.dump(out, open('results/grid_smell.json', 'w'), indent=1)


if __name__ == '__main__':
    {'taste': taste_grid, 'smell': smell_grid}[sys.argv[1]]()
