"""Banana chemistry over time, from kinetics instead of hand-picked 'sweet/yeasty' numbers.

Two scenarios, both starting at eating-ripe (day 0, Phillips et al. 2021 composition):
  * 'intact'  : fruit keeps ripening (starch -> sugars, fruit esters rise), no yeast colonisation
  * 'damaged' : split peel colonised by yeast (1e5 CFU/g) and acetic acid bacteria -> fermentation

Units: pulp concentrations in g/kg fresh weight (FW) unless noted; tissue water 0.75-0.79 kg/kg.
Every parameter is from research/chemistry.json; ones marked ESTIMATE there are estimates here too.
"""

import numpy as np

MW = {'glucose': 180.16, 'fructose': 180.16, 'sucrose': 342.3, 'ethanol': 46.07, 'acetic_acid': 60.05,
      'ethyl_acetate': 88.11, 'acetoin': 88.11, '2-phenylethanol': 122.17, 'isoamyl_alcohol': 88.15,
      'isoamyl_acetate': 130.19, 'isobutyl_acetate': 116.16, 'butyl_acetate': 116.16, 'hexanal': 100.16,
      '2-heptanone': 114.19, 'glycerol': 92.09, 'malic_acid': 134.09, 'citric_acid': 192.12, 'K': 39.1,
      'CO2': 44.01, 'amino_acids': 130.0, 'tannin': 290.0}

# Ripe pulp at day 0 (g/kg FW), Phillips 2021 'ripe'; USDA minerals/protein; Bashmil 2021 tannins.
RIPE = dict(starch=25.2, sucrose=45.6, glucose=53.3, fructose=62.3, water=753.0, K=3.58, protein=10.9,
            malic_acid=3.0, citric_acid=1.0, tannin=0.02, free_amino_acids=1.0)
# Overripe (day 5, intact): Phillips 2021 'overripe'.
OVERRIPE = dict(starch=4.8, sucrose=18.8, glucose=69.6, fructose=67.2, water=786.0)

# Fruit-made volatiles, mg/kg FW (chemistry.json per_compound_defaults ESTIMATE); interpolated ripe -> overripe.
FRUIT_VOLATILES = {  # (ripe, overripe)
    'isoamyl_acetate': (3.0, 5.0), 'isobutyl_acetate': (1.0, 1.5), 'butyl_acetate': (0.5, 0.8),
    'ethyl_acetate': (0.5, 2.0), 'hexanal': (0.5, 0.2), '2-heptanone': (0.1, 0.3),
    'isoamyl_alcohol': (1.5, 3.0), 'ethanol': (300.0, 600.0),
}

# Fermentation parameters.
P = dict(
    mu_max=0.22,        # 1/h, S. cerevisiae anaerobic at 25 C (Verduyn 1990, Q10-corrected)
    Ks=0.1,             # g/L hexose
    Yxs=0.10,           # g dry yeast / g sugar (anaerobic)
    X0=1.5e-3,          # g dry yeast / kg: 1e5 CFU/g x 15 pg (damaged fruit inoculum)
    Xmax=1.5,           # g/kg: 1e8 CFU/g x 15 pg
    q_ferm=0.2,         # g sugar / g yeast / h, non-growth fermentation (ESTIMATE, gives ~2% ethanol by day 7)
    Y_E=0.46,           # g ethanol / g sugar (practical)
    E_max=100.0,        # g/kg ethanol that stops fermentation (10%)
    Z0=1e-4, Zmax=0.2, mu_Z=0.1,   # acetic acid bacteria (relative biomass, 1/h) ESTIMATE
    k_A=0.05,           # g ethanol oxidised per unit Z per h at saturating ethanol ESTIMATE
    K_E=1.0,            # g/kg
    eta_A=0.53 * 1.304, # g acetic / g ethanol (banana vinegar 53% of theoretical)
    y_EA=0.004, y_acetoin=0.003, k_acetoin=0.01, y_PE=0.002, y_IAOH=0.008,  # g / g ethanol made (wine ranges)
    y_glycerol=0.07,    # g / g sugar
    k_starch=0.35 / 24, k_inv=0.2 / 24,  # 1/h fruit enzymes (ESTIMATE)
    protein_frac=0.45,  # yeast protein fraction of dry weight
)


Q10 = 2.0  # biological rate temperature coefficient (chemistry.json S. cerevisiae Q10 ~2 ESTIMATE)
RATE_KEYS = ['mu_max', 'q_ferm', 'mu_Z', 'k_A', 'k_starch', 'k_inv']


def simulate(scenario, days=14, dt_h=0.1, temp_C=25.0):
    """Integrate composition. Returns dict of arrays sampled every hour plus 't_days'.
    Biological and enzymatic rates scale with Q10 = 2 around 25 C."""
    f_T = Q10 ** ((temp_C - 25.0) / 10.0)
    P = dict(globals()['P'])
    for k in RATE_KEYS:
        P[k] = P[k] * f_T
    s = dict(starch=RIPE['starch'], sucrose=RIPE['sucrose'], hexose=RIPE['glucose'] + RIPE['fructose'],
             X=P['X0'] if scenario == 'damaged' else 0.0, Z=P['Z0'] if scenario == 'damaged' else 0.0,
             E=0.0, A=0.0, EA=0.0, acetoin=0.0, PE=0.0, IAOH=0.0, glycerol=0.0, fermented=0.0)
    out = {k: [] for k in list(s) + ['t_days', 'ferm_rate']}
    steps = int(days * 24 / dt_h)
    ferm_rate = 0.0
    for i in range(steps + 1):
        t = i * dt_h
        if i % int(1 / dt_h) == 0:
            for k, v in s.items():
                out[k].append(v)
            out['t_days'].append(t / 24)
            out['ferm_rate'].append(ferm_rate)
        # fruit enzymes
        d_starch = -P['k_starch'] * s['starch']
        d_suc = -0.5 * d_starch - P['k_inv'] * s['sucrose']
        d_hex = -0.5 * d_starch + 1.053 * P['k_inv'] * s['sucrose']
        # yeast
        water_L = 0.77
        hex_gL = s['hexose'] / water_L
        inhib = max(0.0, 1 - s['E'] / P['E_max'])
        mu = P['mu_max'] * hex_gL / (P['Ks'] + hex_gL) * inhib
        dX = mu * s['X'] * (1 - s['X'] / P['Xmax'])
        uptake = (dX / P['Yxs'] + P['q_ferm'] * s['X'] * inhib * hex_gL / (P['Ks'] + hex_gL))
        uptake = min(uptake, s['hexose'] / dt_h)
        dE_made = P['Y_E'] * uptake
        # acetic acid bacteria at the surface
        dZ = P['mu_Z'] * s['Z'] * (1 - s['Z'] / P['Zmax']) * (1 if s['E'] > 0.5 else 0.2)
        oxid = P['k_A'] * s['Z'] / P['Zmax'] * s['E'] / (P['K_E'] + s['E'])
        oxid = min(oxid, s['E'] / dt_h)
        ferm_rate = uptake
        s['starch'] += d_starch * dt_h
        s['sucrose'] += d_suc * dt_h
        s['hexose'] += (d_hex - uptake) * dt_h
        s['X'] += dX * dt_h
        s['Z'] += dZ * dt_h
        s['E'] += (dE_made - oxid) * dt_h
        s['A'] += P['eta_A'] * oxid * dt_h
        s['EA'] += P['y_EA'] * dE_made * dt_h
        s['acetoin'] += (P['y_acetoin'] * dE_made - P['k_acetoin'] * s['acetoin']) * dt_h
        s['PE'] += P['y_PE'] * dE_made * dt_h
        s['IAOH'] += P['y_IAOH'] * dE_made * dt_h
        s['glycerol'] += P['y_glycerol'] * uptake * dt_h
        s['fermented'] += uptake * dt_h
    return {k: np.array(v) for k, v in out.items()}


def ph_from_acids(acids_mM, base_mM):
    """Charge balance: base (K+ salts) + H+ = OH- + sum(dissociated weak acids). Returns pH."""
    pkas = {'malic_acid': [3.40, 5.11], 'citric_acid': [3.13, 4.76, 6.40], 'acetic_acid': [4.76],
            'lactic_acid': [3.86], 'CO2': [6.35]}
    def charge(ph):
        h = 10 ** -ph
        neg = 0.0
        for name, conc in acids_mM.items():
            ks = [10 ** -p for p in pkas[name]]
            # fraction-weighted average charge of a polyprotic acid
            terms = [1.0]
            prod = 1.0
            for k in ks:
                prod *= k / h
                terms.append(prod)
            z = sum(j * terms[j] for j in range(len(terms))) / sum(terms)
            neg += conc * z
        return base_mM + h * 1e3 - neg - (1e-14 / h) * 1e3
    lo, hi = 2.0, 8.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if charge(mid) > 0:  # net positive -> pH must be higher
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


_BASE = None


def composition(scenario, day, temp_C=25.0):
    """Composition of pulp at a given day: taste-relevant concentrations (mM in tissue water)
    and volatile contents (mg/kg FW)."""
    global _BASE
    sim = simulate(scenario, days=max(day, 0.01) + 0.1, temp_C=temp_C)
    i = int(round(day * 24))
    g = {k: v[min(i, len(v) - 1)] for k, v in sim.items()}
    ripeness = min(day * Q10 ** ((temp_C - 25.0) / 10.0) / 5.0, 1.0)  # 0 ripe -> 1 overripe
    water = RIPE['water'] + (OVERRIPE['water'] - RIPE['water']) * ripeness
    L = water / 1000.0  # litres of tissue water per kg
    hexose = g['hexose']
    glu_frac = (RIPE['glucose'] + (OVERRIPE['glucose'] - RIPE['glucose']) * ripeness)
    fru_frac = (RIPE['fructose'] + (OVERRIPE['fructose'] - RIPE['fructose']) * ripeness)
    glucose = hexose * glu_frac / (glu_frac + fru_frac)
    fructose = hexose - glucose
    mM = lambda g_per_kg, mw: g_per_kg / mw / L * 1000
    fruit_vol = {k: (a + (b - a) * ripeness) for k, (a, b) in FRUIT_VOLATILES.items()}
    vol_mg_kg = dict(fruit_vol)
    vol_mg_kg['ethanol'] += g['E'] * 1000
    vol_mg_kg['ethyl_acetate'] += g['EA'] * 1000
    vol_mg_kg['isoamyl_alcohol'] += g['IAOH'] * 1000
    vol_mg_kg['acetoin'] = g['acetoin'] * 1000
    vol_mg_kg['2-phenylethanol'] = g['PE'] * 1000
    vol_mg_kg['acetic_acid'] = g['A'] * 1000
    acids = {'malic_acid': mM(RIPE['malic_acid'], MW['malic_acid']),
             'citric_acid': mM(RIPE['citric_acid'], MW['citric_acid']),
             'acetic_acid': mM(g['A'], MW['acetic_acid'])}
    fermenting = g['ferm_rate'] > 0.05
    co2_mM = 34.0 if fermenting else 3.0  # saturated (1 atm) while fermenting ESTIMATE
    if _BASE is None:  # calibrate buffer base so ripe pulp is pH 4.9 (chemistry.json)
        lo, hi = 0.0, 200.0
        L0 = RIPE['water'] / 1000.0  # calibrate on ripe pulp water, independent of which call comes first
        ripe_acids = {'malic_acid': RIPE['malic_acid'] / MW['malic_acid'] / L0 * 1000,
                      'citric_acid': RIPE['citric_acid'] / MW['citric_acid'] / L0 * 1000, 'acetic_acid': 0.0}
        for _ in range(60):
            mid = (lo + hi) / 2
            if ph_from_acids(ripe_acids, mid) < 4.9:
                lo = mid
            else:
                hi = mid
        _BASE = (lo + hi) / 2
    ph = ph_from_acids({**acids, 'CO2': co2_mM}, _BASE)
    taste = {
        'sucrose_mM': mM(g['sucrose'], MW['sucrose']), 'glucose_mM': mM(glucose, MW['glucose']),
        'fructose_mM': mM(fructose, MW['fructose']), 'ethanol_mM': mM(g['E'] + fruit_vol['ethanol'] / 1000, MW['ethanol']),
        'acetic_acid_mM': acids['acetic_acid'], 'pH': ph, 'K_mM': mM(RIPE['K'], MW['K']), 'Na_mM': 0.6,
        'tannin_mM': mM(RIPE['tannin'], MW['tannin']), 'CO2_mM': co2_mM,
        'yeast_percent_w_v': g['X'] / L / 10,  # g dry yeast per 100 mL tissue water
        'amino_acids_mM': mM(RIPE['free_amino_acids'] + P['protein_frac'] * g['X'] * 0.05, MW['amino_acids']),
        'glycerol_mM': mM(g['glycerol'], MW['glycerol']),
    }
    osm = (taste['sucrose_mM'] + taste['glucose_mM'] + taste['fructose_mM'] + taste['ethanol_mM']
           + taste['acetic_acid_mM'] + 2 * taste['K_mM'] + taste['glycerol_mM'] + sum(acids.values()))
    taste['osmolality_mOsm'] = osm
    extra = {'temp_C': temp_C, 'yeast_g_per_kg': g['X'], 'yeast_protein_g_per_kg': P['protein_frac'] * g['X'],
             'sugar_fermented_g_per_kg': g['fermented'], 'fermentation_rate_g_per_kg_h': g['ferm_rate'],
             'starch_g_per_kg': g['starch'], 'ripeness': ripeness}
    return {'taste': taste, 'volatiles_mg_per_kg': vol_mg_kg, 'extra': extra, 'water_L_per_kg': L}


if __name__ == '__main__':
    for sc, day in [('intact', 0), ('intact', 5), ('damaged', 2), ('damaged', 4), ('damaged', 7), ('damaged', 14)]:
        c = composition(sc, day)
        t, v, e = c['taste'], c['volatiles_mg_per_kg'], c['extra']
        print(f"{sc:7s} day {day:2d} | sugars suc {t['sucrose_mM']:5.0f} glu {t['glucose_mM']:5.0f} fru {t['fructose_mM']:5.0f} mM"
              f" | EtOH {v['ethanol'] / 1e4:5.2f}% acetic {v['acetic_acid'] / 1e4:5.3f}% pH {t['pH']:.2f} osm {t['osmolality_mOsm']:5.0f}"
              f" | yeast {e['yeast_g_per_kg']:.2f} g/kg | EA {v['ethyl_acetate']:5.1f} IAA {v['isoamyl_acetate']:.1f} "
              f"PE {v['2-phenylethanol']:5.1f} acetoin {v['acetoin']:5.1f} mg/kg")
