"""Molecules -> receptor neuron firing rates -> connectome neuron groups.

Olfaction: DoOR 2.0 consensus responses (Muench & Galizia 2016, spontaneous rate subtracted) scaled to
spikes/s, with a Hill concentration dependence around the concentration DoOR/Hallem-style
recordings used (1e-2 dilution; research/receptors.json dose_response). CO2 via Gr21a/Gr63a (glomerulus V).
Taste: Hill fits per gustatory receptor neuron class from tip-recording data (research/receptors.json).
"""

import json

import numpy as np
import pandas as pd

R_GAS = 24.45e-3   # m3/mol at 25 C
ORN_MAX = 290.0    # spikes/s (Olsen 2010)
HILL_N, HILL_K = 0.8, 0.3  # population spike-rate Hill slope ~0.7-0.8 (Si 2019); ref conc ~ 77% of max

RESEARCH = 'research/'
DOOR = 'ref/door/'  # DoOR.data CSVs, fetched by scripts/fetch_door.sh

# compound -> (CAS, solvent used in receptor recordings)
COMPOUNDS = {
    'isoamyl_acetate': ('123-92-2', 'oil'), 'isobutyl_acetate': ('110-19-0', 'oil'),
    'butyl_acetate': ('123-86-4', 'oil'), 'ethyl_acetate': ('141-78-6', 'oil'), 'hexanal': ('66-25-1', 'oil'),
    '2-heptanone': ('110-43-0', 'oil'), 'isoamyl_alcohol': ('123-51-3', 'oil'), 'ethanol': ('64-17-5', 'water'),
    'acetoin': ('513-86-0', 'water'), '2-phenylethanol': ('60-12-8', 'oil'), 'acetic_acid': ('64-19-7', 'water'),
}


def _load_constants():
    chem = json.load(open(RESEARCH + 'chemistry.json'))['volatile_physical_constants']['compounds']
    return chem


class Olfaction:
    def __init__(self, neurons):
        t = neurons.type.fillna('').values
        self.orn_types = sorted(set(x for x in t if x.startswith('ORN_')))
        self.orn_index = {o: np.flatnonzero(t == o) for o in self.orn_types}
        gm = json.load(open(RESEARCH + 'receptors.json'))['olfactory']['glomerulus_map']
        self.rec2orn = {}
        for rec, info in gm.items():
            g = info.get('connectome_ORN') or info.get('glomerulus')
            if not g or '+' in str(g) and rec.startswith('ac'):
                continue
            gl = []
            for part in str(g).replace(' ', '').split('+'):
                if part == 'DL2d/v':
                    gl += ['DL2d', 'DL2v']
                elif part.startswith('VM6'):
                    gl += ['VM6v', 'VM6m', 'VM6l']
                else:
                    gl.append(part.split('(')[0])
            gl = [f'ORN_{x}' for x in gl if f'ORN_{x}' in self.orn_index]
            if gl:
                self.rec2orn[rec] = gl
        # DoOR consensus matrix, spontaneous (SFR) row subtracted, keyed by InChIKey.
        m = pd.read_csv(DOOR + 'door_response_matrix.csv', sep=';', index_col=0)
        sfr = m.loc['SFR']
        self.door = (m.drop(index='SFR') - sfr).clip(lower=0).fillna(0)
        odor = pd.read_csv(DOOR + 'odor.csv', sep=';')
        self.cas2key = dict(zip(odor.CAS, odor.InChIKey))
        self.const = _load_constants()

    def receptor_profile(self, cas):
        key = self.cas2key.get(cas)
        if key is None or key not in self.door.index:
            return None
        return self.door.loc[key]

    def orn_rates_from_profile(self, profile, scale=1.0):
        """profile: Series receptor -> normalised response at reference concentration."""
        rates = {}
        for rec, v in profile.items():
            if rec in self.rec2orn and v > 0:
                for o in self.rec2orn[rec]:
                    rates[o] = rates.get(o, 0.0) + v * ORN_MAX * scale
        return rates

    def reference_ppm(self, name):
        """Air concentration (ppm) at the fly for a 1e-2 dilution, the level DoOR/Hallem responses describe."""
        c = self.const[name]
        cas, solvent = COMPOUNDS[name]
        if solvent == 'water':
            density = {'ethanol': 789, 'acetic_acid': 1049, 'acetoin': 1010}.get(name, 1000)
            aq = 1e-2 * density / c['MW'] * 1000  # mol/m3 in water
            gas = aq * c['Kaw']
            ppm = gas * R_GAS * 1e6
        else:
            ppm = 1e-2 * c['VP_25C_mmHg'] / 760 * 1e6 * 0.04  # oil mole fraction ~0.04 of v/v (MW ~350)
        return ppm * 0.2  # carrier-stream mixing (Hallem 2006)

    def orn_rates(self, ppm_by_compound, co2_excess_ppm=0.0, gains=None):
        """Mixture -> per-ORN-type firing rate (spikes/s)."""
        drive = {}
        for name, ppm in ppm_by_compound.items():
            if ppm <= 0 or name not in COMPOUNDS:
                continue
            prof = self.receptor_profile(COMPOUNDS[name][0])
            if prof is None:
                continue
            c = ppm / self.reference_ppm(name)
            f = c ** HILL_N / (c ** HILL_N + HILL_K ** HILL_N) * (1 + HILL_K ** HILL_N)
            for o, r in self.orn_rates_from_profile(prof, f).items():
                drive[o] = drive.get(o, 0.0) + r
        if co2_excess_ppm > 0:  # Gr21a/Gr63a, behavioural threshold ~ +0.1% (1000 ppm)
            x = co2_excess_ppm / 1000.0
            drive['ORN_V'] = drive.get('ORN_V', 0.0) + 150 * x ** 1.5 / (x ** 1.5 + 1)
        rates = {o: ORN_MAX * np.tanh(r / ORN_MAX) for o, r in drive.items()}
        if gains:
            rates = {o: r * gains.get(o, 1.0) for o, r in rates.items()}
        return rates


def hill(x, ec50, n=1.0):
    x = max(x, 0.0)
    return x ** n / (x ** n + ec50 ** n)


BITTER_HALF_SUPPRESSION = 15.0  # spikes/s ESTIMATE

GRN_TYPES = {
    'sugar': ['LB3b', 'LB3c', 'LB3d'],
    'water': ['LB3a'],
    'bitter': ['LB1a', 'LB1b', 'LB1c', 'LB1d'],
    'ir94e': ['LB1e'],
    'yeast_peg': ['dorsal_tpGRN'],      # amino acid / yeast (Tastekin 2025; Steck 2018)
    'carbonation_peg': ['claw_tpGRN'],  # carbonation Ir56d (Fischler 2007; Tastekin 2025)
}


def taste_rates(t, gains=None):
    """Pulp composition (chemistry.composition()['taste']) -> spikes/s per GRN class."""
    g = dict(gains or {})
    # Hunger changes sugar and bitter GRN *sensitivity* (the concentration needed for a given response;
    # Inagaki 2012/2014), not their maximum firing, so those two gains scale effective concentration.
    s_sugar, s_bitter = g.pop('sugar', 1.0), g.pop('bitter', 1.0)
    # Sugar GRNs: sucrose EC50 ~60 mM, max ~94 spikes/s (Cameron 2010 / Dweck 2022 two-point fit).
    # Labellar L sensilla respond weakly to fructose (Jiao 2008); glucose ~half as effective (ESTIMATE).
    sugar_eq = s_sugar * (t['sucrose_mM'] + 0.5 * t['glucose_mM'] + 0.1 * t['fructose_mM'])
    acid_inhib = 1.0 / (1.0 + (10 ** -t['pH'] / 10 ** -3.5))  # acids suppress sweet GRNs (Charlu 2013) ESTIMATE
    sugar = 94 * hill(sugar_eq, 60) * acid_inhib
    # Water GRNs (ppk28): 12 spikes/s to water, silenced by osmolality (half ~150 mOsm, ESTIMATE).
    water = 12 / (1 + (t['osmolality_mOsm'] / 150) ** 4)
    # Bitter GRNs: tannin (weak agonist, threshold ~0.1-1 mM ESTIMATE), acidity below pH 5 (Charlu 2013),
    # ethanol (aversive taste, Devineni 2009; EC50 ~1 M ESTIMATE), K+ as non-selective high salt (EC50 300 mM).
    acid_drive = max(0.0, 5.0 - t['pH'])  # ~0 at pH 5, grows with acidity
    # Tannin and K+ act through thresholds (tannin 0.1-1 mM; high-salt cells from ~100 mM), so n = 2.
    bitter = 40 * min(1.0, hill(s_bitter * t['tannin_mM'], 0.5, 2.0) + hill(s_bitter * t['ethanol_mM'], 1000)
                      + 0.35 * hill(s_bitter * acid_drive, 1.0) + 0.5 * hill(s_bitter * t['K_mM'], 300, 2.0))
    # Ir94e: weak responses to amino acids / low salt (Shiu 2024; mildly aversive).
    ir94e = 15 * hill(t['amino_acids_mM'], 50)
    # Taste peg yeast neurons: 10% yeast > 500 mM sucrose (Steck 2018); EC50 ~2% w/v yeast (ESTIMATE).
    yeast_peg = 80 * hill(t['yeast_percent_w_v'], 2.0) + 10 * hill(t['amino_acids_mM'], 25)
    # Taste peg carbonation neurons: dissolved CO2, saturation ~34 mM; EC50 10 mM (ESTIMATE, no data).
    carb = 60 * hill(t['CO2_mM'], 10.0, 2.0)
    rates = dict(sugar=sugar, water=water, bitter=bitter, ir94e=ir94e, yeast_peg=yeast_peg, carbonation_peg=carb)
    rates = {k: v * g.get(k, 1.0) for k, v in rates.items()}
    # Bitter-sweet integration at the sugar GRN output: GABAergic presynaptic inhibition via GABA-B receptors on
    # sweet GRN terminals (Chu et al. 2014) plus direct inhibition of sweet GRNs by bitter compounds
    # (Jeong et al. 2013). Not representable in a connectome LIF model, so applied here. Half-suppression at
    # 15 spikes/s bitter GRN firing is an ESTIMATE. Applied after state gains so starved flies' reduced bitter
    # output also releases the brake.
    rates['sugar'] = rates['sugar'] / (1.0 + rates['bitter'] / BITTER_HALF_SUPPRESSION)
    return rates


def to_input_vector(n_neurons, neurons, orn_rates=None, grn_rates=None, spont_orn=0.0):
    t = neurons.type.fillna('').values
    v = np.zeros(n_neurons, np.float32)
    if spont_orn:
        v[pd.Series(t).str.startswith('ORN_').values] = spont_orn
    for o, r in (orn_rates or {}).items():
        v[t == o] += r
    for cls, r in (grn_rates or {}).items():
        v[np.isin(t, GRN_TYPES[cls])] += r
    return v
