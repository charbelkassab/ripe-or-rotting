"""Volatile emission from a banana and transport to a fly at distance r.

Emission (film theory): J_i = h_i * A * C_air,surface,i ; h_i = D_i / delta ; C_air,surface = Kaw_i * C_aq,i
(acetic acid: only the undissociated fraction partitions, f = 1/(1+10^(pH-4.76))).
Transport, two cases:
  * breeze: time-averaged Gaussian plume on the centreline, C = Q / (pi u sy sz), sy = sz = 0.1 r + 0.005 m
  * still air: steady diffusion from a point source, C = Q / (4 pi D r)
Instantaneous plume filaments are far more concentrated than the time average (intermittency);
the fly samples those, so we also report a filament peak = time-average x PEAK_FACTOR (ESTIMATE).
"""

import json

import numpy as np

from chemistry import MW

R_GAS = 24.45e-3          # m3/mol
DELTA = 2e-3              # m, boundary layer thickness (chemistry.json: 1-3 mm)
EXPOSED_AREA = 30e-4      # m2: a split / half-peeled banana exposing ~30 cm2 of pulp (ESTIMATE)
WIND = 0.3                # m/s light indoor breeze
PEAK_FACTOR = 10.0        # filament concentration / time average (ESTIMATE; plumes are highly intermittent)
CONST = json.load(open('research/chemistry.json'))['volatile_physical_constants']['compounds']


def kaw_at(name, temp_C):
    """Air/water partition coefficient at temperature (van 't Hoff on Henry solubility; Kaw = 1/(H R T))."""
    c = CONST[name]
    T = temp_C + 273.15
    dln = c.get('dlnH_d1overT_K') or 6000.0  # typical organic volatile when not tabulated (ESTIMATE)
    return c['Kaw'] * np.exp(-dln * (1 / T - 1 / 298.15)) * 298.15 / T


def emission_mol_s(volatiles_mg_per_kg, water_L_per_kg, pH, temp_C=25.0):
    q = {}
    for name, mg in volatiles_mg_per_kg.items():
        if name not in CONST or mg <= 0:
            continue
        c = CONST[name]
        c_aq = mg / 1000 / MW[name] / water_L_per_kg * 1000  # mol/m3 in tissue water
        kaw = kaw_at(name, temp_C)
        if name == 'acetic_acid':
            kaw = kaw / (1 + 10 ** (pH - 4.76))
        c_surf = kaw * c_aq
        h = c['D_air_cm2_s_best'] * 1e-4 / DELTA
        q[name] = h * EXPOSED_AREA * c_surf
    return q


def co2_emission_mol_s(fermentation_g_per_kg_h, respiration=True, fruit_kg=0.12):
    ml_min = (0.2 if respiration else 0.0) + fermentation_g_per_kg_h * fruit_kg * 272 / 60  # 272 mL CO2 / g sugar
    return ml_min / 1e6 / 60 / R_GAS


def ppm_at(q_mol_s, r_m, mode='breeze'):
    out = {}
    for name, q in q_mol_s.items():
        if mode == 'breeze':
            s = 0.1 * r_m + 0.005
            c = q / (np.pi * WIND * s * s)
        else:
            d = CONST[name]['D_air_cm2_s_best'] * 1e-4 if name in CONST else 1.6e-5
            c = q / (4 * np.pi * d * max(r_m, 0.005))
        out[name] = c * R_GAS * 1e6
    return out
