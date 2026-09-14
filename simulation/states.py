"""Internal states: peripheral sensory gains from the literature + tonic drive to the connectome's own
hunger/satiety neurons. Fold changes from research/receptors.json internal_state_modulation.

This connectome is male, so the strongest known protein-appetite state (mated female, sex peptide)
cannot be represented; 'protein_deprived' uses the taste-peg gain reported for yeast-deprived flies.
"""

import numpy as np

# Peripheral gains (multiplicative on receptor neuron output).
STATES = {
    'fed': dict(grn={}, orn={}, central={'satiety': 15.0}),
    'starved_24h': dict(
        grn={'sugar': 2.0, 'bitter': 0.6},                       # Inagaki 2012/2014, LeDue 2016
        orn={'ORN_DM1': 1.75, 'ORN_DM4': 1.3, 'ORN_DM2': 1.3,   # Root 2011 (sNPF facilitation)
             'ORN_DM5': 0.25,                                     # Ko 2015 (tachykinin), threshold x4
             'ORN_VM2': 0.8, 'ORN_VA3': 0.8},
        central={'hunger': 15.0}),
    'protein_deprived': dict(
        grn={'yeast_peg': 1.5},                                  # Steck 2018 (yeast taste peg gain) ESTIMATE
        orn={}, central={'satiety': 15.0}),
    'starved_and_protein_deprived': dict(
        grn={'sugar': 2.0, 'bitter': 0.6, 'yeast_peg': 1.5},
        orn={'ORN_DM1': 1.75, 'ORN_DM4': 1.3, 'ORN_DM2': 1.3, 'ORN_DM5': 0.25, 'ORN_VM2': 0.8, 'ORN_VA3': 0.8},
        central={'hunger': 15.0}),
}

# Connectome neurons carrying internal state. Tonic Poisson drive (Hz) during the state.
HUNGER_TYPES = ['DNp29', 'NPFL1-I']       # NPF neurons (hunger; Wu 2003)
SATIETY_TYPES = ['IPC', 'Hugin-RG']        # insulin-producing cells, hugin (satiety / feeding suppression)


def central_drive(neurons, state):
    t = neurons.type.fillna('').values
    v = np.zeros(len(t), np.float32)
    for kind, hz in STATES[state]['central'].items():
        types = HUNGER_TYPES if kind == 'hunger' else SATIETY_TYPES
        v[np.isin(t, types)] = hz
    return v
