"""Whole-CNS LIF model for foraging: stabilised MaleCNS with a working olfactory system.

Changes vs the base LIF (Shiu et al. 2024), all needed on this dataset:
  * spike adaptation; dopamine/serotonin/octopamine not fast synapses;
    no Kenyon cell <-> Kenyon cell contacts; sensory neurons are pure inputs
  * synapse gain 0.2 mV/contact (best trade-off found: taste reaches feeding/aversion neurons, odours stay distinct)
  * inhibitory synapses x3 (stand-in for slow GABA-B inhibition, absent from the model)
  * cholinergic antennal-lobe local neurons' chemical output x0.1 (in flies these excitatory LNs
    act mainly through gap junctions; as chemical synapses they formed a self-sustaining loop
    that made every odour look identical)
Inputs are per-neuron Poisson rates (Hz), so each receptor class can be driven by its own response.
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp

V_REST, V_TH = -52.0, -45.0
TAU_M, TAU_SYN, TAU_ADAPT = 20.0, 5.0, 200.0
ADAPT_INC = 3.0
T_REF, DELAY = 2, 2
W_SCALE = 0.2
INH_GAIN = 3.0
ELN_GAIN = 0.1
DT = 1.0
W_INPUT = 68.75

AL_LN = r'^(lLN|il3LN|v2LN|vLN|l2LN|lvLN|mLN)'


def load(path='data/'):
    W = sp.load_npz(path + 'brain.npz').tocsc()
    neurons = pd.read_parquet(path + 'neurons.parquet')
    ann = pd.read_feather(path + 'body-annotations-male-cns-v1.0-minconf-0.5.feather',
                          columns=['bodyId', 'flywireType', 'receptorType', 'entryNerve', 'synonyms'])
    neurons = neurons.merge(ann, on='bodyId', how='left')
    return W, neurons


def stabilise(W, neurons):
    t = neurons.type.fillna('')
    kc = t.str.startswith('KC').values
    modulatory = neurons.nt.isin(['dopamine', 'serotonin', 'octopamine']).values
    sensory = neurons.superclass.fillna('').str.contains('sensory').values
    eln = t.str.match(AL_LN).values & (neurons.nt.values == 'acetylcholine')
    coo = W.tocoo()
    keep = ~(kc[coo.row] & kc[coo.col]) & ~modulatory[coo.col] & ~sensory[coo.row]
    data = coo.data[keep].astype(np.float32)
    col = coo.col[keep]
    data = np.where(eln[col], data * ELN_GAIN, data)
    data = np.where(data < 0, data * INH_GAIN, data)
    return sp.csc_matrix((data, (coo.row[keep], col)), shape=W.shape)


def scramble(W, seed=0):
    rng = np.random.default_rng(seed)
    W = W.tocsc(copy=True)
    W.indices = rng.integers(0, W.shape[0], size=W.nnz).astype(W.indices.dtype)
    W.sum_duplicates()
    return W


class Brain:
    def __init__(self, W, neurons, seed=0):
        self.W = W.tocsc()
        self.n = W.shape[0]
        self.neurons = neurons
        self.rng = np.random.default_rng(seed)
        self.reset()

    def reset(self):
        self.v = np.full(self.n, V_REST, np.float32)
        self.g = np.zeros(self.n, np.float32)
        self.a = np.zeros(self.n, np.float32)
        self.ref = np.zeros(self.n, np.int16)
        self.queue = [np.zeros(0, np.int64) for _ in range(DELAY)]

    def run(self, rates, ms, record=None):
        """rates: dense float array of Poisson input rates in Hz (len n), or None."""
        counts = np.zeros(self.n, np.int32)
        dm = np.float32(DT / TAU_M)
        ds = np.float32(np.exp(-DT / TAU_SYN))
        da = np.float32(np.exp(-DT / TAU_ADAPT))
        if rates is not None:
            idx = np.flatnonzero(rates > 0)
            p = (rates[idx] * DT / 1000.).astype(np.float32)
        else:
            idx = np.zeros(0, np.int64)
        for _ in range(ms):
            delayed = self.queue.pop(0)
            if delayed.size:
                self.g += W_SCALE * np.asarray(self.W[:, delayed].sum(axis=1)).ravel()
            if idx.size:
                self.g[idx[self.rng.random(idx.size) < p]] += W_INPUT
            self.v += ((V_REST - self.v) + self.g - self.a) * dm
            self.g *= ds
            self.a *= da
            refractory = self.ref > 0
            self.v[refractory] = V_REST
            self.ref[refractory] -= 1
            spk = np.flatnonzero(self.v >= V_TH)
            self.v[spk] = V_REST
            self.ref[spk] = T_REF
            self.a[spk] += ADAPT_INC
            counts[spk] += 1
            self.queue.append(spk)
            if record is not None:
                record.append(spk)
        return counts
