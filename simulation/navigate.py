"""Flying flies searching a turbulent odour plume for banana, then tasting on landing.

Stimulus: per-fly instantaneous odour at its position (atmosphere.Plume) -> detection probability from the
connectome olfactory lookup (olf_lookup.py: PN population response above baseline in a 100 ms window).
Reaction:
  * surge upwind after a latency of 190 +/- 75 ms from detection, cast crosswind 450 +/- 165 ms after losing
    it (van Breugel & Dickinson 2014); casting reversals widen over time. These are behavioural parameters,
    not brain outputs: the connectome model did not produce lateralised steering (probe_bilateral.py).
  * optional valence gating: surge only if literature labelled-line glomerular valence of the evoked PN
    pattern is positive (CO2 counted attractive in flight, van Breugel 2018).
  * within 15 cm a surging fly steers to the visible fruit (vision guides odour localisation in free flight, Frye et al. 2003)
  * on landing (within 3 cm), feeding probability from the brain's taste acceptance (run.py taste grid)
    through a logistic calibrated on proboscis extension (feeding_calibration.py); rejecters take off again.
"""

import json

import numpy as np

from atmosphere import Plume

DT = 0.05
T_MAX = 90.0
LAND_R = 0.03
VALENCE_FLIGHT = {'DM1': 1, 'VA2': 1, 'DM2': 1, 'DM4': 1, 'D': -1, 'DA4l': -1, 'DA4m': -1, 'DC3': -1, 'DL1': -1,
                  'DL4': -1, 'DL5': -1, 'V': 1, 'DA2': -1, 'DC4': -1}
TASTE_LABEL = {'ripe': 'ripe (day 0)', 'rotting': 'rotting day 7'}


class Sensory:
    """Detection probability and valence as a function of plume concentration, from the brain lookup."""

    def __init__(self, path='results/olf_lookup.npz'):
        d = np.load(path, allow_pickle=True)
        self.levels = np.log10(d['levels'])
        self.gloms = list(d['gloms'])
        self.d = d
        self.cache = {}

    def table(self, food, temp, state):
        key = (food, temp, state)
        if key in self.cache:
            return self.cache[key]
        olf_state = 'starved_24h' if state.startswith('starved') else 'fed'
        name = f'resp_{food}_{temp}_{olf_state}'
        if name not in self.d:
            name = f'resp_{food}_{temp}_fed'
        resp = self.d[name]
        base = self.d[f'baseline_{olf_state}']
        thr = base[:, 0].mean() + 3 * base[:, 0].std()
        p_detect = (resp[:, :, 0] > thr).mean(1)
        p_detect = np.maximum.accumulate(p_detect)  # monotone in concentration
        glom_act = np.clip(resp[:, :, 2:].mean(1) - base[:, 2:].mean(0), 0, None)
        w = np.array([VALENCE_FLIGHT.get(g, 0) for g in self.gloms], float)
        valence = (glom_act @ w) / (np.abs(w) @ glom_act.T + 1e-9)
        valence = np.where(glom_act.sum(1) > 1, valence, 0.0)
        self.cache[key] = (p_detect, valence, thr)
        return self.cache[key]

    def query(self, food, temp, state, c_unit):
        p, v, _ = self.table(food, temp, state)
        x = np.log10(np.maximum(c_unit, 10 ** self.levels[0] / 10))
        return np.interp(x, self.levels, p, left=0.0, right=p[-1]), np.interp(x, self.levels, v, left=0.0, right=v[-1])


def feeding_probability(food, state, yeast_surface, calib, taste_rows, decision='logistic'):
    """logistic: calibrated on fed PER (55% at 800 mM), steep because model acceptance saturates.
    hyperbolic: p = A / (A + A_800), same calibration point, saturating shape (sensitivity analysis)."""
    enrich = yeast_surface if food == 'rotting' else 1.0
    row = next(r for r in taste_rows if r['food'] == TASTE_LABEL[food] and r['state'] == state
               and r['yeast_enrichment'] == enrich)
    a = row['accept']
    if decision == 'hyperbolic':
        a800 = float(np.interp(np.log10(800), np.log10(calib['conc_mM']), calib['accept_fed']))
        return a / (a + a800)
    return 1.0 / (1.0 + np.exp(-(a - calib['A50']) / calib['slope']))


def simulate(foods, positions, state='fed', temp=25, wind=0.3, n_flies=200, gated=False, yeast_surface=20.0, decision='logistic',
             seed=0, release_x=2.0, release_y=0.3, keep_tracks=0, record_frames=False, t_max=T_MAX, plume_growth=0.01):
    rng = np.random.default_rng(seed)
    sens = Sensory()
    calib = json.load(open('results/feeding_calibration.json'))
    taste_rows = json.load(open('results/grid_taste.json'))
    p_feed = np.array([feeding_probability(f, state, yeast_surface, calib, taste_rows, decision) for f in foods])
    plume = Plume(positions, wind=wind, x_max=release_x + 1.2, seed=seed + 1000, growth=plume_growth)
    for _ in range(int(max(20, 1.3 * (release_x + 1.2) / wind) / DT)):  # spin up until the plume fills the arena
        plume.step(DT)
    src = np.asarray(positions, float)
    N = n_flies
    pos = np.c_[rng.uniform(release_x - 0.05, release_x + 0.05, N), rng.uniform(-release_y, release_y, N)]
    mode = np.zeros(N, int)            # 0 cast, 1 surge, 2 done (fed), 3 lost
    cast_dir = rng.choice([-1, 1], N)
    cast_timer = np.zeros(N)
    cast_period = np.full(N, 0.6)
    surge_at = np.full(N, np.inf)
    cast_at = np.full(N, np.inf)
    last_detect = np.full(N, -np.inf)
    ignore_until = np.zeros((N, len(foods)))
    first_land = np.full(N, -1)
    fed_on = np.full(N, -1)
    t_land = np.full(N, np.nan)
    t_feed = np.full(N, np.nan)
    landings = np.zeros((N, len(foods)), int)
    detections = np.zeros(N, int)
    tracks = [[] for _ in range(keep_tracks)]
    frames = []
    t = 0.0
    while t < t_max and np.any(mode < 2):
        plume.step(DT)
        t += DT
        active = mode < 2
        c = plume.c_unit(pos[active])  # (n_active, n_sources)
        pdet = np.zeros_like(c)
        val = np.zeros_like(c)
        for k, f in enumerate(foods):
            pdet[:, k], val[:, k] = sens.query(f, temp, state, c[:, k])
        best = pdet.argmax(1)
        p100 = pdet[np.arange(len(best)), best]
        p_step = 1 - (1 - p100) ** (DT / 0.1)
        detected = rng.random(len(p_step)) < p_step
        if gated:
            detected &= val[np.arange(len(best)), best] > 0
        idx = np.flatnonzero(active)
        det_idx = idx[detected]
        detections[det_idx] += 1
        last_detect[det_idx] = t
        newly = det_idx[(mode[det_idx] == 0) & ~np.isfinite(surge_at[det_idx])]
        surge_at[newly] = t + np.clip(rng.normal(0.19, 0.075, len(newly)), 0.05, None)
        cast_at[det_idx] = np.inf
        go = idx[(surge_at[idx] <= t) & (mode[idx] == 0)]
        mode[go] = 1
        surge_at[go] = np.inf
        lost = idx[(mode[idx] == 1) & (last_detect[idx] < t) & ~np.isfinite(cast_at[idx])]
        cast_at[lost] = last_detect[lost] + np.clip(rng.normal(0.45, 0.165, len(lost)), 0.1, None)
        back = idx[(mode[idx] == 1) & (cast_at[idx] <= t)]
        mode[back] = 0
        cast_at[back] = np.inf
        cast_period[back] = 0.6
        cast_timer[back] = 0.0
        # movement (ground velocity)
        vel = np.zeros((N, 2))
        s = idx[mode[idx] == 1]
        ang = rng.normal(0, np.deg2rad(15), len(s))
        vel[s, 0] = -0.4 * np.cos(ang)
        vel[s, 1] = 0.4 * np.sin(ang)
        # close range: vision guides odour localisation (Frye et al. 2003) - a surging fly that
        # sees the fruit within 15 cm flies straight at it
        if len(s):
            ds = np.linalg.norm(pos[s][:, None, :] - src[None, :, :], axis=2)
            near = ds.min(1) < 0.15
            if near.any():
                tgt = src[ds.argmin(1)[near]] - pos[s[near]]
                vel[s[near]] = 0.3 * tgt / (np.linalg.norm(tgt, axis=1, keepdims=True) + 1e-9)
        cmask = idx[mode[idx] == 0]
        cast_timer[cmask] += DT
        flip = cmask[cast_timer[cmask] >= cast_period[cmask]]
        cast_dir[flip] *= -1
        cast_timer[flip] = 0
        cast_period[flip] = np.minimum(cast_period[flip] * 1.5, 3.0)
        vel[cmask, 0] = -0.05
        vel[cmask, 1] = 0.3 * cast_dir[cmask]
        pos[idx] += vel[idx] * DT
        # landing
        d = np.linalg.norm(pos[idx][:, None, :] - src[None, :, :], axis=2)
        for j in np.flatnonzero(d.min(1) < LAND_R):
            fi = idx[j]
            k = int(d[j].argmin())
            if t >= ignore_until[fi, k]:
                landings[fi, k] += 1
                if first_land[fi] < 0:
                    first_land[fi] = k
                    t_land[fi] = t
                if rng.random() < p_feed[k]:
                    mode[fi] = 2
                    fed_on[fi] = k
                    t_feed[fi] = t
                else:
                    ignore_until[fi, k] = t + 5.0
                    pos[fi] = src[k] + np.array([0.06, rng.uniform(-0.03, 0.03)])
                    mode[fi] = 0
                    cast_period[fi] = 0.6
        out_of = idx[(pos[idx, 0] > release_x + 1.0) | (pos[idx, 0] < -0.3) | (np.abs(pos[idx, 1]) > max(1.0, 0.35 * release_x))]
        mode[out_of[mode[out_of] < 2]] = 3
        for i in range(keep_tracks):
            tracks[i].append((t, pos[i, 0], pos[i, 1], mode[i]))
        if record_frames and int(round(t / DT)) % 2 == 0:
            frames.append(dict(t=t, puffs=plume.pos.copy(), sig=plume.sigma().copy(), sid=plume.sid.copy(),
                               flies=pos.copy(), mode=mode.copy(), fed_on=fed_on.copy()))
    res = dict(foods=foods, state=state, temp=temp, wind=wind, gated=gated, yeast_surface=yeast_surface, decision=decision,
               n=N, p_feed=p_feed.tolist(),
               first_landing={f: int((first_land == k).sum()) for k, f in enumerate(foods)},
               fed={f: int((fed_on == k).sum()) for k, f in enumerate(foods)},
               never_landed=int((first_land < 0).sum()),
               median_time_to_land=float(np.nanmedian(t_land)) if np.isfinite(t_land).any() else None,
               median_time_to_feed=float(np.nanmedian(t_feed)) if np.isfinite(t_feed).any() else None,
               detections_per_fly=float(detections.mean()))
    return res, tracks, frames
