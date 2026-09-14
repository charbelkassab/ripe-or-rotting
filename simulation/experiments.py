"""All stimulus-reaction experiments -> results/experiments.json (+ results/plume_series.npz, results/tracks.json)."""
import json

import numpy as np

import atmosphere
from atmosphere import Plume
from navigate import Sensory, simulate

TWO = (['ripe', 'rotting'], [(0.0, 0.15), (0.0, -0.15)])
STATES = ['fed', 'starved_24h', 'protein_deprived', 'starved_and_protein_deprived']
SEEDS = list(range(8))
N = 120


def pooled(foods, positions, **kw):
    """Pool plume seeds. With two sources every plume realization is run in both orientations (fruit positions
    swapped), because a single realization's meander can favour one side for tens of seconds."""
    agg = None
    per_run_share = []
    orientations = [positions] if len(positions) < 2 else [positions, positions[::-1]]
    for s in SEEDS:
        for pos in orientations:
            r, _, _ = simulate(foods, pos, n_flies=N, seed=s, **kw)
            if len(foods) == 2:
                tot = sum(r['first_landing'].values())
                per_run_share.append(r['first_landing'][foods[1]] / tot if tot else np.nan)
            if agg is None:
                agg = r
                agg['median_time_to_land'] = [r['median_time_to_land']]
            else:
                agg['n'] += r['n']
                agg['never_landed'] += r['never_landed']
                for key in ['first_landing', 'fed']:
                    for f in foods:
                        agg[key][f] += r[key][f]
                agg['median_time_to_land'].append(r['median_time_to_land'])
    if per_run_share:
        # orientation pairs share a plume; average within pair, SD across plume realizations
        pairs = np.array(per_run_share).reshape(-1, 2).mean(1)
        agg['second_food_first_landing_share_by_plume'] = pairs.tolist()
        agg['share_sd_across_plumes'] = float(np.nanstd(pairs))
    return agg


def detection_vs_distance(duration=60.0, dt=0.05):
    dists = np.array([0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 2.8])
    p = Plume([(0.0, 0.0)], x_max=3.2, seed=7)
    for _ in range(int(20 / dt)):
        p.step(dt)
    series = []
    for _ in range(int(duration / dt)):
        p.step(dt)
        series.append(p.c_unit(np.c_[dists, np.zeros_like(dists)])[:, 0])
    S = np.array(series)
    sens = Sensory()
    out = {}
    for food in ['ripe', 'rotting']:
        for state in ['fed', 'starved_24h']:
            pd_, _ = sens.query(food, 25, state, S)
            p_step = 1 - (1 - pd_) ** (dt / 0.1)
            # probability of at least one detection within a 1 s sampling window
            win = int(1.0 / dt)
            chunks = p_step[: (len(p_step) // win) * win].reshape(-1, win, len(dists))
            p1s = 1 - np.prod(1 - chunks, axis=1)
            out[f'{food}_{state}'] = p1s.mean(0).tolist()
    return dists.tolist(), out, S


def main():
    res = {}
    cal, S = atmosphere.calibrate(duration=200)
    res['plume_calibration'] = [{k: float(v) for k, v in r.items()} for r in cal]
    np.savez('results/plume_series.npz', series=S, distances=np.array([0.1, 0.3, 1.0, 2.0]), dt=0.02)
    print('plume calibration done', flush=True)
    d, det, _ = detection_vs_distance()
    res['detection_vs_distance'] = {'distances_m': d, **det}
    print('detection done', flush=True)
    res['single_source'] = {}
    for food in ['rotting', 'ripe']:
        for state in ['fed', 'starved_24h']:
            r = pooled([food], [(0.0, 0.0)], state=state)
            res['single_source'][f'{food}_{state}'] = r
            print('single', food, state, r['first_landing'], r['fed'], flush=True)
    res['two_choice'] = {}
    for decision in ['logistic', 'hyperbolic']:
      for gated in [False, True]:
        for ys in [20.0, 1.0]:
            for state in STATES:
                r = pooled(*TWO, state=state, gated=gated, yeast_surface=ys, decision=decision)
                key = f"{state}|yeast{int(ys)}|{'gated' if gated else 'detect'}|{decision}"
                res['two_choice'][key] = r
                print('two', key, 'first', r['first_landing'], 'fed', r['fed'], 'never', r['never_landed'], flush=True)
    res['wind'] = {}
    for w in [0.15, 0.3, 0.6]:
        r = pooled(*TWO, state='starved_24h', wind=w, decision='hyperbolic')
        res['wind'][str(w)] = r
        print('wind', w, r['first_landing'], r['fed'], flush=True)
    res['release_distance'] = {}
    for rx in [1.0, 2.0, 3.5, 5.0]:
        r = pooled(*TWO, state='starved_24h', decision='hyperbolic', release_x=rx)
        res['release_distance'][str(rx)] = r
        print('release', rx, r['first_landing'], r['fed'], 'never', r['never_landed'], flush=True)
    res['temperature'] = {}
    for T in [18, 25, 32]:
        r = pooled(*TWO, state='fed', temp=T, decision='hyperbolic')
        res['temperature'][str(T)] = r
        print('temp', T, r['first_landing'], r['fed'], flush=True)
    json.dump(res, open('results/experiments.json', 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    # example tracks for figures
    _, tracks, _ = simulate(*TWO, state='starved_24h', n_flies=40, seed=3, keep_tracks=40, t_max=60)
    json.dump(tracks, open('results/tracks.json', 'w'), default=lambda o: o.item() if hasattr(o, 'item') else str(o))


if __name__ == '__main__':
    main()
