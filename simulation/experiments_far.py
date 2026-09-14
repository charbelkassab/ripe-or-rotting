"""Long-range and control experiments (Results 3.6-3.7) -> merged into results/experiments.json

  python simulation/experiments_far.py

  * release_distance 10 / 15 m in the baseline plume, starved and fed
  * far_field_plume: the same two-choice test with a plume that dilutes as expected (puff growth 0.05)
  * same_food_control: identical fruit on both sides, to check the arena itself has no side bias
"""

import json

import numpy as np

import experiments as E
import navigate as nv

SEEDS = list(range(4))


def dump(res):
    json.dump(res, open('results/experiments.json', 'w'), indent=1,
              default=lambda o: o.item() if hasattr(o, 'item') else str(o))


def share(r, food='rotting'):
    tot = sum(r['first_landing'].values())
    return r['first_landing'][food] / tot if tot else float('nan')


def main():
    E.SEEDS = SEEDS
    res = json.load(open('results/experiments.json'))
    res.pop('separation', None)

    for rx in [10.0, 15.0]:
        for state, key in [('starved_24h', str(rx)), ('fed', f'{rx}_fed')]:
            r = E.pooled(*E.TWO, state=state, decision='hyperbolic', release_x=rx, release_y=0.15 * rx, t_max=240)
            res['release_distance'][key] = r
            print(f'baseline plume, release {rx} m, {state}: rotting first {share(r):.2f} '
                  f'(sd {r["share_sd_across_plumes"]:.2f}), never landed {r["never_landed"]} of {r["n"]}', flush=True)
    dump(res)

    res['far_field_plume'] = {}
    for sep, rx in [(0.3, 2.0), (0.3, 10.0), (4.0, 10.0)]:
        for state in ['fed', 'starved_24h']:
            pos = [(0.0, sep / 2), (0.0, -sep / 2)]
            r = E.pooled(['ripe', 'rotting'], pos, state=state, decision='hyperbolic', release_x=rx,
                         release_y=max(0.3, min(sep, 0.15 * rx)), t_max=240, plume_growth=0.05)
            r['rotting_share'] = share(r)
            res['far_field_plume'][f'{sep}|{rx}|{state}'] = r
            print(f'diluting plume, {sep} m apart, release {rx} m, {state}: rotting first {r["rotting_share"]:.2f} '
                  f'(sd {r["share_sd_across_plumes"]:.2f}), never landed {r["never_landed"]} of {r["n"]}', flush=True)
    dump(res)

    # Same fruit on both sides: '#2' labels share the first fruit's odour and taste tables.
    orig_query = nv.Sensory.query
    nv.Sensory.query = lambda self, f, T, s, c: orig_query(self, f.split('#')[0], T, s, c)
    nv.TASTE_LABEL.update({'ripe#2': 'ripe (day 0)', 'rotting#2': 'rotting day 7'})
    res['same_food_control'] = {}
    for food in ['ripe', 'rotting']:
        r = E.pooled([food, f'{food}#2'], E.TWO[1], state='starved_24h', decision='hyperbolic')
        r['second_share'] = share(r, f'{food}#2')
        res['same_food_control'][food] = r
        print(f'same-food control, {food}: second copy first {r["second_share"]:.2f} '
              f'(sd {r["share_sd_across_plumes"]:.2f})', flush=True)
    nv.Sensory.query = orig_query
    dump(res)


if __name__ == '__main__':
    main()
