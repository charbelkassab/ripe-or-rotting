"""Figures for the paper -> paper/fig*.png"""

import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from atmosphere import Plume
from chemistry import composition
from physics import emission_mol_s

SURFACE, INK, INK2, GRID = '#fbfcfc', '#15201c', '#4f5b57', '#e3e8e6'
BLUE, ORANGE, AQUA, YELLOW, MAGENTA = '#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4'
STATE_COL = {'fed': BLUE, 'starved_24h': ORANGE, 'protein_deprived': AQUA, 'starved_and_protein_deprived': YELLOW}
STATE_LAB = {'fed': 'fed', 'starved_24h': 'starved 24 h', 'protein_deprived': 'protein-deprived',
             'starved_and_protein_deprived': 'starved + protein-deprived'}
plt.rcParams.update({'font.family': ['Helvetica Neue', 'DejaVu Sans'], 'font.size': 9.5, 'axes.edgecolor': GRID,
                     'axes.labelcolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2, 'axes.titlecolor': INK,
                     'axes.titlesize': 10.5, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left',
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.spines.top': False,
                     'axes.spines.right': False, 'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE,
                     'legend.frameon': False, 'lines.linewidth': 2, 'xtick.major.size': 0, 'ytick.major.size': 0,
                     'savefig.facecolor': SURFACE})
os.makedirs('paper', exist_ok=True)


def tag(ax, letter):
    ax.text(-0.14, 1.08, letter, transform=ax.transAxes, fontsize=13, fontweight='bold', color=INK, va='bottom')


def fig_chemistry():
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.4))
    days = np.arange(0, 14.01, 0.5)
    comps = {sc: [composition(sc, d) for d in days] for sc in ['intact', 'damaged']}
    panels = [('Sugars in tissue water (mM)', lambda c: c['taste']['sucrose_mM'] + c['taste']['glucose_mM'] + c['taste']['fructose_mM']),
              ('Ethanol (% w/w)', lambda c: c['volatiles_mg_per_kg']['ethanol'] / 1e4),
              ('Acetic acid (% w/w)', lambda c: c['volatiles_mg_per_kg']['acetic_acid'] / 1e4),
              ('pH (charge balance)', lambda c: c['taste']['pH']),
              ('Yeast biomass (g dry / kg)', lambda c: c['extra']['yeast_g_per_kg'])]
    for ax, (title, f), letter in zip(axes.flat[:5], panels, 'ABCDE'):
        for sc, col, lab in [('intact', BLUE, 'intact, ripening'), ('damaged', ORANGE, 'damaged, yeast-colonised')]:
            ax.plot(days, [f(c) for c in comps[sc]], color=col, label=lab)
        ax.set_title(title)
        ax.set_xlabel('days after eating-ripe')
        tag(ax, letter)
    axes.flat[0].legend(loc='lower left', fontsize=8.5)
    ax = axes.flat[5]
    names = ['ethanol', 'ethyl_acetate', 'isoamyl_alcohol', 'isoamyl_acetate', 'acetic_acid', 'isobutyl_acetate',
             'butyl_acetate', 'hexanal', 'acetoin', '2-heptanone', '2-phenylethanol']
    qs = {}
    for lab, sc, d in [('ripe', 'intact', 0), ('rotting d7', 'damaged', 7)]:
        c = composition(sc, d)
        q = emission_mol_s(c['volatiles_mg_per_kg'], c['water_L_per_kg'], c['taste']['pH'])
        qs[lab] = [q.get(n, 0) * 1e9 for n in names]
    y = np.arange(len(names))
    ax.barh(y + 0.2, np.maximum(qs['ripe'], 1e-4), 0.38, color=BLUE, label='ripe')
    ax.barh(y - 0.2, np.maximum(qs['rotting d7'], 1e-4), 0.38, color=ORANGE, label='rotting day 7')
    ax.set_yticks(y, [n.replace('_', ' ') for n in names], fontsize=8)
    ax.invert_yaxis()
    ax.set_xscale('log')
    ax.set_xlim(1e-2, 1e4)
    ax.set_title('Emission from 30 cm² of pulp (nmol/s)')
    ax.legend(loc='lower right', fontsize=8.5)
    tag(ax, 'F')
    fig.tight_layout()
    fig.savefig('paper/fig2_chemistry.png', dpi=160)


def fig_atmosphere(exp):
    fig = plt.figure(figsize=(11, 6.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.05, 1], hspace=0.45, wspace=0.45)
    ax = fig.add_subplot(gs[0, :])
    p = Plume([(0.0, 0.15), (0.0, -0.15)], x_max=2.4, seed=21)
    for _ in range(int(25 / 0.05)):
        p.step(0.05)
    xs, ys = np.linspace(-0.1, 2.3, 600), np.linspace(-1.0, 1.0, 330)
    field = np.zeros((len(ys), len(xs), 3))
    for k, col in enumerate([(0.16, 0.47, 0.84), (0.92, 0.41, 0.2)]):
        c = p.c_unit_grid(xs, ys, source=k)
        a = np.clip(np.log10(c + 1e-9) / 4.0, 0, 1) ** 1.4
        field += a[..., None] * np.array(col)[None, None]
    bg = np.array([0.985, 0.99, 0.99])
    img = bg * (1 - np.clip(field.sum(-1, keepdims=True), 0, 1)) + field
    ax.imshow(np.clip(img, 0, 1), extent=[xs[0], xs[-1], ys[0], ys[-1]], origin='lower', aspect='auto')
    ax.scatter([0, 0], [0.15, -0.15], s=60, c=[BLUE, ORANGE], edgecolors=INK, zorder=3)
    ax.text(0.04, 0.2, 'ripe', color=INK, fontsize=9)
    ax.text(0.04, -0.25, 'rotting day 7', color=INK, fontsize=9)
    ax.set_title('Instantaneous odour field (wind 0.3 m/s blowing left to right; log concentration)')
    ax.set_xlabel('downwind distance (m)')
    ax.set_ylabel('crosswind (m)')
    ax.grid(False)
    tag(ax, 'A')
    d = np.load('results/plume_series.npz')
    S, dists, dt = d['series'], d['distances'], float(d['dt'])
    ax = fig.add_subplot(gs[1, 0])
    tt = np.arange(len(S)) * dt
    sel = tt < 30
    for j, col in zip([1, 2, 3], [INK, AQUA, MAGENTA]):
        ax.plot(tt[sel], S[sel, j] / S[:, j].mean(), color=col, linewidth=1.2, label=f'{dists[j]:g} m')
    ax.set_yscale('symlog', linthresh=0.1)
    ax.set_title('Concentration ÷ mean at fixed points')
    ax.set_xlabel('time (s)')
    ax.legend(fontsize=8.5, loc='upper right', ncol=3)
    tag(ax, 'B')
    ax = fig.add_subplot(gs[1, 1])
    cal = exp['plume_calibration']
    x = np.arange(len(cal))
    ax.bar(x - 0.2, [c['intermittency'] for c in cal], 0.38, color=AQUA, label='intermittency (fraction of time > 10% of mean)')
    ax.bar(x + 0.2, [c['ratio_to_analytic'] for c in cal], 0.38, color=INK2, label='mean ÷ Gaussian-plume mean')
    ax.set_xticks(x, [f"{c['distance_m']:g} m" for c in cal])
    ax.set_ylim(0, 1.5)
    ax.set_title('Plume statistics')
    ax.legend(fontsize=8, loc='upper left')
    tag(ax, 'C')
    ax = fig.add_subplot(gs[1, 2])
    dr = json.load(open('results/detection_range.json'))
    rows = [r for r in dr['rows'] if r['temp'] == 25 and r['wind'] == 0.3 and r['spread'] == 0.1]
    combos = [('ripe', 'fed'), ('ripe', 'starved_24h'), ('rotting', 'fed'), ('rotting', 'starved_24h')]
    for i, (food, st) in enumerate(combos):
        m = next(r['range_m'] for r in rows if r['food'] == food and r['state'] == st and r['basis'] == 'mean')
        f = next(r['range_m'] for r in rows if r['food'] == food and r['state'] == st and r['basis'] == 'filament_p95')
        col = BLUE if food == 'ripe' else ORANGE
        ax.barh(i, f, color=col, alpha=0.35, zorder=2)
        ax.barh(i, m, color=col, zorder=3)
        ax.text(f + 0.4, i, f'{m:.1f}–{f:.1f} m', va='center', fontsize=8.5, color=INK)
    ax.set_yticks(range(4), [f"{f}, {'fed' if s == 'fed' else 'starved'}" for f, s in combos])
    ax.invert_yaxis()
    ax.set_xlim(0, 30)
    ax.set_title('Detection range (solid: mean, pale: filaments)')
    ax.set_xlabel('distance downwind (m), 0.3 m/s')
    tag(ax, 'D')
    fig.savefig('paper/fig3_atmosphere.png', dpi=160, bbox_inches='tight')


def fig_brain():
    d = np.load('results/olf_lookup.npz', allow_pickle=True)
    probe = json.load(open('results/probe_bilateral.json'))
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), gridspec_kw={'width_ratios': [1.4, 1, 1]})
    ax = axes[0]
    lv = d['levels']
    for food, col in [('ripe', BLUE), ('rotting', ORANGE)]:
        for state, ls in [('fed', '-'), ('starved_24h', '--')]:
            r = d[f'resp_{food}_25_{state}'][:, :, 0]
            ax.plot(lv, r.mean(1), color=col, linestyle=ls, marker='o', markersize=3.5,
                    label=f"{food}, {'fed' if state == 'fed' else 'starved'}")
    base = d['baseline_fed'][:, 0]
    ax.axhline(base.mean() + 3 * base.std(), color=INK2, linewidth=1)
    ax.text(lv[0], base.mean() + 3 * base.std() * 1.25, 'detection threshold (baseline + 3 SD)', color=INK2, fontsize=8)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('plume concentration per unit emission (s/m³)')
    ax.set_title('Projection-neuron spikes in 100 ms')
    ax.legend(fontsize=8, loc='lower right')
    tag(ax, 'A')
    ax = axes[1]
    lat = probe['latency_ms']
    ax.barh(['ORN→PN', 'ORN→lateral horn', 'ORN→descending'], [lat['PN'], lat['LH'], lat['DN']], color=[AQUA, AQUA, INK2])
    for i, v in enumerate([lat['PN'], lat['LH'], lat['DN']]):
        ax.text(v + 0.5, i, f'{v} ms', va='center', color=INK, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, 30)
    ax.set_title('Latency after ORN onset')
    tag(ax, 'B')
    ax = axes[2]
    bl = probe['bilateral']
    labels = [b['stim'] for b in bl]
    ax.bar(np.arange(len(bl)) - 0.2, [b['steer_index'] for b in bl], 0.38, color=INK2, label='DNa01/DNa02 (turning)')
    ax.bar(np.arange(len(bl)) + 0.2, [b['dn_index'] for b in bl], 0.38, color=MAGENTA, label='all descending neurons')
    ax.axhline(0, color=INK, linewidth=0.8)
    ax.set_xticks(np.arange(len(bl)), [l.replace(' ', '\n') for l in labels], fontsize=8)
    ax.set_ylim(-1.05, 1.05)
    ax.set_title('Left − right output index')
    ax.legend(fontsize=8, loc='lower left')
    tag(ax, 'C')
    fig.tight_layout()
    fig.savefig('paper/fig4_brain.png', dpi=160)


def fig_taste():
    cal = json.load(open('results/feeding_calibration.json'))
    g = pd.DataFrame(json.load(open('results/grid_taste.json')))
    fig = plt.figure(figsize=(11, 3.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.6, 1.6], wspace=0.3)
    ax = fig.add_subplot(gs[0])
    ax.plot(cal['conc_mM'], cal['accept_fed'], color=BLUE, marker='o', markersize=4, label='fed')
    ax.plot(cal['conc_mM'], cal['accept_starved'], color=ORANGE, marker='o', markersize=4, label='starved 24 h')
    ax.axhline(cal['A50'], color=INK2, linewidth=1)
    ax.text(10, cal['A50'] + 0.4, 'PER 55% (calibrated at 800 mM, fed)', color=INK2, fontsize=8)
    ax.set_xscale('log')
    ax.set_xlabel('sucrose (mM)')
    ax.set_title('Feeding-neuron spikes, sucrose')
    ax.legend(fontsize=8, loc='lower right')
    tag(ax, 'A')
    foods = ['ripe (day 0)', 'overripe (day 5)', 'rotting day 2', 'rotting day 4', 'rotting day 7', 'rotting day 14']
    for k, (enrich, title) in enumerate([(1.0, 'Pulp (yeast averaged through fruit)'), (20.0, 'Surface yeast colonies (×20)')]):
        ax = fig.add_subplot(gs[k + 1])
        for si, st in enumerate(STATE_COL):
            vals = []
            for food in foods:
                e = 1.0 if food.startswith(('ripe', 'overripe')) else enrich
                vals.append(g[(g.food == food) & (g.state == st) & (g.yeast_enrichment == e)].accept.mean())
            ax.bar(np.arange(len(foods)) + (si - 1.5) * 0.2, vals, 0.17, color=STATE_COL[st], label=STATE_LAB[st], zorder=3)
        ax.set_xticks(np.arange(len(foods)), [f.replace(' (', '\n(').replace('rotting ', 'rotting\n') for f in foods], fontsize=8)
        ax.set_title(title)
        if k == 0:
            ax.legend(fontsize=7.5, loc='upper right')
        tag(ax, 'BC'[k])
    fig.savefig('paper/fig5_taste.png', dpi=160, bbox_inches='tight')


def fig_behaviour(exp):
    tracks = json.load(open('results/tracks.json'))
    fig = plt.figure(figsize=(11, 7.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], hspace=0.45, wspace=0.32)
    ax = fig.add_subplot(gs[0, :2])
    mode_col = {0: '#a9b3b0', 1: AQUA, 2: YELLOW, 3: '#cccccc'}
    for tr in tracks[:30]:
        a = np.array(tr)
        if len(a) < 2:
            continue
        for m in [0, 1]:
            seg = np.where(a[:, 3] == m, a[:, 2], np.nan)
            ax.plot(a[:, 1], seg, color=mode_col[m], linewidth=0.8, alpha=0.9)
    ax.scatter([0, 0], [0.15, -0.15], s=90, c=[BLUE, ORANGE], edgecolors=INK, zorder=4)
    ax.plot([], [], color=mode_col[0], label='casting')
    ax.plot([], [], color=mode_col[1], label='surging (odour detected)')
    ax.set_xlim(-0.1, 2.2)
    ax.set_ylim(-0.8, 0.8)
    ax.set_xlabel('downwind distance (m)')
    ax.set_title('30 starved flies, released at 2 m; ripe (blue) vs rotting (orange)')
    ax.legend(fontsize=8.5, loc='upper right')
    tag(ax, 'A')
    ax = fig.add_subplot(gs[0, 2])
    ss = exp['single_source']
    labels, vals, cols = [], [], []
    for food in ['ripe', 'rotting']:
        for st in ['fed', 'starved_24h']:
            r = ss[f'{food}_{st}']
            labels.append(f"{food}\n{'fed' if st == 'fed' else 'starved'}")
            vals.append(100 * (r['n'] - r['never_landed']) / r['n'])
            cols.append(STATE_COL[st])
    ax.bar(np.arange(4), vals, color=cols, zorder=3)
    for i, v in enumerate(vals):
        ax.text(i, v + 1.5, f'{v:.0f}%', ha='center', fontsize=8.5, color=INK)
    ax.set_xticks(np.arange(4), labels, fontsize=8)
    ax.set_ylim(0, 110)
    ax.set_title('Single source: % reaching fruit')
    tag(ax, 'B')
    tc = exp['two_choice']
    for k, (metric, title) in enumerate([('first_landing', 'First landing on rotting (%)'), ('fed', 'Feeding on rotting (% of feeders)')]):
        ax = fig.add_subplot(gs[1, k])
        variants = [('yeast20|detect|logistic', 'logistic'), ('yeast20|detect|hyperbolic', 'hyperbolic'),
                    ('yeast20|gated|hyperbolic', 'valence-gated'), ('yeast1|detect|hyperbolic', 'no surface yeast')]
        for vi, (suffix, vlab) in enumerate(variants):
            vals = []
            for st in STATE_COL:
                r = tc[f'{st}|{suffix}']
                tot = sum(r[metric].values())
                vals.append(100 * r[metric]['rotting'] / tot if tot else np.nan)
            ax.plot(np.arange(4), vals, marker='o', markersize=5, color=[INK, AQUA, MAGENTA, INK2][vi],
                    linestyle=['-', '-', '--', ':'][vi], label=vlab)
        ax.axhline(50, color=INK2, linewidth=0.8)
        ax.set_xticks(np.arange(4), [{'fed': 'fed', 'starved_24h': 'starved', 'protein_deprived': 'protein-\ndeprived', 'starved_and_protein_deprived': 'both'}[s] for s in STATE_COL], fontsize=8)
        ax.set_ylim(0, 100)
        ax.set_title(title)
        if k == 1:
            ax.legend(fontsize=7.5, loc='lower right')
        tag(ax, 'CD'[k])
    ax = fig.add_subplot(gs[1, 2])
    ff = exp['far_field_plume']
    conds = [('0.3|2.0', '2 m,\n30 cm apart'), ('0.3|10.0', '10 m,\n30 cm apart'), ('4.0|10.0', '10 m,\n4 m apart')]
    for si, (st, col) in enumerate([('fed', BLUE), ('starved_24h', ORANGE)]):
        vals = [100 * ff[f'{c}|{st}']['rotting_share'] for c, _ in conds]
        sds = [100 * ff[f'{c}|{st}']['share_sd_across_plumes'] for c, _ in conds]
        ax.bar(np.arange(3) + (si - 0.5) * 0.36, vals, 0.34, yerr=sds, color=col, label='fed' if st == 'fed' else 'starved 24 h',
               zorder=3, error_kw=dict(ecolor=INK2, lw=1, capsize=2))
    ax.axhline(50, color=INK2, linewidth=0.8)
    ax.set_xticks(range(3), [l for _, l in conds], fontsize=8)
    ax.set_ylim(0, 100)
    ax.set_title('Diluting far-field plume: first landing on rotting (%)')
    ax.legend(fontsize=8, loc='upper left')
    tag(ax, 'E')
    fig.savefig('paper/fig6_behaviour.png', dpi=160, bbox_inches='tight')


if __name__ == '__main__':
    exp = json.load(open('results/experiments.json'))
    fig_chemistry()
    fig_atmosphere(exp)
    fig_brain()
    fig_taste()
    fig_behaviour(exp)
    print('figures written')
