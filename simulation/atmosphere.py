"""Turbulent odour plume: filament (puff) model after Farrell et al. 2002 (Environ Fluid Mech 2:143).

Each source releases puffs at a fixed rate; every puff carries an equal share of the emission.
Puffs are advected by the mean wind plus a large-scale meander (the whole plume swings; Ornstein-Uhlenbeck
crosswind velocity) and a small-scale random velocity (relative dispersion), and grow in size.
The concentration at a point is the sum of 3D Gaussian puffs at the source height.

Concentrations are returned per unit emission (s/m3), so one simulation serves every compound:
C_i(x, t) = Q_i * c_unit(x, t).
"""

import numpy as np


class Plume:
    def __init__(self, sources, wind=0.3, meander_sd=0.15, meander_tau=2.0, eddy_sd=0.10, eddy_tau=0.2,
                 sigma0=0.005, growth=0.01, release_hz=25.0, x_max=3.0, seed=0):
        """sources: list of (x, y) in metres; wind blows toward +x.
        meander_sd: crosswind velocity s.d. of the plume-scale meander (m/s)
        eddy_sd: per-puff random velocity s.d. (m/s); growth: puff sigma growth rate as a fraction of wind speed."""
        self.src = np.asarray(sources, float)
        # turbulent velocities scale with the mean wind (constant turbulence intensity); the defaults are
        # the values calibrated at 0.3 m/s
        k = wind / 0.3
        self.U, self.msd, self.mtau = wind, meander_sd * k, meander_tau
        self.esd, self.etau = eddy_sd * k, eddy_tau
        self.s0, self.g = sigma0, growth
        self.rate, self.xmax = release_hz, x_max
        self.rng = np.random.default_rng(seed)
        self.t = 0.0
        self.v_meander = 0.0
        self.acc = np.zeros(len(self.src))
        self.pos = np.zeros((0, 2))
        self.vel = np.zeros((0, 2))
        self.age = np.zeros(0)
        self.sid = np.zeros(0, int)

    def step(self, dt):
        # large-scale meander (spatially uniform, time-correlated)
        a = np.exp(-dt / self.mtau)
        self.v_meander = a * self.v_meander + self.msd * np.sqrt(1 - a * a) * self.rng.standard_normal()
        # release
        self.acc += self.rate * dt
        new = np.floor(self.acc).astype(int)
        self.acc -= new
        if new.sum():
            sid = np.repeat(np.arange(len(self.src)), new)
            self.pos = np.vstack([self.pos, self.src[sid]])
            self.vel = np.vstack([self.vel, np.zeros((len(sid), 2))])
            self.age = np.concatenate([self.age, self.rng.uniform(0, dt, len(sid))])
            self.sid = np.concatenate([self.sid, sid])
        # small-scale velocity (per puff OU)
        b = np.exp(-dt / self.etau)
        self.vel = b * self.vel + self.esd * np.sqrt(1 - b * b) * self.rng.standard_normal(self.vel.shape)
        self.pos[:, 0] += (self.U + self.vel[:, 0]) * dt
        self.pos[:, 1] += (self.v_meander + self.vel[:, 1]) * dt
        self.age += dt
        keep = self.pos[:, 0] < self.xmax
        self.pos, self.vel, self.age, self.sid = self.pos[keep], self.vel[keep], self.age[keep], self.sid[keep]
        self.t += dt

    def sigma(self):
        return self.s0 + self.g * self.U * self.age

    def c_unit(self, pts):
        """Concentration per unit emission (s/m3) at points (N, 2) from each source -> array (N, n_sources)."""
        pts = np.atleast_2d(pts)
        out = np.zeros((len(pts), len(self.src)))
        if len(self.pos) == 0:
            return out
        s = self.sigma()
        m = 1.0 / self.rate  # each puff carries 1/rate seconds of a 1 mol/s emission
        d2 = ((pts[:, None, :] - self.pos[None, :, :]) ** 2).sum(-1)
        contrib = m / ((2 * np.pi) ** 1.5 * s ** 3)[None, :] * np.exp(-d2 / (2 * s * s)[None, :])
        for k in range(len(self.src)):
            out[:, k] = contrib[:, self.sid == k].sum(1)
        return out

    def c_unit_grid(self, xs, ys, source=None):
        X, Y = np.meshgrid(xs, ys)
        pts = np.c_[X.ravel(), Y.ravel()]
        c = self.c_unit(pts)
        c = c.sum(1) if source is None else c[:, source]
        return c.reshape(X.shape)


def analytic_mean(r, wind=0.3):
    """Time-averaged Gaussian plume centreline per unit emission (physics.py), for calibration."""
    s = 0.1 * r + 0.005
    return 1.0 / (np.pi * wind * s * s)


def calibrate(duration=120.0, dt=0.02, distances=(0.1, 0.3, 1.0, 2.0), seed=0):
    """Time series at fixed downwind points: mean vs analytic, intermittency, peak/mean ratio."""
    p = Plume([(0.0, 0.0)], seed=seed)
    pts = np.array([(d, 0.0) for d in distances])
    series = []
    for i in range(int(duration / dt)):
        p.step(dt)
        if p.t > 15:  # skip spin-up
            series.append(p.c_unit(pts)[:, 0])
    S = np.array(series)
    out = []
    for j, d in enumerate(distances):
        x = S[:, j]
        mean = x.mean()
        thr = 0.1 * mean
        out.append(dict(distance_m=d, mean=mean, analytic=analytic_mean(d), ratio_to_analytic=mean / analytic_mean(d),
                        intermittency=float((x > thr).mean()), p95_over_mean=float(np.percentile(x, 95) / mean),
                        max_over_mean=float(x.max() / mean)))
    return out, S


if __name__ == '__main__':
    res, _ = calibrate()
    for r in res:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
