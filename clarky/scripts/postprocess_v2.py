#!/usr/bin/env python3
"""Summarise a Clark-Y cavitation run: vapour volume, inlet sigma, Cl/Cd, cavity length,
shedding frequency (FFT of vapour volume and Cl), cavity area relative to foil area Sc.

v2 vs v1: adds FFT frequency / Strouhal and Vvap/Sc; defaults match Huang 2013 water (25 C).

Usage: postprocess_v2.py <caseDir> [--tavg 0.03]   (averages over t >= tavg)
Writes <caseDir>/summary.txt and <caseDir>/history.png
"""
import argparse, glob, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

p = argparse.ArgumentParser()
p.add_argument("case")
p.add_argument("--tavg", type=float, default=None, help="start of averaging window [s]")
p.add_argument("--c", type=float, default=0.07)
p.add_argument("--dz", type=float, default=0.001)
p.add_argument("--aoa", type=float, default=8.0)
p.add_argument("--U", type=float, default=10.0)
p.add_argument("--rho", type=float, default=999.19)
p.add_argument("--pSat", type=float, default=3169.0)
p.add_argument("--Sc", type=float, default=0.0809, help="foil area / c^2")
p.add_argument("--alphaCav", type=float, default=0.9, help="alpha.water threshold for cavity")
a = p.parse_args()
pp = os.path.join(a.case, "postProcessing")
q = 0.5 * a.rho * a.U ** 2

def load(pattern):
    files = sorted(glob.glob(os.path.join(pp, pattern)), key=lambda f: float(f.split("/")[-2]))
    d = np.vstack([np.loadtxt(f, comments="#", ndmin=2) for f in files])
    _, i = np.unique(d[:, 0], return_index=True)      # drop restart overlaps
    return d[i]

# vapour volume from integral of alpha.water and domain volume in header
cf = sorted(glob.glob(os.path.join(pp, "cavity/*/volFieldValue.dat")))[0]
Vtot = float([l for l in open(cf) if l.startswith("# Volume")][0].split(":")[1])
cav = load("cavity/*/volFieldValue.dat")
t_v, Vv = cav[:, 0], (Vtot - cav[:, 1]) / (a.c ** 2 * a.dz)   # vapour area / c^2

pin = load("inletP/*/surfaceFieldValue.dat")
t_p, sig = pin[:, 0], (pin[:, 1] - a.pSat) / q

fc = load("forces/*/coefficient.dat")
t_f, Cd, Cl = fc[:, 0], fc[:, 1], fc[:, 4]

# cavity length on suction side from wall alpha (chordwise coordinate in foil frame)
th = np.radians(a.aoa)
def chordwise(x, y):
    xr, yr = x - 0.5 * a.c, y
    return (xr * np.cos(th) - yr * np.sin(th) + 0.5 * a.c) / a.c, (xr * np.sin(th) + yr * np.cos(th)) / a.c

Lc = []
for d in sorted(glob.glob(os.path.join(pp, "foilAlpha/*")), key=lambda s: float(os.path.basename(s))):
    f = os.path.join(d, "alpha.water_foil.raw")
    if not os.path.exists(f):
        continue
    s = np.loadtxt(f, comments="#")
    xc, yc = chordwise(s[:, 0], s[:, 1])
    up = yc > 0.0
    xs, al = xc[up], s[up, 3]
    o = np.argsort(xs)
    xs, al = xs[o], al[o]
    wet = np.where(al >= a.alphaCav)[0]
    # cavity attached at LE: first wet face after the first vapour face
    vap = np.where(al < a.alphaCav)[0]
    if len(vap) == 0:
        L = 0.0
    else:
        after = wet[wet > vap[0]]
        L = xs[after[0]] if len(after) else 1.0
    Lc.append((float(os.path.basename(d)), L))
Lc = np.array(Lc)

t0 = a.tavg if a.tavg is not None else 0.5 * t_v[-1]
def avg(t, y):
    m = t >= t0
    return y[m].mean(), y[m].std()

lines = [f"case            {os.path.abspath(a.case)}",
         f"t_end           {t_v[-1]:.4f} s  ({t_v[-1] * a.U / a.c:.1f} c/U)",
         f"avg window      t >= {t0:.4f} s",
         f"sigma_inlet     {avg(t_p, sig)[0]:.3f} +- {avg(t_p, sig)[1]:.3f}",
         f"Vvap/(c^2 dz)   {avg(t_v, Vv)[0]:.4f} +- {avg(t_v, Vv)[1]:.4f}",
         f"Cl              {avg(t_f, Cl)[0]:.3f} +- {avg(t_f, Cl)[1]:.3f}",
         f"Cd              {avg(t_f, Cd)[0]:.4f} +- {avg(t_f, Cd)[1]:.4f}"]
def peak_freq(t, y):
    m = t >= t0
    if m.sum() < 64:
        return np.nan
    tt = np.linspace(t[m][0], t[m][-1], 4096)
    yy = np.interp(tt, t[m], y[m]); yy -= yy.mean()
    F = np.abs(np.fft.rfft(yy * np.hanning(len(yy))))
    fr = np.fft.rfftfreq(len(tt), tt[1] - tt[0])
    k = 1 + np.argmax(F[1:])
    return fr[k]
fV, fL = peak_freq(t_v, Vv), peak_freq(t_f, Cl)
lines.append(f"Vvap/Sc         {avg(t_v, Vv)[0] / a.Sc:.3f}   (Huang 2013: 0.72 exp, sigma=0.8)")
lines.append(f"f_peak Vvap     {fV:.1f} Hz  St={fV * a.c / a.U:.3f}")
lines.append(f"f_peak Cl       {fL:.1f} Hz  St={fL * a.c / a.U:.3f}   (window {(t_v[-1] - t0) * a.U / a.c:.0f} c/U)")
if len(Lc):
    m = Lc[:, 0] >= t0
    if m.any():
        lines.append(f"L_cav/c         {Lc[m, 1].mean():.3f} +- {Lc[m, 1].std():.3f}  (alpha<{a.alphaCav}, {m.sum()} samples)")
open(os.path.join(a.case, "summary.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))

fig, ax = plt.subplots(4, 1, figsize=(8, 9), sharex=True)
tc = a.U / a.c
ax[0].plot(t_p * tc, sig); ax[0].set_ylabel(r"$\sigma_{inlet}$")
ax[1].plot(t_v * tc, Vv); ax[1].set_ylabel(r"$V_{vap}/c^2$")
ax[2].plot(t_f * tc, Cl, label="Cl"); ax[2].plot(t_f * tc, Cd, label="Cd"); ax[2].legend(); ax[2].set_ylim(-0.5, 2.5)
if len(Lc):
    ax[3].plot(Lc[:, 0] * tc, Lc[:, 1], "o-")
ax[3].set_ylabel(r"$L_{cav}/c$"); ax[3].set_xlabel(r"$tU/c$")
for x in ax:
    x.axvline(t0 * tc, color="gray", ls="--", lw=0.8); x.grid(alpha=0.3)
plt.tight_layout(); plt.savefig(os.path.join(a.case, "history.png"), dpi=110)
