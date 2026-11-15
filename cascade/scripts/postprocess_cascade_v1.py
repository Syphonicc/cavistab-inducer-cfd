#!/usr/bin/env python3
"""Cascade run summary: inlet sigma, cavity area V_c/d^2 (vapour integral and alpha<0.9 region),
cavity length on the suction side (bladeSS, blade-aligned so x = distance from LE), steadiness.

Steadiness metric: std/mean of V_c over the averaging window (quasi-steady K, M need a steady cavity;
kill condition 2 if partial cavities never settle).
Usage: postprocess_cascade_v1.py <caseDir> [--tavg 0.08] [--d 0.05] [--dz 0.001]
Writes <caseDir>/summary.txt and <caseDir>/history.png
"""
import argparse, glob, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

p = argparse.ArgumentParser()
p.add_argument("case")
p.add_argument("--tavg", type=float, default=None)
p.add_argument("--d", type=float, default=0.05)
p.add_argument("--dz", type=float, default=0.001)
p.add_argument("--U", type=float, default=10.0)
p.add_argument("--rho", type=float, default=999.19)
p.add_argument("--pSat", type=float, default=3169.0)
p.add_argument("--alphaCav", type=float, default=0.9)
a = p.parse_args()
pp = os.path.join(a.case, "postProcessing")
q = 0.5 * a.rho * a.U ** 2
A = a.d ** 2 * a.dz                                    # normalising "area" (d^2 per unit span)

def load(pattern):
    files = sorted(glob.glob(os.path.join(pp, pattern)), key=lambda f: float(f.split("/")[-2]))
    dd = np.vstack([np.loadtxt(f, comments="#", ndmin=2) for f in files])
    _, i = np.unique(dd[:, 0], return_index=True)
    return dd[i]

cf = sorted(glob.glob(os.path.join(pp, "cavity/*/volFieldValue.dat")))[0]
Vtot = float([l for l in open(cf) if l.startswith("# Volume")][0].split(":")[1])
cav = load("cavity/*/volFieldValue.dat");      tV, Vv = cav[:, 0], (Vtot - cav[:, 1]) / A
ca = load("cavityArea/*/volFieldValue.dat");   tA, Ac = ca[:, 0], ca[:, 1] / A
pin = load("inletP/*/surfaceFieldValue.dat");  tp, sig = pin[:, 0], (pin[:, 1] - a.pSat) / q

Lc = []
for dd in sorted(glob.glob(os.path.join(pp, "foilAlpha/*")), key=lambda s: float(os.path.basename(s))):
    f = os.path.join(dd, "alpha.water_foil.raw")
    if not os.path.exists(f):
        continue
    s = np.loadtxt(f, comments="#")
    o = np.argsort(s[:, 0]); x, al = s[o, 0], s[o, 3]
    vap = np.where(al < a.alphaCav)[0]
    if len(vap) == 0:
        L = 0.0
    else:
        wet = np.where(al >= a.alphaCav)[0]; after = wet[wet > vap[0]]
        L = x[after[0]] if len(after) else x[-1]
    Lc.append((float(os.path.basename(dd)), L / a.d))
Lc = np.array(Lc)

t0 = a.tavg if a.tavg is not None else 0.5 * tV[-1]
def st(t, y):
    m = t >= t0
    return y[m].mean(), y[m].std()
mV, sV = st(tV, Vv); mA, sA = st(tA, Ac); ms, ss = st(tp, sig)
lines = [f"case             {os.path.abspath(a.case)}",
         f"t_end            {tV[-1]:.4f} s   avg window t >= {t0:.4f} s",
         f"sigma_inlet      {ms:.4f} +- {ss:.4f}",
         f"Vvap/d2          {mV:.5f} +- {sV:.5f}   (std/mean {sV / max(mV, 1e-12):.3f})",
         f"Acav/d2          {mA:.5f} +- {sA:.5f}   (alpha<{a.alphaCav}; std/mean {sA / max(mA, 1e-12):.3f})"]
if len(Lc):
    m = Lc[:, 0] >= t0
    if m.any():
        lines.append(f"Lcav/d           {Lc[m, 1].mean():.3f} +- {Lc[m, 1].std():.3f}   ({m.sum()} samples)")
open(os.path.join(a.case, "summary.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))

fig, ax = plt.subplots(4, 1, figsize=(8, 9), sharex=True)
ax[0].plot(tp, sig); ax[0].set_ylabel(r"$\sigma_{inlet}$"); ax[0].set_ylim(-0.1, max(1.0, 1.2 * ms))
ax[1].plot(tV, Vv, label="vapour integral"); ax[1].plot(tA, Ac, label=r"$\alpha<0.9$ area"); ax[1].legend()
ax[1].set_ylabel(r"$V_c/d^2$")
if len(Lc):
    ax[2].plot(Lc[:, 0], Lc[:, 1], "o-", ms=3)
ax[2].set_ylabel(r"$L_{cav}/d$ (SS)")
fc = load("forces/*/coefficient.dat"); ax[3].plot(fc[:, 0], fc[:, 4], label="Cl"); ax[3].plot(fc[:, 0], fc[:, 1], label="Cd")
ax[3].legend(); ax[3].set_xlabel("t [s]")
for x in ax:
    x.axvline(t0, color="gray", ls="--", lw=0.8); x.grid(alpha=0.3)
plt.tight_layout(); plt.savefig(os.path.join(a.case, "history.png"), dpi=100)
