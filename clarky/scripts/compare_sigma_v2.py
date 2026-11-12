#!/usr/bin/env python3
"""Cl, Cd, cavity area vs sigma_inlet for v5 runs vs DIGITISED Wang (2001) Fig. 16 (8 deg).
v2 vs v1: uses wang2001_fig16_8deg_v1.csv (digitised) instead of eyeballed values; reports the
sigma-shift needed to match CFD Cl to the experimental Cl curve.
Usage (from clarky/): compare_sigma_v2.py -> sigma_sweep_v2.png, sigma_sweep_v2.csv
"""
import glob, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

w = np.loadtxt("wang2001_fig16_8deg_v1.csv", delimiter=",", comments="#", skiprows=5)
rows = []
for f in sorted(glob.glob("runs/v5_reboud_s[0-9]*/summary.txt")):
    if "fine" in f:
        continue
    t = open(f).read()
    g = lambda k: float(re.search(rf"^{re.escape(k)}\s+([-\d.]+)", t, re.M).group(1))
    rows.append((g("sigma_inlet"), g("Cl"), g("Cd"), g("A_cav/Sc")))
rows = np.array(sorted(rows))
np.savetxt("sigma_sweep_v2.csv", rows, delimiter=",", header="sigma_inlet,Cl,Cd,Acav_over_Sc", fmt="%.4f")

print(" sigma  Cl_cfd  Cl_exp   dCl   Cd_cfd  Cd_exp   sigma_exp(same Cl)  shift")
mono = w[w[:, 0] <= 1.2]                       # Cl monotonic in sigma on this branch
for s, cl, cd, a in rows:
    cle, cde = np.interp(s, w[:, 0], w[:, 1]), np.interp(s, w[:, 0], w[:, 2])
    se = np.interp(cl, mono[:, 1], mono[:, 0]) if cl < mono[:, 1].max() else np.nan
    print(f" {s:5.3f}  {cl:6.3f}  {cle:6.3f}  {cl-cle:+6.3f}  {cd:6.4f}  {cde:6.4f}   {se:6.3f}            {s-se:+.3f}")

fig, ax = plt.subplots(1, 3, figsize=(14, 4))
for k, (lab, col) in enumerate([("Cl", 1), ("Cd", 2)]):
    ax[k].plot(w[:, 0], w[:, col], "ks--", mfc="none", label="Wang 2001 (digitised)")
    ax[k].plot(rows[:, 0], rows[:, col], "o-", color="C0", label="CFD v5, 32.6k")
    ax[k].set_xlabel(r"$\sigma$ (upstream)"); ax[k].set_ylabel(lab); ax[k].grid(alpha=0.3); ax[k].legend()
ax[2].plot(rows[:, 0], rows[:, 3], "o-", color="C0", label="CFD v5")
ax[2].plot([0.8], [0.72], "k^", ms=9, label="Huang 2013 exp")
ax[2].set_xlabel(r"$\sigma$"); ax[2].set_ylabel(r"$A_{cav}/S_c$"); ax[2].grid(alpha=0.3); ax[2].legend()
plt.tight_layout(); plt.savefig("sigma_sweep_v2.png", dpi=110)
