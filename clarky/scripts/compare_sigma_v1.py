#!/usr/bin/env python3
"""Cl, Cd, cavity area vs sigma (inlet) for v5 runs vs Wang et al. (2001) Fig. 16 (8 deg).
Experimental points are EYEBALLED from the figure -- replace with digitized values.
Usage: compare_sigma_v1.py  (run from clarky/) -> sigma_sweep_v1.png, sigma_sweep_v1.csv
"""
import glob, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Wang 2001 Fig. 16, alpha = 8 deg, eyeballed (sigma, Cl, Cd)
wang = np.array([[3.0, 1.16, 0.037], [2.0, 1.16, 0.037], [1.4, 1.17, 0.040], [1.2, 1.18, 0.060],
                 [1.0, 1.03, 0.120], [0.8, 0.78, 0.115], [0.55, 0.57, 0.100], [0.28, 0.30, 0.050]])

rows = []
for f in sorted(glob.glob("runs/v5_*/summary.txt")):
    t = open(f).read()
    g = lambda k: float(re.search(rf"^{k}\s+([-\d.]+)", t, re.M).group(1))
    rows.append((g("sigma_inlet"), g("Cl"), g("Cd"), g("A_cav/Sc"), f.split("/")[1]))
rows.sort()
np.savetxt("sigma_sweep_v1.csv", np.array([r[:4] for r in rows]), delimiter=",",
           header="sigma_inlet,Cl,Cd,Acav_over_Sc", fmt="%.4f")

fig, ax = plt.subplots(1, 3, figsize=(14, 4))
s = np.array([r[0] for r in rows])
for k, (lab, col) in enumerate([("Cl", 1), ("Cd", 2)]):
    ax[k].plot(wang[:, 0], wang[:, col], "ks--", mfc="none", label="Wang 2001 (eyeballed)")
    ax[k].plot(s, [r[col] for r in rows], "o-", color="C0", label="CFD v5 (Reboud SST)")
    ax[k].set_xlabel(r"$\sigma$"); ax[k].set_ylabel(lab); ax[k].grid(alpha=0.3); ax[k].legend()
ax[2].plot(s, [r[3] for r in rows], "o-", color="C0", label="CFD v5")
ax[2].plot([0.8], [0.72], "k^", ms=9, label="Huang 2013 exp")
ax[2].set_xlabel(r"$\sigma$"); ax[2].set_ylabel(r"$A_{cav}/S_c$ (10 % vapour)"); ax[2].grid(alpha=0.3); ax[2].legend()
plt.tight_layout(); plt.savefig("sigma_sweep_v1.png", dpi=110)
for r in rows:
    print(f"{r[4]:18s} sigma={r[0]:.3f} Cl={r[1]:.3f} Cd={r[2]:.4f} Acav/Sc={r[3]:.3f}")
