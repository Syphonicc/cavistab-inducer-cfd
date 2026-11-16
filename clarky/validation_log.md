# Clark-Y validation log

## 2026-11-12 — base-flow bug found and fixed (template v5)

**Symptom.** In every interPhaseChangeFoam run (incl. both CloudHPC runs and the first
Reboud run), Cl climbed to ~1.9–2.5 and Cd went negative during the fully wetted hold;
inlet pressure ended up *below* outlet pressure; spurious cavitation at σ = 2–3.

**Isolation (all on the identical 32.6k-cell mesh, same BCs, same initial field):**
| Case | Cl, Cd (~3–20 ms) | Verdict |
|---|---|---|
| simpleFoam, steady SST | 1.18, 0.021 | matches Wang 2001 (≈1.16) |
| pimpleFoam, transient SST / laminar | 1.17, 0.021 / 1.18, 0.016 | OK |
| interPhaseChangeFoam, pure liquid (σ=10) | → 1.9, −0.08 | drifts |
| interFoam (no phase change), SST and laminar | → 1.86, −0.03 | drifts → not the cavitation model |
| interFoam: zeroGradient p on foil, 2 non-orth corr., no relaxation, gauge p, tight tolerances, momentumPredictor yes, PISO, thicker wall cells | all identical drift | ruled out |
| interFoam, laplacian/snGrad `uncorrected` | 1.63, −0.23 | worse |
| **interFoam, laplacian/snGrad `corrected`** | **1.183, 0.0136** | **fixed** |

**Cause.** `cavitatingBullet` tutorial schemes use `limited corrected 0.5`. VoF solvers rebuild
cell velocity from face fluxes (`fvc::reconstruct`); with a partially corrected
non-orthogonal flux on the curved LE cells this produced a saw-tooth near-wall jet
(27/22/25 m/s vs a smooth 17.5 m/s in pimpleFoam) behind the LE → excess circulation.
pimpleFoam uses the cell pressure gradient and is unaffected.

**Consequence.** All earlier interPhaseChangeFoam results (cloud_v1, reboud_s0.8_v3,
local_v2 σ=2.0, v4) are invalid. The Reboud model itself was verified independently
(nut ↓ ×12 in mixture, ×5000 in vapour, unchanged in liquid).

**Other fixes along the way.**
- run_case_v4+: σ set with `sed`, not `foamDictionary -set` (that froze `#eval` pOut).
- v4+: hold at σ = 10 (was 3), slower ramp.
- v5: hold 0.05 s, ramp to target by 0.10 s, run to 0.40 s; Cl during hold now rises
  smoothly 0.84 → 1.12 (Wagner-type start-up), Cd 0.056 → 0.028, no vapour.

Steady SST Cp_min = −2.71 → the model will cavitate a little at σ = 2 (experiment
calls σ = 2 subcavitating; real water sustains slight tension).

## 2026-11-12 — first valid results (template v5, Reboud SST, 32.6k cells)
Averages over t = 0.15–0.40 s (36 c/U); experiment values from reference_data.md.

| | σ = 2.0 | Exp | σ = 0.8 | Exp |
|---|---|---|---|---|
| σ at inlet | 2.005 | — | 0.853 | (σ defined at outlet = 0.80) |
| Cl | 1.177 ± 0.05 | ≈ 1.16 | 0.760 ± 0.20 | ≈ 0.75–0.8 |
| Cd | 0.023 | ≈ 0.037 | 0.139 | ≈ 0.11–0.12 |
| A_cav (10 % vapour) / Sc | 0.02 | subcavitating | 0.98 | 0.72 (Huang) |
| Shedding (FFT Vvap, Cl) | — | — | 24 Hz, St 0.168 | 20–26 Hz, St 0.185 |

FFT resolution 4 Hz (0.25 s window). Open issues: Cd ~15–25 % high and cavity ~35 % large at
σ = 0.8; inlet σ 0.85 vs outlet 0.80 (experiment references upstream pressure) — next run
should target σ_inlet; mesh-independence not yet checked; Wang Cl/Cd values are eyeballed.

## 2026-11-12 — σ sweep (v5): 2.0, 1.4, 1.0, 0.8 → `sigma_sweep_v1.png/.csv`
| σ inlet | Cl | Wang Cl* | Cd | Wang Cd* | A_cav/Sc |
|---|---|---|---|---|---|
| 2.005 | 1.177 | 1.16 | 0.023 | 0.037 | 0.02 |
| 1.458 | 1.055 | 1.17 | 0.112 | 0.040 | 0.34 |
| 1.091 | 0.852 | 1.03 | 0.137 | 0.120 | 0.74 |
| 0.853 | 0.760 | 0.78 | 0.139 | 0.115 | 0.98 |
\*eyeballed from Wang Fig. 16.

**Finding:** the CFD Cl/Cd-vs-σ curve is shifted to higher σ by roughly 0.4–0.5: the model
loses lift / gains drag earlier (at σ = 1.4 it has a large unsteady cavity where the
experiment has a thin stable sheet). Consistent with steady Cp_min = −2.71: the equilibrium
mass-transfer model cavitates as soon as p < pSat, while tunnel water shows inception only
around σ ≈ 1.6–2.0. Agreement at σ = 0.8 is therefore partly coincidental.

**Options:** calibrate Schnerr–Sauer (nuclei n, dNuc, Cv/Cc) against the sheet regime; an
effective tension / pSat offset; check mesh dependence first (LE suction peak is resolution-
sensitive). For K = −dV_c/dp the *slope* of cavity volume vs σ matters — the shift must be
understood before extracting K.

## 2026-11-12 — mesh independence + Schnerr–Sauer calibration at σ = 1.4
**Mesh (CloudHPC, `cloud_v5/fine1.5_s1.4`, 73.7k vs 32.6k cells, 0.40 s, 3.9 h on 4 cores):**
Cl 1.055 → 1.074 (+1.8 %), Cd 0.112 → 0.118 (+5 %), A_cav/Sc 0.34 → 0.38, f(Vvap) 16 Hz both.
→ the σ-shift is not a resolution effect; 32.6k mesh adequate for screening.

**Schnerr–Sauer coefficients (32.6k, run to 0.25 s, avg 0.15–0.25 s = 14 c/U):**
| Run | Cl | Cd | A_cav/Sc |
|---|---|---|---|
| baseline n 1.6e13, Cv 1 | 1.055 | 0.112 | 0.34 |
| Cv 0.1 | 1.086 | 0.121 | 0.40 |
| Cv 0.01 | 1.100 | 0.092 | 0.28 |
| n 1e11 | 1.039 | 0.098 | 0.27 |
| n 1e9 | 1.045 | 0.106 | 0.29 |
| Wang (eyeballed) | ≈ 1.17 | ≈ 0.04 | thin stable sheet |

→ mass-transfer *rate* has little effect (Cd at best 0.09). The error is the inception
threshold (vapour wherever p < pSat; steady Cp_min = −2.71), not the vaporisation rate.
Next candidates: effective tension (pSat offset) / inception criterion; check Wang σ
reference and digitise Fig. 16 before tuning further.

Laptop temps during 4–8-core runs (536 samples): CPU max 87.9 °C, mean 75.8 °C (`runs/temps.csv`).
Earlier run loss at 13:40 was a system suspend, not thermal; runs now wrapped in `systemd-inhibit`.

## 2026-11-13 — Wang Fig. 16 digitised; σ reference checked
- σ reference: Huang 2013 (same group/tunnel) — σ uses the *upstream tunnel pressure*. CFD
  comparisons already use σ_inlet → reference choice explains ≤ 0.05 of the shift.
- Digitised 8° data: `clarky/wang2001_fig16_8deg_v1.csv` (±0.01 Cl, ±0.002 Cd).
  Comparison `scripts/compare_sigma_v2.py` → `sigma_sweep_v2.png/.csv`:

| σ_in | Cl CFD | Cl exp | ΔCl | Cd CFD | Cd exp | σ-shift to match Cl |
|---|---|---|---|---|---|---|
| 0.853 | 0.760 | 0.850 | −0.09 | 0.139 | 0.116 | +0.09 |
| 1.091 | 0.852 | 1.093 | −0.24 | 0.137 | 0.106 | +0.24 |
| 1.458 | 1.055 | 1.166 | −0.11 | 0.112 | 0.040 | +0.41 |
| 2.005 | 1.177 | 1.152 | +0.03 | 0.023 | 0.037 | (flat branch) |

**Revised finding:** σ = 0.8 "agreement" was an artefact of the eyeballed reference (CFD Cl is
0.09 low). The discrepancy is not a constant σ-offset: the experiment keeps Cl ≈ 1.17 (stable
sheet, no lift penalty) down to σ ≈ 1.2 then drops steeply; CFD loses lift gradually from σ ≈ 2.
The model lacks the stable-sheet regime → a single tension/pSat offset cannot fix it across σ.
Candidates: sheet stability / closure (Reboud exponent n, re-entrant jet), other mass-transfer
models (Zwart, Kunz, Merkle), and the inception threshold — test at σ = 1.4 where error is largest.

## 2026-11-15 — model comparison at σ = 1.4 (`run_case_v9.sh`, 2 cores each, avg 0.15–0.30 s)
| Run | σ_in | Cl | Cd | A_cav/Sc |
|---|---|---|---|---|
| baseline Schnerr–Sauer, Reboud n = 10 | 1.458 | 1.055 | 0.112 | 0.34 |
| Kunz (Cc = Cv = 1000, U∞ 10, t∞ 0.007) | 1.485 | 1.107 | 0.126 | 0.47 |
| Reboud n = 3 | 1.436 | 1.024 | 0.100 | 0.26 |
| Merkle (Cc 80, Cv 1e-3) | — | failed: Δt collapsed to 1e-13 s, 14 ms in 18 h | | |
| tension pSat_eff = −16.8 kPa | — | failed: FPE at step 1 (SchnerrSauer uses 0.01·pSat as a floor → invalid for pSat < 0) | | |
| Wang (digitised) | 1.46 | 1.166 | 0.040 | thin stable sheet |

**Go/no-go (Cd < 0.06 and Cl > 1.12): no model passes.** Clark-Y calibration closed.
**Documented limitation:** with interPhaseChangeFoam + k-ω SST (Reboud) and any of the tested
mass-transfer models/coefficients, the stable partial-sheet regime (σ ≈ 1.2–1.6 at 8°) is not
reproduced: cavity too large/unsteady, Cd 2.5–3× too high, Cl 5–12 % low. Inception-threshold
handling (effective tension) remains untested (needs a model that accepts an offset, e.g. a
modified Schnerr–Sauer with p − (pSat − Δp) in the driving term).
Max CPU temperature during the 8-core run: 95.5 °C (alert threshold 97 °C) → use ≤ 6 cores
or better cooling for long runs.

## 2026-11-16 — single-foil quasi-steady compliance from the v5 σ sweep (check "C", no new runs)
| σ_in | V_vap/c² (± std) | A(10 % vapour)/c² |
|---|---|---|
| 0.853 | 0.0590 (0.0144) | 0.0790 |
| 1.091 | 0.0452 (0.0217) | 0.0597 |
| 1.458 | 0.0213 (0.0142) | 0.0275 |
| 2.005 | 0.0010 (0.0003) | 0.0015 |
K_foil = −∂(V_c/c²)/∂σ (finite differences, mid-σ 1.73 / 1.27 / 0.97): 0.037 / 0.065 / 0.058 (V_vap);
0.047 / 0.088 / 0.081 (10 % area). Sign correct (V_c grows monotonically as σ falls), order of
magnitude plausible; low-σ trend unresolved — temporal std of V_c (shedding) is as large as the
differences. → quasi-steady K needs a steady cavity; no kill condition triggered.
