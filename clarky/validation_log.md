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
