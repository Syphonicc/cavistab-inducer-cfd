# cavistab — handover (state as of 2026-11-17)

> **Update 2026-11-17:** the owner invoked kill condition 2; the OpenFOAM cascade route is closed.
> Section 6 below is the pre-decision state, kept for the record. See `docs/FINAL_ASSESSMENT.md`.

Repo: github.com/Syphonicc/cavistab-inducer-cfd
Toolchain: OpenFOAM v2512 (ESI), Python 3 (numpy, scipy, matplotlib), CloudHPC (ESI v2512).

---

## 1. Idea and motivation

**Goal:** instability-aware design of rocket turbopump inducers. Inducers are normally optimised for
suction performance (NPSH), but the failures that matter in test are cavitation instabilities
(rotating cavitation, cavitation surge). Their onset is governed by two cavity properties:

- **cavitation compliance** K = −∂(V_c)/∂p (cavity volume vs inlet pressure / σ)
- **mass-flow gain factor** M = −∂(V_c)/∂(flow rate or incidence)

**Plan:** CFD-trained surrogate geometry → (K, M) → 1-D stability model (Tsujimoto-type) →
optimise geometry for stability margin.

**Novelty check (done at the start):** no published work found that trains a surrogate for K and M
and couples it to a stability model for design. Closest work: MIT/Pisa (Wan, Spakovszky et al. 2021)
identify K, M with a reduced-order model (no ML, no optimisation); steady NPSH surrogates and
DL-based instability classification from video exist but target different questions.
Caveats: Chinese-language literature (CNKI) and industry in-house tools not searchable.

**Phases:** 0 novelty ✅ · 1 Clark-Y validation ✅ (closed with documented limitation) ·
2 cascade K, M ❌ (killed 2026-11-17) · 3 surrogate/optimisation ⬜ · 4 3D inducer ⬜

---

## 2. Phase 1 — Clark-Y hydrofoil validation (`clarky/`)

**Case:** Clark-Y 11.7 %, c = 70 mm, α = 8°, U = 10 m/s, Re = 7×10⁵, tunnel 10c × 2.7c, LE 4c from
inlet, water at 25 °C (Huang et al. 2013 / Wang et al. 2001). 2D, 32.6k-cell structured mesh with a
body-fitted inner layer (`make_blockMesh_v2.py`). Solver: `interPhaseChangeFoam`, Schnerr–Sauer,
k-ω SST with Reboud density correction (own library `clarky/src/kOmegaSSTReboud`).

### Key results (final template v5)
| | CFD | Experiment |
|---|---|---|
| Non-cavitating Cl (σ = 2.0) | 1.177 | 1.152 (Wang Fig. 16, digitised) |
| Cloud-shedding (σ = 0.8) | 24 Hz, St 0.168 | 20–26 Hz, St 0.185 |
| Mesh independence (σ = 1.4, 73.7k vs 32.6k, CloudHPC) | Cl +1.8 %, Cd +5 %, same shedding freq. | — |

σ sweep vs digitised Wang data (`clarky/sigma_sweep_v2.png/.csv`):
| σ_in | Cl CFD | Cl exp | Cd CFD | Cd exp |
|---|---|---|---|---|
| 0.853 | 0.760 | 0.850 | 0.139 | 0.116 |
| 1.091 | 0.852 | 1.093 | 0.137 | 0.106 |
| 1.458 | 1.055 | 1.166 | 0.112 | 0.040 |
| 2.005 | 1.177 | 1.152 | 0.023 | 0.037 |

**Limitation (documented):** the model has no stable partial-sheet regime. The experiment keeps
Cl ≈ 1.17 with a thin stable sheet down to σ ≈ 1.2; CFD loses lift from σ ≈ 2 and has a large
unsteady cavity at σ = 1.4 (Cd 2.5–3× too high). Not fixed by: mesh refinement, Schnerr–Sauer
coefficients (Cv ×0.1, ×0.01; nuclei n 1e11, 1e9), Kunz model, Reboud n = 3. Merkle was unstable
(Δt → 1e-13); an effective-tension test (negative pSat) crashes Schnerr–Sauer (FPE).
Cause: equilibrium inception threshold (vapour wherever p < pSat; steady Cp_min = −2.71) plus
RANS/homogeneous-mixture sheet stability — a known model-fidelity limit.

**Check C (single-foil compliance from the sweep):** K_foil = −∂(V_c/c²)/∂σ ≈ 0.037 / 0.065 / 0.058
at σ 1.73 / 1.27 / 0.97 — correct sign, plausible magnitude, trend unresolved (shedding scatter).

Full history: `clarky/validation_log.md`. Targets: `clarky/reference_data.md`,
`clarky/wang2001_fig16_8deg_v1.csv`.

---

## 3. Phase 2 — 2D cascade K, M (`cascade/`)

**Target:** Otsuka, Tsujimoto, Kamijo & Furuya (1996) J. Fluids Eng. 118:400 — closed-cavity theory
for a flat-plate cascade, β = 15° (stagger 75°), α = 2.5–7.5°, σ/(2α) ≈ 0.5–2, with comparison to
measured K, M (SSME LPOTP, ARIANE V LOX). Definitions and plan: `cascade_notes.md`.

**Case:** one blade passage, blade-aligned coordinates, translational cyclicAMI periodicity, pitch
d = 50 mm, chord 3d. Mesh v1 (`make_cascade_blockMesh_v1.py`, 18.3k cells, zero-thickness plates);
mesh v2 (`make_cascade_blockMesh_v2.py`, 24.6k cells, t = 1 mm, semicircular LE, tapered TE).
Skewed faces only in the inlet/outlet trapezoid blocks.

### Run tally — no run held a steady partial cavity at controlled inlet σ
| Run | Setup | Outcome |
|---|---|---|
| smoke (v1) | outlet pressure set | stable, but σ_inlet = 0.014: cascade raises static p by 0.19–0.38 q |
| smoke2 (v2) | pressure inlet + flowRate outlet | diverged at 6 ms (before any vapour) |
| smoke3 (v3) | outlet-p feedback controller on inlet p | controller oscillated/diverged once a cavity formed (t = 0.093 s) |
| smoke4 (v4) | controller frozen after ramp | stable but choked supercavity (whole suction side, σ_in → 0.014) — real static instability of fixed-outlet-p + fixed-flow operation |
| smoke5 (v5) | inlet fixed U + inlet static p | diverged at cavity onset (t = 0.081 s) |
| A1 (v6) | v5 + slower ramp, Co 0.5, nOuter 3 | diverged at σ ≈ 0.2 |
| A2 (v7) | A1 + finite-thickness rounded-LE blade | diverged at σ ≈ 0.45 |
| CF1 (cf_v1, CloudHPC) | `cavitatingFoam` (compressible barotropic), mapped liquid start | FPE at cavity onset, σ = 1.41, V_vap ≈ 3e-10 m³ |

**Diagnosis:** both ESI solver families fail at cavity onset in the cascade. Likely contributors:
(1) extreme stiffness at phase change (ρ_l/ρ_v ≈ 43,000; barotropic ψ jumps ~16×), (2) my inlet
BC (fixed U + fixed p) is over-constrained — a real system has compliance (tank/lines) to absorb
cavity volume changes, (3) ESI OpenFOAM has only transient, segregated cavitation solvers (no
steady-state or coupled cavitating solver), (4) tight high-solidity passage with AMI.
**Not impossible in general** — cavitating inducers are routinely simulated with commercial codes
(steady cavitating RANS with σ stepping) and in published OpenFOAM inducer studies.

**Kill condition 2 (cascade cannot hold quasi-steady partial cavities) — criteria MET, EXECUTED 2026-11-17 by the owner.**
The fair statement is "this tool chain/setup could not, within the time box", not "infeasible".

---

## 4. Important lessons / bugs found (reusable)

1. **`limited corrected 0.5` (cavitatingBullet tutorial schemes) corrupts VoF solvers** on meshes
   with curved non-orthogonal cells: `fvc::reconstruct` turns the partially corrected flux into a
   spurious near-wall jet at the LE → Cl drifted 1.18 → 1.9 with negative drag, even in pure liquid
   (`interFoam` too; `pimpleFoam` unaffected). Fix: `corrected` laplacian/snGrad. Found via ~12
   controlled tests (`clarky/runs/diag_*`).
2. **`foamDictionary -set` on a dict with `#eval` entries freezes the evaluated values** (the
   "σ = 2.0" cloud run actually ran at σ = 0.8). Use `sed` for parameters.
3. Start the cavitating run from a **high-σ hold** (σ = 10) then ramp; low hold σ let the transient LE
   suction peak cavitate and the cavity persisted.
4. **Cascade with outlet pressure set + fixed flow is statically unstable** for partial cavities
   (growth lowers inlet p → more growth → choke). σ must be referenced/imposed upstream.
5. A per-time-step outlet-p controller becomes unstable once a cavity exists (cavity compliance
   delays the response — a surge-like loop).
6. `cavitatingFoam` impulsive start = numerical water hammer (ρ c ΔU ≈ 15 MPa); start from a mapped
   developed liquid flow.
7. Laptop: wrap long runs in `systemd-inhibit` (a suspend killed runs once); 8 cores reached
   95.5 °C → use ≤ 6 cores; `scripts/monitor_temps_v1.sh` logs temps with a 97 °C alert.
8. Never `pkill -f <pattern>` from a shell whose command line contains the pattern (it kills itself).

---

## 5. Repository map

| Path | Contents |
|---|---|
| `README.md` | project summary |
| `handover.md` | this file |
| `cascade_notes.md` | Otsuka definitions, cascade plan, full cascade setup log |
| `clarky/validation_log.md` | full Phase 1 log (bug diagnosis, sweeps, calibration, check C) |
| `clarky/reference_data.md`, `wang2001_fig16_8deg_v1.csv` | experimental targets (digitised) |
| `clarky/case_template{,_v2..v5}` | Clark-Y templates; **v5 = validated** (corrected schemes, Reboud, σ=10 hold) |
| `clarky/scripts/make_blockMesh_v1/v2.py` | Clark-Y mesh generators (v2 used) |
| `clarky/scripts/run_case_v1..v9.sh` | runners (v7: v5 template; v8: SS overrides; v9: model/pSat/Reboud overrides) |
| `clarky/scripts/postprocess_v1..v3.py`, `compare_sigma_v1/v2.py` | analysis (v3 + compare_v2 current) |
| `clarky/scripts/monitor_temps_v1.sh` | temperature logger with alert |
| `clarky/src/kOmegaSSTReboud` | Reboud-corrected k-ω SST (incompressible; `wmake libso`) |
| `cascade/case_template_v1..v7`, `case_template_cf_v1` | cascade templates (see run tally) |
| `cascade/scripts/make_cascade_blockMesh_v1/v2.py` | cascade meshes (zero-thickness / finite-thickness) |
| `cascade/scripts/run_cascade_*.sh` | cascade runners (cf_v2 supports MAP_FROM/MAP_TIME) |
| `cascade/scripts/postprocess_cascade_v1/v2.py` | cascade analysis (v2 handles alpha.vapour) |
| `papers/` (git-ignored) | Wang 2001, Huang 2013, Otsuka 1996 PDFs + text |

Run data (git-ignored, ~8 GB): `clarky/runs/` (5.1 GB), `cascade/runs/` (2.9 GB). Only
`summary.txt` / `history.png` of key runs are tracked. Cloud results: `clarky/runs/cloud_v1`
(invalid, pre-fix), `clarky/runs/cloud_v5` (fine mesh), `cascade/runs/cloud_cf1` (CF1).

**CloudHPC usage:** 4 jobs (2 × ~51 min, 1 × 3.9 h, 1 × 4.7 h on 4 cores) plus one user-run stale
σ = 2.0 job → roughly 80–90 of 300 vCPU-h if 4 cores bill as 8 vCPU; check the dashboard.

---

## 6. Open decision and next steps

**Decision pending (user):** invoke kill condition 2 for the OpenFOAM cascade route?
Recommended: close the *OpenFOAM cascade route*, keep the *project* open, via one of:

1. **Commercial solver** (Ansys CFX/Fluent — check institute academic licence): steady cavitating
   RANS with σ stepping gives quasi-steady K directly; this is standard practice for inducers.
2. **OpenFOAM BC rethink:** total-pressure inlet + mass-flow outlet (+ explicit upstream
   compliance), and debug why v2 diverged pre-cavitation. Uncertain, days.
3. **Collaborate:** email inducer-dynamics groups (MIT GTL — Spakovszky; Pisa/Alta — Pasini; Osaka/
   JAXA line) with Phase 1 results, the scheme-bug finding and the cascade choking analysis; ask to
   test the K/M-surrogate idea with their solver.
4. **If stopping:** write the assessment note — Phase 1 validation + limitation, scheme bug, cascade
   choking instability, onset stiffness in both solver families, check C.

Papers still needed: Brennen, Meissner, Lo & Hoffman (1982) J. Fluids Eng. 104(4):428 (measured
K, M); Wan, Spakovszky, Manfredi & Pasini (2021) J. Eng. Gas Turbines Power 143(2):021011
(MIT DSpace); Tsujimoto, Kamijo & Yoshida (1993) J. Fluids Eng. 115(1):135.

Hard kill conditions (agreed 2026-11-15; not self-executing — user go-ahead required):
(1) wrong K/M trend vs Otsuka; (2) no quasi-steady partial cavities; (3) no clear answer by
2026-11-29 or ~25 cascade runs.
