# Phase 2 notes — 2D cascade K, M (validation target)

## Otsuka, Tsujimoto, Kamijo & Furuya (1996), J. Fluids Eng. 118(2) — `papers/Otsuka1996_JFE.pdf`
- Linear closed-cavity model of partial cavitation on a **flat-plate cascade** (2D), unsteady.
- Cavitation number σ = (p₁ − p_v)/(½ρU²), U = upstream velocity, α = flow incidence.
- Cavity area V_c (per blade, 2D). Definitions (non-dimensionalised by spacing d):
  - **K = −∂(V_c/d²)/∂σ** (cavitation compliance)
  - **M = −∂(V_c/d²)/∂(α)** — mass-flow gain via incidence (flow-rate) perturbation
  - Complex in general; **quasi-steady (k = 0) values are real** — what steady/averaged CFD gives.
- Baseline geometry: vane angle β = 15° (stagger γ = 75°), spacing d = 2π (cavity length units);
  Fig. 9 comparison uses β = 10.27° ("typical modern inducer").
- Results vs cavity length l and σ/(2α); effects of incidence α and mean cavity length (Figs 3–8).
- **Fig. 9: comparison with measured K, M** (tendencies only — paper warns not quantitative):
  SSME LPOTP (Brennen et al. 1982), ARIANE V LOX inducer (Iacopozzi et al. 1993),
  LE-7 LOX (Shimura 1993). Reduced frequency k* = ωd/U.

## Plan for the CFD cascade (to design)
1. Flat-plate / thin-blade cascade, stagger 75° (β = 15°) — match Otsuka baseline first.
2. Quasi-steady K: run σ ± Δσ at fixed α, K = −ΔV_c/Δσ (normalised by d²).
   Quasi-steady M: run α ± Δα at fixed σ, M = −ΔV_c/Δα.
3. Compare K, M vs σ/(2α) with Otsuka's quasi-steady curves; then β = 10.27° vs Fig. 9 data.
4. Clark-Y finding to carry over: model lacks a stable-sheet regime → expect CFD cavities longer/
   more unsteady at given σ; compare at matched **cavity length** as well as matched σ.

## References still needed (exact citations from Otsuka's reference list)
- Tsujimoto, Kamijo & Yoshida (1993), "A Theoretical Analysis of Rotating Cavitation in
  Inducers," J. Fluids Eng. 115(1):135–141.
- Brennen, Meissner, Lo & Hoffman (1982), "Scale Effects in the Dynamic Transfer Functions for
  Cavitating Inducers," J. Fluids Eng. 104(4):428–433. ← measured K, M (SSME LPOTP)
- Wan, Spakovszky, Manfredi & Pasini (2021), J. Eng. Gas Turbines Power 143(2):021011
  — open access, MIT DSpace handle 1721.1/134548.

## Cascade setup log
- v1 (`case_template_v1`, outlet-imposed pressure): smoke at σ_out 0.2 gave σ_inlet 0.014 — the cascade
  raises static pressure by 0.38 q (non-cavitating) / 0.19 q (heavily cavitating), so inlet σ was
  uncontrolled; cavity still growing at 0.12 s (near-supercavitating, σ/2α ≈ 0.08). Not a valid test.
- v2 (`case_template_v2`, `run_cascade_v2.sh`): Otsuka's conditions — inlet static p_rgh imposed
  (ramped σ 10 → target), inflow direction fixed at α (pressureDirectedInletVelocity), outlet fixed
  volumetric flow U·d·sin(β−α)·dz (flowRateOutletVelocity). BC check: inlet p exact, stable.
- Post-processing: `scripts/postprocess_cascade_v1.py` (σ_in, V_c/d², A_cav/d², L_cav/d, std/mean).
- v2 diverged at 6.15 ms (Cl → 1e33, Δt → 1e-32): pressure inlet + flowRate outlet unstable here.
- v3 (`case_template_v3`, `run_cascade_v3.sh`): v1's stable BCs (fixed U at inlet, p at outlet) +
  codedFixedValue outlet controller driving the inlet area-averaged p_rgh to the hold/ramp target
  (gain 0.05, once per time step, |Δp| ≤ 5 kPa/step). Check: inlet p within 0.5 % of target in 2 ms.
- v3 smoke (resumed from 0.05 s) diverged at t = 0.093 s, ~13 ms after the ramp: with a cavity present
  the per-step controller drove inlet p into 0.5–3.9 MPa oscillations (Δt → 7e-9). Cavity compliance
  delays the outlet→inlet pressure response, so the closed loop over-corrects (a surge-like loop).
- v4 (`case_template_v4`, `run_cascade_v4.sh`): controller only during hold + ramp (gain 0.02/step),
  outlet p frozen for t ≥ tRamp → stable v1-type BCs while the cavity develops; inlet σ then drifts,
  so each run records its actual σ_in (K is a slope through the runs' (σ_in, V_c) points).
- v4 smoke (controller frozen after ramp) ran stably to 0.20 s but drifted to σ_in = 0.014: one cavity
  covering the whole suction side (x 0–0.153 m, up to 36 % of passage height), passage at p ≈ pv —
  choked supercavity. Fixed outlet p + fixed flow is statically unstable for partial cavities
  (growth lowers inlet p → more growth). Experiments set σ upstream (tank) → must impose inlet p.
- v5 (`case_template_v5`, `run_cascade_v5.sh`): inlet fixedValue U (U, α) + inlet static p_rgh
  ramp; outlet zeroGradient U (inletOutlet) and p. 4.9 ms check: stable, inlet p exact.
- v5 smoke diverged at t = 0.0809 s, exactly when the ramp reached σ_in = 0.2 and vapour appeared
  (Δt → 1.7e-10 s). Hold phase was clean (Cl ≈ 0.1, inlet p on target).
- Tally (2026-11-16): 5 cascade smokes — v1 wrong σ (outlet-referenced), v2 diverged (6 ms),
  v3 controller oscillated/diverged, v4 choked supercavity (σ_in 0.014), v5 diverged at cavity onset.
  No run has yet held a partial cavity at controlled inlet σ. Common factor in the stiff failures:
  cavity onset at the zero-thickness sharp LE with an inlet pressure constraint.
- Bounded retry A (agreed 2026-11-16, max 2 runs / 3 days, then kill condition 2 if no steady partial cavity):
  A1 = `case_template_v6` (v5 BCs, ramp to 0.15 s, maxCo/maxAlphaCo 0.5, nOuter 3), zero-thickness mesh;
  A2 = finite-thickness blade with rounded LE (new mesh topology, in progress).
- A2 = `case_template_v7` + `make_cascade_blockMesh_v2.py` (t = 1 mm, semicircular LE r = 0.5 mm,
  TE taper 10 mm; 24.6k cells; nose blocks use 45° bisector edges; bad faces only in end blocks),
  `run_cascade_v7.sh`, 2 cores alongside A1 (6 cores total, temp limit). Both started 2026-11-16.
