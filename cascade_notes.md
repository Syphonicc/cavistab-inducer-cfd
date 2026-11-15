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
