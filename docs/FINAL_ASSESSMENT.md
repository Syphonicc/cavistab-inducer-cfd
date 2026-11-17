# Final assessment — where the project stopped (2026-11-17)

## Decision
**Kill condition 2 invoked by the project owner on 2026-11-17: the OpenFOAM cascade route is closed.**
The condition (agreed 2026-11-15): *the cascade cannot hold quasi-steady partial cavities, so K and M
are undefined.* It was met by 8 runs across both ESI solver families (`interPhaseChangeFoam` ×7
configurations, `cavitatingFoam` ×1). See [`PROJECT_HISTORY.md`](PROJECT_HISTORY.md).

The accurate statement is: **this tool chain and setup could not do it within the time box.** It is
not a proof that the problem is infeasible. Cavitating inducers are routinely simulated with
commercial codes (steady cavitating RANS with σ stepping) and in published OpenFOAM studies.

## What was delivered
- A Clark-Y cavitating-hydrofoil validation (Phase 1), good at σ = 2.0 and for shedding frequency,
  not in the stable-sheet regime (σ ≈ 1.2–1.6), for every model tested.
- A reusable bug finding: `limited corrected 0.5` schemes corrupt VoF solvers on curved
  non-orthogonal meshes.
- A Reboud-corrected k-ω SST library for OpenFOAM v2512 (`clarky/src/kOmegaSSTReboud`).
- A cascade mesh generator (zero-thickness and finite-thickness blades), 8 case templates and
  runners, and an analysis script for cascade V_c, σ_in and cavity length.
- Analysis of why the fixed-outlet-pressure and per-step-controller cascade setups are unstable.

## What was not achieved
- No steady partial cavity in the cascade, so no CFD K or M, no comparison with Otsuka (1996).
- No stable-sheet regime on the hydrofoil, so a trustworthy single-foil K could not be extracted.
- Phase 3 (surrogate and optimisation) and Phase 4 (3D inducer) were never started.

## Open hypotheses (untested)
1. The over-constrained inlet (fixed U + fixed p) is the main cause. Try total-pressure inlet +
   mass-flow outlet with explicit upstream compliance.
2. Onset stiffness is intrinsic to the segregated solvers. A coupled or steady cavitating solver
   would avoid it.
3. The inception threshold (p < pSat with no tension or nuclei limit) is the Phase 1 error. Needs a
   modified Schnerr–Sauer that accepts an offset.

## Possible continuations (none started)
1. Commercial solver (Ansys CFX/Fluent, institute licence): steady cavitating RANS with σ stepping.
2. OpenFOAM boundary-condition rethink (hypothesis 1), days of work, uncertain.
3. Collaboration: groups with working cavitating-inducer tools (MIT GTL, Pisa/Alta, Osaka/JAXA).
4. Theory-only route: Otsuka-type linear closed-cavity model for K, M, validated against published
   measurements (Brennen et al. 1982), skipping CFD.

Papers still needed: Brennen et al. (1982) J. Fluids Eng. 104(4):428; Wan et al. (2021) JEGTP
143(2):021011; Tsujimoto et al. (1993) J. Fluids Eng. 115(1):135.

## Reproducibility notes
- Run data (~8 GB) and the source PDFs are not in this repository. Only `summary.txt` / `history.png`
  of key runs are tracked.
- Experimental targets are digitised from Wang et al. (2001) Fig. 16 (±0.01 Cl, ±0.002 Cd) in
  `clarky/wang2001_fig16_8deg_v1.csv`; see the original paper for authoritative values.
- Hardware: 4-core laptop (≤ 6 cores for thermal reasons), CloudHPC 4-core jobs.
