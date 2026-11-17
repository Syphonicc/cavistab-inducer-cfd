# Project history — every attempt, mistake and outcome

Chronological record, including what failed and why. Detailed per-phase logs:
[`clarky/validation_log.md`](../clarky/validation_log.md), [`cascade_notes.md`](../cascade_notes.md).
Toolchain: OpenFOAM v2512 (ESI), `interPhaseChangeFoam` / `cavitatingFoam`, Python 3, a 4-core laptop,
plus CloudHPC (300 free vCPU-h) for a few long jobs.

## Phase 0 — novelty check (done)
No published work found that trains a surrogate for cavitation compliance K and mass-flow gain
factor M and couples it to a 1-D stability model for design. Closest: Wan, Spakovszky, Manfredi &
Pasini (2021), which identifies K, M with a reduced-order model (no ML, no optimisation).
Not searchable: Chinese-language literature, industry in-house tools.

## Phase 1 — Clark-Y hydrofoil validation (closed 2026-11-15, with a documented limitation)
Case: Clark-Y 11.7 %, c = 70 mm, α = 8°, U = 10 m/s, 2D, 32.6k cells; targets from Wang et al.
(2001) and Huang et al. (2013).

| Date | What was tried | Outcome |
|---|---|---|
| early | `interPhaseChangeFoam`, Schnerr–Sauer, k-ω SST, schemes copied from the `cavitatingBullet` tutorial | Cl drifted to 1.9–2.5, Cd negative, spurious cavitation at σ = 2–3 |
| 11-12 | ~12 controlled isolation tests (simpleFoam, pimpleFoam, interFoam, BC/tolerance/PISO variants) | **Root cause: `limited corrected 0.5` laplacian/snGrad corrupts VoF solvers** on curved non-orthogonal cells (`fvc::reconstruct` → near-wall jet at the LE). `corrected` fixes it (Cl 1.183, Cd 0.0136 in interFoam). |
| 11-12 | Fixed template (v5), Reboud-corrected SST, σ = 10 hold then ramp | σ = 2.0: Cl 1.177 vs 1.152 (exp). σ = 0.8: shedding 24 Hz, St 0.168 vs 0.185 |
| 11-12 | σ sweep 2.0 / 1.4 / 1.0 / 0.8 | Cl, Cd off between σ ≈ 1.1 and 1.5; CFD loses lift gradually from σ ≈ 2, experiment keeps Cl ≈ 1.17 down to σ ≈ 1.2 |
| 11-12 | Mesh independence (73.7k vs 32.6k, CloudHPC) | Cl +1.8 %, Cd +5 %, same shedding frequency → not a resolution effect |
| 11-12 | Schnerr–Sauer coefficients (Cv ×0.1, ×0.01; nuclei 1e11, 1e9) | Cd 0.09–0.12 vs exp 0.04: mass-transfer rate is not the cause |
| 11-13 | Digitised Wang Fig. 16; checked σ reference (upstream pressure) | "Agreement at σ = 0.8" was an artefact of eyeballed reference values (CFD Cl actually 0.09 low). Not a constant σ offset. |
| 11-15 | Kunz, Reboud n = 3, Merkle, effective tension (negative pSat) | Kunz/Reboud-3 no better; **Merkle** collapsed (Δt → 1e-13 s); **tension test** crashed (Schnerr–Sauer floors at 0.01·pSat, invalid for pSat < 0) |
| 11-16 | Check "C": single-foil quasi-steady compliance from the sweep | K_foil ≈ 0.037 / 0.065 / 0.058: correct sign, plausible size, trend unresolved (shedding scatter ≈ the differences) |

Go/no-go (Cd < 0.06 and Cl > 1.12 at σ = 1.4): no model passes. The model has no stable
partial-sheet regime. Likely cause: equilibrium inception threshold (vapour wherever p < pSat; steady
Cp_min = −2.71) plus RANS/homogeneous-mixture sheet stability.

## Phase 2 — 2D cascade quasi-steady K and M (killed 2026-11-17)
Target: Otsuka, Tsujimoto, Kamijo & Furuya (1996), flat-plate cascade, β = 15°, pitch d = 50 mm.
Needs steady partial cavities at several controlled inlet σ and incidences, so that
K = −∂(V_c/d²)/∂σ and M = −∂(V_c/d²)/∂α are defined.

| Run | Setup | Outcome |
|---|---|---|
| smoke (v1) | outlet pressure set | stable, but σ_inlet = 0.014: the cascade raises static p by 0.19–0.38 q, so inlet σ was uncontrolled |
| smoke2 (v2) | pressure inlet + flowRate outlet (Otsuka's conditions) | diverged at 6 ms, before any vapour |
| smoke3 (v3) | outlet-pressure controller steering inlet p | controller oscillated and diverged once a cavity formed (t = 0.093 s): cavity compliance delays the response, a surge-like loop |
| smoke4 (v4) | controller frozen after ramp | stable, but a choked supercavity covering the whole suction side (σ_in → 0.014): fixed outlet p + fixed flow is statically unstable for partial cavities |
| smoke5 (v5) | fixed inlet U + inlet static p ramp | diverged at cavity onset (t = 0.081 s) |
| A1 (v6) | v5 + slower ramp, Co 0.5, nOuter 3 | diverged at σ ≈ 0.2 |
| A2 (v7) | A1 + finite-thickness rounded-LE blade (new mesh) | diverged at σ ≈ 0.45, mid-ramp, which weakens the "ramp-end kink" hypothesis |
| CF1 (cf_v1) | `cavitatingFoam` (compressible barotropic), mapped liquid start, CloudHPC, 4.7 h | FPE at cavity onset, σ_in = 1.41, V_vap ≈ 3e-10 m³ |

Notes on the CF1 setup: an impulsive start crashed within 0.09 ms (numerical water hammer,
ρcΔU ≈ 15 MPa), so it was started from A2's mapped developed liquid flow. Two earlier local CF1
attempts were stopped by the user, one at t = 0.042 s with no checkpoint.

Diagnosis (hypotheses, not proven): (1) extreme stiffness at phase change (ρ_l/ρ_v ≈ 43 000;
barotropic ψ jumps ≈ 16× when vapour appears); (2) my inlet BC (fixed U + fixed p) over-constrains
the problem, a real system has upstream compliance; (3) ESI OpenFOAM only has transient, segregated
cavitation solvers; (4) tight high-solidity passage with cyclicAMI.

## Mistakes and lost effort (honest list)
1. Ran the first `interPhaseChangeFoam` cases, including a CloudHPC job, with the tutorial's
   `limited corrected 0.5` schemes. All of them are invalid (cloud_v1, reboud_s0.8_v3, local_v2, v4).
2. `foamDictionary -set` on a dictionary with `#eval` entries froze the evaluated values: the run
   labelled σ = 2.0 actually ran at σ = 0.8. Use `sed`.
3. Compared against eyeballed experimental values before digitising Wang Fig. 16, which produced a
   false "agreement at σ = 0.8".
4. A low initial σ hold let the start-up suction peak cavitate and the cavity persisted. Start from σ = 10.
5. Cascade v1 imposed outlet pressure, so inlet σ was uncontrolled and the "σ = 0.2" run was really 0.014.
6. A system suspend killed a long run (13:40, 11-12); an 8-core run hit 95.5 °C. Use `systemd-inhibit`, ≤ 6 cores.
7. `pkill -f <pattern>` from a shell whose own command line contained the pattern killed that shell.
8. Starting `cavitatingFoam` impulsively gives a numerical water hammer. Map a developed liquid state.
9. Cascade v1–v5 were tried one at a time before the stiffness at onset was understood. A cheaper
   diagnosis (a single-cell phase-change stiffness test) was not done first.
