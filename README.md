# cavistab

Instability-aware design of rocket turbopump inducers.

Inducers are usually optimised for suction performance (NPSH), while the failures that matter in
testing are cavitation instabilities (rotating cavitation, cavitation surge). Their onset is governed
by two cavity properties: the cavitation compliance **K** and the mass-flow gain factor **M**.
This project builds a CFD-trained surrogate geometry → (K, M), couples it to a 1-D stability model,
and optimises inducer geometry for stability margin.

## Status: stopped 2026-11-17 (OpenFOAM cascade route closed)
- Phase 1 (Clark-Y hydrofoil validation): closed. Good at σ = 2.0 and shedding frequency; no stable
  partial-sheet regime for any model tested. See [`clarky/validation_log.md`](clarky/validation_log.md).
- Phase 2 (2D cascade K, M): **killed.** 8 runs with `interPhaseChangeFoam` and `cavitatingFoam` never
  held a steady partial cavity. See [`cascade_notes.md`](cascade_notes.md).
- Phases 3–4 (surrogate, 3D inducer): never started.

This is a negative-results repository. Read in this order:
1. [`docs/FINAL_ASSESSMENT.md`](docs/FINAL_ASSESSMENT.md) — decision, what was and wasn't achieved, how to continue
2. [`docs/PROJECT_HISTORY.md`](docs/PROJECT_HISTORY.md) — every attempt, every mistake
3. [`handover.md`](handover.md) — full state snapshot and repo map

Reusable: the `limited corrected` scheme bug (VoF solvers on non-orthogonal meshes), the
Reboud-corrected k-ω SST library, the cascade mesh generators.

## Layout
| Path | Contents |
|---|---|
| `clarky/case_template_v*` | OpenFOAM case templates (current: v5) |
| `clarky/scripts/` | mesh generator, run scripts, post-processing (versioned, never overwritten) |
| `clarky/src/kOmegaSSTReboud` | k-ω SST with Reboud density correction (`wmake libso`) |
| `clarky/geometry/clarky.dat` | Clark-Y coordinates (UIUC Airfoil Coordinates Database) |
| `clarky/reference_data.md` | experimental targets |
| `clarky/runs/v5_*/` | tracked: `summary.txt`, `history.png` only |
| `cascade/` | cascade case templates v1–v7 + cavitatingFoam, mesh generators, runners, post-processing |
| `docs/` | final assessment and full project history |

## Toolchain
OpenFOAM v2512 (ESI), `interPhaseChangeFoam`, Schnerr–Sauer, Python 3 (numpy, matplotlib).

## Not included
Run data (~8 GB) and the source papers (copyrighted PDFs) are not in the repository. Experimental
targets are digitised from published figures; see the original papers for authoritative values.
