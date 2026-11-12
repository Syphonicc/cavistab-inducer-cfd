# cavistab

Instability-aware design of rocket turbopump inducers.

Inducers are usually optimised for suction performance (NPSH), while the failures that matter in
testing are cavitation instabilities (rotating cavitation, cavitation surge). Their onset is governed
by two cavity properties: the cavitation compliance **K** and the mass-flow gain factor **M**.
This project builds a CFD-trained surrogate geometry → (K, M), couples it to a 1-D stability model,
and optimises inducer geometry for stability margin.

## Status
Step 1 — validate cavitating CFD on the Clark-Y hydrofoil (Wang et al. 2001; Huang et al. 2013).
See [`clarky/validation_log.md`](clarky/validation_log.md).

## Layout
| Path | Contents |
|---|---|
| `clarky/case_template_v*` | OpenFOAM case templates (current: v5) |
| `clarky/scripts/` | mesh generator, run scripts, post-processing (versioned, never overwritten) |
| `clarky/src/kOmegaSSTReboud` | k-ω SST with Reboud density correction (`wmake libso`) |
| `clarky/geometry/clarky.dat` | Clark-Y coordinates (UIUC Airfoil Coordinates Database) |
| `clarky/reference_data.md` | experimental targets |
| `clarky/runs/v5_*/` | tracked: `summary.txt`, `history.png` only |

## Toolchain
OpenFOAM v2512 (ESI), `interPhaseChangeFoam`, Schnerr–Sauer, Python 3 (numpy, matplotlib).
