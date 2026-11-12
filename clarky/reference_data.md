# Clark-Y validation reference data

Sources: Wang et al. (2001) PAS 37:551–581 (`papers/Wang2001_PAS.pdf`),
Huang et al. (2013) JFE 135:071301 (`papers/Huang2013_JFE.pdf`).

## Setup (Huang 2013, Sec. 2 and 4; Fig. 2)
- Clark-Y 11.7 %, c = 0.07 m, α = 8°, U∞ = 10 m/s, Re = 7×10⁵
- Test section 10c long × 2.7c high × 1c wide; LE 4c from inlet; foil ~1.4c above bottom wall
  (model: centred, i.e. 0.05c offset ignored)
- Water 25 °C: ρl 999.19, μl 1.139e-3, ρv 0.02308, μv 9.8626e-6, pv 3169 Pa
- Inlet Tu 2 %; symmetry top/bottom; outlet pressure set from σ
- Init fully wetted, tunnel pressure lowered progressively to target σ
- Their CFD: CFX, 150k cells, y⁺ ≈ 1, k–ε with DCM/FBDCM eddy-viscosity corrections

## Targets
| Quantity | σ = 2.0 (subcavitating) | σ = 0.8 (sheet/cloud) | Source |
|---|---|---|---|
| Cl | ≈ 1.16 | ≈ 0.75–0.8 | Wang Fig. 16, eyeballed — digitize properly |
| Cd | ≈ 0.037 | ≈ 0.11–0.12 | Wang Fig. 16, eyeballed — digitize properly |
| Shedding frequency | — | 20 Hz (Wang Table 2); St = fc/U = 0.185 → ≈ 26 Hz, T_ref ≈ 40 ms = 5.7 c/U (Huang text) | |
| Mean cavity area (10 % vapour isoline) | — | 0.72 Sc measured, 0.70 Sc predicted | Huang Sec. 5.1 |

Sc (Clark-Y section area) = 0.0809 c² (from UIUC coordinates).

## Caveats
- Huang Table 1 lists the measured mean shedding frequency as 69.2 Hz, inconsistent with the
  same paper's St = 0.185 / T_ref = 40 ms and with Wang's 20 Hz. Treat ≈ 20–26 Hz as the target
  unless resolved.
- Huang: standard k–ε gave 25.5 Hz, DCM 46.5 Hz, FBDCM 72.0 Hz (vs their 69.2 Hz) — i.e.
  unmodified eddy-viscosity models change the shedding dynamics substantially.
- Wang regime map at 8°: sheet cavity σ ≈ 1.40–1.55, cloud σ ≈ 0.8, supercavitation σ ≈ 0.28.
