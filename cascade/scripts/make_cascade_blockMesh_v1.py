#!/usr/bin/env python3
"""blockMeshDict for one passage of a 2D flat-plate cascade (zero-thickness blades), blade-aligned.

Geometry (Otsuka et al. 1996 baseline): vane angle beta = 15 deg (stagger 75 deg), pitch d.
Blade i lies on y = 0, x in [0, Lc]; blade i+1 is blade i shifted by p = d (cos beta, sin beta).
Domain lower boundary = suction side of blade i (flow at +alpha, i.e. upward), upper boundary =
pressure side of blade i+1. Up/downstream of the blades the lower/upper lines are translational
cyclicAMI pairs (separation vector p). Inlet/outlet edges are parallel to p so all interior blocks
are rectangles (orthogonal cells at both leading edges).

Blocks (x-stations, dc = d cos beta, h = d sin beta):
  B1 inlet trapezoid  [Xin .. x1]      B2 [x1 .. 0]         B3 [0 .. dc]   (suction-side LE)
  B4 [dc .. Lc] (passage, blade both sides)                 B5 [Lc .. Lc+dc]
  B6 [Lc+dc .. x6]    B7 outlet trapezoid [x6 .. Xout]
Usage: make_cascade_blockMesh_v1.py <caseDir> [--d 0.05] [--beta 15] [--Lc 3] [--fine 1]
"""
import argparse
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("case")
p.add_argument("--d", type=float, default=0.05, help="pitch [m]")
p.add_argument("--beta", type=float, default=15.0, help="vane angle from cascade axis [deg]")
p.add_argument("--Lc", type=float, default=3.0, help="blade chord [d]")
p.add_argument("--up", type=float, default=2.0, help="inlet distance upstream of LE [d]")
p.add_argument("--down", type=float, default=2.5, help="outlet distance downstream of TE_(i+1) [d]")
p.add_argument("--h1", type=float, default=5e-5, help="first wall cell height [m]")
p.add_argument("--dz", type=float, default=0.001)
p.add_argument("--fine", type=float, default=1.0)
a = p.parse_args()

d, b = a.d, np.radians(a.beta)
dc, h = d * np.cos(b), d * np.sin(b)
Lc = a.Lc * d
Xin, x1 = -a.up * d, -0.9 * d
x6 = Lc + dc + 0.5 * d
Xout = x6 + a.down * d
assert x1 > Xin + dc, "inlet trapezoid not convex"

f = a.fine
N = lambda n: max(4, int(round(n * f)))
Ny = N(60)
n1, n2, n3, n4, n5, n6, n7 = N(15), N(20), N(70), N(120), N(40), N(10), N(30)

# two-sided wall-normal grading across the passage height h: first cell h1/f at both walls
def two_sided(L, h1, n):
    half = n // 2
    lo, hi = 1.0 + 1e-9, 3.0
    for _ in range(200):
        q = 0.5 * (lo + hi)
        if h1 * (q ** half - 1) / (q - 1) > L / 2:
            hi = q
        else:
            lo = q
    r = q ** (half - 1)
    return f"((0.5 0.5 {r:.4f}) (0.5 0.5 {1 / r:.4f}))"
gy = two_sided(h, a.h1 / f, Ny)

# x-gradings (last/first): cluster at LE of blade i (x=0), LE of blade i+1 (x=dc), TEs
g = {1: 0.25, 2: 0.15, 3: "((0.5 0.5 4) (0.5 0.5 0.5))", 4: "((0.3 0.4 3) (0.7 0.6 0.5))",
     5: "((0.5 0.5 2) (0.5 0.5 0.5))", 6: 1.0, 7: 6.0}

lo_x = [Xin, x1, 0.0, dc, Lc, Lc + dc, x6, Xout]           # lower-line stations (0..7)
up_x = [Xin + dc, x1, 0.0, dc, Lc, Lc + dc, x6, Xout + dc]  # upper-line stations (8..15)
V = [(x, 0.0) for x in lo_x] + [(x, h) for x in up_x]
nv = len(V)

def vtx(i, z):
    return f"    ({V[i][0]:.8f} {V[i][1]:.8f} {z})"

def hexb(k, nx, gx):
    vs = (k, k + 1, 8 + k + 1, 8 + k)
    vv = " ".join(map(str, vs)) + " " + " ".join(str(v + nv) for v in vs)
    return f"    hex ({vv}) ({nx} {Ny} 1) edgeGrading ({gx} {gx} {gx} {gx}  {gy} {gy} {gy} {gy}  1 1 1 1)"

blocks = [hexb(k, n, g[k + 1]) for k, n in enumerate([n1, n2, n3, n4, n5, n6, n7])]

def fl(k):  # lower face of block k
    return f"({k} {k + 1} {k + 1 + nv} {k + nv})"
def fu(k):  # upper face of block k
    return f"({8 + k} {8 + k + nv} {8 + k + 1 + nv} {8 + k + 1})"

sep = f"({dc:.8f} {h:.8f} 0)"
sepN = f"({-dc:.8f} {-h:.8f} 0)"
patches = [
    ("inlet", "patch", [f"(0 {8} {8 + nv} {nv})"], ""),
    ("outlet", "patch", [f"(7 {7 + nv} {15 + nv} 15)"], ""),
    ("bladeSS", "wall", [fl(2), fl(3)], ""),            # suction side of blade i
    ("bladePS", "wall", [fu(3), fu(4)], ""),            # pressure side of blade i+1
    ("cycUpLow", "cyclicAMI", [fl(0), fl(1)],
     f"neighbourPatch cycUpUp; transform translational; separationVector {sep};"),
    ("cycUpUp", "cyclicAMI", [fu(0), fu(1), fu(2)],
     f"neighbourPatch cycUpLow; transform translational; separationVector {sepN};"),
    ("cycDnLow", "cyclicAMI", [fl(4), fl(5), fl(6)],
     f"neighbourPatch cycDnUp; transform translational; separationVector {sep};"),
    ("cycDnUp", "cyclicAMI", [fu(5), fu(6)],
     f"neighbourPatch cycDnLow; transform translational; separationVector {sepN};"),
]

ncell = Ny * (n1 + n2 + n3 + n4 + n5 + n6 + n7)
out = ["FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }",
       f"// make_cascade_blockMesh_v1.py: d={d} beta={a.beta} Lc={a.Lc}d h1={a.h1 / f:.2e} fine={f}",
       f"// dc={dc:.6f} h={h:.6f} cells={ncell}", "scale 1;", "vertices", "("]
out += [vtx(i, 0) for i in range(nv)] + [vtx(i, a.dz) for i in range(nv)]
out += [");", "edges ();", "blocks", "("] + blocks + [");", "boundary", "("]
for name, typ, faces, extra in patches:
    out += [f"    {name}", "    {", f"        type {typ};"]
    if extra:
        out += ["        " + e.strip() + ";" for e in extra.rstrip(";").split(";")]
    out += ["        faces", "        ("] + [f"            {fc}" for fc in faces] + ["        );", "    }"]
out += [");", "defaultPatch { name frontAndBack; type empty; }"]
open(f"{a.case}/system/blockMeshDict", "w").write("\n".join(out) + "\n")
print(f"cascade blockMeshDict: {ncell} cells, dc={dc:.4f} h={h:.4f} Lc={Lc:.3f} m, sep={sep}")
