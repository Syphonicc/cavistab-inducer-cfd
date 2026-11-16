#!/usr/bin/env python3
"""blockMeshDict for one passage of a 2D cascade of FINITE-THICKNESS blades (v2).

v2 vs v1: blade thickness t = 2r with a semicircular leading edge (radius r) and a trailing edge
tapered to a point over Lt; v1 used zero-thickness plates (sharp LE -> singular suction peak, stiff
cavity onset). Blade-aligned coordinates: blade i has its mid-plane on y = 0, LE nose point at the
origin, suction half (y >= 0) is the domain's lower wall; blade i+1 = blade i + p, p = d(cos b, sin b),
its pressure half (y <= h) is the upper wall. Up/downstream of the blades the mid-plane lines are
translational cyclicAMI pairs.

At each nose the boundary turns 90 deg; the block edge leaving the nose follows the 45 deg bisector
(two 45 deg corners instead of a degenerate one), and short transition blocks confine the skew.
Lower stations: L0 L1 N P1 P2 Q R1 T0 TE D1 D2 L6 L7 ; upper: U0 U1 Uq Us1 Us2 N' R1u T0u TEu D1u D2u U6 U7
Usage: make_cascade_blockMesh_v2.py <caseDir> [--r 0.0005] [--Lt 0.01] [--fine 1]
"""
import argparse
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("case")
p.add_argument("--d", type=float, default=0.05)
p.add_argument("--beta", type=float, default=15.0)
p.add_argument("--Lc", type=float, default=3.0, help="chord [d]")
p.add_argument("--r", type=float, default=0.0005, help="LE radius = half thickness [m]")
p.add_argument("--Lt", type=float, default=0.01, help="TE taper length [m]")
p.add_argument("--up", type=float, default=2.0)
p.add_argument("--down", type=float, default=2.5)
p.add_argument("--h1", type=float, default=5e-5)
p.add_argument("--dz", type=float, default=0.001)
p.add_argument("--fine", type=float, default=1.0)
a = p.parse_args()

d, b, r, Lt = a.d, np.radians(a.beta), a.r, a.Lt
dc, h = d * np.cos(b), d * np.sin(b)
Lc = a.Lc * d
Xin, x1 = -a.up * d, -0.9 * d
s1, s2 = 0.010, 0.026                       # end of LE transition / start of N' transition
q = dc - (h - r)                            # foot of the 45 deg edge from N' on blade i
s3 = dc + 0.010
x6 = Lc + dc + 0.5 * d
Xout = x6 + a.down * d
assert x1 < -h and s2 < q and s3 < Lc - Lt and Lc < dc + Lc - Lt

lo = [(Xin, 0), (x1, 0), (0, 0), (s1, r), (s2, r), (q, r), (s3, r), (Lc - Lt, r), (Lc, 0),
      (Lc + dc - Lt, 0), (Lc + dc, 0), (x6, 0), (Xout, 0)]
up = [(Xin + dc, h), (x1, h), (-h, h), (s1, h), (s2, h), (dc, h), (s3, h - r), (Lc - Lt, h - r),
      (Lc, h - r), (Lc + dc - Lt, h - r), (Lc + dc, h), (x6, h), (Xout + dc, h)]
nb = len(lo) - 1
V = lo + up
nv = len(V)
U0 = len(lo)                                 # index offset of upper points

f = a.fine
N = lambda n: max(4, int(round(n * f)))
Ny = N(60)
#            B1   B2   B3a  B3b  B3c  B4a  B4b   B4c  B5   B5b  B6   B7
ncx = [N(n) for n in (15, 20, 50, 25, 25, 50, 110, 20, 40, 15, 10, 30)]
gx = [0.25, 0.1, 10.0, 1.0, 0.5, 6.0, "((0.3 0.4 2) (0.7 0.6 0.5))", 0.5, 2.0, 0.5, 1.0, 6.0]

def two_sided(L, h1, n):
    half = n // 2
    lo_, hi_ = 1.0 + 1e-9, 3.0
    for _ in range(200):
        qq = 0.5 * (lo_ + hi_)
        if h1 * (qq ** half - 1) / (qq - 1) > L / 2:
            hi_ = qq
        else:
            lo_ = qq
    rr = qq ** (half - 1)
    return f"((0.5 0.5 {rr:.4f}) (0.5 0.5 {1 / rr:.4f}))"
gy = two_sided(h - 2 * r, a.h1 / f, Ny)

def vtx(i, z):
    return f"    ({V[i][0]:.8f} {V[i][1]:.8f} {z})"

def pl(pts, z):
    return " ".join(f"({x:.8f} {y:.8f} {z})" for x, y in pts)

# nose arcs: blade i suction half (centre (r,0), 180 -> 90 deg), blade i+1 pressure half
th1 = np.radians(np.linspace(180, 90, 24))[1:]
arcLo = [(r + r * np.cos(t), r * np.sin(t)) for t in th1]               # ends at (r, r)
th2 = np.radians(np.linspace(180, 270, 24))[1:]
arcUp = [(dc + r + r * np.cos(t), h + r * np.sin(t)) for t in th2]      # ends at (dc+r, h-r)

edges = []
for k, z in ((0, 0), (nv, a.dz)):
    edges.append(f"    polyLine {2 + k} {3 + k} ({pl(arcLo, z)})")             # N -> P1
    edges.append(f"    polyLine {U0 + 5 + k} {U0 + 6 + k} ({pl(arcUp, z)})")   # N' -> R1u

blocks = []
for k in range(nb):
    vs = (k, k + 1, U0 + k + 1, U0 + k)
    vv = " ".join(map(str, vs)) + " " + " ".join(str(v + nv) for v in vs)
    g = gx[k]
    blocks.append(f"    hex ({vv}) ({ncx[k]} {Ny} 1) edgeGrading ({g} {g} {g} {g}  "
                  f"{gy} {gy} {gy} {gy}  1 1 1 1)  // B{k}")

def fl(k):
    return f"({k} {k + 1} {k + 1 + nv} {k + nv})"
def fu(k):
    return f"({U0 + k} {U0 + k + nv} {U0 + k + 1 + nv} {U0 + k + 1})"

sepP = f"({dc:.8f} {h:.8f} 0)"
sepN = f"({-dc:.8f} {-h:.8f} 0)"
patches = [
    ("inlet", "patch", [f"(0 {U0} {U0 + nv} {nv})"], []),
    ("outlet", "patch", [f"({nb} {nb + nv} {U0 + nb + nv} {U0 + nb})"], []),
    ("bladeSS", "wall", [fl(k) for k in range(2, 8)], []),            # N .. TE (suction half, blade i)
    ("bladePS", "wall", [fu(k) for k in range(5, 10)], []),           # N' .. TE' (pressure half, i+1)
    ("cycUpLow", "cyclicAMI", [fl(0), fl(1)],
     ["neighbourPatch cycUpUp", "transform translational", f"separationVector {sepP}"]),
    ("cycUpUp", "cyclicAMI", [fu(k) for k in range(0, 5)],
     ["neighbourPatch cycUpLow", "transform translational", f"separationVector {sepN}"]),
    ("cycDnLow", "cyclicAMI", [fl(k) for k in range(8, 12)],
     ["neighbourPatch cycDnUp", "transform translational", f"separationVector {sepP}"]),
    ("cycDnUp", "cyclicAMI", [fu(k) for k in range(10, 12)],
     ["neighbourPatch cycDnLow", "transform translational", f"separationVector {sepN}"]),
]

ncell = Ny * sum(ncx)
out = ["FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }",
       f"// make_cascade_blockMesh_v2.py: d={d} beta={a.beta} Lc={a.Lc}d r={r} (t={2 * r}) Lt={Lt} fine={f}",
       f"// dc={dc:.6f} h={h:.6f} cells={ncell}", "scale 1;", "vertices", "("]
out += [vtx(i, 0) for i in range(nv)] + [vtx(i, a.dz) for i in range(nv)]
out += [");", "edges", "("] + edges + [");", "blocks", "("] + blocks + [");", "boundary", "("]
for name, typ, faces, extra in patches:
    out += [f"    {name}", "    {", f"        type {typ};"] + [f"        {e};" for e in extra]
    out += ["        faces", "        ("] + [f"            {fc}" for fc in faces] + ["        );", "    }"]
out += [");", "defaultPatch { name frontAndBack; type empty; }"]
open(f"{a.case}/system/blockMeshDict", "w").write("\n".join(out) + "\n")
print(f"cascade v2 blockMeshDict: {ncell} cells, t={2 * r * 1e3:.2f} mm, dc={dc:.4f} h={h:.4f} Lc={Lc:.3f} m")
