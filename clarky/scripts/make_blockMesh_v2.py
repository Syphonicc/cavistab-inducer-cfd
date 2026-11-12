#!/usr/bin/env python3
"""Generate system/blockMeshDict for a 2D Clark-Y hydrofoil in a water tunnel (v2).

v2 vs v1: a body-fitted inner layer (foil offset along its normals by --delta) wraps the
foil, so near-wall grid lines are orthogonal and the LE C-block no longer folds.

Blocks: B0u/B0l upstream | I1u/I1l LE (split at LE), I2 upper, I3 lower (inner layer) | W1/W2 inner wake |
        O1u/O1l front, O2 top, O3 bottom, O4/O5 outer wake
Normal direction runs foil -> outer everywhere, so shared edges share one grading spec.

Usage: make_blockMesh_v2.py <caseDir> [--dat ...] [--aoa 8] [--fine 1.0]
"""
import argparse
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("case")
p.add_argument("--dat", default="geometry/clarky.dat")
p.add_argument("--chord", type=float, default=0.07)
p.add_argument("--aoa", type=float, default=8.0, help="deg, nose up, about mid-chord")
p.add_argument("--H", type=float, default=0.095, help="half tunnel height [m]")
p.add_argument("--up", type=float, default=3.0, help="inlet distance upstream of LE [c]")
p.add_argument("--down", type=float, default=6.0, help="outlet distance downstream of TE [c]")
p.add_argument("--dz", type=float, default=0.001, help="span thickness [m]")
p.add_argument("--split", type=float, default=0.4, help="x/c of upper/lower block split")
p.add_argument("--delta", type=float, default=0.08, help="inner layer thickness [c]")
p.add_argument("--h1", type=float, default=4e-5, help="first wall cell height [m]")
p.add_argument("--fine", type=float, default=1.0, help="global cell-count multiplier")
a = p.parse_args()
c = a.chord

# ---- read Lednicer-format coordinates: upper LE->TE, lower LE->TE
rows = []
for line in open(a.dat).read().splitlines()[1:]:
    s = line.split()
    if len(s) == 2:
        rows.append((float(s[0]), float(s[1])))
nu, nl = int(rows[0][0]), int(rows[0][1])
up = np.array(rows[1:1 + nu])
lo = np.array(rows[1 + nu:1 + nu + nl])
up[-1, 1] = lo[-1, 1] = 0.0          # close blunt TE (0.12%c) to a sharp point

# ---- scale, rotate about mid-chord (positive aoa = nose up)
th = np.radians(a.aoa)
R = np.array([[np.cos(th), np.sin(th)], [-np.sin(th), np.cos(th)]])
def xf(pts):
    return (((pts - [0.5, 0.0]) @ R.T) + [0.5, 0.0]) * c
upX, loX = xf(up), xf(lo)

# ---- outward normals on the open curve TE(lower) -> LE -> TE(upper)
loop = np.vstack([loX[::-1], upX[1:]])           # clockwise? check sign below
t = np.gradient(loop, axis=0)
t /= np.linalg.norm(t, axis=1)[:, None]
n = np.column_stack([t[:, 1], -t[:, 0]])
if n[nl - 1, 0] > 0:                             # LE normal must point upstream (-x)
    n = -n
nLo, nUp = n[:nl][::-1], np.vstack([n[nl - 1], n[nl:]])
d = a.delta * c
upO, loO = upX + d * nUp, loX + d * nLo

iu = int(np.argmin(abs(up[:, 0] - a.split)))
il = int(np.argmin(abs(lo[:, 0] - a.split)))
LE, TE, U1, L1 = upX[0], upX[-1], upX[iu], loX[il]
U1o, L1o, TEu, TEl = upO[iu], loO[il], upO[-1], loO[-1]

xLE = LE[0]
Xi, Xf, Xo = xLE - a.up * c, xLE - 1.0 * c, TE[0] + a.down * c
H = a.H
V = [(Xi, -H), (Xf, -H), (Xf, H), (Xi, H),                       # 0-3
     tuple(L1), tuple(U1), tuple(TE),                             # 4-6
     tuple(L1o), tuple(U1o), tuple(TEu), tuple(TEl),              # 7-10
     (TEu[0], H), (TEl[0], -H),                                   # 11-12
     (Xo, H), (Xo, TEu[1]), (Xo, TE[1]), (Xo, TEl[1]), (Xo, -H),  # 13-17
     tuple(LE), tuple(upO[0]), (Xf, upO[0][1]), (Xi, upO[0][1])]  # 18-21
nv = len(V)

# ---- cell counts
f = a.fine
Nup, Nc, Nin, Nout, Nu, Nl, Nw = (max(4, int(round(k * f)))
                                  for k in (20, 110, 45, 30, 70, 45, 90))

def arclen(pts):
    return float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))

def geo(L, h1, N):
    """total expansion (last/first) for N cells over L with first cell h1"""
    lo_, hi_ = 1.0 + 1e-9, 3.0
    for _ in range(200):
        q = 0.5 * (lo_ + hi_)
        if h1 * (q ** N - 1) / (q - 1) > L:
            hi_ = q
        else:
            lo_ = q
    return q ** (N - 1), h1 * q ** (N - 1)

Gin, hEdge = geo(d, a.h1 / f, Nin)
Lout = np.mean([np.linalg.norm(np.array(V[2]) - U1o), H - TEu[1], H + TEl[1]])
Gout, _ = geo(Lout, hEdge * 1.2, Nout)

sU, sL = arclen(upX[:iu + 1]), arclen(loX[:il + 1])
fU = sU / (sU + sL)
cU = int(round(Nc * fU))
gLE = 6.0
gTEu, gTEl, gW = 0.5, 0.6, 25.0

# ---- dictionary text
def vtx(i, z):
    x, y = V[i]
    return f"    ({x:.8f} {y:.8f} {z})  // {i + (nv if z else 0)}"

def pl(pts, z):
    return " ".join(f"({x:.8f} {y:.8f} {z})" for x, y in pts)

def interior_wrap(U, L, i_u, i_l, le):
    return np.vstack([U[i_u - 1:0:-1], [le], L[1:i_l]])

edges = []
for k, z in ((0, 0), (nv, a.dz)):
    edges += [
        f"    polyLine {5 + k} {18 + k} ({pl(upX[iu - 1:0:-1], z)})",
        f"    polyLine {18 + k} {4 + k} ({pl(loX[1:il], z)})",
        f"    polyLine {5 + k} {6 + k} ({pl(upX[iu + 1:-1], z)})",
        f"    polyLine {4 + k} {6 + k} ({pl(loX[il + 1:-1], z)})",
        f"    polyLine {8 + k} {19 + k} ({pl(upO[iu - 1:0:-1], z)})",
        f"    polyLine {19 + k} {7 + k} ({pl(loO[1:il], z)})",
        f"    polyLine {8 + k} {9 + k} ({pl(upO[iu + 1:-1], z)})",
        f"    polyLine {7 + k} {10 + k} ({pl(loO[il + 1:-1], z)})",
    ]

def eg(g1, g2):
    # x1 edges: v0v1 v3v2 v7v6 v4v5 ; x2 edges: v0v3 v1v2 v5v6 v4v7 ; z edges
    return f"{g1[0]} {g1[1]} {g1[1]} {g1[0]}  {g2[0]} {g2[1]} {g2[1]} {g2[0]}  1 1 1 1"

def hexb(vs, nn, grading, tag):
    vv = " ".join(str(v) for v in vs) + " " + " ".join(str(v + nv) for v in vs)
    return f"    hex ({vv}) ({nn[0]} {nn[1]} 1) edgeGrading ({grading})  // {tag}"

cL = Nc - cU
blocks = [
    hexb((21, 20, 2, 3), (Nup, cU), eg((0.3, 0.3), (2, 2)), "B0u upstream"),
    hexb((0, 1, 20, 21), (Nup, cL), eg((0.3, 0.3), (0.5, 0.5)), "B0l upstream"),
    # inner layer: normal = foil -> offset
    hexb((5, 8, 19, 18), (Nin, cU), eg((Gin, Gin), (1 / gLE, 1 / gLE)), "I1u LE upper"),
    hexb((18, 19, 7, 4), (Nin, cL), eg((Gin, Gin), (gLE, gLE)), "I1l LE lower"),
    hexb((6, 9, 8, 5), (Nin, Nu), eg((Gin, Gin), (1 / gTEu, 1 / gTEu)), "I2 upper"),
    hexb((4, 7, 10, 6), (Nin, Nl), eg((Gin, Gin), (gTEl, gTEl)), "I3 lower"),
    hexb((6, 15, 14, 9), (Nw, Nin), eg((gW, gW), (Gin, Gin)), "W1 inner wake up"),
    hexb((6, 10, 16, 15), (Nin, Nw), eg((Gin, Gin), (gW, gW)), "W2 inner wake lo"),
    # outer: normal = offset -> tunnel
    hexb((8, 2, 20, 19), (Nout, cU), eg((Gout, Gout), (1 / gLE, 0.5)), "O1u front upper"),
    hexb((19, 20, 1, 7), (Nout, cL), eg((Gout, Gout), (gLE, 2)), "O1l front lower"),
    hexb((8, 9, 11, 2), (Nu, Nout), eg((gTEu, 1), (Gout, Gout)), "O2 top"),
    hexb((7, 1, 12, 10), (Nout, Nl), eg((Gout, Gout), (gTEl, 1)), "O3 bottom"),
    hexb((9, 14, 13, 11), (Nw, Nout), eg((gW, gW), (Gout, Gout)), "O4 wake top"),
    hexb((10, 12, 17, 16), (Nout, Nw), eg((Gout, Gout), (gW, gW)), "O5 wake bottom"),
]

def face(a_, b_):
    return f"({a_} {b_} {b_ + nv} {a_ + nv})"

patches = {
    "inlet": ("patch", [face(0, 21), face(21, 3)]),
    "outlet": ("patch", [face(14, 13), face(15, 14), face(16, 15), face(17, 16)]),
    "topWall": ("symmetryPlane", [face(3, 2), face(2, 11), face(11, 13)]),
    "bottomWall": ("symmetryPlane", [face(0, 1), face(1, 12), face(12, 17)]),
    "hydrofoil": ("wall", [face(5, 18), face(18, 4), face(5, 6), face(4, 6)]),
}

ncell = Nup * Nc + Nin * (Nc + Nu + Nl + 2 * Nw) + Nout * (Nc + Nu + Nl + 2 * Nw)
out = ["FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }",
       f"// make_blockMesh_v2.py: aoa={a.aoa} c={c} H={H} delta={a.delta}c h1={a.h1 / f} fine={f}",
       f"// Gin={Gin:.2f} Gout={Gout:.2f} cells={ncell}",
       "scale 1;", "vertices", "("]
out += [vtx(i, 0) for i in range(nv)] + [vtx(i, a.dz) for i in range(nv)]
out += [");", "edges", "("] + edges + [");", "blocks", "("] + blocks + [");", "boundary", "("]
for name, (typ, faces) in patches.items():
    out += [f"    {name}", "    {", f"        type {typ};", "        faces", "        ("]
    out += [f"            {fc}" for fc in faces] + ["        );", "    }"]
out += [");", "defaultPatch { name frontAndBack; type empty; }"]
open(f"{a.case}/system/blockMeshDict", "w").write("\n".join(out) + "\n")
print(f"blockMeshDict v2: {ncell} cells, Gin={Gin:.1f} Gout={Gout:.1f}, "
      f"LE=({LE[0]:.4f},{LE[1]:.4f}) TE=({TE[0]:.4f},{TE[1]:.4f})")
