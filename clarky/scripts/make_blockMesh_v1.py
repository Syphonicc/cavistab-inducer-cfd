#!/usr/bin/env python3
"""Generate system/blockMeshDict for a 2D Clark-Y hydrofoil in a water tunnel.

Topology (6 blocks, C-wrap around LE, H elsewhere):
  B0 upstream | B1 front (wraps LE) | B2 upper foil | B3 lower foil | B4/B5 wake
Normal direction always runs foil -> outer, so one grading spec fits all shared edges.

Usage: make_blockMesh_v1.py <caseDir> [--aoa 8] [--fine 1.0]
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
    q = (pts - [0.5, 0.0]) @ R.T
    return (q + [0.5, 0.0]) * c

upX, loX = xf(up), xf(lo)
iu = int(np.argmin(abs(up[:, 0] - a.split)))
il = int(np.argmin(abs(lo[:, 0] - a.split)))
LE, TE, U1, L1 = upX[0], upX[-1], upX[iu], loX[il]

def unit(v):
    return v / np.linalg.norm(v)

def normal(pts, i, sign):
    t = unit(pts[i + 1] - pts[i - 1])
    return sign * np.array([-t[1], t[0]])

nU, nL = normal(upX, iu, 1.0), normal(loX, il, -1.0)

xLE = LE[0]
Xi, Xf, Xo = xLE - a.up * c, xLE - 1.0 * c, TE[0] + a.down * c
H = a.H
v2 = {0: (Xi, -H), 1: (Xf, -H), 2: (Xf, H), 3: (Xi, H), 4: tuple(L1), 5: tuple(U1),
      6: tuple(TE), 7: (TE[0], H), 8: (TE[0], -H), 9: (Xo, H), 10: (Xo, TE[1]), 11: (Xo, -H)}
nv = len(v2)

# ---- cell counts
f = a.fine
Nup, Nc, Nn, Nu, Nl, Nw = (max(4, int(round(n * f))) for n in (20, 110, 80, 70, 45, 90))

def arclen(pts):
    return float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))

# normal-direction grading: geometric from h1 on a representative edge length
Lref = np.mean([np.linalg.norm(np.array(v2[5]) - v2[2]), H - TE[1], H + TE[1]])
def expansion(L, h1, N):
    lo_, hi_ = 1.0 + 1e-9, 2.0
    for _ in range(200):
        q = 0.5 * (lo_ + hi_)
        if h1 * (q ** N - 1) / (q - 1) > L:
            hi_ = q
        else:
            lo_ = q
    return q ** (N - 1)
Gn = expansion(Lref, a.h1, Nn)

# LE wrap (U1 -> LE -> L1): cluster at LE
sU, sL = arclen(upX[:iu + 1]), arclen(loX[:il + 1])
fU = sU / (sU + sL)
cU = int(round(Nc * fU))
gLE = 6.0
wrap = f"(({fU:.4f} {cU / Nc:.4f} {1 / gLE}) ({1 - fU:.4f} {(Nc - cU) / Nc:.4f} {gLE}))"
gTEu, gTEl = 0.5, 0.6     # mild clustering toward TE along foil
gW = 25.0                 # wake growth from TE

# ---- write
def vtx(i, z):
    x, y = v2[i]
    return f"    ({x:.8f} {y:.8f} {z})  // {i + (nv if z else 0)}"

def pl(pts, z):
    return " ".join(f"({x:.8f} {y:.8f} {z})" for x, y in pts)

def spline_pts(p0, n, p1):
    # leave the foil along the surface normal, then bend toward the outer vertex
    d = np.linalg.norm(np.array(p1) - p0)
    return [p0 + n * 0.25 * d, 0.5 * (p0 + n * 0.5 * d) + 0.5 * np.array(p1)]

wrapPts = np.vstack([upX[iu - 1:0:-1], [LE], loX[1:il]])   # interior pts U1->LE->L1
upPts = upX[iu + 1:-1]
loPts = loX[il + 1:-1]
sU1 = spline_pts(U1, nU, v2[2])
sL1 = spline_pts(L1, nL, v2[1])

edges = []
for k, z in ((0, 0), (nv, a.dz)):
    edges += [f"    polyLine {5 + k} {4 + k} ({pl(wrapPts, z)})",
              f"    polyLine {5 + k} {6 + k} ({pl(upPts, z)})",
              f"    polyLine {4 + k} {6 + k} ({pl(loPts, z)})",
              f"    spline {5 + k} {2 + k} ({pl(sU1, z)})",
              f"    spline {4 + k} {1 + k} ({pl(sL1, z)})"]

def hexb(vs, n, grading):
    vv = " ".join(str(v) for v in vs) + " " + " ".join(str(v + nv) for v in vs)
    return f"    hex ({vv}) ({n[0]} {n[1]} 1) edgeGrading ({grading})"

def eg(g1, g2):
    # edgeGrading: 4 edges x1 (v0v1, v3v2, v7v6, v4v5), 4 edges x2 (v0v3, v1v2, v5v6, v4v7), 4 edges z
    # args: g1 = (v0v1, v3v2), g2 = (v0v3, v1v2)
    return f"{g1[0]} {g1[1]} {g1[1]} {g1[0]}  {g2[0]} {g2[1]} {g2[1]} {g2[0]}  1 1 1 1"

blocks = [
    # B0 upstream: x1 inlet->Xf, x2 bottom->top (uniform, matches B1 Xf edge)
    hexb((0, 1, 2, 3), (Nup, Nc), eg((0.3, 0.3), (1, 1))),
    # B1 front: v0=U1 v1=F_T v2=F_B v3=L1; x1 normal (foil->outer), x2 U1->LE->L1 / F_T->F_B
    hexb((5, 2, 1, 4), (Nn, Nc), eg((Gn, Gn), (wrap, 1))),
    # B2 upper: v0=U1 v1=TE v2=T2 v3=F_T; x1 along upper foil / top wall, x2 normal
    hexb((5, 6, 7, 2), (Nu, Nn), eg((gTEu, 1), (Gn, Gn))),
    # B3 lower: v0=L1 v1=F_B v2=B2 v3=TE; x1 normal, x2 along lower foil / bottom wall
    hexb((4, 1, 8, 6), (Nn, Nl), eg((Gn, Gn), (gTEl, 1))),
    # B4 wake upper: v0=TE v1=O_W v2=O_T v3=T2; x1 wake, x2 normal
    hexb((6, 10, 9, 7), (Nw, Nn), eg((gW, gW), (Gn, Gn))),
    # B5 wake lower: v0=TE v1=B2 v2=O_B v3=O_W; x1 normal, x2 wake
    hexb((6, 8, 11, 10), (Nn, Nw), eg((Gn, Gn), (gW, gW))),
]

def face(a_, b_):
    return f"({a_} {b_} {b_ + nv} {a_ + nv})"

patches = {
    "inlet": ("patch", [face(0, 3)]),
    "outlet": ("patch", [face(10, 9), face(11, 10)]),
    "topWall": ("symmetryPlane", [face(3, 2), face(2, 7), face(7, 9)]),
    "bottomWall": ("symmetryPlane", [face(0, 1), face(1, 8), face(8, 11)]),
    "hydrofoil": ("wall", [face(5, 4), face(5, 6), face(4, 6)]),
    "frontAndBack": ("empty", [f"({' '.join(str(v) for v in vs)})" for vs in
                               [(0, 3, 2, 1), (1, 2, 5, 4), (5, 2, 7, 6), (4, 6, 8, 1),
                                (6, 7, 9, 10), (6, 10, 11, 8)]] +
                     [f"({' '.join(str(v + nv) for v in vs)})" for vs in
                      [(0, 1, 2, 3), (1, 4, 5, 2), (5, 6, 7, 2), (4, 1, 8, 6),
                       (6, 10, 9, 7), (6, 8, 11, 10)]]),
}

out = ["FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }",
       f"// generated by make_blockMesh_v1.py: aoa={a.aoa} c={c} H={H} h1={a.h1} fine={f}",
       f"// Gn={Gn:.2f}  cells ~ {Nup*Nc + Nn*Nc + Nu*Nn + Nn*Nl + 2*Nw*Nn}",
       "scale 1;", "vertices", "("]
out += [vtx(i, 0) for i in range(nv)] + [vtx(i, a.dz) for i in range(nv)]
out += [");", "edges", "("] + edges + [");", "blocks", "("] + blocks + [");", "boundary", "("]
for name, (typ, faces) in patches.items():
    out += [f"    {name}", "    {", f"        type {typ};", "        faces", "        ("]
    out += [f"            {fc}" for fc in faces] + ["        );", "    }"]
out += [");"]
open(f"{a.case}/system/blockMeshDict", "w").write("\n".join(out) + "\n")
print(f"blockMeshDict written: Gn={Gn:.1f}, ~{Nup*Nc + Nn*Nc + Nu*Nn + Nn*Nl + 2*Nw*Nn} cells, "
      f"LE=({LE[0]:.4f},{LE[1]:.4f}) TE=({TE[0]:.4f},{TE[1]:.4f})")
