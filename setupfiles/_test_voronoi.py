"""Verification for the VoronoiGrid target (both methods)."""
import numpy as np
import time
from sph_interp import Particles, VoronoiGrid, interpolate

rng = np.random.default_rng(0)


def mk_particles(pos, mass=None, rho=None, h=None, vals=None, box=None,
                 periodic=False):
    N = pos.shape[0]
    if mass is None: mass = np.ones(N)
    if rho is None: rho = np.ones(N)
    if h is None: h = np.full(N, 0.5)
    if vals is None: vals = {"rho": rho.copy(), "ones": np.ones(N)}
    return Particles(pos=pos, mass=mass, rho=rho, h=h, values=vals,
                     box=box, periodic=periodic)


print("=== 1. single particle, kernel fully inside a random tessellation ===")
# generators filling [-2,2]^3; particle at origin with h=0.5 (support r<1)
gens = rng.uniform(-2.0, 2.0, size=(400, 3))
grid = VoronoiGrid.from_points(gens, bounds=[[-2, 2]] * 3)
print(f"  Ncell={grid.ncell}, vol_sum={grid.volume.sum():.4f} (box=64)")
p = mk_particles(np.zeros((1, 3)), h=np.array([0.5]))

t0 = time.time()
res = interpolate(p, grid, method="petkova")
m = res.data["_mass"].sum()
print(f"  petkova: total mass = {m:.6f} (expect 1.0), err {100*(m-1):+.3f}%  "
      f"[{time.time()-t0:.2f}s]")

res_s = interpolate(p, grid, method="sph")
ms = (res_s.data["rho"] * grid.volume).sum()
print(f"  sph:     sum(rho*vol) = {ms:.6f} (approx)")

print("\n=== 2. many particles, mass conservation over the tessellation ===")
# uniform random particles inside [-1,1]^3, mesh covers [-2,2]^3 so all kernels
# (h=0.3, support 0.6) stay inside the meshed domain
Np = 300
ppos = rng.uniform(-1.0, 1.0, size=(Np, 3))
pmass = rng.uniform(0.5, 1.5, size=Np)
prho = rng.uniform(0.8, 1.2, size=Np)
pp = mk_particles(ppos, mass=pmass, rho=prho, h=np.full(Np, 0.3),
                  vals={"rho": prho, "ones": np.ones(Np)})
res2 = interpolate(pp, grid, method="petkova")
m2 = res2.data["_mass"].sum()
print(f"  petkova: deposited mass {m2:.5f} vs particle mass {pmass.sum():.5f}"
      f"  err {100*(m2/pmass.sum()-1):+.3f}%")

print("\n=== 3. uniform density recovery (SPH) ===")
# particles on a jittered lattice, constant rho; SPH value of 'ones' ~ 1
lat = np.stack(np.meshgrid(*[np.linspace(-0.9, 0.9, 8)] * 3, indexing="ij"), -1)
lpos = lat.reshape(-1, 3) + rng.uniform(-0.02, 0.02, size=(512, 3))
lp = mk_particles(lpos, h=np.full(512, 0.35),
                  vals={"ones": np.ones(512), "rho": np.ones(512)})
res3 = interpolate(lp, grid, method="sph")
ones = res3.data["ones"]
hit = res3.norm > 0
print(f"  sph 'ones' over hit cells: mean {ones[hit].mean():.4f} "
      f"min {ones[hit].min():.4f} max {ones[hit].max():.4f} (expect ~1)")

print("\n=== 4. from_particles round trip ===")
vg = VoronoiGrid.from_particles(pp)
print(f"  from_particles -> Ncell={vg.ncell} (expect {Np}), "
      f"vol_sum={vg.volume.sum():.4f}")
res4 = interpolate(pp, vg, method="petkova")
m4 = res4.data["_mass"].sum()
print(f"  petkova on particle-tessellation: mass {m4:.5f} vs {pmass.sum():.5f} "
      f"err {100*(m4/pmass.sum()-1):+.3f}%")

print("\nDONE")
