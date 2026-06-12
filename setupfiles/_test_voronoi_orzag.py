import numpy as np, time
from sph_interp import from_tipsy, VoronoiGrid, interpolate

snap = ('../test_cases/orzag/orzag64_N64lattice_3deta05GASOLINE/'
        'orzag64_N64lattice_3deta05GASOLINE.00100')
p = from_tipsy(snap, fields=("rho", "Bmag"), periodic=False)
print(f"N={p.n}, box={p.box}, particle mass total={p.mass.sum():.5f}", flush=True)

half = 0.5 * p.box
t0 = time.time()
vg = VoronoiGrid.from_particles(
    p, bounds=[[-half[a] - 0.05, half[a] + 0.05] for a in range(3)])
print(f"Voronoi build: Ncell={vg.ncell}, vol_sum={vg.volume.sum():.4f} "
      f"(box vol={np.prod(p.box):.4f}) [{time.time()-t0:.1f}s]", flush=True)

t0 = time.time()
res = interpolate(p, vg, method="petkova")
m = res.data["_mass"].sum()
print(f"petkova: deposited mass {m:.5f} vs {p.mass.sum():.5f}  "
      f"err {100*(m/p.mass.sum()-1):+.3f}% [{time.time()-t0:.1f}s]", flush=True)
print(f"  rho range [{res.data['rho'].min():.4f},{res.data['rho'].max():.4f}], "
      f"Bmag finite: {np.isfinite(res.data['Bmag']).all()}", flush=True)

t0 = time.time()
res_s = interpolate(p, vg, method="sph")
print(f"sph: rho range [{res_s.data['rho'].min():.4f},"
      f"{res_s.data['rho'].max():.4f}] [{time.time()-t0:.1f}s]", flush=True)
print("DONE", flush=True)
