import numpy as np
import pytest

from glassgen import resample


def _setup(n=400, seed=0):
    rng = np.random.default_rng(seed)
    box = np.array([2.0, 2.0, 2.0])
    pos = (rng.random((n, 3)) - 0.5) * box
    mass = np.full(n, 3.0)
    aux = {'vel': rng.random((n, 3)), 'u': np.full(n, 2.0),
           'B': rng.random((n, 3)), 'momenti': np.full(n, 1.0),
           'rho': np.full(n, 1.0), 'h': np.full(n, 0.1)}
    return pos, mass, aux, box


def test_split_count_and_mass_conservation():
    pos, mass, aux, box = _setup()
    n = len(pos)
    for factor in (2, 4, 8):
        p2, m2, a2 = resample.split(pos, mass, aux, box=box, factor=factor)
        assert len(p2) == n * factor
        assert np.isclose(m2.sum(), mass.sum())          # total mass conserved
        assert np.allclose(m2, mass[0] / factor)         # per-child = parent/factor


def test_split_rejects_non_power_of_two():
    pos, mass, aux, box = _setup()
    with pytest.raises(ValueError):
        resample.split(pos, mass, aux, box=box, factor=3)


def test_split_field_rules():
    pos, mass, aux, box = _setup()
    p2, m2, a2 = resample.split(pos, mass, aux, box=box, factor=8)
    # extensive fields (u, B, momenti) are halved 3x -> sum conserved
    assert np.isclose(aux['u'].sum(), a2['u'].sum())
    assert np.allclose(aux['B'].sum(0), a2['B'].sum(0))
    assert np.isclose(aux['momenti'].sum(), a2['momenti'].sum())
    # intensive fields inherited -> momentum conserved (vel mass-weighted)
    p0 = (mass[:, None] * aux['vel']).sum(0)
    p1 = (m2[:, None] * a2['vel']).sum(0)
    assert np.allclose(p0, p1)


def test_split_children_local_and_in_box():
    pos, mass, aux, box = _setup()
    p2, m2, a2 = resample.split(pos, mass, aux, box=box, factor=2,
                                nsmooth=64, dist=0.4)
    assert np.isfinite(p2).all()
    assert np.all(np.abs(p2) <= box / 2 + 1e-9)
    # each child pair sits near its parent: min-image separation small (< box/2)
    sep = p2[0::2] - p2[1::2]
    sep -= box * np.rint(sep / box)
    assert np.all(np.linalg.norm(sep, axis=1) < box.min() / 2)


def test_split_custom_halve_fields():
    pos, mass, aux, box = _setup()
    aux2 = {'Bx': aux['B'][:, 0].copy()}
    # by default 'Bx' is inherited (not in HALVE_FIELDS); override to halve it
    _, _, a_inh = resample.split(pos, mass, aux2, box=box, factor=2)
    _, _, a_hlv = resample.split(pos, mass, aux2, box=box, factor=2,
                                 halve_fields={'Bx'})
    assert np.isclose(a_inh['Bx'].sum(), 2 * aux2['Bx'].sum())  # inherited
    assert np.isclose(a_hlv['Bx'].sum(), aux2['Bx'].sum())      # halved


def test_merge_conserves_mass_and_momentum():
    pos, mass, aux, box = _setup(n=600)
    pm, mm, am = resample.merge(pos, mass, aux, box=box, target_factor=2)
    assert len(pm) < len(pos)                       # some pairs merged
    assert np.isclose(mm.sum(), mass.sum())         # mass conserved
    p0 = (mass[:, None] * aux['vel']).sum(0)
    p1 = (mm[:, None] * am['vel']).sum(0)
    assert np.allclose(p0, p1)                       # momentum conserved
    # extensive fields summed back
    assert np.isclose(am['u'].sum(), aux['u'].sum())


def test_open_domain_split():
    # box=None (open) should still split and conserve mass
    rng = np.random.default_rng(1)
    pos = rng.random((300, 3))
    mass = np.full(300, 1.0)
    p2, m2, _ = resample.split(pos, mass, {}, box=None, factor=2)
    assert len(p2) == 600
    assert np.isclose(m2.sum(), mass.sum())
    assert np.isfinite(p2).all()


def test_coarsen_exact_count_and_roundtrip():
    # coarsen() hits count//factor EXACTLY and is the inverse-by-count of split:
    # split(coarsen(x, F), F) restores the original count (1:1).
    pos, mass, aux, box = _setup(n=512)        # divisible by 8
    n = len(pos)
    for factor in (2, 4, 8):
        cp, cm, ca = resample.coarsen(pos, mass, aux, box=box, factor=factor)
        assert len(cp) == n // factor                  # EXACT target count
        assert np.isclose(cm.sum(), mass.sum())        # mass conserved
        assert np.isclose(ca['u'].sum(), aux['u'].sum())   # extensive aux summed
        sp, sm, _ = resample.split(cp, cm, ca, box=box, factor=factor)
        assert len(sp) == n                            # round-trips to input N
        assert np.isclose(sm.sum(), mass.sum())


def test_coarsen_rejects_non_power_of_two():
    pos, mass, aux, box = _setup()
    with pytest.raises(ValueError):
        resample.coarsen(pos, mass, aux, box=box, factor=3)
