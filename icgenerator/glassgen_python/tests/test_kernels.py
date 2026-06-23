import numpy as np

from glassgen.kernels import dw_wc2, w_wc2, W0


def test_unit_volume_integral():
    # int W dV = (1/pi h^3) * 4 pi int w(q) q^2 h^3 dq = 4 int w(q) q^2 dq = 1
    # (the Wendland C2 normalization 21/16 is baked into w_wc2 for exactly this)
    q = np.linspace(0.0, 2.0, 200001)
    w = np.array([w_wc2(x * x) for x in q])
    integral = 4.0 * np.trapz(w * q * q, q)
    assert abs(integral - 1.0) < 1e-6


def test_compact_support():
    assert w_wc2(4.0) == 0.0
    assert w_wc2(9.0) == 0.0
    assert dw_wc2(4.0) == 0.0
    # Wendland C2 central value is 21/16 (used by the self contribution)
    assert w_wc2(0.0) == W0
    assert abs(W0 - 21.0 / 16.0) < 1e-12
    assert dw_wc2(0.0) == 0.0


def test_derivative_matches_finite_difference():
    # dw_wc2(q^2) = (dW/dq)/q
    eps = 1e-6
    for q in [0.3, 0.7, 1.0, 1.3, 1.5, 1.9]:
        dwdq = (w_wc2((q + eps) ** 2) - w_wc2((q - eps) ** 2)) / (2 * eps)
        assert abs(dw_wc2(q * q) - dwdq / q) < 1e-5, q


def test_smooth_single_polynomial():
    # Wendland C2 is a single polynomial on [0, 2] (no piecewise break at
    # q=1): continuous value AND derivative everywhere interior
    for qb in (0.5, 1.0, 1.5):
        assert abs(w_wc2(qb ** 2 - 1e-9) - w_wc2(qb ** 2 + 1e-9)) < 1e-7
        assert abs(dw_wc2(qb ** 2 - 1e-9) - dw_wc2(qb ** 2 + 1e-9)) < 1e-7
