"""Tests de LIBRE sobre trayectorias sintéticas exactas (sintetico.trayectoria, con blur)."""
import numpy as np
import pytest

from pinzas import fisica, libre
from pinzas.sintetico import trayectoria

T = 296.15
R = 1.5e-6                    # PS de 3 µm
ETA = 1.7e-3                  # 20 % glicerina
GAMMA = fisica.gamma_stokes(ETA, R)
D = fisica.difusion(T, GAMMA)  # ~0.085 µm²/s
FPS, W = 35.0, 10e-3
SIG = 0.05 * 90e-9            # 0.05 px a 90 nm/px


def video(semilla, n=2000, nb=5, perdidos=20, jitter=2e-5, frac_nan=0.01, D_=D):
    """Un 'video' de nb partículas libres: tiempos con jitter y cuadros perdidos, ruido, NaN."""
    rng = np.random.default_rng(semilla)
    t = np.arange(n) / FPS + rng.normal(0, jitter, n)
    keep = np.ones(n, bool)
    keep[rng.choice(np.arange(1, n - 1), perdidos, replace=False)] = False
    t = t[keep]
    X = np.stack([trayectoria(t, W, D_, rng=rng)[0] for _ in range(nb)], 1)
    Y = np.stack([trayectoria(t, W, D_, rng=rng)[0] for _ in range(nb)], 1)
    X += rng.normal(0, SIG, X.shape)
    Y += rng.normal(0, SIG, Y.shape)
    X[rng.random(X.shape) < frac_nan] = np.nan
    return t, X, Y


def test_msd_usa_retardos_reales_con_huecos():
    # movimiento balístico x = v t: MSD(tau) = v² tau² exacto si se usan los retardos reales
    dt = 0.01
    t = np.delete(np.arange(200) * dt, [50, 51, 120])        # tres cuadros perdidos
    v = 1e-6
    r = libre.msd(t, v * t, 0 * t, max_lag=10)
    m = np.arange(1, 11)
    assert np.allclose(r['tau'], m * dt)
    assert np.allclose(r['msd'][:, 0], (v * m * dt) ** 2, rtol=1e-9)
    # con índices en lugar de tiempos, el retardo 1 mezclaría pares de 1 y de 3 cuadros
    ingenuo = np.mean(np.diff(v * t) ** 2)
    assert ingenuo > 1.01 * (v * dt) ** 2


def test_ajuste_msd_recupera_modelo_exacto():
    tau = np.arange(1, 11) / FPS
    s2 = 3e-17
    m = libre.modelo_msd(tau, D, s2, W)
    f = libre.ajustar_msd(tau, m, W, n_puntos=5)
    assert f['D'] == pytest.approx(D, rel=1e-9)
    assert f['sigma2'] == pytest.approx(s2, rel=1e-6)


def test_D_dentro_del_2_por_ciento_en_10_videos():
    rel = [libre.analizar(*video(s), W, T, R)['D_medio'] / D - 1 for s in range(10)]
    print(f'\nLIBRE: D/D_verdadero - 1 = {np.mean(rel):+.4f} ± {np.std(rel) / np.sqrt(10):.4f} (10 videos)')
    assert abs(np.mean(rel)) < 0.02


def test_eta_faxen_e_incertidumbre():
    t, X, Y = video(3, nb=8)
    h = 3 * R
    r = libre.analizar(t, X, Y, W, T, R, h=h)
    assert r['eta'] == pytest.approx(ETA, rel=4 * r['err_eta'] / r['eta'])
    assert r['eta_faxen'] == pytest.approx(r['eta'] / fisica.faxen_paralelo(R, h))
    # el error por dispersión entre partículas es del orden del real (D_i con ~5 % cada uno)
    assert 0.005 < r['err_D'] / r['D_medio'] < 0.05
    assert r['msd'].shape == r['msd_ajuste'].shape == (30, 8)


def test_histograma_de_desplazamientos_gaussiano():
    t, X, Y = video(4, nb=10)
    h = libre.analizar(t, X, Y, W, T, R)['hist']
    assert h['p_ks'] > 1e-3
    assert abs(h['curtosis_exceso']) < 4 * h['err_curtosis']
    assert h['cuentas'].sum() == pytest.approx(h['gauss'].sum(), rel=0.01)


def test_correccion_por_resta_de_deriva_comun():
    # DERIVA resta el desplazamiento medio de las nb partículas: cada MSD baja en (1 - 1/nb)
    nb = 4
    sin, con = [], []
    for s in range(10):
        t, X, Y = video(100 + s, nb=nb, frac_nan=0)
        X = X - X.mean(1, keepdims=True)
        Y = Y - Y.mean(1, keepdims=True)
        sin.append(libre.analizar(t, X, Y, W, T, R)['D_medio'] / D)
        con.append(libre.analizar(t, X, Y, W, T, R, n_deriva=nb)['D_medio'] / D)
    assert np.mean(sin) == pytest.approx(1 - 1 / nb, abs=0.03)
    assert np.mean(con) == pytest.approx(1, abs=0.03)


def test_una_sola_particula_error_nan():
    t, X, Y = video(5, nb=1)
    r = libre.analizar(t, X, Y, W, T, R)
    assert np.isfinite(r['D_medio']) and np.isnan(r['err_D'])
    assert r['avisos']
