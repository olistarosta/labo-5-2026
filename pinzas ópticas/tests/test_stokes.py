"""Tests de STOKES con datos sintéticos: OU alrededor de x_eq = gamma v / kappa.

Con la platina a velocidad v, gamma (dx/dt - v) = -kappa (x - x_trampa) + ruido. En régimen
la solución es x_eq + (el mismo OU que sin flujo), x_eq = x_trampa + gamma v / kappa: se genera
así, sumándole x_eq a la trayectoria atrapada de sintetico.trayectoria. Las beads pegadas se
mueven con la platina: x = x0 + v t + ruido.
"""
import numpy as np
import pytest

from pinzas import fisica, stokes
from pinzas.sintetico import trayectoria

T = 296.15
GAMMA = fisica.gamma_stokes(1.7e-3, 1.5e-6) * float(fisica.faxen_paralelo(1.5e-6, 5e-6))
D = fisica.difusion(T, GAMMA)
KX, KY = 3e-6, 5e-6
FPS, W, SIG = 100.0, 5e-3, 4.5e-9
X_TRAMPA, Y_TRAMPA = 12e-6, 8e-6


def segmento(v, rng, dur=20.0, v_escape=None, eje='x'):
    t = np.arange(int(dur * FPS)) / FPS
    vx, vy = (v, 0.0) if eje == 'x' else (0.0, v)
    x = X_TRAMPA + GAMMA * vx / KX + trayectoria(t, W, D, KX, GAMMA, T, 20, rng)[0]
    y = Y_TRAMPA + GAMMA * vy / KY + trayectoria(t, W, D, KY, GAMMA, T, 20, rng)[0]
    if v_escape is not None and abs(v) > v_escape:      # se escapa: el fluido se la lleva
        x = x + vx * t
        y = y + vy * t
    x += rng.normal(0, SIG, len(t))
    y += rng.normal(0, SIG, len(t))
    x0 = np.array([3e-6, 30e-6])
    xp = x0 + np.outer(vx * t, [1, 1]) + rng.normal(0, SIG, (len(t), 2))
    yp = np.array([5e-6, 20e-6]) + np.outer(vy * t, [1, 1]) + rng.normal(0, SIG, (len(t), 2))
    return {'t': t, 'x': x, 'y': y, 'x_pegadas': xp, 'y_pegadas': yp}


def test_velocidad_de_las_pegadas():
    rng = np.random.default_rng(0)
    s = segmento(15e-6, rng)
    v = stokes.velocidad_platina(s['t'], s['x_pegadas'], s['y_pegadas'])
    assert v['vx'] == pytest.approx(15e-6, rel=1e-3)
    assert abs(v['vy']) < 5 * v['err_vy'] + 1e-12


@pytest.mark.parametrize('eje,kappa', [('x', KX), ('y', KY)])
def test_pendiente_gamma_sobre_kappa(eje, kappa):
    rng = np.random.default_rng(1)
    vs = np.array([-20, -10, 0, 10, 20]) * 1e-6
    r = stokes.analizar([segmento(v, rng, eje=eje) for v in vs], gamma=GAMMA,
                        kappa={'x': KX, 'y': KY})
    aj = r['ajuste'][eje]
    print(f'\nSTOKES {eje}: pendiente / (gamma/kappa) = {aj["pendiente"] * kappa / GAMMA:.4f} '
          f'± {aj["err_pendiente"] * kappa / GAMMA:.4f}, chi2_red = {aj["chi2_red"]:.2f}')
    assert aj['pendiente'] == pytest.approx(GAMMA / kappa, rel=0.03)
    assert abs(aj['pendiente'] - GAMMA / kappa) < 3.5 * aj['err_pendiente']
    assert r['kappa_stokes'][eje] == pytest.approx(kappa, rel=0.03)
    assert r['pendiente_esperada'][eje] == pytest.approx(GAMMA / kappa)
    assert aj['ordenada'] == pytest.approx(X_TRAMPA if eje == 'x' else Y_TRAMPA, abs=5e-9)


def test_velocidad_de_escape_y_fuerza_maxima():
    rng = np.random.default_rng(2)
    vs = np.array([0, 10, 20, 30, 40, 50, 60]) * 1e-6
    v_esc = 45e-6
    segs = [segmento(v, rng, dur=5, v_escape=v_esc) for v in vs]
    r = stokes.analizar(segs, gamma=GAMMA, kappa={'x': KX, 'y': KY}, radio_escape=1.5e-6)
    assert list(r['atrapada']) == [True] * 5 + [False] * 2
    e = r['escape']
    assert e['v_abajo'] == pytest.approx(40e-6, rel=1e-3)
    assert e['v_arriba'] == pytest.approx(50e-6, rel=1e-3)
    assert e['F_max'] == pytest.approx(GAMMA * 45e-6, rel=1e-3)
    # la pendiente usa solo los segmentos en que siguió atrapada
    assert r['ajuste']['x']['n'] == 5
    assert r['ajuste']['x']['pendiente'] == pytest.approx(GAMMA / KX, rel=0.05)


def test_escape_casos_borde():
    assert np.isnan(stokes.velocidad_escape([1, 2], [True, True])['v_max'])
    assert np.isnan(stokes.velocidad_escape([1, 2], [False, False])['v_max'])
    r = stokes.velocidad_escape([-3, 1, 2, 4], [False, True, True, False])
    assert (r['v_abajo'], r['v_arriba'], r['v_max']) == (2, 3, 2.5)
