"""Tests de TRAMPA sobre trayectorias sintéticas exactas (OU promediado durante W + ruido)."""
import numpy as np
import pytest

from pinzas import fisica, trampa
from pinzas.fisica import KB, factor_blur_covarianza as B, factor_blur_varianza as S
from pinzas.sintetico import trayectoria

T = 296.15
GAMMA = fisica.gamma_stokes(1.7e-3, 1.5e-6)     # PS 3 µm, 20 % glicerina
D = fisica.difusion(T, GAMMA)
KX, KY = 3e-6, 5e-6                            # 3 y 5 pN/µm: trampa anisótropa
SIG = 0.05 * 90e-9                             # ruido de localización 0.05 px


def atrapada(semilla, fps=100.0, W=5e-3, dur=300.0, kx=KX, ky=KY, perdidos=30, sig=SIG):
    rng = np.random.default_rng(semilla)
    n = int(dur * fps)
    t = np.arange(n) / fps
    keep = np.ones(n, bool)
    keep[rng.choice(np.arange(1, n - 1), perdidos, replace=False)] = False
    t = t[keep]
    x = trayectoria(t, W, D, kx, GAMMA, T, 40, rng)[0] + rng.normal(0, sig, len(t))
    y = trayectoria(t, W, D, ky, GAMMA, T, 40, rng)[0] + rng.normal(0, sig, len(t))
    # bead pegada: mismo ruido, más una deriva lenta que el método de diferencias ignora
    tp = np.arange(5000) / fps
    ruido = trampa.piso_de_ruido(rng.normal(0, sig, 5000) + 2e-9 * tp, rng.normal(0, sig, 5000))
    return t, x, y, ruido


# ---------------------------------------------------------------- fórmulas
def test_factores_de_blur():
    assert S(0.0) == pytest.approx(1.0)
    assert B(0.0) == pytest.approx(1.0)
    a = np.array([1e-5, 1e-4 * 0.999, 1e-4 * 1.001, 0.1, 1.5, 3.0])
    exacto = 2 / a - 2 / a ** 2 * (1 - np.exp(-a))
    assert np.allclose(S(a[3:]), exacto[3:], rtol=1e-12)
    assert np.allclose(S(a[:3]), 1 - a[:3] / 3, rtol=1e-6)              # continuidad de la serie
    assert np.allclose(B(a[3:]), 2 * (np.cosh(a[3:]) - 1) / a[3:] ** 2, rtol=1e-10)
    assert np.allclose(B(a[:3]), 1 + a[:3] ** 2 / 12, rtol=1e-12)
    assert np.all(np.diff(S(np.linspace(0, 5, 50))) < 0)


def test_psd_forma_cerrada_igual_a_suma_de_alias():
    dt, tau_c, v, W = 0.01, 0.016, 1.4e-15, 0.006
    f = np.linspace(0.5, 49.5, 30)
    cerrada = trampa.modelo_psd(f, dt, v, tau_c, 0.0, W, filtro=True)
    suma = trampa.psd_suma_alias(f, dt, v, tau_c, W, n_alias=20000)
    assert np.allclose(cerrada, suma, rtol=2e-4)
    # sin exposición, la Lorentziana con aliasing de BS&F
    assert np.allclose(trampa.modelo_psd(f, dt, v, tau_c, 0.0, 0.0, True),
                       trampa.modelo_psd(f, dt, v, tau_c, 0.0, 0.0, False), rtol=1e-12)


def test_piso_de_ruido_ignora_deriva_lenta():
    rng = np.random.default_rng(0)
    n = 20000
    x = rng.normal(0, SIG, n) + np.linspace(0, 20 * SIG, n)
    r = trampa.piso_de_ruido(x, rng.normal(0, SIG, n))
    assert r['var_x'] == pytest.approx(SIG ** 2, rel=0.04)
    assert r['var_x_directa'] > 10 * SIG ** 2


# ---------------------------------------------------------------- criterios del encargo
@pytest.fixture(scope='module')
def tres_videos():
    out = []
    for s in range(3):
        t, x, y, ruido = atrapada(s)
        out.append(trampa.analizar(t, x, y, 5e-3, T, GAMMA, ruido=ruido, n_boot=150, semilla=s))
    return out


def test_kappa_con_correcciones_dentro_del_5_por_ciento(tres_videos):
    tab = tres_videos[0]['tabla']
    verdad = np.where(tab['eje'] == 'x', KX, KY)
    rel = np.array([r['tabla']['kappa'] / verdad - 1 for r in tres_videos])        # [3, 10]
    z = np.array([(r['tabla']['kappa'] - verdad) / r['tabla']['err_kappa'] for r in tres_videos])
    print('\nTRAMPA con correcciones (media de 3 videos de 300 s, 100 fps, W = 5 ms):')
    for i in range(10):
        print(f"  {tab['metodo'][i]:16s} {tab['eje'][i]}  {rel[:, i].mean():+.4f}   z = "
              + ' '.join(f'{v:+.2f}' for v in z[:, i]))
    assert np.all(np.abs(rel.mean(0)) < 0.05)
    assert np.all(np.isfinite(tab['err_kappa']))
    assert np.all(np.abs(z) < 3.5)                       # cada uno compatible con la verdad...
    for r in tres_videos:                                # ...y los cinco entre sí
        for eje in 'xy':
            assert r['compatibilidad'][eje]['z_max'] < 3.5
    # tau_c de autocorrelación, espectro y MSD
    tc = tab['tau_c']
    ok = np.isfinite(tc)
    assert np.allclose(np.mean([r['tabla']['tau_c'][ok] for r in tres_videos], 0),
                       (GAMMA / verdad)[ok], rtol=0.05)


def test_sin_correccion_de_blur_sesga_como_predice_S():
    # alpha = W / tau_c = 1.5: S = 0.643, kappa sin corregir blur = kappa / S = +56 %
    W, fps = 24e-3, 40.0
    alpha = W * KX / GAMMA
    t, x, _, ruido = atrapada(7, fps=fps, W=W, dur=300, ky=KX)
    var = np.var(x)
    k_sin, _ = trampa.kappa_varianza(var, T, ruido['var_x'], W, GAMMA, corregir_blur=False)
    k_con, a = trampa.kappa_varianza(var, T, ruido['var_x'], W, GAMMA, corregir_blur=True)
    predicho = 1 / float(S(alpha))
    print(f'\nsin blur: kappa/kappa_v = {k_sin / KX:.3f}, predicho 1/S(alpha={alpha:.2f}) = {predicho:.3f};'
          f' con corrección {k_con / KX:.3f}')
    assert k_sin / KX > 1.4                                     # dirección: sobreestima
    assert k_sin / KX == pytest.approx(predicho, rel=0.03)      # magnitud: la de S(alpha)
    assert k_con / KX == pytest.approx(1, abs=0.03)
    assert a == pytest.approx(alpha, rel=0.03)
    # el análisis completo: tabla sin corrección sobreestima en equipartición/Boltzmann como 1/(S + s2/v)
    r = trampa.analizar(t, x, x, W, T, GAMMA, ruido=ruido, n_boot=0)
    esperado = 1 / (float(S(alpha)) + ruido['var_x'] * KX / (KB * T))
    for m in ('equiparticion', 'boltzmann'):
        i = list(r['tabla_sin_correccion']['metodo']).index(m)
        assert r['tabla_sin_correccion']['kappa'][i] / KX == pytest.approx(esperado, rel=0.03)
        assert r['tabla']['kappa'][i] / KX == pytest.approx(1, abs=0.03)


def test_boltzmann_detecta_potencial_no_armonico():
    rng = np.random.default_rng(3)
    s = np.sqrt(KB * T / KX)

    def muestras(beta):              # U = kappa x²/2 + beta x⁴ por rechazo
        x = rng.normal(0, s, 400000)
        acepta = rng.random(len(x)) < np.exp(-beta * x ** 4 / (KB * T))
        return x[acepta][:100000]

    bordes = np.arange(-5, 5.0001, 0.2) * s
    arm = trampa.ajuste_boltzmann(bordes, np.histogram(muestras(0.0), bordes)[0].astype(float), T)
    anarm = trampa.ajuste_boltzmann(bordes, np.histogram(muestras(0.3 * KB * T / s ** 4), bordes)[0]
                                    .astype(float), T)
    assert arm['kappa_raw'] == pytest.approx(KX, rel=0.02)
    assert arm['anarmonia_2s'] < 0.05 and arm['chi2_red'] < 2
    assert anarm['anarmonia_2s'] > 0.1


def test_sin_gamma_da_tau_c_pero_no_kappa():
    t, x, y, ruido = atrapada(11, dur=120)
    r = trampa.analizar(t, x, y, 5e-3, T, None, ruido=ruido, n_boot=0, sin_correccion=False)
    tab = r['tabla']
    i = (tab['metodo'] == 'autocorrelacion') & (tab['eje'] == 'x')
    assert np.isnan(tab['kappa'][i]).all() and np.isfinite(tab['tau_c'][i]).all()
    j = (tab['metodo'] == 'equiparticion') & (tab['eje'] == 'x')
    assert tab['kappa'][j][0] == pytest.approx(KX, rel=0.06)      # blur con tau_c de la autocorrelación
    assert any('gamma' in a for a in r['avisos'])


def test_graficos_completos():
    t, x, y, ruido = atrapada(12, dur=60)
    r = trampa.analizar(t, x, y, 5e-3, T, GAMMA, ruido=ruido, n_boot=0, sin_correccion=False)
    g = r['graficos']['x']
    for c in ('x_hist', 'P', 'U', 'U_ajuste', 'tau', 'C', 'C_ajuste', 'msd', 'msd_ajuste',
              'f', 'psd', 'psd_ajuste'):
        assert len(g[c]) > 5 and np.isfinite(g[c]).any(), c
    assert len(g['f']) == len(g['psd']) == len(g['psd_ajuste'])


def test_pasaaltos_contra_deriva_lenta():
    # deriva de 2 sigma en 300 s más una oscilación de 0.01 Hz; corte 0.05 Hz << f_c = 10 Hz
    t, x, y, ruido = atrapada(13)
    s = np.sqrt(KB * T / KX)
    deriva = 2 * s * t / t[-1] + s * np.sin(2 * np.pi * 0.01 * t)
    sin_f = trampa.analizar(t, x + deriva, y, 5e-3, T, GAMMA, ruido=ruido, n_boot=0,
                            sin_correccion=False)
    con_f = trampa.analizar(t, x + deriva, y, 5e-3, T, GAMMA, ruido=ruido, n_boot=0,
                            sin_correccion=False, pasaaltos_hz=0.05)
    i = con_f['tabla']['eje'] == 'x'
    eq = i & (con_f['tabla']['metodo'] == 'equiparticion')
    assert sin_f['tabla']['kappa'][eq][0] / KX < 0.7                 # la deriva infla la varianza
    assert np.all(np.abs(con_f['tabla']['kappa'][i] / KX - 1) < 0.06)
    assert con_f['parametros']['f_min_hz'] == pytest.approx(0.2)
