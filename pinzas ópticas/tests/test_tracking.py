"""Tests de pinzas.tracking sobre sintéticos (sintetico.dibujar / camara / trayectoria).

Criterios del encargo: bead quieta con sesgo < 0.02 px; ruido estimado de un video de bead
quieta igual al verdadero ± 10 %; `seguir` recupera una OU de rms ~0.4 px; deriva 'media'.
"""
import numpy as np
import pytest

from pinzas import sintetico as sint
from pinzas import tracking as tr
from pinzas.fisica import KB

R = 8.0                 # px: bead de ~1.5 µm a 90 nm/px
NIVEL, CONTRASTE = 600.0, 0.35
METODOS = list(tr.METODOS)


def fondo_de(forma, R=R):
    """Iluminación despareja por defecto de sintetico.dibujar (sin beads)."""
    return sint.dibujar(forma, np.zeros((0, 2)), R)


def preparar(ideal, fondo, rng, nivel=NIVEL):
    """Cuadro con ruido de cámara y fondo restado (lo que hará imagen.py): c / (nivel·fondo) - 1."""
    return sint.camara(ideal, nivel, rng=rng) / (nivel * fondo) - 1.0


# ---------------------------------------------------------------- detectar
def test_detectar_encuentra_las_beads_con_fondo_disparejo():
    rng = np.random.default_rng(0)
    forma = (160, 220)
    verdad = np.array([[40.3, 35.7], [120.6, 50.2], [180.1, 120.9], [60.8, 125.4]])
    crudo = sint.camara(sint.dibujar(forma, verdad, R, contraste=CONTRASTE), NIVEL, rng=rng)
    cand = tr.detectar(crudo.astype(float), R)             # sin restar fondo a propósito
    assert cand.shape == (4, 3)
    for x, y in verdad:
        d = np.hypot(cand[:, 0] - x, cand[:, 1] - y)
        assert d.min() < 1.0
        assert 0.6 * R < cand[np.argmin(d), 2] < 1.5 * R


def test_detectar_cuadro_vacio_no_inventa():
    rng = np.random.default_rng(1)
    forma = (120, 120)
    crudo = sint.camara(fondo_de(forma), NIVEL, rng=rng).astype(float)
    assert tr.detectar(crudo, R).shape == (0, 3)


# ---------------------------------------------------------------- localizadores
@pytest.mark.parametrize('metodo,suav', [('simetria_radial', 0.0), ('simetria_radial', 0.2 * R),
                                         ('gauss2d', 0.0)])
def test_sesgo_sin_ruido_en_grilla_subpixel(metodo, suav):
    """Imagen ideal (super-muestreo 8), posiciones 0..1 px en x e y: solo el sesgo del método."""
    N = 64
    fondo = np.ones((N, N))
    lado = tr.lado_recorte_sugerido(R)
    err = []
    for fx in np.linspace(0, 1, 11):
        for fy in np.linspace(0, 1, 11):
            x0, y0 = 32 + fx, 32 + fy
            img = sint.dibujar((N, N), [[x0, y0]], R, fondo=fondo, super_muestreo=8) - 1
            x, y = tr._localizar_en(img, 32.2, 31.9, lado, metodo, True, suav)
            err.append((x - x0, y - y0))
    assert np.abs(err).max() < 0.02


@pytest.mark.parametrize('metodo', METODOS)
def test_sesgo_con_ruido_bead_quieta(metodo):
    """Grilla 4×4 de posiciones subpíxel, 300 realizaciones cada una, vía `seguir`.

    sesgo = |<x_medido> - x_verdadero| < 0.02 px en cada posición y eje. (El error de estimar
    la media con 300 muestras es ~0.003 px.)
    """
    rng = np.random.default_rng(2)
    N, M = 62, 300
    fondo = fondo_de((N, N))
    for fx in (0.0, 0.25, 0.5, 0.75):
        for fy in (0.0, 0.25, 0.5, 0.75):
            x0, y0 = 31 + fx, 31 + fy
            ideal = sint.dibujar((N, N), [[x0, y0]], R, contraste=CONTRASTE, fondo=fondo)
            cuadros = [preparar(ideal, fondo, rng) for _ in range(M)]
            r = tr.seguir(lambda i: cuadros[i], M, [[x0 + 0.4, y0 - 0.4]], metodo, radio_px=R)
            assert r['valido'].all()
            assert abs(r['x'].mean() - x0) < 0.02
            assert abs(r['y'].mean() - y0) < 0.02


@pytest.mark.parametrize('metodo', METODOS)
def test_ruido_de_un_video_de_bead_pegada_igual_al_verdadero(metodo):
    """Ruido de localización: 'verdadero' = desvío sobre 1000 realizaciones independientes en
    posiciones subpíxel al azar; 'estimado' = desvío temporal de UN video de 1000 cuadros de una
    bead pegada que deriva lento (0.8 px en el video), con seguir + quitar_deriva('recta').
    Tienen que coincidir ± 10 % en cada eje."""
    rng = np.random.default_rng(3)
    N, n = 62, 1000
    fondo = fondo_de((N, N))
    lado = tr.lado_recorte_sugerido(R)
    suav = 0.2 * R if metodo == 'simetria_radial' else 0.0
    err = []
    for _ in range(n):
        x0, y0 = 31 + rng.random(), 31 + rng.random()
        img = preparar(sint.dibujar((N, N), [[x0, y0]], R, contraste=CONTRASTE, fondo=fondo),
                       fondo, rng)
        err.append(np.subtract(tr._localizar_en(img, x0, y0, lado, metodo, True, suav), (x0, y0)))
    verdadero = np.std(err, axis=0)

    xs = 30.2 + 0.8 * np.arange(n) / n
    ys = 31.7 - 0.5 * np.arange(n) / n

    def cuadro(i):
        ideal = sint.dibujar((N, N), [[xs[i], ys[i]]], R, contraste=CONTRASTE, fondo=fondo)
        return preparar(ideal, fondo, np.random.default_rng(1000 + i))

    r = tr.seguir(cuadro, n, [[xs[0], ys[0]]], metodo, radio_px=R)
    xc, yc, _, _ = tr.quitar_deriva(r['x'], r['y'], 'recta')
    estimado = np.array([np.nanstd(xc), np.nanstd(yc)])
    assert np.all(np.abs(estimado / verdadero - 1) < 0.10), (estimado, verdadero)


# ---------------------------------------------------------------- seguir
def test_seguir_recupera_una_trayectoria_OU():
    """Dos beads atrapadas (OU, κ = 3 pN/µm, γ = 5e-8 N·s/m, 100 fps, W = 5 ms): rms ~0.4 px."""
    rng = np.random.default_rng(4)
    m_por_px, T = 90e-9, 296.15
    kappa, gamma, W = 3e-6, 5e-8, 5e-3
    D = KB * T / gamma
    n = 600
    t = np.arange(n) / 100.0
    verdad = np.zeros((n, 2, 2))                         # [cuadro, bead, (x, y)] px
    centros = [(40.3, 40.6), (100.7, 45.2)]
    for b in range(2):
        for e in range(2):
            x_m, _ = sint.trayectoria(t, W, D, kappa=kappa, gamma=gamma, rng=rng)
            verdad[:, b, e] = centros[b][e] + x_m / m_por_px
    forma = (90, 140)
    fondo = fondo_de(forma)
    cuadros = [preparar(sint.dibujar(forma, verdad[i], R, contraste=CONTRASTE, fondo=fondo),
                        fondo, rng) for i in range(n)]
    llamadas = []
    for metodo in METODOS:
        r = tr.seguir(lambda i: cuadros[i], n, verdad[0], metodo, radio_px=R,
                      progreso=lambda i, n: llamadas.append(i))
        assert r['valido'].all()
        senal = np.std(verdad[:, :, 0], axis=0).mean()
        assert 0.3 < senal < 0.6                          # el caso del encargo
        for X, e in ((r['x'], 0), (r['y'], 1)):
            err = X - verdad[:, :, e]
            assert np.sqrt(np.mean(err ** 2)) < 0.25 * senal
            assert abs(np.mean(err)) < 0.02
        assert r['ms_por_bead_cuadro'] > 0
    assert llamadas == list(range(n)) * len(METODOS)


def test_seguir_marca_NaN_si_salta_o_sale_del_cuadro():
    rng = np.random.default_rng(5)
    forma = (70, 120)
    fondo = fondo_de(forma)
    n = 40
    pos = np.zeros((n, 2, 2))
    pos[:, 0] = [30.0, 35.0]
    pos[20:, 0] = [42.0, 35.0]                           # salto de 12 px en el cuadro 20
    pos[:, 1, 0] = 80 + 2.0 * np.arange(n)               # se va por el borde derecho
    pos[:, 1, 1] = 35.0
    cuadros = [preparar(sint.dibujar(forma, pos[i], R, fondo=fondo), fondo, rng)
               for i in range(n)]
    r = tr.seguir(lambda i: cuadros[i], n, pos[0], radio_px=R)
    assert r['valido'][:20, 0].all() and not r['valido'][20:, 0].any()
    sale = 80 + 2.0 * np.arange(n) + r['lado_recorte_px'] // 2 + 1 > forma[1]
    assert not r['valido'][sale, 1].any()
    assert np.isnan(r['x'][~r['valido']]).all()
    assert np.nanmin(r['contraste_z'][r['valido']]) > 20      # las beads de verdad, holgadas


# ---------------------------------------------------------------- deriva
def _libres(nb, n, D_px2, rng, frac_nan=0.0):
    """nb beads libres (random walk, D en px²/cuadro) más NaN al azar."""
    x = np.cumsum(rng.normal(0, np.sqrt(2 * D_px2), (n, nb)), axis=0) + rng.uniform(20, 200, nb)
    y = np.cumsum(rng.normal(0, np.sqrt(2 * D_px2), (n, nb)), axis=0) + rng.uniform(20, 200, nb)
    if frac_nan:
        malo = rng.random((n, nb)) < frac_nan
        x[malo] = np.nan
        y[malo] = np.nan
    return x, y


def test_quitar_deriva_media_recupera_una_deriva_lineal():
    rng = np.random.default_rng(6)
    nb, n, D = 8, 3000, 0.05
    x, y = _libres(nb, n, D, rng, frac_nan=0.05)
    i = np.arange(n)
    vx, vy = 0.01, -0.004                                # px / cuadro: 30 y 12 px en el video
    xc, yc, dx, dy = tr.quitar_deriva(x + vx * i[:, None], y + vy * i[:, None], 'media')
    assert dx.shape == (n,) and dx[0] == 0
    # la deriva estimada = inyectada + promedio de los brownianos de las nb beads
    # (desvío ~ sqrt(2 D t / nb) ~ 6 px al final): la pendiente sale a ~1e-3 px/cuadro
    assert abs(np.polyfit(i, dx, 1)[0] - vx) < 1.5e-3
    assert abs(np.polyfit(i, dy, 1)[0] - vy) < 1.5e-3
    # en los datos corregidos no queda deriva común: el desplazamiento medio por paso es 0
    paso = np.nanmean(np.diff(xc, axis=0), axis=1)
    assert np.nanmax(np.abs(np.cumsum(np.nan_to_num(paso)))) < 1e-9
    assert np.isnan(xc[np.isnan(x)]).all()               # NaN sigue siendo NaN


def test_quitar_deriva_media_le_quita_1_sobre_nb_al_browniano():
    """Documenta el sesgo de 'media' sin referencia: var del paso corregido = 2D (1 - 1/nb)."""
    rng = np.random.default_rng(7)
    nb, n, D = 5, 20000, 0.05
    x, y = _libres(nb, n, D, rng)
    xc, _, _, _ = tr.quitar_deriva(x, y, 'media')
    var = np.var(np.diff(xc, axis=0))
    assert abs(var / (2 * D * (1 - 1 / nb)) - 1) < 0.03


def test_quitar_deriva_con_beads_de_referencia_no_sesga_la_libre():
    rng = np.random.default_rng(8)
    n, D = 20000, 0.05
    i = np.arange(n)
    deriva = 0.002 * i + 3 * np.sin(2 * np.pi * i / 7000)
    libre, _ = _libres(1, n, D, rng)
    pegadas = rng.normal(0, 0.05, (n, 2)) + [50.0, 90.0]  # solo ruido de localización
    x = np.column_stack([libre[:, 0], pegadas]) + deriva[:, None]
    xc, _, dx, _ = tr.quitar_deriva(x, x.copy(), 'media', referencia=[1, 2])
    assert np.std(dx - (deriva - deriva[0])) < 0.05
    var = np.var(np.diff(xc[:, 0]))
    assert abs(var / (2 * D) - 1) < 0.03


def _ou_con_huecos(rng, n=30000, fps=200.0):
    """OU en px (rms 0.4, f_c = 10 Hz) con 2 % de cuadros perdidos (tiempos reales con huecos)."""
    m_por_px, T, kappa, gamma = 90e-9, 296.15, 3e-6, 5e-8
    t = np.arange(n) / fps
    t = t[rng.random(n) > 0.02]
    x_m, _ = sint.trayectoria(t, 2e-3, KB * T / gamma, kappa=kappa, gamma=gamma, rng=rng)
    return t, x_m / m_por_px


def test_quitar_deriva_recta_por_bead():
    rng = np.random.default_rng(9)
    t, x = _ou_con_huecos(rng, n=10000)
    con = x + 0.05 * t + 7.0
    xc, _, dx, _ = tr.quitar_deriva(con[:, None], con[:, None], 'recta', t=t)
    assert dx.shape == (t.size, 1)
    assert abs(np.std(xc) / np.std(x) - 1) < 0.01
    assert abs(np.mean(xc)) < 1e-9


def test_quitar_deriva_pasaaltos_con_huecos_conserva_la_varianza_OU():
    """Deriva lenta (0.01 Hz, 2 px) + OU con f_c = 10 Hz; fc = 0.1 Hz (f_c/100)."""
    rng = np.random.default_rng(10)
    t, x = _ou_con_huecos(rng)
    lenta = 2.0 * np.sin(2 * np.pi * 0.01 * t) + 0.02 * t
    con = (x + lenta)[:, None]
    xc, _, dx, _ = tr.quitar_deriva(con, con, 'pasaaltos', fc_hz=0.1, t=t)
    lejos = (t > t[0] + 10) & (t < t[-1] - 10)            # sin los bordes del suavizado
    assert np.std((dx[:, 0] - lenta)[lejos]) < 0.05
    assert abs(np.std(xc[lejos, 0]) / np.std(x[lejos]) - 1) < 0.03
    with pytest.raises(ValueError):
        tr.quitar_deriva(con, con, 'pasaaltos', t=t)


def test_sigma_t_pasaaltos_da_1_sobre_raiz2_en_fc():
    fc = 0.3
    s = tr.sigma_t_pasaaltos(fc)
    assert np.isclose(np.exp(-2 * np.pi ** 2 * s ** 2 * fc ** 2), 1 / np.sqrt(2))


# ---------------------------------------------------------------- láser y unidades
@pytest.mark.parametrize('metodo', METODOS)
def test_mancha_de_laser_concentrica_no_corre_y_corrida_empuja(metodo):
    """Mancha gaussiana brillante (0.3 del fondo, s = 3 px). Concéntrica: simétrica, no corre
    el centro. A 3 px en +x: llena el núcleo oscuro de ese lado y el centro medido se ALEJA
    (-x), por décimas de px: comparable a la señal (cuánto: TouchDesigner/herramientas/comparar_tracking.py)."""
    N = 64
    fondo = np.ones((N, N))
    yy, xx = np.mgrid[0:N, 0:N]
    x0, y0 = 32.3, 31.6
    lado = tr.lado_recorte_sugerido(R)
    suav = 0.2 * R if metodo == 'simetria_radial' else 0.0
    bead = sint.dibujar((N, N), [[x0, y0]], R, contraste=CONTRASTE, fondo=fondo) - 1
    for d, cond in ((0.0, lambda e: abs(e) < 0.01), (3.0, lambda e: e < -0.05)):
        mancha = 0.3 * np.exp(-((xx - x0 - d) ** 2 + (yy - y0) ** 2) / (2 * 3.0 ** 2))
        x, _ = tr._localizar_en(bead + mancha, x0, y0, lado, metodo, True, suav)
        assert cond(x - x0), (d, x - x0)


def test_px_a_m():
    assert np.allclose(tr.px_a_m([1.0, np.nan, -2.0], 9e-8), [9e-8, np.nan, -1.8e-7],
                       equal_nan=True)
