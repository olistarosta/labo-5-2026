"""imagen: canales Bayer coherentes con a_coordenadas_sensor, demosaico, fondo, vista y retícula."""
import numpy as np
import pytest

from pinzas import imagen, io, sintetico

FASES = ['RGGB', 'GRBG', 'GBRG', 'BGGR']


def _meta(fase):
    return {'camara': {'formato_pixel': 'bayer' if fase else 'mono', 'bayer_fase': fase, 'bits': 10}}


def _mancha(forma, x0, y0, sigma=4.0):
    """1 - 0.5·gaussiana centrada en (x0, y0) px del sensor, evaluada en los centros de píxel."""
    yy, xx = np.mgrid[0:forma[0], 0:forma[1]].astype(float)
    return 1 - 0.5 * np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * sigma ** 2))


def _centroide(img):
    w = 1 - img / img.max()                 # sin umbral: un umbral sesga el centroide ~0.02 px
    yy, xx = np.mgrid[0:img.shape[0], 0:img.shape[1]]
    return np.sum(w * xx) / w.sum(), np.sum(w * yy) / w.sum()


@pytest.mark.parametrize('fase', FASES)
def test_planos_y_coordenadas_exactos(fase):
    """El píxel (v, u) de cada plano es el píxel a_coordenadas_sensor(u, v) del sensor."""
    rng = np.random.default_rng(0)
    rgb = rng.integers(0, 1024, (16, 20, 3)).astype(np.uint16)
    crudo = sintetico.mosaico_bayer(rgb, fase)
    for canal, k in (('rojo', 0), ('azul', 2)):
        p = imagen.preparar(crudo, _meta(fase), canal)
        assert p.shape == (8, 10) and p.dtype == np.float32
        v, u = np.mgrid[0:8, 0:10]
        xs, ys = imagen.a_coordenadas_sensor(u, v, canal, fase)
        np.testing.assert_array_equal(p, rgb[ys.astype(int), xs.astype(int), k])
    # verde: promedio de los dos G; su centro está en (2u + 1/2, 2v + 1/2)
    p = imagen.preparar(crudo, _meta(fase), 'verde')
    (a, b), (c, d) = imagen.posiciones_bayer('G', fase)
    np.testing.assert_allclose(p, 0.5 * (rgb[a::2, b::2, 1].astype(float) + rgb[c::2, d::2, 1]))
    assert imagen.a_coordenadas_sensor(0, 0, 'verde', fase) == (0.5, 0.5)
    # ida y vuelta
    x, y = imagen.a_coordenadas_canal(*imagen.a_coordenadas_sensor(3.3, 1.7, 'rojo', fase, (5, 3, 10, 10)),
                                      'rojo', fase, (5, 3, 10, 10))
    assert (float(x), float(y)) == pytest.approx((3.3, 1.7))


@pytest.mark.parametrize('fase', FASES)
@pytest.mark.parametrize('canal', ['rojo', 'verde', 'azul', 'gris'])
def test_posicion_de_una_mancha_en_cada_canal(fase, canal):
    """Una mancha en (x0, y0) del sensor se recupera en (x0, y0) pasando por cada canal."""
    x0, y0 = 31.37, 26.81
    gris = _mancha((64, 64), x0, y0)
    crudo = sintetico.mosaico_bayer(np.stack([0.8 * gris, gris, 0.6 * gris], -1) * 800, fase)
    p = imagen.preparar(crudo, _meta(fase), canal)
    xc, yc = _centroide(p)
    xs, ys = imagen.a_coordenadas_sensor(xc, yc, canal, fase)
    assert (float(xs), float(ys)) == pytest.approx((x0, y0), abs=0.01)


def test_recorte_con_fase_impar_se_lleva_a_pares():
    x0, y0 = 40.2, 30.6
    gris = _mancha((64, 80), x0, y0)
    crudo = sintetico.mosaico_bayer(np.stack([gris, gris, gris], -1) * 800, 'RGGB')
    rec = (21, 13, 41, 37)                                   # impares: queda (20, 12, 40, 36)
    assert imagen.recorte_efectivo(rec, crudo.shape, bayer=True) == (20, 12, 40, 36)
    p = imagen.preparar(crudo, _meta('RGGB'), 'rojo', recorte=rec)
    assert p.shape == (18, 20)
    xs, ys = imagen.a_coordenadas_sensor(*_centroide(p), 'rojo', 'RGGB', recorte=rec)
    assert (float(xs), float(ys)) == pytest.approx((x0, y0), abs=0.01)
    # mono: recorte exacto, coordenadas identidad + desplazamiento
    pm = imagen.preparar((gris * 800).astype(np.uint16), _meta(None), 'verde', recorte=rec)
    assert pm.shape == (37, 41)
    assert imagen.a_coordenadas_sensor(0, 0, 'verde', None, recorte=rec) == (21.0, 13.0)


def test_demosaico_reproduce_constantes_y_valores_medidos():
    rgb = np.ones((10, 12, 3)) * np.array([300.0, 500.0, 200.0])
    d = imagen.demosaico_bilineal(sintetico.mosaico_bayer(rgb, 'GBRG'), 'GBRG')
    np.testing.assert_allclose(d, rgb, rtol=1e-6)            # también en los bordes
    rng = np.random.default_rng(1)
    rgb = rng.uniform(0, 1000, (10, 12, 3))
    cr = sintetico.mosaico_bayer(rgb, 'RGGB')
    d = imagen.demosaico_bilineal(cr, 'RGGB')
    np.testing.assert_allclose(d[0::2, 0::2, 0], cr[0::2, 0::2], rtol=1e-6)   # R medido queda igual


def test_px_sensor_por_px():
    assert imagen.px_sensor_por_px('verde', _meta('RGGB')) == 2
    assert imagen.px_sensor_por_px('gris', _meta('RGGB')) == 1
    assert imagen.px_sensor_por_px('verde', _meta(None)) == 1


def test_grabacion_bayer_fondo_y_vista(tmp_path):
    beads = [dict(x0_px=40.0, y0_px=30.0, D_m2s=2e-12, radio_px=6.0)]
    car = sintetico.grabar(str(tmp_path), 'bayer', 'libre', n=40, fps=20.0, W=5e-3, beads=beads,
                           forma=(64, 80), bayer='RGGB', prob_perdido=0.0, seed=2)
    g = io.leer_grabacion(car)
    assert g.meta['camara']['formato_pixel'] == 'bayer'
    v = imagen.preparar(g.cuadro(0), g.meta, 'verde')
    assert v.shape == (32, 40)
    f = imagen.fondo_mediana(g, canal='verde', n_max=21)
    p = imagen.preparar(g.cuadro(5), g.meta, 'verde', fondo=f)
    # lejos de la bead el fondo dividido queda ~1 (la iluminación despareja se fue)
    xs, ys = g.meta['verdad']['x_px'][5][0], g.meta['verdad']['y_px'][5][0]
    xc, yc = imagen.a_coordenadas_canal(xs, ys, 'verde', 'RGGB')
    yy, xx = np.mgrid[0:32, 0:40]
    lejos = np.hypot(xx - xc, yy - yc) > 10
    assert abs(np.median(p[lejos]) - 1) < 0.01
    assert np.std(p[lejos]) < 0.06                           # ruido de disparo de ~2·600 cuentas
    pr = imagen.preparar(g.cuadro(5), g.meta, 'verde', fondo=f, modo_fondo='restar')
    assert abs(np.median(pr[lejos])) < 10
    with pytest.raises(ValueError):
        imagen.preparar(g.cuadro(5), g.meta, 'gris', fondo=f)
    # vista rápida
    vr = imagen.vista_rapida(g.cuadro(0), g.meta, max_lado=20)
    assert vr.dtype == np.uint8 and vr.shape == (16, 20, 3)
    vm = imagen.vista_rapida(np.zeros((100, 300), np.uint16) + 7, _meta(None), max_lado=100)
    assert vm.shape == (33, 100) and vm.dtype == np.uint8


@pytest.mark.parametrize('tipo', ['lineas', 'puntos'])
@pytest.mark.parametrize('paso_px,angulo', [(11.3, 0.0), (17.77, 0.06), (8.2, -0.3)])
def test_escala_reticula(tipo, paso_px, angulo):
    """Escala recuperada < 0.5 % (criterio); se imprime el error medido."""
    crudo = sintetico.reticula((256, 300), paso_px, angulo=angulo, tipo=tipo, rng=3,
                               desplazamiento_px=(3.1, 5.4))
    paso_m = 10e-6
    r = imagen.escala_reticula(crudo.astype(float), paso_m)
    verdad = paso_m / paso_px
    err = r['m_por_px'] / verdad - 1
    print(f'\n[reticula {tipo} P={paso_px} ang={angulo}] error {err * 100:+.4f} %, '
          f'err estimado {r["err_m_por_px"] / verdad * 100:.4f} %, P = {r["periodo_px"]}, '
          f'ang = {r["angulo_rad"]:.4f}')
    assert abs(err) < 0.005
    np.testing.assert_allclose(r['periodo_px'], paso_px, rtol=0.005)
    assert r['angulo_rad'] == pytest.approx(angulo, abs=2e-3)
    assert r['err_m_por_px'] < 0.005 * verdad
    assert r['lineas_oscuras']
    assert len(r['lineas_x']) >= 256 // paso_px and r['lineas_x'].shape[1] == 4


def test_lineas_detectadas_caen_sobre_las_lineas():
    paso = 12.5
    crudo = sintetico.reticula((128, 128), paso, tipo='lineas', desplazamiento_px=(4.0, 7.0),
                               ruido_lectura=0.0, rng=1)
    r = imagen.escala_reticula(crudo.astype(float), 1e-5)
    xs = np.sort(r['lineas_x'][:, 0])                        # líneas verticales: x constante
    esperadas = 4.0 + paso * np.arange(-1, 12)
    esperadas = esperadas[(esperadas >= 0) & (esperadas <= 127)]
    np.testing.assert_allclose(xs, esperadas, atol=0.1)


def test_escala_reticula_en_plano_bayer():
    """Retícula vista por el plano verde (media resolución): la escala sigue siendo por px del sensor."""
    paso_px = 13.1
    crudo = sintetico.reticula((256, 256), paso_px, angulo=0.02, bayer='RGGB', rng=4)
    meta = _meta('RGGB')
    v = imagen.preparar(crudo, meta, 'verde')
    r = imagen.escala_reticula(v, 10e-6, px_sensor_por_px=imagen.px_sensor_por_px('verde', meta))
    err = r['m_por_px'] / (10e-6 / paso_px) - 1
    print(f'\n[reticula bayer verde] error {err * 100:+.4f} %')
    assert abs(err) < 0.005
