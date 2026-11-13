"""importar: grabaciones de ThorCam (TIFF multipágina con tags TSI, partido en _0/_1), secuencia
numerada, RGB, OME-TIFF y AVI se importan SIN copiar ni tocar los originales, y la cadena anda
igual que sobre la grabación nuestra con los mismos cuadros."""
import hashlib
import os
import warnings

import numpy as np
import pytest
import tifffile

from pinzas import imagen, importar, io, sintetico, tiempo

MPX = 9e-8
NOMBRE_TC = 'image_2026-11-05T10-00-00.123'          # como los nombra ThorCam: prefijo + sello


def _tags_tsi(contador, ts_ns, roi, bits=10):
    """Tags privados por página como los documenta el Apéndice D del manual de ThorCam."""
    x0, y0, w, h = roi
    v = {32768: bits, 32770: 1, 32771: 1, 32772: x0, 32773: y0, 32774: w, 32775: h,
         32777: 40_000_000, 32778: 1}
    if contador is not None:
        v[32779] = int(contador)
    if ts_ns is not None:
        v[32781], v[32782] = int(ts_ns) >> 32, int(ts_ns) & 0xFFFFFFFF
    return [(k, 'I', 1, val, True) for k, val in sorted(v.items())]


def _escribir_tiff(ruta, cuadros, tags=None, photometric='minisblack'):
    """Multipágina: una página por cuadro (TiffWriter.write en un bucle), con tags por página."""
    with tifffile.TiffWriter(ruta) as tw:
        for k, img in enumerate(cuadros):
            tw.write(img, photometric=photometric, metadata=None, contiguous=False,
                     extratags=tags[k] if tags else ())


def _sintetica(tmp_path, bayer='RGGB', n=240, perdido=0.02, seed=3):
    b = dict(x0_px=32.3, y0_px=30.8, D_m2s=2e-13, kappa_x=3e-6, kappa_y=2e-6, radio_px=6.0)
    c = sintetico.grabar(str(tmp_path / 'orig'), 'orig', 'atrapada', n, 200.0, 1e-3, [b],
                         m_por_px=MPX, forma=(64, 64), bayer=bayer, prob_perdido=perdido, seed=seed)
    m = io.leer_json(os.path.join(c, 'grabacion.json'))
    m['muestra']['diametro_m'] = 3e-6
    io.escribir_grabacion_json(c, m)
    return io.leer_grabacion(c)


def _hash(rutas):
    return {r: hashlib.sha1(open(r, 'rb').read()).hexdigest() for r in rutas}


_BLOQ = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'TouchDesigner', 'bloques')


def _tiempo(carpeta):
    """TIEMPO de una grabación (lo mismo que usan los bloques)."""
    g = io.leer_grabacion(carpeta)
    return tiempo.reconstruir(g.tiempos, tick_hz=(g.meta.get('camara') or {}).get('tick_hz'))


def _tracking(carpeta, radio=3.0):
    """Posiciones de los bloques B0–B6 (sin fondo) en cada cuadro: lista de arrays [k, 3]."""
    import sys
    if _BLOQ not in sys.path:
        sys.path.insert(0, _BLOQ)
    import _lote
    import _td_bloque as tb
    etapas = []
    for n in ('b0_cuadro', 'b1_canal', 'b2_fondo', 'b3_fourier', 'b4_clahe', 'b5_realce', 'b6_deteccion',
              'b7_localizacion'):
        m = _lote.cargar(os.path.join(_BLOQ, n + '.py'))
        pp = {k: v for k, t, v, *_ in m.PARAMETROS if t != 'Pulse'}
        pp.update({'Formato': 'grabacion', 'Carpeta': carpeta, 'Metodo': 'no', 'Radio': radio,
                   'Chico': 1.0, 'Grande': 20.0, 'Umbral': 6.0})      # beads chicas (los defaults son para MP4)
        etapas.append((m, pp, {}))
    ctx = tb.Ctx(carpeta)
    return [tb.correr(i, etapas, ctx)['posiciones'] for i in range(len(ctx.grabacion))]


META = dict(tipo='atrapada', muestra=dict(diametro_m=3e-6, glicerina_pct=20.0),
            camara=dict(modelo='CS165CU', exposicion_us=1000.0, periodo_us=5000.0),
            escala_m_por_px=MPX)


def test_thorcam_multipagina_en_dos_partes_y_tracking(tmp_path, monkeypatch):
    monkeypatch.setenv('PINZAS_DATOS', str(tmp_path / 'datos'))
    g0 = _sintetica(tmp_path)
    n = len(g0)
    src = tmp_path / 'thorcam'
    src.mkdir()
    ti = g0.tiempos
    cuadros = [g0.cuadro(i) for i in range(n)]
    tags = [_tags_tsi(ti['contador_camara'][i], ti['t_camara'][i], (0, 0, 64, 64)) for i in range(n)]
    partes = [src / f'{NOMBRE_TC}_0.tif', src / f'{NOMBRE_TC}_1.tif']
    _escribir_tiff(str(partes[0]), cuadros[:150], tags[:150])
    _escribir_tiff(str(partes[1]), cuadros[150:], tags[150:])
    antes = _hash([str(p) for p in partes])

    car = importar.importar(str(partes[0]), dict(META, etiqueta='tc'))     # junta _0 y _1
    g = io.leer_grabacion(car)
    m = g.meta
    assert os.path.dirname(car) == os.path.join(str(tmp_path / 'datos'), 'GRABACIONES')
    assert not os.path.exists(os.path.join(car, 'cuadros'))              # no copia cuadros
    assert m['fuente']['archivos'] == [p.name for p in partes] and m['fuente']['paginas'] == [150, n - 150]
    assert m['fuente']['programa'].startswith('ThorCam')
    assert m['camara']['formato_pixel'] == 'bayer' and m['camara']['bayer_fase'] == 'RGGB'
    assert m['camara']['bits'] == 10 and m['camara']['roi'] == [0, 0, 64, 64]
    assert m['camara']['binning'] == [1, 1] and m['camara']['pixel_clock_hz'] == 4e7
    assert m['camara']['tick_hz'] == 1e9 and m['tiempo_fuente'] == 'camara'
    assert m['hora_inicio_iso'] == '2026-11-05T10:00:00.123'
    assert m['tipo'] == 'atrapada' and m['muestra']['glicerina_pct'] == 20.0
    # cuadros idénticos, leídos de la fuente (página 149 del _0 y página 0 del _1)
    assert len(g) == n
    for i in (0, 149, 150, n - 1):
        np.testing.assert_array_equal(g.cuadro(i), cuadros[i])
    assert g.tiempos['archivo'][150] == f'{NOMBRE_TC}_1.tif#0'
    # tiempos: contador y timestamp de los tags → todos los perdidos detectados
    inf = m['informe_tiempo']
    perd = g0.meta['verdad']['contadores_perdidos']
    assert len(perd) > 0 and inf['cuadros_perdidos'] == len(perd) and inf['contador'] == 'hardware'
    faltan = np.concatenate([np.arange(a, b + 1) for a, b in inf['posiciones_perdidos']])
    np.testing.assert_array_equal(faltan, perd)
    assert not any('timestamps' in f for f in m['importacion']['faltan'])

    # TIEMPO y los bloques: idénticos a los de la grabación original
    t_imp, t_org = _tiempo(car), _tiempo(g0.carpeta)
    np.testing.assert_allclose(t_imp['t'], t_org['t'], atol=1e-12)
    np.testing.assert_array_equal(t_imp['huecos'], t_org['huecos'])
    tr_imp, tr_org = _tracking(car), _tracking(g0.carpeta)
    for a, b in zip(tr_imp, tr_org):
        np.testing.assert_array_equal(a, b)
    assert np.mean([len(a) > 0 for a in tr_imp]) > 0.95
    g.cerrar()
    assert _hash([str(p) for p in partes]) == antes                      # originales intactos


def test_roi_impar_corre_la_fase_y_avisa_supuesta(tmp_path, monkeypatch):
    monkeypatch.setenv('PINZAS_DATOS', str(tmp_path / 'datos'))
    g0 = _sintetica(tmp_path, bayer='GRBG', n=20, perdido=0.0)
    src = tmp_path / 'a.tif'
    _escribir_tiff(str(src), [g0.cuadro(i) for i in range(20)],
                   [_tags_tsi(i + 1, 1000 * i, (101, 40, 64, 64)) for i in range(20)])
    car = importar.importar(str(src), dict(camara=dict(modelo='CS165CU')))
    m = io.leer_json(os.path.join(car, 'grabacion.json'))
    assert m['camara']['bayer_fase'] == importar.mover_fase('RGGB', 101, 40) == 'GRBG'
    assert m['camara']['bayer_fase_supuesta'] is True
    assert any('SUPUESTA' in f for f in m['importacion']['faltan'])
    assert any('exposicion_us' in f for f in m['importacion']['faltan'])
    assert m['camara']['exposicion_us'] is None                         # no se inventa


def test_secuencia_numerada_sin_tags_es_nominal(tmp_path, monkeypatch):
    monkeypatch.setenv('PINZAS_DATOS', str(tmp_path / 'datos'))
    src = tmp_path / 'seq'
    src.mkdir()
    rng = np.random.default_rng(0)
    cuadros = [rng.integers(0, 1023, (16, 20)).astype(np.uint16) for _ in range(12)]
    for k in rng.permutation(12):                                       # escritos en desorden
        tifffile.imwrite(str(src / f'cuadro_{k + 1}.tif'), cuadros[k], photometric='minisblack')
    car = importar.importar(str(src), dict(camara=dict(formato_pixel='mono', periodo_us=10000.0)))
    g = io.leer_grabacion(car)
    assert g.meta['fuente']['archivos'][:3] == ['cuadro_1.tif', 'cuadro_2.tif', 'cuadro_3.tif']
    np.testing.assert_array_equal(g.cuadro(9), cuadros[9])              # orden natural: 10 después de 9
    assert g.meta['tiempo_fuente'] == 'nominal'
    inf = g.meta['informe_tiempo']
    assert inf['metodo'] == 'nominal' and 'SIN timestamps' in inf['aviso']
    assert inf['cuadros_perdidos'] is None and inf['contador'] == 'ninguno'
    assert 'NO HAY timestamps' in inf['explicacion']
    assert any('timestamps' in f for f in g.meta['importacion']['faltan'])
    assert np.isnan(g.tiempos['contador_camara']).all() and np.isnan(g.tiempos['t_camara']).all()
    r = _tiempo(car)
    np.testing.assert_allclose(r['t'], np.arange(12) * 0.01)
    assert r['cuadros_perdidos'] == 0 and len(r['tramos']) == 1


def test_rgb_verde_es_el_plano_g_y_tracking_igual_al_mono(tmp_path, monkeypatch):
    monkeypatch.setenv('PINZAS_DATOS', str(tmp_path / 'datos'))
    g0 = _sintetica(tmp_path, bayer=None, n=120, perdido=0.0)
    n = len(g0)
    rgb = [np.stack([g0.cuadro(i) // 2, g0.cuadro(i), np.full((64, 64), 7, np.uint16)], axis=-1)
           for i in range(n)]
    src = tmp_path / 'rgb.tif'
    _escribir_tiff(str(src), rgb, [_tags_tsi(i + 1, None, (0, 0, 64, 64)) for i in range(n)],
                   photometric='rgb')
    car = importar.importar(str(src), dict(META, camara=dict(periodo_us=5000.0)))
    g = io.leer_grabacion(car)
    m = g.meta
    assert m['camara']['formato_pixel'] == 'rgb' and m['camara'].get('bayer_fase') is None
    assert any('RGB procesados' in a for a in m['importacion']['avisos'])
    c = g.cuadro(5)
    assert c.shape == (64, 64, 3) and c.dtype == np.uint16
    np.testing.assert_array_equal(imagen.preparar(c, m, 'verde'), g0.cuadro(5).astype(np.float32))
    np.testing.assert_array_equal(imagen.preparar(c, m, 'rojo'), (g0.cuadro(5) // 2).astype(np.float32))
    assert imagen.preparar(c, m, 'verde', recorte=(3, 5, 10, 11)).shape == (11, 10)
    assert imagen.px_sensor_por_px('verde', m) == 1.0
    assert imagen.a_coordenadas_sensor(2.5, 3.0, 'verde', None, (3, 5, 10, 11)) == (5.5, 8.0)
    assert imagen.vista_rapida(c, m).shape == (64, 64, 3)
    # contador sin timestamp: perdidos detectables, t nominal
    assert m['tiempo_fuente'] == 'nominal' and m['informe_tiempo']['contador'] == 'hardware'
    tr, tr0 = _tracking(car, 6.0), _tracking(g0.carpeta, 6.0)        # RGB: el verde es de resolución completa
    for a, b in zip(tr, tr0):
        np.testing.assert_array_equal(a, b)
    assert sum(len(a) for a in tr) > 0


def test_ome_tiff_delta_t_y_exposicion(tmp_path, monkeypatch):
    monkeypatch.setenv('PINZAS_DATOS', str(tmp_path / 'datos'))
    d = np.arange(5 * 8 * 8, dtype=np.uint16).reshape(5, 8, 8)
    dt = [0.0, 0.01, 0.02, 0.04, 0.05]
    tifffile.imwrite(str(tmp_path / 'x.ome.tif'), d, ome=True, photometric='minisblack',
                     metadata=dict(axes='TYX', PhysicalSizeX=3.45,
                                   Plane=dict(DeltaT=dt, ExposureTime=[0.002] * 5)))
    car = importar.importar(str(tmp_path / 'x.ome.tif'), dict(camara=dict(formato_pixel='mono')))
    g = io.leer_grabacion(car)
    assert g.meta['tiempo_fuente'] == 'ome_deltat' and g.meta['camara']['exposicion_us'] == 2000.0
    assert g.meta['fuente']['ome_tamano_px_m'] == pytest.approx(3.45e-6)
    assert g.meta.get('escala_m_por_px') is None                         # no se usa como escala
    r = tiempo.reconstruir(g.tiempos, tick_hz=g.meta['camara']['tick_hz'])
    np.testing.assert_allclose(r['t'], dt, atol=1e-9)
    assert r['informe']['contador'] == 'ninguno'
    np.testing.assert_array_equal(g.cuadro(3), d[3])


def test_avi_8_bits_avisa_y_lee_al_azar(tmp_path, monkeypatch):
    import cv2
    monkeypatch.setenv('PINZAS_DATOS', str(tmp_path / 'datos'))
    ruta = str(tmp_path / 'grab.avi')
    w = cv2.VideoWriter(ruta, cv2.VideoWriter_fourcc(*'MJPG'), 25.0, (48, 32))
    assert w.isOpened()
    yy, xx = np.mgrid[0:32, 0:48]
    for k in range(20):
        f = np.stack([(xx * 5 + 9 * k) % 256, (yy * 7 + 9 * k) % 256, np.full_like(xx, 40)], -1).astype(np.uint8)
        w.write(f)
    w.release()
    car = importar.importar(ruta, dict(tipo='libre'))
    g = io.leer_grabacion(car)
    m = g.meta
    assert len(g) == 20 and m['fuente']['formato'] == 'video' and m['camara']['bits'] == 8
    assert m['camara']['formato_pixel'] == 'rgb'
    assert any('subpíxel' in a for a in m['importacion']['avisos'])
    assert m['tiempo_fuente'] == 'nominal' and m['informe_tiempo']['periodo_s'] == pytest.approx(0.04)
    c15, c3 = g.cuadro(15), g.cuadro(3)
    assert c3.shape == (32, 48, 3) and c3.max() <= 255
    sec = io.LectorFuente(dict(m['fuente']))
    for k in range(4):
        ref = sec.leer('grab.avi', k)
    sec.cerrar()
    np.testing.assert_array_equal(c3, ref)
    g.cerrar()


def test_reconstruir_nominal_y_hardware_sin_host_sin_warnings():
    n = 50
    nan = np.full(n, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        r = tiempo.reconstruir(dict(contador_camara=nan, t_camara=nan, t_host_ns=nan,
                                    t_nominal_s=np.arange(n) * 0.005))
        assert r['informe']['metodo'] == 'nominal' and np.isnan(r['informe']['cuadros_perdidos'])
        np.testing.assert_allclose(r['t'], np.arange(n) * 0.005)
        c = np.r_[np.arange(1, 21), np.arange(23, 53)]                  # faltan 21 y 22
        tc = (c - 1) * 5_000_000.0
        r = tiempo.reconstruir(dict(contador_camara=c, t_camara=tc, t_host_ns=nan), tick_hz=1e9)
        assert r['informe']['metodo'] == 'hardware' and r['cuadros_perdidos'] == 2
        assert r['informe']['posiciones_perdidos'] == [[21, 22]]
        assert np.isnan(r['informe']['ppm_reloj'])
        assert r['informe']['fps_medido'] == pytest.approx(200.0)
        r = tiempo.reconstruir(dict(contador_camara=nan, t_camara=nan, t_host_ns=nan))
        assert r['informe']['metodo'] == 'sin_tiempos' and np.isnan(r['t']).all()


def test_auxiliares():
    assert importar.mover_fase('RGGB', 1, 0) == 'GRBG'
    assert importar.mover_fase('RGGB', 0, 1) == 'GBRG'
    assert importar.mover_fase('RGGB', 1, 1) == 'BGGR'
    assert importar.hora_de_nombre('image_2022-02-17T14-55-13.447_0.tif') == '2022-02-17T14:55:13.447'
    assert importar.hora_de_nombre('cuadro_1.tif') is None
    conf = []
    m = importar.fusionar(dict(a=1, b=None, c=dict(x=2)), dict(a=3, b=4, c=dict(x=2, y=5)), '', conf)
    assert m == dict(a=1, b=4, c=dict(x=2, y=5)) and len(conf) == 1 and conf[0].startswith('a:')
    mos = sintetico.mosaico_bayer(np.stack([np.full((8, 8), v) for v in (300.0, 600.0, 200.0)], -1), 'RGGB')
    assert importar.diagonal_verde(mos)[1] == '01-10'
    assert importar.diagonal_verde(np.full((8, 8), 5.0))[1] == 'parece_mono'
    assert sorted(['a_10.tif', 'a_9.tif', 'a_1.tif'], key=importar.orden_natural) == ['a_1.tif', 'a_9.tif', 'a_10.tif']


def test_serie_thorcam_primera_parte_sin_numero(tmp_path):
    """ThorCam nombra las partes 'x.tif', 'x_1.tif', ..., 'x_10.tif': es UNA serie, en ese orden."""
    for n in ('50G-50W.tif', '50G-50W_1.tif', '50G-50W_2.tif', '50G-50W_10.tif', 'otra.tif'):
        (tmp_path / n).write_bytes(b'')
    esperado = ['50G-50W.tif', '50G-50W_1.tif', '50G-50W_2.tif', '50G-50W_10.tif']
    for desde in ('50G-50W.tif', '50G-50W_2.tif'):
        assert importar.listar_archivos(str(tmp_path / desde))[2] == esperado
    assert importar.listar_archivos(str(tmp_path / '50G-50W.tif'), serie=False)[2] == ['50G-50W.tif']


def test_serie_thorcam_nombre_que_termina_en_numero(tmp_path):
    """'BROWNIANO1.tif', 'BROWNIANO1_1.tif', ...: el 1 es parte del nombre, y es UNA serie."""
    esperado = ['BROWNIANO1.tif', 'BROWNIANO1_1.tif', 'BROWNIANO1_2.tif', 'BROWNIANO1_10.tif']
    for n in esperado:
        (tmp_path / n).write_bytes(b'')
    for desde in ('BROWNIANO1.tif', 'BROWNIANO1_2.tif', ''):
        assert importar.listar_archivos(str(tmp_path / desde))[2] == esperado
    sec = tmp_path / 'sec'
    sec.mkdir()
    for n in ('cuadro001.tif', 'cuadro002.tif', 'cuadro010.tif'):          # secuencia sin '_': sigue siendo una
        (sec / n).write_bytes(b'')
    assert importar.listar_archivos(str(sec))[2] == ['cuadro001.tif', 'cuadro002.tif', 'cuadro010.tif']


def test_script_de_linea_de_comando(tmp_path, monkeypatch, capsys):
    import importlib.util
    monkeypatch.setenv('PINZAS_DATOS', str(tmp_path / 'datos'))
    src = tmp_path / f'{NOMBRE_TC}_0.tif'
    _escribir_tiff(str(src), [np.full((8, 8), i, np.uint16) for i in range(6)],
                   [_tags_tsi(i + 1, 10_000_000 * i, (0, 0, 8, 8)) for i in range(6)])
    ruta = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        'TouchDesigner', 'herramientas', 'importar_thorcam.py')
    spec = importlib.util.spec_from_file_location('importar_thorcam', ruta)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    car = mod.main([str(src), '--tipo', 'atrapada', '--diametro_um', '3', '--glicerina_pct', '20',
                    '--exposicion_us', '2000', '--fps', '100', '--escala_nm_px', '90',
                    '--temperatura_C', '23'])
    out = capsys.readouterr().out
    assert 'Creada:' in out and 'FALTA' in out and 'camara' in out
    m = io.leer_json(os.path.join(car, 'grabacion.json'))
    assert m['muestra']['diametro_m'] == pytest.approx(3e-6)
    assert m['escala_m_por_px'] == pytest.approx(9e-8)
    assert m['condiciones']['temperatura_K'] == pytest.approx(296.15)
    assert m['camara']['exposicion_us'] == 2000.0 and m['camara']['periodo_us'] == pytest.approx(1e4)
    assert m['informe_tiempo']['fps_medido'] == pytest.approx(100.0)


def test_video_descarta_cuadros_repetidos(tmp_path, monkeypatch):
    """Video con fps variable / de pantalla: los cuadros que son copia del anterior se descartan y
    cada cuadro nuevo conserva el tiempo en que aparece (índice original / fps)."""
    import cv2
    monkeypatch.setenv('PINZAS_DATOS', str(tmp_path / 'datos'))
    ruta = str(tmp_path / 'rep.avi')
    w = cv2.VideoWriter(ruta, cv2.VideoWriter_fourcc(*'MJPG'), 10.0, (48, 32))
    yy, xx = np.mgrid[0:32, 0:48]
    contenido = [0, 1, 1, 2, 3, 3, 3, 4, 5, 6]            # el cuadro 2 repite al 1, el 5 y el 6 al 4...
    for k in contenido:
        w.write(np.stack([(xx * 5 + 40 * k) % 256, (yy * 7 + 40 * k) % 256, np.full_like(xx, 40)], -1).astype(np.uint8))
    w.release()
    car = importar.importar(ruta, dict(tipo='libre'))
    g = io.leer_grabacion(car)
    assert len(g) == 7 and g.meta['fuente']['cuadros_repetidos'] == 3
    assert any('REPETIDOS' in a for a in g.meta['importacion']['avisos'])
    np.testing.assert_allclose(_tiempo(car)['t'], np.array([0, 1, 3, 4, 7, 8, 9]) * 0.1, atol=1e-9)
    g.cerrar()
