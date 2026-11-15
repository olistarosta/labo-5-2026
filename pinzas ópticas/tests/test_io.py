"""io: formato de disco de CONTRATO.md, ida y vuelta, y que nada se pise."""
import json
import os

import numpy as np
import pytest
import tifffile

from pinzas import io, sintetico


def _beads():
    return [dict(x0_px=20.3, y0_px=17.8, D_m2s=4e-13, radio_px=5.0),
            dict(x0_px=44.0, y0_px=40.0, D_m2s=1e-13, kappa_x=2e-6, kappa_y=5e-6, radio_px=6.0)]


def test_grabar_y_leer(tmp_path):
    car = sintetico.grabar(str(tmp_path), 'prueba io', 'atrapada', n=80, fps=100.0, W=2e-3,
                           beads=_beads(), forma=(64, 64), prob_perdido=0.03, seed=1)
    assert os.path.basename(car).endswith('_prueba_io')
    g = io.leer_grabacion(car)
    ng = len(g.meta['verdad']['t_inicio_s'])
    assert len(g) == ng == len(os.listdir(os.path.join(car, 'cuadros')))
    assert g.nombre == g.meta['nombre'] == os.path.basename(car)
    # columnas y tipos
    assert list(g.tiempos) == io.COLUMNAS_TIEMPOS
    assert g.tiempos['contador_camara'].dtype == np.int64
    assert g.tiempos['t_host_ns'].dtype == np.int64
    assert g.tiempos['t_camara'].dtype == np.float64
    assert g.tiempos['archivo'][0] == '000000.tif'
    # los perdidos no tienen fila: el contador salta
    saltos = np.diff(g.tiempos['contador_camara'])
    assert np.sum(saltos - 1) == len(g.meta['verdad']['contadores_perdidos'])
    # cuadro crudo uint16 de 10 bits, de solo lectura y cacheado
    c = g.cuadro(3)
    assert c.dtype == np.uint16 and c.shape == (64, 64) and c.max() < 1024
    assert not c.flags.writeable
    assert g.cuadro(3) is c
    np.testing.assert_array_equal(g.cuadro(-1), tifffile.imread(os.path.join(car, 'cuadros', f'{ng - 1:06d}.tif')))
    with pytest.raises(IndexError):
        g.cuadro(ng)
    # verdad: t_camara en ticks del reloj REAL, tick_hz NOMINAL en el json
    assert g.meta['camara']['tick_hz'] == 1e9
    assert g.meta['verdad']['tick_hz_real'] != 1e9
    assert np.shape(g.meta['verdad']['x_px']) == (ng, 2)


def test_cache_lru_chico(tmp_path):
    car = sintetico.grabar(str(tmp_path), 'cache', 'libre', n=20, fps=50.0, W=1e-3, beads=[],
                           forma=(8, 8), prob_perdido=0.0, seed=0)
    g = io.leer_grabacion(car, cache=4)
    for i in range(10):
        g.cuadro(i)
    assert list(g._cache) == [6, 7, 8, 9]


def test_t_camara_vacio_es_nan(tmp_path):
    cols = dict(archivo=['000000.tif', '000001.tif'], contador_camara=[1, 2],
                t_camara=[np.nan, np.nan], t_host_ns=[10, 20], exposicion_us=[100.0, 100.0])
    io.escribir_tiempos_csv(str(tmp_path), cols)
    texto = open(tmp_path / 'tiempos.csv', encoding='utf-8').read().splitlines()
    assert texto[0] == 'archivo,contador_camara,t_camara,t_host_ns,exposicion_us'
    assert texto[1] == '000000.tif,1,,10,100'
    t = io.leer_csv(str(tmp_path / 'tiempos.csv'))
    assert np.isnan(t['t_camara']).all()
    assert t['t_host_ns'].dtype == np.int64


def test_ticks_grandes_exactos(tmp_path):
    tc = np.array([4_321_000_000_000, 4_321_010_000_137, 9_000_000_000_000_001 // 1000])
    io.escribir_tiempos_csv(str(tmp_path), dict(archivo=['a', 'b', 'c'], contador_camara=[1, 2, 3],
                                                t_camara=tc, t_host_ns=[1, 2, 3], exposicion_us=[1, 1, 1]))
    t = io.leer_csv(str(tmp_path / 'tiempos.csv'))
    np.testing.assert_array_equal(t['t_camara'].astype(np.int64), tc)


def test_bloques_multipagina(tmp_path):
    car = tmp_path / 'grab'
    (car / 'cuadros').mkdir(parents=True)
    pilas = [np.arange(5 * 4 * 3, dtype=np.uint16).reshape(5, 4, 3) + 100 * b for b in range(2)]
    for b, p in enumerate(pilas):
        tifffile.imwrite(str(car / 'cuadros' / f'bloque_{b:06d}.tif'), p, photometric='minisblack')
    # sin columna 'archivo' útil: página global i
    io.escribir_tiempos_csv(str(car), dict(archivo=[''] * 10, contador_camara=np.arange(10),
                                           t_camara=np.arange(10), t_host_ns=np.arange(10),
                                           exposicion_us=np.ones(10)))
    g = io.leer_grabacion(str(car))
    assert len(g) == 10
    np.testing.assert_array_equal(g.cuadro(7), pilas[1][2])
    np.testing.assert_array_equal(g.cuadro(0), pilas[0][0])
    # con 'archivo' = bloque#página
    arch = [f'bloque_{i // 5:06d}.tif#{i % 5}' for i in range(10)][::-1]
    io.escribir_tiempos_csv(str(car), dict(archivo=arch, contador_camara=np.arange(10),
                                           t_camara=np.arange(10), t_host_ns=np.arange(10),
                                           exposicion_us=np.ones(10)))
    g = io.leer_grabacion(str(car))
    np.testing.assert_array_equal(g.cuadro(0), pilas[1][4])


def test_nueva_carpeta_grabacion_no_pisa(tmp_path):
    a = io.nueva_carpeta_grabacion('x', raiz=str(tmp_path))
    b = io.nueva_carpeta_grabacion('x', raiz=str(tmp_path))
    c = io.nueva_carpeta_grabacion('x', raiz=str(tmp_path))
    assert len({a, b, c}) == 3
    assert all(os.path.isdir(os.path.join(r, 'cuadros')) for r in (a, b, c))
    nombre = os.path.basename(a)
    assert len(nombre.split('_')[0]) == 10 and len(nombre.split('_')[1]) == 6


def test_analisis_no_pisa_y_json_nativo(tmp_path):
    car = io.nueva_carpeta_grabacion('g', raiz=str(tmp_path))
    a1 = io.nueva_carpeta_analisis(car)
    a2 = io.nueva_carpeta_analisis(io.leer_grabacion(car))
    assert a1 != a2 and os.path.dirname(a1) == os.path.join(car, 'analisis')
    assert os.path.isdir(os.path.join(a1, 'figuras'))
    par = dict(método='MSD', ventana=np.int64(15), escala=np.float32(9e-8), umbral=np.nan,
               lista=np.arange(3), anidado={'κ': (np.float64(1.5), np.bool_(True))})
    res = [dict(bead=0, D=np.float64(4.1e-13), err=np.nan), dict(bead=1, D=2e-13, extra='ñ')]
    rp, rr = io.guardar_analisis(a1, par, res)
    texto = open(rp, encoding='utf-8').read()
    assert 'método' in texto and 'κ' in texto                # ensure_ascii=False
    p = json.loads(texto)
    assert p['ventana'] == 15 and p['umbral'] is None and p['lista'] == [0, 1, 2]
    assert p['anidado']['κ'] == [1.5, True]
    r = io.leer_csv(rr)
    assert list(r) == ['bead', 'D', 'err', 'extra']
    assert r['D'][0] == 4.1e-13 and np.isnan(r['err']).all()
    with pytest.raises(FileExistsError):
        io.guardar_analisis(a1, par, res)
    # dict de columnas y dict de escalares
    io.guardar_analisis(a2, {}, dict(t=np.linspace(0, 1, 5), x=np.arange(5.0)))
    assert len(io.leer_csv(os.path.join(a2, 'resultados.csv'))['t']) == 5
    a3 = io.nueva_carpeta_analisis(car)
    io.guardar_analisis(a3, {}, dict(D=1.0, err_D=0.1))
    assert io.leer_csv(os.path.join(a3, 'resultados.csv'))['err_D'][0] == 0.1


def _video_mjpg(ruta, n=12, forma=(48, 64)):
    """AVI Motion-JPEG chico: cada cuadro tiene su número grabado como nivel de gris."""
    import cv2
    w = cv2.VideoWriter(io.ruta_ascii(os.path.dirname(ruta)) + os.sep + os.path.basename(ruta),
                        cv2.VideoWriter_fourcc(*'MJPG'), 10, forma[::-1])
    assert w.isOpened()
    for i in range(n):
        w.write(np.full(forma + (3,), 10 + 18 * i, np.uint8))
    w.release()


def test_video_en_carpeta_con_acentos(tmp_path):
    """cv2 no abre rutas con acentos en Windows ('pinzas ópticas'): abrir_video tiene que poder."""
    carpeta = tmp_path / 'pinzas ópticas'
    carpeta.mkdir()
    ruta = str(carpeta / 'v.avi')
    _video_mjpg(ruta)
    assert io.ruta_ascii(ruta).isascii()
    cap, _ = io.abrir_video(ruta)
    ok, f = cap.read()
    assert ok and abs(int(f[..., 1].mean()) - 10) <= 3


def test_video_lectura_exacta_en_desorden(tmp_path):
    """El lector devuelve el cuadro pedido aunque se pida en desorden (en MP4 la búsqueda de
    cv2 puede caer en el de al lado; acá se fuerza el camino de leer en orden con fourcc != MJPG)."""
    ruta = str(tmp_path / 'v.avi')
    _video_mjpg(ruta)
    for fourcc in ('MJPG', 'H264'):
        lec = io.LectorFuente(dict(formato='video', carpeta=str(tmp_path), archivos=['v.avi'],
                                   formato_pixel='rgb', fourcc=fourcc))
        for k in (5, 2, 11, 0, 7, 8, 3):
            assert abs(int(lec.leer('v.avi', k)[..., 1].mean()) - (10 + 18 * k)) <= 3, (fourcc, k)
