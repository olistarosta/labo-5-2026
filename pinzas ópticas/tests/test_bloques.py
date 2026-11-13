"""Los bloques de TouchDesigner/bloques/ (los containers de TD) corridos sin TD, como hace el lote."""
import csv
import os
import sys

import numpy as np

from pinzas import fisica, io, sintetico, tracking

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BLOQUES = os.path.join(RAIZ, 'TouchDesigner', 'bloques')
sys.path.insert(0, BLOQUES)
import _lote  # noqa: E402

ARCHIVOS = ['b0_cuadro', 'b1_canal', 'b2_fondo', 'b3_fourier', 'b4_clahe', 'b5_realce', 'b6_deteccion',
            'b7_localizacion', 'b8_trayectorias']


# Los defaults de los bloques están calibrados para los MP4 (beads de ~17 px, 8 bits). Los
# sintéticos son grabaciones propias con beads de ~4 px: estos son sus valores.
SINTETICOS = {'b0_cuadro': {'Formato': 'grabacion'}, 'b3_fourier': {'Chico': 1.0, 'Grande': 20.0},
              'b6_deteccion': {'Umbral': 6.0},
              'b8_trayectorias': {'Maxsalto': 5.0, 'Memoria': 3, 'Largomin': 20}}


def _defectos(nombre):
    m = _lote.cargar(os.path.join(BLOQUES, nombre + '.py'))
    p = {n: v for n, tipo, v, *_ in m.PARAMETROS if tipo != 'Pulse'}
    p.update(SINTETICOS.get(nombre, {}))
    return p


def test_simetria_radial_igual_a_la_validada():
    m = _lote.cargar(os.path.join(BLOQUES, 'b7_localizacion.py'))
    rng = np.random.default_rng(0)
    for _ in range(20):
        x, y = 10 + rng.random(2)
        rec = sintetico.camara(sintetico.dibujar((21, 21), [[x, y]], 4.0), rng=rng).astype(float)
        for s in (0.0, 0.8):
            assert np.allclose(m.simetria_radial(rec, s), tracking.simetria_radial(rec, s), atol=1e-12)


def test_lote_libre_da_csv_con_las_beads(tmp_path, monkeypatch):
    monkeypatch.setenv('PINZAS_DATOS', str(tmp_path / 'datos'))
    MPX, R = 9e-8, 1.5e-6
    rad = 0.5 * R / MPX * 0.5
    D = fisica.difusion(296.15, fisica.gamma_stokes(1.7e-3, R))
    rng = np.random.default_rng(3)
    c = sintetico.grabar(str(tmp_path), 'libre', 'libre', 300, 30.0, 5e-3,
                         [dict(x0_px=float(a), y0_px=float(b), D_m2s=D, radio_px=rad) for a, b in rng.uniform(40, 150, (4, 2))],
                         m_por_px=MPX, forma=(192, 192), seed=4)
    meta = io.leer_json(os.path.join(c, 'grabacion.json'))
    meta['escala_m_por_px'] = MPX
    io.escribir_grabacion_json(c, meta)
    bloques = [dict(nombre=n, archivo=os.path.join(BLOQUES, n + '.py'), p=_defectos(n)) for n in ARCHIVOS]
    bloques[0]['p']['Carpeta'] = c
    bloques[5]['p']['Radio'] = rad
    carpeta, npart, nfil = _lote.procesar(bloques)
    assert npart == 4
    with open(os.path.join(carpeta, 'trayectorias.csv'), encoding='utf-8') as f:
        filas = list(csv.DictReader(f))
    assert len(filas) > 0.95 * 4 * len(io.leer_grabacion(c))
    for f in ('particulas.csv', 'tiempo.csv', 'parametros.json', 'LEEME.txt', 'grabacion.json'):
        assert os.path.exists(os.path.join(carpeta, f))
