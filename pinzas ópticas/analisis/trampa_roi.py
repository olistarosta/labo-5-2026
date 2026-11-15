"""Detecciones alrededor de la trampa en TODOS los cuadros, sin enlazar: lo que hace falta para seguir a
la partícula atrapada durante una patada, cuando se aparta varios px por cuadro y el seguimiento
general (trayectorias.csv, salto máximo 10 px) la pierde.

    python trampa_roi.py "BARRIDO 1" ["LASER NE" ...]        (sin argumentos: todas las del catálogo)
    python trampa_roi.py "BARRIDO 1" --tapar=11              (control: con el núcleo del láser rellenado)

Lee directo los TIFF (sin compresión) de C:\\datos_pinzas\\GRABACIONES: solo una ventana de
±VENTANA px alrededor del láser. Usa los MISMOS pasos que la conversión general, con el código de
los bloques (TouchDesigner/bloques): canal verde, fondo por mediana temporal con el disco de la
trampa rellenado, pasabanda de Fourier, realce, detección (umbral más bajo) y localización por
simetría radial.

Escribe datos/<grabación>/trampa_detecciones.csv: una fila por detección
    cuadro
    x_px, y_px   posición (px del sensor, mismo sistema que trayectorias.csv)
    z            contraste (realce / ruido)
"""
import importlib.util
import os
import re
import sys

import numpy as np
import pandas as pd
from PIL import Image

from pinzas_datos import CATALOGO, DATOS, RAIZ

GRABACIONES = r"C:\datos_pinzas\GRABACIONES"
VENTANA = 128                # medio lado de la ventana (px)
N_FONDO = 61                 # cuadros de la mediana temporal
RADIO_RELLENO = 45.0         # disco de la trampa que se rellena en el fondo (px), como en la conversión
UMBRAL = 8.0                 # la conversión general usa 12; acá interesa no perder a la partícula borrosa
P_REALCE, P_LOCAL = {'Radio': 17.0}, {'Lado': 0, 'Suavizado': 0.2}
BLOQUES = os.path.join(RAIZ, 'TouchDesigner', 'bloques')


def _bloque(nombre):
    spec = importlib.util.spec_from_file_location(nombre, os.path.join(BLOQUES, nombre + '.py'))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


b2, b3, b5, b6, b7 = (_bloque(n) for n in ('b2_fondo', 'b3_fourier', 'b5_realce', 'b6_deteccion', 'b7_localizacion'))


def orden(n):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', os.path.splitext(n)[0])]


def paginas(nombre):
    """[(ruta, desplazamiento en bytes de los píxeles, ancho)] de cada cuadro, en orden."""
    carpeta = os.path.join(GRABACIONES, nombre)
    out = []
    for f in sorted([f for f in os.listdir(carpeta) if f.lower().endswith('.tif')], key=orden):
        ruta = os.path.join(carpeta, f)
        with Image.open(ruta) as im:
            ancho = im.size[0]
            for k in range(im.n_frames):
                im.seek(k)
                out.append((ruta, im.tile[0][2], ancho))
    return out


def ventana_verde(pag, i0, j0):
    ruta, off, ancho = pag
    lado = 2 * VENTANA + 1
    a = np.memmap(ruta, dtype=np.uint8, mode='r', offset=off + i0 * ancho * 3, shape=(lado, ancho, 3))
    return np.asarray(a[:, j0:j0 + lado, 1], np.float32)


def detectar(img, fondo, nucleo=None):
    """Los pasos B2 (dividir) → B3 → B5 → B6 → B7 sobre la ventana. Devuelve [k, 3]: x, y (ventana), z.
    nucleo = (x, y, radio) en px de la ventana: ese disco (el láser) se rellena con una interpolación
    suave antes de detectar, para que el punto oscuro que deja el láser no corra la posición."""
    plano = np.where(fondo > 0, img / np.where(fondo > 0, fondo, 1), 1.0).astype(np.float32)
    if nucleo is not None:
        plano = b2.rellenar_disco(plano, *nucleo)
    trabajo = b3.pasabanda(plano, 1.5, 60.0)
    d = {'img': trabajo, 'trabajo': trabajo}
    d = b5.aplicar(d, P_REALCE, None)
    d = b6.aplicar(d, {'Umbral': UMBRAL, 'Maximo': 20}, None)
    d = b7.aplicar(d, P_LOCAL, None)
    return d['posiciones']


def medir(nombre, tapar_nucleo=0.0):
    """tapar_nucleo > 0: radio (px) del disco del láser que se rellena en cada cuadro; la salida va a
    trampa_detecciones_tapado.csv. Sirve para medir cuánto corre la posición el punto oscuro del láser
    (la cámara le resta rojo al verde: ver PINZAS 4)."""
    tx, ty = CATALOGO[nombre]['trampa_px']
    j0, i0 = int(round(tx)) - VENTANA, int(round(ty)) - VENTANA
    pags = paginas(nombre)
    idx = np.unique(np.linspace(0, len(pags) - 1, N_FONDO).round().astype(int))
    fondo = np.median(np.stack([ventana_verde(pags[i], i0, j0) for i in idx]), axis=0).astype(np.float32)
    fondo = b2.rellenar_disco(fondo, tx - j0, ty - i0, RADIO_RELLENO)
    nucleo = (tx - j0, ty - i0, float(tapar_nucleo)) if tapar_nucleo > 0 else None
    filas = []
    for c, pag in enumerate(pags):
        for x, y, z in detectar(ventana_verde(pag, i0, j0), fondo, nucleo):
            filas.append((c, x + j0, y + i0, z))
        if c % 1000 == 0 and c:
            print('%s: %d de %d' % (nombre, c, len(pags)), flush=True)
    d = pd.DataFrame(filas, columns=['cuadro', 'x_px', 'y_px', 'z'])
    salida = 'trampa_detecciones_tapado.csv' if nucleo else 'trampa_detecciones.csv'
    d.to_csv(os.path.join(DATOS, nombre, salida), index=False, float_format='%.6g')
    print('%s: %d detecciones en %d cuadros (%.2f por cuadro)' % (nombre, len(d), len(pags), len(d) / len(pags)), flush=True)


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    tapar = next((float(a.split('=')[1]) for a in sys.argv[1:] if a.startswith('--tapar=')), 0.0)
    for n in (args or list(CATALOGO)):
        medir(n, tapar)
