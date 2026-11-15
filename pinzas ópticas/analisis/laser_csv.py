"""¿Hay una partícula en la trampa? La partícula atrapada dispersa el láser y se ve una mancha roja.
Mide en cada cuadro cuánto rojo de más hay en la posición del láser y escribe datos/<grabación>/laser.csv.

    python laser_csv.py "BARRIDO 1" ["LASER NE" ...]        (sin argumentos: todas las del catálogo)

Lee directo los TIFF de C:\\datos_pinzas\\GRABACIONES (sin compresión: solo las filas que hacen falta).
Columnas de laser.csv
    cuadro
    rojo        cociente rojo/verde promediado en un disco de 6 px en el láser, menos la mediana del
                mismo cociente en un anillo de 25 a 40 px (≈ 0 sin partícula; > 0 con partícula)
    rojo_x_px, rojo_y_px   centroide de la mancha roja (px del sensor); vacío si no hay mancha
"""
import os
import re
import sys

import numpy as np
import pandas as pd
from PIL import Image

from pinzas_datos import CATALOGO, DATOS

GRABACIONES = r"C:\datos_pinzas\GRABACIONES"
M = 40                                                     # medio lado de la ventana (px)


def orden(n):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', os.path.splitext(n)[0])]


def medir(nombre):
    tx, ty = CATALOGO[nombre]['trampa_px']
    j0, i0 = int(round(tx)) - M, int(round(ty)) - M
    yy, xx = np.mgrid[i0:i0 + 2 * M + 1, j0:j0 + 2 * M + 1]
    r = np.hypot(xx - tx, yy - ty)
    disco, anillo = r <= 6, (r >= 25) & (r <= 40)
    carpeta = os.path.join(GRABACIONES, nombre)
    filas = []
    for f in sorted([f for f in os.listdir(carpeta) if f.lower().endswith('.tif')], key=orden):
        ruta = os.path.join(carpeta, f)
        with Image.open(ruta) as im:
            ancho, alto = im.size
            desplazamientos = []
            for k in range(im.n_frames):
                im.seek(k)
                desplazamientos.append(im.tile[0][2])       # dónde empiezan los píxeles de la página k
        for off in desplazamientos:
            a = np.memmap(ruta, dtype=np.uint8, mode='r', offset=off + i0 * ancho * 3, shape=(2 * M + 1, ancho, 3))
            v = np.asarray(a[:, j0:j0 + 2 * M + 1, :], np.float32)
            q = v[:, :, 0] / (v[:, :, 1] + 1.0)
            exceso = q - np.median(q[anillo])
            rojo = float(exceso[disco].mean())
            w = np.clip(exceso, 0, None) * (r <= 15)
            if rojo > 0.05 and w.sum() > 0:
                cx, cy = float((w * xx).sum() / w.sum()), float((w * yy).sum() / w.sum())
            else:
                cx = cy = np.nan
            filas.append((rojo, cx, cy))
    d = pd.DataFrame(filas, columns=['rojo', 'rojo_x_px', 'rojo_y_px'])
    d.insert(0, 'cuadro', np.arange(len(d)))
    os.makedirs(os.path.join(DATOS, nombre), exist_ok=True)
    d.to_csv(os.path.join(DATOS, nombre, 'laser.csv'), index=False, float_format='%.4g')
    print('%s: %d cuadros | rojo: mediana %.3f, p5 %.3f, p95 %.3f' % (nombre, len(d), d.rojo.median(), d.rojo.quantile(.05), d.rojo.quantile(.95)), flush=True)


if __name__ == '__main__':
    for n in (sys.argv[1:] or list(CATALOGO)):
        medir(n)
