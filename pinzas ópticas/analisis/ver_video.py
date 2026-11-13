"""Arma un MP4 para MIRAR una grabación de ThorCam (TIFF multipágina partido en varios archivos).

Uso:
    python ver_video.py "C:\\datos_pinzas\\GRABACIONES\\BROWNIANO 1"
    python ver_video.py "C:\\datos_pinzas\\GRABACIONES\\BROWNIANO 1" --fps 10 --salida D:\\otro.mp4

Le pasás la carpeta de la grabación (junta todos los .tif en orden numérico: X.tif, X_1.tif,
X_2.tif, ...) o un solo .tif (convierte solo ese). El video queda en
C:\\datos_pinzas\\VISTAS\\<nombre>.mp4, a resolución completa, en H.264 (lo abre cualquier
reproductor de Windows).

fps: los TIFF de ThorCam no lo traen. Si no lo pasás, se estima con las fechas de modificación
de los archivos (cuadros de cada archivo / segundos desde que se cerró el anterior, la mediana);
si no se puede (un solo archivo), 10. Es una estimación: sirve para mirar, no para medir.

El MP4 es solo para mirar (está comprimido): el análisis se hace siempre sobre los TIFF, que
este script solo lee. Necesita numpy, Pillow y ffmpeg en el PATH.
"""
import argparse
import os
import re
import subprocess
import sys

import numpy as np
from PIL import Image


def orden_numerico(nombre):
    """'X.tif' < 'X_1.tif' < 'X_2.tif' < 'X_10.tif' (el orden alfabético pondría _10 antes que _2)."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', os.path.splitext(nombre)[0])]


def estimar_fps(archivos, paginas):
    """Cada archivo se cierra cuando se escribió su último cuadro: entre la fecha de un archivo y
    la del siguiente se escribieron los cuadros del siguiente. Mediana de (cuadros / segundos)
    entre archivos vecinos, para que una pausa del disco no la corra."""
    dt = np.diff([os.path.getmtime(f) for f in archivos])
    tasas = [n / d for n, d in zip(paginas[1:], dt) if d > 0]
    return float(np.median(tasas)) if tasas else None


def main():
    ap = argparse.ArgumentParser(description='TIFF de ThorCam -> MP4 para mirar.')
    ap.add_argument('ruta', help='carpeta de la grabación, o un .tif')
    ap.add_argument('--fps', type=float, help='cuadros por segundo (si falta, se estima)')
    ap.add_argument('--salida', help='ruta del .mp4 (si falta, C:/datos_pinzas/VISTAS/<nombre>.mp4)')
    a = ap.parse_args()

    ruta = os.path.abspath(a.ruta)
    if os.path.isdir(ruta):
        nombres = sorted((f for f in os.listdir(ruta) if f.lower().endswith(('.tif', '.tiff'))),
                         key=orden_numerico)
        archivos = [os.path.join(ruta, f) for f in nombres]
        nombre = os.path.basename(ruta)
    else:
        archivos = [ruta]
        nombre = os.path.splitext(os.path.basename(ruta))[0]
    if not archivos:
        sys.exit(f'no hay ningún .tif en {ruta}')

    salida = a.salida or os.path.join(os.environ.get('PINZAS_DATOS', 'C:/datos_pinzas'),
                                      'VISTAS', nombre + '.mp4')
    salida = os.path.abspath(salida)
    os.makedirs(os.path.dirname(salida), exist_ok=True)

    paginas = []
    for f in archivos:
        with Image.open(f) as im:
            paginas.append(im.n_frames)
    total = sum(paginas)

    if a.fps:
        fps, origen = a.fps, 'el que pasaste'
    else:
        fps = estimar_fps(archivos, paginas)
        origen = 'estimado de las fechas de los archivos' if fps else 'no se pudo estimar: 10 por defecto'
        fps = round(fps, 2) if fps else 10.0

    with Image.open(archivos[0]) as im:
        primero = np.asarray(im)
    alto, ancho = primero.shape[:2]
    color = primero.ndim == 3
    # cámara de más de 8 bits: se bajan a 8 tirando los bits de abajo (los mismos en todo el video)
    correr = max(0, int(primero.max()).bit_length() - 8) if primero.dtype != np.uint8 else 0

    print(f'{len(archivos)} archivos, {total} cuadros de {ancho}x{alto}, '
          f'{fps:g} fps ({origen}) -> {total / fps:.0f} s de video')
    print(f'salida: {salida}')

    parcial = salida[:-4] + '.parcial.mp4'
    ffmpeg = subprocess.Popen(
        ['ffmpeg', '-y', '-loglevel', 'error',
         '-f', 'rawvideo', '-pix_fmt', 'rgb24' if color else 'gray', '-s', f'{ancho}x{alto}',
         '-r', f'{fps:g}', '-i', '-',
         # lados pares y yuv420p: lo que piden los reproductores de Windows; bt709 = colores bien
         '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2,scale=out_color_matrix=bt709,format=yuv420p',
         '-c:v', 'libx264', '-crf', '18', '-preset', 'medium',
         '-colorspace', 'bt709', '-color_primaries', 'bt709', '-color_trc', 'bt709',
         parcial],
        stdin=subprocess.PIPE)

    hechos = 0
    for k, (f, n) in enumerate(zip(archivos, paginas), 1):
        with Image.open(f) as im:                      # de a un cuadro: nunca hay más que eso en memoria
            for p in range(n):
                im.seek(p)
                cuadro = np.asarray(im)
                if cuadro.shape != primero.shape:
                    sys.exit(f'{f}, página {p}: tamaño {cuadro.shape}, esperaba {primero.shape}')
                if cuadro.dtype != np.uint8:
                    cuadro = np.clip(cuadro >> correr, 0, 255).astype(np.uint8)
                ffmpeg.stdin.write((cuadro[..., :3] if color else cuadro).tobytes())
        hechos += n
        print(f'  [{k}/{len(archivos)}] {os.path.basename(f)}: {hechos}/{total} cuadros '
              f'({100 * hechos / total:.0f} %)', flush=True)

    ffmpeg.stdin.close()
    if ffmpeg.wait() != 0:
        sys.exit('ffmpeg falló (el error está arriba); no se generó el video')
    os.replace(parcial, salida)
    print(f'listo: {salida} ({os.path.getsize(salida) / 1e6:.0f} MB)')


if __name__ == '__main__':
    main()
