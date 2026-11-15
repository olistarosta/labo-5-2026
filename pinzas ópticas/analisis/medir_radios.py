"""Radio aparente de cada partícula seguida, medido sobre la imagen (perfil radial).

    "C:/Program Files/Derivative/TouchDesigner/bin/python.exe" analisis/medir_radios.py "datos/<grabación>" [--cuadros 30] [--largo_min 100]

Hay que correrlo con el Python de TouchDesigner (lee el video con tifffile / cv2, que el Python
global no tiene). Lee trayectorias.csv y el video original (parametros.json → grabación importada).

Para cada partícula con al menos --largo_min cuadros medidos sin colisión:
  1. toma hasta --cuadros cuadros repartidos en su trayectoria;
  2. en cada uno recorta la imagen alrededor de (x, y) y calcula el perfil radial: la MEDIANA de la
     intensidad (canal verde) en anillos de 1 px centrados en la posición del tracking, dividida por
     la mediana del anillo exterior (el fondo local) menos 1. Fondo = 0; bead oscura = negativo.
     La mediana en el anillo aguanta una bead vecina o una mancha que tape parte del anillo;
  3. promedia los perfiles de todos los cuadros (la iluminación moteada se promedia porque la
     partícula se mueve) y mide sobre el perfil promedio y sobre el de cada cuadro:
       r_min     radio del mínimo del perfil (el anillo oscuro del borde; 0 si la bead es un disco
                 oscuro sin centro claro)
       r_borde   radio donde el perfil, yendo hacia afuera desde el mínimo, vuelve a la mitad de
                 su profundidad: el borde exterior. ES EL RADIO APARENTE que se informa.
       contraste profundidad del mínimo;  centro: valor del perfil en r = 0 (centro claro > 0:
                 la bead hace de lente, depende del foco).

Escribe:
  datos/<grabación>/radios.csv     una fila por partícula: n_cuadros, r_borde_px (del perfil promedio),
                                   err_r_borde_px (desvío entre cuadros / √n), disp_r_borde_px (desvío
                                   entre cuadros: cuánto cambia con el foco), r_min_px, contraste, centro
  datos/<grabación>/perfiles.csv   el perfil promedio de cada partícula (columna r_px y una por partícula)
  figuras/<grabación>/radios/recortes.png   seis recortes de cada partícula (una fila cada una, repartidos
                                   en el tiempo) con el círculo de r_borde: para MIRAR qué es cada una
                                   (bead sola, agregado, fuera de foco)

OJO: el radio aparente NO es el radio de la bead. Depende del foco (una bead desenfocada se ve más
grande) y de la difracción (~0.3 µm de borde). Sirve para comparar partículas entre sí y para
detectar tamaños distintos, no como medida absoluta mejor que ~0.2 µm.
"""
import argparse
import json
import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)
if RAIZ not in sys.path:
    sys.path.insert(0, RAIZ)
from pinzas import io  # noqa: E402  (agrega TouchDesigner/lib al path en el Python de TD)


def leer_trayectorias(ruta):
    """trayectorias.csv → dict de arrays (cuadro, particula, x_px, y_px, interpolado, colision); vacío = NaN.
    (El Python de TD no trae pandas.)"""
    import csv
    with open(ruta, encoding='utf-8', newline='') as f:
        filas = list(csv.DictReader(f))
    num = lambda c: np.array([float(r[c]) if r[c] != '' else np.nan for r in filas])
    return {c: num(c) for c in ('cuadro', 'particula', 'x_px', 'y_px', 'interpolado', 'colision')}


def escribir_csv(ruta, columnas):
    """dict nombre → array, a CSV (NaN = celda vacía)."""
    import csv
    with open(ruta, 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f)
        w.writerow(columnas)
        for fila in zip(*columnas.values()):
            w.writerow(['' if isinstance(v, float) and not np.isfinite(v) else ('%.5g' % v if isinstance(v, float) else v)
                        for v in fila])


def perfil_radial(img, x, y, r_max, r_fondo):
    """Mediana de img en anillos de 1 px alrededor de (x, y), relativa al fondo: p[r] = I(r)/I_fondo − 1.
    img: recorte 2D; (x, y): centro en px del recorte; r_fondo = (r0, r1): anillo que define el fondo.
    Devuelve p [r_max] (anillo k = radios k … k+1, centrado en k + 0.5) o None si el fondo es 0."""
    yy, xx = np.indices(img.shape)
    r = np.hypot(xx - x, yy - y)
    k = r.astype(int)
    fondo = np.median(img[(r >= r_fondo[0]) & (r < r_fondo[1])])
    if not fondo > 0:
        return None
    p = np.full(r_max, np.nan)
    for i in range(r_max):
        v = img[k == i]
        if v.size:
            p[i] = np.median(v)
    return p / fondo - 1


def medir_perfil(p, r_buscar):
    """r_min, r_borde, contraste y centro de un perfil radial p (ver el encabezado). NaN si no hay borde.
    r_buscar: radio máximo donde se busca el mínimo (para no agarrar una vecina)."""
    r = np.arange(len(p)) + 0.5
    nan = (np.nan, np.nan, np.nan, np.nan)
    if not np.isfinite(p).all():
        return nan
    s = np.convolve(p, [0.25, 0.5, 0.25], mode='same')          # suavizado mínimo (anillos de 1 px ruidosos)
    s[0], s[-1] = p[0], p[-1]
    i = int(np.argmin(s[:int(r_buscar)]))
    pmin = s[i]
    if pmin >= 0:
        return nan
    r_min = r[i]
    if 0 < i < len(s) - 1:                                        # parábola por los tres puntos del mínimo
        den = s[i - 1] - 2 * s[i] + s[i + 1]
        if den > 0:
            r_min = r[i] + 0.5 * (s[i - 1] - s[i + 1]) / den
    if i == 0:
        r_min = 0.0
    fuera = np.flatnonzero(s[i:] >= pmin / 2)
    if not len(fuera) or fuera[0] == 0:
        return nan
    j = i + fuera[0]
    r_borde = r[j - 1] + (pmin / 2 - s[j - 1]) / (s[j] - s[j - 1]) * (r[j] - r[j - 1])
    return r_min, r_borde, -pmin, p[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('datos', help='carpeta datos/<grabación> (con trayectorias.csv y parametros.json)')
    ap.add_argument('--cuadros', type=int, default=30, help='cuadros por partícula')
    ap.add_argument('--largo_min', type=int, default=100, help='cuadros medidos sin colisión para entrar')
    a = ap.parse_args()

    carpeta = os.path.abspath(a.datos)
    nombre = os.path.basename(carpeta)
    with open(os.path.join(carpeta, 'parametros.json'), encoding='utf-8') as f:
        par = json.load(f)
    R0 = float(par['bloques']['b5_realce']['Radio'])              # radio usado para detectar (px)
    g = io.leer_grabacion(par['grabacion'])
    tr = leer_trayectorias(os.path.join(carpeta, 'trayectorias.csv'))
    bien = (tr['colision'] == 0) & (tr['interpolado'] == 0) & np.isfinite(tr['x_px'])

    r_max = int(round(3.0 * R0))
    r_fondo = (2.3 * R0, 3.0 * R0)
    lado = r_max + 2
    pedidos = {}                                                   # cuadro → [(partícula, x, y)]
    ids = []
    for pid in np.unique(tr['particula']).astype(int):
        suyas = np.flatnonzero(bien & (tr['particula'] == pid))
        if len(suyas) < a.largo_min:
            continue
        ids.append(pid)
        sel = suyas[np.unique(np.linspace(0, len(suyas) - 1, a.cuadros).round().astype(int))]
        for c, x, y in zip(tr['cuadro'][sel], tr['x_px'][sel], tr['y_px'][sel]):
            pedidos.setdefault(int(c), []).append((pid, x, y))

    perfiles = {p: [] for p in ids}
    recortes = {}
    for c in sorted(pedidos):                                      # en orden: los videos se leen hacia adelante
        img = np.asarray(g.cuadro(c), np.float32)
        if img.ndim == 3:
            img = img[:, :, 1]                                     # verde, como en la conversión
        for pid, x, y in pedidos[c]:
            j, i = int(round(x)), int(round(y))
            if i - lado < 0 or j - lado < 0 or i + lado + 1 > img.shape[0] or j + lado + 1 > img.shape[1]:
                continue
            rec = img[i - lado:i + lado + 1, j - lado:j + lado + 1]
            p = perfil_radial(rec, x - (j - lado), y - (i - lado), r_max, r_fondo)
            if p is not None and np.isfinite(p).all():
                perfiles[pid].append(p)
                recortes.setdefault(pid, []).append((c, rec.copy(), x - (j - lado), y - (i - lado)))
    g.cerrar()

    filas, medios = [], {}
    for pid in ids:
        P = np.array(perfiles[pid])
        if len(P) < 3:
            continue
        medio = P.mean(0)
        medios[pid] = medio
        r_min, r_borde, contraste, centro = medir_perfil(medio, 1.8 * R0)
        cada = np.array([medir_perfil(p, 1.8 * R0)[1] for p in P])
        cada = cada[np.isfinite(cada)]
        filas.append(dict(particula=pid, n_cuadros=len(P), r_borde_px=r_borde,
                          err_r_borde_px=cada.std(ddof=1) / np.sqrt(len(cada)) if len(cada) > 2 else np.nan,
                          disp_r_borde_px=cada.std(ddof=1) if len(cada) > 2 else np.nan,
                          r_borde_mediana_px=np.median(cada) if len(cada) else np.nan,
                          r_min_px=r_min, contraste=contraste, centro=centro))
    claves = list(filas[0]) if filas else ['particula']
    escribir_csv(os.path.join(carpeta, 'radios.csv'),
                 {c: [int(f[c]) if c in ('particula', 'n_cuadros') else float(f[c]) for f in filas] for c in claves})
    escribir_csv(os.path.join(carpeta, 'perfiles.csv'),
                 dict(r_px=list(np.arange(r_max) + 0.5), **{str(p): [float(x) for x in v] for p, v in medios.items()}))

    # una imagen para mirar: un recorte de cada partícula con su círculo
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig_dir = os.path.join(RAIZ, 'figuras', nombre, 'radios')
    os.makedirs(fig_dir, exist_ok=True)
    n = len(filas)
    col = 6                                                        # recortes por partícula, repartidos en el tiempo
    fig, axs = plt.subplots(max(n, 1), col, figsize=(1.5 * col, 1.6 * max(n, 1)), squeeze=False)
    for ax in axs.ravel():
        ax.axis('off')
    th = np.linspace(0, 2 * np.pi, 100)
    for fila_ax, f in zip(axs, filas):
        todos = recortes[f['particula']]
        elegidos = [todos[k] for k in np.unique(np.linspace(0, len(todos) - 1, col).round().astype(int))]
        for ax, (c, rec, x, y) in zip(fila_ax, elegidos):
            ax.imshow(rec, cmap='gray')
            if np.isfinite(f['r_borde_px']):
                ax.plot(x + f['r_borde_px'] * np.cos(th), y + f['r_borde_px'] * np.sin(th), color='#E76F51', lw=0.8)
            ax.set_title('cuadro %d' % c, fontsize=7)
        fila_ax[0].text(-0.08, 0.5, 'partícula %d\n%.1f px' % (f['particula'], f['r_borde_px']), fontsize=8,
                        ha='right', va='center', transform=fila_ax[0].transAxes)
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, 'recortes.png'), dpi=150, bbox_inches='tight')
    print('  '.join(claves))
    for f in filas:
        print('  '.join('%d' % f[c] if c in ('particula', 'n_cuadros') else '%.2f' % f[c] for c in claves))
    print('%s: %d partículas -> radios.csv, perfiles.csv y %s' % (nombre, n, os.path.join(fig_dir, 'recortes.png')))


if __name__ == '__main__':
    main()
