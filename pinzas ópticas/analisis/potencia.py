"""Barrido de potencia con la partícula atrapada y la platina quieta (LASER 5 y BARRIDO POTENCIAS 2).

Lo usa `PINZAS 4 - barrido de potencia.ipynb`. Tres partes:

1. VIDEO  (`medir_video`): lee los TIFF originales (C:\\datos_pinzas\\GRABACIONES) alrededor del láser y
   mide en cada cuadro (a) cuánta luz roja del láser llega a la cámara, en el núcleo y lejos de él, y
   (b) la posición de la partícula otra vez, con el núcleo del láser TAPADO y en dos canales (verde y
   azul), para chequear las posiciones de `atrapada.csv` y estimar el ruido de localización.
2. SINTÉTICOS (`ou_sintetico`): Ornstein-Uhlenbeck muestreado como la cámara (promedio durante la
   exposición + ruido blanco), para validar los estimadores.
3. ESTIMADORES por escalón (`analizar_tramo`): envuelve `pinzas.trampa` y agrega la covarianza entre
   canales y el ajuste conjunto de la autocovarianza.

Unidades de este módulo: µm, s, pN/µm (1 pN/µm = 1e-6 N/m), kT, salvo donde se dice otra cosa.
"""
import os
import re
import sys

import numpy as np
import pandas as pd

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)
if RAIZ not in sys.path:
    sys.path.insert(0, RAIZ)

from pinzas import trampa as _tr                                    # noqa: E402
from pinzas.fisica import KB, factor_blur_covarianza as B_blur, factor_blur_varianza as S_blur  # noqa: E402
from pinzas_datos import CATALOGO, DATOS, ESCALA_UM                 # noqa: E402

GRABACIONES = r"C:\datos_pinzas\GRABACIONES"
SALIDA = os.path.join(DATOS, 'pinzas')                              # CSV de este análisis
BARRIDOS = ('LASER 5', 'BARRIDO POTENCIAS 2')
# cuadros en que la trampa está vacía y el láser apagado o casi: de ahí sale el fondo de la imagen
CUADROS_VACIOS = {'LASER 5': (2400, 2540), 'BARRIDO POTENCIAS 2': (2760, 3000)}
R_NUCLEO = 11.0          # px: radio que se tapa alrededor del láser al localizar
KT_PN_UM = lambda T: KB * T * 1e18                                  # kT en pN·µm


# ================================================================ 1. VIDEO
def _orden(n):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', os.path.splitext(n)[0])]


def indice_tiff(nombre):
    """Lista de páginas de la grabación: (ruta, desplazamiento de los píxeles, ancho, alto).
    Los TIFF de ThorCam son RGB de 8 bits sin compresión: se leen por memmap sin cargar el archivo."""
    from PIL import Image
    carpeta = os.path.join(GRABACIONES, nombre)
    paginas = []
    for f in sorted([f for f in os.listdir(carpeta) if f.lower().endswith('.tif')], key=_orden):
        ruta = os.path.join(carpeta, f)
        with Image.open(ruta) as im:
            ancho, alto = im.size
            for k in range(im.n_frames):
                im.seek(k)
                paginas.append((ruta, im.tile[0][2], ancho, alto))
    return paginas


def ventana(pagina, i0, j0, alto, ancho_v):
    """Recorte [alto, ancho_v, 3] (uint8) de una página, desde la fila i0 y la columna j0."""
    ruta, off, ancho, _ = pagina
    a = np.memmap(ruta, dtype=np.uint8, mode='r', offset=off + i0 * ancho * 3, shape=(alto, ancho, 3))
    return np.array(a[:, j0:j0 + ancho_v, :])


def simetria_radial_pesos(I, pesos=None):
    """Centro de simetría radial (Parthasarathy 2012), igual que `pinzas.tracking.simetria_radial`,
    pero con un peso extra por vértice (0 = ese gradiente no se usa). Sirve para tapar una zona.
    I: [Ny, Nx]; pesos: [Ny-1, Nx-1] o None. Devuelve (x, y) en px del recorte."""
    from scipy import ndimage
    I = np.asarray(I, float)
    Ny, Nx = I.shape
    xm, ym = np.meshgrid(np.arange(Nx - 1) - (Nx - 2) / 2.0, np.arange(Ny - 1) - (Ny - 2) / 2.0)
    dIdu = I[:-1, 1:] - I[1:, :-1]
    dIdv = I[:-1, :-1] - I[1:, 1:]
    fdu = ndimage.uniform_filter(dIdu, 3, mode='constant', cval=0.0)
    fdv = ndimage.uniform_filter(dIdv, 3, mode='constant', cval=0.0)
    g2 = fdu * fdu + fdv * fdv
    if pesos is not None:
        g2 = g2 * pesos
    with np.errstate(divide='ignore', invalid='ignore'):
        m = -(fdv + fdu) / (fdu - fdv)
    m = np.where(np.isfinite(m), m, 0.0)
    g2 = np.where(np.isfinite(g2), g2, 0.0)
    b = ym - m * xm
    s = g2.sum()
    if not s > 0:
        return np.nan, np.nan
    xcen, ycen = (g2 * xm).sum() / s, (g2 * ym).sum() / s
    w = g2 / np.maximum(np.hypot(xm - xcen, ym - ycen), 1e-6) / (m * m + 1)
    sw, smmw, smw = w.sum(), (m * m * w).sum(), (m * w).sum()
    smbw, sbw = (m * b * w).sum(), (b * w).sum()
    det = smw * smw - smmw * sw
    if det == 0 or not np.isfinite(det):
        return np.nan, np.nan
    return (smbw * sw - smw * sbw) / det + (Nx - 1) / 2.0, (smbw * smw - smmw * sbw) / det + (Ny - 1) / 2.0


def _cruce_rojo(pag, i0, j0, n, nuc3, cuadros):
    """Cuánto baja el verde (y el azul) por cada cuenta de rojo: la cámara entrega RGB ya procesado
    y el rojo del láser se RESTA de los otros canales. Regresión temporal por píxel, G(t) y B(t)
    contra R(t), en los píxeles a menos de 3 px del láser (ahí el rojo cambia mucho de un cuadro a
    otro). Devuelve (c_G, c_B) > 0."""
    V = np.array([ventana(pag[k], i0, j0, n, n).astype(np.float32)[nuc3] for k in cuadros])
    R = V[..., 0] - V[..., 0].mean(0)
    c = []
    for canal in (1, 2):
        Y = V[..., canal] - V[..., canal].mean(0)
        c.append(-float(np.median((R * Y).mean(0) / (R * R).mean(0))))
    return c[0], c[1]


def _localizar(rec, tx, ty, xc0, yc0, x_abs0, y_abs0, suavizado, r_tapa, r_max):
    """Simetría radial del recorte `rec` (ya dividido por el fondo), sin los vértices a menos de
    r_tapa px del láser (tx, ty) ni a más de r_max px del centro de partida (xc0, yc0).
    x_abs0, y_abs0: coordenadas absolutas del píxel [0, 0] del recorte. Devuelve (x, y) absolutos."""
    from scipy import ndimage
    I = ndimage.gaussian_filter(rec, suavizado, mode='nearest') if suavizado else rec
    ny, nx = I.shape
    yv, xv = np.mgrid[0:ny - 1, 0:nx - 1] + 0.5
    xa, ya = xv + x_abs0, yv + y_abs0
    peso = (np.hypot(xa - xc0, ya - yc0) < r_max).astype(float)
    if r_tapa:
        peso *= np.hypot(xa - tx, ya - ty) > r_tapa
    x, y = simetria_radial_pesos(I, peso)
    return x + x_abs0, y + y_abs0


ESTIMADORES = {
    # nombre: (canal, corregir el cruce del rojo, suavizado px, radio tapado en el láser px)
    'c': (1, True, 3.0, 8.0),       # el que se usa: verde corregido + núcleo tapado
    'u': (1, True, 3.4, 0.0),       # verde corregido, sin tapar
    'b': (2, True, 3.0, R_NUCLEO),  # azul corregido + núcleo tapado (ruido de píxel independiente)
    's': (1, False, 3.4, 0.0),      # verde tal cual, sin tapar: lo que hace la conversión a CSV
}


def medir_video(nombre, medio=110, lado=32, r_max=27.0, avisar=True):
    """Relee el video alrededor del láser. Devuelve (DataFrame con una fila por cuadro, dict de extras).

      cuadro, t_s
      rojo_nucleo   R/R_lejos en un disco de 6 px en el láser, menos su valor con la trampa vacía y
                    el láser apagado. R_lejos = media del rojo a más de 100 px (iluminación).
      rojo_halo     lo mismo en el anillo de 25 a 50 px (ahí el rojo vale pocas cuentas: no satura)
      rojo_lejos    R/G a más de 100 px del láser, menos su valor con el láser apagado: luz del láser
                    desparramada en todo el campo
      rojo_max      máximo del canal rojo a menos de 15 px del láser (cuentas, 0-255)
      x<e>_px, y<e>_px   posición de la partícula (px del sensor) por simetría radial sobre el canal
                    dividido por el fondo, con cuatro variantes e (ver ESTIMADORES):
                      c  verde + c_G·(rojo − rojo del fondo), sin los píxeles a < 8 px del láser
                      u  verde corregido igual, sin tapar nada
                      b  azul corregido, sin los píxeles a < 11 px del láser
                      s  verde tal cual, sin tapar (réplica de lo que hace la conversión)
    La posición de partida de cada cuadro es la de `atrapada.csv`; sin ella, queda vacío.
    extras: c_G, c_B (cruce del rojo), fondo [n, n, 3], plantillas para `ruido_por_plantilla`.
    """
    tx, ty = CATALOGO[nombre]['trampa_px']
    a = pd.read_csv(os.path.join(DATOS, nombre, 'atrapada.csv'))
    pag = indice_tiff(nombre)
    i0, j0 = int(round(ty)) - medio, int(round(tx)) - medio
    n = 2 * medio + 1
    yy, xx = np.mgrid[i0:i0 + n, j0:j0 + n]
    r = np.hypot(xx - tx, yy - ty)
    nuc, halo, lejos = r <= 6, (r > 25) & (r <= 50), r > 100
    v0, v1 = CUADROS_VACIOS[nombre]
    vacios = np.array([ventana(pag[k], i0, j0, n, n).astype(np.float32) for k in range(v0, v1)])
    fondo = np.median(vacios, axis=0)
    ref_n = np.mean([v[..., 0][nuc].mean() / v[..., 0][lejos].mean() for v in vacios])
    ref_h = np.mean([v[..., 0][halo].mean() / v[..., 0][lejos].mean() for v in vacios])
    ref_l = np.mean([v[..., 0][lejos].mean() / v[..., 1][lejos].mean() for v in vacios])
    t = a.t_s.values
    c_G, c_B = _cruce_rojo(pag, i0, j0, n, r <= 3, np.flatnonzero((t > 5) & (t < 145))[::3])
    cruce = {1: c_G, 2: c_B}
    xp = a.x_um.values / ESCALA_UM + a.x0_px.values
    yp = a.y_um.values / ESCALA_UM + a.y0_px.values
    filas = []
    for k, p in enumerate(pag):
        v = ventana(p, i0, j0, n, n).astype(np.float32)
        R, G = v[..., 0], v[..., 1]
        f = dict(cuadro=k, t_s=t[k],
                 rojo_nucleo=R[nuc].mean() / R[lejos].mean() - ref_n,
                 rojo_halo=R[halo].mean() / R[lejos].mean() - ref_h,
                 rojo_lejos=R[lejos].mean() / G[lejos].mean() - ref_l,
                 rojo_max=float(R[r <= 15].max()))
        for e in ESTIMADORES:
            f['x%s_px' % e] = f['y%s_px' % e] = np.nan
        if np.isfinite(xp[k]):
            jc, ic = int(round(xp[k])) - j0, int(round(yp[k])) - i0
            if lado <= ic < n - lado and lado <= jc < n - lado:
                sl = (slice(ic - lado, ic + lado + 1), slice(jc - lado, jc + lado + 1))
                exceso = v[sl][..., 0] - fondo[sl][..., 0]
                for e, (canal, corr, suav, tapa) in ESTIMADORES.items():
                    rec = (v[sl][..., canal] + (cruce[canal] * exceso if corr else 0.0)) / fondo[sl][..., canal]
                    f['x%s_px' % e], f['y%s_px' % e] = _localizar(rec, tx, ty, xp[k], yp[k], jc - lado + j0,
                                                                  ic - lado + i0, suav, tapa, r_max)
        filas.append(f)
        if avisar and k % 500 == 0:
            print('  %s: cuadro %d de %d' % (nombre, k, len(pag)), flush=True)
    d = pd.DataFrame(filas)
    extras = dict(c_G=c_G, c_B=c_B, fondo=fondo, vacios=vacios, origen=(i0, j0), trampa=(tx, ty))
    return d, extras


def ruido_por_plantilla(nombre, extras, cuadros, lado=32, r_max=27.0):
    """Ruido de localización que produce el RUIDO DE PÍXEL, medido en la misma grabación.

    Plantilla = promedio de la imagen (verde y azul corregidos del cruce del rojo, divididos por el
    fondo) sobre `cuadros` (un escalón con la partícula atrapada): la partícula sin ruido (apenas
    emborronada por su propio movimiento, < 1 px). A la plantilla se le suma el ruido REAL de la
    cámara en esos mismos píxeles: la diferencia entre dos cuadros seguidos con la trampa vacía,
    dividida por √2 y por el fondo (trae la correlación entre píxeles vecinos del procesado de color
    y no se entera de nada que cambie despacio). Se localiza cada imagen
    plantilla + ruido con los mismos estimadores y se devuelve el desvío estándar [px] de cada uno:
    dict {'xc': .., 'yc': .., 'xu': .., ...}.
    No incluye lo que agregue el láser (su mancha cambia de forma de un cuadro a otro) ni el ruido
    de fotones extra del rojo en el núcleo: es una cota INFERIOR del ruido de localización.
    """
    tx, ty = extras['trampa']
    i0, j0 = extras['origen']
    fondo, vacios = extras['fondo'], extras['vacios']
    pag = indice_tiff(nombre)
    n = fondo.shape[0]
    ic, jc = int(round(ty)) - i0, int(round(tx)) - j0
    sl = (slice(ic - lado, ic + lado + 1), slice(jc - lado, jc + lado + 1))
    acum = {1: 0.0, 2: 0.0}
    for k in cuadros:
        v = ventana(pag[k], i0, j0, n, n).astype(np.float32)[sl]
        exceso = v[..., 0] - fondo[sl][..., 0]
        acum[1] = acum[1] + (v[..., 1] + extras['c_G'] * exceso) / fondo[sl][..., 1]
        acum[2] = acum[2] + (v[..., 2] + extras['c_B'] * exceso) / fondo[sl][..., 2]
    plantilla = {c: acum[c] / len(cuadros) for c in acum}
    res = {e: [] for e in ESTIMADORES if e != 's'}
    dv = (vacios[1:] - vacios[:-1])[:, sl[0], sl[1]] / np.sqrt(2)      # ruido: diferencias de cuadros seguidos
    for d in dv[::2]:
        for e in res:
            canal, _, suav, tapa = ESTIMADORES[e]
            # el ruido del rojo también entra por la corrección del cruce
            ruido = (d[..., canal] + (extras['c_G'] if canal == 1 else extras['c_B']) * d[..., 0]) / fondo[sl][..., canal]
            res[e].append(_localizar(plantilla[canal] + ruido, tx, ty, tx, ty, jc - lado + j0, ic - lado + i0,
                                     suav, tapa, r_max))
    out = {}
    for e, v in res.items():
        v = np.array(v)
        out['x' + e], out['y' + e] = float(np.nanstd(v[:, 0], ddof=1)), float(np.nanstd(v[:, 1], ddof=1))
    out['n'] = len(dv[::2])
    return out


def ruido_objeto_fijo(nombre, extras, centro_px, cuadros, lado=32, r_max=27.0, suavizado=3.0):
    """Chequeo del método de la plantilla con algo que NO se mueve: una mancha de polvo del sensor
    (centro_px, la misma en todas las grabaciones). Devuelve (medido_x, medido_y, plantilla_x,
    plantilla_y) [px]: desvío de su posición cuadro a cuadro, y el que predice plantilla + ruido."""
    pag = indice_tiff(nombre)
    cx, cy = centro_px
    i0, j0 = int(round(cy)) - lado, int(round(cx)) - lado
    n = 2 * lado + 1
    V = np.array([ventana(pag[k], i0, j0, n, n).astype(np.float32)[..., 1] for k in cuadros])
    lejos = np.ones((n, n), bool)
    lejos[8:-8, 8:-8] = False
    V = V / V[:, lejos].mean(1)[:, None, None]                    # saca el parpadeo de la lámpara
    med = np.array([_localizar(v, -1e9, -1e9, cx, cy, j0, i0, suavizado, 0.0, r_max) for v in V])
    plantilla = V.mean(0)
    ruido = (V[1:] - V[:-1]) / np.sqrt(2)
    sim = np.array([_localizar(plantilla + r_, -1e9, -1e9, cx, cy, j0, i0, suavizado, 0.0, r_max) for r_ in ruido[::2]])
    dif = lambda u: float(np.sqrt(np.mean(np.diff(u) ** 2) / 2))
    return dif(med[:, 0]), dif(med[:, 1]), float(sim[:, 0].std(ddof=1)), float(sim[:, 1].std(ddof=1))


# ================================================================ escalones y tabla por cuadro
RETRASO_S = 1.5        # el cambio de potencia se ve en el video ~1.5 s después de cada múltiplo de 30 s
MARGEN_INI_S = 5.0     # se descarta desde el múltiplo de 30 s hasta 5 s después (3.5 s tras el cambio)
MARGEN_FIN_S = 0.5     # y se usa hasta 0.5 s después del múltiplo siguiente (1 s antes del cambio)
N_ATRAPADA = 5         # escalones (1..5: 30 a 7 mW) en que la partícula está atrapada todo el tramo


def escalones():
    """Los escalones del cuaderno con el tramo que se usa de cada uno.
    Columnas: grabacion, n, t_ini_s, t_fin_s (cuaderno), potencia_mW, potencia_texto, nota,
    t_cambio_s (cuándo cambia de verdad: nominal + RETRASO_S), t_ini_usado_s, t_fin_usado_s."""
    cu = pd.read_csv(os.path.join(DATOS, 'cuaderno', 'potencias_cuaderno.csv'))
    cu['t_cambio_s'] = cu.t_ini_s + np.where(cu.n > 1, RETRASO_S, 0.0)
    cu['t_ini_usado_s'] = cu.t_ini_s + MARGEN_INI_S
    cu['t_fin_usado_s'] = cu.t_fin_s + MARGEN_FIN_S
    return cu


def tabla_cuadros(recalcular=False, avisar=True):
    """Tabla por cuadro de las dos grabaciones (`datos/pinzas/potencia_cuadros.csv`) y el ruido de
    píxel por escalón (`potencia_ruido_pixel.csv`). Si ya existen se leen; si no (o recalcular=True)
    se relee el video (unos 5 minutos). Devuelve (cuadros, ruido)."""
    os.makedirs(SALIDA, exist_ok=True)
    r1, r2 = os.path.join(SALIDA, 'potencia_cuadros.csv'), os.path.join(SALIDA, 'potencia_ruido_pixel.csv')
    if not recalcular and os.path.exists(r1) and os.path.exists(r2):
        return pd.read_csv(r1), pd.read_csv(r2)
    esc = escalones()
    tablas, ruidos = [], []
    for g in BARRIDOS:
        if avisar:
            print('releyendo el video de', g, flush=True)
        d, ex = medir_video(g, avisar=False)
        a = pd.read_csv(os.path.join(DATOS, g, 'atrapada.csv'))
        x0, y0 = a.x0_px.values[0], a.y0_px.values[0]
        out = pd.DataFrame(dict(grabacion=g, cuadro=a.cuadro, t_s=a.t_s))
        e = esc[esc.grabacion == g]
        out['escalon'] = np.clip(np.searchsorted(e.t_cambio_s.values, a.t_s.values, side='right'), 1, 8)
        usable = np.zeros(len(a), bool)
        for _, f in e.iterrows():
            usable |= (a.t_s.values >= f.t_ini_usado_s) & (a.t_s.values < f.t_fin_usado_s) & (f.n <= N_ATRAPADA)
        out['usable'] = usable
        out['particula'] = a.particula
        out['x_conv_um'], out['y_conv_um'] = a.x_um, a.y_um
        for q in ESTIMADORES:
            out['x%s_um' % q] = (d['x%s_px' % q].values - x0) * ESCALA_UM
            out['y%s_um' % q] = (d['y%s_px' % q].values - y0) * ESCALA_UM
        for c in ('rojo_nucleo', 'rojo_halo', 'rojo_lejos', 'rojo_max'):
            out[c] = d[c].values
        out['plat_dx_um'], out['plat_dy_um'] = a.plat_dx_um, a.plat_dy_um
        tablas.append(out)
        for _, f in e[e.n <= N_ATRAPADA].iterrows():
            ks = np.flatnonzero((a.t_s.values >= f.t_ini_usado_s) & (a.t_s.values < f.t_fin_usado_s))
            rp = ruido_por_plantilla(g, ex, ks[::2])
            fila = dict(grabacion=g, escalon=int(f.n), c_G=ex['c_G'], c_B=ex['c_B'], n_cuadros_vacios=rp['n'])
            for q in 'cub':
                fila['ruido_x%s_um' % q] = rp['x' + q] * ESCALA_UM
                fila['ruido_y%s_um' % q] = rp['y' + q] * ESCALA_UM
            ruidos.append(fila)
        ks = np.flatnonzero((a.t_s.values >= 5) & (a.t_s.values < 30))
        for nombre_p, centro in (('polvo_1', (970.0, 400.0)), ('polvo_2', (462.5, 155.0))):
            mx, my, px_, py_ = ruido_objeto_fijo(g, ex, centro, ks)
            ruidos.append(dict(grabacion=g, escalon=0, c_G=np.nan, c_B=np.nan, n_cuadros_vacios=len(ks),
                               ruido_xc_um=mx * ESCALA_UM, ruido_yc_um=my * ESCALA_UM,
                               ruido_xu_um=px_ * ESCALA_UM, ruido_yu_um=py_ * ESCALA_UM,
                               ruido_xb_um=np.nan, ruido_yb_um=np.nan))
    cuadros, ruido = pd.concat(tablas, ignore_index=True), pd.DataFrame(ruidos)
    cuadros.to_csv(r1, index=False, float_format='%.5g')
    ruido.to_csv(r2, index=False, float_format='%.5g')
    return cuadros, ruido


# ================================================================ 2. SINTÉTICOS
def ou_sintetico(n, dt, tau_c, var, W, var_ruido, rng=None, deriva=None):
    """Ornstein-Uhlenbeck de varianza `var` [µm²] y tiempo de relajación tau_c [s], muestreado como
    la cámara: n cuadros cada dt, cada uno el PROMEDIO durante la exposición W (0 < W <= dt), más
    ruido blanco de varianza var_ruido. deriva: None o función t -> µm que se suma.
    Devuelve (t, x) [s, µm]. Usa `pinzas.sintetico.trayectoria` (paso exacto, 40 subpasos)."""
    from pinzas.sintetico import trayectoria
    rng = np.random.default_rng(rng)
    t = np.arange(n) * dt
    W = max(W, dt / 400)
    gamma = 1.0
    kappa = gamma / tau_c
    T = var * 1e-12 * kappa / KB                                     # para que kT/kappa = var
    x, _ = trayectoria(t, W, KB * T / gamma, kappa, gamma, T, sub=40, rng=rng)
    x = x * 1e6 + rng.standard_normal(n) * np.sqrt(var_ruido)
    if deriva is not None:
        x = x + deriva(t)
    return t, x


# ================================================================ 3. ESTIMADORES por tramo
def quitar_tendencia(t, x, grado=2):
    """x menos un polinomio de grado `grado` ajustado en el tramo (la deriva lenta del equilibrio)."""
    ok = np.isfinite(x)
    u = (t - t[ok].mean()) / (np.ptp(t[ok]) or 1.0)
    return x - np.polyval(np.polyfit(u[ok], x[ok], grado), u)


def analizar_tramo(t, x_um, W, T, var_ruido_um2, gamma=None, f_min_hz=0.2, grado=2, n_bloques=12,
                   n_boot=200, semilla=0):
    """Todos los estimadores de `pinzas.trampa` sobre un tramo de un eje, en las unidades de acá.

    t [s], x_um [µm] (un escalón, un eje). W: exposición supuesta [s]. T [K].
    var_ruido_um2: varianza del ruido de localización que se resta [µm²].
    gamma [N·s/m] o None: con gamma, la corrección de blur de equipartición/Boltzmann se resuelve
      autoconsistente (alpha = W κ/γ) y autocorrelación/espectro dan también κ = γ/τ_c; sin gamma,
      el blur usa el τ_c de la autocorrelación (malo si no se resuelve).
    Antes se resta un polinomio de grado `grado` (deriva lenta); el espectro se ajusta desde f_min_hz.
    Devuelve dict plano: n, dt, media_um, var_um2 (sin tendencia), var_cruda_um2, rho1 (correlación
    a un cuadro, de la serie sin tendencia), y para cada método m de trampa.METODOS:
    kappa_<m> [pN/µm], err_kappa_<m>, tau_<m> [s], err_tau_<m>; además var_psd_um2 (varianza
    térmica del ajuste del espectro), anarmonia_2s, chi2_boltzmann, asimetria, curtosis_exceso,
    y 'graficos' (lo de trampa, en SI).
    """
    t = np.asarray(t, float)
    x = np.asarray(x_um, float)
    ok = np.isfinite(x)
    t, x = t[ok], x[ok]
    xd = quitar_tendencia(t, x, grado) if grado is not None else x - x.mean()
    dt = float(np.median(np.diff(t)))
    tramos = np.array([[0, len(t) - 1]])
    est = _tr.analizar_eje(t, xd * 1e-6, W, T, gamma, tramos, var_ruido_um2 * 1e-12, True, None,
                           n_bloques, n_boot, semilla, f_min_hz, None)
    out = dict(n=len(x), dt=dt, media_um=float(x.mean()), var_um2=float(xd.var()), var_cruda_um2=float(x.var()),
               rho1=float(np.mean(xd[1:] * xd[:-1]) / xd.var()))
    for m in _tr.METODOS:
        out['kappa_' + m] = est[m]['kappa'] * 1e6
        out['err_kappa_' + m] = est[m]['err_kappa'] * 1e6
        out['tau_' + m] = est[m]['tau_c']
        out['err_tau_' + m] = est[m]['err_tau_c']
    out['var_psd_um2'] = est['espectro']['ajuste']['v'] * 1e12
    out['err_var_psd_um2'] = est['espectro']['ajuste']['err_v'] * 1e12
    out['var_msd_um2'] = est['msd']['ajuste']['v'] * 1e12
    out['amp_acf_um2'] = est['autocorrelacion']['amplitud'] * 1e12
    out['anarmonia_2s'] = est['armonicidad']['anarmonia_2s']
    out['chi2_boltzmann'] = est['armonicidad']['chi2_red']
    z = (xd - xd.mean()) / xd.std()
    out['asimetria'] = float(np.mean(z ** 3))
    out['curtosis_exceso'] = float(np.mean(z ** 4) - 3)
    out['graficos'] = est['graficos']
    out['x_sin_tendencia_um'] = xd
    return out


def autocovarianza(x, max_lag, n_bloques=12, n_boot=300, rng=0):
    """Autocovarianza C(k), k = 0..max_lag, de una serie ya centrada, con error por bootstrap de
    bloques contiguos (los pares se asignan al bloque del primer punto). Devuelve (C, err)."""
    x = np.asarray(x, float)
    n = len(x)
    bl = (np.arange(n) * n_bloques) // n
    S = np.zeros((n_bloques, max_lag + 1))
    N = np.zeros((n_bloques, max_lag + 1))
    for k in range(max_lag + 1):
        p = x[:n - k] * x[k:]
        S[:, k] = np.bincount(bl[:n - k], p, minlength=n_bloques)
        N[:, k] = np.bincount(bl[:n - k], minlength=n_bloques)
    C = S.sum(0) / N.sum(0)
    rng = np.random.default_rng(rng)
    w = rng.multinomial(n_bloques, np.full(n_bloques, 1 / n_bloques), size=n_boot)
    return C, ((w @ S) / (w @ N)).std(0, ddof=1)


def difusion_libre(t, x, y, lags=(1, 2, 3, 4), trozo_s=6.0, n_boot=400, rng=0):
    """D [µm²/s] de un tramo de trayectoria casi libre, por eje: pendiente/2 de la recta ajustada a
    la VARIANZA del desplazamiento (no su media cuadrática: así no entra la velocidad media de
    arrastre) contra el retardo, en los retardos `lags`. El tramo se parte en trozos de trozo_s y
    la varianza se toma dentro de cada trozo. La ordenada da 2·(ruido² − D·W/3).
    Devuelve dict: Dx, Dy, err_Dx, err_Dy (bootstrap de trozos), ord_x, ord_y [µm²], n."""
    t, x, y = (np.asarray(v, float) for v in (t, x, y))
    dt = np.median(np.diff(t))
    trozos = np.floor((t - t[0]) / trozo_s).astype(int)
    ids = np.unique(trozos)
    res = {}
    for eje, v in (('x', x), ('y', y)):
        S = np.zeros((len(ids), len(lags)))        # suma de (d - media del trozo)²
        N = np.zeros((len(ids), len(lags)))
        for a, i in enumerate(ids):
            m = trozos == i
            tv, vv = t[m], v[m]
            for b, L in enumerate(lags):
                d = vv[L:] - vv[:-L]
                bien = np.isfinite(d) & (np.abs((tv[L:] - tv[:-L]) / dt - L) < 0.3)
                if bien.sum() > 5:
                    S[a, b] = np.sum((d[bien] - d[bien].mean()) ** 2)
                    N[a, b] = bien.sum() - 1
        tau = np.array(lags) * dt
        ajuste = lambda w: np.polyfit(tau, (w @ S) / (w @ N), 1)
        p = ajuste(np.ones(len(ids)))
        g = np.random.default_rng(rng)
        W_ = g.multinomial(len(ids), np.full(len(ids), 1 / len(ids)), size=n_boot).astype(float)
        P = np.array([ajuste(w) for w in W_])
        res['D' + eje], res['err_D' + eje] = p[0] / 2, P[:, 0].std(ddof=1) / 2
        res['ord_' + eje] = p[1]
    res['n'] = int(np.isfinite(x).sum())
    return res


def flujo_libres(nombre, t_max=150.0, largo_min=300, lejos_px=120.0, t_min=0.0):
    """Velocidad media de arrastre de las partículas LIBRES de la grabación (recta x(t), y(t) de
    cada trayectoria con al menos largo_min cuadros entre t_min y t_max y a más de lejos_px del láser).
    Devuelve DataFrame: particula, n, vx_um_s, vy_um_s, x_px, y_px (posición media)."""
    tr = pd.read_csv(os.path.join(DATOS, nombre, 'trayectorias.csv'))
    tx, ty = CATALOGO[nombre]['trampa_px']
    tr = tr[tr.x_px.notna() & (tr.t_s >= t_min) & (tr.t_s < t_max) & (tr.colision == 0)]
    filas = []
    for pid, d in tr.groupby('particula'):
        if len(d) < largo_min or np.hypot(d.x_px - tx, d.y_px - ty).min() < lejos_px:
            continue
        filas.append(dict(particula=int(pid), n=len(d), vx_um_s=np.polyfit(d.t_s, d.x_um, 1)[0],
                          vy_um_s=np.polyfit(d.t_s, d.y_um, 1)[0], x_px=d.x_px.mean(), y_px=d.y_px.mean()))
    return pd.DataFrame(filas)


def radio_aparente(nombre, cuadros, x_px, y_px, corregir_rojo=False, c_G=0.0):
    """Radio aparente [px] de la partícula con la definición de `medir_radios.py` (r_borde: donde el
    perfil radial, yendo hacia afuera desde el anillo oscuro, recupera la mitad de su profundidad),
    sobre el verde dividido por el fondo. El mínimo se busca entre 8 y 32 px (así no se toma el punto
    oscuro que deja el láser en el centro). Devuelve (r_borde del perfil promedio, desvío entre
    cuadros, r del mínimo, profundidad, perfil promedio [40])."""
    from medir_radios import perfil_radial
    pag = indice_tiff(nombre)
    tx, ty = CATALOGO[nombre]['trampa_px']
    medio = 110
    i0, j0 = int(round(ty)) - medio, int(round(tx)) - medio
    n = 2 * medio + 1
    v0, v1 = CUADROS_VACIOS[nombre]
    fondo = np.median([ventana(pag[k], i0, j0, n, n).astype(np.float32) for k in range(v0, v1, 3)], axis=0)
    perfiles = []
    for k, x, y in zip(cuadros, x_px, y_px):
        if not np.isfinite(x) or min(x - j0, y - i0, j0 + n - x, i0 + n - y) < 52:
            continue
        v = ventana(pag[k], i0, j0, n, n).astype(np.float32)
        img = (v[..., 1] + (c_G * (v[..., 0] - fondo[..., 0]) if corregir_rojo else 0.0)) / fondo[..., 1]
        jc, ic = int(round(x)) - j0, int(round(y)) - i0
        p = perfil_radial(img[ic - 50:ic + 51, jc - 50:jc + 51], x - j0 - (jc - 50), y - i0 - (ic - 50), 40, (40, 50))
        if p is not None and np.isfinite(p).all():
            perfiles.append(p)
    P = np.array(perfiles)

    def borde(p):
        s = np.convolve(p, [0.25, 0.5, 0.25], mode='same')
        i = 8 + int(np.argmin(s[8:32]))
        fuera = np.flatnonzero(s[i:] >= s[i] / 2)
        if s[i] >= 0 or not len(fuera) or fuera[0] == 0:
            return np.nan, np.nan, np.nan
        j = i + fuera[0]
        return (j - 0.5) + (s[i] / 2 - s[j - 1]) / (s[j] - s[j - 1]), i + 0.5, -s[i]
    rb, rm, prof = borde(P.mean(0))
    cada = np.array([borde(p)[0] for p in P])
    return rb, float(np.nanstd(cada, ddof=1)), rm, prof, P.mean(0)


def verosimilitud_exposicion(series, W_rel, var_ruido, f_min_hz=0.2):
    """¿Qué exposición prefieren los datos? Para cada W = W_rel·dt se ajusta el espectro de cada serie
    (v y τ_c libres, ruido fijo) con el modelo de `trampa.modelo_psd` y se suma −ln L (Whittle).
    series: lista de (t, x_um ya sin tendencia); var_ruido: lista de varianzas [µm²].
    Devuelve array [len(W_rel)] con 2·(−ln L) menos su mínimo (se lee como un χ²: +1 = 1 sigma)."""
    datos = []
    for (t, x), s2 in zip(series, var_ruido):
        seg = [(0, len(t) - 1)]
        f, P, dts = _tr.periodogramas(t, x * 1e-6, seg)
        m = f >= f_min_hz
        datos.append((f[m], P[m], dts[m], s2 * 1e-12, float(np.var(x)) * 1e-12, float(np.median(np.diff(t)))))
    tot = np.zeros(len(W_rel))
    for i, wr in enumerate(W_rel):
        for f, P, dts, s2, var, dt in datos:
            a = _tr.ajuste_espectro(f, P, dts, wr * dt, max(var - s2, 0.1 * var), 2 * dt, s2, True)
            Pm = _tr.modelo_psd(f, dts, a['v'], a['tau_c'], s2, wr * dt, True)
            tot[i] += 2 * np.sum(np.log(Pm) + P / Pm)
    return tot - tot.min()


def validar_sinteticos(casos, n_rep=60, T=296.15, gamma=None, semilla=1, **kw):
    """Corre `analizar_tramo` sobre series sintéticas y compara con la verdad.
    casos: lista de dict con n, dt, tau_c [s], var [µm²], W [s], var_ruido [µm²] y, opcional,
    W_supuesto (la exposición que se le dice al análisis; por defecto la verdadera) y
    'deriva' (amplitud [µm] de una deriva lenta parabólica que se suma).
    Devuelve DataFrame: una fila por caso y estimador con verdad, media, sesgo_% y dispersion_%."""
    rng = np.random.default_rng(semilla)
    filas = []
    for i, c in enumerate(casos):
        kT = KT_PN_UM(T)
        kappa = kT / c['var']
        g = kappa * 1e-6 * c['tau_c'] if gamma is None else gamma     # N·s/m coherente con el caso
        res = []
        for _ in range(n_rep):
            der = None
            if c.get('deriva'):
                a1, a2 = rng.standard_normal(2) * c['deriva']
                der = (lambda t, a1=a1, a2=a2, L=c['n'] * c['dt']: a1 * (t / L - 0.5) + a2 * (t / L - 0.5) ** 2 * 4)
            t, x = ou_sintetico(c['n'], c['dt'], c['tau_c'], c['var'], c['W'], c['var_ruido'], rng, der)
            a = analizar_tramo(t, x, c.get('W_supuesto', c['W']), T, c['var_ruido'], g, n_boot=0, **kw)
            res.append(a)
        for nombre, verdad in [('kappa_' + m, kappa) for m in _tr.METODOS] + \
                              [('tau_' + m, c['tau_c']) for m in ('autocorrelacion', 'espectro', 'msd')]:
            v = np.array([r[nombre] for r in res], float)
            ok = np.isfinite(v)
            filas.append(dict(caso=i, descripcion=c.get('nombre', ''), estimador=nombre, verdad=verdad,
                              mediana=np.median(v[ok]) if ok.any() else np.nan,
                              sesgo_pct=100 * (np.median(v[ok]) / verdad - 1) if ok.any() else np.nan,
                              dispersion_pct=100 * 0.7413 * np.subtract(*np.percentile(v[ok], [75, 25])) / verdad if ok.any() else np.nan,
                              fallos_pct=100 * (1 - ok.mean())))
    return pd.DataFrame(filas)


def exposicion_por_estirado(nombre, cuadro_ini, n_cuadros, direccion, ancho=7):
    """Exposición por el estirado (motion blur) de las partículas durante una patada de la platina.

    En cada cuadro de la patada se toma el perfil de la imagen (verde / fondo, sin lo que cambia
    despacio) a lo largo de la dirección de movimiento `direccion` = (ux, uy), promediado en una
    franja de 2·ancho + 1 px, partiendo de la posición de cada partícula aislada en `cuadro_ini`
    (todavía quieta). El tramo oscuro del perfil es el trazo de la partícula: su ancho a media
    profundidad es v·W (si es mayor que el diámetro) y su centro es la posición a mitad de la
    exposición, así que v·dt = (centro[k+1] − centro[k−1]) / 2. Cociente = W/dt.
    Devuelve DataFrame: cuadro, particula (índice), v_dt_px, largo_px, diametro_px, W_sobre_dt.
    Solo quedan los trazos de más de 1.6 diámetros con v·dt > 40 px. Es una medición ruidosa
    (los trazos de partículas vecinas se pisan): mirar la mediana.
    """
    from scipy import ndimage
    pag = indice_tiff(nombre)
    cuadro = lambda k: ventana(pag[k], 0, 0, pag[k][3], pag[k][2])[..., 1].astype(np.float32)
    fondo = np.median([cuadro(i) for i in np.linspace(0, len(pag) - 1, 25).astype(int)], axis=0)
    tr = pd.read_csv(os.path.join(DATOS, nombre, 'trayectorias.csv'))
    tx, ty = CATALOGO[nombre]['trampa_px']
    P = tr[(tr.cuadro == cuadro_ini) & tr.x_px.notna()][['x_px', 'y_px']].values
    dist = np.hypot(P[:, None, 0] - P[None, :, 0], P[:, None, 1] - P[None, :, 1]) + np.eye(len(P)) * 1e9
    P = P[(dist.min(1) > 90) & (np.hypot(P[:, 0] - tx, P[:, 1] - ty) > 80)]
    u = np.array(direccion, float) / np.hypot(*direccion)
    perp = np.array([-u[1], u[0]])
    s = np.arange(-60, 1300)
    res = {i: [] for i in range(len(P))}
    for k in range(cuadro_ini, cuadro_ini + n_cuadros):
        g = cuadro(k) / fondo
        g = ndimage.gaussian_filter(g - ndimage.gaussian_filter(g, 60), 1.5)
        for i, (x, y) in enumerate(P):
            prof = np.zeros(len(s))
            for a in range(-ancho, ancho + 1):
                prof += ndimage.map_coordinates(g, [y + s * u[1] + a * perp[1], x + s * u[0] + a * perp[0]],
                                                order=1, mode='constant', cval=np.nan)
            prof /= 2 * ancho + 1
            if np.isfinite(prof).sum() < 100 or np.nanmin(prof) > -0.012:
                res[i].append((np.nan, np.nan))
                continue
            j = int(np.nanargmin(prof))
            a0 = a1 = j
            while a0 > 0 and prof[a0 - 1] < prof[j] / 2:
                a0 -= 1
            while a1 < len(s) - 1 and prof[a1 + 1] < prof[j] / 2:
                a1 += 1
            res[i].append((0.5 * (s[a0] + s[a1]), s[a1] - s[a0] + 1))
    filas = []
    for i, r in res.items():
        c, L = np.array(r).T
        for q in range(1, n_cuadros - 1):
            vdt = (c[q + 1] - c[q - 1]) / 2
            if np.isfinite(vdt) and L[q] > 1.6 * L[0] and vdt > 40:
                filas.append(dict(cuadro=cuadro_ini + q, particula=i, v_dt_px=vdt, largo_px=L[q], diametro_px=L[0],
                                  W_sobre_dt=L[q] / vdt))
    return pd.DataFrame(filas)


def ruido_de_referencia(ruido, grabacion, eje, estimador='c'):
    """Varianza [µm²] del ruido de píxel que se resta: promedio de los escalones 1 y 2 (ahí la
    plantilla está menos emborronada por el movimiento de la propia partícula)."""
    r = ruido[(ruido.grabacion == grabacion) & ruido.escalon.isin([1, 2])]
    return float((r['ruido_%s%s_um' % (eje, estimador)] ** 2).mean())


def analizar_escalones(cuadros, ruido, gamma, T=296.15, W_rel=(0.0, 0.5, 1.0), estimador='c',
                       err_rel_ruido=0.2, n_boot=200, f_min_hz=0.2):
    """Corre `analizar_tramo` en cada grabación, escalón atrapado (1..N_ATRAPADA), eje y exposición
    supuesta W = W_rel·dt. Devuelve (tabla, detalles):
      tabla: una fila por (grabacion, escalon, eje, W_rel) con todo lo de `analizar_tramo` más
             potencia_mW, t_ini_usado_s, t_fin_usado_s, ruido_um2, f_c_<m> [Hz] = 1/(2π τ_<m>),
             err_kappa_ruido (cuánto cambia κ de equipartición si el ruido restado cambia ±20 %),
             rojo_nucleo, rojo_halo, rojo_lejos (medias del tramo), sd_rojo_nucleo, corr_xy.
      detalles[(grabacion, escalon, eje, W_rel)] = dict con 'graficos' y la serie sin tendencia.
    Con W = 0 no se usa gamma (no hay blur que corregir)."""
    esc = escalones()
    filas, det = [], {}
    for g in BARRIDOS:
        d = cuadros[cuadros.grabacion == g]
        for _, e in esc[(esc.grabacion == g) & (esc.n <= N_ATRAPADA)].iterrows():
            m = ((d.t_s >= e.t_ini_usado_s) & (d.t_s < e.t_fin_usado_s)).values
            t = d.t_s.values[m]
            dt = float(np.median(np.diff(t)))
            xy = {eje: d['%s%s_um' % (eje, estimador)].values[m] for eje in 'xy'}
            rxy = float(np.corrcoef(quitar_tendencia(t, xy['x']), quitar_tendencia(t, xy['y']))[0, 1])
            for eje in 'xy':
                s2 = ruido_de_referencia(ruido, g, eje, estimador)
                for wr in W_rel:
                    a = analizar_tramo(t, xy[eje], wr * dt, T, s2, gamma if wr > 0 else None,
                                       f_min_hz=f_min_hz, n_boot=n_boot, semilla=int(e.n) * 10 + (eje == 'y'))
                    det[(g, int(e.n), eje, wr)] = dict(graficos=a.pop('graficos'), x=a.pop('x_sin_tendencia_um'), t=t)
                    fila = dict(grabacion=g, escalon=int(e.n), potencia_mW=e.potencia_mW, eje=eje, W_rel=wr,
                                W_s=wr * dt, t_ini_usado_s=e.t_ini_usado_s, t_fin_usado_s=e.t_fin_usado_s,
                                ruido_um2=s2)
                    fila.update(a)
                    for met in ('autocorrelacion', 'espectro', 'msd'):
                        fila['f_c_' + met] = 1 / (2 * np.pi * a['tau_' + met]) if a['tau_' + met] > 0 else np.nan
                    vt = a['var_um2'] - s2
                    fila['err_kappa_ruido'] = a['kappa_equiparticion'] * err_rel_ruido * s2 / vt if vt > 0 else np.nan
                    for c in ('rojo_nucleo', 'rojo_halo', 'rojo_lejos'):
                        fila[c] = float(d[c].values[m].mean())
                    fila['sd_rojo_nucleo'] = float(d.rojo_nucleo.values[m].std())
                    fila['corr_xy'] = rxy
                    filas.append(fila)
    return pd.DataFrame(filas), det


def nulo_armonico(n, dt, tau_c, var, W, var_ruido, n_rep=300, rng=0, grado=2):
    """Asimetría y curtosis (exceso) de series sintéticas de un pozo ARMÓNICO con los mismos
    parámetros y el mismo largo, tratadas igual (sin tendencia): contra eso se compara lo medido.
    Devuelve (asimetrias [n_rep], curtosis [n_rep])."""
    rng = np.random.default_rng(rng)
    A, K = [], []
    for _ in range(n_rep):
        t, x = ou_sintetico(n, dt, tau_c, var, W, var_ruido, rng)
        z = quitar_tendencia(t, x, grado)
        z = (z - z.mean()) / z.std()
        A.append(np.mean(z ** 3))
        K.append(np.mean(z ** 4) - 3)
    return np.array(A), np.array(K)


def ajuste_lineal(x, y, err):
    """Recta y = a + b x por cuadrados mínimos pesados. Devuelve dict a, b, err_a, err_b, chi2, gl,
    x0 (= −a/b, dónde cruza el cero) y err_x0."""
    x, y, err = (np.asarray(v, float) for v in (x, y, err))
    w = 1 / err ** 2
    X = np.vstack([np.ones_like(x), x]).T
    cov = np.linalg.inv(X.T @ (X * w[:, None]))
    a, b = cov @ (X.T @ (w * y))
    chi2 = float(np.sum(w * (y - a - b * x) ** 2))
    x0 = -a / b
    err_x0 = abs(x0) * np.sqrt(cov[0, 0] / a ** 2 + cov[1, 1] / b ** 2 - 2 * cov[0, 1] / (a * b))
    return dict(a=a, b=b, err_a=np.sqrt(cov[0, 0]), err_b=np.sqrt(cov[1, 1]), chi2=chi2, gl=len(x) - 2,
                x0=x0, err_x0=err_x0)
