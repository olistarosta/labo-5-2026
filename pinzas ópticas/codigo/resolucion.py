"""Resolución del video, precisión con la que se mide el centro de una partícula, y forma del contorno.

Lo usa `resolución y precisión.ipynb`. Cuatro partes:

1. ÓPTICA    `borde_reticula`: qué tan nítida es la imagen, con los bordes de las líneas de la retícula.
2. RUIDO     `ruido_de_pixel`: cuánto cambia un píxel de un cuadro al siguiente sin que cambie la escena.
3. CENTRO    `quietas` + `precision_pares`: partículas pegadas al vidrio, que no se mueven: lo que
             "se mueve" su centro medido es el error de medición.
             `plantilla` + `prueba_corrimientos`: se corre la imagen de una partícula real una fracción
             conocida de píxel, se le suma ruido real de la cámara y se la localiza con el mismo
             código de la conversión: da el error sistemático y el aleatorio.
             `cota_teorica`: lo mejor que se podría lograr con ese ruido y esa imagen.
4. FORMA     `contorno`: radio del borde en función del ángulo → cuánto se aparta de un círculo.

Lee los TIFF originales (C:\\datos_pinzas\\GRABACIONES) con `trampa_roi.paginas` y localiza con el código de
los bloques de TouchDesigner (`trampa_roi.detectar`): el mismo de la conversión a CSV.
Unidades: px del sensor, salvo donde dice nm (1 px = ESCALA_UM µm).
"""
import os

import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage

from pinzas_datos import CATALOGO, DATOS, ESCALA_UM, RAIZ, cargar
import trampa_roi as tr

NM_PX = ESCALA_UM * 1e3
MEDIO = tr.VENTANA                                   # medio lado de la ventana de localización (px)


# ---------------------------------------------------------------- 1. óptica: bordes de la retícula
def borde_reticula(periodo_px=105.956, filas=(600, 1000), paso=0.25):
    """Perfil promedio de una línea de la retícula, sobremuestreado: como las líneas están inclinadas,
    cada fila corta el borde en una fracción distinta de píxel, y juntando 400 filas el perfil queda
    cada 0.25 px. Devuelve dict:
      u [px] (posición respecto del centro de la línea), perfil (intensidad relativa a la iluminación),
      pendiente [1/px] (derivada del perfil, suavizada),
      fwhm_borde_px   ancho a media altura de la derivada en cada borde, promedio de los dos lados:
                      qué tan ancho queda en la imagen un escalón perfectamente abrupto
      fwhm_por_lado, ancho_linea_px (distancia entre los dos bordes), contraste (profundidad de la línea)."""
    a = np.asarray(Image.open(os.path.join(RAIZ, 'CALIBRACION 0,01mm.tif')).convert('RGB'), float)[:, :, 1]
    x = np.arange(a.shape[1])
    # dónde está el centro de las líneas en cada fila: fase del patrón en franjas de 20 filas → recta
    cent, fase = [], []
    for y0 in range(filas[0], filas[1], 20):
        p = a[y0:y0 + 20].mean(0)
        p = p - ndimage.uniform_filter1d(p, 201, mode='nearest')
        fase.append(np.angle(np.sum(-p * np.exp(-2j * np.pi * x / periodo_px))))        # −p: las líneas son oscuras
        cent.append(y0 + 10)
    corr = np.polyfit(cent, -np.unwrap(fase) * periodo_px / (2 * np.pi), 1)              # centro(y) = corr[0]·y + corr[1]
    us, vs = [], []
    for y in range(filas[0], filas[1]):
        fila = a[y] / ndimage.uniform_filter1d(a[y], 301, mode='nearest')               # sin la iluminación lenta
        us.append((x - np.polyval(corr, y) + periodo_px / 2) % periodo_px - periodo_px / 2)
        vs.append(fila)
    us, vs = np.concatenate(us), np.concatenate(vs)
    bordes = np.arange(-periodo_px / 2, periodo_px / 2 + paso, paso)
    k = np.digitize(us, bordes) - 1
    ok = (k >= 0) & (k < len(bordes) - 1)
    n = np.bincount(k[ok], minlength=len(bordes) - 1)
    perfil = np.bincount(k[ok], vs[ok], len(bordes) - 1) / np.maximum(n, 1)
    u = 0.5 * (bordes[1:] + bordes[:-1])
    pend = np.gradient(ndimage.uniform_filter1d(perfil, 5), u)                           # suavizado de 1.25 px
    fwhm, pos = [], []
    for lado in (-1, 1):
        d = np.where(u * lado > 0, pend * lado, 0.0)                                     # borde izquierdo: baja; derecho: sube
        i = int(np.argmax(d))
        j0 = i
        while j0 > 0 and d[j0 - 1] > d[i] / 2:
            j0 -= 1
        j1 = i
        while j1 < len(d) - 1 and d[j1 + 1] > d[i] / 2:
            j1 += 1
        fwhm.append((j1 - j0 + 1) * paso)
        pos.append(u[i])
    centro = 0.5 * (pos[0] + pos[1])
    return dict(u=u - centro, perfil=perfil, pendiente=pend, fwhm_borde_px=float(np.mean(fwhm)), fwhm_por_lado=fwhm,
                ancho_linea_px=float(pos[1] - pos[0]), contraste=float(1 - perfil.min() / np.percentile(perfil, 95)))


# ---------------------------------------------------------------- cuadros y fondo
def cuadro_verde(pag):
    """Un cuadro entero (canal verde, float32)."""
    ruta, off, ancho = pag
    a = np.memmap(ruta, dtype=np.uint8, mode='r', offset=off, shape=(1024, ancho, 3))
    return np.asarray(a[:, :, 1], np.float32)


def fondo_mediana(pags, n=41):
    """Mediana temporal de n cuadros repartidos en la grabación: iluminación y lo que no se mueve nunca."""
    idx = np.unique(np.linspace(0, len(pags) - 1, n).round().astype(int))
    f = np.median(np.stack([cuadro_verde(pags[i]) for i in idx]), axis=0)
    return np.where(f > 0, f, 1).astype(np.float32)


def recorte(img, x, y, m):
    """Recorte de (2m+1)² px centrado en el píxel más cercano a (x, y). None si se sale de la imagen."""
    j, i = int(round(x)), int(round(y))
    if i - m < 0 or j - m < 0 or i + m + 1 > img.shape[0] or j + m + 1 > img.shape[1]:
        return None
    return img[i - m:i + m + 1, j - m:j + m + 1]


# ---------------------------------------------------------------- 3a. partículas quietas
def quietas(nombre, n_min=150, desvio_max_px=0.6, lejos_trampa_px=80.0, margen=5):
    """Partículas que no se mueven mientras la platina está quieta (pegadas al vidrio). Una fila por
    trayectoria: particula, n, cuadro_ini, cuadro_fin, x_px, y_px, sx_px, sy_px (desvío del centro
    medido, con la recta de deriva restada), z (contraste)."""
    d = cargar(nombre)
    t, pl = d['tr'], d['platina']
    paso = np.hypot(pl.dx_px, pl.dy_px).values
    quieta = np.convolve(paso > 1.5, np.ones(2 * margen + 1), 'same') == 0
    t = t[t.x_px.notna() & (t.colision == 0) & quieta[t.cuadro.values]]
    tx, ty = d['trampa_px']
    filas = []
    for pid, g in t.groupby('particula'):
        if len(g) < n_min or np.hypot(g.x_px.mean() - tx, g.y_px.mean() - ty) < lejos_trampa_px:
            continue
        # un solo tramo quieto: el más largo de cuadros seguidos
        c = g.cuadro.values
        cortes = np.flatnonzero(np.diff(c) > 3)
        ini = np.r_[0, cortes + 1]
        fin = np.r_[cortes, len(c) - 1]
        k = int(np.argmax(fin - ini))
        s = g.iloc[ini[k]:fin[k] + 1]
        if len(s) < n_min:
            continue
        rx = s.x_px.values - np.polyval(np.polyfit(s.cuadro, s.x_px, 1), s.cuadro)
        ry = s.y_px.values - np.polyval(np.polyfit(s.cuadro, s.y_px, 1), s.cuadro)
        if rx.std() > desvio_max_px or ry.std() > desvio_max_px:
            continue
        filas.append(dict(grabacion=nombre, particula=int(pid), n=len(s), cuadro_ini=int(s.cuadro.min()), cuadro_fin=int(s.cuadro.max()),
                          x_px=s.x_px.mean(), y_px=s.y_px.mean(), sx_px=rx.std(), sy_px=ry.std(), z=s.z.median()))
    return pd.DataFrame(filas)


def precision_pares(nombre, tabla=None):
    """Precisión real del centro con pares de partículas quietas que conviven en el mismo tramo.
    Si las dos están pegadas al vidrio, la distancia entre ellas no cambia: lo que cambia la
    distancia MEDIDA es el error de las dos mediciones. sigma de una = desvío de la diferencia / raíz
    de 2. La diferencia cancela además lo que mueve a las dos juntas (vibración o deriva de la platina).
    Devuelve una fila por par: a, b, n, separacion_px, sigma_x_px, sigma_y_px (de una partícula, por
    pares), comun_x_px, comun_y_px (movimiento que comparten: raíz de la covarianza)."""
    q = quietas(nombre) if tabla is None else tabla
    t = cargar(nombre)['tr']
    filas = []
    for i in range(len(q)):
        for j in range(i + 1, len(q)):
            a, b = q.iloc[i], q.iloc[j]
            c0, c1 = max(a.cuadro_ini, b.cuadro_ini), min(a.cuadro_fin, b.cuadro_fin)
            if c1 - c0 < 100:
                continue
            ta = t[(t.particula == a.particula) & (t.cuadro >= c0) & (t.cuadro <= c1) & (t.colision == 0)].set_index('cuadro')
            tb = t[(t.particula == b.particula) & (t.cuadro >= c0) & (t.cuadro <= c1) & (t.colision == 0)].set_index('cuadro')
            m = ta.join(tb, lsuffix='_a', rsuffix='_b', how='inner').dropna(subset=['x_px_a', 'x_px_b'])
            if len(m) < 100:
                continue
            fila = dict(grabacion=nombre, a=int(a.particula), b=int(b.particula), n=len(m),
                        separacion_px=float(np.hypot(a.x_px - b.x_px, a.y_px - b.y_px)))
            for eje in 'xy':
                ca, cb = m[eje + '_px_a'].values, m[eje + '_px_b'].values
                quitar = lambda v: v - np.polyval(np.polyfit(m.index, v, 1), m.index)
                ra, rb = quitar(ca), quitar(cb)
                fila['sigma_%s_px' % eje] = float(np.std(ra - rb) / np.sqrt(2))
                cov = float(np.mean(ra * rb))
                fila['comun_%s_px' % eje] = float(np.sign(cov) * np.sqrt(abs(cov)))
            filas.append(fila)
    return pd.DataFrame(filas)


# ---------------------------------------------------------------- 2 y 3b. plantilla, ruido y prueba con corrimientos
def plantilla(pags, fondo, fila, m=40, n_cuadros=120):
    """Imagen promedio de una partícula quieta (cuadros divididos por el fondo; 1 = fondo), de
    (2m+1)² px, y cuántos cuadros se promediaron. `fila`: una fila de `quietas`."""
    cuadros = np.unique(np.linspace(fila.cuadro_ini, fila.cuadro_fin, n_cuadros).round().astype(int))
    pila = []
    for c in cuadros:
        r = recorte(cuadro_verde(pags[c]) / fondo, fila.x_px, fila.y_px, m)
        if r is not None:
            pila.append(r)
    return np.mean(pila, axis=0), len(pila)


def recortes_vacios(nombre, pags, fondo, cuantos=60, lejos_px=150.0, semilla=0):
    """Recortes de (2·MEDIO+1)² px de zonas SIN partículas seguidas (cuadros divididos por el fondo):
    el ruido real de la cámara, con lo que quede de iluminación despareja. Devuelve [cuantos, L, L]
    y el nivel medio del fondo en cuentas."""
    d = cargar(nombre)
    t = d['tr'][d['tr'].x_px.notna()]
    tx, ty = d['trampa_px']
    rng = np.random.default_rng(semilla)
    por_cuadro = {c: g[['x_px', 'y_px']].values for c, g in t.groupby('cuadro')}
    salida, niveles, intentos = [], [], 0
    while len(salida) < cuantos and intentos < 40 * cuantos:
        intentos += 1
        c = int(rng.integers(0, len(pags)))
        x = float(rng.uniform(MEDIO + 5, 1280 - MEDIO - 5))
        y = float(rng.uniform(MEDIO + 5, 1024 - MEDIO - 5))
        pts = por_cuadro.get(c)
        if pts is None or np.hypot(x - tx, y - ty) < lejos_px + MEDIO:
            continue
        if np.min(np.hypot(pts[:, 0] - x, pts[:, 1] - y)) < lejos_px:
            continue
        salida.append(recorte(cuadro_verde(pags[c]) / fondo, x, y, MEDIO))
        niveles.append(float(recorte(fondo, x, y, MEDIO).mean()))
    return np.stack(salida), float(np.mean(niveles))


def ruido_de_pixel(pags, fondo, cuadros, zona=(slice(100, 400), slice(700, 1100))):
    """Ruido de un píxel entre un cuadro y el siguiente (desvío robusto de la diferencia / raíz de 2),
    en una zona de la imagen, relativo al fondo y en cuentas. También la correlación del ruido entre
    píxeles vecinos (0 = cada píxel por su cuenta)."""
    rel, cuentas, vecinos = [], [], []
    for c in cuadros:
        a, b = cuadro_verde(pags[c])[zona], cuadro_verde(pags[c + 1])[zona]
        d = (b - a) / np.sqrt(2)
        mad = lambda v: 1.4826 * np.median(np.abs(v - np.median(v)))
        cuentas.append(mad(d))
        rel.append(mad(d / fondo[zona]))
        dd = d - d.mean()
        vecinos.append(float(np.mean(dd[:, 1:] * dd[:, :-1]) / np.mean(dd * dd)))
    return dict(relativo=float(np.mean(rel)), cuentas=float(np.mean(cuentas)), fondo_cuentas=float(fondo[zona].mean()),
                correlacion_vecinos=float(np.mean(vecinos)))


def poner(vacio, plant, dx, dy):
    """Recorte vacío + la partícula de la plantilla con su centro corrido (dx, dy) px respecto del centro
    del recorte (interpolación con splines cúbicos)."""
    m = (plant.shape[0] - 1) // 2
    yy, xx = np.mgrid[-m:m + 1, -m:m + 1]
    borde = np.clip((m - np.hypot(xx, yy)) / 6.0, 0, 1)                    # la plantilla se apaga suave en el borde
    p = ndimage.shift((plant - 1.0) * borde, (dy, dx), order=3, mode='constant', cval=0.0)
    img = vacio.copy()
    c = MEDIO
    img[c - m:c + m + 1, c - m:c + m + 1] += p
    return img


def localizar(img):
    """Centro de la partícula más cercana al centro del recorte, con los pasos de la conversión
    (Fourier → realce → detección → simetría radial). Devuelve (dx, dy) px respecto del centro, o NaN."""
    pos = tr.detectar(img, np.ones_like(img))
    if len(pos) == 0:
        return np.nan, np.nan
    r = np.hypot(pos[:, 0] - MEDIO, pos[:, 1] - MEDIO)
    k = int(np.argmin(r))
    return (pos[k, 0] - MEDIO, pos[k, 1] - MEDIO) if r[k] < 12 else (np.nan, np.nan)


def prueba_corrimientos(plant, vacios, corrimientos, eje='x'):
    """Para cada corrimiento conocido (px) pone la plantilla corrida sobre cada recorte vacío y la
    localiza. Devuelve tabla: corrimiento_px, medido_px (media), error_px (media − corrimiento, con el
    corrimiento 0 como referencia), desvio_px (dispersión entre recortes), n."""
    base = np.array([localizar(poner(v, plant, 0.0, 0.0)) for v in vacios])
    k = 0 if eje == 'x' else 1
    filas = []
    for s in corrimientos:
        med = np.array([localizar(poner(v, plant, s if k == 0 else 0.0, s if k == 1 else 0.0)) for v in vacios])
        dif = med[:, k] - base[:, k]                       # respecto de la misma realización sin correr: saca el sesgo de la plantilla
        ok = np.isfinite(dif)
        filas.append(dict(corrimiento_px=s, medido_px=float(np.mean(dif[ok])), error_px=float(np.mean(dif[ok]) - s),
                          desvio_px=float(np.std(med[ok, k] - np.mean(med[ok, k]))), n=int(ok.sum())))
    return pd.DataFrame(filas), base


def cota_teorica(plant, ruido_rel):
    """Lo mejor que se puede medir el centro de una imagen conocida con ruido blanco de desvío
    `ruido_rel` por píxel (cota de Cramér-Rao): sigma = ruido / raíz(suma de (dI/dx)²). En px."""
    gy, gx = np.gradient(plant)
    return float(ruido_rel / np.sqrt(np.sum(gx ** 2))), float(ruido_rel / np.sqrt(np.sum(gy ** 2)))


# ---------------------------------------------------------------- 4. forma del contorno
def centro_de_simetria(plant):
    """Centro de la plantilla por simetría radial (el mismo algoritmo), en px del recorte."""
    return tr.b7.simetria_radial(plant, 3.0)


def contorno(plant, n_angulos=36, r_max=32.0, r_desde=8.0):
    """Radio del borde de la partícula en función del ángulo. En cada sector angular se toma el perfil
    radial (interpolado), se busca el anillo oscuro (mínimo entre r_desde y r_max) y, hacia afuera, el
    radio donde el perfil recupera la mitad de la profundidad: la misma definición que `medir_radios`.
    Devuelve dict: angulo [rad], r_px [n_angulos], r0 (radio medio), y para m = 1..4 la amplitud a_m
    [px] y la fase del armónico m de r(θ):  r(θ) = r0 + Σ a_m cos(m(θ − fase_m)).
      a_1: error al elegir el centro;  a_2: elipse (cociente de ejes ≈ 1 + 2·a_2/r0);  a_3, a_4: otras."""
    cx, cy = centro_de_simetria(plant)
    ang = np.arange(n_angulos) * 2 * np.pi / n_angulos
    rr = np.arange(0.0, r_max + 0.25, 0.25)
    medio_sector = np.pi / n_angulos
    r_borde = np.full(n_angulos, np.nan)
    for k, a in enumerate(ang):
        perf = np.zeros(len(rr))
        for da in np.linspace(-medio_sector, medio_sector, 7):              # promedio dentro del sector
            xs, ys = cx + rr * np.cos(a + da), cy + rr * np.sin(a + da)
            perf += ndimage.map_coordinates(plant, [ys, xs], order=1, mode='nearest')
        perf = perf / 7 - 1.0
        s = (rr >= r_desde)
        i = np.flatnonzero(s)[int(np.argmin(perf[s]))]
        prof = perf[i]
        for j in range(i, len(rr) - 1):
            if perf[j + 1] > prof / 2:
                r_borde[k] = rr[j] + (prof / 2 - perf[j]) / (perf[j + 1] - perf[j]) * 0.25
                break
    ok = np.isfinite(r_borde)
    out = dict(angulo=ang, r_px=r_borde, r0=float(np.nanmean(r_borde)), centro=(cx, cy), validos=int(ok.sum()))
    if ok.all():
        F = np.fft.rfft(r_borde) / n_angulos
        for m_ in (1, 2, 3, 4):
            out['a%d' % m_] = float(2 * np.abs(F[m_]))
            out['fase%d' % m_] = float(-np.angle(F[m_]) / m_)
    return out


def plantilla_redonda(plant):
    """La misma partícula hecha perfectamente redonda: promedio en ángulo de su perfil radial."""
    cx, cy = centro_de_simetria(plant)
    m = (plant.shape[0] - 1) // 2
    yy, xx = np.mgrid[0:2 * m + 1, 0:2 * m + 1]
    r = np.hypot(xx - cx, yy - cy)
    rr = np.arange(0, m * 1.5, 0.5)
    perf = np.array([plant[(r >= a) & (r < a + 0.5)].mean() if ((r >= a) & (r < a + 0.5)).any() else 1.0 for a in rr])
    return np.interp(r, rr + 0.25, perf)


def pilas_vacias(nombre, pags, fondo, cuadro_ini, cuadro_fin, cuantos=12, n_cuadros=40, lejos_px=150.0, semilla=0):
    """Para separar el error que cambia de un cuadro al siguiente del que queda fijo en cada lugar:
    `cuantos` lugares sin partículas seguidas durante todo el tramo [cuadro_ini, cuadro_fin] (platina
    quieta) y, en cada uno, `n_cuadros` cuadros repartidos en el tramo. Devuelve [cuantos, n, L, L]
    (cuadros divididos por el fondo) y los centros [(x, y)]."""
    d = cargar(nombre)
    t = d['tr']
    t = t[t.x_px.notna() & (t.cuadro >= cuadro_ini) & (t.cuadro <= cuadro_fin)]
    pts = t[['x_px', 'y_px']].values
    tx, ty = d['trampa_px']
    rng = np.random.default_rng(semilla)
    centros, intentos = [], 0
    while len(centros) < cuantos and intentos < 5000:
        intentos += 1
        x = float(rng.uniform(MEDIO + 5, 1280 - MEDIO - 5))
        y = float(rng.uniform(MEDIO + 5, 1024 - MEDIO - 5))
        if np.hypot(x - tx, y - ty) < lejos_px + 60:
            continue
        if len(pts) and np.min(np.hypot(pts[:, 0] - x, pts[:, 1] - y)) < lejos_px:
            continue
        if any(np.hypot(x - a, y - b) < 60 for a, b in centros):
            continue
        centros.append((x, y))
    cuadros = np.unique(np.linspace(cuadro_ini, cuadro_fin, n_cuadros).round().astype(int))
    pila = np.zeros((len(centros), len(cuadros), 2 * MEDIO + 1, 2 * MEDIO + 1), np.float32)
    for k, c in enumerate(cuadros):
        plano = cuadro_verde(pags[c]) / fondo
        for i, (x, y) in enumerate(centros):
            pila[i, k] = recorte(plano, x, y, MEDIO)
    return pila, centros


def separar_errores(plant, pila):
    """Pone la plantilla (sin correr) en cada cuadro de cada lugar y la localiza. Devuelve dict:
      temporal_x_px, temporal_y_px   desvío entre cuadros en un mismo lugar (mediana de los lugares):
                                     el error que cambia de un cuadro al siguiente
      fijo_x_px, fijo_y_px           desvío entre lugares de la posición media: el corrimiento que
                                     mete el fondo de cada lugar, y que no cambia con el tiempo
      posiciones [lugares, cuadros, 2]"""
    K, M = pila.shape[:2]
    pos = np.full((K, M, 2), np.nan)
    for i in range(K):
        for k in range(M):
            pos[i, k] = localizar(poner(pila[i, k], plant, 0.0, 0.0))
    return dict(temporal_x_px=float(np.nanmedian(np.nanstd(pos[:, :, 0], axis=1))), temporal_y_px=float(np.nanmedian(np.nanstd(pos[:, :, 1], axis=1))),
                fijo_x_px=float(np.nanstd(np.nanmean(pos[:, :, 0], axis=1))), fijo_y_px=float(np.nanstd(np.nanmean(pos[:, :, 1], axis=1))),
                posiciones=pos, perdidas=int(np.isnan(pos[:, :, 0]).sum()))


def localizar_por_plantilla(img, plant, buscar=6):
    """Otra forma de medir el centro, para comparar: correlación de la imagen con la plantilla de la
    propia partícula (filtro adaptado) y máximo con precisión de subpíxel (parábola). Con ruido blanco
    es casi lo mejor posible, pero necesita conocer de antemano cómo se ve la partícula.
    Devuelve (dx, dy) px respecto del centro del recorte."""
    from scipy.signal import fftconvolve
    m = (plant.shape[0] - 1) // 2
    yy, xx = np.mgrid[-m:m + 1, -m:m + 1]
    borde = np.clip((m - np.hypot(xx, yy)) / 6.0, 0, 1)
    nucleo = (plant - 1.0) * borde
    nucleo = nucleo - nucleo.mean()
    c = MEDIO
    zona = img[c - m - buscar:c + m + buscar + 1, c - m - buscar:c + m + buscar + 1] - 1.0
    cc = fftconvolve(zona, nucleo[::-1, ::-1], mode='valid')                # [2·buscar+1]²
    i, j = np.unravel_index(np.argmax(cc), cc.shape)
    if not (0 < i < cc.shape[0] - 1 and 0 < j < cc.shape[1] - 1):
        return np.nan, np.nan
    sub = lambda a, b, d: 0.5 * (a - d) / (a - 2 * b + d)
    return (j - buscar + sub(cc[i, j - 1], cc[i, j], cc[i, j + 1]), i - buscar + sub(cc[i - 1, j], cc[i, j], cc[i + 1, j]))
