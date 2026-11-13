"""De cuadro crudo a imagen de trabajo: canal Bayer o demosaico, fondo, recorte; vista; escala.

Coordenadas (CONTRATO.md): x = columna, y = fila, el centro del píxel (i, j) está en (x=j, y=i).

Bayer. `bayer_fase` (ej. 'RGGB') nombra el color de los píxeles (0,0), (0,1), (1,0), (1,1) del
cuadro GUARDADO (si la cámara recorta un ROI que empieza en fila o columna impar, la fase del
archivo cambia: el grabador tiene que guardar la del archivo). Extraer un canal da un "plano"
de media resolución: el píxel (v, u) del plano rojo es el píxel (2v + di, 2u + dj) del sensor,
donde (di, dj) es la posición de R en la celda 2×2. Por eso, una posición (x, y) medida en el
plano corresponde en el sensor a

    x_sensor = 2·x + dj,   y_sensor = 2·y + di                    (rojo, azul)
    x_sensor = 2·x + 1/2,  y_sensor = 2·y + 1/2                   (verde)

El verde promedia los dos G de la celda, que están en diagonal: (0,1) y (1,0) para RGGB/BGGR,
(0,0) y (1,1) para GRBG/GBRG; en las cuatro fases el centro del promedio es (1/2, 1/2). El
promedio de dos muestras a ±(1/2, -1/2) (o ±(1/2, 1/2)) es una convolución con un núcleo
simétrico: no corre el centro de un objeto simétrico (una bead), solo lo suaviza un poco.
`a_coordenadas_sensor` hace esta cuenta. La escala `escala_m_por_px` de grabacion.json es por
píxel del SENSOR: las posiciones medidas en un plano se pasan al sensor ANTES de multiplicar.

RGB (formato_pixel = 'rgb'): cuadros [alto, ancho, 3] ya procesados por otro programa (por
ejemplo ThorCam con Image Type sRGB / Linear sRGB, o un AVI; ver pinzas.importar). Cada canal
está a resolución COMPLETA (interpolado por el demosaico del programa): 'verde' = el plano G,
'rojo' = R, 'azul' = B, 'gris' = (R + G + B) / 3; px_sensor_por_px = 1 y a_coordenadas_sensor
es la identidad (más el recorte).
"""
import numpy as np

_POS = [(0, 0), (0, 1), (1, 0), (1, 1)]
_LETRA = {'rojo': 'R', 'verde': 'G', 'azul': 'B'}


# ---------------------------------------------------------------- metadatos y geometría
def _camara(meta):
    """Acepta el meta completo de grabacion.json o directamente su dict 'camara'."""
    meta = meta or {}
    return meta.get('camara', meta)


def formato_pixel(meta):
    """'bayer' | 'mono' | 'rgb' (camara.formato_pixel de grabacion.json; sin dato = 'mono')."""
    return str(_camara(meta).get('formato_pixel') or 'mono').lower()


def es_bayer(meta):
    return formato_pixel(meta) == 'bayer'


def es_rgb(meta):
    return formato_pixel(meta) == 'rgb'


def fase_bayer(meta):
    f = _camara(meta).get('bayer_fase') or 'RGGB'
    return str(f).upper()


def posiciones_bayer(letra, fase):
    """Posiciones (di, dj) de la celda 2×2 donde está el color `letra` ('R', 'G', 'B')."""
    fase = fase.upper()
    if sorted(fase) != ['B', 'G', 'G', 'R']:
        raise ValueError(f'fase Bayer no válida: {fase}')
    return [p for p, l in zip(_POS, fase) if l == letra]


def px_sensor_por_px(canal, meta):
    """Píxeles del sensor por píxel de la imagen que devuelve `preparar`: 2 en un plano Bayer,
    si no 1 (mono, RGB ya demosaicado, 'gris', 'crudo')."""
    return 2.0 if es_bayer(meta) and canal in _LETRA else 1.0


def recorte_efectivo(recorte, forma, bayer=False):
    """Recorte (x, y, ancho, alto) en px del sensor, recortado a la imagen y, si es Bayer,
    con x, y, ancho y alto llevados a pares (hacia abajo) para no cambiar la fase.
    None → imagen completa (también llevada a pares si es Bayer). Devuelve (x, y, ancho, alto)."""
    alto_img, ancho_img = forma[:2]
    x, y, a, h = (0, 0, ancho_img, alto_img) if recorte is None else [int(v) for v in recorte]
    x, y = max(0, x), max(0, y)
    a, h = min(a, ancho_img - x), min(h, alto_img - y)
    if bayer:
        x, y = x - x % 2, y - y % 2
        a, h = min(a, ancho_img - x), min(h, alto_img - y)
        a, h = a - a % 2, h - h % 2
    return x, y, max(a, 0), max(h, 0)


def a_coordenadas_sensor(x, y, canal, fase, recorte=None):
    """Posición (x, y) en la imagen de `preparar` → píxeles del sensor (del cuadro guardado).

    canal: 'rojo' | 'verde' | 'azul' (plano Bayer, ver fórmulas del módulo) o 'gris' | 'crudo'
           | mono (identidad). fase: 'RGGB', ... o None (mono → identidad).
    recorte: el (x0, y0, ancho, alto) que se pasó a `preparar` (se suma el recorte EFECTIVO, ya
           llevado a pares). Acepta escalares o arrays. Devuelve (x_sensor, y_sensor).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if fase is not None and canal in _LETRA:
        pos = posiciones_bayer(_LETRA[canal], fase)
        di = np.mean([p[0] for p in pos])
        dj = np.mean([p[1] for p in pos])
        x, y = 2 * x + dj, 2 * y + di
    if recorte is not None:
        x0, y0 = int(recorte[0]), int(recorte[1])
        if fase is not None:
            x0, y0 = x0 - x0 % 2, y0 - y0 % 2
        x, y = x + max(x0, 0), y + max(y0, 0)
    return x, y


def a_coordenadas_canal(x_sensor, y_sensor, canal, fase, recorte=None):
    """Inversa de `a_coordenadas_sensor`: px del sensor → posición en la imagen de `preparar`."""
    x = np.asarray(x_sensor, float)
    y = np.asarray(y_sensor, float)
    if recorte is not None:
        x0, y0 = int(recorte[0]), int(recorte[1])
        if fase is not None:
            x0, y0 = x0 - x0 % 2, y0 - y0 % 2
        x, y = x - max(x0, 0), y - max(y0, 0)
    if fase is not None and canal in _LETRA:
        pos = posiciones_bayer(_LETRA[canal], fase)
        di = np.mean([p[0] for p in pos])
        dj = np.mean([p[1] for p in pos])
        x, y = (x - dj) / 2, (y - di) / 2
    return x, y


# ---------------------------------------------------------------- canales
def extraer_canal(crudo, canal, fase):
    """Plano de media resolución del canal ('rojo', 'verde', 'azul') de un mosaico Bayer → float32.

    Verde = promedio de los dos G de cada celda 2×2. Si el alto o el ancho son impares se
    descarta la última fila/columna.
    """
    c = np.asarray(crudo)
    c = c[: c.shape[0] // 2 * 2, : c.shape[1] // 2 * 2]
    planos = [c[di::2, dj::2].astype(np.float32) for di, dj in posiciones_bayer(_LETRA[canal], fase)]
    return planos[0] if len(planos) == 1 else 0.5 * (planos[0] + planos[1])


def demosaico_bilineal(crudo, fase):
    """Demosaico bilineal → float32 [alto, ancho, 3] (R, G, B).

    Convolución normalizada por canal: C = (K * (m·I)) / (K * m), con m la máscara del canal y
    K = [[1,2,1],[2,4,2],[1,2,1]] para R y B, K = [[0,1,0],[1,4,1],[0,1,0]] para G. En el píxel
    del propio color devuelve el valor medido; en los demás, el promedio de los vecinos más
    cercanos de ese color. En los bordes usa solo los vecinos que existen (no inventa).
    """
    from scipy.ndimage import convolve
    img = np.asarray(crudo, dtype=np.float64)
    alto, ancho = img.shape
    kRB = np.array([[1, 2, 1], [2, 4, 2], [1, 2, 1]], float)
    kG = np.array([[0, 1, 0], [1, 4, 1], [0, 1, 0]], float)
    out = np.empty((alto, ancho, 3), np.float32)
    for k, (letra, K) in enumerate([('R', kRB), ('G', kG), ('B', kRB)]):
        m = np.zeros((alto, ancho))
        for di, dj in posiciones_bayer(letra, fase):
            m[di::2, dj::2] = 1.0
        num = convolve(img * m, K, mode='constant', cval=0.0)
        den = convolve(m, K, mode='constant', cval=0.0)
        out[..., k] = num / den
    return out


def preparar(crudo, meta, canal='verde', fondo=None, recorte=None, modo_fondo='dividir'):
    """Cuadro crudo uint16 → imagen de trabajo float32.

    crudo : [alto, ancho] uint16 del sensor (mosaico Bayer si formato_pixel == 'bayer').
    meta  : grabacion.json (o su 'camara'): usa formato_pixel y bayer_fase.
    canal : 'verde' | 'rojo' | 'azul' → plano Bayer de media resolución (ver módulo);
            'gris' → demosaico bilineal y promedio (R + G + B) / 3, resolución completa;
            'crudo' → el mosaico tal cual, en float. En mono se ignora (imagen tal cual).
    recorte : (x, y, ancho, alto) en px del SENSOR, se aplica antes de extraer el canal; en
            Bayer se lleva a pares (recorte_efectivo) para no cambiar la fase.
    fondo : imagen float de la misma forma que la salida (ver fondo_mediana) o None.
    modo_fondo : 'dividir' → I / fondo (corrección de campo plano: iluminación despareja
            multiplicativa; fondo <= 0 → NaN) | 'restar' → I - fondo.
            Ojo: si la bead casi no se mueve (atrapada, pegada), la mediana de cuadros la
            contiene y dividir/restar la borra. Para esas, un fondo con la bead fuera de campo.
    RGB [alto, ancho, 3] (formato_pixel == 'rgb'): 'verde'/'rojo'/'azul' = ese plano a resolución
            completa; 'gris' y 'crudo' = (R + G + B) / 3.
    Devuelve float32 [alto', ancho'].
    """
    bayer = es_bayer(meta)
    x, y, a, h = recorte_efectivo(recorte, np.shape(crudo), bayer)
    img = np.asarray(crudo)[y:y + h, x:x + a]
    if img.ndim == 3:
        out = _canal_rgb(img, canal)
    elif canal == 'crudo' or not bayer:
        out = img.astype(np.float32)
    elif canal == 'gris':
        out = demosaico_bilineal(img, fase_bayer(meta)).mean(axis=2, dtype=np.float32)
    elif canal in _LETRA:
        out = extraer_canal(img, canal, fase_bayer(meta))
    else:
        raise ValueError(f"canal desconocido: {canal!r} (verde, rojo, azul, gris, crudo)")
    if fondo is not None:
        f = np.asarray(fondo, np.float32)
        if f.shape != out.shape:
            raise ValueError(f'fondo {f.shape} y cuadro {out.shape} no coinciden')
        if modo_fondo == 'restar':
            out = out - f
        else:
            with np.errstate(divide='ignore', invalid='ignore'):
                out = np.where(f > 0, out / np.where(f > 0, f, 1), np.nan).astype(np.float32)
    return out


def _canal_rgb(img, canal):
    """Plano de un cuadro RGB [alto, ancho, 3 o 4] → float32 [alto, ancho] (ver `preparar`)."""
    if canal in _LETRA:
        return img[..., ('rojo', 'verde', 'azul').index(canal)].astype(np.float32)
    if canal in ('gris', 'crudo'):
        return img[..., :3].mean(axis=2, dtype=np.float32)
    raise ValueError(f"canal desconocido: {canal!r} (verde, rojo, azul, gris, crudo)")


def fondo_mediana(grabacion, indices=None, canal='verde', recorte=None, n_max=51):
    """Fondo = mediana píxel a píxel de cuadros preparados (sin fondo) de una grabación.

    grabacion : objeto io.Grabacion. indices: cuadros a usar; None = n_max repartidos en toda
    la grabación. Devuelve float32 con la forma de preparar(...).
    """
    n = len(grabacion)
    if indices is None:
        indices = np.unique(np.linspace(0, n - 1, min(n_max, n)).round().astype(int))
    pila = np.stack([preparar(grabacion.cuadro(i), grabacion.meta, canal, recorte=recorte)
                     for i in indices])
    return np.median(pila, axis=0).astype(np.float32)


# ---------------------------------------------------------------- vista
def _reducir(img, paso):
    """Promedio por bloques paso×paso (descarta el borde que no completa un bloque)."""
    if paso <= 1:
        return img.astype(np.float32)
    h, w = img.shape[0] // paso * paso, img.shape[1] // paso * paso
    r = img[:h, :w].astype(np.float32).reshape(h // paso, paso, w // paso, paso, *img.shape[2:])
    return r.mean(axis=(1, 3))


def vista_rapida(crudo, meta, max_lado=512, estirar=True):
    """Imagen chica en uint8 SOLO para mostrar (nunca para analizar).

    Bayer → RGB [h, w, 3] armado con los planos R, G (promedio) y B; RGB → tal cual (sin alfa);
    mono → gris [h, w].
    Se reduce por promedio de bloques hasta que el lado mayor sea <= max_lado.
    estirar=True: escala lineal entre los percentiles 0.5 y 99.5; False: entre 0 y 2^bits - 1.
    """
    if es_bayer(meta):
        f = fase_bayer(meta)
        img = np.stack([extraer_canal(crudo, c, f) for c in ('rojo', 'verde', 'azul')], axis=-1)
    else:
        img = np.asarray(crudo)
        img = img[..., :3] if img.ndim == 3 else img
    paso = int(np.ceil(max(img.shape[:2]) / float(max_lado)))
    r = _reducir(img, paso)
    if estirar:
        lo, hi = np.percentile(r, [0.5, 99.5])
    else:
        lo, hi = 0.0, 2.0 ** int(_camara(meta).get('bits', 10)) - 1
    hi = hi if hi > lo else lo + 1
    return np.clip((r - lo) * (255.0 / (hi - lo)), 0, 255).astype(np.uint8)


# ---------------------------------------------------------------- escala con retícula
def _dtft(fw, k, yy, xx):
    """F(k) = sum_r fw(r) exp(-i 2π k·r), k = (kx, ky) en ciclos/px (separable: O(alto·ancho))."""
    ex = np.exp(-2j * np.pi * k[0] * xx)
    ey = np.exp(-2j * np.pi * k[1] * yy)
    return ey @ (fw @ ex)


def _refinar(fw, k0, yy, xx, paso_bin):
    """Máximo de |F(k)| cerca de k0 (Nelder-Mead): frecuencia subpíxel del pico."""
    from scipy.optimize import minimize
    r = minimize(lambda k: -abs(_dtft(fw, k, yy, xx)), np.asarray(k0, float), method='Nelder-Mead',
                 options=dict(xatol=1e-8, fatol=1e-12, maxiter=400,
                              initial_simplex=[k0, k0 + [paso_bin / 2, 0], k0 + [0, paso_bin / 2]]))
    return r.x, _dtft(fw, r.x, yy, xx)


def _ventaneada(img):
    f = np.asarray(img, float)
    f = np.where(np.isfinite(f), f, np.nanmean(f))
    f = f - f.mean()
    return f * np.hanning(f.shape[0])[:, None] * np.hanning(f.shape[1])[None, :]


def _picos(fw, kmin, n_periodos_min=3):
    """Los dos vectores de onda fundamentales de la retícula (pico grueso en la FFT con relleno ×2).

    k1 = máximo del espectro de potencia fuera de |k| < kmin; k2 = máximo entre las direcciones
    a 90° ± 30° de k1. Si a k/2 o k/3 hay un pico > 0.3 del de k, k era un armónico y se baja
    al subarmónico. Devuelve ([k1, k2], paso_bin) con k en ciclos/px (k2 puede ser None).
    """
    alto, ancho = fw.shape
    P = np.abs(np.fft.fft2(fw, s=(2 * alto, 2 * ancho))) ** 2
    fy = np.fft.fftfreq(2 * alto)[:, None]
    fx = np.fft.fftfreq(2 * ancho)[None, :]
    kmod = np.hypot(fx, fy)
    P = np.where((kmod >= kmin) & (fx >= 0), P, 0.0)       # medio plano: k y -k son el mismo pico

    def mejor(mask):
        Q = np.where(mask, P, 0.0)
        if Q.max() <= 0:
            return None
        i, j = np.unravel_index(np.argmax(Q), Q.shape)
        k = np.array([fx[0, j], fy[i, 0]])
        for _ in range(3):                                  # bajar a la fundamental
            for m in (3, 2):
                ks = k / m
                if np.hypot(*ks) < kmin:
                    continue
                d = np.hypot(fx - ks[0], fy - ks[1]) <= 1.5 / (2 * min(alto, ancho))
                if np.any(d) and np.max(np.where(d, P, 0)) > 0.3 * Q[i, j]:
                    i2, j2 = np.unravel_index(np.argmax(np.where(d, P, 0)), P.shape)
                    k = np.array([fx[0, j2], fy[i2, 0]])
                    break
            else:
                break
        return k

    k1 = mejor(np.ones_like(P, bool))
    ang = np.arctan2(fy, fx) - np.arctan2(k1[1], k1[0])
    perp = np.abs(np.cos(ang)) < np.sin(np.deg2rad(30))
    k2 = mejor(perp)
    if k2 is not None and np.sqrt(P[np.unravel_index(np.argmax(np.where(perp, P, 0)), P.shape)]) \
            < 0.05 * np.sqrt(P.max()):
        k2 = None                                           # no hay segunda familia de líneas
    return [k1, k2], 1.0 / (2 * max(alto, ancho))


def _segmentos(k, fase, forma):
    """Líneas k·r + fase/2π = m (m entero) recortadas a la imagen → [nl, 4] (x0, y0, x1, y1) px."""
    alto, ancho = forma
    k = np.asarray(k, float)
    kk = k @ k
    d = np.array([-k[1], k[0]]) / np.sqrt(kk)               # dirección de la línea
    esquinas = np.array([[0, 0], [ancho - 1, 0], [0, alto - 1], [ancho - 1, alto - 1]], float)
    proy = esquinas @ k + fase / (2 * np.pi)
    segs = []
    L = np.hypot(alto, ancho)
    for m in range(int(np.floor(proy.min())), int(np.ceil(proy.max())) + 1):
        r0 = (m - fase / (2 * np.pi)) * k / kk
        t0, t1 = -L, L                                      # Liang-Barsky contra [0,ancho-1]x[0,alto-1]
        for p, q in ((-d[0], r0[0]), (d[0], ancho - 1 - r0[0]), (-d[1], r0[1]), (d[1], alto - 1 - r0[1])):
            if abs(p) < 1e-12:
                if q < 0:
                    t0, t1 = 1, 0
                continue
            t = q / p
            if p < 0:
                t0 = max(t0, t)
            else:
                t1 = min(t1, t)
        if t1 > t0:
            segs.append([*(r0 + t0 * d), *(r0 + t1 * d)])
    return np.array(segs, float).reshape(-1, 4)


def escala_reticula(img, paso_m, px_sensor_por_px=1.0, n_periodos_min=3):
    """Escala [m/px del sensor] a partir de una imagen de retícula (líneas o puntos) de período paso_m.

    img : 2D float (por ejemplo preparar(cuadro, meta, 'verde')). paso_m: período de la retícula [m].
    px_sensor_por_px : píxeles del sensor por píxel de img (px_sensor_por_px(canal, meta): 2 en
          un plano Bayer). Así la escala devuelta es siempre por píxel del SENSOR.
    Método: se resta la media, se aplica una ventana de Hann 2D y se busca el pico fundamental
    de |FFT|² para cada familia de líneas (k1 y la perpendicular k2, en ciclos/px). Cada pico se
    refina maximizando |F(k)| = |sum f(r) w(r) e^{-i2π k·r}| en k continuo (estimador de máxima
    verosimilitud de la frecuencia de una sinusoide). Período P = 1/|k| (no depende del giro de
    la retícula); m_por_px = paso_m / (P · px_sensor_por_px).
    Error: se repite el refinamiento en los cuatro cuadrantes; err_P = desvío de los cuatro / 2.
    Es una cota conservadora (cada cuadrante tiene la mitad de períodos por lado).
    Devuelve dict: m_por_px, err_m_por_px, fuente='reticula', periodo_px [2] (x, y; en px de
    img), err_periodo_px [2], m_por_px_xy [2], diferencia_xy_rel (P_y/P_x - 1), angulo_rad
    (giro de la familia x), k [2, 2] ciclos/px, fase_rad [2], lineas_x, lineas_y ([nl, 4]
    segmentos (x0, y0, x1, y1) en px de img, sobre las líneas detectadas: para dibujar),
    lineas_oscuras (bool). Lo que no se pudo medir queda NaN.
    """
    f = np.asarray(img, float)
    alto, ancho = f.shape
    yy, xx = np.arange(alto, dtype=float), np.arange(ancho, dtype=float)
    fw = _ventaneada(f)
    kmin = n_periodos_min / min(alto, ancho)
    ks, paso_bin = _picos(fw, kmin)
    oscuras = bool(np.mean((f - f.mean()) ** 3) < 0)        # líneas oscuras sobre fondo claro

    fams = []
    for k0 in ks:
        if k0 is None:
            fams.append(None)
            continue
        k, F = _refinar(fw, k0, yy, xx, paso_bin)
        k = k if k[0] >= 0 else -k                          # representante con kx >= 0
        F = _dtft(fw, k, yy, xx)
        fase = float(np.angle(F)) + (np.pi if oscuras else 0.0)
        Pq = []
        for i0, j0 in ((0, 0), (0, ancho // 2), (alto // 2, 0), (alto // 2, ancho // 2)):
            sub = f[i0:i0 + alto // 2, j0:j0 + ancho // 2]
            kq, _ = _refinar(_ventaneada(sub), k, yy[:alto // 2], xx[:ancho // 2], 2 * paso_bin)
            Pq.append(1.0 / np.hypot(*kq))
        fams.append(dict(k=k, P=1.0 / np.hypot(*k), eP=float(np.std(Pq, ddof=1) / 2), fase=fase))

    vivas = [fm for fm in fams if fm is not None]
    vivas.sort(key=lambda fm: -abs(fm['k'][0]) / np.hypot(*fm['k']))   # primero la familia "x"
    fx = vivas[0] if vivas else None
    fy = vivas[1] if len(vivas) > 1 else None
    nan2 = np.array([np.nan, np.nan])

    def g(fm, clave, defecto=np.nan):
        return fm[clave] if fm is not None else defecto

    P = np.array([g(fx, 'P'), g(fy, 'P')])
    eP = np.array([g(fx, 'eP'), g(fy, 'eP')])
    m_xy = paso_m / (P * px_sensor_por_px)
    e_xy = m_xy * eP / P
    ok = np.isfinite(m_xy)
    m = float(np.mean(m_xy[ok])) if ok.any() else np.nan
    e = float(np.sqrt(np.sum(e_xy[ok] ** 2)) / ok.sum()) if ok.any() else np.nan
    return dict(m_por_px=m, err_m_por_px=e, fuente='reticula', paso_m=float(paso_m),
                periodo_px=P, err_periodo_px=eP, m_por_px_xy=m_xy, err_m_por_px_xy=e_xy,
                diferencia_xy_rel=float(P[1] / P[0] - 1) if ok.all() else np.nan,
                angulo_rad=float(np.arctan2(fx['k'][1], fx['k'][0])) if fx else np.nan,
                k=np.array([g(fx, 'k', nan2), g(fy, 'k', nan2)]),
                fase_rad=np.array([g(fx, 'fase'), g(fy, 'fase')]),
                lineas_x=_segmentos(fx['k'], fx['fase'], f.shape) if fx else np.zeros((0, 4)),
                lineas_y=_segmentos(fy['k'], fy['fase'], f.shape) if fy else np.zeros((0, 4)),
                lineas_oscuras=oscuras, px_sensor_por_px=float(px_sensor_por_px))
