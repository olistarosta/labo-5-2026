"""Tracking subpíxel de beads de campo claro: detección, localización, seguimiento y deriva.

Convenciones (CONTRATO.md): x = columna, y = fila, el centro del píxel (i, j) está en
(x=j, y=i). Todo en px salvo `px_a_m`. Lo no confiable es NaN: no se interpola.

Flujo:  detectar (candidatos en un cuadro)  ->  seguir (posición subpíxel por cuadro con
`simetria_radial` o `gauss2d`)  ->  quitar_deriva  ->  px_a_m.
"""
import time
import warnings

import numpy as np
from scipy import ndimage
from scipy.optimize import least_squares

METODOS = ('simetria_radial', 'gauss2d')


# ---------------------------------------------------------------- detección
def realzar(img, radio_px):
    """Imagen filtrada en la que cada bead de campo claro es un máximo positivo.

    img: float 2D (cuentas, o fracción del fondo; cualquier escala). radio_px: radio de la bead.

    Por qué así: una bead de campo claro tiene el CENTRO OSCURO (aprox. gaussiano, ancho
    ~0.6 R) y un anillo claro; el fondo es disparejo pero lento. Un filtro de banda
    diferencia-de-gaussianas sobre la imagen INVERTIDA, con s1 = 0.6 R (el ancho del núcleo,
    filtro adaptado) y s2 = 2.5 s1 (quita el fondo lento), deja un pico positivo en el centro
    de cada bead y ~0 en el fondo. Es la aproximación barata del Laplaciano de gaussiana
    (LoG) a la escala de la bead. Antes se divide por un fondo muy suavizado (s = 4 R) para
    que el contraste no dependa de la iluminación local.
    Devuelve float 2D del mismo tamaño (adimensional).
    """
    img = np.asarray(img, float)
    fondo = ndimage.gaussian_filter(img, 4.0 * radio_px, mode='reflect')
    escala = np.nanmedian(np.abs(fondo))
    if not np.isfinite(escala) or escala == 0:
        escala = 1.0
    plano = -(img - fondo) / (np.abs(fondo) + 1e-3 * escala)    # invertida: centro oscuro -> pico
    s1 = max(0.6 * radio_px, 0.8)
    return ndimage.gaussian_filter(plano, s1) - ndimage.gaussian_filter(plano, 2.5 * s1)


def detectar(img, radio_px, umbral=6.0, max_beads=50, separacion_min_px=None):
    """Candidatos a bead en un cuadro de campo claro.

    img        : float 2D (canal elegido; con o sin fondo restado: el filtro lo tolera).
    radio_px   : radio nominal de la bead [px].
    umbral     : el pico filtrado tiene que superar `umbral` veces el ruido robusto de la imagen
                 filtrada (1.4826 · MAD).
    max_beads  : se devuelven a lo sumo los `max_beads` más contrastados.
    separacion_min_px : distancia mínima entre candidatos (default 2 R).

    Máximos locales de `realzar` (ver ahí el porqué del filtro), descartando los que están a
    menos de 1.5 R del borde. Posición refinada al centro de masa del pico filtrado en una
    ventana de ±R/2 (precisión ~0.1 px: sirve para elegir, NO para medir; eso lo hace `seguir`).
    El radio se estima por la escala que maximiza el LoG normalizado (s² ∇²G) en {0.7, 1,
    1.4}·0.6 R, dividida por 0.6: es aproximado.

    Devuelve [k, 3] float (x, y, radio) en px, ordenado de mayor a menor contraste.
    """
    img = np.asarray(img, float)
    f = realzar(img, radio_px)
    ruido = 1.4826 * np.nanmedian(np.abs(f - np.nanmedian(f)))
    sep = separacion_min_px or 2.0 * radio_px
    tam = max(3, int(2 * np.floor(sep / 2) + 1))
    pico = (f == ndimage.maximum_filter(f, size=tam, mode='nearest')) & (f > umbral * ruido)
    borde = int(np.ceil(1.5 * radio_px))
    pico[:borde, :] = pico[-borde:, :] = False
    pico[:, :borde] = pico[:, -borde:] = False
    ii, jj = np.nonzero(pico)
    if ii.size == 0:
        return np.zeros((0, 3))
    orden = np.argsort(f[ii, jj])[::-1][:max_beads]
    ii, jj = ii[orden], jj[orden]
    # radio por la escala del LoG normalizado
    escalas = np.array([0.7, 1.0, 1.4]) * max(0.6 * radio_px, 0.8)
    plano = -(img - ndimage.gaussian_filter(img, 4.0 * radio_px))
    resp = np.array([-(s * s) * ndimage.gaussian_laplace(plano, s)[ii, jj] for s in escalas])
    radios = escalas[np.argmax(resp, axis=0)] / 0.6
    h = max(1, int(round(radio_px / 2)))
    out = []
    for i, j, r in zip(ii, jj, radios):
        sub = f[i - h:i + h + 1, j - h:j + h + 1]
        w = np.clip(sub - 0.5 * f[i, j], 0, None)
        di, dj = np.mgrid[-h:h + 1, -h:h + 1]
        sw = w.sum()
        out.append((j + (w * dj).sum() / sw, i + (w * di).sum() / sw, r))
    return np.array(out, float)


# ---------------------------------------------------------------- localizadores
def _recortar(img, x_px, y_px, lado_px):
    """Recorte cuadrado de lado impar centrado en el píxel más cercano a (x, y).

    Devuelve (recorte, j0, i0) con (j0, i0) la esquina superior izquierda en la imagen, o
    (None, j0, i0) si el recorte se sale de la imagen.
    """
    h = lado_px // 2
    jc, ic = int(np.floor(x_px + 0.5)), int(np.floor(y_px + 0.5))
    j0, i0 = jc - h, ic - h
    if j0 < 0 or i0 < 0 or j0 + lado_px > img.shape[1] or i0 + lado_px > img.shape[0]:
        return None, j0, i0
    return img[i0:i0 + lado_px, j0:j0 + lado_px], j0, i0


def simetria_radial(recorte, suavizado_px=0.0):
    """Centro de simetría radial de un recorte (Parthasarathy, Nat. Methods 9, 724, 2012).

    recorte     : float 2D [Ny, Nx] (Nx, Ny >= 4), la bead aproximadamente centrada. Sirve igual
                  para campo claro (centro oscuro + anillo) que para fluorescencia: solo usa que
                  el perfil es radialmente simétrico, no su signo.
    suavizado_px: 0 (default) = exactamente el algoritmo del paper. > 0 = antes se suaviza el
                  recorte con una gaussiana de ese ancho [px] (bordes 'nearest'). NO está en el
                  paper: lo agregamos porque con beads grandes (R >~ 8 px) el gradiente de la
                  bead por píxel es chico (~contraste/R) y el ruido de cada vértice domina; el
                  promedio 3×3 del paper no alcanza. Suavizar con una gaussiana preserva la
                  simetría radial, así que no corre el centro. Medido: s ~ 0.2 R deja el ruido
                  igual al de gauss2d (ver TouchDesigner/herramientas/comparar_tracking.py).
    Devuelve (x, y) en px en coordenadas del recorte (centro del píxel [i, j] en (x=j, y=i)),
    o (nan, nan) si el sistema es singular.

    Algoritmo (fiel al radialcenter.m del material suplementario del paper):
      1. Gradientes en la grilla rotada 45°, definidos en los vértices entre 4 píxeles:
         dI/du = I[i, j+1] - I[i+1, j],   dI/dv = I[i, j] - I[i+1, j+1].
      2. Suavizado de cada gradiente con un promedio 3×3 (conv2 'same', bordes con ceros).
      3. Cada vértice k define la recta que pasa por (x_k, y_k) con la pendiente del gradiente
         m_k = -(dv + du) / (du - dv)  (= g_y / g_x en la grilla original).
      4. El centro (xc, yc) minimiza Σ w_k d_k², con d_k la distancia del centro a la recta k
         y peso w_k = |∇I|²_k / dist(k, centroide de |∇I|²): solución cerrada de mínimos
         cuadrados (lsradialcenterfit).
    Detalles numéricos del original que se respetan: m = NaN -> m sin suavizar -> 0; m = ±inf
    -> 10·max(m finito). Única diferencia: si un vértice cae exactamente sobre el centroide
    (distancia 0, peso infinito en el original) se usa distancia 1e-6 px.
    """
    I = np.asarray(recorte, float)
    if suavizado_px and suavizado_px > 0:
        I = ndimage.gaussian_filter(I, suavizado_px, mode='nearest')
    Ny, Nx = I.shape
    xm = np.arange(Nx - 1) - (Nx - 2) / 2.0          # coordenadas de los vértices, centradas
    ym = np.arange(Ny - 1) - (Ny - 2) / 2.0
    xm, ym = np.meshgrid(xm, ym)
    dIdu = I[:-1, 1:] - I[1:, :-1]
    dIdv = I[:-1, :-1] - I[1:, 1:]
    fdu = ndimage.uniform_filter(dIdu, 3, mode='constant', cval=0.0)
    fdv = ndimage.uniform_filter(dIdv, 3, mode='constant', cval=0.0)
    dImag2 = fdu * fdu + fdv * fdv
    with np.errstate(divide='ignore', invalid='ignore'):
        m = -(fdv + fdu) / (fdu - fdv)
        sin_suav = (dIdv + dIdu) / (dIdu - dIdv)     # tal cual en el original (signo incluido)
    m = np.where(np.isnan(m), sin_suav, m)
    m = np.where(np.isnan(m), 0.0, m)
    inf = np.isinf(m)
    if inf.any():
        fin = m[~inf]
        m = np.where(inf, 10 * (fin.max() if fin.size else 1.0), m)
    b = ym - m * xm
    sdI2 = dImag2.sum()
    if not sdI2 > 0:
        return np.nan, np.nan
    xcen = (dImag2 * xm).sum() / sdI2
    ycen = (dImag2 * ym).sum() / sdI2
    dist = np.sqrt((xm - xcen) ** 2 + (ym - ycen) ** 2)
    w = dImag2 / np.maximum(dist, 1e-6)
    # lsradialcenterfit
    wm2p1 = w / (m * m + 1)
    sw = wm2p1.sum()
    smmw = (m * m * wm2p1).sum()
    smw = (m * wm2p1).sum()
    smbw = (m * b * wm2p1).sum()
    sbw = (b * wm2p1).sum()
    det = smw * smw - smmw * sw
    if det == 0 or not np.isfinite(det):
        return np.nan, np.nan
    xc = (smbw * sw - smw * sbw) / det
    yc = (smbw * smw - smmw * sbw) / det
    # de coordenadas centradas a coordenadas del recorte (0-based, centro del píxel = entero)
    return xc + (Nx - 1) / 2.0, yc + (Ny - 1) / 2.0


def _modelo_gauss(p, xx, yy):
    x0, y0, s, A, c = p
    e = np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * s * s))
    return A * e + c, e


def _residuo_gauss(p, xx, yy, z):
    return _modelo_gauss(p, xx, yy)[0] - z


def _jac_gauss(p, xx, yy, z):
    x0, y0, s, A, c = p
    _, e = _modelo_gauss(p, xx, yy)
    dx, dy = xx - x0, yy - y0
    Ae = A * e
    return np.column_stack([Ae * dx / (s * s), Ae * dy / (s * s),
                            Ae * (dx * dx + dy * dy) / s ** 3, e, np.ones_like(e)])


def gauss2d(recorte, sigma0_px=None, radio_ventana_px=None):
    """Ajuste gaussiano 2D con offset al NÚCLEO OSCURO de una bead de campo claro.

    recorte          : float 2D, la bead aproximadamente centrada.
    sigma0_px        : ancho inicial [px]; default lado/6.
    radio_ventana_px : solo entran al ajuste los píxeles a menos de esta distancia del centro
                       del recorte (default: todo el recorte).

    Modelo: z = A exp(-((x-x0)² + (y-y0)²) / 2 s²) + c ajustado a z = -recorte (imagen
    INVERTIDA: el centro oscuro pasa a ser un pico, A > 0). Por qué invertida y no el perfil
    con anillo: el núcleo es lo que tiene más señal y es casi gaussiano; el anillo es
    radialmente simétrico, así que no corre el centro, solo agranda s y baja c (es un error de
    modelo simétrico). Ajustar el perfil completo agrega parámetros (radio y ancho del anillo)
    que dependen de la óptica y del desenfoque. La integración sobre el área del píxel tampoco
    corre el centro (convolución simétrica). Mínimos cuadrados no pesados, Levenberg-Marquardt
    con jacobiano analítico (scipy.optimize.least_squares, method='lm').

    Devuelve (x, y, s, A, c): (x, y) px en coordenadas del recorte, s px, A y c en unidades de
    la imagen invertida. Todo NaN si no converge, si A <= 0, si s cae fuera de [0.3, lado] px
    o si el centro cae fuera del recorte.
    """
    z = -np.asarray(recorte, float)
    Ny, Nx = z.shape
    yy, xx = np.mgrid[0:Ny, 0:Nx].astype(float)
    xcen, ycen = (Nx - 1) / 2.0, (Ny - 1) / 2.0
    if radio_ventana_px is not None:
        dentro = (xx - xcen) ** 2 + (yy - ycen) ** 2 <= radio_ventana_px ** 2
        xx, yy, z = xx[dentro], yy[dentro], z[dentro]
    else:
        xx, yy, z = xx.ravel(), yy.ravel(), z.ravel()
    malo = (np.nan,) * 5
    if not np.all(np.isfinite(z)):
        return malo
    c0 = np.percentile(z, 20)
    w = np.clip(z - c0, 0, None)
    sw = w.sum()
    x0, y0 = ((w * xx).sum() / sw, (w * yy).sum() / sw) if sw > 0 else (xcen, ycen)
    s0 = sigma0_px or min(Nx, Ny) / 6.0
    p0 = [x0, y0, s0, max(z.max() - c0, 1e-12), c0]
    try:
        r = least_squares(_residuo_gauss, p0, jac=_jac_gauss, args=(xx, yy, z), method='lm',
                          xtol=1e-8, ftol=1e-10)
    except Exception:
        return malo
    x, y, s, A, c = r.x
    s = abs(s)
    if (not r.success or A <= 0 or not 0.3 <= s <= max(Nx, Ny)
            or not (-0.5 <= x <= Nx - 0.5 and -0.5 <= y <= Ny - 0.5)):
        return malo
    return x, y, s, A, c


def localizar(recorte, metodo='simetria_radial', suavizado_px=0.0):
    """(x, y) px en coordenadas del recorte con el método pedido (ver METODOS).

    suavizado_px solo lo usa 'simetria_radial' (ver ahí).
    """
    if metodo == 'simetria_radial':
        return simetria_radial(recorte, suavizado_px)
    if metodo == 'gauss2d':
        return gauss2d(recorte)[:2]
    raise ValueError(f'método desconocido: {metodo!r}; opciones {METODOS}')


# ---------------------------------------------------------------- seguimiento
def lado_recorte_sugerido(radio_px):
    """Lado impar del recorte: ~3 R (núcleo + anillo + margen), al menos 9 px."""
    return int(max(9, 2 * int(np.ceil(1.5 * radio_px)) + 1))


def contraste_z(recorte, x_px, y_px, radio_px):
    """Cuánto más oscuro que el borde del recorte es el núcleo de la bead, en sigmas.

    recorte: float 2D; (x_px, y_px): centro hallado, en coordenadas del recorte; radio_px: R.
    z = (mediana(borde) - media(núcleo)) / (sigma_borde / sqrt(n_núcleo)), con núcleo = píxeles
    a menos de R/2 del centro, borde = píxeles a más de 1.6 R (afuera del anillo) y sigma_borde =
    1.4826 · MAD del borde. Es el z-score del núcleo oscuro de campo claro: una bead real da
    decenas; un recorte de puro ruido, ~1. Si el borde tiene < 8 píxeles se usa el anillo
    exterior del recorte (última fila/columna). Devuelve float (NaN si no se puede calcular).
    """
    rec = np.asarray(recorte, float)
    yy, xx = np.indices(rec.shape)
    r = np.hypot(xx - x_px, yy - y_px)
    nucleo = rec[r <= 0.5 * radio_px]
    borde = rec[r >= 1.6 * radio_px]
    if borde.size < 8:
        borde = np.concatenate([rec[0], rec[-1], rec[1:-1, 0], rec[1:-1, -1]])
    if nucleo.size == 0:
        return np.nan
    med = np.median(borde)
    sig = 1.4826 * np.median(np.abs(borde - med))
    if not sig > 0:
        return np.nan
    return (med - nucleo.mean()) / (sig / np.sqrt(nucleo.size))


def seguir(obtener_cuadro, n, elegidas_px, metodo='simetria_radial', lado_recorte_px=None,
           max_salto_px=None, radio_px=None, suavizado_px=None, z_min=8.0, recentrar=True,
           progreso=None):
    """Posición subpíxel de cada bead elegida en cada cuadro.

    obtener_cuadro : función i -> float 2D, el cuadro i ya preparado (canal elegido, fondo
                     restado). Este módulo no sabe de canales: eso es de imagen.py.
    n              : cantidad de cuadros (i = 0 .. n-1).
    elegidas_px    : [nb, 2] (x, y) px, posición en el cuadro 0 (por ejemplo de `detectar`).
    metodo         : 'simetria_radial' | 'gauss2d'.
    lado_recorte_px: lado (impar) del recorte; default `lado_recorte_sugerido(radio_px)`.
                     Hace falta `radio_px` o `lado_recorte_px`.
    max_salto_px   : salto máximo entre la posición nueva y la última válida; más que eso = NaN
                     (default lado/4). La búsqueda sigue alrededor de la última válida.
    radio_px       : radio de la bead [px]; define los defaults de lado y suavizado.
    suavizado_px   : solo 'simetria_radial': suavizado previo (ver simetria_radial). None =
                     0.2·radio_px si se dio radio_px, si no 0 (paper puro).
    z_min          : contraste mínimo del núcleo oscuro (ver contraste_z) para aceptar la
                     posición; descarta recortes donde ya no hay bead (se fue, se desenfocó).
                     None = no se chequea. Necesita radio_px (si no, se usa lado/3).
                     Medido en sintéticos (R = 8 px): bead con señal alta z ~ 50, media ~ 25,
                     baja ~ 10; recorte de puro ruido z <~ 5. Default 8: con señal muy baja
                     bajarlo (y mirar contraste_z).
    recentrar      : si la posición cae a más de medio píxel del centro del recorte, se recorta
                     de nuevo centrado en ella y se vuelve a localizar (una sola vez).
    progreso       : callback opcional progreso(i, n), llamado después de cada cuadro.

    Una posición es NaN (valido=False) si el recorte sale del cuadro, si el localizador falla,
    si el resultado cae fuera del recorte, si salta más de max_salto_px o si el contraste z <
    z_min. No se interpola.

    Devuelve dict(x=[n, nb], y=[n, nb] px, valido=[n, nb] bool, contraste_z=[n, nb] (NaN donde no
    se llegó a calcular), metodo=str, lado_recorte_px=int, max_salto_px=float,
    suavizado_px=float, ms_por_bead_cuadro=float).
    """
    elegidas = np.atleast_2d(np.asarray(elegidas_px, float))[:, :2]
    nb = elegidas.shape[0]
    if lado_recorte_px is None:
        if radio_px is None:
            raise ValueError('hace falta radio_px o lado_recorte_px')
        lado_recorte_px = lado_recorte_sugerido(radio_px)
    lado = int(lado_recorte_px) | 1                         # impar
    salto = lado / 4.0 if max_salto_px is None else float(max_salto_px)
    if metodo not in METODOS:
        raise ValueError(f'método desconocido: {metodo!r}; opciones {METODOS}')
    if suavizado_px is None:
        suavizado_px = 0.2 * radio_px if radio_px is not None else 0.0
    if metodo != 'simetria_radial':
        suavizado_px = 0.0
    radio_z = radio_px if radio_px is not None else lado / 3.0
    X = np.full((n, nb), np.nan)
    Y = np.full((n, nb), np.nan)
    Z = np.full((n, nb), np.nan)
    ultima = elegidas.copy()
    t0 = time.perf_counter()
    for i in range(n):
        img = np.asarray(obtener_cuadro(i), float)
        for b in range(nb):
            xb, yb = _localizar_en(img, ultima[b, 0], ultima[b, 1], lado, metodo, recentrar,
                                   suavizado_px)
            if not (np.isfinite(xb) and np.hypot(xb - ultima[b, 0], yb - ultima[b, 1]) <= salto):
                continue
            if z_min is not None:
                rec, j0, i0 = _recortar(img, xb, yb, lado)
                if rec is None:
                    continue
                Z[i, b] = contraste_z(rec, xb - j0, yb - i0, radio_z)
                if not Z[i, b] >= z_min:
                    continue
            X[i, b], Y[i, b] = xb, yb
            ultima[b] = xb, yb
        if progreso is not None:
            progreso(i, n)
    ms = 1e3 * (time.perf_counter() - t0) / max(1, n * nb)
    return dict(x=X, y=Y, valido=np.isfinite(X), contraste_z=Z, metodo=metodo,
                lado_recorte_px=lado,
                max_salto_px=salto, suavizado_px=float(suavizado_px), ms_por_bead_cuadro=ms)


def _localizar_en(img, x_px, y_px, lado, metodo, recentrar=True, suavizado_px=0.0):
    """Posición (x, y) px en la imagen, buscando en un recorte alrededor de (x_px, y_px)."""
    for _ in range(2 if recentrar else 1):
        rec, j0, i0 = _recortar(img, x_px, y_px, lado)
        if rec is None:
            return np.nan, np.nan
        xr, yr = localizar(rec, metodo, suavizado_px)
        if not (np.isfinite(xr) and np.isfinite(yr)
                and -0.5 <= xr <= lado - 0.5 and -0.5 <= yr <= lado - 0.5):
            return np.nan, np.nan
        x_nuevo, y_nuevo = xr + j0, yr + i0
        c = (lado - 1) / 2.0
        if max(abs(xr - c), abs(yr - c)) <= 0.5:
            break
        x_px, y_px = x_nuevo, y_nuevo
    return x_nuevo, y_nuevo


# ---------------------------------------------------------------- deriva
def _deriva_media(z):
    """Deriva común [n] a partir de los desplazamientos cuadro a cuadro promediados entre beads.

    z: [n, nb] con NaN. El incremento entre dos cuadros consecutivos CON DATOS se promedia sobre
    las beads válidas en ambos; la deriva es la suma acumulada, 0 en el primer cuadro con datos.
    Cuadros sin ninguna bead válida -> NaN. Si entre dos cuadros con datos no hay ninguna bead
    en común, el incremento no se puede medir: se toma 0 y se avisa (warning).
    """
    n = z.shape[0]
    deriva = np.full(n, np.nan)
    con_datos = np.nonzero(np.isfinite(z).any(axis=1))[0]
    if con_datos.size == 0:
        return deriva
    inc = z[con_datos[1:]] - z[con_datos[:-1]]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        dm = np.nanmean(inc, axis=1)
    sin = ~np.isfinite(dm)
    if sin.any():
        warnings.warn(f'quitar_deriva: {sin.sum()} pasos sin ninguna bead en común; '
                      'incremento de deriva tomado como 0 ahí')
        dm[sin] = 0.0
    deriva[con_datos] = np.concatenate([[0.0], np.cumsum(dm)])
    return deriva


def _deriva_recta(z, t):
    """Recta por mínimos cuadrados en t de cada columna de z ([n, nb] con NaN) -> [n, nb]."""
    out = np.full(z.shape, np.nan)
    for b in range(z.shape[1]):
        ok = np.isfinite(z[:, b]) & np.isfinite(t)
        if ok.sum() < 2:
            continue
        p = np.polyfit(t[ok] - t[ok].mean(), z[ok, b], 1)
        out[:, b] = np.polyval(p, t - t[ok].mean())
    return out


def sigma_t_pasaaltos(fc_hz):
    """Ancho temporal [s] del núcleo gaussiano cuya transferencia vale 1/√2 en fc_hz.

    H(f) = exp(-2 π² σ² f²) = 1/√2  ->  σ = sqrt(ln 2) / (2 π fc).
    """
    return np.sqrt(np.log(2.0)) / (2 * np.pi * fc_hz)


def _deriva_lenta(z, t, fc_hz):
    """Suavizado gaussiano muy lento de cada columna de z con tiempos reales t (con huecos).

    Promedio pesado de Nadaraya-Watson con núcleo gaussiano de ancho σ = sigma_t_pasaaltos(fc):
    s(t) = Σ K(t - t_k) z_k / Σ K(t - t_k) sobre las muestras VÁLIDAS (NaN no pesa). Para que sea
    rápido se ubica cada muestra en una grilla uniforme de paso = mediana de dt (los huecos
    quedan vacíos, no se rellenan); el error de esa cuantización es < dt/2, despreciable frente a
    σ. Cerca de los extremos el promedio es asimétrico y sigue peor una deriva con pendiente.
    Devuelve [n, nb]; NaN donde no hay ninguna muestra válida dentro de ~4σ.
    """
    dt = np.nanmedian(np.diff(t))
    k = np.round((t - t[0]) / dt).astype(np.int64)
    m = int(k[-1]) + 1
    sig = sigma_t_pasaaltos(fc_hz) / dt
    out = np.full(z.shape, np.nan)
    for b in range(z.shape[1]):
        ok = np.isfinite(z[:, b])
        num = np.zeros(m)
        den = np.zeros(m)
        np.add.at(num, k[ok], z[ok, b])
        np.add.at(den, k[ok], 1.0)
        num = ndimage.gaussian_filter1d(num, sig, mode='constant', truncate=4.0)
        den = ndimage.gaussian_filter1d(den, sig, mode='constant', truncate=4.0)
        with np.errstate(invalid='ignore', divide='ignore'):
            s = num / den
        s[den < 1e-6 * den.max()] = np.nan
        out[:, b] = s[k]
    return out


def quitar_deriva(x, y, modo='media', fc_hz=None, t=None, referencia=None):
    """Estima y resta la deriva de la muestra.

    x, y      : [n, nb] posiciones (px o m: salen en las mismas unidades), NaN = perdido.
    modo      : 'media'    — deriva común = suma acumulada del desplazamiento cuadro a cuadro
                             promediado sobre las beads (ver _deriva_media). Para browniano libre
                             con varias beads (encargo, LIBRE).
                'recta'    — se resta un ajuste lineal en el tiempo (TRAMPA).
                'pasaaltos'— se resta un suavizado gaussiano muy lento con corte fc_hz (TRAMPA):
                             queda un pasa-altos. fc_hz tiene que estar MUY por debajo de f_c de
                             la trampa (p. ej. f_c/100): lo que pase por debajo de fc se pierde.
    fc_hz     : frecuencia de corte [Hz] para 'pasaaltos' (obligatoria ahí). Definición en
                sigma_t_pasaaltos: la transferencia del suavizado vale 1/√2 en fc.
    t         : [n] s, tiempos reales de cada cuadro (con huecos). Obligatorio en 'pasaaltos';
                en 'recta' si es None se usa el índice del cuadro.
    referencia: índices de beads de referencia (pegadas al vidrio, se mueven con la muestra).
                Si se da, la deriva sale SOLO de ellas (con el modo pedido aplicado a su
                promedio), es común a todas y se resta a todas las beads.

    Devuelve (x_corr, y_corr, deriva_x, deriva_y). deriva_* tiene forma [n] si es común
    ('media', o cualquier modo con referencia) y [n, nb] si es por bead ('recta' o
    'pasaaltos' sin referencia). x_corr = x - deriva (broadcast). Con 'recta' y 'pasaaltos' sin
    referencia la deriva incluye la posición media, así que x_corr queda centrado en 0 (posición
    relativa al equilibrio); con 'media' la deriva vale 0 en el primer cuadro con datos.

    OJO 'media' sin referencia: cada bead también está en el promedio, así que se le resta
    1/nb de su propio movimiento. Para beads independientes con el mismo D, el MSD corregido
    vale (1 - 1/nb) del verdadero (más la deriva). Ver informe: no se corrige acá.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if x.ndim == 1:
        x, y = x[:, None], y[:, None]
    n = x.shape[0]
    if modo not in ('media', 'recta', 'pasaaltos'):
        raise ValueError(f"modo desconocido: {modo!r}; opciones 'media', 'recta', 'pasaaltos'")
    if modo == 'pasaaltos':
        if fc_hz is None or t is None:
            raise ValueError("'pasaaltos' necesita fc_hz y t")
    tt = np.arange(n, dtype=float) if t is None else np.asarray(t, float)

    def una(z):
        if referencia is not None:
            ref = z[:, list(np.atleast_1d(referencia))]
            if modo == 'media':
                return _deriva_media(ref)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                # cada referencia relativa a su propia media, así un hueco en una no da un salto
                zr = np.nanmean(ref - np.nanmean(ref, axis=0), axis=1)[:, None]
            d =_deriva_recta(zr, tt) if modo == 'recta' else _deriva_lenta(zr, tt, fc_hz)
            return d[:, 0]
        if modo == 'media':
            return _deriva_media(z)
        if modo == 'recta':
            return _deriva_recta(z, tt)
        return _deriva_lenta(z, tt, fc_hz)

    dx, dy = una(x), una(y)
    bx = dx[:, None] if dx.ndim == 1 else dx
    by = dy[:, None] if dy.ndim == 1 else dy
    return x - bx, y - by, dx, dy


# ---------------------------------------------------------------- unidades
def px_a_m(x_px, m_por_px):
    """px -> m: x_m = x_px · m_por_px (m_por_px de ESCALA). NaN sigue siendo NaN."""
    return np.asarray(x_px, float) * float(m_por_px)
