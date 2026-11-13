"""Detectar en cada cuadro y ENLAZAR detecciones en trayectorias (beads libres, Stokes).

Para la bead atrapada alcanza `tracking.seguir` (una bead casi quieta, se la busca alrededor de
la última posición). Con beads libres o pegadas al vidrio que cruzan el campo con la platina
hay muchas, entran y salen, a veces no se ven unos cuadros y se cruzan. Entonces:

    detectar_todas   en CADA cuadro: candidatos (realce rápido) + posición subpíxel
                     (tracking.localizar en un recorte) + contraste (tracking.contraste_z)
    enlazar          cuadro a cuadro, asignación óptima detección ↔ trayectoria
                     (scipy.optimize.linear_sum_assignment), con memoria para huecos
    postprocesar     cortes por huecos largos, colisiones, interpolación de huecos cortos
                     (marcada), descarte de cortas, ids 0..P-1
    seguir_todas     las tres en fila
    desde_seguir     lleva el resultado de tracking.seguir al mismo formato y lo postprocesa

Formato LARGO: una fila por (partícula, cuadro), arrays [m] alineados: cuadro, particula, x, y
(px del plano de trabajo, x = columna, y = fila, CONTRATO.md), z (contraste), interpolado,
colision, hueco_previo. Lo no medido es NaN o no tiene fila; lo interpolado está marcado.
"""
import time
from collections import deque

import numpy as np
from scipy import ndimage
from scipy.optimize import linear_sum_assignment
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from . import tracking

# parámetro: (defecto, unidad, qué hace) — los muestra la interfaz (pasos.criterios_tracking)
CRITERIOS = {
    'max_salto_px': (5.0, 'px', 'distancia máxima entre un cuadro y el siguiente para que sea la misma partícula'),
    'memoria_cuadros': (3, 'cuadros', 'cuántos cuadros seguidos puede faltar una partícula sin cortar su trayectoria'),
    'interpolar_max_cuadros': (0, 'cuadros', 'huecos de hasta tantos cuadros se rellenan por interpolación lineal (marcados); 0 = nunca'),
    'largo_min_cuadros': (20, 'cuadros', 'trayectorias más cortas se descartan'),
    'dist_colision_px': (0.0, 'px', 'si dos partículas están más cerca que esto hay colisión; 0 = 2.5 × radio'),
    'accion_colision': ('marcar', '', 'marcar (se guarda igual, con colision = 1) · descartar (esos cuadros quedan NaN) · cortar (la trayectoria termina y sigue con otro id)'),
    'z_min': (8.0, '', 'contraste mínimo de una detección (tracking.contraste_z)'),
    'prediccion_cuadros': (10, 'cuadros', 'la posición esperada se extrapola con la velocidad media de los últimos tantos cuadros (Stokes, deriva); 0 = se busca alrededor de la última posición'),
}
ACCIONES = ('marcar', 'descartar', 'cortar')


def criterios_completos(criterios=None):
    """CRITERIOS con los valores de `criterios` encima (las claves que falten = defecto)."""
    c = {k: d[0] for k, d in CRITERIOS.items()}
    c.update({k: v for k, v in (criterios or {}).items() if v is not None})
    if c['accion_colision'] not in ACCIONES:
        raise ValueError(f"accion_colision desconocida: {c['accion_colision']!r}; opciones {ACCIONES}")
    return c


def _tiempos(t, n):
    """t [n] s con NaN → (tt [n] sin NaN, dt mediano). Sin t (o < 2 finitos) = índice del cuadro.

    Los NaN (tiempo no confiable) se rellenan interpolando en el índice: solo se usan para escalar
    la búsqueda y para interpolar huecos, nunca se exportan.
    """
    idx = np.arange(n, dtype=float)
    if t is None:
        return idx, 1.0
    t = np.asarray(t, float).ravel()
    if t.size < n:
        raise ValueError(f't tiene {t.size} tiempos y hay detecciones hasta el cuadro {n - 1}')
    t = t[:n]
    ok = np.isfinite(t)
    if ok.sum() < 2:
        return idx, 1.0
    tt = np.where(ok, t, np.interp(idx, idx[ok], t[ok]))
    d = np.diff(tt)
    d = d[d > 0]
    return tt, (float(np.median(d)) if d.size else 1.0)


# ---------------------------------------------------------------- detección en cada cuadro
def _binado(radio_px):
    """Factor de binado del realce: el núcleo (s1 = 0.6 R) queda con >= 1.5 px binados."""
    return max(1, int(max(0.6 * radio_px, 0.8) // 1.5))


def realzar_rapido(img, radio_px):
    """Lo mismo que tracking.realzar (filtro adaptado al núcleo oscuro), pero barato.

    img: float 2D (cuentas, cuadro dividido por el fondo, o con el fondo RESTADO: fondo ~ 0).
    Devuelve (f, b): f float32 [alto // b, ancho // b] (pico positivo en cada bead), b = binado.

    Por qué más barato: se bina b×b (promedio) con b = floor(0.6 R / 1.5), así el núcleo
    gaussiano de la bead (s1 = 0.6 R) sigue teniendo >= 1.5 px y el filtro adaptado no pierde
    señal, pero se filtran b² veces menos píxeles; y no se estima el radio por escala (LoG en 3
    escalas de tracking.detectar, lo más caro). El pico solo sirve para arrancar: la posición
    la da `tracking.localizar` a resolución completa.
    Diferencia con tracking.realzar: se divide por el fondo suavizado solo si es claramente
    positivo (min > 0.25 · mediana). Con un cuadro de fondo restado (~0) dividir amplifica el
    ruido donde el fondo cruza 0 y aparecen falsos picos.
    """
    img = np.asarray(img, np.float32)
    if not np.all(np.isfinite(img)):
        img = np.where(np.isfinite(img), img, np.float32(np.nanmedian(img)))
    b = _binado(radio_px)
    if b > 1:
        h, w = img.shape[0] // b, img.shape[1] // b
        img = img[:h * b, :w * b].reshape(h, b, w, b).mean(axis=(1, 3), dtype=np.float32)
    fondo = ndimage.gaussian_filter(img, 4.0 * radio_px / b, mode='reflect')
    med = float(np.median(fondo))
    plano = -(img - fondo)
    if med > 0 and float(fondo.min()) > 0.25 * med:
        plano /= fondo
    s1 = max(0.6 * radio_px, 0.8) / b
    f = ndimage.gaussian_filter(plano, s1) - ndimage.gaussian_filter(plano, 2.5 * s1)
    return f, b


def candidatos(img, radio_px, umbral=6.0, max_beads=200):
    """Posiciones gruesas (x, y) px de las beads de un cuadro: [k, 2], de mayor a menor pico.

    Máximos locales de `realzar_rapido` (vecindad 2 R) que superan `umbral` veces el ruido
    robusto del realce (1.4826 · MAD), como tracking.detectar. Subpíxel por centro de masa 3×3
    en la grilla binada (~b/4 px): alcanza para centrar el recorte de `localizar`.
    """
    f, b = realzar_rapido(img, radio_px)
    if min(f.shape) < 3:
        return np.zeros((0, 2))
    ruido = 1.4826 * float(np.median(np.abs(f - np.median(f))))
    tam = max(3, int(2 * np.floor(radio_px / b) + 1))
    pico = (f == ndimage.maximum_filter(f, size=tam, mode='nearest')) & (f > umbral * ruido)
    pico[0, :] = pico[-1, :] = pico[:, 0] = pico[:, -1] = False
    ii, jj = np.nonzero(pico)
    if ii.size == 0:
        return np.zeros((0, 2))
    orden = np.argsort(f[ii, jj])[::-1][:max_beads]
    ii, jj = ii[orden], jj[orden]
    di, dj = np.mgrid[-1:2, -1:2]
    sub = np.stack([f[ii + a, jj + c] for a, c in zip(di.ravel(), dj.ravel())], axis=1)
    w = np.clip(sub - 0.5 * f[ii, jj][:, None], 0, None)
    sw = w.sum(axis=1)
    xb = jj + (w * dj.ravel()).sum(axis=1) / sw
    yb = ii + (w * di.ravel()).sum(axis=1) / sw
    return np.stack([b * xb + (b - 1) / 2.0, b * yb + (b - 1) / 2.0], axis=1)


def detectar_todas(obtener_cuadro, n, radio_px, umbral=6.0, max_beads=200, metodo='simetria_radial',
                   lado_recorte_px=None, z_min=8.0, progreso=None):
    """Detecciones subpíxel en cada cuadro, en formato largo.

    obtener_cuadro : función i -> float 2D, cuadro i preparado (como para tracking.seguir).
    n              : cuadros (0 .. n-1).
    radio_px       : radio de la bead [px del plano de trabajo].
    umbral         : pico del realce >= umbral · ruido robusto (como tracking.detectar).
    max_beads      : a lo sumo tantos candidatos por cuadro (los más contrastados).
    metodo         : localizador subpíxel, 'simetria_radial' | 'gauss2d' (tracking.METODOS);
                     simetria_radial con suavizado 0.2 R, como tracking.seguir.
    lado_recorte_px: recorte del localizador; None = tracking.lado_recorte_sugerido(radio_px).
    z_min          : contraste mínimo (tracking.contraste_z) para aceptar la detección; None = no
                     se mira. Por qué: los máximos del realce que no son beads (ruido, polvo,
                     una bead que se desenfoca) dan z ~ 1-5; una bead real, decenas.
    progreso       : progreso(fraccion) opcional, después de cada cuadro.

    En cada cuadro: candidatos (realce binado, ver realzar_rapido) → tracking._localizar_en (con
    recentrado) → contraste_z en el recorte centrado. Se descarta la detección si el recorte sale
    del cuadro (bead entrando o saliendo: medio recorte afuera sesga la posición), si el
    localizador falla o se va a más de R del candidato, o si z < z_min. Si dos candidatos
    convergen a menos de R/2 se queda el de mayor z.

    Devuelve dict: cuadro [m] int, x, y [m] px, z [m], n_por_cuadro [n] int, más radio_px,
    umbral, metodo, lado_recorte_px, z_min, ms_por_cuadro (detección + localización).
    """
    if metodo not in tracking.METODOS:
        raise ValueError(f'método desconocido: {metodo!r}; opciones {tracking.METODOS}')
    R = float(radio_px)
    lado = int(lado_recorte_px or tracking.lado_recorte_sugerido(R)) | 1
    suav = 0.2 * R if metodo == 'simetria_radial' else 0.0
    C, X, Y, Z = [], [], [], []
    n_por = np.zeros(n, int)
    t0 = time.perf_counter()
    for i in range(n):
        img = np.asarray(obtener_cuadro(i), float)
        cand = candidatos(img, R, umbral, max_beads)
        res = []
        for xc, yc in cand:
            xb, yb = tracking._localizar_en(img, xc, yc, lado, metodo, True, suav)
            if not (np.isfinite(xb) and np.isfinite(yb)) or np.hypot(xb - xc, yb - yc) > R:
                continue
            rec, j0, i0 = tracking._recortar(img, xb, yb, lado)
            if rec is None:
                continue
            z = tracking.contraste_z(rec, xb - j0, yb - i0, R)
            if z_min is not None and not z >= z_min:
                continue
            res.append((xb, yb, z))
        if res:
            res = np.array(res, float)
            res = res[np.argsort(-np.nan_to_num(res[:, 2], nan=-np.inf))]
            keep = []
            for k in range(len(res)):
                if all(np.hypot(*(res[k, :2] - res[q, :2])) >= 0.5 * R for q in keep):
                    keep.append(k)
            res = res[keep]
            C.append(np.full(len(res), i))
            X.append(res[:, 0])
            Y.append(res[:, 1])
            Z.append(res[:, 2])
            n_por[i] = len(res)
        if progreso is not None:
            progreso((i + 1) / n)
    ms = 1e3 * (time.perf_counter() - t0) / max(n, 1)
    cat = (lambda L, dt: np.concatenate(L).astype(dt) if L else np.zeros(0, dt))
    return dict(cuadro=cat(C, int), x=cat(X, float), y=cat(Y, float), z=cat(Z, float),
                n_por_cuadro=n_por, radio_px=R, umbral=float(umbral), metodo=metodo,
                lado_recorte_px=lado, z_min=z_min, ms_por_cuadro=ms)


# ---------------------------------------------------------------- enlace
_GRANDE = 1e6


def _asignar(costo, ok):
    """Asignación óptima trayectorias (filas) ↔ detecciones (columnas) con 'no asignar' posible.

    costo [nT, nD] = (d / d_max)² (<= 1 donde ok). Dejar una trayectoria o una detección sin
    pareja cuesta 1 cada una, así que un par permitido (costo <= 1 < 2) siempre conviene salvo
    que compita con otro: se minimiza la suma de desplazamientos² normalizados (máxima
    verosimilitud para pasos gaussianos, como trackpy). Matriz aumentada (Jaqaman et al., Nat.
    Methods 5, 695, 2008) resuelta por componentes conexas del grafo de pares permitidos (los
    problemas reales son ralos: casi todo es 1 × 1). Devuelve lista de (fila, columna).
    """
    nT, nD = ok.shape
    filas, cols = np.nonzero(ok)
    if filas.size == 0:
        return []
    g = coo_matrix((np.ones(filas.size), (filas, nT + cols)), shape=(nT + nD, nT + nD))
    _, etiqueta = connected_components(g, directed=False)
    pares = []
    for e in np.unique(etiqueta[filas]):
        tf = np.nonzero(etiqueta[:nT] == e)[0]
        dc = np.nonzero(etiqueta[nT:] == e)[0]
        if tf.size == 1 and dc.size == 1:
            pares.append((tf[0], dc[0]))
            continue
        a, b = tf.size, dc.size
        sub_ok = ok[np.ix_(tf, dc)]
        A = np.zeros((a + b, b + a))
        A[:a, :b] = np.where(sub_ok, costo[np.ix_(tf, dc)], _GRANDE)
        A[:a, b:] = _GRANDE
        A[:a, b:][np.arange(a), np.arange(a)] = 1.0
        A[a:, :b] = _GRANDE
        A[a:, :b][np.arange(b), np.arange(b)] = 1.0
        r, c = linear_sum_assignment(A)
        for i, j in zip(r, c):
            if i < a and j < b and sub_ok[i, j]:
                pares.append((tf[i], dc[j]))
    return pares


def enlazar(det, t, max_salto_px, memoria_cuadros, prediccion_cuadros=None):
    """Id de trayectoria de cada detección (enlace cuadro a cuadro).

    det             : dict con cuadro [m] int, x [m], y [m] px (p. ej. de detectar_todas). Filas con
                      x o y NaN quedan sin enlazar.
    t               : [n] s, tiempo real de cada cuadro (TIEMPO['t']; NaN = no confiable, se
                      interpola en el índice) o None = índice del cuadro.
    max_salto_px    : d_max [px] entre la posición esperada y la detección, para cuadros
                      consecutivos.
    memoria_cuadros : una trayectoria que no aparece en memoria_cuadros + 1 cuadros seguidos se
                      cierra (no se le asigna nada más).
    prediccion_cuadros : None = CRITERIOS. La posición esperada es la última + v·Δt, con v = la
                      velocidad media de los últimos `prediccion_cuadros` pasos de esa trayectoria
                      (una recién nacida usa la mediana de las demás). 0 = la última posición.
                      Por qué: con Stokes (platina a velocidad constante) o deriva, la bead se
                      corre v·Δt por cuadro; sin predicción max_salto tendría que cubrir eso y
                      se vuelve ambiguo con beads cerca, y en un hueco de k cuadros el corrimiento
                      es k·v·Δt. Para movimiento browniano la velocidad media de k pasos ~ 0 y
                      agrega solo σ²/k a la varianza del paso (despreciable con k = 10).

    En cada cuadro, entre las trayectorias abiertas y las detecciones, el par (k, j) está
    PERMITIDO si d_kj <= max_salto_px · sqrt(Δt / dt), con Δt el tiempo real desde la última
    detección de k y dt el período mediano. Por qué la raíz: en difusión el desplazamiento crece
    como sqrt(Δt), así un hueco de k cuadros (o un cuadro perdido de la cámara) agranda la
    búsqueda como sqrt(k + 1). Entre los permitidos se elige la asignación que minimiza
    Σ (d / d_max)² (ver _asignar). Las detecciones sin pareja empiezan trayectorias nuevas.

    Devuelve particula [m] int (mismo orden que det; −1 = no enlazada). Ids en orden de
    aparición.
    """
    pred_k = CRITERIOS['prediccion_cuadros'][0] if prediccion_cuadros is None else int(prediccion_cuadros)
    memoria = int(memoria_cuadros)
    salto = float(max_salto_px)
    cuadro = np.asarray(det['cuadro'], int)
    x = np.asarray(det['x'], float)
    y = np.asarray(det['y'], float)
    particula = np.full(cuadro.size, -1, int)
    validas = np.isfinite(x) & np.isfinite(y) & (cuadro >= 0)
    if not validas.any():
        return particula
    n = int(cuadro[validas].max()) + 1
    if t is not None:
        n = max(n, len(np.asarray(t).ravel()))
    tt, dt = _tiempos(t, n)
    orden = np.nonzero(validas)[0]
    orden = orden[np.argsort(cuadro[orden], kind='stable')]
    cortes = np.nonzero(np.diff(cuadro[orden]))[0] + 1
    abiertas = []                       # ids
    hist = []                           # por id: deque de (t, x, y)
    ultimo = []                         # por id: último cuadro
    for grupo in np.split(orden, cortes):
        f = int(cuadro[grupo[0]])
        tf = tt[f]
        abiertas = [k for k in abiertas if f - ultimo[k] - 1 <= memoria]
        P = np.column_stack([x[grupo], y[grupo]])
        asignada = np.zeros(len(grupo), bool)
        if abiertas:
            vs = []
            for k in abiertas:
                h = hist[k]
                vs.append(((h[-1][1] - h[0][1]) / (h[-1][0] - h[0][0]),
                           (h[-1][2] - h[0][2]) / (h[-1][0] - h[0][0]))
                          if pred_k > 0 and len(h) >= 2 and h[-1][0] > h[0][0] else None)
            conv = [v for v in vs if v is not None]
            vmed = np.median(np.array(conv), axis=0) if conv else np.zeros(2)
            Q = np.empty((len(abiertas), 2))
            gate = np.empty(len(abiertas))
            for q, k in enumerate(abiertas):
                tl, xl, yl = hist[k][-1]
                v = vs[q] if vs[q] is not None else (vmed if pred_k > 0 else (0.0, 0.0))
                Q[q] = xl + v[0] * (tf - tl), yl + v[1] * (tf - tl)
                gate[q] = salto * np.sqrt(max(1.0, (tf - tl) / dt))
            d2 = ((Q[:, None, :] - P[None, :, :]) ** 2).sum(axis=2) / (gate[:, None] ** 2)
            for q, j in _asignar(d2, d2 <= 1.0):
                k = abiertas[q]
                particula[grupo[j]] = k
                hist[k].append((tf, P[j, 0], P[j, 1]))
                ultimo[k] = f
                asignada[j] = True
        for j in np.nonzero(~asignada)[0]:
            k = len(hist)
            hist.append(deque([(tf, P[j, 0], P[j, 1])], maxlen=max(pred_k, 1) + 1))
            ultimo.append(f)
            abiertas.append(k)
            particula[grupo[j]] = k
    return particula


# ---------------------------------------------------------------- postproceso
def postprocesar(cuadro, particula, x, y, z, t, radio_px, criterios):
    """Trayectorias enlazadas → trayectorias limpias en formato largo.

    cuadro, particula, x, y, z : [m] (formato largo; z puede ser None). particula < 0 o x/y NaN
                 = no se usan.
    t          : [n] s tiempos reales de los cuadros (o None = índice). Para interpolar en t real.
    radio_px   : radio de la bead [px] (define dist_colision_px = 0).
    criterios  : dict (ver CRITERIOS; las claves que falten = defecto).

    Pasos, en este orden, y por qué:
      1. Huecos > memoria_cuadros cortan la trayectoria (sigue con otro id): después de tanto
         tiempo sin verla no hay garantía de que sea la misma (en `enlazar` ya no pasa; en
         `desde_seguir` sí).
      2. Colisiones: en cada cuadro, dos partículas a menos de dist_colision_px (0 = 2.5 R). Por
         qué 2.5 R: el recorte del localizador tiene medio lado ~1.5 R y el anillo de la vecina
         llega a ~R de su centro; más cerca que 1.5 R + R la vecina entra al recorte y corre la
         posición (y el enlace puede intercambiar ids). En los huecos ≤ memoria la posición se
         estima interpolando, así se detecta también la colisión en la que una de las dos dejó
         de verse (las dos se funden en una sola detección al cruzarse). Si la fusión dura más
         que la memoria se reconoce porque una trayectoria termina pegada a otra y luego otra
         empieza pegada a la misma (ver _fusiones): esos cuadros intermedios también son
         colisión.
         accion_colision: 'marcar' → colision = True y se guarda igual; 'descartar' → x, y, z =
         NaN en esos cuadros (la fila queda, colision = True); 'cortar' → esos cuadros se sacan
         y la trayectoria sigue con otro id después de la colisión.
      3. Interpolación: huecos de 1..interpolar_max_cuadros cuadros se rellenan con la recta
         entre las dos mediciones vecinas, en t real (así un cuadro perdido de la cámara no
         deforma la recta); interpolado = True, z = NaN. Con 'descartar'/'cortar' no se
         interpola hacia, desde ni a través de una colisión.
      4. Trayectorias con menos de largo_min_cuadros posiciones MEDIDAS (finitas, no
         interpoladas) se descartan: no alcanzan para un MSD ni una velocidad.
      5. Ids renumerados 0..P-1 en orden de (id original, cuadro).

    Devuelve dict (formato largo, ordenado por particula y cuadro): cuadro [m] int, particula
    [m] int, x, y [m] px, z [m], interpolado [m] bool, colision [m] bool, hueco_previo [m] int
    (cuadros que faltan entre esta fila y la anterior de la misma partícula; 0 en la primera),
    particula_origen [m] int (id de entrada), resumen = dict de arrays por partícula
    (particula, particula_origen, primer_cuadro, ultimo_cuadro, n_cuadros (filas con posición),
    n_medidos, n_interpolados, n_colision, x_medio, y_medio (de lo medido), desplazamiento_total_px
    (|r_último − r_primero|)) + contadores globales: detecciones (filas de entrada con posición),
    no_enlazadas, descartadas_cortas (trayectorias), colisiones (pares partícula-cuadro en
    colisión, antes de descartar cortas), cortes (= cortes_hueco + cortes_colision),
    cortes_hueco, cortes_colision; y criterios = los usados.
    """
    c = criterios_completos(criterios)
    memoria = int(c['memoria_cuadros'])
    interp_max = int(c['interpolar_max_cuadros'])
    largo_min = int(c['largo_min_cuadros'])
    dcol = float(c['dist_colision_px']) or 2.5 * float(radio_px)
    accion = c['accion_colision']

    cuadro = np.asarray(cuadro, int).ravel()
    particula = np.asarray(particula, int).ravel()
    x = np.asarray(x, float).ravel()
    y = np.asarray(y, float).ravel()
    z = np.full(x.size, np.nan) if z is None else np.asarray(z, float).ravel()
    con_pos = np.isfinite(x) & np.isfinite(y)
    usar = con_pos & (particula >= 0) & (cuadro >= 0)
    cont = dict(detecciones=int(con_pos.sum()), no_enlazadas=int((con_pos & (particula < 0)).sum()))
    sel = np.nonzero(usar)[0]
    sel = sel[np.lexsort((cuadro[sel], particula[sel]))]
    if sel.size:                                  # una sola fila por (partícula, cuadro)
        dup = np.r_[False, (np.diff(particula[sel]) == 0) & (np.diff(cuadro[sel]) == 0)]
        sel = sel[~dup]
    C, P, X, Y, Zz = cuadro[sel], particula[sel], x[sel], y[sel], z[sel]
    n = int(C.max()) + 1 if C.size else 0
    if t is not None:
        n = max(n, len(np.asarray(t).ravel()))
    tt, _ = _tiempos(t, n)

    # 1. segmentos: cambia la partícula o hueco > memoria
    nuevo = np.r_[True, (np.diff(P) != 0) | (np.diff(C) > memoria + 1)] if C.size else np.zeros(0, bool)
    cortes_hueco = int((nuevo[1:] & (np.diff(P) == 0)).sum()) if C.size else 0
    seg = np.cumsum(nuevo) - 1

    # 2. representación densa: cada segmento en todos sus cuadros; los huecos (≤ memoria) con la
    #    recta en t real entre las dos mediciones vecinas
    fin_seg = np.r_[nuevo[1:], True] if C.size else np.zeros(0, bool)
    g = np.where(fin_seg, 0, np.r_[np.diff(C), 1] - 1) if C.size else np.zeros(0, int)
    rep = 1 + g
    k = np.repeat(np.arange(C.size), rep)                     # medición de la que sale la fila
    off = np.arange(k.size) - np.repeat(np.cumsum(rep) - rep, rep)
    dC = C[k] + off
    medido = off == 0
    k2 = np.minimum(k + 1, max(C.size - 1, 0))
    with np.errstate(invalid='ignore', divide='ignore'):
        fr = np.where(medido, 0.0, (tt[dC] - tt[C[k]]) / (tt[C[k2]] - tt[C[k]]))
    dX = X[k] + fr * (X[k2] - X[k])
    dY = Y[k] + fr * (Y[k2] - Y[k])
    dZ = np.where(medido, Zz[k], np.nan)
    dS = seg[k]
    hueco_len = g[k]                                          # largo del hueco de cada fila virtual

    # colisiones en cada cuadro (árbol en (x, y, cuadro·L): solo se tocan puntos del mismo cuadro)
    col = np.zeros(k.size, bool)
    if k.size > 1:
        L = 10.0 * (np.nanmax(np.abs(np.r_[dX, dY])) + dcol + 1.0)
        pts = np.column_stack([dX, dY, dC * L])
        arbol = cKDTree(pts)
        pares = arbol.query_pairs(dcol, output_type='ndarray')
        if len(pares):
            pares = pares[np.hypot(dX[pares[:, 0]] - dX[pares[:, 1]], dY[pares[:, 0]] - dY[pares[:, 1]]) < dcol]
            col[pares.ravel()] = True
        _fusiones(col, arbol, pts, dS, dC, dcol, n, 10 * (memoria + 1))
    cont['colisiones'] = int(col.sum())

    # 3. qué filas virtuales se rellenan
    hueco_col = np.zeros(C.size, bool)                        # el hueco que sigue a la medición k toca colisión
    if accion != 'marcar' and k.size:
        np.logical_or.at(hueco_col, k[~medido], col[~medido])
        mcol = np.zeros(C.size, bool)
        mcol[k[medido]] = col[medido]
        hueco_col |= mcol | mcol[np.minimum(np.arange(C.size) + 1, C.size - 1)]
    rellenar = ~medido & (hueco_len <= interp_max) & ~hueco_col[k]
    fila = medido | rellenar
    if accion == 'cortar':
        caida = col          # medidas o estimadas en un hueco: la trayectoria se parte ahí
        nuevo_s = np.r_[True, (dS[1:] != dS[:-1]) | caida[:-1]] if k.size else np.zeros(0, bool)
        dS2 = np.cumsum(nuevo_s) - 1
        fila &= ~caida
        cortes_col = len(np.unique(dS2[fila])) - len(np.unique(dS[fila]))
        dS = dS2
    else:
        cortes_col = 0
        if accion == 'descartar':
            dX, dY, dZ = dX.copy(), dY.copy(), dZ.copy()
            dX[medido & col] = dY[medido & col] = dZ[medido & col] = np.nan
    cont.update(cortes_hueco=cortes_hueco, cortes_colision=int(cortes_col),
                cortes=cortes_hueco + int(cortes_col))

    oC, oS, oX, oY, oZ = dC[fila], dS[fila], dX[fila], dY[fila], dZ[fila]
    oI, oCol, oP = ~medido[fila], col[fila], P[k[fila]]

    # 4. cortas (cuenta lo medido con posición)
    ids, inv = np.unique(oS, return_inverse=True)
    n_med = np.bincount(inv, weights=(~oI & np.isfinite(oX)).astype(float), minlength=ids.size)
    larga = n_med >= largo_min
    cont['descartadas_cortas'] = int((~larga).sum())
    q = larga[inv]
    oC, oS, oX, oY, oZ, oI, oCol, oP = oC[q], inv[q], oX[q], oY[q], oZ[q], oI[q], oCol[q], oP[q]

    # 5. ids 0..P-1
    _, oS = np.unique(oS, return_inverse=True)
    oS = oS.astype(int)
    hueco = np.r_[0, np.where(np.diff(oS) == 0, np.diff(oC) - 1, 0)].astype(int) if oC.size else np.zeros(0, int)
    res = _resumen(oS, oP, oC, oX, oY, oI, oCol)
    res.update(cont)
    return dict(cuadro=oC.astype(int), particula=oS, x=oX, y=oY, z=oZ, interpolado=oI.astype(bool),
                colision=oCol.astype(bool), hueco_previo=hueco, particula_origen=oP.astype(int),
                resumen=res, criterios=c, dist_colision_px=dcol)


def _fusiones(col, arbol, pts, dS, dC, dcol, n, max_cuadros):
    """Marca (en `col`, in place) los cuadros en que dos beads se FUNDIERON en una detección.

    Al cruzarse, dos beads a menos de ~1.5 R dan un solo pico: una trayectoria (B) termina pegada
    a otra (A) en el cuadro e y, cuando se separan, otra trayectoria (B', a veces con el id de A
    intercambiado) empieza pegada a A en el cuadro f. Si la fusión dura más que la memoria, B y
    B' son segmentos distintos y en e < cuadro < f solo existe A, así que el chequeo de pares no
    ve nada. Regla: si un segmento termina (antes del último cuadro) a menos de dcol de A en e,
    y otro empieza a menos de dcol de A en f, con e < f <= e + max_cuadros (el primero que
    aparezca), las filas de A en e < cuadro < f son colisión (posición = el promedio de dos
    beads). max_cuadros = 10 (memoria + 1): acota el falso positivo de una bead que se va de
    foco junto a A y otra que aparece mucho después.
    """
    nS = int(dS.max()) + 1
    a = np.searchsorted(dS, np.arange(nS))
    b = np.r_[a[1:], dS.size] - 1
    ini, fin = dC[a], dC[b]

    def vecinos(idx):
        out = []
        for q, lst in zip(idx, arbol.query_ball_point(pts[idx], dcol)):
            out.append([dS[j] for j in lst if dS[j] != dS[q]
                        and np.hypot(*(pts[j, :2] - pts[q, :2])) < dcol])
        return out

    fines = np.nonzero(fin < n - 1)[0]
    inicios = np.nonzero(ini > 0)[0]
    if not fines.size or not inicios.size:
        return
    empiezan = {}                                     # ocultadora -> [(f, segmento que empieza)]
    for s, occ in zip(inicios, vecinos(a[inicios])):
        for o in occ:
            empiezan.setdefault(o, []).append((ini[s], s))
    for s, occ in zip(fines, vecinos(b[fines])):
        e = fin[s]
        for o in occ:
            fs = [f for f, s2 in empiezan.get(o, []) if e < f <= e + max_cuadros and s2 != s]
            if fs:
                f = min(fs)
                lo, hi = max(e + 1, ini[o]), min(f - 1, fin[o])
                if lo <= hi:
                    col[a[o] + lo - ini[o]:a[o] + hi - ini[o] + 1] = True


def _resumen(S, Pori, C, X, Y, I, Col):
    """Arrays por partícula (S ya es 0..P-1, filas ordenadas por S y cuadro)."""
    nP = int(S.max()) + 1 if S.size else 0
    if nP == 0:
        vac = lambda dt: np.zeros(0, dt)
        return dict(particula=vac(int), particula_origen=vac(int), primer_cuadro=vac(int),
                    ultimo_cuadro=vac(int), n_cuadros=vac(int), n_medidos=vac(int),
                    n_interpolados=vac(int), n_colision=vac(int), x_medio=vac(float),
                    y_medio=vac(float), desplazamiento_total_px=vac(float))
    ini = np.searchsorted(S, np.arange(nP))
    fin = np.r_[ini[1:], S.size] - 1
    fin_ok = np.isfinite(X)
    med = fin_ok & ~I
    cnt = lambda m: np.bincount(S, weights=m.astype(float), minlength=nP).astype(int)
    with np.errstate(invalid='ignore', divide='ignore'):
        xm = np.bincount(S[med], X[med], nP) / np.bincount(S[med], minlength=nP)
        ym = np.bincount(S[med], Y[med], nP) / np.bincount(S[med], minlength=nP)
    desp = np.full(nP, np.nan)
    for p in range(nP):
        idx = np.nonzero(fin_ok[ini[p]:fin[p] + 1])[0] + ini[p]
        if idx.size:
            desp[p] = np.hypot(X[idx[-1]] - X[idx[0]], Y[idx[-1]] - Y[idx[0]])
    return dict(particula=np.arange(nP), particula_origen=Pori[ini].astype(int),
                primer_cuadro=C[ini].astype(int), ultimo_cuadro=C[fin].astype(int),
                n_cuadros=cnt(fin_ok), n_medidos=cnt(med), n_interpolados=cnt(I),
                n_colision=cnt(Col), x_medio=xm, y_medio=ym, desplazamiento_total_px=desp)


# ---------------------------------------------------------------- de punta a punta
def seguir_todas(obtener_cuadro, n, t, radio_px, criterios, umbral=6.0, metodo='simetria_radial',
                 lado_recorte_px=None, progreso=None, max_beads=200):
    """detectar_todas + enlazar + postprocesar (beads libres, Stokes).

    obtener_cuadro, n, radio_px, umbral, metodo, lado_recorte_px, max_beads: ver detectar_todas.
    t: [n] s (TIEMPO['t']) o None. criterios: ver CRITERIOS (z_min sale de acá).
    progreso(fraccion): detección 0 → 0.95, enlace y postproceso hasta 1.
    Devuelve lo de postprocesar, más deteccion = dict(n_detecciones, n_por_cuadro, ms_por_cuadro,
    lado_recorte_px, metodo) y ms_enlace (enlace + postproceso, total).
    """
    c = criterios_completos(criterios)
    prog = (lambda f: progreso(0.95 * f)) if progreso is not None else None
    det = detectar_todas(obtener_cuadro, n, radio_px, umbral, max_beads, metodo, lado_recorte_px,
                         c['z_min'], prog)
    t0 = time.perf_counter()
    p = enlazar(det, t, c['max_salto_px'], c['memoria_cuadros'], c['prediccion_cuadros'])
    r = postprocesar(det['cuadro'], p, det['x'], det['y'], det['z'], t, radio_px, c)
    r['ms_enlace'] = 1e3 * (time.perf_counter() - t0)
    r['deteccion'] = dict(n_detecciones=int(det['cuadro'].size), n_por_cuadro=det['n_por_cuadro'],
                          ms_por_cuadro=det['ms_por_cuadro'], lado_recorte_px=det['lado_recorte_px'],
                          metodo=metodo)
    if progreso is not None:
        progreso(1.0)
    return r


def desde_seguir(r_seguir, t, radio_px, criterios):
    """Resultado de tracking.seguir ([n, nb] con NaN) → formato largo postprocesado.

    r_seguir : dict de tracking.seguir (x, y [n, nb] px; contraste_z [n, nb] opcional).
    t        : [n] s o None. radio_px: radio [px]. criterios: ver CRITERIOS (z_min, max_salto y
               prediccion no se usan acá: ya los aplicó seguir).
    La bead b es la partícula b; sus cuadros con NaN son huecos: ≤ memoria no cortan (y se
    interpolan si interpolar_max_cuadros lo pide), > memoria cortan (sigue con otro id).
    Colisiones, cortas e ids como en postprocesar. Sin huecos ni colisiones y con la bead más
    larga que largo_min_cuadros, x/y de la partícula b son EXACTAMENTE r_seguir['x'][:, b].
    """
    X = np.asarray(r_seguir['x'], float)
    Y = np.asarray(r_seguir['y'], float)
    if X.ndim == 1:                                           # una sola bead
        X, Y = X[:, None], Y[:, None]
    n, nb = X.shape
    Z = np.asarray(r_seguir.get('contraste_z', np.full((n, nb), np.nan)), float).reshape(n, nb)
    b, i = np.meshgrid(np.arange(nb), np.arange(n))           # [n, nb]
    ok = np.isfinite(X) & np.isfinite(Y)
    orden = np.argsort(b[ok], kind='stable')
    r = postprocesar(i[ok][orden], b[ok][orden], X[ok][orden], Y[ok][orden], Z[ok][orden],
                     t if t is not None else None, radio_px, criterios)
    return r
