"""Datos sintéticos de parámetros conocidos, para validar todo el análisis.

Tres capas, de abajo hacia arriba:
  1. dinámica   — trayectorias exactas (difusión libre u Ornstein-Uhlenbeck) promediadas
                  durante la exposición W (así aparece el blur).
  2. imagen     — una bead de campo claro dibujada en una posición subpíxel, fondo despar
                  ejo, ruido de disparo y de lectura, 10 bits, mosaico Bayer opcional.
  3. grabación  — `grabar(...)` escribe una carpeta con el MISMO formato que una real
                  (ver CONTRATO.md), con artefactos de tiempo y la verdad en grabacion.json.
"""
import numpy as np
from scipy.signal import lfilter

from .fisica import KB


# ---------------------------------------------------------------- 1. dinámica
def trayectoria(t_inicio, W, D, kappa=0.0, gamma=None, T=None, sub=40, rng=None, x0=None):
    """Trayectoria en UN eje promediada durante la exposición de cada cuadro.

    t_inicio : [n] s, instante en que empieza la exposición de cada cuadro (creciente,
               puede tener huecos: cuadros perdidos o jitter).
    W        : exposición [s].
    D        : coeficiente de difusión [m²/s] (= kB T / gamma).
    kappa    : rigidez [N/m]; 0 = difusión libre.
    gamma    : arrastre [N·s/m]; hace falta si kappa > 0 (tau_c = gamma / kappa).
    sub      : pasos finos por exposición.
    x0       : posición inicial [m]; None = 0 si libre, equilibrio sorteado si atrapada.

    Integra con la solución EXACTA del paso (no Euler) en una grilla fina uniforme de paso
    h = W / sub: libre  x += sqrt(2 D h) xi;   OU  x = a x + sqrt(kB T/kappa (1 - a²)) xi,
    con a = exp(-h / tau_c). Cada cuadro devuelve el promedio de x sobre [t_i, t_i + W].

    Devuelve (x_medida [n], x_instantanea [n]) en m: la instantánea es x(t_i + W/2).
    """
    rng = np.random.default_rng(rng)
    t_inicio = np.asarray(t_inicio, float)
    h = W / sub
    k0 = np.round((t_inicio - t_inicio[0]) / h).astype(np.int64)    # índice fino de cada inicio
    nfino = int(k0[-1]) + sub + 1
    xi = rng.standard_normal(nfino)
    if kappa > 0:
        if gamma is None:
            raise ValueError('con kappa > 0 hace falta gamma')
        tau = gamma / kappa
        a = np.exp(-h / tau)
        var_eq = D * gamma / kappa            # = kB T / kappa
        s = np.sqrt(var_eq * (1 - a * a))
        inicio = rng.standard_normal() * np.sqrt(var_eq) if x0 is None else x0
        x, _ = lfilter([1.0], [1.0, -a], s * xi, zi=[a * inicio])
    else:
        inicio = 0.0 if x0 is None else x0
        x = inicio + np.cumsum(np.sqrt(2 * D * h) * xi)
    c = np.concatenate([[0.0], np.cumsum(x)])
    # promedio de x[k0 .. k0+sub-1] (sub muestras dentro de la exposición)
    medida = (c[k0 + sub] - c[k0]) / sub
    inst = x[k0 + sub // 2]
    return medida, inst


# ---------------------------------------------------------------- 2. imagen
def perfil_bead(r, radio_px, ancho_anillo_px=None):
    """Perfil radial de una bead de campo claro (adimensional, fondo = 0).

    Centro oscuro y un anillo brillante en el borde: -exp(-r²/2s²) + 0.6 exp(-(r-R)²/2w²).
    No es óptica de Mie: alcanza para tener la forma (centro + anillos) que ve la cámara.
    """
    s = 0.6 * radio_px
    w = ancho_anillo_px or max(0.8, 0.25 * radio_px)
    return -np.exp(-r ** 2 / (2 * s * s)) + 0.6 * np.exp(-(r - radio_px) ** 2 / (2 * w * w))


def dibujar(forma, posiciones_px, radio_px, contraste=0.35, fondo=None, super_muestreo=4):
    """Imagen ideal (sin ruido) en unidades de "fracción del fondo", con beads en posiciones subpíxel.

    forma: (alto, ancho). posiciones_px: [nb, 2] (x, y) en px, convención de CONTRATO.md.
    fondo: imagen de iluminación (1 = media) o None para iluminación despareja por defecto.
    Integra cada píxel con super_muestreo² sub-muestras (la cámara integra el área del píxel).
    Devuelve float64 [alto, ancho], ~1 en el fondo.
    """
    alto, ancho = forma
    if fondo is None:
        yy, xx = np.mgrid[0:alto, 0:ancho]
        fondo = 1 + 0.08 * np.cos(2 * np.pi * xx / (2.3 * ancho)) * np.cos(2 * np.pi * yy / (1.7 * alto)) \
            + 0.04 * (xx / ancho - 0.5)
    img = np.array(fondo, float, copy=True)
    m = super_muestreo
    off = (np.arange(m) + 0.5) / m - 0.5
    L = int(np.ceil(radio_px * 2.5)) + 2
    for x, y in np.atleast_2d(posiciones_px):
        if not (np.isfinite(x) and np.isfinite(y)):
            continue
        j0, j1 = max(0, int(x) - L), min(ancho, int(x) + L + 1)
        i0, i1 = max(0, int(y) - L), min(alto, int(y) + L + 1)
        if j0 >= j1 or i0 >= i1:
            continue
        ii, jj = np.mgrid[i0:i1, j0:j1]
        acc = np.zeros(ii.shape)
        for dy in off:
            for dx in off:
                acc += perfil_bead(np.hypot(jj + dx - x, ii + dy - y), radio_px)
        img[i0:i1, j0:j1] += contraste * fondo[i0:i1, j0:j1] * acc / (m * m)
    return img


def camara(ideal, nivel=600.0, ruido_lectura=2.5, bits=10, rng=None):
    """De imagen ideal a cuadro crudo: ruido de disparo (Poisson) + lectura (gauss), cuantizado.

    ideal: float, ~1 en el fondo. nivel: cuentas del fondo (en ADU, ganancia 1 e-/ADU).
    Devuelve uint16 en [0, 2^bits - 1].
    """
    rng = np.random.default_rng(rng)
    e = rng.poisson(np.clip(ideal, 0, None) * nivel).astype(float)
    e += rng.normal(0, ruido_lectura, e.shape)
    return np.clip(np.round(e), 0, 2 ** bits - 1).astype(np.uint16)


def mosaico_bayer(rgb, fase='RGGB'):
    """Muestrea una imagen [alto, ancho, 3] con un patrón Bayer -> [alto, ancho] (un canal por píxel)."""
    canal = {'R': 0, 'G': 1, 'B': 2}
    out = np.empty(rgb.shape[:2], rgb.dtype)
    for k, (di, dj) in enumerate([(0, 0), (0, 1), (1, 0), (1, 1)]):
        out[di::2, dj::2] = rgb[di::2, dj::2, canal[fase[k]]]
    return out


# ---------------------------------------------------------------- 3. grabación
TICK_HZ_NOMINAL = 1e9            # lo que "informa el SDK" en el sintético (el reloj real difiere)
GANANCIAS_RGB = (0.85, 1.0, 0.6)  # respuesta relativa R, G, B a la lámpara (solo para Bayer)


def tiempos_camara(n, fps, W, tick_hz_real=1.0000137e9, tick_hz_nominal=TICK_HZ_NOMINAL,
                   jitter_host_s=2e-3, prob_perdido=0.002, latencia_s=3e-3, prob_cola=0.002,
                   cola_s=(0.01, 0.05), jitter_camara_s=1e-6, procesar_s=20e-6, contador0=1,
                   rng=None):
    """Relojes de una grabación: lo que la cámara y el host anotarían por cuadro.

    n   : cuadros que produce la cámara (contador contador0 .. contador0 + n - 1).
    fps : tasa pedida; la cámara la arma con SU reloj: período = round(tick_hz_nominal / fps)
          ticks, que en segundos verdaderos dura período / tick_hz_real (no 1/fps exacto).
    Inicio de exposición verdadero: t_k = k·período/tick_hz_real + N(0, jitter_camara_s).
    t_camara_k = t0 + round(t_k · tick_hz_real)  (ticks del reloj REAL; se supone que marca el
    inicio de la exposición — a confirmar con la Zelux; un corrimiento fijo no cambia nada).
    Host: el cuadro está listo en t_k + W + latencia_s y llega con una demora extra
    exponencial de media jitter_host_s (siempre >= 0: el host no puede recibir antes); con
    probabilidad prob_cola se suma una cola U(cola_s). La entrega es en orden (FIFO) y cada
    cuadro tarda procesar_s: llegada_k = max(listo_k + demora_k, llegada_{k-1} + procesar_s).
    Perdidos: cada cuadro (salvo el primero y el último) se pierde con probabilidad prob_perdido;
    un perdido no se guarda (no tiene fila).
    Devuelve dict: columnas (contador_camara, t_camara, t_host_ns, exposicion_us) de los
    GUARDADOS, t_verdad [n_g] s desde el inicio de exposición del primer guardado,
    contadores_perdidos [n_p], periodo_real_s.
    """
    rng = np.random.default_rng(rng)
    k = np.arange(n)
    P_ticks = int(round(tick_hz_nominal / fps))
    t = k * (P_ticks / tick_hz_real) + rng.normal(0, jitter_camara_s, n)
    t -= t[0]
    t0_ticks = 4_321_000_000_000
    t_cam = t0_ticks + np.round(t * tick_hz_real).astype(np.int64)
    demora = rng.exponential(jitter_host_s, n) if jitter_host_s > 0 else np.zeros(n)
    cola = rng.random(n) < prob_cola
    demora[cola] += rng.uniform(cola_s[0], cola_s[1], cola.sum())
    listo = t + W + latencia_s + demora
    llegada = np.maximum.accumulate(listo - procesar_s * k) + procesar_s * k
    t_host = np.round((123456.789 + llegada) * 1e9).astype(np.int64)
    perdido = rng.random(n) < prob_perdido
    perdido[[0, -1]] = False
    g = ~perdido
    contador = contador0 + k
    return dict(contador_camara=contador[g], t_camara=t_cam[g], t_host_ns=t_host[g],
                exposicion_us=np.full(g.sum(), W * 1e6), t_verdad=t[g] - t[g][0],
                contadores_perdidos=contador[perdido], periodo_real_s=P_ticks / tick_hz_real)


def _iluminacion(forma):
    """La misma iluminación despareja que usa `dibujar` por defecto (1 = media)."""
    alto, ancho = forma
    yy, xx = np.mgrid[0:alto, 0:ancho]
    return 1 + 0.08 * np.cos(2 * np.pi * xx / (2.3 * ancho)) * np.cos(2 * np.pi * yy / (1.7 * alto)) \
        + 0.04 * (xx / ancho - 0.5)


def _a_crudo(ideal, bayer, nivel, ruido_lectura, rng):
    """Imagen ideal (gris, ~1 en el fondo) → cuadro crudo uint16 de 10 bits (mosaico si bayer)."""
    if bayer:
        ideal = mosaico_bayer(ideal[..., None] * np.asarray(GANANCIAS_RGB)[None, None, :], bayer)
    return camara(ideal, nivel=nivel, ruido_lectura=ruido_lectura, rng=rng)


def grabar(carpeta_raiz, etiqueta, tipo, n, fps, W, beads, m_por_px=9e-8, forma=(128, 128),
           bayer=None, tick_hz_real=1.0000137e9, jitter_host_s=2e-3, prob_perdido=0.002,
           deriva_m_s=(0.0, 0.0), nivel=600.0, ruido_lectura=2.5, seed=None, T=None,
           contraste=0.35, con_timestamp=True, **kw_tiempos):
    """Escribe una grabación sintética con el formato de CONTRATO.md y devuelve su carpeta.

    carpeta_raiz : dónde crear <AAAA-MM-DD_HHMMSS>_<etiqueta>/ (en los tests, un tmp_path).
    tipo  : 'libre' | 'atrapada' | 'stokes' | 'pegada' | 'reticula' (va a grabacion.json).
    n     : cuadros que produce la cámara; se guardan ~n·(1 - prob_perdido) (ver tiempos_camara).
    fps   : tasa pedida [1/s]. W: exposición [s].
    beads : lista de dicts: x0_px, y0_px (px del sensor), D_m2s [m²/s], kappa_x, kappa_y [N/m]
            (0 = libre), gamma [N·s/m] (None = kB T / D), radio_px, pegada (bool: no se mueve).
    m_por_px : escala [m/px del sensor]. forma: (alto, ancho) px. bayer: None (mono) o 'RGGB'...
    tick_hz_real : frecuencia verdadera del reloj de la cámara; en grabacion.json va la nominal
            (1e9), como haría un SDK que informa el valor de catálogo.
    jitter_host_s : media de la demora exponencial del host [s]. prob_perdido: por cuadro.
    deriva_m_s : (vx, vy) [m/s], deriva lineal común a todas las beads (platina).
    nivel, ruido_lectura : cuentas del fondo y ruido de lectura [ADU]. T: temperatura [K].
    con_timestamp : False → t_camara vacío (para probar el método sin timestamp de hardware).
    kw_tiempos : se pasan a tiempos_camara (latencia_s, prob_cola, cola_s, ...).

    Dinámica por eje con `trayectoria` (solución exacta, promediada sobre W), muestreada en los
    instantes verdaderos de los cuadros GUARDADOS (los perdidos no se dibujan, pero el tiempo
    corre igual). Ojo: `trayectoria` lleva cada inicio a la grilla fina de paso W/40.
    La deriva se suma promediada sobre la exposición: v·(t + W/2).
    `verdad` en grabacion.json: los parámetros, t_inicio_s [n_g], x_px, y_px [n_g, nb]
    (promedio durante la exposición), x_inst_px, y_inst_px (instantánea en t + W/2) y los
    contadores perdidos.
    """
    import os
    from datetime import datetime
    from . import io
    from .fisica import T_AMBIENTE

    rng = np.random.default_rng(seed)
    T = T_AMBIENTE if T is None else T
    tm = tiempos_camara(n, fps, W, tick_hz_real=tick_hz_real, jitter_host_s=jitter_host_s,
                        prob_perdido=prob_perdido, rng=rng, **kw_tiempos)
    t = tm['t_verdad']
    ng, nb = len(t), len(beads)
    X, Y = np.zeros((ng, nb)), np.zeros((ng, nb))
    Xi, Yi = np.zeros((ng, nb)), np.zeros((ng, nb))
    beads_v = []
    for b, bd in enumerate(beads):
        bd = dict(bd)
        D = float(bd.get('D_m2s', 0.0))
        bd.setdefault('kappa_x', 0.0)
        bd.setdefault('kappa_y', 0.0)
        bd.setdefault('pegada', False)
        if bd.get('gamma') is None and D > 0:
            bd['gamma'] = KB * T / D
        for eje, kap, M, Mi in (('x', bd['kappa_x'], X, Xi), ('y', bd['kappa_y'], Y, Yi)):
            if bd['pegada'] or D <= 0:
                med = inst = np.zeros(ng)
            else:
                med, inst = trayectoria(t, W, D, kappa=kap, gamma=bd['gamma'], rng=rng)
            v = deriva_m_s[0] if eje == 'x' else deriva_m_s[1]
            p0 = bd['x0_px'] if eje == 'x' else bd['y0_px']
            M[:, b] = p0 + (med + v * (t + W / 2)) / m_por_px
            Mi[:, b] = p0 + (inst + v * (t + W / 2)) / m_por_px
        beads_v.append(bd)

    carpeta = io.nueva_carpeta_grabacion(etiqueta, raiz=carpeta_raiz)
    F = _iluminacion(forma)
    radios = sorted({float(bd.get('radio_px', 6.0)) for bd in beads_v})
    archivos = []
    for j in range(ng):
        ideal = F.copy()
        for r in radios:
            sel = [b for b, bd in enumerate(beads_v) if float(bd.get('radio_px', 6.0)) == r]
            pos = np.stack([X[j, sel], Y[j, sel]], axis=1)
            ideal += dibujar(forma, pos, r, contraste=contraste, fondo=F) - F
        archivos.append(io.escribir_cuadro(carpeta, j, _a_crudo(ideal, bayer, nivel, ruido_lectura, rng)))

    t_cam = tm['t_camara'] if con_timestamp else np.full(ng, np.nan)
    io.escribir_tiempos_csv(carpeta, dict(archivo=archivos, contador_camara=tm['contador_camara'],
                                          t_camara=t_cam, t_host_ns=tm['t_host_ns'],
                                          exposicion_us=tm['exposicion_us']))
    alto, ancho = forma
    meta = dict(
        version=1, nombre=os.path.basename(carpeta), tipo=tipo,
        hora_inicio_iso=datetime.now().isoformat(timespec='milliseconds'),
        muestra=dict(material='sintético', diametro_m=None, glicerina_pct=None, preparada=None),
        condiciones=dict(temperatura_K=T, laser=None, profundidad_m=None,
                         velocidad_platina_m_s=None),
        camara=dict(modelo='sintetico', serie='', exposicion_us=W * 1e6, periodo_us=1e6 / fps,
                    roi=[0, 0, ancho, alto], binning=[1, 1], ganancia=0, bits=10,
                    formato_pixel='bayer' if bayer else 'mono', bayer_fase=bayer,
                    tick_hz=TICK_HZ_NOMINAL if con_timestamp else None),
        escala_m_por_px=m_por_px, notas='sintético (pinzas.sintetico.grabar)',
        verdad=dict(seed=seed, n_camara=n, fps_pedido=fps, periodo_real_s=tm['periodo_real_s'],
                    tick_hz_real=tick_hz_real, W_s=W, jitter_host_s=jitter_host_s,
                    prob_perdido=prob_perdido, contadores_perdidos=tm['contadores_perdidos'],
                    T_K=T, m_por_px=m_por_px, deriva_m_s=list(deriva_m_s), beads=beads_v,
                    nivel=nivel, ruido_lectura=ruido_lectura, con_timestamp=con_timestamp,
                    t_inicio_s=t, x_px=X, y_px=Y, x_inst_px=Xi, y_inst_px=Yi))
    io.escribir_grabacion_json(carpeta, meta)
    return carpeta


def reticula(forma, paso_px, angulo=0.0, tipo='lineas', ancho_px=1.5, contraste=0.5,
             desplazamiento_px=(0.0, 0.0), nivel=600.0, ruido_lectura=2.5, bayer=None,
             super_muestreo=4, rng=None):
    """Cuadro crudo uint16 de una retícula de período paso_px (px del sensor), girada `angulo` rad.

    tipo: 'lineas' (rejilla de líneas oscuras en x y en y) o 'puntos' (puntos oscuros en los
    cruces). Perfil de cada línea: gaussiana de ancho a media altura ancho_px, periódica:
    g(u) = exp(-d²/2s²), d = distancia a la línea más cercana, s = ancho_px / 2.3548.
    Líneas: I = ilum·(1 - c·(g_u + g_v - g_u g_v)); puntos: I = ilum·(1 - c·g_u g_v).
    Se integra el área de cada píxel con super_muestreo² sub-muestras, iluminación despareja
    como `dibujar`, y pasa por `camara` (y por `mosaico_bayer` si bayer).
    """
    rng = np.random.default_rng(rng)
    alto, ancho = forma
    s = ancho_px / 2.3548
    m = super_muestreo
    off = (np.arange(m) + 0.5) / m - 0.5
    yy, xx = np.mgrid[0:alto, 0:ancho].astype(float)
    c, si = np.cos(angulo), np.sin(angulo)
    acc = np.zeros(forma)
    for dy in off:
        for dx in off:
            x = xx + dx - desplazamiento_px[0]
            y = yy + dy - desplazamiento_px[1]
            u, v = x * c + y * si, -x * si + y * c
            gu = np.exp(-(((u + paso_px / 2) % paso_px) - paso_px / 2) ** 2 / (2 * s * s))
            gv = np.exp(-(((v + paso_px / 2) % paso_px) - paso_px / 2) ** 2 / (2 * s * s))
            acc += (gu + gv - gu * gv) if tipo == 'lineas' else gu * gv
    ideal = _iluminacion(forma) * (1 - contraste * acc / (m * m))
    return _a_crudo(ideal, bayer, nivel, ruido_lectura, rng)
