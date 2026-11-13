"""Análisis y visualización de las trayectorias que exporta TouchDesigner (PROCESAR → CSV).

Uso (desde un notebook de esta carpeta):

    from trayectorias import *
    datos = cargar(ultima_salida())          # o cargar(r"C:\\datos_pinzas\\...\\csv\\<fecha>")
    t, X, Y, ids = matrices(datos)           # [cuadros] y [cuadros, partículas], NaN = no medido
    graficar_trayectorias(datos)

Unidades: si la salida trae escala (parámetro Escala de B8 o grabacion.json), las posiciones van
en µm; si no, en px. Todas las funciones devuelven también la unidad que usaron.

Física del MSD (2D, todos los orígenes, con los tiempos reales): pinzas/libre.py.
MSD(τ) = 4 D τ + 4 σ²   (σ = ruido de localización por eje; sin blur porque el MP4 no trae la
exposición). Con deriva v: MSD = 4 D τ + v² τ² + 4σ² → primero se resta la deriva común.
"""
import glob
import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)                                   # pinzas ópticas/
if RAIZ not in sys.path:
    sys.path.insert(0, RAIZ)
from pinzas import libre  # noqa: E402

KB = 1.380649e-23

AZUL, NARANJA, VERDE, ROJO, GRIS = '#2A6FB5', '#E08A2E', '#2CA02C', '#C0392B', '#7f7f7f'


def usar_estilo(fontsize=12):
    """El estilo de los notebooks de fluido/analisis (serif, fondo claro)."""
    plt.style.use('seaborn-v0_8')
    plt.rcParams.update({'font.family': 'serif', 'font.size': fontsize, 'axes.labelsize': fontsize,
                         'axes.titlesize': fontsize, 'xtick.labelsize': fontsize - 1,
                         'ytick.labelsize': fontsize - 1, 'legend.fontsize': fontsize - 2,
                         'mathtext.fontset': 'cm', 'legend.frameon': True, 'figure.dpi': 110})


# ---------------------------------------------------------------- lectura
def ultima_salida(raiz=None):
    """La carpeta de CSV más reciente (busca trayectorias.csv debajo de C:/datos_pinzas)."""
    raiz = raiz or os.environ.get('PINZAS_DATOS', 'C:/datos_pinzas')
    todas = glob.glob(os.path.join(raiz, '**', 'trayectorias.csv'), recursive=True)
    if not todas:
        raise FileNotFoundError(f'no hay ningún trayectorias.csv debajo de {raiz}: ¿corriste PROCESAR?')
    return os.path.dirname(max(todas, key=os.path.getmtime))


def cargar(carpeta):
    """Lee una salida de PROCESAR. Devuelve un dict con:
    tr (trayectorias.csv), particulas, tiempo (DataFrames), meta (grabacion.json),
    parametros (parametros.json), escala_um (µm por px, o None), carpeta, cuadro0 (imagen o None)."""
    leer = lambda n: pd.read_csv(os.path.join(carpeta, n))
    with open(os.path.join(carpeta, 'grabacion.json'), encoding='utf-8') as f:
        meta = json.load(f)
    with open(os.path.join(carpeta, 'parametros.json'), encoding='utf-8') as f:
        par = json.load(f)
    m = par.get('escala_m_por_px')
    escala = m * 1e6 if m is not None and np.isfinite(m) and m > 0 else None
    img = None
    if os.path.exists(os.path.join(carpeta, 'cuadro0.png')):
        img = plt.imread(os.path.join(carpeta, 'cuadro0.png'))
    return dict(tr=leer('trayectorias.csv'), particulas=leer('particulas.csv'), tiempo=leer('tiempo.csv'),
                meta=meta, parametros=par, escala_um=escala, carpeta=carpeta, cuadro0=img)


def matrices(datos, escala_um=None, largo_min=0):
    """Formato ancho: t [n] s, X, Y [n, nb] (µm si hay escala, si no px), ids [nb], unidad.
    Los cuadros sin medir (huecos, colisiones descartadas) quedan NaN. Los interpolados se usan
    como medidos (están marcados en tr.interpolado). largo_min: partículas con menos cuadros
    medidos se dejan afuera."""
    tr, ti = datos['tr'], datos['tiempo']
    esc = escala_um if escala_um is not None else datos['escala_um']
    n = len(ti)
    ids = [p for p, g in tr.groupby('particula') if g['x_px'].notna().sum() >= largo_min]
    X = np.full((n, len(ids)), np.nan)
    Y = np.full((n, len(ids)), np.nan)
    for k, p in enumerate(ids):
        g = tr[tr.particula == p]
        X[g.cuadro.values, k] = g.x_px.values
        Y[g.cuadro.values, k] = g.y_px.values
    if esc:
        X, Y = X * esc, Y * esc
    return ti.t_s.values.astype(float), X, Y, np.array(ids), ('µm' if esc else 'px')


# ---------------------------------------------------------------- deriva
def deriva(X, Y, minimo=2):
    """Deriva común (flujo, platina): promedio, sobre las partículas medidas en dos cuadros
    seguidos, de su desplazamiento entre esos cuadros; acumulado desde el cuadro 0.
    Devuelve dx, dy [n] y cuántas partículas entraron en cada paso.
    El browniano de cada partícula se promedia: el error de la deriva por paso es ~σ_B/√N.
    Con menos de `minimo` partículas en un paso, ese paso se interpola."""
    paso_x, paso_y = np.diff(X, axis=0), np.diff(Y, axis=0)
    ok = np.isfinite(paso_x) & np.isfinite(paso_y)
    cuantas = ok.sum(1)
    with np.errstate(invalid='ignore'):
        vx = np.where(cuantas >= minimo, np.nansum(np.where(ok, paso_x, 0), 1) / np.maximum(cuantas, 1), np.nan)
        vy = np.where(cuantas >= minimo, np.nansum(np.where(ok, paso_y, 0), 1) / np.maximum(cuantas, 1), np.nan)
    i = np.arange(len(vx))
    for v in (vx, vy):
        bien = np.isfinite(v)
        if bien.any() and not bien.all():
            v[~bien] = np.interp(i[~bien], i[bien], v[bien])
        v[~np.isfinite(v)] = 0.0
    return np.r_[0, np.cumsum(vx)], np.r_[0, np.cumsum(vy)], np.r_[0, cuantas]


def sin_deriva(X, Y):
    """Resta la deriva COMÚN (la misma para todas las partículas)."""
    dx, dy, _ = deriva(X, Y)
    return X - dx[:, None], Y - dy[:, None]


def velocidades(t, X, Y):
    """Velocidad media de cada partícula: pendiente de la recta x(t), y(t) [unidad/s] → vx, vy [nb]."""
    vx, vy = np.full(X.shape[1], np.nan), np.full(X.shape[1], np.nan)
    for k in range(X.shape[1]):
        ok = np.isfinite(X[:, k]) & np.isfinite(Y[:, k])
        if ok.sum() >= 3:
            vx[k] = np.polyfit(t[ok], X[ok, k], 1)[0]
            vy[k] = np.polyfit(t[ok], Y[ok, k], 1)[0]
    return vx, vy


def sin_deriva_por_particula(t, X, Y):
    """Resta a cada partícula SU deriva: la recta x(t), y(t) ajustada a su propia trayectoria.

    Hace falta cuando las partículas no derivan todas igual (flujo con perfil: beads a distinta
    altura van a distinta velocidad); ahí restar la deriva común deja v_i² τ² en el MSD.
    Costo: restar la velocidad media de un browniano de duración T achica el MSD en ~τ/T (y algo
    más cerca de τ ~ T). Con el ajuste de D en los primeros retardos (τ ≲ 0.3 s) y trayectorias
    de más de ~10 s, el sesgo es < 3 %; por eso conviene un LARGO_MIN alto."""
    vx, vy = velocidades(t, X, Y)
    Xs, Ys = X.copy(), Y.copy()
    for k in range(X.shape[1]):
        ok = np.isfinite(X[:, k])
        if np.isfinite(vx[k]) and ok.any():
            t0 = t[ok].mean()
            Xs[:, k] = X[:, k] - vx[k] * (t - t0)
            Ys[:, k] = Y[:, k] - vy[k] * (t - t0)
    return Xs, Ys


# ---------------------------------------------------------------- MSD y D
def msd(t, X, Y, max_lag=30):
    """MSD por partícula (pinzas.libre.msd): tau [L], msd [L, nb], n_pares [L, nb]."""
    return libre.msd(t, X, Y, max_lag=max_lag)


def msd_conjunto(m):
    """MSD promedio de todas las partículas, pesado por la cantidad de pares, y su error
    (desvío entre partículas / √N)."""
    N = m['n_pares']
    with np.errstate(invalid='ignore', divide='ignore'):
        prom = np.nansum(np.nan_to_num(m['msd']) * N, 1) / N.sum(1)
        cuantas = np.isfinite(m['msd']).sum(1)
        err = np.nanstd(m['msd'], 1) / np.sqrt(np.maximum(cuantas, 1))
    return prom, err


def columnas(m, sel):
    """El MSD de solo algunas partículas (sel: máscara o índices)."""
    return dict(m, **{k: m[k][:, sel] for k in ('msd', 'msd_x', 'msd_y', 'n_pares')})


def pendiente_loglog(m, desde=4):
    """Pendiente α de log MSD vs log τ entre el retardo `desde` y el último, por partícula.
    α ≈ 1: difusión libre. α < 0.5: el MSD se aplana → partícula CONFINADA (atrapada en la pinza,
    o pegada). α > 1: queda deriva."""
    i, j = desde, len(m['tau']) - 1
    with np.errstate(invalid='ignore', divide='ignore'):
        return (np.log(m['msd'][j]) - np.log(m['msd'][i])) / (np.log(m['tau'][j]) - np.log(m['tau'][i]))


def coeficientes_D(m, n_puntos=4):
    """D [unidad²/s] y σ por partícula: ajuste lineal de los primeros `n_puntos` retardos
    (pinzas.libre.ajustar_msd, sin blur: W = 0)."""
    filas = []
    for k in range(m['msd'].shape[1]):
        a = libre.ajustar_msd(m['tau'], m['msd'][:, k], 0.0, n_puntos=n_puntos)
        filas.append(dict(D=a['D'], sigma=a['sigma'], ordenada=a['ordenada'], n_pares_1=int(m['n_pares'][0, k])))
    return pd.DataFrame(filas)


def D_teorico(diametro_um, T=296.15, eta=0.93e-3):
    """Stokes-Einstein: D = kB T / (6 π η R) [µm²/s] (agua a 23 °C por defecto)."""
    return KB * T / (6 * np.pi * eta * diametro_um * 0.5e-6) * 1e12


# ---------------------------------------------------------------- gráficos
def _colores(n):
    return plt.cm.tab20(np.linspace(0, 1, max(n, 1)) % 1)


def graficar_trayectorias(datos, X=None, Y=None, ids=None, ax=None, fondo=True, titulo=None):
    """Trayectorias sobre el primer cuadro del video (en px, como la imagen)."""
    if X is None:
        _, X, Y, ids, _ = matrices(datos, escala_um=0)
    ax = ax or plt.subplots(figsize=(9, 7))[1]
    if fondo and datos.get('cuadro0') is not None:
        ax.imshow(datos['cuadro0'], cmap='gray', zorder=0)
    for k, c in zip(range(X.shape[1]), _colores(X.shape[1])):
        ok = np.isfinite(X[:, k])
        if not ok.any():
            continue
        ax.plot(X[:, k], Y[:, k], '-', color=c, lw=1.2)
        i0 = np.flatnonzero(ok)[0]
        ax.plot(X[i0, k], Y[i0, k], 'o', color=c, ms=4)
        ax.annotate(str(ids[k]), (X[i0, k], Y[i0, k]), color=c, fontsize=8, xytext=(4, 4), textcoords='offset points')
    ax.set_xlabel('x [px]')
    ax.set_ylabel('y [px]')
    ax.set_aspect('equal')
    if datos.get('cuadro0') is not None:
        h, w = datos['cuadro0'].shape[:2]
        ax.set_xlim(0, w)
        ax.set_ylim(h, 0)
    else:
        ax.invert_yaxis()
    ax.grid(False)
    ax.set_title(titulo or '%d trayectorias (o = donde empieza cada una)' % X.shape[1])
    return ax


def graficar_xy_t(t, X, Y, ids, unidad, axs=None):
    """x(t) e y(t) de cada partícula, relativas a su primera posición."""
    axs = axs if axs is not None else plt.subplots(2, 1, figsize=(9, 6), sharex=True)[1]
    for k, c in zip(range(X.shape[1]), _colores(X.shape[1])):
        ok = np.isfinite(X[:, k])
        if ok.any():
            i0 = np.flatnonzero(ok)[0]
            axs[0].plot(t, X[:, k] - X[i0, k], color=c, lw=1)
            axs[1].plot(t, Y[:, k] - Y[i0, k], color=c, lw=1)
    axs[0].set_ylabel('x − x₀ [%s]' % unidad)
    axs[1].set_ylabel('y − y₀ [%s]' % unidad)
    axs[1].set_xlabel('t [s]')
    return axs


def graficar_deriva(t, X, Y, unidad, ax=None):
    dx, dy, n = deriva(X, Y)
    ax = ax or plt.subplots(figsize=(9, 3.5))[1]
    ax.plot(t, dx, color=AZUL, label='x')
    ax.plot(t, dy, color=NARANJA, label='y')
    v = np.polyfit(t, dx, 1)[0], np.polyfit(t, dy, 1)[0]
    ax.set_title('deriva común: vx = %.3g, vy = %.3g %s/s (promedio de %.1f partículas por paso)'
                 % (v[0], v[1], unidad, np.mean(n[1:])))
    ax.set_xlabel('t [s]')
    ax.set_ylabel('deriva [%s]' % unidad)
    ax.legend()
    return ax, v


def graficar_velocidades(datos, t, Xp, Yp, ax=None):
    """Velocidad media de cada partícula (flecha, en px/s) dibujada en su posición media, sobre
    el primer cuadro. Si las flechas son distintas, la deriva no es común a todas."""
    vx, vy = velocidades(t, Xp, Yp)
    xm, ym = np.nanmean(Xp, 0), np.nanmean(Yp, 0)
    ax = ax or plt.subplots(figsize=(9, 7))[1]
    if datos.get('cuadro0') is not None:
        ax.imshow(datos['cuadro0'], cmap='gray')
    rap = np.hypot(vx, vy)
    q = ax.quiver(xm, ym, vx, vy, rap, cmap='plasma', angles='xy', scale_units='xy', scale=0.1, width=0.004)
    plt.colorbar(q, ax=ax, label='|v| [px/s]', shrink=0.8)
    ax.set_aspect('equal')
    if datos.get('cuadro0') is None:
        ax.invert_yaxis()
    ax.grid(False)
    ax.set_xlabel('x [px]')
    ax.set_ylabel('y [px]')
    ax.set_title('velocidad media de cada partícula (flecha = 10 s de recorrido): |v| de %.0f a %.0f px/s'
                 % (np.nanmin(rap), np.nanmax(rap)))
    return ax, (vx, vy)


def graficar_msd(m, unidad, ax=None, D=None, sigma=None, etiqueta='', color=AZUL, por_particula=True):
    """MSD de cada partícula (fino) y el promedio (grueso, con error), en log-log."""
    ax = ax or plt.subplots(figsize=(7, 5))[1]
    if por_particula:
        for k in range(m['msd'].shape[1]):
            ax.plot(m['tau'], m['msd'][:, k], '-', color=color, alpha=0.18, lw=0.8)
    prom, err = msd_conjunto(m)
    ax.errorbar(m['tau'], prom, err, fmt='o', color=color, ms=4, capsize=2, label='promedio ' + etiqueta)
    if D is not None:
        tau = m['tau'][np.isfinite(m['tau'])]
        ax.plot(tau, 4 * D * tau + 4 * (sigma or 0) ** 2, '--', color=color, lw=1.2,
                label='4Dτ + 4σ²: D = %.3g %s²/s' % (D, unidad))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('τ [s]')
    ax.set_ylabel('MSD [%s²]' % unidad)
    ax.legend()
    return ax


def graficar_histograma(t, X, Y, unidad, retardo=1, ax=None):
    """Desplazamientos Δx, Δy a `retardo` cuadros vs la gaussiana de la misma varianza."""
    dt0 = libre.periodo_nominal(t)
    d = libre.desplazamientos(t, X, Y, retardo, dt0)
    ax = ax or plt.subplots(figsize=(7, 4.5))[1]
    s = np.std(d)
    ax.hist(d, bins=51, density=True, color=AZUL, alpha=0.6, label='Δx y Δy (%d pares)' % len(d))
    u = np.linspace(-4 * s, 4 * s, 300)
    ax.plot(u, np.exp(-u ** 2 / (2 * s ** 2)) / np.sqrt(2 * np.pi * s ** 2), color=NARANJA,
            label='gaussiana, σ = %.3g %s' % (s, unidad))
    ax.set_yscale('log')
    ax.set_ylim(1e-3 / s, None)
    ax.set_xlabel('desplazamiento en %d cuadro(s) [%s]' % (retardo, unidad))
    ax.set_ylabel('densidad')
    ax.legend()
    return ax


# ================================================================ comparación de formas de sacar la deriva
# (lo usa "BROWNIANO comparación.ipynb"). Todas las funciones aceptan huecos (NaN) y tiempos reales.
def tiempo_regular(t, salto_grande=2.5):
    """Tiempos de un video cuyo contenedor repite cuadros (MP4 a 30 cuadros/s de una cámara más lenta).

    El importador descarta los cuadros repetidos y a cada cuadro nuevo le deja el tiempo del contenedor,
    que está cuantizado: los intervalos salen de 1 o 2 períodos del contenedor aunque la cámara vaya
    pareja. Acá se supone que entre dos saltos grandes (> salto_grande × el intervalo más corto) la
    cámara fue pareja: dentro de cada tramo t = a + índice × período, con el período medio de todos los
    intervalos normales. Los saltos grandes (cuadros que faltan) se conservan.
    Devuelve t_regular [n] y el período [s]."""
    t = np.asarray(t, float)
    dt = np.diff(t)
    corto = np.min(dt[dt > 0])
    normal = dt <= salto_grande * corto
    periodo = dt[normal].mean()
    cortes = np.flatnonzero(~normal) + 1
    tr_ = np.empty_like(t)
    for i0, i1 in zip(np.r_[0, cortes], np.r_[cortes, len(t)]):
        i = np.arange(i0, i1)
        tr_[i] = np.mean(t[i] - (i - i0) * periodo) + (i - i0) * periodo
    return tr_, float(periodo)


def tendencia_local(t, A, ventana_s):
    """Tendencia lenta de cada columna de A [n, nb]: en cada instante, la recta ajustada a los puntos
    medidos de esa partícula que están a menos de ventana_s/2, evaluada en ese instante.

    Es una media móvil que no se sesga en los bordes ni con huecos: lejos de los bordes y sin huecos
    coincide con la media móvil centrada de ancho ventana_s; una deriva de velocidad constante la
    devuelve exacta en toda la trayectoria. Con menos de 3 puntos en la ventana usa su promedio."""
    t = np.asarray(t, float)
    T = np.full(A.shape, np.nan)
    for k in range(A.shape[1]):
        ok = np.isfinite(A[:, k])
        if ok.sum() < 2:
            continue
        tv = t[ok] - t[ok].mean()
        xv = A[ok, k]
        lo = np.searchsorted(tv, tv - ventana_s / 2, 'left')
        hi = np.searchsorted(tv, tv + ventana_s / 2, 'right')
        S1, St, Stt, Sx, Stx = (np.r_[0, np.cumsum(v)] for v in (np.ones_like(tv), tv, tv * tv, xv, tv * xv))
        n, st, stt, sx, stx = (S[hi] - S[lo] for S in (S1, St, Stt, Sx, Stx))
        den = n * stt - st * st
        with np.errstate(invalid='ignore', divide='ignore'):
            b = np.where((n >= 3) & (den > 1e-12), (n * stx - st * sx) / den, 0.0)
        T[ok, k] = (sx - b * st) / n + b * tv
    return T


def sin_deriva_media_movil(t, X, Y, ventana_s):
    """Pasa-altos por partícula: a cada trayectoria se le resta su tendencia local (`tendencia_local`,
    media móvil de ancho ventana_s). Saca cualquier deriva que cambie más lento que la ventana, sea o
    no común a todas; el costo es que también saca la parte lenta del browniano (ver `simular`)."""
    return X - tendencia_local(t, X, ventana_s), Y - tendencia_local(t, Y, ventana_s)


def inverso_n_deriva(X, Y):
    """<1/N> por partícula: promedio, sobre sus pasos, de 1/(cuántas partículas entraron en la deriva
    común de ese paso). Restar la deriva común (que incluye a la propia partícula) achica su MSD en
    (1 − <1/N>) si todas tienen el mismo D: D = D_medido / (1 − <1/N>)."""
    ok = np.isfinite(np.diff(X, axis=0)) & np.isfinite(np.diff(Y, axis=0))
    cuantas = ok.sum(1)
    inv = 1.0 / np.maximum(cuantas, 1)
    return np.array([inv[ok[:, k]].mean() if ok[:, k].any() else np.nan for k in range(X.shape[1])])


def msd_pares(t, X, Y, max_lag=30, solape_min=200):
    """MSD de la SEPARACIÓN entre dos partículas, para cada par con al menos `solape_min` cuadros en común.

    r_ij = r_i − r_j no tiene la deriva común (aunque cambie en el tiempo) y no se le resta nada, así
    que no pierde browniano: MSD_ij(τ) = 4 (D_i + D_j) τ + 4 (σ_i² + σ_j²). No cancela la deriva que
    es distinta en cada partícula: queda (v_i − v_j)² τ².
    Devuelve el dict de `msd` (una columna por par) y pares [np, 2] (índices de columnas de X)."""
    pares = [(i, j) for i in range(X.shape[1]) for j in range(i + 1, X.shape[1])
             if (np.isfinite(X[:, i]) & np.isfinite(X[:, j])).sum() >= solape_min]
    if not pares:
        return None, np.zeros((0, 2), int)
    pares = np.array(pares)
    m = libre.msd(t, X[:, pares[:, 0]] - X[:, pares[:, 1]], Y[:, pares[:, 0]] - Y[:, pares[:, 1]], max_lag=max_lag)
    return m, pares


def D_de_pares(m, pares, nb, n_puntos=4):
    """D del conjunto a partir de los pares: promedio de (D_i + D_j)/2 pesado por la cantidad de datos
    de cada par, y su error por jackknife sobre las PARTÍCULAS (los pares no son independientes: cada
    partícula está en varios). También D por partícula, resolviendo D_i + D_j = S_ij por cuadrados
    mínimos (NaN si el sistema no lo determina).
    Devuelve dict: D, err_D, ordenada (por partícula), S [np], D_i [nb], n_particulas, n_pares."""
    aj = [libre.ajustar_msd(m['tau'], m['msd'][:, k], 0.0, n_puntos=n_puntos) for k in range(len(pares))]
    S = np.array([a['D'] for a in aj])
    a0 = np.array([a['ordenada'] for a in aj])
    w = m['n_pares'][0].astype(float)
    bien = np.isfinite(S)
    media = lambda sel: np.sum(S[sel] * w[sel]) / np.sum(w[sel]) / 2
    usadas = np.unique(pares[bien])
    jk = np.array([media(bien & (pares != p).all(1)) for p in usadas if (bien & (pares != p).all(1)).any()])
    err = np.sqrt((len(jk) - 1) / len(jk) * np.sum((jk - jk.mean()) ** 2)) if len(jk) > 1 else np.nan
    A = np.zeros((bien.sum(), nb))
    A[np.arange(bien.sum()), pares[bien, 0]] = 1
    A[np.arange(bien.sum()), pares[bien, 1]] = 1
    sw = np.sqrt(w[bien])
    Di = np.full(nb, np.nan)
    cols = np.flatnonzero(A.any(0))
    if len(cols) and np.linalg.matrix_rank(A[:, cols]) == len(cols):
        Di[cols] = np.linalg.lstsq(A[:, cols] * sw[:, None], S[bien] * sw, rcond=None)[0]
    return dict(D=media(bien), err_D=err, ordenada=np.sum(a0[bien] * w[bien]) / np.sum(w[bien]) / 2, S=S, D_i=Di,
                n_particulas=len(usadas), n_pares=int(bien.sum()))


def ajustar_msd_deriva(tau, msd_, n_puntos=None):
    """Ajuste MSD(τ) = 4 D τ + v² τ² + 4 σ²  (browniano + deriva de velocidad constante v).

    Cuadrados mínimos lineales en (1, τ, τ²) sobre los primeros n_puntos retardos (None = todos), con
    peso 1/τ (el error del MSD crece con el retardo). No resta nada a las trayectorias: el browniano
    queda entero. v² puede salir negativo por estadística (se informa igual; v = NaN).
    Devuelve dict: D, v2, v, sigma2, ordenada."""
    tau, msd_ = np.asarray(tau, float), np.asarray(msd_, float)
    ok = np.isfinite(tau) & np.isfinite(msd_)
    tv, mv = tau[ok][:n_puntos], msd_[ok][:n_puntos]
    if len(tv) < 4:
        return dict(D=np.nan, v2=np.nan, v=np.nan, sigma2=np.nan, ordenada=np.nan)
    c2, c1, c0 = np.polyfit(tv, mv, 2, w=1 / tv)
    return dict(D=c1 / 4, v2=c2, v=np.sqrt(c2) if c2 >= 0 else np.nan, sigma2=c0 / 4, ordenada=c0)


def correlacion_pasos(t, X, Y):
    """Correlación entre un paso y el siguiente, por partícula (promedio de x e y), usando solo los
    pasos de un período nominal seguidos. Devuelve r [nb] y cuántos pares de pasos entraron [nb].
    Browniano puro: 0. Ruido de localización: negativa. Exposición larga o deriva: positiva."""
    uno = np.abs(np.diff(t) / libre.periodo_nominal(t) - 1) < 0.25
    r, n = np.full(X.shape[1], np.nan), np.zeros(X.shape[1], int)
    for k in range(X.shape[1]):
        rr = []
        for A in (X, Y):
            d = np.where(uno, np.diff(A[:, k]), np.nan)
            ok = np.isfinite(d[:-1]) & np.isfinite(d[1:])
            if ok.sum() > 10:
                rr.append(np.corrcoef(d[:-1][ok], d[1:][ok])[0, 1])
                n[k] = ok.sum()
        if rr:
            r[k] = np.mean(rr)
    return r, n


def pasos_por_particula(t, X, Y):
    """Pasos de un período de cada partícula: lista de (dx, dy), solo los pasos medidos."""
    uno = np.abs(np.diff(t) / libre.periodo_nominal(t) - 1) < 0.25
    sal = []
    for k in range(X.shape[1]):
        dx, dy = np.diff(X[:, k]), np.diff(Y[:, k])
        ok = uno & np.isfinite(dx) & np.isfinite(dy)
        sal.append((dx[ok], dy[ok]))
    return sal


def resumen_pasos(t, X, Y):
    """Números del histograma de desplazamientos a un período, juntando las partículas después de
    dividir los pasos de cada una (y de cada eje) por su desvío estándar:
    media_sobre_sigma (media del paso / desvío, promedio cuadrático entre partículas y ejes: mide la
    deriva que queda), curtosis (en exceso; 0 = gaussiana) y err_curtosis, n_pasos,
    anisotropia (σ_y/σ_x, media entre partículas) y err_anisotropia."""
    from scipy import stats
    z, ms, an = [], [], []
    for dx, dy in pasos_por_particula(t, X, Y):
        if len(dx) < 20:
            continue
        for d in (dx, dy):
            z.append((d - d.mean()) / d.std())
            ms.append(d.mean() / d.std())
        an.append(dy.std() / dx.std())
    if not z:
        return dict(media_sobre_sigma=np.nan, curtosis=np.nan, err_curtosis=np.nan, n_pasos=0,
                    anisotropia=np.nan, err_anisotropia=np.nan)
    z = np.concatenate(z)
    return dict(media_sobre_sigma=float(np.sqrt(np.mean(np.square(ms)))), curtosis=float(stats.kurtosis(z)),
                err_curtosis=float(np.sqrt(24 / len(z))), n_pasos=len(z), anisotropia=float(np.mean(an)),
                err_anisotropia=float(np.std(an) / np.sqrt(len(an))))


def simular(t, mascara, D, sigma=0.0, v=None, deriva_comun=None, rng=None):
    """Browniano 2D sintético con los MISMOS tiempos y los mismos huecos que una medición.

    t [n] s; mascara [n, nb] (True = medido); D escalar o [nb] (unidad²/s); sigma: ruido de
    localización por eje (escalar o [nb]); v [nb, 2]: velocidad constante de cada partícula;
    deriva_comun [n, 2]: un corrimiento igual para todas. Devuelve X, Y [n, nb] con NaN donde
    mascara es False. Sirve para medir cuánto browniano se lleva cada forma de sacar la deriva."""
    rng = rng or np.random.default_rng()
    n, nb = mascara.shape
    D = np.broadcast_to(np.asarray(D, float), (nb,))
    t = np.asarray(t, float)
    dt = np.diff(t)
    sal = []
    for eje in range(2):
        pasos = rng.standard_normal((n - 1, nb)) * np.sqrt(2 * D[None, :] * dt[:, None])
        A = np.vstack([np.zeros((1, nb)), np.cumsum(pasos, 0)])
        if v is not None:
            A = A + np.asarray(v)[None, :, eje] * (t - t[0])[:, None]
        if deriva_comun is not None:
            A = A + np.asarray(deriva_comun)[:, eje][:, None]
        A = A + rng.standard_normal((n, nb)) * np.broadcast_to(np.asarray(sigma, float), (nb,))[None, :]
        sal.append(np.where(mascara, A, np.nan))
    return sal[0], sal[1]


METODOS_DERIVA = ('sin restar', 'deriva común', 'recta por partícula', 'media móvil', 'pares', 'ajuste con deriva')


def analizar_deriva(t, X, Y, metodo, ventana_s=None, max_lag=30, n_ajuste=4, solape_min=200, corregir=None):
    """Aplica una forma de sacar la deriva y devuelve los números que se comparan.

    metodo: uno de METODOS_DERIVA.
      'sin restar'           las posiciones como salen del tracking
      'deriva común'         resta el promedio de los pasos de todas (`sin_deriva`); D se corrige por
                             1/(1 − <1/N>) (`inverso_n_deriva`)
      'recta por partícula'  `sin_deriva_por_particula`
      'media móvil'          `sin_deriva_media_movil` con ventana_s
      'pares'                MSD de la separación de cada par (`msd_pares`, `D_de_pares`)
      'ajuste con deriva'    sin restar nada, MSD de cada partícula ajustado con 4Dτ + v²τ² + 4σ²
                             sobre los max_lag retardos (`ajustar_msd_deriva`)
    corregir: factor por el que se divide D (el que da la simulación gemela: D recuperado / D verdadero).
    Devuelve dict:
      D_i [nb]      D por partícula (en 'pares' sale del sistema D_i + D_j = S_ij; NaN si no se puede resolver)
      D, err_D, n   media de D_i y desvío entre partículas / √n. En 'pares', si el sistema no se puede
                    resolver, D es el del conjunto de pares (D_pares: promedio de (D_i + D_j)/2 pesado
                    por pesos_pares) y err_D su jackknife sobre partículas; D_es_del_conjunto lo avisa
      D_conjunto    D del ajuste del MSD promedio (pesado por la cantidad de datos de cada partícula)
      ordenada      la del MSD promedio (por partícula)
      alfa, err_alfa  pendiente log-log del MSD promedio en todos los retardos (en 'ajuste con deriva',
                    del MSD menos el término v²τ² ajustado) y su error por jackknife sobre partículas (o pares)
      tau, msd, err_msd   el MSD promedio y su error
      corr, err_corr      correlación entre pasos seguidos (media pesada)
      los de `resumen_pasos`, y Xs, Ys (las posiciones corregidas)."""
    if metodo not in METODOS_DERIVA:
        raise ValueError('metodo: uno de %s' % (METODOS_DERIVA,))
    f = 1.0 if corregir is None else float(corregir)
    nb = X.shape[1]
    if metodo == 'deriva común':
        Xs, Ys = sin_deriva(X, Y)
    elif metodo == 'recta por partícula':
        Xs, Ys = sin_deriva_por_particula(t, X, Y)
    elif metodo == 'media móvil':
        Xs, Ys = sin_deriva_media_movil(t, X, Y, ventana_s)
    else:
        Xs, Ys = X, Y
    r = dict(metodo=metodo, Xs=Xs, Ys=Ys)
    if metodo == 'pares':
        m, pares = msd_pares(t, X, Y, max_lag, solape_min)
        if m is None:
            return dict(r, D=np.nan, err_D=np.nan, n=0, D_i=np.full(nb, np.nan), D_conjunto=np.nan, alfa=np.nan,
                        err_alfa=np.nan, ordenada=np.nan, corr=np.nan, err_corr=np.nan, tau=None, msd=None, err_msd=None,
                        media_sobre_sigma=np.nan, curtosis=np.nan, err_curtosis=np.nan, n_pasos=0,
                        anisotropia=np.nan, err_anisotropia=np.nan)
        dp = D_de_pares(m, pares, nb, n_ajuste)
        prom, err = msd_conjunto(m)
        prom, err = prom / 2, err / 2                              # por partícula: MSD_ij / 2
        Xr, Yr = X[:, pares[:, 0]] - X[:, pares[:, 1]], Y[:, pares[:, 0]] - Y[:, pares[:, 1]]
        c, nc = correlacion_pasos(t, Xr, Yr)
        Di = dp['D_i'] / f
        resuelto = np.isfinite(Di).sum() == dp['n_particulas'] and dp['n_particulas'] > 1
        if resuelto:
            Dm, eD = float(np.nanmean(Di)), float(np.nanstd(Di, ddof=1) / np.sqrt(dp['n_particulas']))
        else:
            Dm, eD = dp['D'] / f, dp['err_D'] / f
        r.update(D_i=Di, D=Dm, err_D=eD, n=dp['n_particulas'], n_pares=dp['n_pares'], D_pares=dp['D'] / f,
                 err_D_pares=dp['err_D'] / f, D_es_del_conjunto=not resuelto, pares=pares,
                 pesos_pares=m['n_pares'][0].astype(float), **resumen_pasos(t, Xr, Yr))
        g = 1.0
    else:
        m = msd(t, Xs, Ys, max_lag)
        prom, err = msd_conjunto(m)
        if metodo == 'ajuste con deriva':
            aj = [ajustar_msd_deriva(m['tau'], m['msd'][:, k]) for k in range(nb)]
            Di = np.array([a['D'] for a in aj])
            r['v2_i'] = np.array([a['v2'] for a in aj])
        else:
            Di = coeficientes_D(m, n_puntos=n_ajuste).D.values
        g = 1.0
        if metodo == 'deriva común':
            inv = inverso_n_deriva(X, Y)
            Di = Di / (1 - inv)
            w = m['n_pares'][0]
            g = 1 - float(np.nansum(inv * w) / np.sum(w[np.isfinite(inv)]))      # para el MSD promedio
        Di = Di / f
        bien = np.isfinite(Di)
        c, nc = correlacion_pasos(t, Xs, Ys)
        r.update(D_i=Di, D=float(np.mean(Di[bien])) if bien.any() else np.nan,
                 err_D=float(np.std(Di[bien], ddof=1) / np.sqrt(bien.sum())) if bien.sum() > 1 else np.nan,
                 n=int(bien.sum()), **resumen_pasos(t, Xs, Ys))
    tau = m['tau']
    if metodo == 'ajuste con deriva':
        a = ajustar_msd_deriva(tau, prom)
        r['v2_conjunto'] = a['v2']
    else:
        a = libre.ajustar_msd(tau, prom, 0.0, n_puntos=n_ajuste)

    def pendiente(p):
        if metodo == 'ajuste con deriva':
            p = p - ajustar_msd_deriva(tau, p)['v2'] * tau ** 2
        ok = np.isfinite(tau) & np.isfinite(p) & (p > 0)
        return float(np.polyfit(np.log(tau[ok]), np.log(p[ok]), 1)[0]) if ok.sum() > 2 else np.nan

    N, Mk = m['n_pares'], np.nan_to_num(m['msd'])
    escala_ = 0.5 if metodo == 'pares' else 1.0
    jk = []
    for k in range(Mk.shape[1]):                       # jackknife: sacar una partícula (o un par) por vez
        sel = np.arange(Mk.shape[1]) != k
        with np.errstate(invalid='ignore', divide='ignore'):
            jk.append(pendiente(escala_ * (Mk[:, sel] * N[:, sel]).sum(1) / N[:, sel].sum(1)))
    jk = np.array(jk)
    jk = jk[np.isfinite(jk)]
    bien = np.isfinite(c)
    r.update(tau=tau, msd=prom, err_msd=err, D_conjunto=a['D'] / f / g, ordenada=a['ordenada'], alfa=pendiente(prom),
             err_alfa=float(np.sqrt((len(jk) - 1) / len(jk) * np.sum((jk - jk.mean()) ** 2))) if len(jk) > 1 else np.nan,
             corr=float(np.average(c[bien], weights=nc[bien])) if bien.any() else np.nan,
             err_corr=float(1 / np.sqrt(nc[bien].sum() / 2)) if bien.any() else np.nan)
    return r
