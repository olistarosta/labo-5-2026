"""LIBRE — browniano libre: MSD con los tiempos reales -> D, sigma y eta.

Física (ver el encargo, "Browniano libre"):
  MSD(tau) = < |r(t+tau) - r(t)|² >_t        promediado sobre TODOS los orígenes t
  MSD(tau) = 4 D (tau - W/3) + 4 sigma²        Savin & Doyle, Biophys. J. 88, 623 (2005)
             (2 ejes; W = exposición; sigma = ruido de localización POR EJE; vale para tau >= W)
  eta = kB T / (6 pi R D)                      (con Faxén si se conoce h)

NUNCA se usa el promedio acumulado de la ec. (30) del manual del kit: con r² = 4 D t ese
promedio crece como 2 D t y da la mitad de la pendiente (eta al doble).

Convenciones (CONTRATO.md): SI, t [n] s reales (con huecos), x, y [n, nb] m, NaN = no medido.
Los retardos se agrupan por retardo REAL: el par (i, j) cae en el retardo m = round((t_j - t_i)/dt0),
con dt0 el período nominal. Un par que cruza un cuadro perdido usa su retardo verdadero.
"""
import numpy as np

from .fisica import KB, eta_desde_D, faxen_paralelo


# ---------------------------------------------------------------- tiempos
def periodo_nominal(t, tramos=None):
    """Período de cuadro dt0 [s]: mediana de los dt dentro de los tramos sin huecos.

    t: [n] s. tramos: [m, 2] (i0, i1) inclusive, o None (mediana de todos los dt; un hueco
    aislado no la mueve).
    """
    t = np.asarray(t, float)
    dt = np.diff(t)
    if tramos is None or len(tramos) == 0:
        return float(np.median(dt))
    sel = [dt[i0:i1] for i0, i1 in np.asarray(tramos, int) if i1 > i0]
    return float(np.median(np.concatenate(sel))) if sel else float(np.median(dt))


def tramos_desde_t(t, tolerancia=1.5):
    """Tramos sin huecos deducidos de los tiempos, si TIEMPO no los dio.

    Corta donde dt > tolerancia * mediana(dt). Devuelve [m, 2] (i0, i1) inclusive.
    """
    t = np.asarray(t, float)
    dt = np.diff(t)
    cortes = np.flatnonzero(dt > tolerancia * np.median(dt))
    i0 = np.concatenate([[0], cortes + 1])
    i1 = np.concatenate([cortes, [len(t) - 1]])
    return np.stack([i0, i1], axis=1)


def pares_por_retardo(t, X, max_lag, dt0, bloque=None, n_bloques=1):
    """Sumas sobre todos los pares (i, j>i), agrupadas por retardo real y (opcional) por bloque.

    t: [n] s.  X: [n, k] (k series, p. ej. partículas o un eje), NaN = no medido.
    max_lag: retardo máximo en períodos nominales (m = 1..max_lag).  dt0: período nominal [s].
    bloque: [n] int (bloque de cada cuadro; el par pertenece al bloque de i) o None.

    El par (i, j) cae en m = round((t_j - t_i)/dt0): los retardos son los REALES, no j - i.
    Devuelve dict de arrays [n_bloques, max_lag + 1, k] (el índice m = 0 queda vacío):
      N (pares), s_tau (suma de retardos reales), s_d2 (suma de (X_j - X_i)²),
      s_ab (suma de X_i X_j), s_a (suma de X_i), s_b (suma de X_j).
    """
    t = np.asarray(t, float)
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]
    n, k = X.shape
    L = int(max_lag)
    bl = np.zeros(n, int) if bloque is None else np.asarray(bloque, int)
    forma = (n_bloques * (L + 1),)
    claves = ('N', 's_tau', 's_d2', 's_ab', 's_a', 's_b')
    out = {c: np.zeros((n_bloques, L + 1, k)) for c in claves}
    fin = np.isfinite(X)
    for d in range(1, min(L, n - 1) + 1):        # índice j - i: un retardo real m >= j - i
        i = np.arange(n - d)
        j = i + d
        m = np.rint((t[j] - t[i]) / dt0).astype(np.int64)
        ok_t = (m >= 1) & (m <= L)
        if not ok_t.any():
            continue
        for c in range(k):
            ok = ok_t & fin[i, c] & fin[j, c]
            if not ok.any():
                continue
            ii, jj = i[ok], j[ok]
            idx = bl[ii] * (L + 1) + m[ok]
            a, b = X[ii, c], X[jj, c]
            for clave, w in (('N', None), ('s_tau', t[jj] - t[ii]), ('s_d2', (b - a) ** 2),
                             ('s_ab', a * b), ('s_a', a), ('s_b', b)):
                out[clave][:, :, c] += np.bincount(idx, w, minlength=forma[0]).reshape(n_bloques, L + 1)
    return out


# ---------------------------------------------------------------- MSD
def msd(t, x, y, tramos=None, max_lag=30, dt0=None):
    """MSD 2D por partícula promediado sobre todos los orígenes, con retardos reales.

    t: [n] s. x, y: [n, nb] m (o [n] para una partícula). NaN = no medido; si falta x o y
    el cuadro no se usa. tramos: solo para estimar dt0 (los retardos ya salen de t).
    max_lag: retardo máximo en períodos nominales.

    Devuelve dict:
      tau [L] s (retardo real medio de cada grupo, NaN si no hubo pares),
      msd [L, nb] m², msd_x, msd_y [L, nb] m², n_pares [L, nb], dt0 s.
    """
    x = np.atleast_2d(np.asarray(x, float).T).T
    y = np.atleast_2d(np.asarray(y, float).T).T
    malo = ~(np.isfinite(x) & np.isfinite(y))
    x = np.where(malo, np.nan, x)
    y = np.where(malo, np.nan, y)
    if dt0 is None:
        dt0 = periodo_nominal(t, tramos)
    nb = x.shape[1]
    p = pares_por_retardo(t, np.concatenate([x, y], axis=1), max_lag, dt0)
    N = p['N'][0, 1:, :nb]
    with np.errstate(invalid='ignore', divide='ignore'):
        mx = p['s_d2'][0, 1:, :nb] / N
        my = p['s_d2'][0, 1:, nb:] / N
        tau = p['s_tau'][0, 1:, :nb].sum(1) / N.sum(1)
    mx[N == 0] = np.nan
    my[N == 0] = np.nan
    return {'tau': tau, 'msd': mx + my, 'msd_x': mx, 'msd_y': my, 'n_pares': N, 'dt0': dt0}


def n_puntos_michalet(x_reducido):
    """Cuántos retardos ajustar, según el ruido reducido x = |sigma_ef²| / (D dt0).

    Michalet, PRE 82, 041914 (2010): con x << 1 lo óptimo es ajustar 2 puntos; el número
    crece con x. Los cortes de acá salen de una simulación propia (N = 1000 cuadros, 300
    repeticiones, ajuste lineal sin pesos) que reproduce esa tendencia:
      x < 0.5 -> 2;  x < 2 -> 3;  x < 10 -> 6;  si no -> 10.
    sigma_ef² = sigma² - D W / 3 (el blur actúa como un ruido negativo en la ordenada).
    """
    if not np.isfinite(x_reducido):
        return 4
    x = abs(x_reducido)
    return 2 if x < 0.5 else 3 if x < 2 else 6 if x < 10 else 10


def ajustar_msd(tau, msd_, W, n_puntos=None, dt0=None):
    """Ajuste lineal MSD(tau) = 4 D (tau - W/3) + 4 sigma² sobre los primeros retardos.

    tau: [L] s; msd_: [L] m² (2D); W: exposición [s].
    n_puntos: retardos a usar; None = elige con `n_puntos_michalet` a partir de un ajuste
    preliminar de 4 puntos (hace falta dt0, si no se toma tau[0]).
    Mínimos cuadrados sin pesos (los puntos están correlacionados: el error de D NO sale de
    acá sino de la dispersión entre partículas).

    Devuelve dict: D [m²/s], sigma2 [m²] por eje (puede salir negativo por estadística:
    se informa igual), sigma [m] (NaN si sigma2 < 0), n_puntos, ordenada [m²], pendiente [m²/s].
    """
    tau = np.asarray(tau, float)
    msd_ = np.asarray(msd_, float)
    ok = np.isfinite(tau) & np.isfinite(msd_)
    tv, mv = tau[ok], msd_[ok]
    nan = {'D': np.nan, 'sigma2': np.nan, 'sigma': np.nan, 'n_puntos': 0,
           'ordenada': np.nan, 'pendiente': np.nan}
    if len(tv) < 2:
        return nan
    if n_puntos is None:
        b, a = np.polyfit(tv[:4], mv[:4], 1)
        D0 = b / 4
        s2_ef = a / 4                     # = sigma² - D W / 3
        paso = dt0 if dt0 is not None else tv[0]
        n_puntos = n_puntos_michalet(s2_ef / (D0 * paso)) if D0 > 0 else 4
    n_puntos = int(min(max(n_puntos, 2), len(tv)))
    b, a = np.polyfit(tv[:n_puntos], mv[:n_puntos], 1)
    D = b / 4
    sigma2 = a / 4 + D * W / 3
    return {'D': D, 'sigma2': sigma2, 'sigma': np.sqrt(sigma2) if sigma2 >= 0 else np.nan,
            'n_puntos': n_puntos, 'ordenada': a, 'pendiente': b}


def modelo_msd(tau, D, sigma2, W):
    """MSD(tau) = 4 D (tau - W/3) + 4 sigma² [m²] (Savin & Doyle 2005), para graficar."""
    return 4 * D * (np.asarray(tau, float) - W / 3) + 4 * sigma2


# ---------------------------------------------------------------- histograma
def desplazamientos(t, x, y, retardo, dt0):
    """Desplazamientos Δx, Δy [m] de todos los pares con retardo real `retardo` períodos.

    t: [n] s; x, y: [n, nb] m. Junta todas las partículas y los dos ejes en un vector.
    """
    x = np.atleast_2d(np.asarray(x, float).T).T
    y = np.atleast_2d(np.asarray(y, float).T).T
    t = np.asarray(t, float)
    out = []
    for d in range(1, retardo + 1):
        i = np.arange(len(t) - d)
        j = i + d
        sel = np.rint((t[j] - t[i]) / dt0).astype(int) == retardo
        for X in (x, y):
            dX = X[j[sel]] - X[i[sel]]
            out.append(dX[np.isfinite(dX)])
    return np.concatenate(out) if out else np.zeros(0)


def chequeo_gaussiano(d, n_bins=41):
    """Histograma de desplazamientos y su comparación con una gaussiana de la misma varianza.

    d: [k] m. Devuelve dict: bordes [n_bins+1] m, cuentas [n_bins], gauss [n_bins] (cuentas
    esperadas con media 0 y varianza <d²>), curtosis_exceso (0 si gaussiano; error ~ sqrt(24/k)),
    p_ks (Kolmogorov-Smirnov contra esa normal; supone pares independientes, que a retardo 1
    casi se cumple: el ruido de localización correlaciona levemente pasos vecinos).
    """
    from scipy import stats
    d = np.asarray(d, float)
    d = d[np.isfinite(d)]
    if len(d) < 10:
        return {'bordes': np.zeros(0), 'cuentas': np.zeros(0), 'gauss': np.zeros(0),
                'curtosis_exceso': np.nan, 'err_curtosis': np.nan, 'p_ks': np.nan, 'n': len(d)}
    s = np.sqrt(np.mean(d ** 2))
    bordes = np.linspace(-4 * s, 4 * s, n_bins + 1)
    cuentas, _ = np.histogram(d, bordes)
    gauss = len(d) * np.diff(stats.norm.cdf(bordes, 0, s))
    return {'bordes': bordes, 'cuentas': cuentas, 'gauss': gauss,
            'curtosis_exceso': float(stats.kurtosis(d)), 'err_curtosis': float(np.sqrt(24 / len(d))),
            'p_ks': float(stats.kstest(d / s, 'norm').pvalue), 'n': len(d)}


# ---------------------------------------------------------------- análisis completo
def analizar(t, x, y, W, T, R, h=None, tramos=None, max_lag=30, n_puntos=None,
             retardo_hist=1, n_deriva=None):
    """Análisis LIBRE completo: MSD -> D por partícula y promedio, sigma, eta, histograma.

    t: [n] s reales. x, y: [n, nb] m, ya sin deriva (DERIVA). W: exposición [s]. T: temperatura [K].
    R: radio de las beads [m]. h: distancia del centro a la pared [m] o None (sin Faxén).
    tramos: [m, 2] de TIEMPO (solo para dt0). max_lag: retardos a calcular (períodos).
    n_puntos: retardos del ajuste (None = regla de Michalet, ver `n_puntos_michalet`).
    retardo_hist: retardo (en períodos) del histograma de desplazamientos.
    n_deriva: si DERIVA restó el desplazamiento medio de n_deriva partículas (incluida esta),
      cada MSD queda multiplicado por (1 - 1/n_deriva) (partículas independientes e iguales);
      se corrige multiplicando por n_deriva / (n_deriva - 1). None = no se corrige.

    Devuelve dict:
      D [nb], sigma [nb], sigma2 [nb], n_puntos [nb]      por partícula
      D_medio, err_D (std entre partículas / sqrt(nb); NaN si nb < 2), nb_validas
      eta, err_eta [Pa·s] sin Faxén;  eta_faxen, err_eta_faxen (NaN si h es None)
      tau [L] s, msd [L, nb] m², msd_ajuste [L, nb] m², msd_medio [L], dt0
      hist: dict de `chequeo_gaussiano` al retardo `retardo_hist`
      avisos: lista de str
    """
    avisos = []
    r = msd(t, x, y, tramos, max_lag)
    M = r['msd']
    if n_deriva is not None and n_deriva >= 2:
        M = M * n_deriva / (n_deriva - 1)
    nb = M.shape[1]
    fits = [ajustar_msd(r['tau'], M[:, k], W, n_puntos, r['dt0']) for k in range(nb)]
    D = np.array([f['D'] for f in fits])
    validas = np.isfinite(D) & (D > 0)
    if (~validas).any():
        avisos.append(f'{int((~validas).sum())} partícula(s) sin D válido (NaN)')
    Dv = D[validas]
    D_medio = float(np.mean(Dv)) if len(Dv) else np.nan
    err_D = float(np.std(Dv, ddof=1) / np.sqrt(len(Dv))) if len(Dv) >= 2 else np.nan
    if len(Dv) < 2:
        avisos.append('menos de 2 partículas: el error de D (dispersión entre partículas) es NaN')
    eta = eta_desde_D(D_medio, T, R)
    err_eta = eta * err_D / D_medio
    if h is not None:
        f = float(faxen_paralelo(R, h))
        eta_f, err_eta_f = eta / f, err_eta / f
    else:
        eta_f = err_eta_f = np.nan
    ajuste = np.stack([modelo_msd(r['tau'], fi['D'], fi['sigma2'], W) for fi in fits], axis=1)
    d = desplazamientos(t, x, y, retardo_hist, r['dt0'])
    return {
        'D': D, 'sigma': np.array([f['sigma'] for f in fits]),
        'sigma2': np.array([f['sigma2'] for f in fits]),
        'n_puntos': np.array([f['n_puntos'] for f in fits]),
        'D_medio': D_medio, 'err_D': err_D, 'nb_validas': int(validas.sum()),
        'eta': float(eta), 'err_eta': float(err_eta), 'eta_faxen': float(eta_f),
        'err_eta_faxen': float(err_eta_f),
        'tau': r['tau'], 'msd': M, 'msd_ajuste': ajuste, 'msd_medio': np.nanmean(M, axis=1),
        'dt0': r['dt0'], 'hist': chequeo_gaussiano(d), 'avisos': avisos,
    }
