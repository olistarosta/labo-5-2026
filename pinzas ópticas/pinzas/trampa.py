"""TRAMPA — rigidez kappa de la trampa por cinco métodos, por eje, con correcciones.

Modelo (encargo, "Partícula atrapada"), por eje:
    gamma dx/dt = -kappa x + sqrt(2 gamma kB T) xi(t),   tau_c = gamma/kappa,  f_c = kappa/(2 pi gamma)
    v = <x²> = kB T / kappa                                (varianza instantánea)

Lo que mide la cámara en el cuadro i es la posición PROMEDIADA durante la exposición W más el
ruido de localización (varianza s2, blanco). Con alpha = W / tau_c y W <= dt (exposición <=
período), las fórmulas exactas para ese proceso son (ver fisica.factor_blur_*):
    Var(x_i)                    = v S(alpha) + s2                     Wong & Halvorsen 2006
    Cov(x_i, x_j), t_j - t_i = tau >= W  = v B(alpha) e^(-tau/tau_c)   (el ruido no entra)
    MSD(tau >= W)               = 2 s2 + 2 v [S(alpha) - B(alpha) e^(-tau/tau_c)]
    PSD del muestreo (paso dt)  = dt [s2 + v S + v B 2(r cos th - r²)/(1 - 2 r cos th + r²)],
                                  r = e^(-dt/tau_c), th = 2 pi f dt.
La última es la suma cerrada de los alias de la Lorentziana filtrada por la exposición,
sum_n L(f + n/dt) sinc²(pi (f + n/dt) W)  (Berg-Sørensen & Flyvbjerg, RSI 75, 594, 2004, con
el filtro de exposición en lugar del filtro del fotodiodo); sin exposición (S = B = 1) es la
Lorentziana con aliasing de BS&F: dt [s2 + v (1 - r²)/(1 - 2 r cos th + r²)].

Los cinco métodos:
  equiparticion   kappa = kB T S(alpha) / (Var - s2), alpha = W kappa / gamma  (autoconsistente)
  boltzmann       U(x) = -kB T ln P(x), parábola -> Var_B = kB T / U''; misma corrección
  autocorrelacion C(tau) = A e^(-tau/tau_c) con tau >= dt (sin tau = 0) -> kappa = gamma / tau_c
  espectro        máxima verosimilitud (Whittle) del periodograma de cada tramo sin huecos con
                  la PSD de arriba -> tau_c -> kappa = gamma / tau_c
  msd             ajuste del MSD atrapado de arriba -> kappa = kB T / v  (y tau_c)

"Sin corrección" (para ver qué hacen las correcciones): equipartición y Boltzmann sin restar
ruido ni corregir blur; MSD con (2 kB T/kappa)(1 - e^(-tau/tau_c)); espectro con la Lorentziana
con aliasing sin filtro de exposición ni ruido; la autocorrelación no tiene corrección.

Incertidumbres: bootstrap por bloques (los cuadros se parten en n_bloques bloques contiguos;
cada réplica sortea bloques con reposición y recalcula TODO desde las sumas de cada bloque:
momentos, histograma y sumas por retardo, así no se inventan pares en las uniones). Para el
espectro: matriz de Fisher de la verosimilitud de Whittle (analítica). La incertidumbre de gamma
NO está incluida (es sistemática y común a autocorrelación y espectro).
"""
import numpy as np

from .fisica import KB, factor_blur_covarianza as B_blur, factor_blur_varianza as S_blur
from .libre import pares_por_retardo, periodo_nominal, tramos_desde_t

METODOS = ('equiparticion', 'boltzmann', 'autocorrelacion', 'espectro', 'msd')


# ---------------------------------------------------------------- ruido y deriva
def piso_de_ruido(x_pegada, y_pegada, metodo='diferencias'):
    """Varianza del ruido de localización por eje, de una (o varias) beads pegadas al vidrio.

    x_pegada, y_pegada: [n] o [n, k] m (mismo ROI, exposición e iluminación que la medición).
    metodo: 'diferencias' -> s2 = <(x_{i+1} - x_i)²> / 2  (ruido blanco; no lo infla la deriva lenta)
            'varianza'    -> s2 = Var(x)  (literal; incluye la deriva que haya quedado)
    Con varias beads se promedian las varianzas.
    Devuelve dict: var_x, var_y [m²] (las del método elegido), sigma_x, sigma_y [m],
    var_x_dif, var_y_dif, var_x_directa, var_y_directa, n, metodo.
    Si var_directa >> var_dif hay deriva o ruido no blanco (vibración): mirarlo.
    """
    out = {'metodo': metodo}
    for eje, X in (('x', x_pegada), ('y', y_pegada)):
        X = np.asarray(X, float)
        X = X[:, None] if X.ndim == 1 else X
        dif = np.nanmean(np.nanmean(np.diff(X, axis=0) ** 2, axis=0)) / 2
        directa = np.nanmean(np.nanvar(X, axis=0))
        out[f'var_{eje}_dif'] = float(dif)
        out[f'var_{eje}_directa'] = float(directa)
        out[f'var_{eje}'] = float(dif if metodo == 'diferencias' else directa)
        out[f'sigma_{eje}'] = float(np.sqrt(out[f'var_{eje}']))
        out['n'] = int(np.isfinite(X).sum())
    return out


def segmentos_continuos(x, tramos, n_min=2):
    """Sub-tramos (i0, i1) inclusive sin huecos de tiempo NI NaN, de largo >= n_min."""
    fin = np.isfinite(np.asarray(x, float))
    seg = []
    for i0, i1 in np.asarray(tramos, int):
        ok = fin[i0:i1 + 1]
        borde = np.diff(np.concatenate([[0], ok.astype(int), [0]]))
        for a, b in zip(np.flatnonzero(borde == 1), np.flatnonzero(borde == -1)):
            if b - a >= n_min:
                seg.append((i0 + a, i0 + b - 1))
    return seg


def pasaaltos(t, x, tramos, f_corte):
    """Pasaaltos opcional contra la deriva: x - (x suavizada con una gaussiana EN TIEMPO REAL).

    La línea de base es el promedio pesado (Nadaraya-Watson) con núcleo gaussiano de ancho
    s_t = 1 / (2 pi f_corte):  base(t_i) = sum_j K(t_i - t_j) x_j / sum_j K(t_i - t_j),
    sobre los cuadros que existen (los perdidos y los NaN no pesan; no se interpola nada, y en
    los bordes la normalización evita transitorios). Se calcula en la grilla de período dt0.
    Respuesta en potencia G(f) = (1 - exp(-(f/f_corte)²/2))². f_corte [Hz] tiene que estar MUY
    por debajo de f_c: a una Lorentziana le saca una fracción de la varianza ~ 1.03 f_corte / f_c
    (integral de 1 - G contra la Lorentziana, f_corte << f_c) que NO se corrige: con
    f_corte = f_c / 100 es 1 %. El espectro se ajusta desde 4 f_corte (ahí G = 0.9993).
    t: [n] s; x: [n] m; tramos: [m, 2] (solo para dt0). Devuelve x filtrada [n] m (NaN donde era NaN).
    """
    from scipy.signal import fftconvolve
    t = np.asarray(t, float)
    x = np.asarray(x, float)
    dt = periodo_nominal(t, tramos)
    k = np.rint((t - t[0]) / dt).astype(np.int64)
    fin = np.isfinite(x)
    w = np.zeros(k[-1] + 1)
    xw = np.zeros(k[-1] + 1)
    np.add.at(w, k[fin], 1.0)
    np.add.at(xw, k[fin], x[fin])
    st = 1.0 / (2 * np.pi * f_corte) / dt                       # ancho en cuadros
    u = np.arange(-int(np.ceil(5 * st)), int(np.ceil(5 * st)) + 1)
    nucleo = np.exp(-0.5 * (u / st) ** 2)
    with np.errstate(invalid='ignore', divide='ignore'):
        base = fftconvolve(xw, nucleo, 'same') / fftconvolve(w, nucleo, 'same')
    return np.where(fin, x - base[k], np.nan)


# ---------------------------------------------------------------- estadísticos por bloque
def estadisticos_por_bloque(t, x, dt0, max_lag, n_bloques, bordes):
    """Sumas suficientes de cada bloque contiguo de cuadros, para el bootstrap por bloques.

    t: [n] s; x: [n] m (centrada); bordes: bordes del histograma [m].
    Devuelve dict de arrays con primer eje = bloque: n, s1, s2 [K]; hist [K, nbins];
    N, s_tau, s_d2, s_ab, s_a, s_b [K, max_lag + 1] (ver libre.pares_por_retardo).
    """
    x = np.asarray(x, float)
    n = len(x)
    bl = (np.arange(n) * n_bloques) // n
    fin = np.isfinite(x)
    b, xf = bl[fin], x[fin]
    nb = len(bordes) - 1
    k = np.clip(np.searchsorted(bordes, xf, side='right') - 1, -1, nb)
    dentro = (k >= 0) & (k < nb)
    st = {'n': np.bincount(b, minlength=n_bloques).astype(float),
          's1': np.bincount(b, xf, minlength=n_bloques),
          's2': np.bincount(b, xf * xf, minlength=n_bloques),
          'hist': np.bincount(b[dentro] * nb + k[dentro], minlength=n_bloques * nb)
                    .reshape(n_bloques, nb).astype(float)}
    p = pares_por_retardo(t, x, max_lag, dt0, bl, n_bloques)
    st.update({c: v[:, :, 0] for c, v in p.items()})
    return st


def _agregar(st, w=None):
    """Suma (ponderada por w [K], el número de veces que la réplica tomó cada bloque)."""
    if w is None:
        return {c: v.sum(0) for c, v in st.items()}
    return {c: np.tensordot(w, v, axes=(0, 0)) for c, v in st.items()}


def _momentos(a):
    """De las sumas agregadas a: media, varianza, tau [L], C(tau) [L], MSD(tau) [L] (m = 1..L)."""
    mu = a['s1'] / a['n']
    var = a['s2'] / a['n'] - mu * mu
    N = a['N'][1:]
    with np.errstate(invalid='ignore', divide='ignore'):
        tau = a['s_tau'][1:] / N
        C = (a['s_ab'][1:] - mu * (a['s_a'][1:] + a['s_b'][1:])) / N + mu * mu
        M = a['s_d2'][1:] / N
    return mu, var, tau, C, M


# ---------------------------------------------------------------- estimadores
def kappa_varianza(var, T, var_ruido=0.0, W=0.0, gamma=None, tau_c=None, corregir_blur=True):
    """kappa [N/m] de una varianza medida (equipartición, o la de Boltzmann), con correcciones.

    var: varianza medida [m²]; var_ruido: varianza del ruido a restar [m²] (0 = no se resta).
    Blur (Wong & Halvorsen): var - var_ruido = (kB T / kappa) S(alpha), alpha = W / tau_c.
      - con gamma: tau_c = gamma / kappa, se resuelve kappa autoconsistente (raíz única: el lado
        kappa (var - s2) crece y kB T S(W kappa / gamma) decrece).
      - sin gamma pero con tau_c (p. ej. de la autocorrelación): kappa = kB T S(W/tau_c) / (var - s2).
      - corregir_blur=False o W = 0: kappa = kB T / (var - var_ruido).
    Devuelve (kappa, alpha); NaN si var - var_ruido <= 0.
    """
    from scipy.optimize import brentq
    vp = var - (var_ruido or 0.0)
    if not np.isfinite(vp) or vp <= 0:
        return np.nan, np.nan
    k0 = KB * T / vp
    if not corregir_blur or not W:
        return k0, np.nan
    if gamma is not None:
        g = lambda k: k * vp - KB * T * float(S_blur(W * k / gamma))
        k = brentq(g, k0 * 1e-9, k0, xtol=k0 * 1e-12, rtol=1e-12)
        return k, W * k / gamma
    if tau_c is not None and np.isfinite(tau_c) and tau_c > 0:
        a = W / tau_c
        return k0 * float(S_blur(a)), a
    return np.nan, np.nan


def _poisson_polinomio(u, n, grado, iter_max=60):
    """Ajuste de máxima verosimilitud de Poisson de ln(mu_k) = polinomio(u_k) (IRLS).

    u: [nb] abscisas escaladas (adimensionales); n: [nb] cuentas (los ceros también entran).
    Devuelve coeficientes [grado+1] de mayor a menor potencia (como np.polyval) o None.
    """
    X = np.vander(u, grado + 1)
    beta = np.polyfit(u, np.log(n + 0.5), grado)
    for _ in range(iter_max):
        eta = X @ beta
        mu = np.exp(np.clip(eta, -700, 700))
        z = eta + (n - mu) / mu
        sw = np.sqrt(mu)
        nuevo, *_ = np.linalg.lstsq(X * sw[:, None], z * sw, rcond=None)
        if not np.all(np.isfinite(nuevo)):
            return None
        if np.max(np.abs(nuevo - beta)) < 1e-10:
            return nuevo
        beta = nuevo
    return beta


def ajuste_boltzmann(bordes, cuentas, T, sheppard=True):
    """Potencial U(x) = -kB T ln P(x) del histograma y ajuste parabólico.

    bordes [nb+1] m (uniformes), cuentas [nb]. U = c0 + c1 x + c2 x²  ->  kappa_B = 2 c2 [N/m].
    El ajuste es de máxima verosimilitud de Poisson sobre TODAS las cuentas (ln mu_k =
    -U(x_k)/kB T + cte), no cuadrados mínimos sobre ln(cuentas): el ln de pocas cuentas está
    sesgado y los pesos 1/cuentas correlacionan con los datos (con cuadrados mínimos medimos
    un sesgo de -0.5 a -1.3 % en kappa).
    sheppard: el histograma ensancha la varianza en w²/12 (corrección de Sheppard); se devuelve
      var_B = kB T / kappa_raw - w²/12 (con w = sigma/5 es 0.33 %).
    Armonicidad: se ajusta también U4 = a0 + a1 x + a2 x² + a3 x³ + a4 x⁴ (Poisson) y se mide
    cuánto pesan los términos no cuadráticos a ±2 sigma de la parábola (donde U - Umin = 2 kB T):
    anarmonia_2s = max_± |a3 x³ + a4 x⁴| / (2 kB T)  (x medido desde la media); además el
    chi2_red de Pearson de la parábola sobre los bins con >= 5 cuentas esperadas.
    Devuelve dict: kappa_raw [N/m], var_B [m²], x0 [m], chi2_red, anarmonia_2s,
    centros [nb] m, U [nb] J (NaN sin cuentas), err_U = kB T/sqrt(n), U_ajuste [nb] J, usado [nb].
    """
    centros = 0.5 * (bordes[1:] + bordes[:-1])
    w = np.diff(bordes)
    N = cuentas.sum()
    with np.errstate(divide='ignore', invalid='ignore'):
        U = np.where(cuentas > 0, -KB * T * np.log(cuentas / (N * w)), np.nan)
        eU = KB * T / np.sqrt(cuentas)
    out = {'kappa_raw': np.nan, 'var_B': np.nan, 'x0': np.nan, 'chi2_red': np.nan,
           'anarmonia_2s': np.nan, 'centros': centros, 'U': U, 'err_U': eU,
           'U_ajuste': np.full_like(centros, np.nan), 'usado': cuentas > 0}
    if (cuentas > 0).sum() < 4 or N < 20:
        return out
    m0 = np.sum(cuentas * centros) / N
    esc = np.sqrt(np.sum(cuentas * (centros - m0) ** 2) / N)
    u = (centros - m0) / esc
    c = _poisson_polinomio(u, cuentas, 2)
    if c is None or c[0] >= 0:
        return out
    kappa = -2 * KB * T * c[0] / esc ** 2
    x0 = m0 - c[1] / (2 * c[0]) * esc
    mu = np.exp(np.polyval(c, u))
    ok = mu >= 5
    chi2 = np.sum((cuentas[ok] - mu[ok]) ** 2 / mu[ok]) / max(ok.sum() - 3, 1)
    var_B = KB * T / kappa - (np.mean(w) ** 2 / 12 if sheppard else 0.0)
    out.update(kappa_raw=kappa, var_B=var_B, x0=x0, chi2_red=float(chi2),
               U_ajuste=-KB * T * np.log(mu / (N * w)))
    c4 = _poisson_polinomio(u, cuentas, 4)
    if c4 is not None:
        s = np.sqrt(KB * T / kappa)
        us = (np.array([x0 - 2 * s, x0 + 2 * s]) - m0) / esc
        no_cuad = c4[0] * us ** 4 + c4[1] * us ** 3            # términos x⁴ y x³ [kB T]
        out['anarmonia_2s'] = float(np.max(np.abs(no_cuad)) / 2)
    return out


def ajuste_autocorrelacion(tau, C, dt0, tau_c0, n_max=None):
    """Ajuste C(tau) = A e^(-tau/tau_c) sobre los retardos tau >= dt0 (tau = 0 EXCLUIDO).

    tau, C: [L] desde el primer retardo. Usa los retardos hasta 3 tau_c (mínimo 3 puntos,
    máximo n_max); itera dos veces porque el rango depende de tau_c.
    Para tau >= W el blur solo cambia la amplitud (A = v B(alpha)), no tau_c.
    Devuelve (tau_c [s], A [m²], n_puntos). NaN si el ajuste no converge.
    """
    from scipy.optimize import curve_fit
    ok = np.isfinite(tau) & np.isfinite(C)
    tau, C = tau[ok], C[ok]
    n_max = len(tau) if n_max is None else min(n_max, len(tau))
    tc = tau_c0 if np.isfinite(tau_c0) and tau_c0 > 0 else dt0
    A = np.nan
    n = 0
    for _ in range(2):
        n = int(np.clip(np.ceil(3 * tc / dt0), 3, n_max))
        if n < 2:
            return np.nan, np.nan, n
        esc = abs(C[0]) if C[0] != 0 else 1.0
        try:
            p, _ = curve_fit(lambda s, a, lt: a * np.exp(-s / np.exp(lt)), tau[:n], C[:n] / esc,
                             p0=[C[0] / esc * np.exp(tau[0] / tc), np.log(tc)], maxfev=2000)
        except (RuntimeError, ValueError):
            return np.nan, np.nan, n
        tc, A = float(np.exp(p[1])), float(p[0] * esc)
        if not np.isfinite(tc) or tc > 1e3 * dt0 or tc < 1e-3 * dt0:
            return np.nan, np.nan, n
    return tc, A, n


def modelo_msd_atrapado(tau, v, tau_c, var_ruido, W, corregido=True):
    """MSD de un eje [m²]: corregido 2 s2 + 2 v [S(a) - B(a) e^(-tau/tau_c)], a = W/tau_c;
    sin corrección (2 kB T/kappa)(1 - e^(-tau/tau_c)) = 2 v (1 - e^(-tau/tau_c)).  Vale para tau >= W."""
    e = np.exp(-np.asarray(tau, float) / tau_c)
    if not corregido:
        return 2 * v * (1 - e)
    a = W / tau_c
    return 2 * var_ruido + 2 * v * (S_blur(a) - B_blur(a) * e)


def ajuste_msd_atrapado(tau, M, W, dt0, v0, tau_c0, var_ruido=None, corregido=True, n_max=None):
    """Ajuste por cuadrados mínimos del MSD atrapado de un eje sobre los retardos hasta 5 tau_c.

    tau, M: [L] (desde el primer retardo). var_ruido: s2 fijo [m²], o None = s2 libre (solo
    corregido). Parámetros ln v, ln tau_c (+ s2). Devuelve dict v, tau_c, var_ruido, n_puntos.
    """
    from scipy.optimize import least_squares
    ok = np.isfinite(tau) & np.isfinite(M)
    tau, M = tau[ok], M[ok]
    n_max = len(tau) if n_max is None else min(n_max, len(tau))
    nan = {'v': np.nan, 'tau_c': np.nan, 'var_ruido': np.nan, 'n_puntos': 0}
    if not (np.isfinite(v0) and v0 > 0):
        return nan
    tc0 = tau_c0 if np.isfinite(tau_c0) and tau_c0 > 0 else dt0
    n = int(np.clip(np.ceil(5 * tc0 / dt0), 4, n_max))
    if n < 3:
        return nan
    libre = corregido and var_ruido is None
    s_fijo = 0.0 if not corregido else (var_ruido or 0.0)

    def res(p):
        s2 = p[2] * v0 if libre else s_fijo
        return (modelo_msd_atrapado(tau[:n], np.exp(p[0]), np.exp(p[1]), s2, W, corregido) - M[:n]) / v0

    p0 = [np.log(v0), np.log(tc0)] + ([0.01] if libre else [])
    try:
        r = least_squares(res, p0, method='lm' if n > len(p0) else 'trf')
    except (ValueError, RuntimeError):
        return nan
    v, tc = float(np.exp(r.x[0])), float(np.exp(r.x[1]))
    if tc > 1e3 * dt0 or tc < 1e-3 * dt0:
        return nan
    return {'v': v, 'tau_c': tc, 'var_ruido': float(r.x[2] * v0) if libre else s_fijo, 'n_puntos': n}


# ---------------------------------------------------------------- espectro
def periodogramas(t, x, segmentos, n_min=64):
    """Periodograma de cada segmento continuo (sin huecos ni NaN), con su propio dt.

    P_k = (dt/N) |sum_j x_j e^(-2 pi i j k/N)|²  [m²/Hz], f_k = k / (N dt), k = 1 .. ceil(N/2) - 1
    (sin f = 0 ni Nyquist). x se centra en cada segmento. Devuelve (f, P, dt) concatenados.
    """
    fs, Ps, dts = [], [], []
    for i0, i1 in segmentos:
        n = i1 - i0 + 1
        if n < n_min:
            continue
        s = x[i0:i1 + 1] - np.mean(x[i0:i1 + 1])
        dt = (t[i1] - t[i0]) / (n - 1)
        F = np.fft.rfft(s)
        k = np.arange(1, (n - 1) // 2 + 1)
        fs.append(k / (n * dt))
        Ps.append(dt / n * np.abs(F[k]) ** 2)
        dts.append(np.full(len(k), dt))
    if not fs:
        return np.zeros(0), np.zeros(0), np.zeros(0)
    return np.concatenate(fs), np.concatenate(Ps), np.concatenate(dts)


def modelo_psd(f, dt, v, tau_c, var_ruido=0.0, W=0.0, filtro=True):
    """PSD esperada del periodograma de las posiciones muestreadas [m²/Hz] (expresión cerrada).

    filtro=True: Lorentziana con aliasing y filtro de exposición (ver docstring del módulo),
      dt [s2 + v S(a) + v B(a) 2 (r cos th - r²) / (1 - 2 r cos th + r²)],  a = W/tau_c.
    filtro=False: BS&F con aliasing, muestreo instantáneo: dt [s2 + v (1 - r²)/(1 - 2 r cos th + r²)].
    r = e^(-dt/tau_c), th = 2 pi f dt. f, dt: arrays [Hz], [s]; v [m²]; tau_c [s]; s2 [m²].
    """
    r = np.exp(-dt / tau_c)
    c = np.cos(2 * np.pi * f * dt)
    den = 1 - 2 * r * c + r * r
    if filtro:
        a = W / tau_c
        return dt * (var_ruido + v * S_blur(a) + v * B_blur(a) * 2 * (r * c - r * r) / den)
    return dt * (var_ruido + v * (1 - r * r) / den)


def psd_suma_alias(f, dt, v, tau_c, W, n_alias=2000):
    """Misma PSD (sin ruido) por la suma explícita de alias, para verificar la forma cerrada:
    sum_n L(f + n/dt) sinc²(pi (f + n/dt) W),  L(f) = 2 v tau_c / (1 + (2 pi f tau_c)²) la
    densidad BILATERAL del OU [m²/Hz] (integra a v). Lento: solo para tests.
    Todas las PSD de este módulo son bilaterales: sum_k P_k / (N dt) sobre -N/2..N/2 = varianza.
    """
    n = np.arange(-n_alias, n_alias + 1)[:, None]
    fn = np.asarray(f, float)[None, :] + n / dt
    L = 2 * v * tau_c / (1 + (2 * np.pi * fn * tau_c) ** 2)
    return np.sum(L * np.sinc(fn * W) ** 2, axis=0)


def ajuste_espectro(f, P, dt, W, v0, tau_c0, var_ruido=None, filtro=True, f_min=None, f_max=None):
    """Máxima verosimilitud (Whittle) del periodograma con `modelo_psd`.

    -ln L = sum_k [ln P_mod(f_k) + P_k / P_mod(f_k)]  (P_k ~ P_mod chi²_2 / 2, independientes).
    Es el estimador de máxima verosimilitud de BS&F para bloques de 1 frecuencia; vale para
    tramos de largos distintos porque cada punto lleva su dt.
    Parámetros ln v, ln tau_c y, si var_ruido es None y filtro=True, ln s2 (libre).
    Sin corrección (filtro=False) el ruido no se modela (s2 = 0).
    Incertidumbre: Fisher I_ab = sum_k d_a ln P_mod d_b ln P_mod -> cov = I^-1.
    Devuelve dict v, tau_c, var_ruido, err_tau_c, err_v, n_frec.
    """
    from scipy.optimize import minimize
    sel = np.isfinite(P) & (P > 0)
    if f_min is not None:
        sel &= f >= f_min
    if f_max is not None:
        sel &= f <= f_max
    f, P, dt = f[sel], P[sel], dt[sel]
    nan = {'v': np.nan, 'tau_c': np.nan, 'var_ruido': np.nan, 'err_tau_c': np.nan,
           'err_v': np.nan, 'n_frec': int(len(f))}
    if len(f) < 10 or not (np.isfinite(v0) and v0 > 0):
        return nan
    libre = filtro and var_ruido is None
    s_fijo = (var_ruido or 0.0) if filtro else 0.0
    tc0 = tau_c0 if np.isfinite(tau_c0) and tau_c0 > 0 else np.median(dt)

    def modelo(p):
        s2 = np.exp(p[2]) if libre else s_fijo
        return modelo_psd(f, dt, np.exp(p[0]), np.exp(p[1]), s2, W, filtro)

    def nll(p):
        m = modelo(p)
        if not np.all(np.isfinite(m)) or np.any(m <= 0):
            return np.inf
        return np.sum(np.log(m) + P / m) / len(f)

    p0 = [np.log(v0), np.log(tc0)] + ([np.log(0.01 * v0)] if libre else [])
    r = minimize(nll, p0, method='Nelder-Mead',
                 options={'xatol': 1e-7, 'fatol': 1e-12, 'maxiter': 4000, 'maxfev': 8000})
    p = r.x
    # Fisher numérica en los parámetros logarítmicos
    g = []
    for i in range(len(p)):
        h = np.zeros(len(p))
        h[i] = 1e-5
        g.append((np.log(modelo(p + h)) - np.log(modelo(p - h))) / 2e-5)
    g = np.array(g)
    try:
        cov = np.linalg.inv(g @ g.T)
    except np.linalg.LinAlgError:
        cov = np.full((len(p), len(p)), np.nan)
    v, tc = float(np.exp(p[0])), float(np.exp(p[1]))
    return {'v': v, 'tau_c': tc, 'var_ruido': float(np.exp(p[2])) if libre else s_fijo,
            'err_tau_c': tc * float(np.sqrt(cov[1, 1])), 'err_v': v * float(np.sqrt(cov[0, 0])),
            'n_frec': int(len(f))}


def psd_en_bloques(f, P, Pm, n_bins=40):
    """Promedio del periodograma (y del modelo) en bins logarítmicos de frecuencia, para graficar.
    Devuelve (f_media, P_media, Pm_media, n_por_bin)."""
    ok = (f > 0) & np.isfinite(P)
    f, P, Pm = f[ok], P[ok], Pm[ok]
    if len(f) == 0:
        return (np.zeros(0),) * 4
    b = np.geomspace(f.min(), f.max() * 1.0001, n_bins + 1)
    k = np.clip(np.searchsorted(b, f, side='right') - 1, 0, n_bins - 1)
    n = np.bincount(k, minlength=n_bins)
    m = n > 0
    prom = lambda w: np.bincount(k, w, minlength=n_bins)[m] / n[m]
    return prom(f), prom(P), prom(Pm), n[m]


# ---------------------------------------------------------------- un eje
def _estimar(a, ctx, corregido):
    """Los cuatro métodos basados en momentos/sumas (todos menos el espectro), para un juego
    de sumas agregadas `a` (los datos completos o una réplica bootstrap). Un ajuste que no
    converge da NaN (sin warnings)."""
    import warnings
    with warnings.catch_warnings(), np.errstate(all='ignore'):
        warnings.simplefilter('ignore')
        return _estimar_sin_avisos(a, ctx, corregido)


def _estimar_sin_avisos(a, ctx, corregido):
    T, W, gamma, s2, dt0 = ctx['T'], ctx['W'], ctx['gamma'], ctx['var_ruido'], ctx['dt0']
    s2_resta = (s2 or 0.0) if corregido else 0.0
    _, var, tau, C, M = _momentos(a)
    out = {}
    tc_acf, A, n_acf = ajuste_autocorrelacion(tau, C, dt0, ctx['tau_c0'], ctx['max_lag'])
    out['autocorrelacion'] = {'kappa': gamma / tc_acf if gamma else np.nan, 'tau_c': tc_acf,
                              'amplitud': A, 'n_puntos': n_acf}
    k_eq, al = kappa_varianza(var, T, s2_resta, W, gamma, tc_acf, corregir_blur=corregido)
    out['equiparticion'] = {'kappa': k_eq, 'tau_c': np.nan, 'alpha': al, 'var': var}
    bz = ajuste_boltzmann(ctx['bordes'], a['hist'], T)
    k_bz, _ = kappa_varianza(bz['var_B'], T, s2_resta, W, gamma, tc_acf, corregir_blur=corregido)
    out['boltzmann'] = {'kappa': k_bz, 'tau_c': np.nan, 'ajuste': bz}
    v0 = max(var - s2_resta, 1e-3 * var)
    ms = ajuste_msd_atrapado(tau, M, W, dt0, v0, tc_acf if np.isfinite(tc_acf) else ctx['tau_c0'],
                             s2 if corregido else 0.0, corregido, ctx['max_lag'])
    out['msd'] = {'kappa': KB * T / ms['v'] if ms['v'] > 0 else np.nan, 'tau_c': ms['tau_c'],
                  'ajuste': ms}
    out['_momentos'] = (var, tau, C, M)
    return out


def analizar_eje(t, x, W, T, gamma, tramos, var_ruido=None, corregido=True, max_lag=None,
                 n_bloques=20, n_boot=200, rng=None, f_min_hz=None, f_max_hz=None):
    """Los cinco métodos en UN eje. x: [n] m ya centrada (y filtrada si corresponde).

    Devuelve dict: metodo -> {kappa, err_kappa, tau_c, err_tau_c, ...} y 'graficos'.
    """
    rng = np.random.default_rng(rng)
    t = np.asarray(t, float)
    x = np.asarray(x, float)
    dt0 = periodo_nominal(t, tramos)
    fin = np.isfinite(x)
    sd = np.nanstd(x)
    # tau_c preliminar con la correlación a 1 cuadro (solo para elegir rangos)
    seg = segmentos_continuos(x, tramos, 2)
    num = sum(np.sum(x[i0:i1] * x[i0 + 1:i1 + 1]) for i0, i1 in seg)
    den = sum(np.sum(x[i0:i1] ** 2) for i0, i1 in seg)
    r1 = np.clip(num / den, 0.02, 0.995) if den > 0 else 0.5
    tau_c0 = -dt0 / np.log(r1)
    if max_lag is None:
        max_lag = int(np.clip(np.ceil(8 * tau_c0 / dt0), 10, 200))
    bordes = np.nanmean(x) + sd * np.arange(-5, 5.0001, 0.2)       # bins de sigma/5
    ctx = {'T': T, 'W': W, 'gamma': gamma, 'var_ruido': var_ruido, 'dt0': dt0,
           'tau_c0': tau_c0, 'max_lag': max_lag, 'bordes': bordes}
    st = estadisticos_por_bloque(t, x, dt0, max_lag, n_bloques, bordes)
    est = _estimar(_agregar(st), ctx, corregido)
    var, tau, C, M = est.pop('_momentos')
    # espectro
    f, P, dts = periodogramas(t, x, seg)
    tc_ini = est['autocorrelacion']['tau_c']
    tc_ini = tc_ini if np.isfinite(tc_ini) else tau_c0
    s2_esp = var_ruido if corregido else 0.0
    v0 = max(var - (s2_esp or 0.0), 1e-3 * var)
    with np.errstate(all='ignore'):
        es = ajuste_espectro(f, P, dts, W, v0, tc_ini, s2_esp, filtro=corregido,
                             f_min=f_min_hz, f_max=f_max_hz)
    est['espectro'] = {'kappa': gamma / es['tau_c'] if gamma else np.nan, 'tau_c': es['tau_c'],
                       'err_tau_c': es['err_tau_c'], 'ajuste': es}
    # bootstrap por bloques
    reps = {m: {'kappa': [], 'tau_c': []} for m in METODOS if m != 'espectro'}
    for _ in range(int(n_boot)):
        w = rng.multinomial(n_bloques, np.full(n_bloques, 1.0 / n_bloques)).astype(float)
        e = _estimar(_agregar(st, w), ctx, corregido)
        for m in reps:
            reps[m]['kappa'].append(e[m]['kappa'])
            reps[m]['tau_c'].append(e[m]['tau_c'])
    for m in METODOS:
        d = est[m]
        if m == 'espectro':
            d['err_kappa'] = d['kappa'] * d['err_tau_c'] / d['tau_c'] if gamma else np.nan
        else:
            sd_ = lambda v: float(np.nanstd(v, ddof=1)) if np.isfinite(v).sum() > 2 else np.nan
            d['err_kappa'] = sd_(np.array(reps[m]['kappa'])) if n_boot else np.nan
            d['err_tau_c'] = sd_(np.array(reps[m]['tau_c'])) if n_boot else np.nan
        d['f_c'] = 1 / (2 * np.pi * d['tau_c']) if np.isfinite(d['tau_c']) else np.nan
        d['err_f_c'] = d['f_c'] * d['err_tau_c'] / d['tau_c'] if np.isfinite(d['tau_c']) else np.nan
    # datos para graficar
    bz = est['boltzmann']['ajuste']
    N = st['hist'].sum()
    Pm = (modelo_psd(f, dts, es['v'], es['tau_c'], es['var_ruido'], W, corregido)
          if np.isfinite(es['tau_c']) else np.full_like(f, np.nan))
    fb, Pb, Pmb, nfb = psd_en_bloques(f, P, Pm)
    tc_a, A = est['autocorrelacion']['tau_c'], est['autocorrelacion']['amplitud']
    ms = est['msd']['ajuste']
    graf = {
        'x_hist': bz['centros'], 'P': st['hist'].sum(0) / (N * np.diff(bordes)) if N else bz['centros'] * np.nan,
        'U': bz['U'], 'err_U': bz['err_U'], 'U_ajuste': bz['U_ajuste'], 'U_usado': bz['usado'],
        'tau': tau, 'C': C, 'C0': var,
        'C_ajuste': A * np.exp(-tau / tc_a) if np.isfinite(tc_a) else tau * np.nan,
        'msd': M, 'msd_ajuste': (modelo_msd_atrapado(tau, ms['v'], ms['tau_c'], ms['var_ruido'], W, corregido)
                                 if np.isfinite(ms['tau_c']) else tau * np.nan),
        'f': fb, 'psd': Pb, 'psd_ajuste': Pmb, 'psd_n': nfb,
    }
    est['graficos'] = graf
    est['dt0'] = dt0
    est['max_lag'] = max_lag
    est['var'] = var
    est['armonicidad'] = {'chi2_red': bz['chi2_red'], 'anarmonia_2s': bz['anarmonia_2s']}
    return est


# ---------------------------------------------------------------- análisis completo
def _tabla(por_eje):
    cols = ('kappa', 'err_kappa', 'tau_c', 'err_tau_c', 'f_c', 'err_f_c')
    tab = {'metodo': [], 'eje': []}
    tab.update({c: [] for c in cols})
    for eje, est in por_eje.items():
        for m in METODOS:
            tab['metodo'].append(m)
            tab['eje'].append(eje)
            for c in cols:
                tab[c].append(float(est[m].get(c, np.nan)))
    return {c: np.array(v) for c, v in tab.items()}


def compatibilidad(tabla, eje):
    """Media pesada de kappa entre métodos de un eje, chi² reducido y máximo |z| por pares.
    Ojo: equipartición, Boltzmann y MSD usan la misma varianza (están correlacionados)."""
    sel = (tabla['eje'] == eje) & np.isfinite(tabla['kappa']) & np.isfinite(tabla['err_kappa']) \
        & (tabla['err_kappa'] > 0)
    k, e = tabla['kappa'][sel], tabla['err_kappa'][sel]
    if len(k) < 2:
        return {'kappa_medio': np.nan, 'err_kappa_medio': np.nan, 'chi2_red': np.nan, 'z_max': np.nan}
    w = 1 / e ** 2
    km = np.sum(w * k) / np.sum(w)
    z = np.abs(k[:, None] - k[None, :]) / np.sqrt(e[:, None] ** 2 + e[None, :] ** 2)
    return {'kappa_medio': float(km), 'err_kappa_medio': float(1 / np.sqrt(np.sum(w))),
            'chi2_red': float(np.sum(w * (k - km) ** 2) / (len(k) - 1)), 'z_max': float(np.max(z))}


def analizar(t, x, y, W, T, gamma, tramos=None, ruido=None, corregir=True, pasaaltos_hz=None,
             max_lag=None, n_bloques=20, n_boot=200, semilla=0, f_min_hz=None, f_max_hz=None,
             sin_correccion=True):
    """Análisis TRAMPA completo de una bead atrapada, por eje.

    t: [n] s reales (TIEMPO). x, y: [n] m de la bead atrapada, ya sin deriva (NaN = perdido).
    W: exposición [s]. T: temperatura [K].
    gamma: arrastre [N·s/m] = gamma_stokes(eta, R) * faxen_paralelo(R, h); None -> los métodos
      que dan tau_c no dan kappa, y el blur de equipartición/Boltzmann usa tau_c de la autocorrelación.
    tramos: [m, 2] sin huecos (TIEMPO); None = se deducen de t.
    ruido: dict de `piso_de_ruido` (usa var_x, var_y); None = no se resta en equipartición ni
      Boltzmann (aviso) y MSD y espectro ajustan el ruido como parámetro libre.
    corregir: aplicar ruido y blur (True) o no.
    pasaaltos_hz: corte del pasaaltos contra deriva (None = no; ver `pasaaltos`). Si se usa y no
      se da f_min_hz, el espectro se ajusta desde 4 * pasaaltos_hz.
    n_bloques, n_boot, semilla: bootstrap por bloques (n_boot = 0 -> errores NaN, más rápido).
    sin_correccion: calcular también la tabla sin ninguna corrección.

    Devuelve dict:
      tabla:  {metodo [10], eje [10], kappa N/m, err_kappa, tau_c s, err_tau_c, f_c Hz, err_f_c}
      tabla_sin_correccion: igual (o None)
      compatibilidad: {'x': {...}, 'y': {...}} (ver `compatibilidad`)
      alpha: {'x', 'y'} = W / tau_c (tau_c del espectro), S_alpha
      armonicidad: {'x': {chi2_red, anarmonia_2s}, 'y': ...}
      graficos: {'x': {x_hist, P, U, U_ajuste, tau, C, C_ajuste, msd, msd_ajuste, f, psd,
                       psd_ajuste}, 'y': ...}   (SI: m, 1/m, J, s, m², Hz, m²/Hz)
      ruido_usado, parametros, avisos
    """
    t = np.asarray(t, float)
    if tramos is None:
        tramos = tramos_desde_t(t)
    avisos = []
    rng = np.random.default_rng(semilla)
    if f_min_hz is None and pasaaltos_hz:
        f_min_hz = 4 * pasaaltos_hz
    por_eje, sin, graf, alpha, S_a, armon, ruido_usado = {}, {}, {}, {}, {}, {}, {}
    for eje, X in (('x', x), ('y', y)):
        X = np.asarray(X, float)
        X = X - np.nanmean(X)
        if pasaaltos_hz:
            X = pasaaltos(t, X, tramos, pasaaltos_hz)
        s2 = None if ruido is None else float(ruido[f'var_{eje}'])
        ruido_usado[eje] = s2
        semilla_eje = rng.integers(2 ** 32)
        est = analizar_eje(t, X, W, T, gamma, tramos, s2, corregir, max_lag, n_bloques, n_boot,
                           semilla_eje, f_min_hz, f_max_hz)
        por_eje[eje] = est
        graf[eje] = est['graficos']
        armon[eje] = est['armonicidad']
        tc = est['espectro']['tau_c']
        alpha[eje] = W / tc if np.isfinite(tc) else np.nan
        S_a[eje] = float(S_blur(alpha[eje])) if np.isfinite(alpha[eje]) else np.nan
        if sin_correccion:
            sin[eje] = analizar_eje(t, X, W, T, gamma, tramos, None, False, max_lag, n_bloques,
                                    n_boot, semilla_eje, f_min_hz, f_max_hz)
        fc = 1 / (2 * np.pi * tc) if np.isfinite(tc) else np.nan
        if np.isfinite(fc) and 1 / est['dt0'] < 2 * fc:
            avisos.append(f'{eje}: fps = {1 / est["dt0"]:.0f} < 2 f_c = {2 * fc:.1f} Hz: '
                          'autocorrelación y espectro poco confiables')
        if np.isfinite(alpha[eje]) and alpha[eje] > 3:
            avisos.append(f'{eje}: alpha = W/tau_c = {alpha[eje]:.2f} > 3: la corrección de blur domina')
        if s2 is not None and s2 > 0.3 * est['var']:
            avisos.append(f'{eje}: el ruido es {100 * s2 / est["var"]:.0f} % de la varianza medida')
        if pasaaltos_hz and np.isfinite(fc) and pasaaltos_hz > fc / 20:
            avisos.append(f'{eje}: pasaaltos {pasaaltos_hz} Hz no está muy por debajo de f_c = {fc:.1f} Hz')
        if np.isfinite(armon[eje]['anarmonia_2s']) and armon[eje]['anarmonia_2s'] > 0.1:
            avisos.append(f'{eje}: el potencial se aparta de la parábola '
                          f'({100 * armon[eje]["anarmonia_2s"]:.0f} % a 2 sigma)')
    if ruido is None and corregir:
        avisos.append('sin bead pegada: no se restó ruido en equipartición ni Boltzmann')
    if gamma is None:
        avisos.append('sin gamma: autocorrelación y espectro dan tau_c pero no kappa')
    tabla = _tabla(por_eje)
    return {
        'tabla': tabla,
        'tabla_sin_correccion': _tabla(sin) if sin_correccion else None,
        'compatibilidad': {e: compatibilidad(tabla, e) for e in ('x', 'y')},
        'alpha': alpha, 'S_alpha': S_a, 'armonicidad': armon, 'graficos': graf,
        'ruido_usado': ruido_usado,
        'parametros': {'W': W, 'T': T, 'gamma': gamma, 'corregir': corregir,
                       'pasaaltos_hz': pasaaltos_hz, 'n_bloques': n_bloques, 'n_boot': n_boot,
                       'semilla': semilla, 'f_min_hz': f_min_hz, 'f_max_hz': f_max_hz,
                       'max_lag': por_eje['x']['max_lag'], 'dt0': por_eje['x']['dt0']},
        'avisos': avisos,
    }
