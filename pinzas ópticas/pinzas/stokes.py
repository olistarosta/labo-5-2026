"""STOKES — arrastre de Stokes: corrimiento de equilibrio vs velocidad y velocidad de escape.

Física. La cámara y la trampa están fijas; la platina (y con ella el vidrio y el fluido, que se
mueve como un sólido entre los dos vidrios) se mueve a velocidad v. Las beads pegadas al vidrio
se mueven a v en la imagen. Sobre la bead atrapada, por eje:
    gamma_par (dx/dt - v) = -kappa (x - x_trampa) + ruido térmico
=>  x(t) = x_eq + OU,   x_eq = x_trampa + gamma_par v / kappa        (en régimen, t >> tau_c)
O sea: la bead se corre en el sentido de la velocidad del fluido, Δx = (gamma_par/kappa) v, y
fluctúa alrededor de x_eq igual que sin flujo. La pendiente de Δx vs v da gamma/kappa [s] =
tau_c, que se compara con TRAMPA (tau_c de autocorrelación/espectro o gamma/kappa).
Escape: si la fuerza de arrastre supera la máxima de la trampa, la bead se va:
    F_max = gamma_par v_max.

IMPORTANTE: las posiciones de la bead atrapada tienen que estar en el sistema de la CÁMARA
(trampa fija). Si DERIVA restó el movimiento de las beads pegadas, se restó justamente la
platina y el corrimiento desaparece. La velocidad sale de las beads pegadas en la misma imagen.

Convenciones: SI (m, s, N/m), NaN = no medido. Cada "segmento" es un tramo a velocidad
constante (sin las aceleraciones), con su t, la bead atrapada y las pegadas.
"""
import numpy as np


def velocidad_platina(t, x_pegadas, y_pegadas):
    """Velocidad de la platina [m/s] por ajuste lineal de la posición media de las beads pegadas.

    t: [n] s. x_pegadas, y_pegadas: [n] o [n, k] m (k beads pegadas; se promedian por cuadro).
    Devuelve dict vx, vy, err_vx, err_vy [m/s] (error de cuadrados mínimos con la dispersión
    de los residuos; si hay vibraciones correlacionadas lo subestima).
    """
    t = np.asarray(t, float)
    out = {}
    for eje, X in (('x', x_pegadas), ('y', y_pegadas)):
        X = np.asarray(X, float)
        X = X[:, None] if X.ndim == 1 else X
        with np.errstate(invalid='ignore'):
            m = np.nanmean(X, axis=1) if np.isfinite(X).any() else np.full(len(t), np.nan)
        ok = np.isfinite(m)
        if ok.sum() < 3:
            out[f'v{eje}'], out[f'err_v{eje}'] = np.nan, np.nan
            continue
        tt = t[ok] - t[ok].mean()
        b, a = np.polyfit(tt, m[ok], 1)
        res = m[ok] - (a + b * tt)
        out[f'v{eje}'] = float(b)
        out[f'err_v{eje}'] = float(np.sqrt(np.sum(res ** 2) / (ok.sum() - 2) / np.sum(tt ** 2)))
    return out


def corrimiento(x, n_bloques=10):
    """Posición media de la bead atrapada en un segmento y su error [m].

    x: [n] m. El error sale de las medias de n_bloques bloques contiguos (las muestras están
    correlacionadas en tau_c): err = std(medias) / sqrt(n_bloques).
    Devuelve (media, err).
    """
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < n_bloques * 2:
        return np.nan, np.nan
    medias = np.array([b.mean() for b in np.array_split(x, n_bloques)])
    return float(x.mean()), float(medias.std(ddof=1) / np.sqrt(n_bloques))


def ajustar_corrimiento(v, dx, err_dx):
    """Recta pesada dx = a + b v  ->  b = gamma_par / kappa [s].

    v: [k] m/s; dx: [k] m (posición media de cada segmento, con cualquier origen: la ordenada
    a absorbe la posición de la trampa); err_dx: [k] m. El error de v se ignora (es mucho menor).
    Devuelve dict pendiente [s], err_pendiente, ordenada [m], err_ordenada, chi2_red, n.
    """
    v, dx, e = (np.asarray(a, float) for a in (v, dx, err_dx))
    ok = np.isfinite(v) & np.isfinite(dx) & np.isfinite(e) & (e > 0)
    nan = {'pendiente': np.nan, 'err_pendiente': np.nan, 'ordenada': np.nan,
           'err_ordenada': np.nan, 'chi2_red': np.nan, 'n': int(ok.sum())}
    if ok.sum() < 2 or np.ptp(v[ok]) == 0:
        return nan
    w = 1 / e[ok] ** 2
    X = np.stack([np.ones(ok.sum()), v[ok]], 1)
    cov = np.linalg.inv(X.T @ (X * w[:, None]))
    a, b = cov @ (X.T @ (w * dx[ok]))
    chi2 = np.sum(w * (dx[ok] - a - b * v[ok]) ** 2) / max(ok.sum() - 2, 1)
    return {'pendiente': float(b), 'err_pendiente': float(np.sqrt(cov[1, 1])), 'ordenada': float(a),
            'err_ordenada': float(np.sqrt(cov[0, 0])), 'chi2_red': float(chi2), 'n': int(ok.sum())}


def sigue_atrapada(x, y, x_ref, y_ref, radio, frac_min=0.95):
    """¿La bead siguió en la trampa durante el segmento?

    x, y: [n] m; (x_ref, y_ref): posición de la trampa [m]; radio: distancia máxima [m] (p. ej.
    unas veces el radio de la bead). Atrapada si al menos frac_min de los cuadros tiene tracking
    válido y está a menos de `radio`. Un tracking perdido (NaN) cuenta como fuera.
    """
    d = np.hypot(np.asarray(x, float) - x_ref, np.asarray(y, float) - y_ref)
    return bool(np.mean(np.isfinite(d) & (d < radio)) >= frac_min)


def velocidad_escape(v, atrapada):
    """Velocidad de escape [m/s] a partir de segmentos a distintas |v|.

    v: [k] m/s (se usa |v|); atrapada: [k] bool. v_max queda entre la mayor |v| con la bead
    atrapada y la menor |v| en que escapó (mayor que la anterior). Devuelve dict v_max (punto
    medio), err_v_max (medio intervalo, uniforme -> /sqrt(3) NO aplicado: es una cota), v_abajo,
    v_arriba. NaN si nunca escapó o nunca estuvo atrapada.
    """
    v = np.abs(np.asarray(v, float))
    a = np.asarray(atrapada, bool)
    if not a.any() or a.all():
        return {'v_max': np.nan, 'err_v_max': np.nan,
                'v_abajo': float(v[a].max()) if a.any() else np.nan,
                'v_arriba': float(v[~a].min()) if (~a).any() else np.nan}
    abajo = v[a].max()
    sobre = v[~a & (v > abajo)]
    arriba = sobre.min() if len(sobre) else np.nan
    return {'v_max': float((abajo + arriba) / 2), 'err_v_max': float((arriba - abajo) / 2),
            'v_abajo': float(abajo), 'v_arriba': float(arriba)}


def analizar(segmentos, gamma=None, kappa=None, radio_escape=None, n_bloques=10):
    """Análisis STOKES completo.

    segmentos: lista de dicts, uno por tramo a velocidad constante, con
      t [n] s, x, y [n] m (bead atrapada, sistema de la cámara), x_pegadas, y_pegadas [n] o [n, k] m.
      Incluir uno con la platina quieta (v = 0) mejora la ordenada, pero no es obligatorio.
    gamma: gamma_par [N·s/m] (con Faxén) para F_max y para kappa_stokes = gamma / pendiente.
    kappa: {'x': kx, 'y': ky} [N/m] de TRAMPA, para comparar pendiente con gamma/kappa.
    radio_escape: [m] distancia a la trampa para decidir escape (None = no se evalúa).
      La trampa se toma como la posición media del segmento de menor |v|.

    Devuelve dict:
      v_x, v_y, err_v_x, err_v_y [k] m/s  (de las beads pegadas)
      dx, err_dx, dy, err_dy [k] m        (posición media de la bead atrapada)
      atrapada [k] bool (o None)
      ajuste: {'x': ajustar_corrimiento(...), 'y': ...}  pendiente = gamma/kappa [s]
              (solo con los segmentos en que siguió atrapada, si se evaluó el escape)
      kappa_stokes: {'x', 'y'} = gamma / pendiente [N/m] (NaN sin gamma)
      pendiente_esperada: {'x', 'y'} = gamma / kappa de TRAMPA [s] (NaN si falta)
      escape: dict de velocidad_escape (usa |v| = hypot(vx, vy)); F_max, err_F_max [N]
      avisos
    """
    avisos = []
    k = len(segmentos)
    cols = {c: np.full(k, np.nan) for c in ('v_x', 'v_y', 'err_v_x', 'err_v_y', 'dx', 'err_dx', 'dy', 'err_dy')}
    for i, s in enumerate(segmentos):
        vel = velocidad_platina(s['t'], s['x_pegadas'], s['y_pegadas'])
        cols['v_x'][i], cols['v_y'][i] = vel['vx'], vel['vy']
        cols['err_v_x'][i], cols['err_v_y'][i] = vel['err_vx'], vel['err_vy']
        cols['dx'][i], cols['err_dx'][i] = corrimiento(s['x'], n_bloques)
        cols['dy'][i], cols['err_dy'][i] = corrimiento(s['y'], n_bloques)
    atrapada, esc = None, {'v_max': np.nan, 'err_v_max': np.nan}
    usar = np.ones(k, bool)
    if radio_escape is not None:
        vmod = np.hypot(cols['v_x'], cols['v_y'])
        ref = segmentos[int(np.nanargmin(vmod))]
        xr, yr = np.nanmean(ref['x']), np.nanmean(ref['y'])
        atrapada = np.array([sigue_atrapada(s['x'], s['y'], xr, yr, radio_escape) for s in segmentos])
        esc = velocidad_escape(vmod, atrapada)
        usar = atrapada
        if atrapada.all():
            avisos.append('la bead nunca escapó: F_max > gamma * v máxima probada')
    ajuste, k_st, esperada = {}, {}, {}
    for eje in ('x', 'y'):
        aj = ajustar_corrimiento(cols[f'v_{eje}'][usar], cols[f'd{eje}'][usar], cols[f'err_d{eje}'][usar])
        ajuste[eje] = aj
        k_st[eje] = gamma / aj['pendiente'] if gamma and aj['pendiente'] > 0 else np.nan
        esperada[eje] = gamma / kappa[eje] if (gamma and kappa and eje in kappa) else np.nan
        if aj['n'] < 3:
            avisos.append(f'{eje}: menos de 3 velocidades distintas: pendiente sin control')
    esc['F_max'] = gamma * esc['v_max'] if gamma else np.nan
    esc['err_F_max'] = gamma * esc['err_v_max'] if gamma else np.nan
    out = dict(cols)
    out.update(atrapada=atrapada, ajuste=ajuste, kappa_stokes=k_st, pendiente_esperada=esperada,
               escape=esc, avisos=avisos)
    return out
