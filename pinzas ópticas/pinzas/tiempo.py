"""Reconstrucción del tiempo de cada cuadro y detección de cuadros perdidos.

Reglas (LECCIONES.md, "El tiempo"):
  - El tiempo sale de la cámara; el reloj del host (`t_host_ns`) es control o respaldo.
    Nunca de los fps nominales ni del reloj de TD. Única excepción: una grabación IMPORTADA
    de otro programa que no guardó ningún reloj (pinzas.importar); ahí metodo = 'nominal' y
    el informe lo dice en 'aviso'.
  - Los cuadros perdidos salen SOLO del contador de hardware (salto > 1), nunca de umbrales
    sobre dt.
  - Un contador que se repite o retrocede se reporta tal cual; no se repara ni se inventa.
  - El período se estima con todo el registro ((t_ult - t_pri) / (c_ult - c_pri)), nunca con
    la mediana de dt (en fluidos la mediana eligió un valor cuantizado).

Entrada: el dict `tiempos` de io.leer_grabacion (columnas de tiempos.csv): contador_camara,
t_camara (ticks crudos, NaN si no hay), t_host_ns.
"""
import numpy as np


# ---------------------------------------------------------------- contador
def analizar_contador(contador):
    """Saltos del contador de hardware.

    contador: [n] enteros, uno por cuadro guardado.
    Devuelve dict:
      saltos      [n-1] c[i+1] - c[i]
      huecos      índices i con salto > 1 (faltan salto-1 cuadros entre el guardado i y el i+1)
      perdidos    [len(huecos)] cuadros que faltan en cada hueco
      anomalias   índices i con salto <= 0 (contador repetido o que retrocede)
      patron_0_2  cuántas veces aparece un salto 0 seguido de un salto 2: la firma de un
                  contador leído antes de actualizarse (pasó en fluidos). Se reporta, no se repara.
    """
    c = np.asarray(contador, dtype=np.int64)
    s = np.diff(c)
    huecos = np.flatnonzero(s > 1)
    anom = np.flatnonzero(s <= 0)
    p02 = int(np.sum((s[:-1] == 0) & (s[1:] == 2))) if len(s) > 1 else 0
    return dict(saltos=s, huecos=huecos, perdidos=s[huecos] - 1, anomalias=anom, patron_0_2=p02)


def tramos_sin_huecos(n, cortes):
    """Tramos [m, 2] (i0, i1 inclusive) de 0..n-1 cortando después de cada índice de `cortes`."""
    if n == 0:
        return np.zeros((0, 2), dtype=np.int64)
    cortes = np.unique(np.asarray(cortes, dtype=np.int64))
    i0 = np.concatenate([[0], cortes + 1])
    i1 = np.concatenate([cortes, [n - 1]])
    return np.stack([i0, i1], axis=1).astype(np.int64)


# ---------------------------------------------------------------- ajustes
def envolvente_inferior(x, y):
    """Recta y = a + b x que queda POR DEBAJO de todos los puntos y minimiza la suma de residuos.

    Problema: min sum_i (y_i - a - b x_i)  sujeto a  y_i >= a + b x_i.
    Es el estimador de deriva de reloj de Moon, Skelly y Towsley (IEEE INFOCOM 1999), pensado
    para retardos de red: el retardo del host es siempre >= un mínimo y tiene colas largas solo
    hacia arriba, así que el borde inferior de la nube define la recta. La solución es la arista
    de la envolvente convexa inferior que contiene a x = media(x) (maximizar a + b·media(x)).
    Es inmune a las colas hacia arriba; un punto por DEBAJO del resto (host antes que la
    cámara: imposible si el timestamp es bueno) sí lo mueve.

    x, y: [n] (n >= 2, x no todos iguales). Devuelve (a, b) o (nan, nan).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 2 or np.ptp(x) == 0:
        return np.nan, np.nan
    orden = np.argsort(x, kind='stable')
    x0, esc = x[orden[0]], np.ptp(x)
    u = (x[orden] - x0) / esc                       # en [0, 1]: productos cruzados bien escalados
    v = y[orden] - y[orden[0]]
    h = []                                           # cadena monótona de Andrew, solo inferior
    for k in range(len(u)):
        while len(h) >= 2:
            i, j = h[-2], h[-1]
            if (u[j] - u[i]) * (v[k] - v[i]) - (v[j] - v[i]) * (u[k] - u[i]) <= 0:
                h.pop()
            else:
                break
        h.append(k)
    h = np.array(h)
    um = u.mean()
    e = int(np.clip(np.searchsorted(u[h], um, side='right') - 1, 0, len(h) - 2))
    i, j = h[e], h[e + 1]
    if u[j] == u[i]:
        return np.nan, np.nan
    b_u = (v[j] - v[i]) / (u[j] - u[i])
    a_u = v[i] - b_u * u[i]
    b = b_u / esc
    a = a_u + y[orden[0]] - b * x0
    return a, b


def recta_local(x, y, ventana=15):
    """Valor ajustado de una recta de cuadrados mínimos local (método de fluidos).

    Para cada i se ajusta y = a + b x con los 2·ventana+1 puntos i-ventana .. i+ventana (en los
    bordes la ventana se corre para quedar adentro, siempre con el mismo número de puntos) y
    se devuelve ŷ_i = ȳ + b (x_i - x̄), con b = sum (x-x̄)(y-ȳ) / sum (x-x̄)².
    x, y: [n]. Devuelve ŷ [n] (NaN si n < 2).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    if n < 2:
        return np.full(n, np.nan)
    L = min(2 * int(ventana) + 1, n)
    ini = np.clip(np.arange(n) - int(ventana), 0, n - L)
    idx = ini[:, None] + np.arange(L)[None, :]
    X, Y = x[idx], y[idx]
    xm, ym = X.mean(axis=1), Y.mean(axis=1)
    dx = X - xm[:, None]
    sxx = np.sum(dx * dx, axis=1)
    b = np.where(sxx > 0, np.sum(dx * (Y - ym[:, None]), axis=1) / np.where(sxx > 0, sxx, 1), np.nan)
    return ym + b * (x - xm)


# ---------------------------------------------------------------- reconstrucción
def _estadistica_dt(t, saltos):
    """fps (todo el registro, por contador), media/desvío de dt entre cuadros consecutivos
    (salto de contador == 1) y dt máximo entre guardados (incluye huecos)."""
    out = dict(fps_medido=np.nan, periodo_s=np.nan, dt_medio_s=np.nan, dt_desvio_s=np.nan,
               dt_max_s=np.nan)
    if len(t) < 2:
        return out
    dt = np.diff(t)
    ok = np.isfinite(t)
    if ok.sum() >= 2:
        i0, i1 = np.flatnonzero(ok)[[0, -1]]
        dc = np.sum(saltos[i0:i1])                   # cuadros de cámara entre el primero y el último
        if dc > 0 and t[i1] > t[i0]:
            out['periodo_s'] = (t[i1] - t[i0]) / dc
            out['fps_medido'] = 1.0 / out['periodo_s']
    uno = (saltos == 1) & np.isfinite(dt)
    if uno.any():
        out['dt_medio_s'] = float(np.mean(dt[uno]))
        out['dt_desvio_s'] = float(np.std(dt[uno], ddof=1)) if uno.sum() > 1 else np.nan
    if np.isfinite(dt).any():
        out['dt_max_s'] = float(np.nanmax(dt))
    return out


def _con_hardware(tc, th, tick_hz):
    """t [s] desde el primer cuadro a partir de los ticks de hardware.

    Ajuste t_host = a + b·t_camara (envolvente_inferior) sobre toda la grabación:
    b = período del tick medido con el reloj del host → tick_hz_ajustado = 1/b.
    t_i = (t_camara_i - t_camara_0) / tick_hz_usado, con tick_hz_usado = tick_hz si viene,
    si no tick_hz_ajustado. Los cuadros sin t_camara quedan NaN.
    """
    ok = np.isfinite(tc)
    i0 = int(np.flatnonzero(ok)[0])
    x = tc - tc[i0]
    a, b = envolvente_inferior(x[ok], th[ok])        # (nan, nan) si no hay reloj del host
    tick_aj = 1.0 / b if np.isfinite(b) and b > 0 else np.nan
    tick_usado = float(tick_hz) if tick_hz else tick_aj
    t = x / tick_usado
    res = th - (a + b * x)                           # >= 0: demora del host sobre la mínima
    t_aj = x / tick_aj
    info = dict(metodo='hardware', ajuste='envolvente_inferior (Moon et al. 1999)',
                tick_hz_usado=tick_usado, tick_hz_ajustado=tick_aj,
                t_camara_faltantes=int((~ok).sum()),
                t_camara_no_creciente=int(np.sum(np.diff(tc[ok]) <= 0)),
                residuo_max_host_s=_nan_o(np.nanmax, res),
                residuo_mediana_host_s=_nan_o(np.nanmedian, res),
                residuo_p99_host_s=_nan_o(lambda v: np.nanpercentile(v, 99), res),
                diferencia_max_vs_ajuste_s=_nan_o(np.nanmax, np.abs(t - t_aj)))
    if not np.isfinite(th).any():
        info['ajuste'] = 'sin reloj del host: tick_hz nominal, sin ajuste ni ppm'
    return t, info


def _nan_o(f, v):
    """f(v) como float, o NaN si v no tiene ningún valor finito (sin RuntimeWarning)."""
    v = np.asarray(v, float)
    return float(f(v)) if np.isfinite(v).any() else np.nan


AVISO_NOMINAL = ('SIN timestamps reales: t = t_nominal_s (índice / fps nominal). No se miden el '
                 'período real, el jitter ni la deriva del reloj, y NO se detectan cuadros perdidos.')
AVISO_SIN_CONTADOR = ('SIN contador de hardware: los cuadros perdidos NO se pueden detectar '
                      '(cuadros_perdidos = NaN en el informe; 0 en el resultado = ninguno detectado).')


def _nominal(tn):
    """t [s] desde el primer cuadro a partir de la columna t_nominal_s (fps nominal, NO medido)."""
    ok = np.isfinite(tn)
    t = tn - tn[np.flatnonzero(ok)[0]]
    info = dict(metodo='nominal', ajuste='ninguno: t_nominal_s escrito al importar (fps nominal)',
                tick_hz_usado=np.nan, tick_hz_ajustado=np.nan, residuo_max_host_s=np.nan,
                residuo_mediana_host_s=np.nan, residuo_p99_host_s=np.nan,
                diferencia_max_vs_ajuste_s=np.nan, aviso=AVISO_NOMINAL)
    return t, info


def _host_s(th_ns, n):
    """t_host_ns [n] (int, o float con NaN si falta) → s desde el primer valor finito; None → NaN."""
    if th_ns is None:
        return np.full(n, np.nan)
    th_ns = np.asarray(th_ns)
    if n and th_ns.dtype.kind in 'iu':
        return (th_ns - th_ns[0]).astype(float) * 1e-9
    th = th_ns.astype(float)
    ok = np.isfinite(th)
    return (th - th[np.flatnonzero(ok)[0]]) * 1e-9 if ok.any() else th


def _contador(c, n):
    """Columna contador_camara → (int64 [n], hay_contador). Vacía o con NaN = no hay contador:
    se usa el índice 0..n-1 SOLO para las cuentas internas (no detecta perdidos)."""
    c = np.asarray(c)
    if c.dtype.kind in 'iu':
        return c.astype(np.int64), True
    cf = c.astype(float)
    if n and np.isfinite(cf).all():
        return cf.astype(np.int64), True
    return np.arange(n, dtype=np.int64), n == 0


def _sin_hardware(c, th, ventana, ajuste):
    """t [s] desde el primer cuadro sin timestamp de hardware (método de fluidos).

    ajuste='local' : recta_local de t_host contra contador_camara en ±ventana cuadros;
                     t_i = ŷ_i - ŷ_0 (la demora media del host se cancela al restar).
    ajuste='global': una sola recta por envolvente_inferior (supone período de cámara fijo).
    El residuo es t_host - ajuste.
    """
    x = (c - c[0]).astype(float)
    if ajuste == 'global':
        a, b = envolvente_inferior(x, th)
        yh = a + b * x
        nombre = 'envolvente_inferior global contra contador'
    else:
        yh = recta_local(x, th, ventana)
        nombre = f'recta local de cuadrados mínimos, ±{int(ventana)} cuadros'
    t = yh - yh[0]
    res = th - yh
    info = dict(metodo='host', ajuste=nombre, ventana=int(ventana),
                tick_hz_usado=np.nan, tick_hz_ajustado=np.nan,
                residuo_max_host_s=float(np.nanmax(np.abs(res))),
                residuo_mediana_host_s=float(np.nanmedian(np.abs(res))),
                residuo_p99_host_s=float(np.nanpercentile(np.abs(res), 99)),
                diferencia_max_vs_ajuste_s=np.nan)
    return t, info


def reconstruir(tiempos, tick_hz=None, ventana=15, ajuste_sin_hw='local'):
    """Tiempo de cada cuadro guardado, huecos y el informe de tiempo.

    tiempos : dict con 'contador_camara' [n], 't_camara' [n] (ticks; NaN = sin timestamp) y
              't_host_ns' [n] (time.perf_counter_ns al recibir). Columnas de tiempos.csv.
    tick_hz : frecuencia del reloj de la cámara [Hz] informada por el SDK. None = usar la
              ajustada contra el host (ver abajo).
    ventana : semiancho [cuadros] de la recta local cuando no hay timestamp de hardware.
    ajuste_sin_hw : 'local' (método de fluidos, por defecto) o 'global' (ver _sin_hardware).

    Con timestamp de hardware: t_i = (t_camara_i - t_camara_0) / tick_hz_usado. Además se
    ajusta t_host = a + b·t_camara por la envolvente inferior (robusta a jitter y colas del
    host): tick_hz_ajustado = 1/b y ppm_reloj = (tick_hz_ajustado / tick_hz - 1)·1e6 (cuánto
    difiere el reloj de la cámara del nominal, medido con el reloj del host).
    Sin timestamp de hardware (t_camara todo NaN): recta de t_host contra contador_camara.
    Sin timestamp ni reloj del host (grabaciones importadas de otro programa, ver
    pinzas.importar): si hay columna t_nominal_s, metodo = 'nominal' (t = t_nominal_s, con
    informe['aviso'] = AVISO_NOMINAL: no son tiempos medidos); si no, t = NaN ('sin_tiempos').
    contador_camara vacío (NaN) = no hay contador de hardware: no se detectan perdidos
    (informe['cuadros_perdidos'] = NaN, informe['contador'] = 'ninguno', AVISO_SIN_CONTADOR).

    Devuelve dict:
      t         [n] s desde el primer cuadro guardado (NaN = no confiable)
      dt        [n-1] s
      contador  [n] int64 (el índice 0..n-1 si no hay contador de hardware)
      huecos    índices i con cuadros perdidos entre el guardado i y el i+1 (por contador)
      tramos    [m, 2] (i0, i1 inclusive) sin huecos NI anomalías de contador
      cuadros_perdidos  total DETECTADO por contador (int; 0 si no hay contador: ver informe)
      informe   dict JSON-serializable: fps_medido, periodo_s, dt_medio_s, dt_desvio_s
                (entre cuadros con salto de contador 1), dt_max_s (incluye huecos),
                cuadros_perdidos, posiciones_perdidos ([[desde, hasta]] valores de contador
                que faltan), anomalias_contador, patron_0_2, residuo_max_host_s,
                residuo_mediana_host_s, residuo_p99_host_s, metodo, ajuste, tick_hz_usado,
                tick_hz_ajustado, ppm_reloj, diferencia_max_vs_ajuste_s, n_cuadros, contador
                ('hardware' | 'ninguno') y, si corresponde, aviso.
    """
    n = len(tiempos['contador_camara'])
    c, hay_contador = _contador(tiempos['contador_camara'], n)
    th = _host_s(tiempos.get('t_host_ns'), n)
    tc = np.asarray(tiempos.get('t_camara', np.full(n, np.nan)), dtype=float)
    tn = np.asarray(tiempos.get('t_nominal_s', np.full(n, np.nan)), dtype=float)
    ca = analizar_contador(c)

    if n == 0:
        t, info = np.zeros(0), dict(metodo='vacio')
    elif np.isfinite(tc).any():
        t, info = _con_hardware(tc, th, tick_hz)
    elif np.isfinite(th).any():
        t, info = _sin_hardware(c, th, ventana, ajuste_sin_hw)
    elif np.isfinite(tn).any():
        t, info = _nominal(tn)
    else:
        t, info = np.full(n, np.nan), dict(metodo='sin_tiempos', ajuste='ninguno',
                                           aviso='no hay t_camara, t_host_ns ni t_nominal_s')

    nominal = float(tick_hz) if tick_hz else np.nan
    info['ppm_reloj'] = (info.get('tick_hz_ajustado', np.nan) / nominal - 1) * 1e6 \
        if np.isfinite(nominal) else np.nan
    info.update(_estadistica_dt(t, ca['saltos']))
    info['n_cuadros'] = int(n)
    info['cuadros_perdidos'] = int(ca['perdidos'].sum()) if hay_contador else np.nan
    info['posiciones_perdidos'] = [[int(c[i] + 1), int(c[i + 1] - 1)] for i in ca['huecos']]
    info['anomalias_contador'] = [[int(i), int(ca['saltos'][i])] for i in ca['anomalias']]
    info['patron_0_2'] = ca['patron_0_2']
    info['contador'] = 'hardware' if hay_contador else 'ninguno'
    if not hay_contador:
        info['aviso'] = ' '.join(a for a in (info.get('aviso', ''), AVISO_SIN_CONTADOR) if a)
    info = {k: (float(v) if isinstance(v, (float, np.floating)) else v) for k, v in info.items()}

    cortes = np.union1d(ca['huecos'], ca['anomalias'])
    return dict(t=t, dt=np.diff(t), contador=c, huecos=ca['huecos'],
                tramos=tramos_sin_huecos(n, cortes), cuadros_perdidos=int(ca['perdidos'].sum()),
                informe=info)


# ---------------------------------------------------------------- en vivo
def estadistica_en_vivo(contadores, t_camara=None, t_host_ns=None, tick_hz=None, n_ultimos=100):
    """Números livianos para el monitor del grabador (no reemplazan a `reconstruir`).

    contadores : [n] contador de hardware de los cuadros recibidos hasta ahora.
    t_camara   : [n] ticks (o None / NaN); t_host_ns: [n] ns del host (o None).
    tick_hz    : frecuencia nominal del reloj de la cámara [Hz].
    fps = (c[-1] - c[-N]) / (t[-1] - t[-N]) sobre los últimos N = n_ultimos cuadros (cuenta los
    perdidos por contador, así que es la tasa de la CÁMARA); se usa t_camara/tick_hz si está,
    si no t_host. perdidos: total acumulado por contador.
    Devuelve dict: fps, fps_host, fuente, perdidos, anomalias_contador, n, dt_max_s (ventana).
    """
    c = np.asarray(contadores, dtype=np.int64)
    n = len(c)
    out = dict(fps=np.nan, fps_host=np.nan, fuente='', perdidos=0, anomalias_contador=0,
               n=int(n), dt_max_s=np.nan)
    if n < 2:
        return out
    s = np.diff(c)
    out['perdidos'] = int(np.sum(s[s > 1] - 1))
    out['anomalias_contador'] = int(np.sum(s <= 0))
    k = max(0, n - int(n_ultimos))
    dc = c[-1] - c[k]

    def _fps(t):
        dtt = t[-1] - t[k]
        return dc / dtt if dtt > 0 else np.nan

    th = None
    if t_host_ns is not None:
        th = (np.asarray(t_host_ns[k:]) - np.asarray(t_host_ns)[k]).astype(float) * 1e-9
        th = np.concatenate([np.zeros(k), th])          # índices alineados con c
        out['fps_host'] = _fps(th)
    tc = None if t_camara is None else np.asarray(t_camara, dtype=float)
    if tc is not None and tick_hz and np.all(np.isfinite(tc[[k, -1]])):
        tt = (tc - tc[k]) / float(tick_hz)
        out['fps'], out['fuente'] = _fps(tt), 'camara'
        out['dt_max_s'] = float(np.nanmax(np.diff(tt[k:])))
    elif th is not None:
        out['fps'], out['fuente'] = out['fps_host'], 'host'
        out['dt_max_s'] = float(np.max(np.diff(th[k:])))
    return out
