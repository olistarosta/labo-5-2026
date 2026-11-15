"""Fuerza de la trampa por arrastre: qué le pasa a la partícula atrapada cuando se mueve la platina.

    from pinzas_fuerza import *
    s = seguir('BARRIDO 1')              # la partícula atrapada cuadro a cuadro, también cuando se aparta
    f = intervalos('BARRIDO 1')          # una fila por par de cuadros consecutivos alrededor de cada patada
    guardar_fuerza('BARRIDO 1')          # escribe datos/BARRIDO 1/seguida.csv y fuerza_intervalos.csv

La física (número de Reynolds ~1e-4: sin inercia). Sobre la partícula actúan la trampa y el arrastre
del fluido, que se mueve con la platina a velocidad v_f:

    0 = F_trampa(x) + γ (v_f − dx/dt)        →        F_trampa(x) / γ = dx/dt − v_f

Entre dos cuadros consecutivos se miden las dos cosas: cuánto se movió la partícula (Δx, de
trampa_detecciones.csv) y cuánto se corrió la imagen (s, de platina.csv). Entonces

    u = (Δx − s) / dt        [µm/s]

es la fuerza de la trampa dividida γ, promediada en ese intervalo, y se asigna a la posición media
de la partícula en el intervalo. u tiene unidades de velocidad: es la velocidad con la que la
trampa sola movería a la partícula en fluido quieto. Para pasar a pN hace falta γ (ver PINZAS 3).

Límite: con 10 a 12.5 cuadros por segundo y patadas de 0.2 a 0.5 s, la partícula recorre varios µm
por intervalo, así que u es un promedio sobre un tramo del pozo, no la fuerza en un punto.
"""
import os

import numpy as np
import pandas as pd

from pinzas_datos import CATALOGO, DATOS, ESCALA_UM, cargar, patadas, atrapada, emparejar_cuaderno

SALTO_MAX_PX = 50.0          # lo más que puede moverse la partícula seguida de un cuadro al siguiente
ANCLA_PX = 8.0               # una detección a menos de esto del equilibrio y con mancha roja ES la atrapada
LEJOS_PX = 110.0             # si se aleja más que esto del equilibrio, se deja de seguir
ROJO_MIN = 0.3


def seguir(nombre, tabla=None, detecciones='trampa_detecciones.csv'):
    """La partícula atrapada cuadro a cuadro, usando las detecciones sin enlazar de la ventana de la
    trampa (trampa_detecciones.csv). Regla:
      1. si hay una detección a menos de ANCLA_PX del equilibrio y hay mancha roja, es la atrapada;
      2. si no, se sigue a la que venía (la más cercana a su posición anterior, a menos de SALTO_MAX_PX);
      3. si no hay ninguna, o se alejó más de LEJOS_PX, se la da por perdida.
    Devuelve una fila por cuadro: cuadro, t_s, x_um, y_um (respecto del equilibrio), z,
    estado ('anclada' = regla 1, 'seguida' = regla 2, '' = no hay), rojo, plat_dx_um, plat_dy_um, plat_corr."""
    d = cargar(nombre)
    a = atrapada(nombre) if tabla is None else tabla
    x0, y0 = float(a.x0_px.iloc[0]), float(a.y0_px.iloc[0])
    det = pd.read_csv(os.path.join(DATOS, nombre, detecciones))
    por_cuadro = {c: g[['x_px', 'y_px', 'z']].values for c, g in det.groupby('cuadro')}
    n = len(a)
    rojo = a.rojo.values
    x, y, z = np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan)
    estado = np.array([''] * n, dtype=object)
    prev = None
    for c in range(n):
        cand = por_cuadro.get(c)
        if cand is None:
            prev = None
            continue
        r0 = np.hypot(cand[:, 0] - x0, cand[:, 1] - y0)
        k = int(np.argmin(r0))
        if r0[k] < ANCLA_PX and rojo[c] > ROJO_MIN:
            x[c], y[c], z[c], estado[c] = cand[k, 0], cand[k, 1], cand[k, 2], 'anclada'
            prev = (cand[k, 0], cand[k, 1])
            continue
        if prev is None:
            continue
        rp = np.hypot(cand[:, 0] - prev[0], cand[:, 1] - prev[1])
        k = int(np.argmin(rp))
        if rp[k] > SALTO_MAX_PX or r0[k] > LEJOS_PX:
            prev = None
            continue
        x[c], y[c], z[c], estado[c] = cand[k, 0], cand[k, 1], cand[k, 2], 'seguida'
        prev = (cand[k, 0], cand[k, 1])
    out = pd.DataFrame(dict(cuadro=np.arange(n), t_s=a.t_s.values, x_um=(x - x0) * ESCALA_UM, y_um=(y - y0) * ESCALA_UM,
                            z=z, estado=estado, rojo=rojo))
    for c in ('plat_dx_um', 'plat_dy_um', 'plat_corr'):
        out[c] = a[c].values
    return out


def intervalos(nombre, antes=3, despues=15, seguida=None, tabla_patadas=None):
    """Una fila por par de cuadros consecutivos (i, i+1) alrededor de cada patada confiable que
    empezó con una partícula atrapada.

    Columnas
      grabacion, patada (n de patadas.csv), cuadro (i), t_s, dt_s (tiempo entre los dos cuadros)
      e_x, e_y          dirección de la patada (versor del desplazamiento total de la imagen)
      xi_um             posición media de la partícula en el intervalo, a lo largo de la patada
      eta_um            ídem, perpendicular
      xi_ini_um, xi_fin_um   posición a lo largo en el cuadro i y en el i+1
      vf_um_s           velocidad del fluido en el intervalo, a lo largo (corrimiento / dt)
      u_par_um_s        (Δxi − corrimiento) / dt: fuerza de la trampa / γ a lo largo (ver arriba)
      u_perp_um_s       ídem perpendicular
      fase              'antes', 'durante' (la platina se mueve) o 'despues'
      estado_ini, estado_fin   'anclada' o 'seguida' en cada extremo
    """
    s = seguir(nombre) if seguida is None else seguida
    p = patadas(nombre) if tabla_patadas is None else tabla_patadas
    dt = float(np.median(np.diff(s.t_s.values)))
    filas = []
    for r in p.itertuples():
        if not (r.confiable and r.antes):
            continue
        norma = np.hypot(r.dx_um, r.dy_um)
        ex, ey = r.dx_um / norma, r.dy_um / norma
        for i in range(max(r.cuadro_ini - 1 - antes, 0), min(r.cuadro_fin + despues, len(s) - 1)):
            a, b = s.iloc[i], s.iloc[i + 1]
            if not (np.isfinite(a.x_um) and np.isfinite(b.x_um)):
                continue
            xi_a, xi_b = a.x_um * ex + a.y_um * ey, b.x_um * ex + b.y_um * ey
            et_a, et_b = -a.x_um * ey + a.y_um * ex, -b.x_um * ey + b.y_um * ex
            sp = b.plat_dx_um * ex + b.plat_dy_um * ey                # corrimiento entre i e i+1
            sn = -b.plat_dx_um * ey + b.plat_dy_um * ex
            fase = 'antes' if i + 1 < r.cuadro_ini else ('durante' if i + 1 <= r.cuadro_fin else 'despues')
            filas.append(dict(grabacion=nombre, patada=r.n, cuadro=i, t_s=a.t_s, dt_s=dt, e_x=ex, e_y=ey,
                              xi_um=0.5 * (xi_a + xi_b), eta_um=0.5 * (et_a + et_b), xi_ini_um=xi_a, xi_fin_um=xi_b,
                              vf_um_s=sp / dt, u_par_um_s=(xi_b - xi_a - sp) / dt, u_perp_um_s=(et_b - et_a - sn) / dt,
                              fase=fase, estado_ini=a.estado, estado_fin=b.estado))
    return pd.DataFrame(filas)


def resumen_patadas(nombre, seguida=None, tabla_patadas=None):
    """patadas.csv con lo que hizo la partícula, medido con el seguimiento de la ventana de la trampa:
      xi_max_um         mayor apartamiento a lo largo de la patada, durante ella y hasta 1 s después
      volvio            ¿está de nuevo anclada entre 1 y 3 s después?
      t_vuelta_s        cuánto tardó en volver a estar anclada desde el fin de la patada (NaN si no volvió)
      resultado         'no se aparto' (xi_max < 0.3 µm), 'se aparto y volvio', 'escapo', o
                        'sin particula' (no había ninguna atrapada antes)"""
    s = seguir(nombre) if seguida is None else seguida
    p = (patadas(nombre) if tabla_patadas is None else tabla_patadas).copy()
    dt = float(np.median(np.diff(s.t_s.values)))
    f1 = int(round(1 / dt))
    xi_max, volvio, t_vuelta, res = [], [], [], []
    for r in p.itertuples():
        norma = np.hypot(r.dx_um, r.dy_um)
        ex, ey = r.dx_um / norma, r.dy_um / norma
        w = s.iloc[max(r.cuadro_ini - 1, 0):r.cuadro_fin + f1]
        xi = (w.x_um * ex + w.y_um * ey).dropna()
        xm = float(xi.loc[xi.abs().idxmax()]) if len(xi) else np.nan
        post = s.iloc[r.cuadro_fin + f1:r.cuadro_fin + 3 * f1]
        v = bool(len(post) and (post.estado == 'anclada').mean() > 0.5)
        desp = s.iloc[r.cuadro_fin + 1:r.cuadro_fin + 3 * f1]
        anc = desp.index[(desp.estado == 'anclada').values]
        tv = (anc[0] - r.cuadro_fin) * dt if v and len(anc) else np.nan
        xi_max.append(xm), volvio.append(v), t_vuelta.append(tv)
        res.append('sin particula' if not r.antes else ('escapo' if not v else
                                                        ('no se aparto' if abs(xm) < 0.3 else 'se aparto y volvio')))
    p['xi_max_um'], p['volvio'], p['t_vuelta_s'], p['resultado'] = xi_max, volvio, t_vuelta, res
    return p


# ---------------------------------------------------------------- modelo del pozo
def fuerza_modelo(xi, k, w):
    """Fuerza de un pozo gaussiano dividida γ [µm/s]:  u(ξ) = −k ξ exp(−ξ² / 2w²).
    Potencial U(ξ) = −γ k w² exp(−ξ²/2w²): cerca del centro es un resorte de constante κ = γ k;
    la fuerza es máxima en ξ = w y vale γ k w e^(−1/2); la profundidad del pozo es γ k w².
    k [1/s] = κ/γ = 1/τ_c (la inversa del tiempo de relajación);  w [µm] = ancho del pozo."""
    xi = np.asarray(xi, float)
    return -k * xi * np.exp(-xi ** 2 / (2 * w ** 2))


def avanzar(xi0, vf, dt, k, w, pasos=40):
    """Dónde termina la partícula después de dt si arranca en xi0 y el fluido va a velocidad vf
    (constante en el intervalo):  dξ/dt = vf + u(ξ).  Runge-Kutta 4. Todo vectorizado."""
    xi = np.array(xi0, float)
    vf, h = np.asarray(vf, float), np.asarray(dt, float) / pasos
    g = lambda x: vf + fuerza_modelo(x, k, w)
    for _ in range(pasos):
        k1 = g(xi)
        k2 = g(xi + 0.5 * h * k1)
        k3 = g(xi + 0.5 * h * k2)
        k4 = g(xi + h * k3)
        xi = xi + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return xi


def ajustar_pozo(f, k0=35.0, w0=1.2, escala_um=0.3, n_boot=200, semilla=0):
    """Ajusta k y w del pozo gaussiano a los intervalos de `f` (tabla de `intervalos`): para cada
    intervalo se integra el movimiento desde su posición inicial y se compara dónde debería terminar
    la partícula con dónde terminó. Así el ajuste no supone que la fuerza es constante en el
    intervalo (la partícula recorre varios µm). Pérdida robusta (soft-L1, escala_um): un intervalo
    mal seguido no manda. Error: bootstrap remuestreando PATADAS enteras (los intervalos de una
    misma patada no son independientes); se informa el medio ancho del intervalo 16–84 %.
    Sesgo medido con patadas sintéticas (PINZAS 3, sección 3): k sale hasta 5 % por debajo del
    verdadero, w hasta 4 % y u_max entre 4 y 8 %, porque la posición de cada cuadro es un promedio
    durante la exposición y la velocidad del fluido se toma constante dentro del intervalo.
    Devuelve dict: k, w, err_k, err_w, u_max (= k w e^-1/2), err_u_max, k_w2 (= k w², profundidad/γ),
    err_k_w2, n_intervalos, n_patadas, residuo_um (desvío robusto de los residuos)."""
    from scipy.optimize import least_squares
    f = f[np.isfinite(f.xi_ini_um) & np.isfinite(f.xi_fin_um)]
    datos = f[['xi_ini_um', 'vf_um_s', 'dt_s', 'xi_fin_um']].values
    grupos = (f.grabacion.astype(str) + '#' + f.patada.astype(str)).values

    def ajustar(d):
        res = lambda p: avanzar(d[:, 0], d[:, 1], d[:, 2], np.exp(p[0]), np.exp(p[1])) - d[:, 3]
        return least_squares(res, [np.log(k0), np.log(w0)], loss='soft_l1', f_scale=escala_um,
                             bounds=([np.log(1.0), np.log(0.2)], [np.log(500.0), np.log(6.0)]))      # k en 1/s, w en µm

    r = ajustar(datos)
    k, w = np.exp(r.x)
    rng = np.random.default_rng(semilla)
    unicos = np.unique(grupos)
    indices = {g: np.flatnonzero(grupos == g) for g in unicos}
    bk, bw = [], []
    for _ in range(n_boot if len(unicos) >= 3 else 0):
        sel = np.concatenate([indices[g] for g in rng.choice(unicos, len(unicos))])
        try:
            rb = ajustar(datos[sel])
            bk.append(np.exp(rb.x[0])), bw.append(np.exp(rb.x[1]))
        except Exception:
            pass
    bk, bw = np.array(bk), np.array(bw)
    # medio ancho del intervalo 16–84 %: como un desvío estándar, pero no lo arruina una réplica que diverge
    err = lambda v: float(0.5 * (np.percentile(v, 84) - np.percentile(v, 16))) if len(v) > 10 else np.nan
    res = avanzar(datos[:, 0], datos[:, 1], datos[:, 2], k, w) - datos[:, 3]
    return dict(k=float(k), w=float(w), err_k=err(bk), err_w=err(bw),
                u_max=float(k * w * np.exp(-0.5)), err_u_max=err(bk * bw * np.exp(-0.5)),
                k_w2=float(k * w ** 2), err_k_w2=err(bk * bw ** 2), n_intervalos=len(datos), n_patadas=len(unicos),
                residuo_um=float(1.4826 * np.median(np.abs(res - np.median(res)))))


def guardar_fuerza(nombre):
    """Escribe datos/<nombre>/seguida.csv, fuerza_intervalos.csv y patadas.csv (con el resultado).
    Devuelve (seguida, intervalos, patadas)."""
    a = atrapada(nombre)
    p0 = patadas(nombre, tabla=a)
    s = seguir(nombre, tabla=a)
    f = intervalos(nombre, seguida=s, tabla_patadas=p0)
    p = emparejar_cuaderno(resumen_patadas(nombre, seguida=s, tabla_patadas=p0), nombre)
    c = os.path.join(DATOS, nombre)
    a.to_csv(os.path.join(c, 'atrapada.csv'), index=False, float_format='%.5g')
    s.to_csv(os.path.join(c, 'seguida.csv'), index=False, float_format='%.5g')
    f.to_csv(os.path.join(c, 'fuerza_intervalos.csv'), index=False, float_format='%.5g')
    p.to_csv(os.path.join(c, 'patadas.csv'), index=False, float_format='%.5g')
    return s, f, p
