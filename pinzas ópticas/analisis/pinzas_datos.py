"""Grabaciones con la pinza prendida: de trayectorias.csv + platina.csv a tablas listas para analizar.

    from pinzas_datos import *
    a = atrapada('BARRIDO 1')        # una fila por cuadro: la partícula que está en la trampa
    p = patadas('BARRIDO 1')         # una fila por movimiento de la platina
    guardar_tablas('BARRIDO 1')      # escribe datos/BARRIDO 1/atrapada.csv y patadas.csv

Convenciones
- Posiciones en µm con la escala de la retícula (ESCALA_UM). x hacia la derecha, y hacia ABAJO
  (como la imagen).
- La trampa no se mueve en la imagen: se mueve la platina, y con ella TODO lo demás. `platina.csv`
  (TouchDesigner/herramientas/platina.py) da cuánto se corre la imagen de un cuadro al siguiente.
  Ese corrimiento es la velocidad del fluido respecto de la trampa: el arrastre empuja a la
  partícula atrapada en ESE sentido.
- Los valores del motor que figuran en el cuaderno NO se usan para calcular nada: solo para
  emparejar cada patada con su renglón del cuaderno.
"""
import os

import numpy as np
import pandas as pd

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)                                   # pinzas ópticas/
DATOS = os.path.join(RAIZ, 'datos')

ESCALA_UM = 0.09461          # µm por px (retícula de 0.01 mm; ver PINZAS 1)
ERR_ESCALA_UM = 0.0002

# nombre → posición del láser en la imagen (px), cuadros por segundo (por fechas de archivo) y qué es
CATALOGO = {
    'LASER NE':            dict(trampa_px=(496.9, 477.2), fps=10.022, tipo='patadas', eje='x',
                                que='dos patadas (a = 49 y 100 µm/s², paso 200 µm); sin tiempos en el cuaderno'),
    'BARRIDO 1':           dict(trampa_px=(496.9, 477.2), fps=10.021, tipo='patadas', eje='x',
                                que='barrido de aceleración: 15 patadas, paso 100 µm, a = 100 a 420 µm/s²'),
    'MUEVO A MANO':        dict(trampa_px=(477.5, 485.3), fps=10.021, tipo='patadas', eje='x',
                                que='platina movida a mano'),
    'BARRIDO 2':           dict(trampa_px=(477.5, 485.3), fps=10.019, tipo='patadas', eje='x',
                                que='3 patadas, paso 250 µm, a = 450 µm/s²'),
    'LASER 5':             dict(trampa_px=(477.5, 485.3), fps=13.128, tipo='potencia', eje=None,
                                que='barrido de potencia 1 (30 a 0 mW de a 30 s), partícula quieta'),
    'BARRIDO POTENCIAS 2': dict(trampa_px=(477.5, 485.3), fps=12.501, tipo='potencia', eje=None,
                                que='barrido de potencia 2 (30 a < 3 mW de a 30 s), partícula quieta'),
    'BARRIDO STEPSIZE':    dict(trampa_px=(477.5, 485.3), fps=12.502, tipo='patadas', eje='x',
                                que='barrido de paso: 7 patadas, a = 450 µm/s², paso 25 y 50 µm'),
    'BARRIDO STEPSIZE 2':  dict(trampa_px=(477.5, 485.3), fps=12.503, tipo='patadas', eje='x',
                                que='barrido de paso: 13 patadas, a = 250 µm/s², paso 25 a 150 µm'),
    'BARRIDO STEPSIZE 3':  dict(trampa_px=(477.5, 485.3), fps=12.501, tipo='patadas', eje='y',
                                que='no figura en el cuaderno ("ahora nos vamos a mover verticalmente")'),
}


# ---------------------------------------------------------------- lectura
def cargar(nombre):
    """Lo que hay en datos/<nombre>/: tr (trayectorias), tiempo, platina (o None), particulas."""
    c = os.path.join(DATOS, nombre)
    leer = lambda n: pd.read_csv(os.path.join(c, n)) if os.path.exists(os.path.join(c, n)) else None
    return dict(nombre=nombre, carpeta=c, tr=leer('trayectorias.csv'), tiempo=leer('tiempo.csv'),
                platina=leer('platina.csv'), particulas=leer('particulas.csv'), laser=leer('laser.csv'),
                **CATALOGO.get(nombre, {}))


def velocidad_platina(platina):
    """Velocidad de la imagen CENTRADA en cada cuadro [µm/s]: promedio del corrimiento que llega al
    cuadro y del que sale. (El corrimiento i mide entre los cuadros i-1 e i; la posición de la
    partícula en el cuadro i es un promedio durante su exposición.) Devuelve vx, vy [n]."""
    t = platina.t_s.values
    dt = np.median(np.diff(t))
    dx, dy = platina.dx_px.values * ESCALA_UM, platina.dy_px.values * ESCALA_UM
    sig = lambda d: np.r_[d[1:], 0.0]
    return (dx + sig(dx)) / (2 * dt), (dy + sig(dy)) / (2 * dt)


# ---------------------------------------------------------------- la partícula atrapada
def atrapada(nombre, radio_px=45.0, quieta_px=1.5, margen=5, rojo_min=0.1, lejos_px=4.0):
    """Una fila por cuadro con la partícula que está en la trampa (la más cercana al láser a menos
    de radio_px) y el movimiento de la platina en ese cuadro.

    Columnas
      cuadro, t_s
      particula        id de la trayectoria (cambia si el seguimiento se corta); vacío = no hay
      n_candidatas     cuántas partículas seguidas hay a menos de radio_px del láser
      x_um, y_um       posición respecto de la posición de equilibrio de la trampa (x0_px, y0_px)
      r_um             distancia a la posición de equilibrio
      z, colision      contraste de la detección y marca de colisión del seguimiento
      rojo             cuánto rojo dispersa lo que hay en el láser (laser.csv); ≈ 0 = trampa vacía
      en_trampa        True si esa partícula está atrapada (ver abajo)
      plat_dx_um, plat_dy_um   corrimiento de la imagen entre el cuadro anterior y este
      plat_corr        altura del pico de correlación de ese corrimiento
      vx_um_s, vy_um_s velocidad de la imagen centrada en el cuadro
      quieta           True si la platina no se movió en los `margen` cuadros de alrededor
      x0_px, y0_px     posición de equilibrio (px): mediana de los cuadros atrapada y quieta

    en_trampa: hay una partícula seguida cerca del láser, hay mancha roja (rojo > rojo_min; en los
    barridos de potencia no se pide, porque ahí el rojo baja con la potencia) y su trayectoria no
    descansa en otro lado (si la mediana de un tramo con la platina quieta queda a más de lejos_px
    del equilibrio, es una partícula pegada al vidrio que quedó cerca del láser, no la atrapada).
    """
    d = cargar(nombre)
    tr, pl = d['tr'], d['platina']
    tx, ty = d['trampa_px']
    n = len(d['tiempo'])
    tr = tr[tr.x_px.notna()].copy()
    tr['dist'] = np.hypot(tr.x_px - tx, tr.y_px - ty)
    cerca = tr[tr.dist < radio_px]
    eleg = cerca.loc[cerca.groupby('cuadro').dist.idxmin()].set_index('cuadro')
    out = pd.DataFrame(dict(cuadro=np.arange(n), t_s=d['tiempo'].t_s.values))
    out['particula'] = eleg.particula.reindex(out.cuadro).values
    out['n_candidatas'] = cerca.groupby('cuadro').size().reindex(out.cuadro, fill_value=0).values
    xp, yp = eleg.x_px.reindex(out.cuadro).values, eleg.y_px.reindex(out.cuadro).values
    out['z'] = eleg.z.reindex(out.cuadro).values
    out['colision'] = eleg.colision.reindex(out.cuadro).values
    if pl is not None:
        out['plat_dx_um'] = pl.dx_px.values * ESCALA_UM
        out['plat_dy_um'] = pl.dy_px.values * ESCALA_UM
        out['plat_corr'] = pl['corr'].values
        out['vx_um_s'], out['vy_um_s'] = velocidad_platina(pl)
        paso = np.hypot(pl.dx_px.values, pl.dy_px.values)
        mov = paso > quieta_px
        cerca_mov = np.convolve(mov, np.ones(2 * margen + 1), 'same') > 0
        out['quieta'] = ~cerca_mov
    else:
        out['quieta'] = True
    out['rojo'] = d['laser'].rojo.values[:n] if d['laser'] is not None else np.nan
    hay = out.particula.notna().values.copy()
    if d['laser'] is not None and d.get('tipo') != 'potencia':
        hay = hay & (out.rojo.values > rojo_min)
    quieta = out.quieta.values
    for _ in range(2):                                    # equilibrio → sacar tramos que descansan lejos → repetir
        ref = hay & quieta if (hay & quieta).sum() >= 20 else hay
        x0, y0 = (np.nanmedian(xp[ref]), np.nanmedian(yp[ref])) if ref.any() else (tx, ty)
        for pid in pd.unique(out.particula[hay].dropna()):
            s = hay & (out.particula == pid).values
            sq = s & quieta
            if sq.sum() >= 15 and np.hypot(np.nanmedian(xp[sq]) - x0, np.nanmedian(yp[sq]) - y0) > lejos_px:
                hay[s] = False
    out['en_trampa'] = hay
    out['x_um'], out['y_um'] = (xp - x0) * ESCALA_UM, (yp - y0) * ESCALA_UM
    out['r_um'] = np.hypot(out.x_um, out.y_um)
    out['x0_px'], out['y0_px'] = x0, y0
    cols = ['cuadro', 't_s', 'particula', 'n_candidatas', 'x_um', 'y_um', 'r_um', 'z', 'colision', 'rojo', 'en_trampa']
    cols += [c for c in ('plat_dx_um', 'plat_dy_um', 'plat_corr', 'vx_um_s', 'vy_um_s') if c in out]
    return out[cols + ['quieta', 'x0_px', 'y0_px']]


# ---------------------------------------------------------------- patadas de la platina
def salto_por_particulas(tr, c0, c1, tol_px=3.0, max_px=900.0):
    """Corrimiento de la imagen entre los cuadros c0 y c1 medido con las PARTÍCULAS seguidas, sin
    usar la correlación de imágenes: el corrimiento más votado entre todos los pares (partícula en
    c0, partícula en c1). Devuelve (dx, dy [px], votos, n0, n1); NaN si no hay acuerdo."""
    a = tr[(tr.cuadro == c0) & tr.x_px.notna()][['x_px', 'y_px']].values
    b = tr[(tr.cuadro == c1) & tr.x_px.notna()][['x_px', 'y_px']].values
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan, 0, len(a), len(b)
    dif = (b[None, :, :] - a[:, None, :]).reshape(-1, 2)
    dif = dif[np.hypot(dif[:, 0], dif[:, 1]) < max_px]
    if len(dif) == 0:
        return np.nan, np.nan, 0, len(a), len(b)
    votos = (np.hypot(dif[:, None, 0] - dif[None, :, 0], dif[:, None, 1] - dif[None, :, 1]) < tol_px).sum(1)
    k = int(np.argmax(votos))
    cerca = np.hypot(dif[:, 0] - dif[k, 0], dif[:, 1] - dif[k, 1]) < tol_px
    if votos[k] < 2:
        return np.nan, np.nan, int(votos[k]), len(a), len(b)
    return float(dif[cerca, 0].mean()), float(dif[cerca, 1].mean()), int(votos[k]), len(a), len(b)


def patadas(nombre, umbral_px=1.5, union=4, minimo_px=8.0, tabla=None):
    """Una fila por movimiento de la platina (tramo de cuadros seguidos en que la imagen se corre
    más de umbral_px por cuadro; se unen tramos separados por hasta `union` cuadros; se descartan
    los que suman menos de minimo_px).

    Columnas
      n, cuadro_ini, cuadro_fin, t_ini_s, t_fin_s, duracion_s
      dx_um, dy_um          desplazamiento total de la imagen (suma de los corrimientos)
      dx_part_um, dy_part_um, votos   lo mismo medido con las partículas seguidas antes y después
      v_pico_um_s           mayor velocidad centrada |v| durante la patada
      v_media_um_s          desplazamiento / duración
      vx_pico_um_s, vy_pico_um_s   componentes de la velocidad en el cuadro del pico
      corr_min              menor pico de correlación durante la patada
      x_max_um, y_max_um    mayor apartamiento de la partícula atrapada (con signo) durante la
                            patada y hasta 3 cuadros después, en cada eje
      r_max_um              mayor distancia al equilibrio en ese lapso
      antes, despues        ¿había partícula atrapada 1 s antes y entre 1 y 3 s después?
      misma                 ¿la que quedó después es la misma trayectoria que antes?
      confiable             False si la correlación fue mala (corr_min < 0.2) y las partículas
                            seguidas no confirman el desplazamiento: no es un movimiento de platina
    """
    d = cargar(nombre)
    pl, tr = d['platina'], d['tr']
    a = atrapada(nombre) if tabla is None else tabla
    t = pl.t_s.values
    dt = np.median(np.diff(t))
    paso = np.hypot(pl.dx_px.values, pl.dy_px.values)
    idx = np.flatnonzero(paso > umbral_px)
    tramos = []
    if len(idx):
        ini = ant = idx[0]
        for i in idx[1:]:
            if i - ant > union:
                tramos.append((ini, ant))
                ini = i
            ant = i
        tramos.append((ini, ant))
    vx, vy = velocidad_platina(pl)
    v = np.hypot(vx, vy)
    filas = []
    for c0, c1 in tramos:
        s = slice(c0, c1 + 1)
        dx, dy = pl.dx_px.values[s].sum(), pl.dy_px.values[s].sum()
        if np.hypot(dx, dy) < minimo_px:
            continue
        k = c0 - 1 + int(np.argmax(v[c0 - 1:c1 + 1]))             # el pico centrado puede caer un cuadro antes
        px, py, votos, _, _ = salto_por_particulas(tr, max(c0 - 2, 0), min(c1 + 2, len(t) - 1))
        lapso = a.iloc[max(c0 - 1, 0):c1 + 4]
        f1 = int(round(1 / dt))
        prev = a.iloc[max(c0 - f1, 0):max(c0 - 1, 1)]
        post = a.iloc[c1 + f1:c1 + 3 * f1]
        antes = prev.en_trampa.mean() > 0.5 if len(prev) else False
        desp = post.en_trampa.mean() > 0.5 if len(post) else False
        misma = bool(antes and desp and prev.particula.mode().iloc[0] == post.particula.mode().iloc[0])

        def extremo(col):
            w = lapso[col].dropna()
            return float(w.loc[w.abs().idxmax()]) if len(w) else np.nan
        filas.append(dict(cuadro_ini=c0, cuadro_fin=c1, t_ini_s=t[c0], t_fin_s=t[c1], duracion_s=(c1 - c0 + 1) * dt,
                          dx_um=dx * ESCALA_UM, dy_um=dy * ESCALA_UM,
                          dx_part_um=px * ESCALA_UM, dy_part_um=py * ESCALA_UM, votos=votos,
                          v_pico_um_s=v[k], v_media_um_s=np.hypot(dx, dy) * ESCALA_UM / ((c1 - c0 + 1) * dt),
                          vx_pico_um_s=vx[k], vy_pico_um_s=vy[k], corr_min=pl['corr'].values[s].min(),
                          x_max_um=extremo('x_um'), y_max_um=extremo('y_um'),
                          r_max_um=float(lapso.r_um.max()) if lapso.r_um.notna().any() else np.nan,
                          antes=bool(antes), despues=bool(desp), misma=misma,
                          confiable=bool(pl['corr'].values[s].min() >= 0.2 or
                                         (np.isfinite(px) and np.hypot(px - dx, py - dy) < 0.2 * np.hypot(dx, dy) + 3))))
    out = pd.DataFrame(filas)
    if len(out):
        out.insert(0, 'n', np.arange(1, len(out) + 1))
    return out


def emparejar_cuaderno(p, nombre, tolerancia_s=4.0):
    """Agrega a la tabla de patadas el renglón del cuaderno más cercano en tiempo (si hay uno a
    menos de tolerancia_s): n_cuaderno, a_nominal_um_s2, paso_nominal_um, resultado_cuaderno."""
    cu = pd.read_csv(os.path.join(DATOS, 'cuaderno', 'patadas_cuaderno.csv'))
    cu = cu[(cu.grabacion == nombre) & cu.t_cuaderno_s.notna()]
    p = p.copy()
    for c in ('n_cuaderno', 't_cuaderno_s', 'a_nominal_um_s2', 'paso_nominal_um', 'resultado_cuaderno'):
        p[c] = np.nan if c != 'resultado_cuaderno' else None
    # todos los pares (renglón, movimiento) a menos de la tolerancia, del más cercano al más lejano:
    # cada renglón y cada movimiento se usan una sola vez
    pares = sorted((abs(p.t_ini_s[k] - r.t_cuaderno_s), i, k) for i, r in cu.iterrows() for k in p.index
                   if abs(p.t_ini_s[k] - r.t_cuaderno_s) <= tolerancia_s)
    usados_i, usados_k = set(), set()
    for _, i, k in pares:
        if i in usados_i or k in usados_k:
            continue
        usados_i.add(i), usados_k.add(k)
        r = cu.loc[i]
        p.loc[k, ['n_cuaderno', 't_cuaderno_s', 'a_nominal_um_s2', 'paso_nominal_um']] = [
            r.n, r.t_cuaderno_s, r.a_nominal_um_s2, r.paso_nominal_um]
        p.loc[k, 'resultado_cuaderno'] = r.resultado_cuaderno if isinstance(r.resultado_cuaderno, str) else None
    return p


def guardar_tablas(nombre):
    """Escribe datos/<nombre>/atrapada.csv y, si hay platina.csv, patadas.csv. Devuelve (a, p)."""
    a = atrapada(nombre)
    a.to_csv(os.path.join(DATOS, nombre, 'atrapada.csv'), index=False, float_format='%.5g')
    p = None
    if cargar(nombre)['platina'] is not None:
        p = emparejar_cuaderno(patadas(nombre, tabla=a), nombre)
        p.to_csv(os.path.join(DATOS, nombre, 'patadas.csv'), index=False, float_format='%.5g')
    return a, p
