"""Movimiento browniano de las partículas LIBRES de todas las grabaciones, con y sin láser.

    from browniano_todas import *
    for nombre in grabaciones():
        if faltan(nombre): continue
        g = cargar(nombre)                 # posiciones [µm] solo en los cuadros que se pueden usar
        tabla = por_particula(g)           # una fila por partícula: D, radio aparente, distancia al láser…

Lo usa `1 - movimiento browniano/todas las grabaciones por radio.ipynb`.

Ideas
- Cuadros que se usan. Con láser la platina da patadas: se usan solo los TRAMOS en que está quieta
  (`tramos_quietos`, a partir de platina.csv) y nunca se calcula un desplazamiento que cruce de un
  tramo a otro. También se sacan los cuadros con colisión y las pausas de escritura del disco.
- D sin que entre el arrastre. Cada tramo se parte en TROZOS de unos segundos. Para cada retardo
  (1, 2, 3, 4 cuadros) se calcula la VARIANZA de los desplazamientos alrededor de su promedio dentro
  del trozo: el promedio es el arrastre del líquido, y al restarlo no entra. La varianza crece en
  línea recta con el retardo: V(τ) = 2·D·τ + 2·ruido². D sale de la pendiente (el ruido de
  localización solo corre la recta hacia arriba). Es la misma idea que `potencia.difusion_libre`,
  con la corrección exacta por haber restado el promedio (`_suma_y_grados`): sin ella D sale unos
  7 % de menos con trozos de 6 s.
- El error de D de cada partícula es por remuestreo de los trozos (bootstrap): se repite la cuenta
  muchas veces eligiendo trozos al azar y se mira cuánto cambia D.

Unidades: µm, s, µm²/s. x hacia la derecha, y hacia abajo (como la imagen).
"""
import os

import numpy as np
import pandas as pd

from pinzas_datos import CATALOGO, DATOS, ESCALA_UM

# ---------------------------------------------------------------- qué grabaciones hay
SIN_LASER = ('BROWNIANO 1', 'BROWNIANO 2')
DIA_19 = tuple('BARRIDO POTENCIAS %d' % i for i in range(3, 9))          # las del 19/11, por si el catálogo no trae el día
# pausas de escritura del disco (datos/catalogo.csv): cada archivo TIFF trae 167 cuadros; en estos
# archivos faltan cuadros y no se sabe cuáles. Se descartan enteros: (primer cuadro, último cuadro).
PAUSAS = {'BROWNIANO 1': [(26 * 167, 28 * 167 - 1)],
          'BARRIDO 1': [(17 * 167, 18 * 167 - 1)],
          'BARRIDO STEPSIZE 3': [(34 * 167, 35 * 167 - 1)]}
NECESARIOS = ('trayectorias.csv', 'tiempo.csv', 'radios.csv')
NECESARIOS_PATADAS = ('platina.csv',)        # donde la platina da patadas, sin platina.csv no se puede saber cuándo

LAGS_D = (1, 2, 3, 4)        # retardos (en cuadros) de los que sale D
TROZO_D_S = 6.0              # largo de los trozos para D [s]
LAGS_MSD = tuple(range(1, 21))
TROZO_MSD_S = 15.0           # trozos más largos para mirar el MSD hasta 20 cuadros


def grabaciones():
    """Las dos sin láser y todas las del CATALOGO (con láser), en ese orden."""
    return list(SIN_LASER) + list(CATALOGO)


def dia(nombre):
    """Día de la grabación, 'dd/mm': de la columna `dia` de datos/catalogo.csv; si no está ahí,
    19/11 para las de DIA_19 y 12/11 para el resto."""
    ruta = os.path.join(DATOS, 'catalogo.csv')
    if os.path.exists(ruta):
        cat = pd.read_csv(ruta)
        if 'dia' in cat:
            d = cat.loc[cat.grabacion == nombre, 'dia']
            if len(d) and isinstance(d.iloc[0], str):
                return d.iloc[0][:5]
    return '19/11' if nombre in DIA_19 else '12/11'


def muestra(nombre):
    """Grupo de la grabación: su día, salvo BROWNIANO 1, que es otra mezcla."""
    if nombre == 'BROWNIANO 1':
        return '12/11, muestra de BROWNIANO 1'
    return dia(nombre)


def apagado_s(nombre):
    """Segundo en que se apagó el láser en esa grabación (primer escalón con potencia_mW = 0 en
    datos/cuaderno/potencias_cuaderno.csv), o None si no se apagó o no figura."""
    ruta = os.path.join(DATOS, 'cuaderno', 'potencias_cuaderno.csv')
    if not os.path.exists(ruta):
        return None
    cu = pd.read_csv(ruta)
    cu = cu[(cu.grabacion == nombre) & (cu.potencia_mW == 0)]
    return float(cu.t_ini_s.min()) if len(cu) else None


def faltan(nombre):
    """Archivos de datos/<nombre>/ que hacen falta y todavía no están (lista vacía = se puede analizar).
    En los barridos de potencia la platina no se toca: platina.csv sirve para verificarlo, pero no
    es obligatorio (si falta, se usan todos los cuadros y `cargar` deja platina = None)."""
    pedir = NECESARIOS + (NECESARIOS_PATADAS if CATALOGO.get(nombre, {}).get('tipo') == 'patadas' else ())
    return [a for a in pedir if not os.path.exists(os.path.join(DATOS, nombre, a))]


# ---------------------------------------------------------------- cuadros que se pueden usar
def tramos_quietos(platina, umbral_px=1.5, union=4, minimo_px=8.0, salto_px=3.0, margen=5):
    """A partir de platina.csv (cuánto se corre toda la imagen de un cuadro al siguiente) dice en
    qué cuadros la platina está quieta. Una PATADA es una seguidilla de cuadros en que la imagen
    se corre más de umbral_px por cuadro (se unen las separadas por hasta `union` cuadros) y que
    recorre en total más de minimo_px, o que tiene algún paso de más de salto_px. Es el criterio de
    `pinzas_datos.patadas`, pero más estricto: acá se suma el camino recorrido y no el desplazamiento
    neto (un ida y vuelta también es un movimiento de la platina) y un solo paso grande alcanza (con
    la platina quieta el ruido de la medición del corrimiento no pasa de 3 px).
    Se descartan además `margen` cuadros antes y después de cada una.
    Devuelve (movida [n] bool, lista de patadas (cuadro_ini, cuadro_fin, dx_px, dy_px))."""
    dx, dy = platina.dx_px.fillna(0).values, platina.dy_px.fillna(0).values
    idx = np.flatnonzero(np.hypot(dx, dy) > umbral_px)
    grupos = []
    if len(idx):
        ini = ant = idx[0]
        for i in idx[1:]:
            if i - ant > union:
                grupos.append((ini, ant))
                ini = i
            ant = i
        grupos.append((ini, ant))
    movida = np.zeros(len(dx), bool)
    pat = []
    for c0, c1 in grupos:
        sx, sy = dx[c0:c1 + 1].sum(), dy[c0:c1 + 1].sum()
        paso = np.hypot(dx[c0:c1 + 1], dy[c0:c1 + 1])
        if paso.sum() >= minimo_px or paso.max() >= salto_px:
            pat.append((int(c0), int(c1), float(sx), float(sy)))
            movida[max(c0 - margen, 0):c1 + margen + 1] = True
    return movida, pat


def _numerar_tramos(usable):
    """usable [n] bool → tramo [n] int: 0, 1, 2… para cada seguidilla de cuadros usables; -1 = no se usa."""
    tramo = np.full(len(usable), -1)
    ini = np.flatnonzero(usable & ~np.r_[False, usable[:-1]])
    fin = np.flatnonzero(usable & ~np.r_[usable[1:], False])
    for k, (a, b) in enumerate(zip(ini, fin)):
        tramo[a:b + 1] = k
    return tramo


def _trozos(tramo, n_trozo):
    """Parte cada tramo en trozos iguales de alrededor de n_trozo cuadros → trozo [n] int (-1 = no se usa)."""
    trozo = np.full(len(tramo), -1)
    k = 0
    for tr in np.unique(tramo[tramo >= 0]):
        idx = np.flatnonzero(tramo == tr)
        partes = max(1, int(round(len(idx) / n_trozo)))
        for parte in np.array_split(idx, partes):
            trozo[parte] = k
            k += 1
    return trozo


# ---------------------------------------------------------------- lectura
def cargar(nombre):
    """Lee una grabación y deja las posiciones solo en los cuadros que se pueden usar.

    Devuelve dict:
      nombre, muestra, dia ('dd/mm'), laser (bool), dt [s], n (cuadros), trampa_px (o None)
      ids [nb]            las trayectorias con radio medido (radios.csv: las de 100 cuadros o más)
      X, Y [n, nb] µm     NaN = no medido, con colisión, platina moviéndose o pausa de escritura
      Xc, Yc [n, nb] px   todas las posiciones medidas (para la distancia al láser)
      tramo [n]           número de tramo quieto de cada cuadro (-1 = no se usa)
      patadas             lista (cuadro_ini, cuadro_fin, dx_px, dy_px)
      radios              radios.csv con índice `particula`
      atrapadas           ids que figuran en la trampa en atrapada.csv (vacío si no está el archivo)
      platina             platina.csv o None
      cuenta              dict: cuántas trayectorias y posiciones había y cuántas se sacaron
    """
    c = os.path.join(DATOS, nombre)
    leer = lambda a: pd.read_csv(os.path.join(c, a)) if os.path.exists(os.path.join(c, a)) else None
    tr, tiempo, radios = leer('trayectorias.csv'), leer('tiempo.csv'), leer('radios.csv')
    platina, atr = leer('platina.csv'), leer('atrapada.csv')
    laser = nombre in CATALOGO
    n = len(tiempo)
    dt = float(tiempo.t_s.values[-1] - tiempo.t_s.values[0]) / (n - 1)       # de punta a punta: t_s trae solo 6 cifras
    pat = []
    if platina is not None:
        platina = platina.set_index('cuadro').reindex(range(n)).reset_index()      # una fila por cuadro, sí o sí
        movida, pat = tramos_quietos(platina)
        usable = ~movida
    else:
        usable = np.ones(n, bool)
    for a, b in PAUSAS.get(nombre, []):
        usable[a:b + 1] = False
    tramo = _numerar_tramos(usable)

    tr = tr[tr.x_px.notna()]
    radios = radios.set_index('particula')
    ids = np.array(sorted(set(radios.index) & set(tr.particula.unique())))
    col = {p: k for k, p in enumerate(ids)}
    tr = tr[tr.particula.isin(col)]
    k = tr.particula.map(col).values
    f = tr.cuadro.values
    Xc = np.full((n, len(ids)), np.nan)
    Yc = np.full((n, len(ids)), np.nan)
    Xc[f, k], Yc[f, k] = tr.x_px.values, tr.y_px.values
    bien = np.zeros((n, len(ids)), bool)
    bien[f, k] = tr.colision.values == 0
    bien &= usable[:, None]
    X = np.where(bien, Xc * ESCALA_UM, np.nan)
    Y = np.where(bien, Yc * ESCALA_UM, np.nan)
    atrapadas = set()
    if atr is not None:
        atrapadas = set(atr[atr.en_trampa & atr.particula.notna()].particula.astype(int))
    medidas = np.isfinite(Xc)
    cuenta = dict(trayectorias=int(pd.read_csv(os.path.join(c, 'particulas.csv')).shape[0]), con_radio=len(ids),
                  posiciones=int(medidas.sum()), con_colision=int((medidas & (~bien) & usable[:, None]).sum()),
                  con_platina_movida=int((medidas & ~usable[:, None]).sum()),
                  cuadros_usables=int(usable.sum()))
    return dict(nombre=nombre, muestra=muestra(nombre), dia=dia(nombre), laser=laser, dt=dt, n=n,
                trampa_px=CATALOGO[nombre]['trampa_px'] if laser else None, ids=ids, X=X, Y=Y, Xc=Xc, Yc=Yc,
                tramo=tramo, patadas=pat, radios=radios.reindex(ids), atrapadas=atrapadas, platina=platina,
                tiene_atrapada=atr is not None, cuenta=cuenta)


# ---------------------------------------------------------------- el estimador de D
def _suma_y_grados(v, trozo, lags):
    """Para una serie v [n] (NaN = no medido) y cada retardo L: en cada trozo, la suma de
    (d − promedio del trozo)² de los desplazamientos d = v[i+L] − v[i] que caben enteros en el trozo,
    y sus «grados de libertad»: por cuánto hay que dividir esa suma para que dé, en promedio, la
    varianza verdadera del desplazamiento de un browniano. Sin huecos son (cantidad − L): restar el
    promedio se lleva más que 1 porque los desplazamientos de retardo L se solapan.
    H es lo mismo para el error de localización (que no se solapa igual): en promedio
    S = 2·D·dt·L·G + 2·ruido²·H.
    Devuelve S, G, H [trozos, len(lags)] y los desplazamientos centrados de retardo lags[0] con su trozo."""
    nt = int(trozo.max()) + 1 if (trozo >= 0).any() else 0
    S = np.zeros((nt, len(lags)))
    G = np.zeros((nt, len(lags)))
    H = np.zeros((nt, len(lags)))
    primero = None
    for b, L in enumerate(lags):
        d = v[L:] - v[:-L]
        ok = np.isfinite(d) & (trozo[L:] == trozo[:-L]) & (trozo[:-L] >= 0)
        if not ok.any():
            continue
        tz = trozo[:-L][ok]
        dd = d[ok]
        cant = np.bincount(tz, minlength=nt).astype(float)
        media = np.bincount(tz, dd, minlength=nt) / np.maximum(cant, 1)
        cen = dd - media[tz]
        S[:, b] = np.bincount(tz, cen ** 2, minlength=nt)
        # solapamiento: pares de desplazamientos del mismo trozo que empiezan a k cuadros uno del otro
        sol = np.zeros(nt)
        for kk in range(1, L):
            par = ok[:-kk] & ok[kk:] & (trozo[:-L][:-kk] == trozo[:-L][kk:])
            sol += (L - kk) * np.bincount(trozo[:-L][:-kk][par], minlength=nt)
        if len(ok) > L:                           # pares de desplazamientos que empiezan justo a L cuadros
            par = ok[:-L] & ok[L:] & (trozo[:-L][:-L] == trozo[:-L][L:])
            justo = np.bincount(trozo[:-L][:-L][par], minlength=nt)
        else:
            justo = np.zeros(nt)
        with np.errstate(invalid='ignore', divide='ignore'):
            G[:, b] = np.where(cant > 0, cant - 1 - 2 * sol / (cant * L), 0)
            H[:, b] = np.where(cant > 0, cant - 1 + justo / cant, 0)
        poco = G[:, b] < 5                      # trozos con muy pocos desplazamientos: no se usan
        S[poco, b], G[poco, b], H[poco, b] = 0, 0, 0
        if b == 0:
            primero = (cen[~poco[tz]], tz[~poco[tz]])
    return S, G, H, primero


def varianza_por_retardo(x, y, tramo, dt, lags=LAGS_D, trozo_s=TROZO_D_S):
    """Varianza del desplazamiento (alrededor del promedio de cada trozo) para cada retardo, por eje.
    Devuelve dict: tau [L] s, Sx, Gx, Hx, Sy, Gy, Hy [trozos, L] (sumas y grados de libertad por trozo),
    Vx, Vy [L] µm² (juntando los trozos), n_trozos (los que aportan al primer retardo)."""
    trozo = _trozos(tramo, max(int(round(trozo_s / dt)), 2 * max(lags) + 6))
    Sx, Gx, Hx, px = _suma_y_grados(x, trozo, lags)
    Sy, Gy, Hy, py = _suma_y_grados(y, trozo, lags)
    with np.errstate(invalid='ignore', divide='ignore'):
        Vx, Vy = Sx.sum(0) / Gx.sum(0), Sy.sum(0) / Gy.sum(0)
    return dict(tau=np.array(lags) * dt, Sx=Sx, Gx=Gx, Hx=Hx, Sy=Sy, Gy=Gy, Hy=Hy, Vx=Vx, Vy=Vy,
                n_trozos=int((Gx[:, 0] > 0).sum()) if len(Gx) else 0, pasos_x=px, pasos_y=py)


def _recta(tau, S, G, H, pesos=None):
    """Pendiente y ordenada de la recta V(τ) = pendiente·τ + ordenada, con V = Σ S / Σ G (trozos
    pesados por `pesos`). La ordenada va multiplicada por Σ H / Σ G (casi 1): es lo que le toca al
    error de localización después de restar el promedio de cada trozo. Mínimos cuadrados sin pesos."""
    w = np.ones(len(S)) if pesos is None else pesos
    with np.errstate(invalid='ignore', divide='ignore'):
        V, h = (w @ S) / (w @ G), (w @ H) / (w @ G)
    if not (np.all(np.isfinite(V)) and np.all(np.isfinite(h))):
        return np.nan, np.nan
    p = np.linalg.lstsq(np.c_[tau, h], V, rcond=None)[0]
    return p[0], p[1]


def difusion(x, y, tramo, dt, lags=LAGS_D, trozo_s=TROZO_D_S, n_boot=300, rng=0):
    """D de una partícula [µm²/s]: mitad de la pendiente de la varianza del desplazamiento contra el
    retardo, por eje, y el promedio de los dos ejes.
    Devuelve dict: D, err_D, Dx, Dy, err_Dx, err_Dy, ord_x, ord_y [µm²] (ordenada = 2·ruido² menos
    un poco por la exposición), n_trozos, pasos (cuántos desplazamientos de un cuadro entraron)."""
    v = varianza_por_retardo(x, y, tramo, dt, lags, trozo_s)
    nt = len(v['Sx'])
    res = dict(n_trozos=v['n_trozos'], pasos=int(v['Gx'][:, 0].sum() + v['n_trozos']) if nt else 0)
    if nt == 0 or v['n_trozos'] == 0:
        return dict(res, D=np.nan, err_D=np.nan, Dx=np.nan, Dy=np.nan, err_Dx=np.nan, err_Dy=np.nan,
                    ord_x=np.nan, ord_y=np.nan)
    g = np.random.default_rng(rng)
    boot = {}
    for eje in 'xy':
        S, G, H = v['S' + eje], v['G' + eje], v['H' + eje]
        m, o = _recta(v['tau'], S, G, H)
        res['D' + eje], res['ord_' + eje] = m / 2, o
        usados = np.flatnonzero(G[:, 0] > 0)
        W = np.zeros((n_boot, nt))
        W[:, usados] = g.multinomial(len(usados), np.full(len(usados), 1 / len(usados)), size=n_boot)
        boot[eje] = np.array([_recta(v['tau'], S, G, H, w)[0] / 2 for w in W])
        res['err_D' + eje] = np.nanstd(boot[eje], ddof=1)
    res['D'] = 0.5 * (res['Dx'] + res['Dy'])
    res['err_D'] = np.nanstd(0.5 * (boot['x'] + boot['y']), ddof=1)
    return res


def pendiente_loglog(tau, V, desde=3):
    """Pendiente de log V contra log τ (ajuste por mínimos cuadrados desde el retardo `desde`).
    1 = browniano libre; cerca de 0 = no se aleja (pegada o encerrada)."""
    ok = np.isfinite(V) & (V > 0) & (np.arange(len(V)) >= desde)
    if ok.sum() < 4:
        return np.nan
    return float(np.polyfit(np.log(tau[ok]), np.log(V[ok]), 1)[0])


# ---------------------------------------------------------------- una fila por partícula
COLUMNAS = ['grabacion', 'muestra', 'particula', 'cuadros', 'cuadros_medidos', 'dist_laser_um', 'en_atrapada',
            'R_ap_um', 'err_R_ap_um', 'disp_R_ap_um', 'contraste', 'D_um2_s', 'err_D_um2_s', 'Dx_um2_s', 'Dy_um2_s',
            'err_Dx_um2_s', 'err_Dy_um2_s', 'ordenada_um2', 'ruido_nm', 'alfa', 'v_um_s', 'x_px', 'y_px']


def por_particula(g, largo_min=200, n_boot=300):
    """Tabla con una fila por trayectoria con radio medido. Columnas:
      grabacion, muestra, particula
      cuadros            posiciones que se usan (sin colisión, con la platina quieta)
      cuadros_medidos    posiciones medidas en total
      dist_laser_um      menor distancia al láser en toda la trayectoria (vacío sin láser)
      en_atrapada        figura como atrapada en atrapada.csv
      R_ap_um, err_R_ap_um, disp_R_ap_um, contraste     de radios.csv (radio aparente, su error, cuánto
                         cambia entre cuadros, y qué tan oscura se ve)
      D_um2_s, err_D_um2_s, Dx_um2_s, Dy_um2_s, err_Dx_um2_s, err_Dy_um2_s
      ordenada_um2       ordenada de la recta de la varianza, promedio de x e y (2·ruido² menos un poco por la exposición)
      ruido_nm           error de localización por eje que sale de la ordenada (vacío si la ordenada es negativa)
      alfa               pendiente log-log de la varianza contra el retardo (1 = browniano)
      v_um_s             velocidad media de arrastre (recta por tramo, promedio pesado)
      x_px, y_px         posición media
    Las que tienen menos de largo_min cuadros usables no se analizan (D vacío)."""
    filas = []
    r = g['radios']
    for k, p in enumerate(g['ids']):
        x, y = g['X'][:, k], g['Y'][:, k]
        ok = np.isfinite(x)
        fila = dict(grabacion=g['nombre'], muestra=g['muestra'], particula=int(p), cuadros=int(ok.sum()),
                    cuadros_medidos=int(np.isfinite(g['Xc'][:, k]).sum()),
                    dist_laser_um=np.nan, en_atrapada=int(p) in g['atrapadas'],
                    R_ap_um=r.r_borde_px.values[k] * ESCALA_UM, err_R_ap_um=r.err_r_borde_px.values[k] * ESCALA_UM,
                    disp_R_ap_um=r.disp_r_borde_px.values[k] * ESCALA_UM, contraste=r.contraste.values[k],
                    x_px=np.nanmean(g['Xc'][:, k]), y_px=np.nanmean(g['Yc'][:, k]))
        if g['trampa_px'] is not None:
            tx, ty = g['trampa_px']
            fila['dist_laser_um'] = float(np.nanmin(np.hypot(g['Xc'][:, k] - tx, g['Yc'][:, k] - ty))) * ESCALA_UM
        if ok.sum() >= largo_min:
            d = difusion(x, y, g['tramo'], g['dt'], n_boot=n_boot, rng=int(p))
            m = varianza_por_retardo(x, y, g['tramo'], g['dt'], LAGS_MSD, TROZO_MSD_S)
            orden = 0.5 * (d['ord_x'] + d['ord_y'])
            fila.update(D_um2_s=d['D'], err_D_um2_s=d['err_D'], Dx_um2_s=d['Dx'], Dy_um2_s=d['Dy'],
                        err_Dx_um2_s=d['err_Dx'], err_Dy_um2_s=d['err_Dy'],
                        ordenada_um2=orden, ruido_nm=1e3 * np.sqrt(orden / 2) if orden > 0 else np.nan,
                        alfa=pendiente_loglog(m['tau'], m['Vx'] + m['Vy']), v_um_s=_arrastre(x, y, g['tramo'], g['dt']))
        filas.append(fila)
    return pd.DataFrame(filas).reindex(columns=COLUMNAS)


def soltadas(g, espera_s=3.0, largo_min=200, n_boot=300):
    """La partícula que estaba atrapada, una vez que se apaga el láser: queda libre y se puede medir
    su D como a cualquier otra. Usa solo los cuadros desde espera_s después del apagado
    (`apagado_s`) de las trayectorias que figuran como atrapadas en atrapada.csv.
    Devuelve una tabla: grabacion, muestra, particula, cuadros, D_um2_s, err_D_um2_s, alfa
    (vacía si el láser no se apagó o no quedan largo_min cuadros)."""
    t0 = apagado_s(g['nombre'])
    filas = []
    if t0 is not None:
        c0 = int(np.ceil((t0 + espera_s) / g['dt']))
        for k, p in enumerate(g['ids']):
            if int(p) not in g['atrapadas']:
                continue
            x, y = g['X'][:, k].copy(), g['Y'][:, k].copy()
            x[:c0], y[:c0] = np.nan, np.nan
            if np.isfinite(x).sum() < largo_min:
                continue
            d = difusion(x, y, g['tramo'], g['dt'], n_boot=n_boot, rng=int(p))
            m = varianza_por_retardo(x, y, g['tramo'], g['dt'], LAGS_MSD, TROZO_MSD_S)
            filas.append(dict(grabacion=g['nombre'], muestra=g['muestra'], particula=int(p), cuadros=int(np.isfinite(x).sum()),
                              D_um2_s=d['D'], err_D_um2_s=d['err_D'], alfa=pendiente_loglog(m['tau'], m['Vx'] + m['Vy'])))
    return pd.DataFrame(filas, columns=['grabacion', 'muestra', 'particula', 'cuadros', 'D_um2_s', 'err_D_um2_s', 'alfa'])


def _arrastre(x, y, tramo, dt):
    """Velocidad media de arrastre [µm/s]: recta x(t), y(t) en cada tramo, promedio pesado por cuadros."""
    vx, vy, w = [], [], []
    t = np.arange(len(x)) * dt
    for tr in np.unique(tramo[(tramo >= 0) & np.isfinite(x)]):
        s = (tramo == tr) & np.isfinite(x)
        if s.sum() >= 30:
            vx.append(np.polyfit(t[s], x[s], 1)[0])
            vy.append(np.polyfit(t[s], y[s], 1)[0])
            w.append(s.sum())
    if not w:
        return np.nan
    return float(np.hypot(np.average(vx, weights=w), np.average(vy, weights=w)))


# ---------------------------------------------------------------- juntar partículas
def juntar_msd(g, ids, lags=LAGS_MSD, trozo_s=TROZO_MSD_S):
    """Varianza del desplazamiento contra el retardo de cada una de las partículas `ids` de la
    grabación g (x e y sumados: es el MSD en 2D sin el arrastre). Devuelve tau [L], V [L, len(ids)]."""
    col = {int(p): k for k, p in enumerate(g['ids'])}
    V = np.full((len(lags), len(ids)), np.nan)
    for j, p in enumerate(ids):
        k = col[int(p)]
        m = varianza_por_retardo(g['X'][:, k], g['Y'][:, k], g['tramo'], g['dt'], lags, trozo_s)
        V[:, j] = m['Vx'] + m['Vy']
    return np.array(lags) * g['dt'], V


def pasos(g, ids, trozo_s=TROZO_D_S):
    """Pasos de un cuadro de las partículas `ids`, ya sin el arrastre (centrados en cada trozo).
    Devuelve una lista con un dict por partícula: dx, dy [µm] (centrados), los pares de pasos
    seguidos (par_x, par_y: a, b) de cada eje para la correlación, y zt_x, zt_y: los pasos divididos
    por el desvío de SU trozo (solo trozos con 40 pasos o más) con nt_x, nt_y: cuántos pasos tiene
    cada uno de esos trozos."""
    col = {int(p): k for k, p in enumerate(g['ids'])}
    trozo = _trozos(g['tramo'], max(int(round(trozo_s / g['dt'])), 14))
    sal = []
    for p in ids:
        k = col[int(p)]
        uno = {}
        for eje, v in (('x', g['X'][:, k]), ('y', g['Y'][:, k])):
            d = v[1:] - v[:-1]
            ok = np.isfinite(d) & (trozo[1:] == trozo[:-1]) & (trozo[:-1] >= 0)
            tz = trozo[:-1]
            nt = int(trozo.max()) + 1
            cant = np.bincount(tz[ok], minlength=nt)
            media = np.bincount(tz[ok], d[ok], minlength=nt) / np.maximum(cant, 1)
            cen = np.where(ok & (cant[np.clip(tz, 0, None)] >= 10), d - media[np.clip(tz, 0, None)], np.nan)
            uno['d' + eje] = cen[np.isfinite(cen)]
            bien = np.isfinite(cen)
            cb = np.bincount(tz[bien], minlength=nt)
            sd = np.sqrt(np.bincount(tz[bien], cen[bien] ** 2, minlength=nt) / np.maximum(cb - 1, 1))
            largo = bien & (cb[np.clip(tz, 0, None)] >= 40)
            uno['zt_' + eje] = cen[largo] / sd[tz[largo]]
            uno['nt_' + eje] = cb[np.unique(tz[largo])]
            seg = np.isfinite(cen[:-1]) & np.isfinite(cen[1:]) & (tz[:-1] == tz[1:])
            uno['par_' + eje] = (cen[:-1][seg], cen[1:][seg])
        sal.append(uno)
    return sal


def movimiento_comun(g, ids, lag=1, trozo_s=TROZO_D_S):
    """Correlación media entre los desplazamientos SIMULTÁNEOS de dos partículas distintas (retardo
    `lag`, sin el arrastre, cada uno dividido por el desvío de su partícula). Si la imagen entera
    temblara (platina, mesa), todas las partículas se moverían juntas y daría mayor que 0.
    Devuelve (correlación, error, cantidad de pares)."""
    col = {int(p): k for k, p in enumerate(g['ids'])}
    trozo = _trozos(g['tramo'], max(int(round(trozo_s / g['dt'])), 2 * lag + 12))
    nt = int(trozo.max()) + 1 if (trozo >= 0).any() else 0
    tz = trozo[:-lag]
    num, den = 0.0, 0.0
    for eje in ('X', 'Y'):
        Z = []
        for p in ids:
            v = g[eje][:, col[int(p)]]
            d = v[lag:] - v[:-lag]
            ok = np.isfinite(d) & (trozo[lag:] == tz) & (tz >= 0)
            cant = np.bincount(tz[ok], minlength=nt)
            media = np.bincount(tz[ok], d[ok], minlength=nt) / np.maximum(cant, 1)
            cen = np.where(ok, d - media[np.clip(tz, 0, None)], np.nan)
            Z.append(cen / np.nanstd(cen) if np.isfinite(cen).sum() > 20 else np.full(len(cen), np.nan))
        if not Z:
            continue
        Z = np.array(Z).T                                   # [cuadros, partículas]
        m = np.isfinite(Z).sum(1)
        s, s2 = np.nansum(Z, 1), np.nansum(Z ** 2, 1)
        pares = m * (m - 1)
        hay = pares > 0
        num += (s[hay] ** 2 - s2[hay]).sum()
        den += pares[hay].sum()
    if den == 0:
        return np.nan, np.nan, 0
    return float(num / den), float(1 / np.sqrt(den / 2 / max(lag, 1))), int(den / 2)


# ---------------------------------------------------------------- simulación para probar el estimador
def simular_como(g, k, D, ruido_um=0.05, v_um_s=0.2, rng=None):
    """Browniano simulado con D conocido, con los MISMOS cuadros medidos, huecos y tramos que la
    partícula k de la grabación g, más un arrastre constante y un error de localización.
    Devuelve x, y [n] µm (NaN donde la partícula real no tiene dato)."""
    rng = rng or np.random.default_rng()
    n, dt = g['n'], g['dt']
    sal = []
    for eje in range(2):
        camino = np.cumsum(rng.standard_normal(n) * np.sqrt(2 * D * dt))
        camino += v_um_s * np.cos(0.7 + eje) * np.arange(n) * dt + rng.standard_normal(n) * ruido_um
        sal.append(np.where(np.isfinite(g['X'][:, k]), camino, np.nan))
    return sal[0], sal[1]
