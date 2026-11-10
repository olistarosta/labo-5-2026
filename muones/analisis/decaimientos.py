"""Herramientas para buscar decaimientos de muones en las mediciones del DRS4.

Los archivos `run_<fecha>_<hora>_idx<N>_selected.txt` que escribe `drs_exam` traen,
por evento, 1024 muestras (~1,7 us a 0,6 GS/s) de los tres centelladores. Un muon que
frena y decae deja el pulso de entrada (~210 ns, el que dispara) y, mas tarde, el pulso
del electron. Aca se hace lo pesado, para que el notebook `decaimientos.ipynb` quede
corto:

    indexar()            recorre todos los archivos una sola vez y guarda, por evento,
                         los pulsos encontrados en cada canal (cache en `_indice/`).
    tabla_eventos()      arma las variables fisicas (pulso de entrada, pulso tardio, dt).
    armar_filtros()      los cortes, cada uno por separado, para ir sumandolos de a uno.
    galeria() / explorar()  muchos eventos a la vez, con los pulsos marcados.
    ajustar_vida_media() ajuste de maxima verosimilitud sin binear.
    ajustar_dos_exponenciales()  lo mismo con dos componentes (mu+ / mu-).
    energia_depositada() cuanta energia deja un muon en una barra (Bethe-Bloch).
"""
from __future__ import annotations

import re
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d
from scipy.optimize import minimize

try:
    from numba import njit
except ImportError:  # sin numba anda igual, pero el indexado tarda mucho mas
    def njit(*args, **kwargs):
        return (lambda f: f) if not (args and callable(args[0])) else args[0]

RAIZ_DATOS = Path(__file__).resolve().parent.parent / "data"
CANALES = ("u1", "u2", "u3")
COLORES = ("#4C72B0", "#E76F51", "#2A9D8F")     # paleta de ESTILO_FIGURAS.md
COLOR_DATOS, COLOR_AJUSTE = "#2A9D8F", "#264653"

# --- deteccion de pulsos (se aplica una sola vez, al indexar) -----------------
VERSION_INDICE = 3
N_SUAVE = 6        # muestras (~10 ns): promedia el ruido periodico de ~100 MHz del DRS4
N_FONDO = 96       # muestras (~160 ns): linea de base local, sigue la deriva lenta
UMBRAL_DET = 20.0  # mV por debajo de la linea de base local; los cortes finos van despues
K_PULSOS = 6       # pulsos guardados por canal y evento (los de mayor amplitud)
MARGEN = 15        # muestras (~25 ns) ignoradas en cada borde: ahi el fondo local no vale

_RE_EVENTO = re.compile(rb"Event #(-?\d+)\s+ts=(\S+) (\S+)")
_RE_ARCHIVO = re.compile(r"^run_(\d{8}_\d{6})_idx(\d+)_selected\.txt$")


# ============================================================== lectura =======
def leer_archivo(ruta):
    """Lee todos los eventos de un archivo.

    Devuelve dict con num (N,), ts (N,) datetime64[ms], offset (N,), t (N,S) y v (N,S,3).
    Los eventos cortados (archivo interrumpido) se descartan.
    """
    crudo = Path(ruta).read_bytes()
    inicios = [m.start() for m in re.finditer(rb"Event #", crudo)]
    num, ts, off, datos = [], [], [], []
    for a, b in zip(inicios, inicios[1:] + [len(crudo)]):
        partes = crudo[a:b].split(b"\n", 2)
        m = _RE_EVENTO.match(partes[0])
        if m is None or len(partes) < 3:
            continue
        try:
            fila = np.array(partes[2].split(), dtype=np.float32)
        except ValueError:
            continue
        if fila.size != 4 * 1024:
            continue
        num.append(int(m.group(1)))
        ts.append((m.group(2) + b"T" + m.group(3)).decode())
        off.append(a)
        datos.append(fila.reshape(1024, 4))
    if not datos:
        return dict(num=np.zeros(0, np.int64), ts=np.zeros(0, "datetime64[ms]"),
                    offset=np.zeros(0, np.int64), t=np.zeros((0, 1024), np.float32),
                    v=np.zeros((0, 1024, 3), np.float32))
    datos = np.stack(datos)
    return dict(num=np.array(num, np.int64), ts=np.array(ts, "datetime64[ms]"),
                offset=np.array(off, np.int64), t=datos[:, :, 0], v=datos[:, :, 1:])


def leer_evento(ruta, offset):
    """Lee un solo evento a partir de su posicion en el archivo: (t, v[S,3])."""
    with open(ruta, "rb") as f:
        f.seek(int(offset))
        trozo = f.read(80_000)
    fin = trozo.find(b"Event #", 7)
    if fin > 0:
        trozo = trozo[:fin]
    fila = np.array(trozo.split(b"\n", 2)[2].split(), dtype=np.float32)
    fila = fila[: fila.size // 4 * 4].reshape(-1, 4)
    return fila[:, 0], fila[:, 1:]


# ============================================================== pulsos ========
def senales(v):
    """Devuelve (x, xs, y): sin linea de base, suavizada, y suavizada menos fondo local.

    v tiene el tiempo en el eje -2 y los canales en el -1.
    """
    x = v - np.median(v, axis=-2, keepdims=True)
    xs = uniform_filter1d(x, N_SUAVE, axis=-2, mode="nearest")
    y = xs - uniform_filter1d(x, N_FONDO, axis=-2, mode="nearest")
    return x, xs, y


@njit(cache=True)
def _buscar_pulsos(y, xs, x, t, umbral, K, margen):
    """Pulsos = tramos con y < -umbral. Guarda los K de mayor amplitud, ordenados en t."""
    N, S = y.shape
    tini = np.full((N, K), np.nan, np.float32)   # cruce del umbral
    tmin = np.full((N, K), np.nan, np.float32)   # instante del minimo
    amp = np.full((N, K), np.nan, np.float32)    # profundidad en xs (suavizada) [mV]
    ampc = np.full((N, K), np.nan, np.float32)   # profundidad en x (cruda) [mV]
    ancho = np.full((N, K), np.nan, np.float32)  # tiempo bajo el umbral [ns]
    carga = np.full((N, K), np.nan, np.float32)  # area [mV ns]
    npul = np.zeros(N, np.int16)
    for i in range(N):
        dt = (t[i, S - 1] - t[i, 0]) / (S - 1)
        n = 0
        j = margen
        fin = S - margen
        while j < fin:
            if y[i, j] < -umbral:
                j0 = j
                a = 0.0
                ac = 0.0
                jm = j
                q = 0.0
                while j < fin and y[i, j] < -umbral:
                    if xs[i, j] < a:
                        a = xs[i, j]
                        jm = j
                    if x[i, j] < ac:
                        ac = x[i, j]
                    q += x[i, j]
                    j += 1
                npul[i] += 1
                # lugar entre los K mas grandes
                if n < K:
                    k = n
                    n += 1
                else:
                    k = 0
                    for m in range(1, K):
                        if amp[i, m] < amp[i, k]:
                            k = m
                    if -a <= amp[i, k]:
                        continue
                tini[i, k] = t[i, j0]
                tmin[i, k] = t[i, jm]
                amp[i, k] = -a
                ampc[i, k] = -ac
                ancho[i, k] = (j - j0) * dt
                carga[i, k] = -q * dt
            else:
                j += 1
        # ordenar por tiempo (insercion; K es chico)
        for p in range(1, n):
            q_ = p
            while q_ > 0 and tini[i, q_ - 1] > tini[i, q_]:
                for arr in (tini, tmin, amp, ampc, ancho, carga):
                    tmp = arr[i, q_]
                    arr[i, q_] = arr[i, q_ - 1]
                    arr[i, q_ - 1] = tmp
                q_ -= 1
    return tini, tmin, amp, ampc, ancho, carga, npul


def extraer_pulsos(t, v):
    """Caracteriza los pulsos de N eventos. t (N,S), v (N,S,3)."""
    x, xs, y = senales(v)
    N = v.shape[0]
    sal = dict(
        base=np.median(v, axis=1).astype(np.float32),
        sigma=(1.4826 * np.median(np.abs(x), axis=1)).astype(np.float32),
        vmax=xs.max(axis=1).astype(np.float32),        # mayor excursion positiva (suavizada)
        vmax_cruda=x.max(axis=1).astype(np.float32),
        npul=np.zeros((N, 3), np.int16),
    )
    for nombre in ("tini", "tmin", "amp", "ampc", "ancho", "carga"):
        sal[nombre] = np.full((N, 3, K_PULSOS), np.nan, np.float32)
    for c in range(3):
        r = _buscar_pulsos(np.ascontiguousarray(y[:, :, c]), np.ascontiguousarray(xs[:, :, c]),
                           np.ascontiguousarray(x[:, :, c]), np.ascontiguousarray(t),
                           UMBRAL_DET, K_PULSOS, MARGEN)
        for nombre, arr in zip(("tini", "tmin", "amp", "ampc", "ancho", "carga"), r[:6]):
            sal[nombre][:, c, :] = arr
        sal["npul"][:, c] = r[6]
    return sal


# ============================================================== indice ========
def listar_runs(raiz=RAIZ_DATOS):
    return sorted(p.name for p in Path(raiz).glob("run_*") if p.is_dir())


def listar_archivos(carpeta_run):
    """Archivos de una corrida, ordenados por idx: [(idx, ruta), ...]."""
    arch = []
    for p in Path(carpeta_run).iterdir():
        m = _RE_ARCHIVO.match(p.name)
        if m:
            arch.append((int(m.group(2)), p))
    return sorted(arch)


def _indexar_archivo(arg):
    idx, ruta = arg
    d = leer_archivo(ruta)
    f = extraer_pulsos(d["t"], d["v"])
    f.update(num=d["num"], ts=d["ts"], offset=d["offset"],
             idx=np.full(len(d["num"]), idx, np.int32))
    return f


def indexar(raiz=RAIZ_DATOS, runs=None, procesos=None, rehacer=False, verbose=True):
    """Arma (o reutiliza) el indice de cada corrida en `<raiz>/_indice/<run>.npz`."""
    raiz = Path(raiz)
    (raiz / "_indice").mkdir(exist_ok=True)
    for run in runs or listar_runs(raiz):
        destino = raiz / "_indice" / f"{run}.npz"
        archivos = listar_archivos(raiz / run)
        if destino.exists() and not rehacer:
            with np.load(destino) as z:
                if int(z["version"]) == VERSION_INDICE and int(z["n_archivos"]) == len(archivos):
                    continue
        if verbose:
            print(f"indexando {run}: {len(archivos)} archivos...", flush=True)
        with ProcessPoolExecutor(procesos) as ex:
            partes = list(ex.map(_indexar_archivo, archivos, chunksize=8))
        junto = {k: np.concatenate([p[k] for p in partes]) for k in partes[0]}
        np.savez_compressed(destino, version=VERSION_INDICE, n_archivos=len(archivos), **junto)
        if verbose:
            print(f"   {len(junto['num'])} eventos", flush=True)
    return cargar_indice(raiz, runs)


@dataclass
class Eventos:
    """Indice de todos los eventos (un elemento por evento en cada array)."""
    raiz: Path
    runs: list            # nombres de las corridas
    run: np.ndarray       # indice dentro de `runs`
    datos: dict           # arrays del indice (idx, num, ts, offset, tini, amp, ...)

    def __len__(self):
        return len(self.run)

    def __getattr__(self, nombre):
        try:
            return self.__dict__["datos"][nombre]
        except KeyError:
            raise AttributeError(nombre) from None

    def ruta(self, i):
        run = self.runs[self.run[i]]
        return self.raiz / run / f"{run}_idx{self.idx[i]}_selected.txt"

    def onda(self, i):
        """(t, v[S,3]) del evento i, leido del archivo original."""
        return leer_evento(self.ruta(i), self.offset[i])

    def etiqueta(self, i):
        return f"{self.runs[self.run[i]][4:12]} idx{self.idx[i]} #{self.num[i]}"


def cargar_indice(raiz=RAIZ_DATOS, runs=None):
    raiz = Path(raiz)
    runs = list(runs or listar_runs(raiz))
    partes, cual = [], []
    for k, run in enumerate(runs):
        with np.load(raiz / "_indice" / f"{run}.npz") as z:
            partes.append({n: z[n] for n in z.files if n not in ("version", "n_archivos")})
        cual.append(np.full(len(partes[-1]["num"]), k, np.int16))
    datos = {n: np.concatenate([p[n] for p in partes]) for n in partes[0]}
    return Eventos(raiz=raiz, runs=runs, run=np.concatenate(cual), datos=datos)


# ============================================================== variables =====
def tabla_eventos(ev, ventana_entrada=(180.0, 260.0), dt_min=40.0, umbral_tardio=30.0,
                  coincidencia=15.0):
    """Variables por evento para cortar.

    Pulso de entrada: el mayor de cada canal dentro de `ventana_entrada` [ns].
    Pulso tardio: el primero, en cualquier canal, que empieza mas de `dt_min` ns despues
    de la entrada y supera `umbral_tardio` mV. dt se mide contra la entrada de u1 (el
    canal que dispara).
    """
    N = len(ev)
    ii = np.arange(N)
    tini, amp = ev.tini, ev.amp
    en_ventana = (tini >= ventana_entrada[0]) & (tini <= ventana_entrada[1])
    a_ent = np.where(en_ventana, amp, -np.inf)
    k_ent = a_ent.argmax(axis=2)
    hay_ent = np.isfinite(a_ent.max(axis=2))
    tab = pd.DataFrame({"run": pd.Categorical.from_codes(ev.run, ev.runs), "idx": ev.idx,
                        "num": ev.num, "ts": ev.ts})
    for c, nom in enumerate(CANALES):
        tab[f"a_ent_{nom}"] = np.where(hay_ent[:, c], amp[ii, c, k_ent[:, c]], 0.0)
        tab[f"t_ent_{nom}"] = np.where(hay_ent[:, c], tini[ii, c, k_ent[:, c]], np.nan)
    t0 = tab["t_ent_u1"].to_numpy()
    t0 = np.where(np.isnan(t0), np.nanmedian(t0), t0)      # sin entrada en u1: referencia tipica
    tab["t0"] = t0
    tab["n_ent"] = hay_ent.sum(axis=1)

    # pulsos tardios (despues de la entrada) y tempranos (antes), por canal
    tarde = (tini > t0[:, None, None] + dt_min) & (amp >= umbral_tardio)
    antes = (tini < ventana_entrada[0] - 20.0) & (amp >= umbral_tardio)
    tab["n_antes"] = antes.sum(axis=(1, 2))
    a_antes = np.where(antes, amp, 0.0).reshape(N, -1)
    tab["a_antes"] = a_antes.max(axis=1)                     # el mayor pulso previo a la entrada
    tab["ancho_antes"] = ev.ancho.reshape(N, -1)[ii, a_antes.argmax(axis=1)]
    t_tarde = np.where(tarde, tini, np.inf)
    k1 = t_tarde.argmin(axis=2)                              # primer tardio de cada canal
    t1 = t_tarde.min(axis=2)
    c1 = t1.argmin(axis=1)                                   # canal con el primer tardio
    t_tar = t1[ii, c1]
    hay = np.isfinite(t_tar)
    kk = k1[ii, c1]
    tab["hay_tardio"] = hay
    tab["canal_tardio"] = np.where(hay, c1 + 1, 0)
    tab["t_tardio"] = np.where(hay, t_tar, np.nan)
    tab["dt"] = tab["t_tardio"] - t0
    for nombre, arr in (("a_tardio", ev.amp), ("ac_tardio", ev.ampc), ("ancho_tardio", ev.ancho),
                        ("q_tardio", ev.carga)):
        tab[nombre] = np.where(hay, arr[ii, c1, kk], np.nan)
    # el pulso tardio suele verse en mas de un centellador: mismos tiempos, +-coincidencia
    junto = tarde & (np.abs(tini - t_tar[:, None, None]) <= coincidencia)
    tab["canales_tardio"] = junto.any(axis=2).sum(axis=1)
    for c, nom in enumerate(CANALES):
        tab[f"a_tar_{nom}"] = np.where(junto[:, c], amp[:, c], 0.0).max(axis=1)
    # topologia: en que centelladores hubo entrada / pulso tardio (ej. "110" = u1 y u2)
    cod = np.array([np.binary_repr(k, 3) for k in range(8)])
    tab["patron_ent"] = cod[hay_ent @ np.array([4, 2, 1])]
    tab["patron_tar"] = cod[junto.any(axis=2) @ np.array([4, 2, 1])]
    ultimo = np.where(hay_ent.any(axis=1), 2 - hay_ent[:, ::-1].argmax(axis=1), -1)
    tab["mismo_canal"] = (junto.any(axis=2) & hay_ent).any(axis=1)   # entrada y tardio en un mismo canal
    tab["frena_en"] = ultimo + 1                             # ultimo centellador con entrada
    tab["tardio_en_freno"] = junto.any(axis=2)[ii, np.maximum(ultimo, 0)] & (ultimo >= 0)
    # otros pulsos tardios en otro instante (tercer pulso, rafagas, pospulsos)
    tab["n_extra"] = (tarde & ~junto).sum(axis=(1, 2))
    tab["a_extra"] = np.where(tarde & ~junto, amp, 0.0).max(axis=(1, 2))
    # calidad de la traza
    tab["vmax"] = ev.vmax.max(axis=1)
    tab["sigma"] = ev.sigma.max(axis=1)
    tab["npul"] = ev.npul.sum(axis=1)
    return tab


def armar_filtros(tab, *, vmax_max=45.0, sigma_max=32.0, a_ent_min=60.0, a_tardio_min=60.0,
                  dt_rango=(100.0, 1400.0), ancho_rango=(10.0, 30.0), a_extra_min=40.0):
    """Devuelve un OrderedDict nombre -> mascara booleana. Cada corte es independiente:
    la seleccion final es el AND de los que se quieran aplicar (ver `combinar`).

    vmax_max, sigma_max  traza limpia: sin excursiones positivas (rafagas de ruido) ni
                         ruido excesivo [mV]
    a_ent_min            amplitud minima del pulso de entrada en u1 [mV]
    a_tardio_min         amplitud minima del pulso tardio [mV]
    dt_rango             ventana de dt aceptada [ns]
    ancho_rango          tiempo bajo umbral del pulso tardio [ns]
    a_extra_min          amplitud a partir de la cual un tercer pulso invalida el evento [mV]
    """
    f = OrderedDict()
    f["traza limpia"] = ((tab["vmax"] < vmax_max) & (tab["sigma"] < sigma_max)).to_numpy()
    f["entrada en u1"] = (tab["a_ent_u1"] >= a_ent_min).to_numpy()
    f["hay pulso tardio"] = tab["hay_tardio"].to_numpy()
    f["amplitud tardio"] = (tab["a_tardio"] >= a_tardio_min).to_numpy()
    f["dt en rango"] = tab["dt"].between(*dt_rango).to_numpy()
    f["ancho tardio"] = tab["ancho_tardio"].between(*ancho_rango).to_numpy()
    f["sin pulsos extra"] = (tab["a_extra"] < a_extra_min).to_numpy()
    f["sin pulsos previos"] = (tab["a_antes"] < a_extra_min).to_numpy()
    f["dos pulsos en un canal"] = tab["mismo_canal"].to_numpy()
    f["tardio donde freno"] = tab["tardio_en_freno"].to_numpy()
    return f


def flujo_de_cortes(filtros):
    """Tabla: cuantos eventos pasan cada corte solo, y acumulando en orden."""
    filas, acum = [], None
    for nombre, m in filtros.items():
        acum = m.copy() if acum is None else acum & m
        filas.append({"corte": nombre, "pasa solo": int(m.sum()), "acumulado": int(acum.sum())})
    df = pd.DataFrame(filas).set_index("corte")
    df["% del anterior"] = (100 * df["acumulado"] / df["acumulado"].shift(1, fill_value=len(m))).round(1)
    return df


def combinar(filtros, activos=None):
    """AND de los filtros nombrados en `activos` (todos si es None)."""
    nombres = list(filtros) if activos is None else list(activos)
    m = np.ones(len(next(iter(filtros.values()))), bool)
    for n in nombres:
        m &= filtros[n]
    return m


def rechazados(filtros, nombre, activos=None):
    """Eventos que pasan todos los filtros activos menos `nombre`, y que `nombre` rechaza.
    Es lo que hay que mirar para saber si ese corte esta tirando decaimientos."""
    activos = list(filtros) if activos is None else list(activos)
    return combinar(filtros, [n for n in activos if n != nombre]) & ~filtros[nombre]


def fondo_accidental(ev, tab, base, a_tardio_min, ancho_rango, dt_rango, ventana=(25.0, 160.0),
                     solo_freno=False):
    """Cota al fondo de coincidencias casuales, a partir de los pulsos ANTERIORES al disparo.

    Un pulso casual cae con la misma probabilidad antes o despues de la entrada. Se cuentan
    los eventos de `base` (mascara) con un pulso en `ventana` [ns] que cumpla los mismos
    cortes de amplitud y ancho que el pulso tardio, y se escala al ancho de `dt_rango`.
    Es una cota superior: parte de esos pulsos previos son decaimientos "al reves" (el
    muon entro sin disparar y el que dispara es el electron).
    """
    previo = ((ev.tini >= ventana[0]) & (ev.tini < ventana[1]) & (ev.amp >= a_tardio_min)
              & (ev.ancho >= ancho_rango[0]) & (ev.ancho <= ancho_rango[1]))
    if solo_freno:
        c = np.maximum(tab["frena_en"].to_numpy() - 1, 0)
        previo = previo[np.arange(len(ev)), c].any(axis=1)
    else:
        previo = previo.any(axis=(1, 2))
    n = int((previo & base).sum())
    escala = (dt_rango[1] - dt_rango[0]) / (ventana[1] - ventana[0])
    return dict(n_previos=n, esperado=n * escala, error=np.sqrt(n) * escala,
                mascara=previo & base)


# ============================================================== galeria =======
def galeria(ev, indices, tab=None, columnas=6, zoom=None, suavizar=True, alto=1.9,
            numerar=True):
    """Dibuja muchos eventos juntos. `indices` son posiciones en `ev` (hasta ~60).

    zoom: None (traza completa), 'tardio' (alrededor del pulso tardio; necesita `tab`)
          o (tmin, tmax) en ns.
    Lineas: gris punteada = entrada (t0); negra = pulso tardio elegido.
    Es un grafico de diagnostico: usa la paleta de ESTILO_FIGURAS.md pero con letra chica
    y una etiqueta por panel (corrida, archivo, evento), que hace falta para identificarlos.
    """
    import matplotlib.pyplot as plt

    indices = np.atleast_1d(indices)
    n = len(indices)
    if n == 0:
        fig, ax = plt.subplots(figsize=(6, 1.2))
        ax.text(0.5, 0.5, "ningun evento en esta seleccion", ha="center", va="center")
        ax.axis("off")
        return fig
    filas = int(np.ceil(n / columnas))
    fig, axs = plt.subplots(filas, columnas, figsize=(3.1 * columnas, alto * filas),
                            squeeze=False, sharex=zoom != "tardio")
    for ax in axs.flat[n:]:
        ax.axis("off")
    for k, (ax, i) in enumerate(zip(axs.flat, indices)):
        t, v = ev.onda(i)
        x, xs, _ = senales(v)
        for c in range(3):
            ax.plot(t, xs[:, c] if suavizar else x[:, c], color=COLORES[c], lw=0.7)
        txt = ev.etiqueta(i)
        if tab is not None:
            fila = tab.iloc[i]
            ax.axvline(fila["t0"], color="0.5", lw=0.6, ls=":")
            if fila["hay_tardio"]:
                ax.axvline(fila["t_tardio"], color="k", lw=0.6, alpha=0.6)
                txt += f"\ndt={fila['dt']:.0f} ns  A={fila['a_tardio']:.0f} mV  u{fila['canal_tardio']}"
                if zoom == "tardio":
                    ax.set_xlim(fila["t_tardio"] - 60, fila["t_tardio"] + 100)
        if isinstance(zoom, tuple):
            ax.set_xlim(*zoom)
        elif zoom is None:
            ax.set_xlim(t[0], t[-1])
        ax.set_title(txt, fontsize=7, pad=2)
        ax.tick_params(labelsize=6, pad=1)
        ax.axhline(0, color="0.8", lw=0.4)
        if numerar:
            ax.text(0.02, 0.05, str(k + 1), transform=ax.transAxes, fontsize=9, weight="bold",
                    color="0.25")
    for c in range(3):
        axs[0, 0].plot([], [], color=COLORES[c], label=CANALES[c])
    axs[0, 0].legend(fontsize=6, loc="lower right", ncol=3, handlelength=1, frameon=False)
    fig.tight_layout(pad=0.4)
    return fig


def ver_evento(ev, i, tab=None):
    """Un evento en grande: traza cruda (clara), suavizada (oscura) y los pulsos detectados."""
    import matplotlib.pyplot as plt

    t, v = ev.onda(i)
    x, xs, _ = senales(v)
    fig, axs = plt.subplots(3, 1, figsize=(11, 6.5), sharex=True)
    for c, ax in enumerate(axs):
        ax.plot(t, x[:, c], color=COLORES[c], lw=0.6, alpha=0.35)
        ax.plot(t, xs[:, c], color=COLORES[c], lw=1.1, label=CANALES[c])
        for tp, ap in zip(ev.tini[i, c], ev.amp[i, c]):
            if np.isfinite(tp):
                ax.annotate(f"{ap:.0f} mV", (tp, -ap), xytext=(4, -2), textcoords="offset points",
                            fontsize=8)
                ax.plot(tp, -ap, "v", color="k", ms=4)
        ax.set_ylabel(f"{CANALES[c]} [mV]")
        if tab is not None:
            ax.axvline(tab["t0"].iloc[i], color="0.5", lw=0.7, ls=":")
            if tab["hay_tardio"].iloc[i]:
                ax.axvline(tab["t_tardio"].iloc[i], color="k", lw=0.7, alpha=0.5)
    axs[-1].set_xlabel("Tiempo [ns]")
    fig.tight_layout()
    return fig


def cargar_vetos(ruta):
    ruta = Path(ruta)
    if not ruta.exists():
        return pd.DataFrame(columns=["run", "idx", "num"])
    return pd.read_csv(ruta)


def mascara_vetos(tab, vetos):
    """True para los eventos que NO estan vetados a mano."""
    if len(vetos) == 0:
        return np.ones(len(tab), bool)
    clave = pd.MultiIndex.from_frame(tab[["run", "idx", "num"]].astype({"run": str}))
    malos = pd.MultiIndex.from_frame(vetos[["run", "idx", "num"]].astype({"run": str}))
    return ~clave.isin(malos)


def explorar(ev, tab, filtros, activos=None, por_pagina=24, ruta_vetos=None, semilla=0):
    """Visor interactivo: tildar filtros de a uno y ver que eventos quedan o se van.

    - Casillas: que filtros estan aplicados (al empezar, los de `activos`).
    - "Mostrar": los que pasan todos los filtros tildados, o los que rechaza un filtro
      en particular habiendo pasado los demas (para ver si ese corte tira decaimientos).
    - "Corrida": todas o una sola.
    - "Vetar": anota a mano eventos que no son decaimientos (numero de panel) en `ruta_vetos`.
    """
    import ipywidgets as w
    import matplotlib.pyplot as plt
    from IPython.display import display

    nombres = list(filtros)
    activos = nombres if activos is None else activos
    casillas = [w.Checkbox(value=n in activos, description=n, indent=False,
                           layout=w.Layout(width="190px")) for n in nombres]
    mostrar = w.Dropdown(options=["pasan los filtros tildados"] + [f"rechazados por: {n}" for n in nombres],
                         description="Mostrar", layout=w.Layout(width="360px"))
    corrida = w.Dropdown(options=[("todas", -1)] + [(r[4:], k) for k, r in enumerate(ev.runs)],
                         description="Corrida", layout=w.Layout(width="260px"))
    orden = w.Dropdown(options=["al azar", "dt creciente", "amplitud tardio creciente", "cronologico"],
                       description="Orden", layout=w.Layout(width="290px"))
    zoom = w.Dropdown(options=[("traza completa", None), ("zoom al pulso tardio", "tardio"),
                               ("zoom a la entrada", (150.0, 400.0))],
                      description="Vista", layout=w.Layout(width="290px"))
    suave = w.Checkbox(value=True, description="suavizado", indent=False, layout=w.Layout(width="110px"))
    n_pag = w.Dropdown(options=[12, 24, 48], value=por_pagina, description="Por pag.",
                       layout=w.Layout(width="170px"))
    pagina = w.BoundedIntText(value=1, min=1, max=1, description="Pagina", layout=w.Layout(width="170px"))
    ant, sig = w.Button(description="< anterior"), w.Button(description="siguiente >")
    vetar_txt = w.Text(placeholder="paneles a vetar, ej: 3, 7, 12", layout=w.Layout(width="260px"))
    vetar_btn = w.Button(description="Vetar", button_style="warning", layout=w.Layout(width="80px"))
    info = w.HTML()
    salida = w.Output()
    estado = {"sel": np.zeros(0, int), "visibles": np.zeros(0, int)}

    def seleccion():
        activos = [n for n, c in zip(nombres, casillas) if c.value]
        if mostrar.value.startswith("rechazados por: "):
            m = rechazados(filtros, mostrar.value.split(": ", 1)[1], activos)
        else:
            m = combinar(filtros, activos)
        if ruta_vetos is not None:
            m &= mascara_vetos(tab, cargar_vetos(ruta_vetos))
        if corrida.value >= 0:
            m &= ev.run == corrida.value
        sel = np.flatnonzero(m)
        if orden.value == "al azar":
            sel = np.random.default_rng(semilla).permutation(sel)
        elif orden.value == "dt creciente":
            sel = sel[np.argsort(tab["dt"].to_numpy()[sel], kind="stable")]
        elif orden.value == "amplitud tardio creciente":
            sel = sel[np.argsort(tab["a_tardio"].to_numpy()[sel], kind="stable")]
        return sel

    def redibujar(recalcular=True):
        if recalcular:
            estado["sel"] = seleccion()
        sel = estado["sel"]
        paginas = max(1, int(np.ceil(len(sel) / n_pag.value)))
        pagina.max = paginas
        pagina.value = min(pagina.value, paginas)
        vis = sel[(pagina.value - 1) * n_pag.value: pagina.value * n_pag.value]
        estado["visibles"] = vis
        info.value = (f"<b>{len(sel):,}</b> eventos de {len(tab):,} "
                      f"({100 * len(sel) / max(1, len(tab)):.3f} %) &nbsp;|&nbsp; pagina {pagina.value} de {paginas}")
        with salida:
            salida.clear_output(wait=True)
            fig = galeria(ev, vis, tab, zoom=zoom.value, suavizar=suave.value)
            display(fig)
            plt.close(fig)

    silencio = {"on": False}

    def al_cambiar(_):
        if silencio["on"]:
            return
        silencio["on"] = True
        pagina.value = 1
        silencio["on"] = False
        redibujar()

    def al_paginar(_):
        if not silencio["on"]:
            redibujar(recalcular=False)

    for c in casillas + [mostrar, corrida, orden]:
        c.observe(al_cambiar, names="value")
    for c in (zoom, suave, n_pag, pagina):
        c.observe(al_paginar, names="value")
    ant.on_click(lambda _: setattr(pagina, "value", max(1, pagina.value - 1)))
    sig.on_click(lambda _: setattr(pagina, "value", min(pagina.max, pagina.value + 1)))

    def al_vetar(_):
        if ruta_vetos is None:
            return
        try:
            cuales = [int(s) for s in re.split(r"[,\s]+", vetar_txt.value.strip()) if s]
        except ValueError:
            return
        elegidos = [estado["visibles"][k - 1] for k in cuales if 1 <= k <= len(estado["visibles"])]
        if not elegidos:
            return
        nuevos = tab.iloc[elegidos][["run", "idx", "num"]].astype({"run": str})
        todos = pd.concat([cargar_vetos(ruta_vetos), nuevos]).drop_duplicates()
        todos.to_csv(ruta_vetos, index=False)
        vetar_txt.value = ""
        redibujar()

    vetar_btn.on_click(al_vetar)
    filas = [w.HBox(casillas, layout=w.Layout(flex_flow="row wrap")),
             w.HBox([mostrar, corrida, orden, zoom, suave]),
             w.HBox([ant, sig, pagina, n_pag] + ([vetar_txt, vetar_btn] if ruta_vetos else [])),
             info, salida]
    panel = w.VBox(filas)
    display(panel)
    redibujar()


# ============================================================== ajuste ========
def _nll(tau, f, dt, a, b):
    if tau <= 0 or not 0 <= f <= 1:
        return 1e300
    norma = tau * (np.exp(-a / tau) - np.exp(-b / tau))
    return -np.log(f * np.exp(-dt / tau) / norma + (1 - f) / (b - a)).sum()


def ajustar_vida_media(dt, rango, frac_fondo=0.0):
    """Maxima verosimilitud sin binear de  f*exp(-dt/tau) + (1-f)*plano  en `rango` [ns].

    frac_fondo: fraccion de fondo plano (1-f). Un numero la deja fija (0 = exponencial
    pura); None la deja libre. Devuelve dict con tau, sus errores (donde -lnL sube 1/2,
    perfilando f si esta libre), f y n.
    """
    a, b = rango
    dt = np.asarray(dt, float)
    dt = dt[(dt >= a) & (dt <= b)]
    libre = frac_fondo is None
    opc = dict(xatol=1e-4, fatol=1e-7, maxiter=4000)

    def perfil(tau_, f0=0.9):
        """-lnL minimizada en f (si esta libre) para un tau dado."""
        if not libre:
            return _nll(tau_, 1.0 - frac_fondo, dt, a, b), 1.0 - frac_fondo
        r = minimize(lambda q: _nll(tau_, q[0], dt, a, b), [f0], method="Nelder-Mead", options=opc)
        return r.fun, r.x[0]

    mejor = None
    for tau0 in (1000.0, 2200.0, 4000.0):
        r = minimize(lambda q: perfil(q[0])[0], [tau0], method="Nelder-Mead", options=opc)
        if mejor is None or r.fun < mejor.fun:
            mejor = r
    tau = mejor.x[0]
    f = perfil(tau)[1]

    def cruce(signo):
        paso = 0.02 * tau
        x0, x1 = tau, tau + signo * paso
        for _ in range(200):
            if x1 <= 1.0:
                return np.nan
            if perfil(x1, f)[0] - mejor.fun >= 0.5:
                break
            x0, x1 = x1, x1 + signo * paso
            paso *= 1.2
        else:
            return np.inf
        for _ in range(30):
            xm = 0.5 * (x0 + x1)
            if perfil(xm, f)[0] - mejor.fun >= 0.5:
                x1 = xm
            else:
                x0 = xm
        return abs(0.5 * (x0 + x1) - tau)

    return dict(tau=tau, err_menos=cruce(-1), err_mas=cruce(+1), f=f, n=len(dt), rango=(a, b),
                fondo_libre=libre, nll=mejor.fun)


TAU_MU, TAU_MU_C, TAU_MU_FE = 2197.0, 2026.0, 206.0   # ns: muon libre, mu- en carbono, mu- en hierro


def _exp_norm(dt, tau, a, b):
    return np.exp(-dt / tau) / (tau * (np.exp(-a / tau) - np.exp(-b / tau)))


def ajustar_dos_exponenciales(dt, rango, tau_lenta=TAU_MU, tau_rapida=TAU_MU_FE, puntos=101):
    """Maxima verosimilitud sin binear de  f*exp(-dt/tau_rapida) + (1-f)*exp(-dt/tau_lenta),
    cada exponencial normalizada en `rango` [ns]. f es la fraccion de los eventos de la
    ventana que esta en la componente rapida. Una tau en None queda libre.

    Devuelve f con su intervalo de 68 % (donde el perfil de -lnL sube 1/2), las taus, y
    `sigmas` = sqrt(2*[lnL(mejor) - lnL(f=0)]): cuanto mejora el ajuste al agregar la
    componente rapida. `perfil` = (f, -lnL - minimo) para graficar.
    """
    a, b = rango
    dt = np.asarray(dt, float)
    dt = dt[(dt >= a) & (dt <= b)]
    libres = [n for n, v in (("tau_lenta", tau_lenta), ("tau_rapida", tau_rapida)) if v is None]
    opc = dict(xatol=1e-2, fatol=1e-8, maxiter=4000)

    def taus(q):
        p = dict(tau_lenta=tau_lenta, tau_rapida=tau_rapida)
        p.update(zip(libres, q))
        return p["tau_lenta"], p["tau_rapida"]

    def nll(f, q):
        tl, tr = taus(q)
        if not 0 <= f <= 1 or (libres and not (20 < tr < 0.8 * tl and tl < 2e4)):
            return 1e300
        return -np.log(f * _exp_norm(dt, tr, a, b) + (1 - f) * _exp_norm(dt, tl, a, b)).sum()

    def perfil(f, semillas):
        if not libres:
            return nll(f, ()), ()
        mejor = None
        for q0 in semillas:
            r = minimize(lambda q: nll(f, q), q0, method="Nelder-Mead", options=opc)
            if mejor is None or r.fun < mejor.fun:
                mejor = r
        return mejor.fun, tuple(mejor.x)

    inicios = {"tau_lenta": (1500.0, 2200.0, 3500.0), "tau_rapida": (80.0, 200.0, 600.0)}
    semillas = [()]
    for n in libres:
        semillas = [q + (x,) for q in semillas for x in inicios[n]]
    fs = np.linspace(0, 1, puntos)
    valores, sol = [], []
    for f in fs:                      # barrido en f, arrancando cada punto desde el anterior
        v, q = perfil(f, semillas + sol[-1:])
        valores.append(v)
        sol.append(q)
    valores = np.array(valores)
    k = int(valores.argmin())
    # afinar f alrededor del minimo del barrido
    lo, hi = fs[max(k - 1, 0)], fs[min(k + 1, puntos - 1)]
    finos = np.linspace(lo, hi, 41)
    vf = np.array([perfil(f, [sol[k]])[0] for f in finos])
    f_mejor, minimo = finos[vf.argmin()], vf.min()
    q_mejor = perfil(f_mejor, [sol[k]])[1]
    todo_f, todo_v = np.concatenate([fs, finos]), np.concatenate([valores, vf]) - minimo
    orden = np.argsort(todo_f)
    todo_f, todo_v = todo_f[orden], todo_v[orden]
    dentro = todo_f[todo_v <= 0.5]
    tl, tr = taus(q_mejor)
    return dict(n=len(dt), f=f_mejor, f_min=dentro.min(), f_max=dentro.max(), tau_lenta=tl, tau_rapida=tr,
                sigmas=np.sqrt(max(0.0, 2 * (valores[0] - minimo))), nll=minimo, rango=(a, b),
                perfil=(todo_f, todo_v))


def graficar_ajuste(dt, res, bins=26, ax=None):
    """Histograma de dt con la curva ajustada, con el estilo de ESTILO_FIGURAS.md.
    Los valores del ajuste no van en la figura: estan en `res`."""
    import matplotlib.pyplot as plt

    a, b = res["rango"]
    dt = np.asarray(dt, float)
    dt = dt[(dt >= a) & (dt <= b)]
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 4.5))
    cuentas, bordes = np.histogram(dt, bins=bins, range=(a, b))
    centro, anchob = 0.5 * (bordes[1:] + bordes[:-1]), bordes[1] - bordes[0]
    tau, f = res["tau"], res["f"]
    xx = np.linspace(a, b, 400)
    norma = tau * (np.exp(-a / tau) - np.exp(-b / tau))
    ax.plot(xx / 1000, len(dt) * anchob * (f * np.exp(-xx / tau) / norma + (1 - f) / (b - a)),
            color=COLOR_AJUSTE, lw=2, zorder=2,
            label="Ajuste exponencial" + (" + fondo" if f < 1 else ""))
    ax.errorbar(centro / 1000, cuentas, yerr=np.sqrt(cuentas), fmt="o", ms=7, capsize=3,
                elinewidth=1.6, capthick=1.6, mew=1.6, color=COLOR_DATOS, mfc=COLOR_DATOS,
                mec=COLOR_DATOS, zorder=3, label="Decaimientos")
    ax.set_xlabel(r"Tiempo entre pulsos $\Delta t$ [$\mu$s]")
    ax.set_ylabel(f"Eventos / {anchob:.0f} ns")
    ax.set_ylim(0, None)
    ax.legend(loc="upper right", facecolor="white", framealpha=0.9)
    ax.figure.tight_layout()
    return ax.figure



# ============================================================== Bethe-Bloch ===
M_MU, M_E = 105.658, 0.510999            # MeV
# centellador plastico (base de viniltolueno): Z/A, I [MeV], densidad [g/cm3] y
# parametros de Sternheimer del efecto de densidad (tablas del PDG)
PLASTICO = dict(ZA=0.54141, I=64.7e-6, rho=1.032, a=0.1610, k=3.2393, x0=0.1464, x1=2.4855, C=3.1997)
ESPESOR_BARRA = 4.5                      # cm


def bethe_bloch(T, mat=PLASTICO):
    """Perdida media de energia de un muon de energia cinetica T [MeV], en MeV cm2/g."""
    T = np.asarray(T, float)
    gamma = 1 + T / M_MU
    bg = np.sqrt(gamma ** 2 - 1)
    beta2 = 1 - 1 / gamma ** 2
    tmax = 2 * M_E * bg ** 2 / (1 + 2 * gamma * M_E / M_MU + (M_E / M_MU) ** 2)
    x = np.log10(bg)
    delta = np.where(x >= mat["x0"], 2 * np.log(10) * x - mat["C"]
                     + mat["a"] * np.clip(mat["x1"] - x, 0, None) ** mat["k"], 0.0)
    return 0.307075 * mat["ZA"] / beta2 * (0.5 * np.log(2 * M_E * bg ** 2 * tmax / mat["I"] ** 2)
                                           - beta2 - delta / 2)


_T_TABLA = np.geomspace(0.5, 2e5, 4000)
_inv = 1 / bethe_bloch(_T_TABLA)
_R_TABLA = np.concatenate([[0], np.cumsum(0.5 * (_inv[1:] + _inv[:-1]) * np.diff(_T_TABLA))]) + 0.25 * _inv[0]


def alcance(T):
    """Cuanto plastico atraviesa un muon de energia cinetica T [MeV] antes de frenar [g/cm2]."""
    return np.interp(T, _T_TABLA, _R_TABLA)


def energia_para_alcance(R):
    """Inversa de `alcance`: energia cinetica [MeV] de un muon que frena tras R g/cm2 de plastico."""
    return np.interp(R, _R_TABLA, _T_TABLA)


def energia_depositada(T, espesor=ESPESOR_BARRA):
    """Energia media [MeV] que deja un muon de energia cinetica T al cruzar `espesor` cm de
    centellador (toda su energia si frena adentro)."""
    d = espesor * PLASTICO["rho"]
    R = alcance(T)
    return np.where(R > d, T - np.interp(np.clip(R - d, 0, None), _R_TABLA, _T_TABLA), T)


def energia_desde_deposito(dE, espesor=ESPESOR_BARRA):
    """Inversa de `energia_depositada` en la rama de baja energia (entre la energia minima
    para cruzar la barra y el minimo de ionizacion): energia cinetica [MeV] de un muon
    lento que deja dE [MeV]."""
    lenta = (alcance(_T_TABLA) > espesor * PLASTICO["rho"]) & (_T_TABLA < 300)
    T = _T_TABLA[lenta]
    return np.interp(-np.asarray(dE, float), -energia_depositada(T, espesor), T)
