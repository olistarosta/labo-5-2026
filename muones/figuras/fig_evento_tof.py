"""Forma de onda de un muón en el montaje de tiempo de vuelo (tres detectores apilados).

Lee cada evento de FIGURAS de la corrida RUN (grupo_7_jueves/medidas/tof/<RUN>_idx*_ondas.txt)
y grafica la tensión de los tres canales en función del tiempo, con el nivel del trigger.

Cada figura tiene la ventana completa del DRS4 a la izquierda y el zoom en XLIM a la derecha,
con el pico de cada pulso marcado (vértice de una parábola en 5 muestras, t_pico).

Uso:  python fig_evento_tof.py   -> guarda evento_tof.{pdf,png} y evento_tof_v2.{pdf,png}
"""
from pathlib import Path
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

AQUI = Path(__file__).resolve().parent
DIR_TOF = AQUI.parent / "grupo_7_jueves" / "medidas" / "tof"
RUN = "tof_d61.2cm_20261029_194928"
THR = 120          # nivel del trigger [mV] (el mismo que el corte de filtrado_muones.ipynb)
XLIM = (80, 130)   # [ns]

CANALES = (1, 2, 3)
NOMBRE = {1: "Módulo 1", 2: "Módulo 2", 3: "Módulo 3"}     # CH1 arriba, CH2 medio, CH3 abajo
COLOR = {1: "#2A9D8F", 2: "#E76F51", 3: "#4C72B0"}
MARCA = {1: "o", 2: "D", 3: "s"}
GRIS = "#264653"


def leer_evento(run, evento):
    """Forma de onda {canal: (t [ns], V [mV])} y fila del CSV de eventos."""
    patron = re.compile(rf"^Event #{evento}\b[^\n]*\n(.*?)(?=^Event #|\Z)", re.M | re.S)
    for ruta in sorted(DIR_TOF.glob(f"{run}_idx*_ondas.txt")):
        m = patron.search(ruta.read_text())
        if m:
            arr = np.array(m.group(1).split(), dtype=float).reshape(-1, 2 * len(CANALES))
            csv = ruta.with_name(ruta.name.replace("_ondas.txt", "_eventos.csv"))
            df = pd.read_csv(csv, comment="#")
            fila = df.loc[df["evento"] == evento].iloc[0]
            return {ch: (arr[:, 2 * k], arr[:, 2 * k + 1]) for k, ch in enumerate(CANALES)}, fila
    raise ValueError(f"no se encontró el evento {evento} en {run}")


def t_pico(t, v, m=2):
    """Tiempo del mínimo de v: vértice de una parábola ajustada a 2m+1 muestras alrededor."""
    ip = int(np.argmin(v[5:-5])) + 5           # sin las 5 muestras de cada punta (artefactos)
    k = slice(ip - m, ip + m + 1)
    a, b, _ = np.polyfit(t[k] - t[ip], v[k], 2)
    return t[ip] - b / (2 * a) if a > 0 else t[ip]


def _estilo():
    plt.style.use("seaborn-v0_8")
    fs = 14
    plt.rcParams.update({"font.family": "serif", "font.size": fs, "axes.labelsize": fs,
                         "xtick.labelsize": fs, "ytick.labelsize": fs, "legend.fontsize": fs,
                         "mathtext.default": "regular", "mathtext.fontset": "cm",
                         "legend.frameon": True})


def _dibujar(ax, ondas, fila, xlim, ms, marca_t="cfd"):
    """marca_t: 'cfd' (t_CFD del CSV), 'pico' (mínimo de cada pulso) o None."""
    for ch in CANALES:
        t, v = ondas[ch]
        k = (t >= xlim[0]) & (t <= xlim[1])
        ax.plot(t[k], v[k], MARCA[ch], ms=ms, color=COLOR[ch], mec=COLOR[ch], zorder=3,
                label=NOMBRE[ch])
        if marca_t == "cfd":
            ax.axvline(fila[f"t_cfd{ch}_ns"], color=COLOR[ch], ls="--", lw=1.4, zorder=2)
        elif marca_t == "pico":
            ax.axvline(t_pico(t, v), color=COLOR[ch], ls="--", lw=1.4, zorder=2)
    ax.axhline(-THR, color=GRIS, ls="--", lw=1.2, zorder=2)
    ax.set_xlim(*xlim)
    ax.set_xlabel("Tiempo $t$ [ns]")


def _leyenda(marca_t="cfd"):
    h = [Line2D([], [], ls="", marker=MARCA[ch], ms=6, color=COLOR[ch]) for ch in CANALES]
    h += [Line2D([], [], color=GRIS, ls="--", lw=1.2),
          Line2D([], [], color="#8a8a8a", ls="--", lw=1.4)]
    l = [NOMBRE[ch] for ch in CANALES] + ["Nivel del trigger",
                                          r"$t_{CFD}$ (50 %)" if marca_t == "cfd" else "Pico del pulso"]
    return h, l


def _guardar(fig, archivo):
    if archivo is not None:
        fig.savefig(f"{archivo}.pdf", bbox_inches="tight")
        fig.savefig(f"{archivo}.png", bbox_inches="tight", dpi=300)


def figura_completa(ondas, fila, ylim, yticks, archivo=None):
    """Ventana completa del DRS4 a la izquierda y zoom en XLIM a la derecha (mismo eje y),
    con el tiempo del pico de cada pulso marcado en el zoom."""
    _estilo()
    fig, (a, b) = plt.subplots(1, 2, figsize=(13, 4.5), sharey=True,
                               gridspec_kw=dict(width_ratios=[1.5, 1]))
    t_fin = max(ondas[ch][0][-1] for ch in CANALES)
    _dibujar(a, ondas, fila, (0, t_fin), ms=2, marca_t=None)
    a.axvspan(*XLIM, color=GRIS, alpha=0.12, lw=0, zorder=1)        # ventana del zoom
    a.set_xticks(np.arange(0, t_fin + 1, 25))
    a.set_ylabel("Tensión [mV]")

    _dibujar(b, ondas, fila, XLIM, ms=3.5, marca_t="pico")
    b.set_xticks(np.arange(XLIM[0], XLIM[1] + 1, 10))
    b.set_ylim(*ylim)
    b.set_yticks(yticks)

    h, l = _leyenda("pico")
    a.legend(h, l, loc="lower right", fontsize=12, facecolor="white", framealpha=0.9)
    fig.tight_layout()
    _guardar(fig, archivo)
    return fig


# (evento, archivo, ylim, yticks, función)
FIGURAS = [
    (11, "evento_tof", (-420, 100), np.arange(-400, 101, 100), figura_completa),
    (8, "evento_tof_v2", (-250, 80), np.arange(-250, 51, 50), figura_completa),
]

if __name__ == "__main__":
    for evento, archivo, ylim, yticks, fn in FIGURAS:
        ondas, fila = leer_evento(RUN, evento)
        print(f"{archivo}: {RUN}, evento {evento} ({fila['timestamp']}): "
              f"A = {-fila['amp1_mV']:.0f}, {-fila['amp2_mV']:.0f}, {-fila['amp3_mV']:.0f} mV, "
              f"Δt31 = {fila['dt31_ns']:.2f} ns, Δt21 = {fila['dt21_ns']:.2f} ns")
        fn(ondas, fila, ylim, yticks, AQUI / archivo)
