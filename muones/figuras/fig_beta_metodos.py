"""Distribución de β = d13 / (c Δt31) evento por evento para dos criterios de tiempo.

Datos: grupo_7_jueves/medidas/tof/*_eventos.csv y *_ondas.txt (montaje de velocidades).
Selección (igual que analisis/filtrado_muones.ipynb, THR_SEL = 120 mV): el pulso baja de
-THR_SEL en los tres canales, ningún canal sube de +THR_SEL (interferencia) y los tres
tienen t_CFD en el CSV.

Tiempo de cada canal, sobre s(t) = -(V - base), con base del CSV:
  50 %    cruce de s = A/2 en el flanco de subida (CFD 50 %, el de drs_tof)
  umbral  cruce de s = THR_SEL en el flanco de subida
Los cruces se buscan caminando hacia atrás desde el pico e interpolando entre muestras.

Resultado de cada criterio: β = d13 / (c ⟨Δt31⟩), con el error de la media de Δt31 y el
de d13 propagados. Δt31 se promedia sobre todos los muones (también Δt31 ≤ 0); el
histograma solo puede mostrar los que tienen Δt31 > 0. Los β >= BETA_MAX van juntos en una
barra rayada al final del eje.

Offsets entre canales: no medidos, no restados.

Uso:  python fig_beta_metodos.py   -> imprime los resultados y guarda beta_metodos.{pdf,png}
      (a) 50 %, (b) umbral, con el eje y compartido
"""
from pathlib import Path
import re

import numpy as np
import pandas as pd

AQUI = Path(__file__).resolve().parent
DIR_TOF = AQUI.parent / "grupo_7_jueves" / "medidas" / "tof"

THR_SEL = 120.0                      # umbral de selección y de tiempo [mV]
BORDE = 5                            # muestras descartadas en cada punta de la ventana
C_M_NS = 0.299792458                 # velocidad de la luz [m/ns]
D12, SD12 = 0.277, 0.001             # [m], entre caras superiores
D23, SD23 = 0.325, 0.001
D13, SD13 = D12 + D23, np.hypot(SD12, SD23)
CANALES = (1, 2, 3)

_RE_EVENTO = re.compile(rb"^Event #(\d+)[^\n]*\n", re.M)


# ----------------------------------------------------------------- datos
def _leer_ondas(ruta):
    datos = ruta.read_bytes()
    ms = list(_RE_EVENTO.finditer(datos))
    ondas = {}
    for k, m in enumerate(ms):
        fin = ms[k + 1].start() if k + 1 < len(ms) else len(datos)
        ondas[int(m.group(1))] = np.array(datos[m.end():fin].split(), dtype=float).reshape(-1, 6)
    return ondas


def _cruce_flanco(t, s, ip, nivel):
    if not (0 < nivel <= s[ip]):
        return np.nan
    j = ip
    while j > 0 and s[j] >= nivel:
        j -= 1
    if s[j] >= nivel:
        return np.nan
    return t[j] + (nivel - s[j]) * (t[j + 1] - t[j]) / (s[j + 1] - s[j])


METODOS = {
    "50 %": lambda t, s, ip: _cruce_flanco(t, s, ip, s[ip] / 2),
    "Umbral": lambda t, s, ip: _cruce_flanco(t, s, ip, THR_SEL),
}


def cargar():
    """Δt31 [ns] de cada muón seleccionado, una columna por criterio."""
    filas = []
    for csv in sorted(DIR_TOF.glob("*_eventos.csv")):
        df = pd.read_csv(csv, comment="#", on_bad_lines="skip")
        ondas = _leer_ondas(csv.with_name(csv.name.replace("_eventos.csv", "_ondas.txt")))
        for _, f in df.iterrows():
            arr = ondas.get(int(f["evento"]))
            if arr is None or any(pd.isna(f[f"t_cfd{ch}_ns"]) for ch in CANALES):
                continue
            s = {ch: -(arr[:, 2 * k + 1] - f[f"base{ch}_mV"]) for k, ch in enumerate(CANALES)}
            v = {ch: s[ch][BORDE:-BORDE] for ch in CANALES}
            if any(v[ch].min() < -THR_SEL for ch in CANALES):        # sube de +THR_SEL: ruido
                continue
            if any(v[ch].max() < THR_SEL for ch in CANALES):         # no baja de -THR_SEL
                continue
            t = {ch: arr[:, 2 * k] for k, ch in enumerate(CANALES)}
            ip = {ch: BORDE + int(np.argmax(v[ch])) for ch in CANALES}
            fila = {"archivo": csv.name, "evento": int(f["evento"])}
            for nombre, fn in METODOS.items():
                fila[nombre] = fn(t[3], s[3], ip[3]) - fn(t[1], s[1], ip[1])
            filas.append(fila)
    return pd.DataFrame(filas)


def resultado(dt):
    """β = d13 / (c ⟨Δt31⟩) y su error (media de Δt31 y d13)."""
    x = dt.dropna()
    m, em = x.mean(), x.std(ddof=1) / np.sqrt(len(x))
    beta = D13 / (C_M_NS * m)
    return dict(N=len(x), dt=m, err_dt=em, beta=beta, err_beta=beta * np.hypot(SD13 / D13, em / m),
                n_neg=int((x <= 0).sum()))


# ----------------------------------------------------------------- figura
COLORES = {"50 %": "#2A9D8F", "Umbral": "#E76F51"}
BETA_MAX = 4                         # los β >= BETA_MAX van juntos en una barra al final
BINS = np.arange(0, BETA_MAX + 0.0001, 0.1)
DESBORDE = (BETA_MAX + 0.2, 0.3)     # barra de los β >= BETA_MAX: (centro, ancho)
YMAX = 11


def figura(dts, archivo=None):
    """Un panel por criterio, eje y compartido: histograma de β por evento (los β >= BETA_MAX en
    una barra rayada al final), con β(⟨Δt31⟩) ± error sombreado y su valor central."""
    import matplotlib
    if archivo is not None:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    plt.style.use("seaborn-v0_8")
    fs = 14
    plt.rcParams.update({"font.family": "serif", "font.size": fs, "axes.labelsize": fs,
                         "xtick.labelsize": fs, "ytick.labelsize": fs, "legend.fontsize": fs,
                         "mathtext.default": "regular", "mathtext.fontset": "cm",
                         "legend.frameon": True})
    GRIS = "#264653"

    fig, axs = plt.subplots(1, len(dts), figsize=(10, 4.5), sharey=True)
    for ax, letra, (nombre, dt) in zip(axs, "ab", dts.items()):
        r = resultado(dt)
        x = dt.dropna()
        beta = D13 / (C_M_NS * x[x > 0])
        ax.hist(beta[beta < BETA_MAX], bins=BINS, color=COLORES[nombre], alpha=0.85, edgecolor="white",
                lw=0.5, zorder=2)
        ax.bar(DESBORDE[0], int((beta >= BETA_MAX).sum()), width=DESBORDE[1], color="white",
               edgecolor=COLORES[nombre], hatch="///", lw=1.2, zorder=2)
        ax.axvspan(r["beta"] - r["err_beta"], r["beta"] + r["err_beta"], color=GRIS,
                   alpha=0.25, lw=0, zorder=3)                  # por delante de las barras
        ax.axvline(r["beta"], color=GRIS, ls=":", lw=2, zorder=4)
        ax.text(0.03, 0.97, f"({letra})", transform=ax.transAxes, ha="left", va="top", fontsize=fs)
        ax.set_xlim(0, DESBORDE[0] + DESBORDE[1])
        ticks = list(np.arange(0, BETA_MAX, 0.5))
        ax.set_xticks(ticks + [DESBORDE[0]],
                      [f"{t:g}" for t in ticks] + [rf"$\geq${BETA_MAX}"])
        ax.set_xlabel(r"Velocidad $\beta = v/c$")
    axs[0].set_ylim(0, YMAX)
    axs[0].set_yticks(np.arange(0, YMAX, 2))                    # cuentas enteras
    axs[0].set_ylabel("Cantidad de muones")
    axs[-1].legend([Patch(color=GRIS, alpha=0.25, lw=0), Line2D([], [], color=GRIS, ls=":", lw=2)],
                   [r"$\beta$ con su incerteza", r"Valor central"],
                   loc="upper right", fontsize=12, facecolor="white", framealpha=0.9)
    fig.tight_layout()
    if archivo is not None:
        fig.savefig(f"{archivo}.pdf", bbox_inches="tight")
        fig.savefig(f"{archivo}.png", bbox_inches="tight", dpi=300)
    return fig


if __name__ == "__main__":
    df = cargar()
    dts = {m: df[m] for m in METODOS}
    print(f"{len(df)} muones (THR_SEL = {THR_SEL:g} mV), d13 = {100 * D13:.1f} ± {100 * SD13:.2f} cm")
    for m, dt in dts.items():
        r = resultado(dt)
        x = dt.dropna()
        fuera = int((D13 / (C_M_NS * x[x > 0]) >= BETA_MAX).sum())
        print(f"  {m:7s} N = {r['N']}  ⟨Δt31⟩ = {r['dt']:.3f} ± {r['err_dt']:.3f} ns  "
              f"β = {r['beta']:.3f} ± {r['err_beta']:.3f}   (Δt31 ≤ 0: {r['n_neg']}, β ≥ {BETA_MAX}: {fuera})")
    figura(dts, AQUI / "beta_metodos")
