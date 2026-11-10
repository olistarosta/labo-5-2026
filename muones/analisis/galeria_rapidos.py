"""Muestra, uno por uno, los decaimientos rapidos (dt < 200 ns) de dos clases de eventos.

En las dos clases el muon dio pulso en u1 y u2 y no en u3 (freno en u2 o en el hierro):
  - electron solo en u3: el pulso tardio esta en un centellador que no vio al muon;
  - electron solo en u2: el pulso tardio esta en el mismo centellador que el pulso de entrada.

Guarda en muones/salida_decaimientos/:
  galeria_rapidos_electron_en_u3.pdf, galeria_rapidos_electron_en_u2.pdf  (todas las trazas, 30 por pagina)
  comparacion_dt_u2_u3.pdf                                                (dt de las dos clases)

Uso (desde muones/analisis):  python galeria_rapidos.py
"""
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import decaimientos as dec

SALIDA = Path(__file__).resolve().parent.parent / "salida_decaimientos"
VENTANA = (60.0, 1400.0)       # dt aceptado [ns]
RAPIDO = 200.0                 # "rapido" = dt menor a esto [ns]
FILTROS = ["traza limpia", "entrada en u1", "hay pulso tardio", "amplitud tardio", "dt en rango",
           "ancho tardio", "sin pulsos extra", "sin pulsos previos"]

ev = dec.cargar_indice()
tab = dec.tabla_eventos(ev)
sel = dec.combinar(dec.armar_filtros(tab, dt_rango=VENTANA), FILTROS)
dt = tab["dt"].to_numpy()
entrada, tardio = tab["patron_ent"].to_numpy(), tab["patron_tar"].to_numpy()
clases = {"electron_en_u3": sel & (entrada == "110") & (tardio == "001"),
          "electron_en_u2": sel & (entrada == "110") & (tardio == "010")}

# ---- galerias: todas las trazas rapidas de cada clase, ordenadas por dt
for nombre, m in clases.items():
    idx = np.flatnonzero(m & (dt < RAPIDO))
    idx = idx[np.argsort(dt[idx])]
    with PdfPages(SALIDA / f"galeria_rapidos_{nombre}.pdf") as pdf:
        for k in range(0, len(idx), 30):
            fig = dec.galeria(ev, idx[k:k + 30], tab, columnas=5, zoom=(120.0, 480.0), suavizar=False, alto=2.2)
            pdf.savefig(fig)
            plt.close(fig)
    print(f"{nombre}: {int(m.sum())} eventos, {len(idx)} con dt < {RAPIDO:.0f} ns ({100 * len(idx) / m.sum():.1f} ± "
          f"{100 * np.sqrt(len(idx)) / m.sum():.1f} %)")

# ---- dt de las dos clases, con el estilo de ESTILO_FIGURAS.md
plt.style.use("seaborn-v0_8")
fs = 14
plt.rcParams.update({"font.family": "serif", "font.size": fs, "axes.labelsize": fs, "xtick.labelsize": fs,
                     "ytick.labelsize": fs, "legend.fontsize": fs, "mathtext.default": "regular",
                     "mathtext.fontset": "cm", "legend.frameon": True})
bordes = np.array([60, 200, 400, 600, 800, 1000, 1200, 1400])
centro, por_100 = 0.5 * (bordes[1:] + bordes[:-1]), np.diff(bordes) / 100
fig, ax = plt.subplots(figsize=(8, 4.5))
for (nombre, m), color, marca, corrimiento, etiqueta in zip(clases.items(), ("#E76F51", "#2A9D8F"), ("D", "o"), (+12, -12),
                                                            ("Electrón solo en el centellador 3", "Electrón solo en el centellador 2")):
    n = np.histogram(dt[m], bordes)[0]
    ax.errorbar(centro + corrimiento, 100 * n / por_100 / m.sum(), yerr=100 * np.sqrt(n) / por_100 / m.sum(), fmt=marca, ms=7,
                capsize=3, elinewidth=1.6, capthick=1.6, mew=1.6, color=color, mfc=color, mec=color, zorder=3, label=etiqueta)
ax.set_xlabel(r"Tiempo entre pulsos $\Delta t$ [ns]")
ax.set_ylabel("Eventos de la clase cada 100 ns [%]")
ax.set_ylim(0, None)
ax.legend(loc="lower left", facecolor="white", framealpha=0.9, fontsize=12)
fig.tight_layout()
fig.savefig(SALIDA / "comparacion_dt_u2_u3.pdf", bbox_inches="tight")
fig.savefig(SALIDA / "comparacion_dt_u2_u3.png", bbox_inches="tight", dpi=300)
