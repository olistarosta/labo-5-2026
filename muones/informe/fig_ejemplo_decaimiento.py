"""Figura de ejemplo para el informe: las tres trazas de un decaimiento.

Uso (desde muones/informe):  python fig_ejemplo_decaimiento.py
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

MUONES = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(MUONES / "analisis"))
import decaimientos as dec

plt.style.use("seaborn-v0_8")
fs = 12
plt.rcParams.update({"font.family": "serif", "font.size": fs, "axes.labelsize": fs, "xtick.labelsize": fs,
                     "ytick.labelsize": fs, "mathtext.default": "regular", "mathtext.fontset": "cm"})

ev = dec.cargar_indice()
tab = dec.tabla_eventos(ev)
# evento del 11/9, archivo idx3, disparo 1247: entrada en los tres centelladores y electron en u3
i = int(np.flatnonzero((tab["run"] == "run_20260911_135116") & (tab["idx"] == 3) & (tab["num"] == 1247))[0])
fig = dec.ver_evento(ev, i, tab)
fig.set_size_inches(8, 5.2)
fig.tight_layout()
fig.savefig(MUONES / "salida_decaimientos" / "reporte_ejemplo_decaimiento.pdf", bbox_inches="tight")
fig.savefig(MUONES / "salida_decaimientos" / "reporte_ejemplo_decaimiento.png", bbox_inches="tight", dpi=300)
print(ev.etiqueta(i), "dt =", round(float(tab["dt"].iloc[i])), "ns")
