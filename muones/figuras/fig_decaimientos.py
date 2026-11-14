"""Figuras y números de la sección de vida media del informe (decaimientos de muones).

Todo con la misma ventana de tiempo entre pulsos, de 60 a 1400 ns:
  1. todos los decaimientos, ajuste de una exponencial;
  2. los que frenan en el hierro (entrada en u1 y u2), una y dos exponenciales;
  3. los que frenan en el centellador 3 (entrada en los tres), una y dos exponenciales.

El análisis está en analisis/dos_casos.py; acá solo se hacen las figuras y se imprimen
los números que van en el texto y en la tabla (informe/seccion_decaimientos.tex).

Uso:  python fig_decaimientos.py     (tarda unos 2 minutos, por las muestras simuladas)
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent / "analisis"))
import decaimientos as dec
import dos_casos as dc

DATOS, SERIE_2, AJUSTE = "#2A9D8F", "#E76F51", "#264653"
ARCHIVOS = {"Grupo 1: frena en el hierro": "decaimientos_hierro",
            "Grupo 2: cruza el hierro y frena en el centellador 3": "decaimientos_centellador3"}


def figura(x, p1, p2=None, archivo=None):
    """Histograma del tiempo entre pulsos con el ajuste de una exponencial y, si se pasa p2, de dos."""
    bordes = np.linspace(dc.A, dc.B, dc.N_BINS + 1)
    cuentas, ancho = np.histogram(x, bordes)[0], bordes[1] - bordes[0]
    t = np.linspace(dc.A, dc.B, 500)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(t / 1000, len(x) * ancho * dc.dens(t, np.asarray(p1)), color=AJUSTE, lw=2, zorder=2,
            label="Ajuste exponencial" if p2 is None else "Una exponencial")
    if p2 is not None:
        ax.plot(t / 1000, len(x) * ancho * dc.dens(t, p2), color=SERIE_2, lw=2, ls="--", zorder=2,
                label="Dos exponenciales")
    ax.errorbar(0.5 * (bordes[1:] + bordes[:-1]) / 1000, cuentas, yerr=np.sqrt(cuentas), fmt="o", ms=7,
                capsize=3, elinewidth=1.6, capthick=1.6, mew=1.6, color=DATOS, mfc=DATOS, mec=DATOS,
                zorder=3, label="Decaimientos")
    ax.set_xlabel(r"Tiempo entre pulsos $\Delta t$ [$\mu$s]")
    ax.set_ylabel("Número de eventos")          # por intervalo; el ancho va en el pie de figura
    ax.set_xlim(0, 1.45)
    ax.set_ylim(0, None)
    ax.legend(loc="upper right", facecolor="white", framealpha=0.9)
    fig.tight_layout()
    if archivo is not None:
        fig.savefig(AQUI / f"{archivo}.pdf", bbox_inches="tight")
        fig.savefig(AQUI / f"{archivo}.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def main():
    dc.estilo()
    grupos = dc.cargar_grupos()

    # 1. todos los decaimientos
    ev = dec.cargar_indice()
    tab = dec.tabla_eventos(ev)
    sel = dec.combinar(dec.armar_filtros(tab, dt_rango=(dc.A, dc.B)), dc.FILTROS)
    x = tab.loc[sel, "dt"].to_numpy()
    una = dec.ajustar_vida_media(x, (dc.A, dc.B))
    c, g, p = dc.chi2_binado(x, (una["tau"],))
    print(f"{len(tab)} eventos registrados; ventana de {dc.A:.0f} a {dc.B:.0f} ns")
    print(f"TODOS: {una['n']} decaimientos, tau = {una['tau']:.0f} ns (+{una['err_mas']:.0f} / -{una['err_menos']:.0f}),"
          f" chi2 = {c:.1f} / {g} (p = {p:.2f})")
    print(f"   suma de los dos grupos: {sum(len(v) for v in grupos.values())}")
    figura(x, (una["tau"],), archivo="decaimientos_todos")

    # 2 y 3. los dos grupos
    for nombre, x in grupos.items():
        r = dc.analizar(x)
        dc.imprimir(nombre, r)
        figura(x, r["p1"], r["p2"], archivo=ARCHIVOS[nombre])
    print("\nfiguras guardadas en", AQUI)


if __name__ == "__main__":
    main()
