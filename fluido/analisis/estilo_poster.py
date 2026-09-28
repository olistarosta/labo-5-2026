"""Estilo de los gráficos del póster: fondo oscuro como los nodos de TouchDesigner.

Uso, en la primera celda:

    from estilo_poster import *
    usar_estilo(fontsize=25)

Los colores salen de poster/poster_fluidos_A0.pdf. La fuente es Inter (la del póster; se baja de
github.com/rsms/inter). Si no está instalada, se usa Segoe UI.
"""
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

FONDO = "#3a3a3a"     # gris de los nodos del póster: el casi negro de los colormaps se distingue
PANEL = "#2f2f2f"     # leyendas y recuadros
TEXTO = "#ececec"
TEXTO_2 = "#b4b4b4"   # texto secundario, lineas de referencia
BORDE = "#7a7a7a"
GRILLA = "#4a4a4a"
BANDA = "#5f636b"     # bandas de incerteza
FLECHAS = "#05080f"   # flechas sobre el mapa de |v| (el extremo oscuro del colormap)

NARANJA = "#e38b2c"   # las palabras clave del póster
AZUL = "#5fb4f5"      # el celeste #9fd0ff de los valores, más saturado para que no se lea gris
VIOLETA = "#a79fdc"   # los conectores entre nodos

# |v|: secuencial "hielo", de casi negro a casi blanco pasando por el azul y el celeste:
# usa todo el rango de luminancia. Las flechas van oscuras (FONDO) encima.
CMAP_RAPIDEZ = LinearSegmentedColormap.from_list(
    "rapidez", ["#05080f", "#122449", "#1d4a8f", "#2a78c4", "#4fb0e6", "#9fdcf2", "#effaff"])
# vorticidad: divergente con el centro oscuro (0 se funde con el fondo)
# naranja = antihorario (> 0), celeste = horario (< 0)
CMAP_VORTICIDAD = LinearSegmentedColormap.from_list(
    "vorticidad", ["#cfeaff", AZUL, "#1d4f7a", "#161b24", "#7a4412", NARANJA, "#ffd9a8"])
for _cmap in (CMAP_RAPIDEZ, CMAP_VORTICIDAD):   # registrados, para poder usarlos por nombre
    if _cmap.name not in matplotlib.colormaps:
        matplotlib.colormaps.register(_cmap)


def usar_estilo(fontsize=25):
    plt.style.use("default")
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Inter", "Segoe UI", "DejaVu Sans"],
        "mathtext.fontset": "custom",   # las fórmulas también en Inter
        "mathtext.rm": "Inter",
        "mathtext.it": "Inter:italic",
        "mathtext.bf": "Inter:bold",
        "mathtext.default": "regular",
        "font.size": fontsize,
        "axes.labelsize": fontsize,
        "axes.titlesize": fontsize,
        "axes.titleweight": "semibold",
        "axes.titlepad": 14,
        "figure.titlesize": fontsize,
        "figure.titleweight": "semibold",
        "xtick.labelsize": fontsize - 2,
        "ytick.labelsize": fontsize - 2,
        "legend.fontsize": fontsize - 4,

        "figure.facecolor": FONDO,
        "savefig.facecolor": FONDO,
        "axes.facecolor": FONDO,
        "axes.edgecolor": BORDE,
        "axes.linewidth": 1.2,
        "axes.labelcolor": TEXTO,
        "axes.titlecolor": TEXTO,
        "text.color": TEXTO,
        "xtick.color": TEXTO_2,
        "ytick.color": TEXTO_2,
        "xtick.labelcolor": TEXTO,
        "ytick.labelcolor": TEXTO,
        "axes.grid": True,
        "grid.color": GRILLA,
        "grid.linewidth": 1,
        "axes.spines.top": False,
        "axes.spines.right": False,

        "legend.frameon": True,
        "legend.facecolor": PANEL,
        "legend.edgecolor": BORDE,
        "legend.framealpha": 0.9,
        "legend.labelcolor": TEXTO,

        "axes.prop_cycle": plt.cycler(color=[AZUL, NARANJA, VIOLETA]),
        "image.cmap": "rapidez",
    })
