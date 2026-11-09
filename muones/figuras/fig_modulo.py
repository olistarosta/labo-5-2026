"""Figura 1: interior de un módulo centellador (corte longitudinal).

Se muestra la mitad trasera del módulo cortado por un plano vertical que pasa
por su eje: así se ve que el centellador va desde la tapa del fondo hasta el
PMT y que el aluminio y la goma espuma lo recubren en todo su largo.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from iso import Camera, Scene, PAL, bezier, shade, LATEX_RC, tex, half_cylinder_x

plt.rcParams.update(LATEX_RC)

cam = Camera(az=22, el=26, light=(0.35, -0.8, 0.9))
sc = Scene(cam)

# --------------------------------------------------------------- geometría
# (sección transversal exagerada respecto del largo para que se vean las capas)
L, W = 80.0, 12.0          # largo y lado exterior de la caja
t, f, e = 0.8, 0.75, 0.5  # pared de fibrofácil, goma espuma, aluminio
cap = 3.0                  # tapas de madera maciza
yc = zc = W / 2            # plano de corte y eje del módulo
x_sc1 = cap + 46.0         # el centellador va de la tapa del fondo hasta acá
CUT = 10                   # capa de dibujo para las caras de corte
k_cut = [0]


def half_box(o, s, color, cut_color=None):
    """Caja cortada en y = yc: caras externas en capa 0, corte en capa CUT."""
    o, s = np.asarray(o, float), np.asarray(s, float)
    sc.box(o, s, color, layer=0, faces=["top", "back", "left", "right", "bottom"])
    x0, z0 = o[0], o[2]
    x1, z1 = o[0] + s[0], o[2] + s[2]
    k_cut[0] += 1
    sc.face([(x0, yc, z0), (x1, yc, z0), (x1, yc, z1), (x0, yc, z1)],
            shade(color, 0.08) if cut_color is None else cut_color,
            layer=CUT + 0.01 * k_cut[0], light=False, lw=0.5)


Lin = L - 2 * cap
# caja de fibrofácil y tapas
half_box((cap, yc, W - t), (Lin, W / 2, t), PAL["wood"])                 # techo
half_box((cap, yc, 0), (Lin, W / 2, t), PAL["wood"])                     # piso
sc.box((cap, W - t, t), (Lin, t, W - 2 * t), PAL["wood"], layer=0)       # pared trasera
half_box((0, yc, 0), (cap, W / 2, W), PAL["wood_cap"])                   # tapa del fondo
half_box((L - cap, yc, 0), (cap, W / 2, W), PAL["wood_cap"])             # tapa de salida
# goma espuma en todo el largo (techo, piso y fondo)
half_box((cap, yc, W - t - f), (Lin, W / 2 - t, f), PAL["foam"])
half_box((cap, yc, t), (Lin, W / 2 - t, f), PAL["foam"])
sc.box((cap, W - t - f, t + f), (Lin, f, W - 2 * t - 2 * f), PAL["foam"], layer=0)
# aluminio envolviendo todo el centellador (incluida la cara del fondo)
zi0, zi1 = t + f, W - t - f
half_box((cap, yc, zi0), (x_sc1 - cap, W / 2 - t - f, zi1 - zi0), PAL["foil"],
         cut_color=shade(PAL["foil"], 0.1))
half_box((cap + e, yc, zi0 + e), (x_sc1 - cap - e, W / 2 - t - f - e, zi1 - zi0 - 2 * e),
         PAL["scint"], cut_color=shade(PAL["scint"], 0.25))

# PMT (ventana + tubo) y base con el preamplificador
r = 3.4
half_cylinder_x(sc, x_sc1, x_sc1 + 0.6, yc, zc, r, PAL["glass"], layer=0,
                cut_layer=CUT + 0.5, cut_color=shade(PAL["glass"], 0.15))
half_cylinder_x(sc, x_sc1 + 0.6, 63.0, yc, zc, r, PAL["pmt"], layer=0,
                cut_layer=CUT + 0.5, cut_color="#454B54")
half_cylinder_x(sc, 63.0, 64.6, yc, zc, r + 0.45, PAL["base_ring"], layer=0,
                cut_layer=CUT + 0.5)
half_cylinder_x(sc, 64.6, 75.5, yc, zc, r + 0.1, PAL["base"], layer=0,
                cut_layer=CUT + 0.5)

# cables: salen de la base, atraviesan la tapa y siguen afuera
ends = {}
for name, dz, spread in [("hv", 2.0, 1.0), ("lv", 0.0, 0.0), ("sig", -2.0, -1.0)]:
    p0 = (75.5, yc, zc + dz)
    p1 = (L, yc, zc + dz)
    out = bezier(p1, (L + 7, yc, zc + dz), (L + 10, yc + 2 * spread, zc - 5 + 2.2 * dz),
                 (L + 17, yc + 3 * spread, zc - 10 + 3.2 * dz))
    sc.line3(np.vstack([p0, out]), PAL[name], lw=1.9, layer=CUT + 1)
    ends[name] = out[-1]

# mitad delantera removida: sólo un contorno muy suave
sc.box((0, 0, 0), (L, yc, W), "#FFFFFF", layer=CUT + 2, alpha=0.0,
       edge="#B9BEC5", lw=0.5, faces=["front", "top", "right"])

fig, ax = plt.subplots(figsize=(7.2, 2.7))
sc.draw(ax)
P = lambda v: cam.p(np.asarray(v, float))[0]

# muón
x_mu = 27.0
top, bot = np.array([x_mu, yc, W + 12]), np.array([x_mu, yc, -9])
ax.annotate("", xy=P(bot), xytext=P(top), zorder=60,
            arrowprops=dict(arrowstyle="-|>,head_width=0.25,head_length=0.5",
                            color=PAL["muon"], lw=1.6, shrinkA=0, shrinkB=0))
mx, my = P(top)
ax.text(mx - 1.3, my - 1.0, r"$\mu$", color=PAL["muon"], fontsize=14, ha="right",
        va="top", zorder=61)

# fotón: línea ondulada desde el punto de paso hasta la ventana del PMT
hit = np.array([x_mu, yc, zc + 1.2])
win = np.array([x_sc1 - 0.6, yc, zc + 0.3])
a, b = P(hit), P(win)
d = b - a
Ld = np.linalg.norm(d)
u, nrm = d / Ld, np.array([-d[1], d[0]]) / Ld
s = np.linspace(0, 1, 400)
taper = np.clip(s / 0.06, 0, 1) * np.clip((1 - s) / 0.10, 0, 1)
wave = a + np.outer(s * (Ld - 1.4), u) + np.outer(0.75 * np.sin(2 * np.pi * 7.5 * s) * taper, nrm)
ax.plot(wave[:, 0], wave[:, 1], color=PAL["photon"], lw=1.3, zorder=62,
        solid_capstyle="round")
ax.annotate("", xy=b, xytext=b - 1.6 * u, zorder=62,
            arrowprops=dict(arrowstyle="-|>,head_width=0.22,head_length=0.45",
                            color=PAL["photon"], lw=1.3, shrinkA=0, shrinkB=0))
ax.scatter(*a, s=26, color=PAL["photon"], ec="white", lw=0.7, zorder=63)

for name, p in ends.items():
    ax.scatter(*P(p), s=14, color=PAL[name], zorder=64, ec="white", lw=0.6)


# ---------------------------------------------------------------- etiquetas
def label(anchor, text, dx, dy, ha="center", va="bottom"):
    x, y = P(anchor)
    ax.plot([x, x + dx], [y, y + dy], color=PAL["ink_soft"], lw=0.6, zorder=70)
    ax.scatter([x], [y], s=6, color=PAL["ink"], zorder=71)
    off = 0.6 if ha == "left" else -0.6 if ha == "right" else 0
    voff = 0.5 if va == "bottom" else -0.5 if va == "top" else 0
    ax.text(x + dx + off, y + dy + voff, tex(text), ha=ha, va=va,
            color=PAL["ink"], zorder=72)


label((1.5, yc, W * 0.62), "Tapa de madera", -4, 9, ha="right")
label((12, W * 0.75, W), "Caja de fibrofácil", -2, 8)
label((44, yc, W - t - f / 2), "Gomaespuma", 0, 10.5)
label((71, yc, zc + 2.6), "Base + preamplificador", 0, 13.5)
label((12, yc, zi0 + e / 2), "Papel aluminio", -3, -9, va="top")
label((40, yc, zc - 2.2), "Centellador", 0, -11.5, va="top")
label((57, yc, zc - 2.0), "PMT", 0, -11.5, va="top")
for name, txt in [("hv", "HV"), ("lv", "±5 V"), ("sig", "Señal")]:
    x, y = P(ends[name])
    ax.text(x + 1.5, y, tex(txt), ha="left", va="center", color=PAL[name], zorder=72)

ax.set_aspect("equal")
ax.axis("off")
ax.autoscale_view()
ax.margins(0.02)
for ext, kw in [("pdf", {}), ("svg", {}), ("png", dict(dpi=400))]:
    fig.savefig(f"fig_modulo.{ext}", bbox_inches="tight", pad_inches=0.03, **kw)
print("ok")
