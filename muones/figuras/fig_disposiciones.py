"""Disposiciones experimentales, una figura por archivo (misma escala en las tres):

  disp_apilados     módulos apilados con hierro entre 2 y 3 + adquisición
  disp_inclinados   en soportes, inclinados un ángulo θ (barrido angular)
  disp_velocidades  en soportes con barras horizontales, cruzados a 60°,
                    con los centros de las tres barras sobre la misma vertical

Uso:  python fig_disposiciones.py            -> las tres
      python fig_disposiciones.py inclinados -> sólo una
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from iso import Camera, Scene, PAL, bezier, rot_z, rot_x, LATEX_RC, tex

plt.rcParams.update(LATEX_RC)
FS = 1.7               # factor de tamaño de letra (para que se lea con las tres en fila)
plt.rcParams["font.size"] = 10 * FS

LM, WM = 80.0, 6.0     # largo y lado de un módulo (sección algo exagerada)
CAP = 3.0
X_SC = 26.0            # centro del centellador medido desde la tapa del fondo
BASE = (14, 20, 2.2)   # pie de soporte
SCALE = 50.0           # unidades de dibujo por pulgada (igual en las tres figuras)
TX = SCALE / 24.0 * FS # cuánto más grande es el texto respecto del dibujo (para usarlas en fila)


# ------------------------------------------------------------- utilidades
def module(sc, R=np.eye(3), t=(0, 0, 0), layer=1):
    """Módulo con origen en la esquina (tapa del fondo, frente, abajo)."""
    T = (R, t)
    side = ["top", "bottom", "front", "back"]
    # sólo caras exteriores: las caras de contacto tapa-cuerpo no se dibujan
    v = sc.box((CAP, 0, 0), (LM - 2 * CAP, WM, WM), PAL["wood"], layer=layer, R=T,
               faces=side)
    d = sc.cam.depth(v.mean(0))[0]
    sc.box((0, 0, 0), (CAP, WM, WM), PAL["wood_cap"], layer=layer, R=T,
           depth=d + 0.01, faces=side + ["left"])
    sc.box((LM - CAP, 0, 0), (CAP, WM, WM), PAL["wood_cap"], layer=layer, R=T,
           depth=d + 0.01, faces=side + ["right"])


def shadow(sc, pts2, a=0.09):
    P = np.c_[np.asarray(pts2, float), np.zeros(len(pts2))]
    sc.items.append((-1, -1e9, "poly",
                     dict(xy=sc.cam.p(P), fc=(0, 0, 0), ec="none", lw=0,
                          alpha=a, ls="-")))


def badge(ax, xy, n):
    ax.add_patch(Circle(xy, 2.3 * TX, fc="white", ec=PAL["ink"], lw=0.8, zorder=80))
    ax.text(xy[0], xy[1] - 0.15, tex(str(n), bf=True), ha="center", va="center",
            fontsize=9 * FS, color=PAL["ink"], zorder=81)


def label(ax, xy, text, dx, dy, ha="left", va="center"):
    x, y = xy
    ax.plot([x, x + dx], [y, y + dy], color=PAL["ink_soft"], lw=0.6, zorder=90)
    ax.scatter([x], [y], s=6, color=PAL["ink"], zorder=91)
    off = 0.8 if ha == "left" else -0.8 if ha == "right" else 0
    voff = 0.6 if va == "bottom" else -0.6 if va == "top" else 0
    ax.text(x + dx + off, y + dy + voff, tex(text), ha=ha, va=va, zorder=92,
            color=PAL["ink"])


def muon_arrow(ax, cam, p_top, p_bot, side=1, lab_at=None):
    a, b = cam.p(p_top)[0], cam.p(p_bot)[0]
    ax.annotate("", xy=b, xytext=a, zorder=70,
                arrowprops=dict(arrowstyle="-|>,head_width=0.22,head_length=0.45",
                                color=PAL["muon"], lw=1.4, shrinkA=0, shrinkB=0))
    q = a if lab_at is None else cam.p(lab_at)[0]
    ax.text(q[0] + 1.6 * side, q[1] - 0.5, r"$\mu$", color=PAL["muon"],
            fontsize=13 * FS, ha="left" if side > 0 else "right", va="top", zorder=71)


def clamp(sc, W, x, z, y_from, y_to, layer):
    """Agarradera: nuez sobre la varilla + brazo hasta el módulo."""
    R, Wf = W
    sc.box((0, 0, 0), (3.0, 3.0, 3.0), PAL["clamp"], layer=layer + 0.01,
           R=(R, Wf((x - 1.5, y_to - 1.5, z - 1.5))))
    sc.line3([Wf((x, y_from, z)), Wf((x, y_to, z))], PAL["clamp"], lw=1.8, layer=layer)


def center_mark(sc, z, layer):
    """Hexágono magenta en la tapa superior: zona común a los tres centelladores."""
    rh = (WM * 0.75 / 2) / np.cos(np.radians(30))
    hx = [(rh * np.cos(np.radians(60 * k)), rh * np.sin(np.radians(60 * k)), z + WM + 0.05)
          for k in range(6)]
    sc.face(hx, PAL["muon"], layer=layer, alpha=0.85, light=False, edge="none", lw=0)


# ======================================================== apilados + DAQ
def apilados(ax):
    cam = Camera(az=26, el=24, light=(0.4, -0.8, 0.9))
    sc = Scene(cam)
    z_iron = 1.6
    zs = {3: 0.0, 2: WM + z_iron, 1: 2 * WM + z_iron}
    shadow(sc, [(-1, -2.5), (LM + 1, -2.5), (LM + 1, WM + 2.5), (-1, WM + 2.5)])
    module(sc, t=(0, 0, zs[3]), layer=1.00)
    sc.box((0, -0.6, WM), (52, WM + 1.2, z_iron), PAL["iron"], layer=1.01)  # hierro
    module(sc, t=(0, 0, zs[2]), layer=1.02)
    module(sc, t=(0, 0, zs[1]), layer=1.03)

    # fuentes, DRS4 y Raspberry Pi
    shadow(sc, [(95, 19), (123, 19), (123, 39), (95, 39)])
    sc.box((96, 20, 0), (26, 18, 18), PAL["cab"], layer=1)
    for k in range(3):
        c = np.array([101 + 8 * k, 20, 12.5])
        sc.cylinder(c, c - (0, 1.0, 0), 1.3, "#D9DDE2", layer=1.2, n=24)
        c2 = np.array([101 + 8 * k, 20, 6.5])
        sc.cylinder(c2, c2 - (0, 0.6, 0), 0.8, PAL["ink"], layer=1.2, n=20)
    shadow(sc, [(95.5, -6.5), (113.5, -6.5), (113.5, 4.5), (95.5, 4.5)])
    sc.box((96, -6, 0), (17, 10, 3.2), PAL["drs"], layer=1)
    sc.cylinder((113, -1, 1.6), (114.2, -1, 1.6), 0.7, PAL["ink_soft"], layer=1.2, n=16)
    shadow(sc, [(124, -9), (132.5, -9), (132.5, -3.4), (124, -3.4)])
    sc.box((124, -9, 0.6), (8.5, 5.6, 0.4), PAL["rpi"], layer=1)
    sc.box((126.8, -7.6, 1.0), (2.4, 2.4, 0.5), "#22262B", layer=1.3)
    sc.box((130.6, -8.6, 1.0), (1.8, 2.0, 1.4), "#C9CED6", layer=1.3)
    sc.box((130.6, -6.0, 1.0), (1.8, 2.0, 1.4), "#C9CED6", layer=1.3)
    sc.box((124.6, -8.7, 1.0), (5.0, 0.6, 0.7), "#22262B", layer=1.3)

    for n, z in zs.items():
        pe = np.array([LM, WM * 0.5, z + WM * 0.5])
        sig = bezier(pe + (0, -1, -0.8), pe + (8, -1.5, -1), (90, -1 + n, 2.0 + 0.4 * n),
                     (96.2, -1.2 + 2.2 * (n - 2), 1.6))
        sc.line3(sig, PAL["sig"], lw=1.3, layer=0.9, depth=-1e6)
        pw = bezier(pe + (0, 1, 0.8), pe + (6, 4, 3), (92, 18, 10 + n),
                    (101 + 8 * (n - 1), 20, 9.5))
        sc.line3(pw, PAL["lv"], lw=1.3, layer=0.8, depth=-1e6)
    usb = bezier((114.2, -1, 1.6), (119, -1, 1.6), (121, -5, 0.8), (124, -6.2, 0.8))
    sc.line3(usb, "#868E96", lw=1.2, layer=1.25)

    sc.draw(ax)
    top = np.array([X_SC, WM / 2, 3 * WM + z_iron + 24])
    muon_arrow(ax, cam, top, np.array([X_SC, WM / 2, -7]), side=-1)
    P = lambda v: cam.p(np.asarray(v, float))[0]
    p2 = P((-0.2, 0, zs[2] + WM / 2))
    for n, z in zs.items():
        pe = P((-0.2, 0, z + WM / 2))
        pb = np.array([p2[0] - 9.0 * TX, p2[1] + (2 - n) * 2.6 * 2.3 * TX])
        ax.plot([pe[0], pb[0]], [pe[1], pb[1]], color=PAL["ink_soft"], lw=0.6, zorder=79)
        badge(ax, pb, n)
    label(ax, P((44, -0.6, WM + z_iron / 2)), "Hierro", 3, -13, ha="center", va="top")
    label(ax, P((118, 29, 18)), "Fuentes HV y ±5 V", 0, 10, ha="center", va="bottom")
    label(ax, P((104.5, -6, 1.6)), "DRS4", -7, -6, ha="center", va="top")
    label(ax, P((128, -9, 0.8)), "Raspberry Pi", 4, -7, ha="left", va="top")


# ===================================================== inclinados θ
def inclinados(ax, TH=30.0):
    cam = Camera(az=-70, el=15, light=(-0.3, -0.7, 0.9))
    sc = Scene(cam)
    y_rod = WM + 2.4                      # varillas del lado hacia el que se inclina
    y_bc = y_rod - 3                      # centro de las bases
    piv = np.array([0, y_bc + BASE[1] / 2, 0])   # giran sobre el borde trasero
    R = rot_x(-TH)                        # inclina hacia +y
    Wt = lambda p: (np.asarray(p, float) - piv) @ R.T + piv
    zb = {3: 7.0, 2: 7.0 + 26.1, 1: 7.0 + 26.1 + 29.5}
    xs_st = (14.0, 66.0)
    H_rod = 80.0

    # pie de apoyo (vertical) con travesaño y contrapeso
    c_rod = Wt((0, y_rod, 66.0))
    y3, z3 = c_rod[1] + 1.8, c_rod[2]
    x3 = 40.0
    shadow(sc, [(x3 - 8, y3 - 11), (x3 + 8, y3 - 11), (x3 + 8, y3 + 11), (x3 - 8, y3 + 11)])
    sc.box((x3 - 7, y3 - 10, 0), BASE, PAL["stand"], layer=0)
    sc.box((x3 - 5.5, y3 + 1.5, 2.2), (11, 8, 7), PAL["weight"], layer=0.05)
    sc.line3([(x3, y3, 2.2), (x3, y3, z3 + 12)], PAL["rod"], lw=2.2, layer=0.1)
    sc.line3([(2, y3 - 0.9, z3), (78, y3 - 0.9, z3)], PAL["rod"], lw=2.2, layer=0.15)
    sc.box((x3 - 1.6, y3 - 2.2, z3 - 1.6), (3.2, 3.2, 3.2), PAL["clamp"], layer=0.2)

    # estructura inclinada
    for xs in xs_st:
        o = np.array([xs - BASE[0] / 2, y_bc - BASE[1] / 2, 0])
        sc.box((0, 0, 0), BASE, PAL["stand"], layer=0.5, R=(R, Wt(o)))
        sc.line3([Wt((xs, y_rod, 2.2)), Wt((xs, y_rod, H_rod))], PAL["rod"], lw=2.2,
                 layer=0.6)
    for n, z in zb.items():
        for xs in xs_st:
            clamp(sc, (R, Wt), xs, z + WM / 2, WM, y_rod, layer=0.7)
        module(sc, R=R, t=Wt((0, 0, z)), layer=1 + z / 1000)

    # cono de aceptancia: todas las rectas que cruzan los centelladores 1 y 3,
    # en el plano perpendicular a las barras que pasa por el centro del centellador.
    # Va detrás de los módulos para que quede prolijo.
    w0, w1 = WM * 0.08, WM * 0.92            # ancho del centellador dentro de la caja
    za, zb_ = zb[1] + WM, zb[3]              # cara superior de 1, inferior de 3
    D = za - zb_
    a0, a1 = np.array([w0, za]), np.array([w1, za])
    b0, b1 = np.array([w0, zb_]), np.array([w1, zb_])
    Ht, Hb = 20.0, 4.0                       # cuánto sobresale arriba y abajo
    tl, tr = a0 + (a0 - b1) * Ht / D, a1 + (a1 - b0) * Ht / D
    br, bl = b1 + (b1 - a0) * Hb / D, b0 + (b0 - a1) * Hb / D
    W3 = lambda yz: Wt((X_SC, yz[0], yz[1]))
    cone = np.array([W3(p) for p in [tl, tr, a1, b1, br, bl, b0, a0]])
    sc.items.append((0.95, 0.0, "poly", dict(xy=cam.p(cone), fc=PAL["muon"], ec="none",
                                             lw=0, alpha=0.16, ls="-")))
    for p, q in [(tl, br), (tr, bl)]:
        sc.line3([W3(p), W3(q)], PAL["muon"], lw=0.8, layer=0.96, ls=(0, (3, 2)),
                 alpha=0.75)

    sc.draw(ax)
    P = lambda v: cam.p(np.asarray(v, float))[0]
    ax_d = R @ np.array([0, 0, 1.0])
    c_top = Wt((X_SC, WM / 2, zb[1] + WM / 2))
    c_bot = Wt((X_SC, WM / 2, zb[3] + WM / 2))
    vtx = c_top + ax_d * 30                  # vértice del ángulo, por encima del cono
    muon_arrow(ax, cam, vtx + ax_d * 14, c_bot - ax_d * 10, side=-1)
    ez, ey = np.array([0, 0, 1.0]), np.array([0, 1.0, 0])
    ax.plot(*np.c_[P(vtx), P(vtx + ez * 26)], color=PAL["ink_soft"], lw=0.8,
            ls=(0, (3, 2)), zorder=69)
    rr = 17.0
    arc = np.array([vtx + rr * (np.cos(a) * ez + np.sin(a) * ey)
                    for a in np.radians(np.linspace(0, TH, 30))])
    ax.plot(*cam.p(arc).T, color=PAL["ink"], lw=0.9, zorder=69)
    am = vtx + (rr + 5) * (np.cos(np.radians(TH / 2)) * ez + np.sin(np.radians(TH / 2)) * ey)
    ax.text(*P(am), r"$\theta$", ha="center", va="center", fontsize=13,
            color=PAL["ink"], zorder=72)

    for n, z in zb.items():
        badge(ax, P(Wt((LM, 0, z + WM / 2))) + (3.6 * TX, 0), n)
    label(ax, P(Wt((xs_st[1], y_rod, H_rod - 2))), "Soporte", 5, 4, ha="left", va="bottom")
    label(ax, P((x3, y3, z3 * 0.55)), "Pie de apoyo", -8, 0, ha="right", va="center")
    label(ax, P((x3, y3 + 9.5, 6)), "Contrapeso", -4, -7, ha="right", va="top")


# ======================================== velocidades (cruzados a 60°)
def velocidades(ax):
    cam = Camera(az=-30, el=33, light=(0.3, -0.8, 0.9))
    sc = Scene(cam)
    zc = {3: 8.0, 2: 8.0 + 32.5, 1: 8.0 + 32.5 + 27.7}   # base de cada módulo
    phis = {1: 30.0, 2: 150.0, 3: 90.0}                   # 60° entre sí
    XR, H_rod, r_rung = 46.0, 84.0, 0.8

    # dos pies con varilla vertical (como en la foto)
    for xs in (-XR, XR):
        shadow(sc, [(xs - 8, -11), (xs + 8, -11), (xs + 8, 11), (xs - 8, 11)])
        sc.box((xs - 7, -10, 0), BASE, PAL["stand"], layer=0.2)
        sc.line3([(xs, 0, 2.2), (xs, 0, H_rod)], PAL["rod"], lw=2.2,
                 layer=2.0 if xs < 0 else 0.3)

    for n, z in zc.items():
        zr = z - r_rung
        # barra horizontal entre las dos varillas, pasa por debajo del centro
        sc.line3([(-XR, 0, zr), (XR, 0, zr)], PAL["rod"], lw=2.2,
                 layer=1 + (z - 1) / 1000)
        for xs in (-XR, XR):
            sc.box((xs - 1.6, -1.6, zr - 1.6), (3.2, 3.2, 3.2), PAL["clamp"],
                   layer=2.1 if xs < 0 else 1 + (z - 1) / 1000 + 0.0001)
        # módulo apoyado con su centro sobre el eje vertical
        Rz = rot_z(phis[n])
        t = np.array([0, 0, z]) - Rz @ np.array([LM / 2, WM / 2, 0])
        module(sc, R=Rz, t=t, layer=1 + z / 1000)
        center_mark(sc, z, layer=1 + z / 1000 + 0.0005)

    sc.draw(ax)
    P = lambda v: cam.p(np.asarray(v, float))[0]
    muon_arrow(ax, cam, np.array([0, 0, zc[1] + WM + 30]), np.array([0, 0, -4]), side=-1)
    for n, z in zc.items():
        Rz = rot_z(phis[n])
        sgn = 1 if n == 1 else -1
        end = Rz @ np.array([sgn * LM / 2, 0, 0]) + (0, 0, z + WM / 2)
        d = P(end) - P((0, 0, z + WM / 2))
        badge(ax, P(end) + 3.6 * TX * d / np.linalg.norm(d), n)
    label(ax, P((XR, 0, H_rod - 3)), "Soporte", 4, 4, ha="left", va="bottom")


# ------------------------------------------------------------- salida
FIGS = dict(apilados=apilados, inclinados=inclinados, velocidades=velocidades)


def render(name):
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1])
    FIGS[name](ax)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.autoscale_view()
    ax.margins(0.02)
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    fig.set_size_inches((x1 - x0) / SCALE, (y1 - y0) / SCALE)
    for ext, kw in [("pdf", {}), ("svg", {}), ("png", dict(dpi=400))]:
        fig.savefig(f"disp_{name}.{ext}", bbox_inches="tight", pad_inches=0.03, **kw)
    plt.close(fig)


if __name__ == "__main__":
    for name in (sys.argv[1:] or FIGS):
        render(name)
        print("ok", name)
