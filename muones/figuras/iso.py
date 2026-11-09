"""Mini motor de dibujo axonométrico para matplotlib.

Proyección ortográfica con cámara (azimut, elevación), sombreado plano por
normal de cara y orden del pintor. Pensado para esquemas limpios de
laboratorio: cajas, cilindros, varillas y curvas en 3D, todo vectorial.
"""
import numpy as np
from matplotlib.patches import Polygon
from matplotlib.colors import to_rgb

# ---------------------------------------------------------------- paleta
PAL = dict(
    ink="#2B2F36",        # texto y líneas de referencia
    ink_soft="#6B7280",
    wood="#E3C08D",       # fibrofácil
    wood_cap="#C99A5E",   # tapas de madera maciza
    foam="#3B4048",       # goma espuma negra
    foil="#C8CDD4",       # papel aluminio
    scint="#7FD6E8",      # centellador acrílico
    pmt="#2E333B",        # cuerpo PMT
    glass="#5FA8D3",
    base_ring="#1E9C8A",  # base UCLA (anillo verde)
    base="#C3C9D1",       # cuerpo plateado de la base
    hv="#E8590C",         # cable alta tensión
    lv="#495057",         # cable ±5 V
    sig="#1C7ED6",        # cable de señal
    muon="#D6336C",
    photon="#F2B705",
    iron="#8A94A3",
    stand="#3A3F46",
    rod="#555B64",
    clamp="#23272D",
    drs="#BCC3CC",
    sma="#D4A72C",
    rpi="#2F9E57",
    cab="#5E6E5B",        # gabinete de fuentes
    weight="#4D545E",
)


def shade(color, f):
    """f>0 aclara hacia blanco, f<0 oscurece hacia negro."""
    c = np.array(to_rgb(color))
    return tuple(c + (1 - c) * f) if f >= 0 else tuple(c * (1 + f))


class Camera:
    def __init__(self, az=-35, el=28, light=(-0.45, -0.75, 0.9)):
        self.az, self.el = np.radians(az), np.radians(el)
        ca, sa = np.cos(self.az), np.sin(self.az)
        ce, se = np.cos(self.el), np.sin(self.el)
        # base de la cámara en coordenadas de mundo
        self.right = np.array([ca, sa, 0.0])
        fwd = np.array([-sa, ca, 0.0])               # hacia el fondo
        self.up = fwd * se + np.array([0, 0, 1]) * ce
        self.toward = -fwd * ce + np.array([0, 0, 1]) * se   # hacia la cámara
        L = np.array(light, float)
        self.light = L / np.linalg.norm(L)

    def p(self, pts):
        pts = np.atleast_2d(pts)
        return np.c_[pts @ self.right, pts @ self.up]

    def depth(self, pts):
        return np.atleast_2d(pts) @ self.toward


class Scene:
    """Acumula caras y las dibuja con el algoritmo del pintor.

    Cada objeto agrega caras con un 'layer' (orden grueso, manual) y la
    profundidad del centroide decide dentro del mismo layer.
    """

    def __init__(self, cam):
        self.cam = cam
        self.items = []   # (layer, depth, kind, payload)

    # ------------------------------------------------------------ primitivas
    def face(self, pts, color, layer=0, alpha=1.0, edge=None, lw=0.5,
             cull=True, light=True, depth=None, ls="-"):
        pts = np.asarray(pts, float)
        n = np.cross(pts[1] - pts[0], pts[2] - pts[0])
        nn = np.linalg.norm(n)
        if nn == 0:
            return
        n /= nn
        if cull and n @ self.cam.toward <= 1e-9:
            return
        if light:
            lam = max(0.0, n @ self.cam.light)
            fc = shade(color, -0.30 + 0.42 * lam)
        else:
            fc = to_rgb(color)
        ec = shade(fc, -0.28) if edge is None else edge
        d = self.cam.depth(pts.mean(0))[0] if depth is None else depth
        self.items.append((layer, d, "poly",
                           dict(xy=self.cam.p(pts), fc=fc, ec=ec, lw=lw,
                                alpha=alpha, ls=ls)))

    def box(self, o, s, color, layer=0, alpha=1.0, edge=None, lw=0.5,
            R=None, faces="all", cull=True, depth=None):
        """Caja de origen o (esquina mínima, local) y tamaño s.
        R: matriz 3x3 de rotación y traslación opcional vía (R, t)."""
        x, y, z = s
        v = np.array([[0, 0, 0], [x, 0, 0], [x, y, 0], [0, y, 0],
                      [0, 0, z], [x, 0, z], [x, y, z], [0, y, z]], float)
        v = v + np.asarray(o, float)
        if R is not None:
            Rm, t = R
            v = v @ np.asarray(Rm).T + np.asarray(t)
        F = dict(bottom=[0, 3, 2, 1], top=[4, 5, 6, 7], front=[0, 1, 5, 4],
                 back=[2, 3, 7, 6], left=[0, 4, 7, 3], right=[1, 2, 6, 5])
        keys = F.keys() if faces == "all" else faces
        dd = self.cam.depth(v.mean(0))[0] if depth is None else depth
        for k in keys:
            self.face(v[F[k]], color, layer, alpha, edge, lw, cull=cull,
                      depth=dd)
        return v

    def cylinder(self, c0, c1, r, color, layer=0, n=40, caps=(True, True),
                 alpha=1.0, cap_color=None, depth=None):
        c0, c1 = np.asarray(c0, float), np.asarray(c1, float)
        ax = c1 - c0
        ax /= np.linalg.norm(ax)
        a = np.cross(ax, [0, 0, 1]) if abs(ax[2]) < 0.9 else np.cross(ax, [1, 0, 0])
        a /= np.linalg.norm(a)
        b = np.cross(ax, a)
        th = np.linspace(0, 2 * np.pi, n, endpoint=False)
        ring = np.cos(th)[:, None] * a + np.sin(th)[:, None] * b
        P0, P1 = c0 + r * ring, c1 + r * ring
        dd = self.cam.depth((c0 + c1) / 2)[0] if depth is None else depth
        for i in range(n):
            j = (i + 1) % n
            q = [P0[i], P0[j], P1[j], P1[i]]
            self.face(q, color, layer, alpha, edge="face", lw=0.4, depth=dd)
        cc = color if cap_color is None else cap_color
        if caps[0]:
            self.face(P0[::-1], cc, layer, alpha, lw=0.5, depth=dd + 1e-3)
        if caps[1]:
            self.face(P1, cc, layer, alpha, lw=0.5, depth=dd + 1e-3)

    def line3(self, pts, color, lw=1.2, layer=0, ls="-", depth=None,
              alpha=1.0, capstyle="round", zshift=0.0):
        pts = np.atleast_2d(np.asarray(pts, float))
        d = self.cam.depth(pts.mean(0))[0] if depth is None else depth
        self.items.append((layer, d + zshift, "line",
                           dict(xy=self.cam.p(pts), c=color, lw=lw, ls=ls,
                                alpha=alpha, cap=capstyle)))

    # ------------------------------------------------------------ render
    def draw(self, ax):
        self.items.sort(key=lambda t: (t[0], t[1]))
        for k, (layer, d, kind, P) in enumerate(self.items):
            z = 1 + k * 1e-4
            if kind == "poly":
                ec = P["fc"] if P["ec"] == "face" else P["ec"]
                ax.add_patch(Polygon(P["xy"], closed=True, fc=P["fc"], ec=ec,
                                     lw=P["lw"], alpha=P["alpha"], zorder=z,
                                     joinstyle="round", ls=P["ls"]))
            else:
                ax.plot(P["xy"][:, 0], P["xy"][:, 1], color=P["c"], lw=P["lw"],
                        ls=P["ls"], alpha=P["alpha"], zorder=z,
                        solid_capstyle=P["cap"], dash_capstyle=P["cap"])
        self.items = []


def bezier(p0, p1, p2, p3, n=60):
    t = np.linspace(0, 1, n)[:, None]
    p0, p1, p2, p3 = map(lambda v: np.asarray(v, float), (p0, p1, p2, p3))
    return ((1 - t) ** 3) * p0 + 3 * ((1 - t) ** 2) * t * p1 + \
        3 * (1 - t) * t ** 2 * p2 + t ** 3 * p3


def rot_z(deg):
    a = np.radians(deg)
    return np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])


def rot_x(deg):
    a = np.radians(deg)
    return np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])


# ---------------------------------------------------------------- tipografía
# Computer Modern (la fuente de LaTeX) vía mathtext: no hace falta LaTeX
# instalado y los acentos/ñ se arman como en LaTeX.
LATEX_RC = {
    "mathtext.fontset": "cm", "font.family": "serif", "font.serif": ["cmr10"],
    "axes.formatter.use_mathtext": True, "font.size": 10,
    "pdf.fonttype": 42, "svg.fonttype": "path",
}
_ACC = {"á": r"\acute{a}", "é": r"\acute{e}", "í": r"\acute{\imath}",
        "ó": r"\acute{o}", "ú": r"\acute{u}", "ñ": r"\tilde{n}", "±": r"{\pm}",
        " ": r"\ ", "μ": r"\mu", "θ": r"\theta", "+": r"{+}"}


def tex(s, bf=False):
    body = "".join(_ACC.get(c, c) for c in s)
    return rf"$\mathbf{{{body}}}$" if bf else rf"$\mathrm{{{body}}}$"


def half_cylinder_x(sc, x0, x1, yc, zc, r, color, layer=0, cut_layer=None,
                    cut_color=None, n=40):
    """Medio cilindro de eje x, sólo la mitad y >= yc (corte en y = yc).
    Las caras curvas van en `layer`; la cara de corte en `cut_layer`."""
    th = np.linspace(-np.pi / 2, np.pi / 2, n)
    ring = np.c_[np.zeros(n), np.cos(th), np.sin(th)] * r      # y = r cos >= 0
    P0 = np.array([x0, yc, zc]) + ring
    P1 = np.array([x1, yc, zc]) + ring
    for i in range(n - 1):
        sc.face([P0[i], P1[i], P1[i + 1], P0[i + 1]], color, layer, edge="face", lw=0.3)
    sc.face(P1[::-1], color, layer, lw=0.4)                      # tapa x1
    sc.face(P0, color, layer, lw=0.4)                            # tapa x0
    cc = shade(color, 0.12) if cut_color is None else cut_color
    sc.face([(x0, yc, zc - r), (x1, yc, zc - r), (x1, yc, zc + r), (x0, yc, zc + r)],
            cc, layer if cut_layer is None else cut_layer, light=False, lw=0.5)
