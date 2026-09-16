# -*- coding: utf-8 -*-
"""
PIV por FFT con deformacion de ventana, multipasada. Mismo metodo que PIVlab
("Multipass FFT window deformation", https://www.pivlab.de/wiki/1-quickstart/):

  1. Preprocesamiento: CLAHE (ecualizacion adaptativa de contraste).
  2. Correlacion cruzada por FFT en ventanas de interrogacion con 50 % de
     solapamiento, en varias pasadas de ventana decreciente (128 -> 64 -> 32).
     El desplazamiento de cada pasada deforma las imagenes para la siguiente,
     que solo tiene que medir el residuo.
  3. Pico sub-pixel por ajuste gaussiano de 3 puntos en cada eje.
  4. Validacion: filtro de desvio estandar (n = 8) y filtro de mediana local
     (umbral 3). Los vectores descartados se reemplazan por interpolacion y
     quedan marcados.

Trabaja en DESPLAZAMIENTOS en pixeles entre dos imagenes. Pasar a velocidad
(escala en mm/px y paso temporal) es trabajo de quien llama.

Necesita numpy y cv2 (TouchDesigner trae los dos).
"""

import cv2
import numpy as np


# ------------------------------------------------------------ preprocesamiento

def preprocesar(img, clahe=True, ventana_clahe=64):
    """Imagen en gris -> float32 0..1, con CLAHE si se pide.

    CLAHE ecualiza el contraste por zonas de `ventana_clahe` px: iguala destellos
    en zonas bien y mal iluminadas, que es lo que mas estabiliza la correlacion.
    """
    a = np.asarray(img)
    if a.ndim == 3:
        a = a[..., 0]
    if a.dtype != np.uint8:
        a = (np.clip(a, 0.0, 1.0) * 255.0).astype(np.uint8)
    if clahe:
        h, w = a.shape
        tiles = (max(1, round(w / ventana_clahe)), max(1, round(h / ventana_clahe)))
        a = cv2.createCLAHE(clipLimit=2.0, tileGridSize=tiles).apply(a)
    return a.astype(np.float32) / 255.0


# -------------------------------------------------------------- utilidades

def _bloques(img, win, step):
    v = np.lib.stride_tricks.sliding_window_view(img, (win, win))
    return v[::step, ::step]


def _grilla(n, win, step):
    return (np.arange(0, n - win + 1, step) + (win - 1) * 0.5).astype(np.float32)


def _interp_campo(f, gy, gx, ty, tx):
    """Interpolacion bilineal de un campo de una grilla regular a otra."""
    def pesos(src, tgt):
        i = np.clip(np.searchsorted(src, tgt) - 1, 0, len(src) - 2)
        w = np.clip((tgt - src[i]) / (src[i + 1] - src[i]), 0.0, 1.0)
        return i, w.astype(np.float32)
    iy, wy = pesos(gy, ty)
    ix, wx = pesos(gx, tx)
    a = f[iy] * (1 - wy)[:, None] + f[iy + 1] * wy[:, None]
    return (a[:, ix] * (1 - wx)[None, :] + a[:, ix + 1] * wx[None, :]).astype(np.float32)


def _vecinos(f):
    p = np.pad(f, 1, mode='edge')
    return np.stack([p[dy:dy + f.shape[0], dx:dx + f.shape[1]]
                     for dy in range(3) for dx in range(3) if (dy, dx) != (1, 1)])


def _mediana_local(u, v, umbral=3.0, eps=0.1):
    """Test de mediana normalizada (Westerweel & Scarano 2005), el de PIVlab."""
    malo = np.zeros(u.shape, bool)
    for f in (u, v):
        nb = _vecinos(f)
        med = np.median(nb, axis=0)
        res = np.median(np.abs(nb - med), axis=0)
        malo |= np.abs(f - med) / (res + eps) > umbral
    return malo


def _filtro_desvio(u, v, malo, n=8.0):
    """Descarta vectores a mas de n desvios estandar del promedio del cuadro."""
    ok = ~malo
    if ok.sum() < 3:
        return malo
    out = malo.copy()
    for f in (u, v):
        m, s = f[ok].mean(), f[ok].std()
        out |= np.abs(f - m) > n * s + 1e-6
    return out


def _rellenar(u, v, malo, dentro):
    """Reemplaza los vectores descartados por el promedio de sus vecinos buenos,
    repitiendo hasta cubrir huecos de mas de un nodo."""
    u, v = u.copy(), v.copy()
    faltan = malo & dentro
    bueno = dentro & ~malo
    for _ in range(50):
        if not faltan.any():
            break
        for f in (u, v):
            g = np.where(bueno, f, 0.0)
            s = _vecinos(g).sum(axis=0)
            c = _vecinos(bueno.astype(np.float32)).sum(axis=0)
            nuevo = faltan & (c > 0)
            f[nuevo] = s[nuevo] / c[nuevo]
        llenos = faltan & (_vecinos(bueno.astype(np.float32)).sum(axis=0) > 0)
        bueno = bueno | llenos
        faltan = faltan & ~llenos
    return u, v


def _suavizar(f):
    p = np.pad(f, 1, mode='edge')
    a = p[:-2] + 2 * p[1:-1] + p[2:]
    return ((a[:, :-2] + 2 * a[:, 1:-1] + a[:, 2:]) / 16.0).astype(np.float32)


# ------------------------------------------------------------------- el motor

class PIV(object):
    """ventanas: tamanos de ventana por pasada, de grande a chica.
    roi: (cx, cy, radio) en px, o None para toda la imagen. Solo se reportan
    nodos cuya ventana final queda entera adentro del circulo."""

    def __init__(self, shape, ventanas=(128, 64, 32), roi=None):
        self.h, self.w = int(shape[0]), int(shape[1])
        self.ventanas = [int(v) for v in ventanas]
        win = self.ventanas[-1]
        self.paso = win // 2
        self.gy = _grilla(self.h, win, self.paso)
        self.gx = _grilla(self.w, win, self.paso)
        self.X, self.Y = np.meshgrid(self.gx, self.gy)
        if roi and roi[2] > 0:
            self.dentro = np.hypot(self.X - roi[0], self.Y - roi[1]) <= roi[2] - win / 2.0
        else:
            self.dentro = np.ones(self.X.shape, bool)
        px = np.arange(self.w, dtype=np.float32)
        py = np.arange(self.h, dtype=np.float32)
        self._px, self._py = px, py
        self._PX, self._PY = np.meshgrid(px, py)

    def _correlacionar(self, A, B, win):
        step = win // 2
        ba = np.array(_bloques(A, win, step), np.float32)
        bb = np.array(_bloques(B, win, step), np.float32)
        ba -= ba.mean(axis=(-2, -1), keepdims=True)
        bb -= bb.mean(axis=(-2, -1), keepdims=True)
        energia = np.sqrt((ba * ba).sum(axis=(-2, -1)) * (bb * bb).sum(axis=(-2, -1)))

        c = np.fft.irfft2(np.fft.rfft2(ba).conj() * np.fft.rfft2(bb), s=(win, win))
        # Solo interesa |desplazamiento| <= win/4 (regla del cuarto). En la
        # salida sin centrar de la FFT esa zona son las cuatro esquinas.
        lim = win // 4
        k = np.r_[win - lim:win, 0:lim + 1]
        sub = c[..., k, :][..., :, k].astype(np.float32)
        # La FFT correlaciona ciclicamente: en un corrimiento d solo se superponen
        # (win - |d|) pixeles por eje. Sin compensarlo, el pico se sesga hacia cero.
        d = np.arange(-lim, lim + 1)
        tri = win / (win - np.abs(d)).astype(np.float32)
        sub *= np.outer(tri, tri)
        sub /= np.maximum(energia, 1e-12)[..., None, None]

        n = 2 * lim + 1
        plano = sub.reshape(sub.shape[0], sub.shape[1], -1)
        idx = plano.argmax(-1)
        iy, ix = idx // n, idx % n
        borde = (iy == 0) | (iy == n - 1) | (ix == 0) | (ix == n - 1)
        iyc, ixc = np.clip(iy, 1, n - 2), np.clip(ix, 1, n - 2)
        I, J = np.ogrid[:sub.shape[0], :sub.shape[1]]

        def gauss3(a0, a1, a2):
            l0, l1, l2 = (np.log(np.maximum(a, 1e-9)) for a in (a0, a1, a2))
            den = 2.0 * (l0 - 2.0 * l1 + l2)
            return np.clip(np.where(np.abs(den) > 1e-12, (l0 - l2) / np.where(den == 0, 1, den), 0), -1, 1)

        dy = gauss3(sub[I, J, iyc - 1, ixc], sub[I, J, iyc, ixc], sub[I, J, iyc + 1, ixc])
        dx = gauss3(sub[I, J, iyc, ixc - 1], sub[I, J, iyc, ixc], sub[I, J, iyc, ixc + 1])
        du = (ix - lim + np.where(borde, 0, dx)).astype(np.float32)
        dv = (iy - lim + np.where(borde, 0, dy)).astype(np.float32)
        # Ventana sin destellos: correlacion nula, el pico no significa nada.
        vacia = energia < 0.02 * max(float(np.median(energia)), 1e-12)
        return du, dv, vacia

    def par(self, a, b):
        """Desplazamiento de b respecto de a, en px, en la grilla final.

        Devuelve (u, v, interpolado): interpolado marca los vectores que la
        validacion descarto y se reemplazaron por sus vecinos.
        """
        u = v = None
        gy0 = gx0 = None
        for i, win in enumerate(self.ventanas):
            step = win // 2
            gy, gx = _grilla(self.h, win, step), _grilla(self.w, win, step)
            if u is None:
                A, B = a, b
            else:
                # Deformacion simetrica: cada imagen se mueve medio desplazamiento,
                # en sentidos opuestos. La interpolacion NO es un detalle: medido
                # sobre un corrimiento conocido de 0.5 px, la bicubica de OpenCV
                # deja un sesgo de -0.057 px y la bilineal +0.029 px; Lanczos
                # (8x8) deja 0.001 px. PIVlab usa spline por el mismo motivo.
                uf = _interp_campo(u, gy0, gx0, self._py, self._px) * 0.5
                vf = _interp_campo(v, gy0, gx0, self._py, self._px) * 0.5
                A = cv2.remap(a, self._PX - uf, self._PY - vf, cv2.INTER_LANCZOS4,
                              borderMode=cv2.BORDER_REPLICATE)
                B = cv2.remap(b, self._PX + uf, self._PY + vf, cv2.INTER_LANCZOS4,
                              borderMode=cv2.BORDER_REPLICATE)
                u = _interp_campo(u, gy0, gx0, gy, gx)
                v = _interp_campo(v, gy0, gx0, gy, gx)
            du, dv, vacia = self._correlacionar(A, B, win)
            u = du if u is None else u + du
            v = dv if v is None else v + dv
            ultima = i == len(self.ventanas) - 1
            malo = _mediana_local(u, v) | vacia
            if ultima:
                # El promedio y el desvio del filtro se toman solo adentro del ROI:
                # afuera no hay glitter y los vectores de ahi no son del flujo.
                malo = _filtro_desvio(u, v, malo | ~self.dentro) & self.dentro
                u, v = _rellenar(u, v, malo, self.dentro)
                return u, v, malo
            todo = np.ones(u.shape, bool)
            u, v = _rellenar(u, v, malo, todo)
            u, v = _suavizar(u), _suavizar(v)
            gy0, gx0 = gy, gx
