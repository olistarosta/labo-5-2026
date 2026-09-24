# -*- coding: utf-8 -*-
"""
PIV por FFT con deformacion de ventana, multipasada. Mismo metodo que PIVlab
("Multipass FFT window deformation", https://www.pivlab.de/wiki/1-quickstart/):

  1. Preprocesamiento: CLAHE (ecualizacion adaptativa de contraste) o lineal
     (la imagen tal cual, sin cuantizar ni recortar).
  2. Correlacion cruzada por FFT en ventanas de interrogacion con 50 % de
     solapamiento, en varias pasadas de ventana decreciente (128 -> 64 -> 32).
     El desplazamiento de cada pasada deforma las imagenes para la siguiente,
     que solo tiene que medir el residuo.
  3. Pico sub-pixel por ajuste gaussiano de 3 puntos en cada eje.
  4. Validacion: ventana sin textura, relacion entre el pico y el segundo pico
     (calidad de la correlacion), filtro de mediana local (umbral 3) y filtro de
     desvio estandar (n = 8).

Los vectores que no pasan la validacion se DESCARTAN (quedan NaN): donde no hay
informacion no se inventa un vector. `rellenar=True` recupera el comportamiento
viejo (reemplazarlos por el promedio de sus vecinos). Adentro del multipasada si
se rellenan, porque la pasada siguiente necesita un predictor en toda la grilla;
eso no sale al resultado, que solo reporta lo medido.

Trabaja en DESPLAZAMIENTOS en pixeles entre dos imagenes. Pasar a velocidad
(escala en mm/px y paso temporal) es trabajo de quien llama.

Necesita numpy y cv2 (TouchDesigner trae los dos).
"""

from concurrent.futures import ThreadPoolExecutor

import cv2
import numpy as np


# ------------------------------------------------------------ preprocesamiento

def preprocesar(img, modo='clahe', ventana_clahe=64):
    """Imagen (gris 0..1 o uint8) -> float32 lista para correlacionar.

    'clahe': como PIVlab. Escala al percentil 99.9, cuantiza a 8 bits y ecualiza
    el contraste por zonas de `ventana_clahe` px: iguala destellos en zonas bien y
    mal iluminadas, que es lo que mas estabiliza la correlacion.
    'lineal': la imagen tal cual (32 bits). No recorta los destellos mas brillantes
    ni cuantiza; la correlacion usa toda la informacion de la proyeccion de color.
    """
    a = np.asarray(img)
    if a.ndim == 3:
        a = a[..., 0]
    if modo == 'lineal':
        return np.ascontiguousarray(a, np.float32)
    if a.dtype != np.uint8:
        tope = max(float(np.percentile(a, 99.9)), 1e-6)
        a = np.rint(np.clip(a / tope, 0.0, 1.0) * 255.0).astype(np.uint8)
    h, w = a.shape
    zonas = (max(1, round(w / ventana_clahe)), max(1, round(h / ventana_clahe)))
    a = cv2.createCLAHE(clipLimit=2.0, tileGridSize=zonas).apply(a)
    return a.astype(np.float32) / 255.0


# -------------------------------------------------------------- utilidades

_POOL = [None]


def _hilos(n=8):
    """Las FFT se reparten en hilos: pocketfft (numpy) suelta el GIL mientras calcula."""
    if _POOL[0] is None:
        _POOL[0] = ThreadPoolExecutor(max_workers=n)
    return _POOL[0]


def _fft_correlacion(ba, bb, win, nh=8):
    """irfft2(conj(FFT(a)) * FFT(b)) de todas las ventanas, en nh hilos."""
    def parte(sl):
        return np.fft.irfft2(np.fft.rfft2(ba[sl]).conj() * np.fft.rfft2(bb[sl]), s=(win, win))
    n = ba.shape[0]
    if n < 64 or nh <= 1:
        return parte(slice(None))
    cortes = [slice(i * n // nh, (i + 1) * n // nh) for i in range(nh)]
    return np.concatenate(list(_hilos(nh).map(parte, cortes)), axis=0)


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
    nodos cuya ventana final queda entera adentro del circulo, y solo se
    correlacionan las ventanas que tocan el ROI (afuera no hay nada que medir).
    pico_min: relacion minima entre el pico de correlacion y el segundo pico.
    rellenar: True = los vectores descartados se reemplazan por sus vecinos."""

    def __init__(self, shape, ventanas=(128, 64, 32), roi=None, pico_min=1.0,
                 rellenar=False, hilos=8):
        self.h, self.w = int(shape[0]), int(shape[1])
        self.ventanas = [int(v) for v in ventanas]
        self.pico_min = float(pico_min)
        self.rellenar = bool(rellenar)
        self.hilos = int(hilos)
        self.roi = roi if (roi and roi[2] > 0) else None
        win = self.ventanas[-1]
        self.paso = win // 2
        self.gy = _grilla(self.h, win, self.paso)
        self.gx = _grilla(self.w, win, self.paso)
        self.X, self.Y = np.meshgrid(self.gx, self.gy)
        if self.roi:
            self.dentro = np.hypot(self.X - roi[0], self.Y - roi[1]) <= roi[2] - win / 2.0
        else:
            self.dentro = np.ones(self.X.shape, bool)
        px = np.arange(self.w, dtype=np.float32)
        py = np.arange(self.h, dtype=np.float32)
        self._px, self._py = px, py
        self._PX, self._PY = np.meshgrid(px, py)

    def cerca_del_roi(self, win):
        """Ventanas que vale la pena correlacionar: las que tocan el ROI, con un
        margen de dos ventanas (4 nodos) para que el test de mediana, el relleno
        y el suavizado de las pasadas intermedias vean lo mismo que si se
        correlacionara toda la imagen. Medido: adentro del ROI los vectores dan
        identicos (< 1e-3 px) y afuera no hay flujo que medir."""
        step = win // 2
        gy, gx = _grilla(self.h, win, step), _grilla(self.w, win, step)
        X, Y = np.meshgrid(gx, gy)
        if not self.roi:
            return np.ones(X.shape, bool)
        return np.hypot(X - self.roi[0], Y - self.roi[1]) <= self.roi[2] + 2 * win

    def _correlacionar(self, A, B, win, sel=None):
        """Devuelve (du, dv, malo, calidad) en la grilla de esa ventana.
        malo = ventana sin textura o pico de correlacion poco definido."""
        step = win // 2
        va, vb = _bloques(A, win, step), _bloques(B, win, step)
        ny, nx = va.shape[:2]
        if sel is None:
            sel = np.ones((ny, nx), bool)
        iy, ix = np.nonzero(sel)
        gy, gx = _grilla(self.h, win, step), _grilla(self.w, win, step)
        if self.roi:
            nucleo = np.hypot(gx[ix] - self.roi[0], gy[iy] - self.roi[1]) <= self.roi[2]
        else:
            nucleo = np.ones(len(iy), bool)
        ba = np.array(va[iy, ix], np.float32)          # solo las ventanas elegidas
        bb = np.array(vb[iy, ix], np.float32)
        ba -= ba.mean(axis=(-2, -1), keepdims=True)
        bb -= bb.mean(axis=(-2, -1), keepdims=True)
        energia = np.sqrt((ba * ba).sum(axis=(-2, -1)) * (bb * bb).sum(axis=(-2, -1)))

        c = _fft_correlacion(ba, bb, win, self.hilos)
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
        plano = sub.reshape(sub.shape[0], -1)
        idx = plano.argmax(-1)
        fil = np.arange(plano.shape[0])
        alto = plano[fil, idx]
        iy0, ix0 = idx // n, idx % n
        # Segundo pico: el maximo afuera de un entorno de 5x5 del primero. Si los
        # dos son parecidos, la correlacion no eligio nada y el vector no vale.
        off = np.arange(-2, 3)
        yy = np.clip(iy0[:, None] + off, 0, n - 1)
        xx = np.clip(ix0[:, None] + off, 0, n - 1)
        tapa = plano.copy()
        np.put_along_axis(tapa, (yy[:, :, None] * n + xx[:, None, :]).reshape(len(fil), -1),
                          -np.inf, axis=1)
        segundo = tapa.max(-1)
        calidad = np.where(alto > 0, alto / np.maximum(segundo, 1e-12), 0.0).astype(np.float32)

        borde = (iy0 == 0) | (iy0 == n - 1) | (ix0 == 0) | (ix0 == n - 1)
        iyc, ixc = np.clip(iy0, 1, n - 2), np.clip(ix0, 1, n - 2)

        def gauss3(a0, a1, a2):
            l0, l1, l2 = (np.log(np.maximum(a, 1e-9)) for a in (a0, a1, a2))
            den = 2.0 * (l0 - 2.0 * l1 + l2)
            return np.clip(np.where(np.abs(den) > 1e-12, (l0 - l2) / np.where(den == 0, 1, den), 0), -1, 1)

        s = sub
        dy = gauss3(s[fil, iyc - 1, ixc], s[fil, iyc, ixc], s[fil, iyc + 1, ixc])
        dx = gauss3(s[fil, iyc, ixc - 1], s[fil, iyc, ixc], s[fil, iyc, ixc + 1])
        du_s = (ix0 - lim + np.where(borde, 0, dx)).astype(np.float32)
        dv_s = (iy0 - lim + np.where(borde, 0, dy)).astype(np.float32)
        # Ventana sin destellos: la correlacion es ruido, el pico no significa nada.
        # La referencia es la mediana de las ventanas de ADENTRO del ROI, para que
        # no dependa de cuanta imagen vacia se haya correlacionado alrededor.
        tipica = float(np.median(energia[nucleo])) if nucleo.any() else float(np.median(energia))
        vacia = energia < 0.02 * max(tipica, 1e-12)
        malo_s = vacia | borde

        du = np.zeros((ny, nx), np.float32)
        dv = np.zeros((ny, nx), np.float32)
        cal = np.zeros((ny, nx), np.float32)
        malo = np.ones((ny, nx), bool)                 # lo que no se calculo es "sin dato"
        du[iy, ix], dv[iy, ix], cal[iy, ix], malo[iy, ix] = du_s, dv_s, calidad, malo_s
        return du, dv, malo, cal

    def par(self, a, b):
        """Desplazamiento de b respecto de a, en px, en la grilla final.

        Devuelve (u, v, descartado). Los vectores descartados quedan NaN, salvo
        que se haya pedido rellenar (ahi valen el promedio de sus vecinos buenos).
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
            du, dv, sin_dato, calidad = self._correlacionar(A, B, win, self.cerca_del_roi(win))
            u = du if u is None else u + du
            v = dv if v is None else v + dv
            malo = _mediana_local(u, v) | sin_dato
            if i == len(self.ventanas) - 1:
                # La calidad del pico solo decide en la ultima pasada, que es la
                # que se reporta. Aplicarla antes empeora: el hueco se rellena con
                # los vecinos y ese predictor peor arrastra a las pasadas que siguen.
                malo |= calidad < self.pico_min
                # El promedio y el desvio del filtro se toman solo adentro del ROI:
                # afuera no hay glitter y los vectores de ahi no son del flujo.
                malo = _filtro_desvio(u, v, malo | ~self.dentro) & self.dentro
                if self.rellenar:
                    u, v = _rellenar(u, v, malo, self.dentro)
                else:
                    u, v = u.copy(), v.copy()
                    u[malo], v[malo] = np.nan, np.nan
                return u, v, malo
            # Entre pasadas si se rellena: la que sigue necesita un predictor en
            # toda la grilla para deformar. Eso no sale al resultado.
            u, v = _rellenar(u, v, malo, np.ones(u.shape, bool))
            u, v = _suavizar(u), _suavizar(v)
            gy0, gx0 = gy, gx
