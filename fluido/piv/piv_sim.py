# -*- coding: utf-8 -*-
"""
Simulador de glitter en un vortice, con campo de velocidades CONOCIDO.

Existe para una sola cosa: poder decir cuanto se equivoca el PIV. Todo lo que
mide el motor se contrasta contra el campo analitico que este modulo inyecta.

Vortice de Rankine
------------------
    w(r) = w0                  si r <= Rc     (nucleo en rotacion rigida)
    w(r) = w0 * Rc^2 / r^2     si r >  Rc     (vortice libre, v_th ~ 1/r)

La velocidad angular w(r) no depende del angulo, asi que la trayectoria exacta
de cada destello es una rotacion pura de angulo w(r)*dt alrededor del centro,
con r constante. No hay integrador numerico ni error de paso: la posicion del
cuadro siguiente es exacta. Con deriva radial activada se integra con RK4.

Lo que el PIV realmente mide
----------------------------
No es la velocidad instantanea sino la CUERDA: (p(t+dt) - p(t)) / dt. En un
flujo curvo las dos difieren, poco pero de forma sistematica: para una rotacion
de angulo th = w*dt la cuerda tiene modulo v*sin(th/2)/(th/2) y esta girada th/2.
A 30 fps y w = 2 rad/s eso es 0.015 % en modulo y 1.9 grados en direccion.
`chord_velocity()` devuelve esa cantidad exacta, que es el patron contra el cual
hay que juzgar al PIV; `true_velocity()` da la instantanea, para poder reportar
las dos cosas por separado y no confundir sesgo del metodo con fisica.
"""

import numpy as np


class GlitterSim:
    """Siembra de destellos advectados por un vortice de Rankine."""

    def __init__(self, h=480, w=640, n=1200, omega0=2.2, rcore=85.0,
                 center=None, inflow=0.0, seed=0, sigma=1.2,
                 twinkle=0.8, dropout=0.05, noise=0.012, vignette=0.35,
                 rmax=None, profile='rankine'):
        self.h, self.w = int(h), int(w)
        self.omega0 = float(omega0)
        self.rcore = float(rcore)
        self.cx, self.cy = center if center else (w * 0.5, h * 0.5)
        self.inflow = float(inflow)
        self.profile = profile
        self.sigma = float(sigma)
        self.twinkle = float(twinkle)
        self.dropout = float(dropout)
        self.noise = float(noise)
        self.rng = np.random.default_rng(seed)
        self.rmax = float(rmax) if rmax else min(self.cx, self.cy, w - self.cx, h - self.cy) * 0.98

        # Sembrado uniforme EN AREA (r = rmax*sqrt(u)), no uniforme en radio:
        # si no, se amontonan en el centro y las ventanas de afuera quedan vacias.
        self.r = self.rmax * np.sqrt(self.rng.random(n))
        self.th = self.rng.uniform(0.0, 2.0 * np.pi, n)
        self.bright = 0.55 + 0.45 * self.rng.random(n)
        self.size = self.sigma * (0.8 + 0.4 * self.rng.random(n))

        # Iluminacion despareja: el pasaaltos del motor tiene que comersela.
        yy, xx = np.mgrid[0:self.h, 0:self.w].astype(np.float32)
        rr = np.hypot(xx - self.w * 0.5, yy - self.h * 0.5) / (0.5 * np.hypot(self.h, self.w))
        self.illum = (1.0 - vignette * rr ** 2).astype(np.float32)

        self._stamp_r = int(np.ceil(3.0 * self.sigma * 1.2))

    # -- campo -------------------------------------------------------------

    def omega(self, r):
        r = np.asarray(r, np.float64)
        if self.profile == 'solid':
            # Rotacion rigida: campo de velocidad LINEAL en x, y. La deformacion
            # de ventana lo representa exactamente, asi que el error que quede
            # es del metodo y no de la resolucion espacial. Es el patron de
            # referencia para separar una cosa de la otra.
            return np.full(r.shape, self.omega0)
        if self.profile == 'lamboseen':
            # Vortice de Lamb-Oseen: mismo aire que Rankine pero SIN el quiebre
            # en r = rcore, que es donde una ventana de interrogacion promedia
            # dos pendientes distintas y sesga.
            rr = np.maximum(r, 1e-9)
            return (self.omega0 * self.rcore ** 2 / rr ** 2) * (
                1.0 - np.exp(-rr ** 2 / self.rcore ** 2)) / (1.0 - np.exp(-1.0))
        return np.where(r <= self.rcore, self.omega0,
                        self.omega0 * self.rcore ** 2 / np.maximum(r, 1e-9) ** 2)

    def true_velocity(self, x, y):
        """Velocidad instantanea en px/s. y crece hacia arriba."""
        dx = np.asarray(x, np.float64) - self.cx
        dy = np.asarray(y, np.float64) - self.cy
        r = np.hypot(dx, dy)
        w = self.omega(r)
        u = -w * dy + self.inflow * (-dx / np.maximum(r, 1e-9))
        v = w * dx + self.inflow * (-dy / np.maximum(r, 1e-9))
        return u, v

    def chord_velocity(self, x, y, dt, midpoint=True):
        """Lo que el PIV deberia medir: cuerda / dt, exacta para la rotacion.

        midpoint decide A QUE POSICION se le atribuye la cuerda, y no es un
        detalle: cambia el patron en un 0.07 % a 30 fps, que es del orden del
        error propio del motor.

        Con deformacion simetrica (las dos imagenes se deforman media diferencia
        cada una, en sentidos opuestos) la correlacion queda centrada en el
        PUNTO MEDIO del recorrido, asi que el vector medido en el nodo (x, y)
        corresponde a la particula cuyo punto medio cae ahi. Esa cuerda es
        exactamente tangencial y vale 2*r*sin(w*dt/2). Con midpoint=False se
        devuelve la convencion de punto inicial, cuya componente tangencial es
        r*sin(w*dt) -- util solo para comparar contra codigos que no deforman.
        """
        dx = np.asarray(x, np.float64) - self.cx
        dy = np.asarray(y, np.float64) - self.cy
        r = np.hypot(dx, dy)
        if self.inflow == 0.0:
            th = self.omega(r) * dt
            if midpoint:
                # Gira media vuelta para atras y media para adelante.
                c, s = np.cos(th * 0.5), np.sin(th * 0.5)
                ax, ay = c * dx + s * dy, -s * dx + c * dy
                nx, ny = c * dx - s * dy, s * dx + c * dy
                return (nx - ax) / dt, (ny - ay) / dt
            c, s = np.cos(th), np.sin(th)
            nx = c * dx - s * dy
            ny = s * dx + c * dy
        else:
            if midpoint:
                ax, ay = self._rk4(dx, dy, -0.5 * dt)
                nx, ny = self._rk4(dx, dy, 0.5 * dt)
                return (nx - ax) / dt, (ny - ay) / dt
            nx, ny = self._rk4(dx, dy, dt)
        return (nx - dx) / dt, (ny - dy) / dt

    def _rk4(self, dx, dy, dt):
        def f(ax, ay):
            r = np.hypot(ax, ay)
            w = self.omega(r)
            ir = 1.0 / np.maximum(r, 1e-9)
            return -w * ay - self.inflow * ax * ir, w * ax - self.inflow * ay * ir
        k1x, k1y = f(dx, dy)
        k2x, k2y = f(dx + 0.5 * dt * k1x, dy + 0.5 * dt * k1y)
        k3x, k3y = f(dx + 0.5 * dt * k2x, dy + 0.5 * dt * k2y)
        k4x, k4y = f(dx + dt * k3x, dy + dt * k3y)
        return (dx + dt / 6.0 * (k1x + 2 * k2x + 2 * k3x + k4x),
                dy + dt / 6.0 * (k1y + 2 * k2y + 2 * k3y + k4y))

    # -- evolucion ---------------------------------------------------------

    def step(self, dt):
        """Avanza los destellos dt segundos. Exacto para inflow = 0."""
        if self.inflow == 0.0:
            self.th = self.th + self.omega(self.r) * dt
        else:
            dx, dy = self.r * np.cos(self.th), self.r * np.sin(self.th)
            nx, ny = self._rk4(dx, dy, dt)
            self.r, self.th = np.hypot(nx, ny), np.arctan2(ny, nx)
            self._respawn()

    def _respawn(self):
        """Los que se van por el desague o por el borde vuelven a sembrarse,
        uniformes en area, para que la densidad no se desmorone."""
        gone = (self.r < 3.0) | (self.r > self.rmax)
        k = int(gone.sum())
        if k:
            self.r[gone] = self.rmax * np.sqrt(self.rng.random(k))
            self.th[gone] = self.rng.uniform(0.0, 2.0 * np.pi, k)

    def positions(self):
        return self.cx + self.r * np.cos(self.th), self.cy + self.r * np.sin(self.th)

    # -- render ------------------------------------------------------------

    def render(self, x=None, y=None, seed_noise=True):
        """Imagen float32 (h, w) en 0..1."""
        if x is None:
            x, y = self.positions()
        n = len(x)
        amp = self.bright.copy()
        if self.twinkle > 0.0:
            # El glitter no brilla parejo: cada lamina refleja segun como quedo
            # orientada, y eso cambia cuadro a cuadro. Es justamente el ruido
            # que la correlacion de conjunto promedia.
            amp *= (1.0 - self.twinkle) + self.twinkle * self.rng.random(n) ** 2
        if self.dropout > 0.0:
            amp *= (self.rng.random(n) > self.dropout)

        img = np.zeros((self.h, self.w), np.float32)
        R = self._stamp_r
        ox = np.arange(-R, R + 1)
        gx = np.round(x).astype(np.int64)
        gy = np.round(y).astype(np.int64)
        fx = (x - gx)[:, None]
        fy = (y - gy)[:, None]
        s2 = 2.0 * (self.size ** 2)[:, None]

        wx = np.exp(-(ox[None, :] - fx) ** 2 / s2)
        wy = np.exp(-(ox[None, :] - fy) ** 2 / s2)
        stamp = (amp[:, None, None] * wy[:, :, None] * wx[:, None, :]).astype(np.float32)

        iy = (gy[:, None] + ox[None, :])
        ix = (gx[:, None] + ox[None, :])
        ok_y = (iy >= 0) & (iy < self.h)
        ok_x = (ix >= 0) & (ix < self.w)
        np.clip(iy, 0, self.h - 1, out=iy)
        np.clip(ix, 0, self.w - 1, out=ix)
        stamp *= (ok_y[:, :, None] & ok_x[:, None, :])

        np.add.at(img, (iy[:, :, None], ix[:, None, :]), stamp)

        img *= self.illum
        if self.noise > 0.0 and seed_noise:
            img += self.rng.normal(0.0, self.noise, img.shape).astype(np.float32)
        return np.clip(img, 0.0, 1.0)

    def render_shifted(self, dx, dy):
        """Render con todos los destellos corridos (dx, dy) exactos. Para medir
        el sesgo sub-pixel sin que se mezcle con el campo del vortice."""
        x, y = self.positions()
        return self.render(x + dx, y + dy)
