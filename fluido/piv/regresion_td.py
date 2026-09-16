# -*- coding: utf-8 -*-
"""
Regresion: el analisis de /analizar en TouchDesigner contra el motor de referencia.

    python regresion_td.py <carpeta_grabacion> <carpeta_del_analisis>/parametros.json [pares]

Para cada camara del JSON:
1. Vector a vector: los mismos cuadros por piv_core.PIV (la referencia validada)
   contra el CSV que escribio TD. Prueba que la red -- proyeccion de color, CLAHE,
   deformacion GLSL, muestreos de grilla, feedback, orden de cocinado -- hace la misma cuenta.
2. Si la camara grabo la simulacion:
   - w(r) del CSV contra el vortice inyectado (medido desde el centro seguido y en el
     sistema del vortice, restando la velocidad del centro);
   - el centro seguido (_centro.csv) contra el centro real que guardo la grabacion.
"""
import json
import os
import sys

import cv2
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from piv_core import PIV, preprocesar      # noqa: E402


def direccion(fondo, objetivo, otro, separable_min=0.25):
    """Igual que contraste_color.direccion en TD: 1 en objetivo, 0 en fondo y en otro."""
    d = objetivo - fondo
    if d @ d < 1e-6:
        d = np.ones(3)
    s = otro - fondo
    w = d.copy()
    if s @ s > 1e-6:
        perp = d - (d @ s) / (s @ s) * s
        if np.linalg.norm(perp) >= separable_min * np.linalg.norm(d):
            w = perp
    return w / (d @ w)


def imagen_piv(ruta, colores, clahe):
    """Cuadro -> lo que entra al bucle en TD: proyeccion del glitter, escala p99.9, CLAHE."""
    im = cv2.imread(ruta, cv2.IMREAD_UNCHANGED)
    if im.ndim == 2:
        rgb = np.repeat(im[..., None], 3, axis=2)
    else:
        rgb = cv2.cvtColor(im[..., :3], cv2.COLOR_BGR2RGB)
    rgb = np.flipud(rgb).astype(np.float32) / 255.0      # TD: fila 0 ABAJO
    f, g, p = (np.array(colores[k]) for k in ('fondo', 'glitter', 'particula'))
    v = np.maximum((rgb - f) @ direccion(f, g, p), 0.0)
    tope = max(float(np.percentile(v, 99.9)), 1e-6)
    gris = np.rint(np.clip(v / tope, 0, 1) * 255).astype(np.uint8)
    return preprocesar(gris, clahe=clahe)


grab, ruta_json = sys.argv[1], sys.argv[2]
npares = int(sys.argv[3]) if len(sys.argv) > 3 else 20
p = json.load(open(ruta_json, encoding='utf-8'))
salida = os.path.dirname(ruta_json)
ventanas = p['piv']['ventanas_px']
clahe = p['piv']['preprocesamiento'] == 'CLAHE'
BRILLO = {'fondo': [0, 0, 0], 'glitter': [1, 1, 1], 'particula': [0, 0, 0]}   # analisis sin colores

for L, cam in p['camaras'].items():
    sub = os.path.join(grab, 'camara_' + L)
    if not os.path.isdir(sub):
        sub = grab
    cal = cam['calibracion']
    s = cal['mm_por_px'] or 1.0
    dt = cam['tiempo']['dt_por_cuadro_s']
    roi = (cal['centro_vortice_px'][0], cal['centro_vortice_px'][1], cal['radio_roi_px'])
    colores = cam.get('colores', BRILLO)
    tiempos = np.loadtxt(os.path.join(sub, 'tiempos.csv'), delimiter=',', skiprows=1, ndmin=2)
    sep = int(p['piv'].get('separacion_cuadros', 1))       # pares (k, k + sep)
    nint = np.maximum(np.rint(np.diff(tiempos[:, 1]) / dt), 1).astype(int)
    # Mismo paso temporal que TD: diferencia de horas (simulación) o dt fijo por los intervalos
    # entre los dos cuadros (cámara).
    cum = np.r_[0, np.cumsum(nint)]
    if cam['tiempo'].get('reloj_exacto'):
        dt_par = tiempos[sep:, 1] - tiempos[:-sep, 1]
    else:
        dt_par = (cum[sep:] - cum[:-sep]) * dt
    csv = np.loadtxt(os.path.join(salida, cam['csv']), delimiter=',', skiprows=1)
    cuadros = sorted(os.listdir(os.path.join(sub, 'cuadros')))
    leer = lambda k: imagen_piv(os.path.join(sub, 'cuadros', cuadros[k]), colores, clahe)

    h, w = leer(0).shape
    ref = PIV((h, w), ventanas=ventanas, roi=roi)
    difs, marcas = [], []
    for k in range(min(npares, cam['pares_analizados'])):
        u, v, interp = ref.par(leer(k), leer(k + sep))
        m = ref.dentro
        filas = csv[csv[:, 0] == k]
        dtp = dt_par[k]
        du_td, dv_td = filas[:, 6] / s * dtp, filas[:, 7] / s * dtp
        it_td = filas[:, 8] > 0.5
        ok = ~interp[m] & ~it_td
        difs.append(np.hypot(du_td - u[m], dv_td - v[m])[ok])
        marcas.append(np.mean(interp[m] != it_td))
    dd = np.concatenate(difs)
    print('== camara %s: TD contra piv_core, %d pares, %d vectores ==' % (L, len(difs), dd.size))
    print('   |d_TD - d_ref| px: mediana %.5f  p95 %.5f  maxima %.5f   marca de interpolado distinta: %.2f %%'
          % (np.median(dd), np.percentile(dd, 95), dd.max(), 100 * np.mean(marcas)))

    info = p.get('grabacion_info', {})
    sim = info.get('camaras', {}).get(L, info).get('simulacion')
    if not sim:
        continue

    # Centro de cada par: el seguido si hay _centro.csv; si no, el marcado.
    k_par = csv[:, 0].astype(int)
    if 'csv_centro' in cam:
        cen = np.loadtxt(os.path.join(salida, cam['csv_centro']), delimiter=',', skiprows=1, ndmin=2)
        cx_par, cy_par, t_par = cen[:, 2], cen[:, 3], cen[:, 1]
        uc, vc = np.gradient(cx_par, t_par), np.gradient(cy_par, t_par)
        cx, cy = cx_par[k_par], cy_par[k_par]
        ucv, vcv = uc[k_par], vc[k_par]
        if tiempos.shape[1] >= 4:
            n = len(cen)
            real_x = 0.5 * (tiempos[:n, 2] + tiempos[sep:n + sep, 2])
            real_y = 0.5 * (tiempos[:n, 3] + tiempos[sep:n + sep, 3])
            e = np.hypot(cx_par - real_x, cy_par - real_y)
            seg = cen[:, 6] > 0.5
            print('   centro seguido contra el real: %d de %d pares con particula, error mediana %.3f px, '
                  'p95 %.3f px, max %.3f px' % (seg.sum(), n, np.median(e[seg]), np.percentile(e[seg], 95), e[seg].max()))
    else:
        cx, cy = cal['centro_vortice_px']
        ucv = vcv = 0.0

    w0, rc = sim['omega_nucleo_rad_s'], sim['radio_nucleo_px']
    x, y = csv[:, 2] - cx, csv[:, 3] - cy
    uu, vv = csv[:, 6] / s - ucv, csv[:, 7] / s - vcv
    r = np.hypot(x, y)
    om = (-uu * y + vv * x) / r ** 2
    vr = (uu * x + vv * y) / r
    print('   w(r) contra el Rankine inyectado (w0=%.2f, Rc=%.0f px), %d pares:' %
          (w0, rc, len(np.unique(csv[:, 0]))))
    for lo, hi in ((20, 60), (60, 110), (110, 160), (160, 260)):
        k = (r >= lo) & (r < hi)
        wt = np.where(r[k] <= rc, w0, w0 * rc * rc / r[k] ** 2)
        dts = dt * sep                                      # la cuerda depende del intervalo del par
        real = np.mean(2 * np.sin(wt * dts / 2) / dts)
        print('     r %3d-%3d px   w %.4f  real %.4f  error %+6.2f %%   v_r %+5.2f px/s' %
              (lo, hi, om[k].mean(), real, 100 * (om[k].mean() / real - 1), vr[k].mean()))
