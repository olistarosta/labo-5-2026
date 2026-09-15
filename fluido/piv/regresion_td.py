# -*- coding: utf-8 -*-
"""
Regresion: el analisis de /analizar en TouchDesigner contra el motor de referencia.

    <python de TouchDesigner> regresion_td.py <carpeta_grabacion> <..._parametros.json> [pares]

Para cada camara del JSON:
1. Vector a vector: los mismos cuadros por piv_core.PIV (la referencia validada)
   contra el CSV que escribio TD. Prueba que la red -- deformacion GLSL, muestreos
   de grilla, feedback, orden de cocinado -- hace exactamente la misma cuenta.
2. Si la camara grabo la simulacion: w(r) del CSV contra el vortice inyectado.
"""
import json
import os
import sys

import cv2
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from piv_core import PIV, preprocesar      # noqa: E402

grab, ruta_json = sys.argv[1], sys.argv[2]
npares = int(sys.argv[3]) if len(sys.argv) > 3 else 20
p = json.load(open(ruta_json, encoding='utf-8'))
salida = os.path.dirname(ruta_json)
ventanas = p['piv']['ventanas_px']
clahe = p['piv']['preprocesamiento'] == 'CLAHE'

for L, cam in p['camaras'].items():
    sub = os.path.join(grab, 'camara_' + L)
    if not os.path.isdir(sub):
        sub = grab
    cal = cam['calibracion']
    s = cal['mm_por_px'] or 1.0
    dt = cam['tiempo']['dt_por_cuadro_s']
    roi = (cal['centro_vortice_px'][0], cal['centro_vortice_px'][1], cal['radio_roi_px'])
    t = np.loadtxt(os.path.join(sub, 'tiempos.csv'), delimiter=',', skiprows=1, ndmin=2)[:, 1]
    nint = np.maximum(np.rint(np.diff(t) / dt), 1).astype(int)
    csv = np.loadtxt(os.path.join(salida, cam['csv']), delimiter=',', skiprows=1)
    cuadros = sorted(os.listdir(os.path.join(sub, 'cuadros')))

    def leer(k):
        im = cv2.imread(os.path.join(sub, 'cuadros', cuadros[k]), cv2.IMREAD_UNCHANGED)
        gris = im[..., 2] if im.ndim == 3 else im
        # TD entrega las imagenes con la fila 0 ABAJO; el archivo la tiene arriba.
        return preprocesar(np.flipud(gris), clahe=clahe)

    h, w = leer(0).shape
    ref = PIV((h, w), ventanas=ventanas, roi=roi)
    difs, marcas = [], []
    for k in range(min(npares, cam['pares_analizados'])):
        u, v, interp = ref.par(leer(k), leer(k + 1))
        m = ref.dentro
        filas = csv[csv[:, 0] == k]
        dtp = nint[k] * dt
        du_td, dv_td = filas[:, 6] / s * dtp, filas[:, 7] / s * dtp
        it_td = filas[:, 8] > 0.5
        ok = ~interp[m] & ~it_td
        difs.append(np.hypot(du_td - u[m], dv_td - v[m])[ok])
        marcas.append(np.mean(interp[m] != it_td))
    dd = np.concatenate(difs)
    print('== camara %s: TD contra piv_core, %d pares, %d vectores ==' % (L, len(difs), dd.size))
    print('   |d_TD - d_ref| px: mediana %.5f  p95 %.5f  maxima %.5f   marca de interpolado distinta: %.2f %%'
          % (np.median(dd), np.percentile(dd, 95), dd.max(), 100 * np.mean(marcas)))

    sim = p.get('grabacion_info', {}).get('camaras', {}).get(L, {}).get('simulacion')
    if sim:
        w0, rc = sim['omega_nucleo_rad_s'], sim['radio_nucleo_px']
        cx, cy = sim['centro_px']
        x, y = csv[:, 2] - cx, csv[:, 3] - cy
        uu, vv = csv[:, 6] / s, csv[:, 7] / s
        r = np.hypot(x, y)
        om = (-uu * y + vv * x) / r ** 2
        vr = (uu * x + vv * y) / r
        print('   w(r) contra el Rankine inyectado (w0=%.2f, Rc=%.0f px), %d pares:' %
              (w0, rc, len(np.unique(csv[:, 0]))))
        for lo, hi in ((20, 60), (60, 110), (110, 160), (160, 260)):
            k = (r >= lo) & (r < hi)
            wt = np.where(r[k] <= rc, w0, w0 * rc * rc / r[k] ** 2)
            real = np.mean(2 * np.sin(wt * dt / 2) / dt)
            print('     r %3d-%3d px   w %.4f  real %.4f  error %+6.2f %%   v_r %+5.2f px/s' %
                  (lo, hi, om[k].mean(), real, 100 * (om[k].mean() / real - 1), vr[k].mean()))
