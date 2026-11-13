"""Importar grabaciones hechas con otro programa (ThorCam / ThorImageCAM de Thorlabs, o cualquier
TIFF, secuencia de TIFF o video) a nuestro formato (CONTRATO.md), SIN copiar los cuadros.

    carpeta = importar(ruta, meta_usuario)      # ruta: archivo o carpeta

Crea una carpeta de grabación NUEVA (io.nueva_carpeta_grabacion) con:
  - grabacion.json: lo que se pudo leer de los archivos + lo que venga en `meta_usuario`, y
    `fuente` = {formato, carpeta, archivos, ...} apuntando a los ORIGINALES (nunca se modifican;
    si se mueven, se cambia fuente.carpeta). `tiempo_fuente` dice de dónde salen los tiempos.
    `importacion` lista lo encontrado, lo que falta, los conflictos y los avisos.
  - tiempos.csv: una fila por cuadro (columnas de CONTRATO.md + t_nominal_s).
io.Grabacion lee después los cuadros de la fuente (io.LectorFuente).

Qué guarda ThorCam (manual ITN000493-D02, rev. V 2022: secciones 3.3.3, 3.4.6, 3.5.3-3.5.4,
3.7 y Apéndice D):
  - Record / Timed Series de cámaras científicas (Zelux, Kiralux, Quantalux): TIFF MULTIPÁGINA
    sin compresión, 16 bits por muestra, o AVI Motion-JPEG de 8 bits (24 bits RGB) escalado
    como se ve en pantalla. Color: con Image Type 'Unprocessed' el TIFF es el mosaico Bayer
    crudo; con los tipos RGB (sRGB, Linear sRGB) es RGB procesado.
  - Tags privados por PÁGINA (TSI_*): 32768 bits significativos, 32770/32771 binning x/y,
    32772-32775 ROI (origen x, y; píxeles x, y), 32777 pixel clock [Hz], 32778 taps,
    32779 NÚMERO DE CUADRO, 32780 tipo de imagen de polarización, 32781/32782 TIMESTAMP
    RELATIVO [ns] (64 bits = alto << 32 | bajo). La exposición y la ganancia NO están en la
    lista documentada: hay que cargarlas a mano (se guardan igual todos los tags >= 32768).
  - Nombres: <prefijo><sello de tiempo>; las grabaciones largas se parten en varios archivos
    <base>_0.tif, <base>_1.tif ... (ej. 'image_2022-02-17T14-55-13.447_0.tif', visto en
    github.com/hacknus/TEMPusVoLa, que lee 32781/32782 de cada página).
ThorImageCAM (sucesor): según Thorlabs, MP4, TIFF 'clásico' multipágina u OME-TIFF; no hay
documentación pública de sus tags. Se leen los TSI_* si están y, en OME-TIFF, Plane DeltaT y
ExposureTime del OME-XML (origen del DeltaT no confirmado).
"""
import os
import re
import xml.etree.ElementTree as ET
from datetime import datetime

import numpy as np

from . import io, rutas, tiempo

rutas.preparar_sys_path()

VERSION_IMPORTADOR = 4          # 2: los videos descartan los cuadros repetidos; 3: x.tif + x_1.tif... es una serie;
#                                 4: también si x termina en número (BROWNIANO1.tif + BROWNIANO1_1.tif...)
EXT_TIFF = ('.tif', '.tiff')
EXT_VIDEO = ('.avi', '.mp4', '.mov', '.mkv', '.wmv')

TAGS_TSI = {32768: 'TSI_SIGNIFICANT_BITS_INT', 32770: 'TSI_BIN_X_INT', 32771: 'TSI_BIN_Y_INT',
            32772: 'TSI_ROI_ORIGIN_X_INT', 32773: 'TSI_ROI_ORIGIN_Y_INT',
            32774: 'TSI_ROI_PIXELS_X_INT', 32775: 'TSI_ROI_PIXELS_Y_INT',
            32777: 'TSI_PIXEL_CLOCK_HZ_INT', 32778: 'TSI_NUM_TAPS_INT',
            32779: 'TSI_FRAME_NUMBER_INT', 32780: 'TSI_POLARIZATION_IMAGE_TYPE',
            32781: 'TSI_TIME_STAMP_RELATIVE_NS_HIGH_BITS',
            32782: 'TSI_TIME_STAMP_RELATIVE_NS_LOW_BITS'}
TAG_CUADRO, TAG_TS_ALTO, TAG_TS_BAJO = 32779, 32781, 32782

AVISO_VIDEO = ('VIDEO: 8 bits por canal y compresión con pérdida (ThorCam: Motion-JPEG escalado '
               'como en pantalla). Sirve para mirar o para un tracking grueso; NO para subpíxel '
               'fino, piso de ruido ni fotometría. Sin timestamps ni contador por cuadro.')
AVISO_SRGB = ('Cuadros RGB procesados por el programa (demosaico + balance de blancos; en sRGB '
              'además gamma): no son cuentas crudas del sensor.')


# ---------------------------------------------------------------- archivos
def orden_natural(nombre):
    """Clave de orden 'humano': img_2 < img_10."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', nombre)]


def _prefijo(nombre, vecinos=()):
    """Nombre sin extensión ni número final: agrupa las partes o los cuadros de una misma serie.
    ThorCam nombra las partes 'x.tif', 'x_1.tif', 'x_2.tif'...: el '_' también se saca.
    vecinos: los demás archivos de la carpeta. Si el nombre no termina en '_N' y entre ellos hay
    un 'nombre_N', es la primera parte de una serie de ThorCam y su número final es parte del
    nombre (BROWNIANO1.tif + BROWNIANO1_1.tif...): no se saca."""
    base = os.path.splitext(nombre)[0]
    sin_parte = re.sub(r'_\d+$', '', base)
    if sin_parte != base:
        return sin_parte
    if any(re.fullmatch(re.escape(base) + r'_\d+', os.path.splitext(v)[0]) for v in vecinos):
        return base
    return re.sub(r'\d+$', '', base)


def listar_archivos(ruta, serie=True):
    """Archivo o carpeta → (formato, carpeta, [nombres en orden]).

    formato: 'tiff' | 'video'. Archivo TIFF con serie=True: él y sus hermanos de la misma serie
    (mismo prefijo + número final: las partes _0, _1 ... de una grabación larga de ThorCam o
    una secuencia numerada). Carpeta: todos los TIFF si son de UNA sola serie, si no el único
    video; si hay varias series o varios videos → ValueError (pasá un archivo).
    """
    ruta = os.path.abspath(ruta)
    if os.path.isfile(ruta):
        carpeta, nombre = os.path.split(ruta)
        ext = os.path.splitext(nombre)[1].lower()
        if ext in EXT_VIDEO:
            return 'video', carpeta, [nombre]
        if ext not in EXT_TIFF:
            raise ValueError(f'no sé leer {nombre} (TIFF: {EXT_TIFF}; video: {EXT_VIDEO})')
        if not serie:
            return 'tiff', carpeta, [nombre]
        tifs = [f for f in os.listdir(carpeta) if f.lower().endswith(EXT_TIFF)]
        p = _prefijo(nombre, tifs)
        hermanos = [f for f in tifs if _prefijo(f, tifs) == p]
        return 'tiff', carpeta, sorted(hermanos, key=orden_natural)
    if not os.path.isdir(ruta):
        raise FileNotFoundError(ruta)
    todos = sorted(os.listdir(ruta), key=orden_natural)
    tifs = [f for f in todos if f.lower().endswith(EXT_TIFF)]
    if tifs:
        series = sorted({_prefijo(f, tifs) for f in tifs})
        if len(series) > 1:
            raise ValueError(f'{ruta} tiene varias series de TIFF ({series[:6]}): pasá un archivo')
        return 'tiff', ruta, tifs
    videos = [f for f in todos if f.lower().endswith(EXT_VIDEO)]
    if len(videos) == 1:
        return 'video', ruta, videos
    raise ValueError(f'{ruta}: {len(videos)} videos y ningún TIFF: pasá un archivo')


def hora_de_nombre(nombre):
    """Sello de ThorCam en el nombre ('..._2022-02-17T14-55-13.447_0.tif') → ISO, o None."""
    m = re.search(r'(\d{4}-\d{2}-\d{2})T(\d{2})-(\d{2})-(\d{2})(\.\d+)?', nombre)
    return f'{m.group(1)}T{m.group(2)}:{m.group(3)}:{m.group(4)}{m.group(5) or ""}' if m else None


# ---------------------------------------------------------------- TIFF
def _valor(tag):
    """Valor de un tag de tifffile → int/float/str/lista nativa (bytes → str)."""
    v = tag.value
    if isinstance(v, bytes):
        return v.decode('latin-1').rstrip('\x00')
    if isinstance(v, (tuple, list, np.ndarray)):
        v = np.asarray(v).tolist()
        return v[0] if len(v) == 1 else v
    if hasattr(v, 'name') and hasattr(v, 'value'):     # enums
        return v.name
    return v.item() if isinstance(v, np.generic) else v


def _entero(tags, codigo):
    t = tags.get(codigo)
    if t is None:
        return None
    v = _valor(t)
    return int(v) if isinstance(v, (int, float)) else None


def timestamp_ns(tags):
    """TSI_TIME_STAMP_RELATIVE_NS: (alto << 32) | bajo [ns], o None si falta alguno de los dos."""
    hi, lo = _entero(tags, TAG_TS_ALTO), _entero(tags, TAG_TS_BAJO)
    if hi is None or lo is None:
        return None
    return ((hi & 0xFFFFFFFF) << 32) | (lo & 0xFFFFFFFF)


def tags_privados(tags, max_largo=64):
    """Tags >= 32768 con valores chicos → {código: valor}. Los blobs largos (ImageJ, etc.) no."""
    out = {}
    for t in tags.values():
        if t.code < 32768:
            continue
        v = _valor(t)
        if isinstance(v, list) and len(v) > max_largo or isinstance(v, str) and len(v) > 2000:
            continue
        out[int(t.code)] = v
    return out


def info_pagina(page):
    """Primera página → dict: forma, dtype, muestras, fotométrica, compresión, bits, tags estándar
    (Software, DateTime, Make, Model, ImageDescription recortada) y tags_privados."""
    tg = page.tags

    def std(nombre):
        t = tg.get(nombre)
        return None if t is None else _valor(t)
    desc = std('ImageDescription')
    return dict(forma=list(page.shape), dtype=str(page.dtype), muestras=int(page.samplesperpixel),
                fotometrica=getattr(page.photometric, 'name', str(page.photometric)),
                planar=getattr(page.planarconfig, 'name', str(page.planarconfig)),
                compresion=getattr(page.compression, 'name', str(page.compression)),
                bits_por_muestra=_valor(tg['BitsPerSample']) if 'BitsPerSample' in tg else None,
                software=std('Software'), fecha_tiff=std('DateTime'), marca=std('Make'),
                modelo=std('Model'), descripcion=desc[:2000] if isinstance(desc, str) else None,
                tags_privados=tags_privados(tg))


def leer_tiff(carpeta, nombres):
    """Recorre TODAS las páginas (solo los tags, no los píxeles) de los TIFF, en orden.

    Devuelve dict: archivo [n] ('nombre#página'), cuadro [n] (TSI_FRAME_NUMBER o None),
    ts_ns [n] (timestamp relativo [ns] o None), paginas [por archivo], info (de la primera
    página, ver info_pagina), ome_xml (str o None), saltadas (páginas de otra forma: miniaturas).
    """
    import tifffile
    out = dict(archivo=[], cuadro=[], ts_ns=[], paginas=[], info=None, ome_xml=None, saltadas=0)
    forma = None
    for nom in nombres:
        k = 0
        with tifffile.TiffFile(os.path.join(carpeta, nom)) as tf:
            tf.pages.useframes = False                       # TiffPage con todos sus tags
            if out['ome_xml'] is None and tf.is_ome:
                out['ome_xml'] = tf.ome_metadata
            for p, page in enumerate(tf.pages):
                if forma is None:
                    forma, out['info'] = page.shape, info_pagina(page)
                if page.shape != forma:
                    out['saltadas'] += 1
                    continue
                out['archivo'].append(f'{nom}#{p}')
                out['cuadro'].append(_entero(page.tags, TAG_CUADRO))
                out['ts_ns'].append(timestamp_ns(page.tags))
                k += 1
        out['paginas'].append(k)
    return out


def camara_tsi(tp):
    """Tags TSI (dict código → valor) → campos de `camara` (solo los que están)."""
    c = {}
    if 32768 in tp:
        c['bits'] = int(tp[32768])
    if 32770 in tp and 32771 in tp:
        c['binning'] = [int(tp[32770]), int(tp[32771])]
    if all(k in tp for k in (32772, 32773, 32774, 32775)):
        c['roi'] = [int(tp[k]) for k in (32772, 32773, 32774, 32775)]
    if 32777 in tp:
        c['pixel_clock_hz'] = float(tp[32777])
    if 32778 in tp:
        c['taps'] = int(tp[32778])
    if 32780 in tp:
        c['tipo_polarizacion'] = {0: 'Intensity', 1: 'DoLP', 2: 'Azimuth', 3: 'Unprocessed',
                                  4: 'QuadView'}.get(int(tp[32780]), int(tp[32780]))
    return c


_A_SEGUNDOS = {'s': 1.0, 'ms': 1e-3, 'µs': 1e-6, 'us': 1e-6, 'ns': 1e-9, 'min': 60.0}
_A_METROS = {'m': 1.0, 'mm': 1e-3, 'µm': 1e-6, 'um': 1e-6, 'nm': 1e-9}


def leer_ome(xml):
    """OME-XML → dict: delta_t_s [planos] (NaN si falta), exposicion_s [planos], tamano_px_m
    (PhysicalSizeX: depende de cómo esté configurado el programa, NO se usa como escala),
    modelo, serie. Lo que no está queda None/NaN."""
    raiz = ET.fromstring(xml)
    px = raiz.find('.//{*}Pixels')
    planos = [] if px is None else px.findall('{*}Plane')

    def seg(pl, clave):
        v = pl.get(clave)
        return float(v) * _A_SEGUNDOS.get(pl.get(clave + 'Unit', 's'), np.nan) if v else np.nan
    tam = None
    if px is not None and px.get('PhysicalSizeX'):
        tam = float(px.get('PhysicalSizeX')) * _A_METROS.get(px.get('PhysicalSizeXUnit', 'µm'), np.nan)
    det = raiz.find('.//{*}Detector')
    orden = sorted(range(len(planos)), key=lambda i: int(planos[i].get('TheT', i)))
    return dict(delta_t_s=np.array([seg(planos[i], 'DeltaT') for i in orden]),
                exposicion_s=np.array([seg(planos[i], 'ExposureTime') for i in orden]),
                tamano_px_m=tam, modelo=None if det is None else det.get('Model'),
                serie=None if det is None else det.get('SerialNumber'))


# ---------------------------------------------------------------- video
REPETIDO_NIVEL = 6        # un píxel "cambió" si difiere en más de 6 niveles (de 255) del cuadro anterior
REPETIDO_FRACCION = 1e-4  # un cuadro es REPETIDO si cambió menos del 0.01 % de los píxeles


def leer_video(ruta):
    """Video → dict: n (contado leyendo, exacto), fps_contenedor (nominal), forma, fourcc, mono
    (True si el primer cuadro tiene B = G = R), cv2_api (backend de cv2 que lo lee) y repetido [n]
    (bool: el cuadro es una copia del anterior).

    Cuadros repetidos: un video grabado de pantalla, o con fps variable, repite el último cuadro
    cuando la cámara no entregó uno nuevo (medido en un MP4 de los grupos anteriores: 109 de 828,
    y en ellos el desplazamiento de las beads es 0). Si se dejaran, el tiempo de cada cuadro
    estaría mal y el MSD tendría un pico falso en Δ = 0. Criterio: menos de REPETIDO_FRACCION de los
    píxeles (canal verde) cambia más de REPETIDO_NIVEL respecto del cuadro anterior; en ese video
    separa limpio (repetidos: paso de las beads 0, p95 0.23 px; nuevos: p5 0.58 px)."""
    import cv2
    cap, api = io.abrir_video(ruta)
    fps = float(cap.get(cv2.CAP_PROP_FPS) or np.nan)
    cc = int(cap.get(cv2.CAP_PROP_FOURCC))
    ok, f0 = cap.read()
    n = int(ok)
    repetido = [False] if ok else []
    prev = f0[..., 1].astype(np.int16) if ok and f0.ndim == 3 else (f0.astype(np.int16) if ok else None)
    while ok:
        ok, f = cap.read()
        if not ok:
            break
        g = (f[..., 1] if f.ndim == 3 else f).astype(np.int16)
        repetido.append(bool((np.abs(g - prev) > REPETIDO_NIVEL).mean() < REPETIDO_FRACCION))
        prev = g
        n += 1
    cap.release()
    if not n:
        raise IOError(f'{ruta}: no tiene cuadros legibles')
    mono = f0.ndim == 2 or bool(np.array_equal(f0[..., 0], f0[..., 1]) and np.array_equal(f0[..., 1], f0[..., 2]))
    return dict(n=n, fps_contenedor=fps if fps > 0 else np.nan, forma=list(f0.shape[:2]), cv2_api=int(api),
                fourcc=''.join(chr((cc >> 8 * i) & 0xFF) for i in range(4)), mono=mono,
                repetido=np.array(repetido, bool))


# ---------------------------------------------------------------- Bayer
def mover_fase(fase, dx, dy):
    """Fase Bayer del sensor → fase de una imagen cuyo píxel (0,0) es el (dy, dx) del sensor."""
    f = str(fase).upper()
    return ''.join(f[2 * ((i + dy) % 2) + (j + dx) % 2] for i in (0, 1) for j in (0, 1))


def diagonal_verde(img):
    """Medias de las 4 subredes 2×2 de un mosaico → (medias [4], diagonal del verde).

    El verde son los dos píxeles de la celda con respuesta parecida: '01-10' (RGGB, BGGR) o
    '00-11' (GRBG, GBRG). Diferencia relativa entre las 4 medias < 1 % → 'parece_mono'.
    Es un diagnóstico (avisa si la fase supuesta no cierra), no una medición de la fase.
    """
    c = np.asarray(img, float)
    c = c[: c.shape[0] // 2 * 2, : c.shape[1] // 2 * 2]
    m = np.array([c[i::2, j::2].mean() for i, j in ((0, 0), (0, 1), (1, 0), (1, 1))])
    if np.ptp(m) < 0.01 * max(m.mean(), 1e-12):
        return m, 'parece_mono'
    return m, '01-10' if abs(m[1] - m[2]) < abs(m[0] - m[3]) else '00-11'


def _es_color(modelo):
    """Modelos color de Thorlabs: terminan en C o CU (CS165CU, CS505CU, DCC1645C...)."""
    return bool(modelo) and bool(re.search(r'(CU|C)$', str(modelo).strip().upper()))


# ---------------------------------------------------------------- tiempos
def columnas_tiempos(archivo, cuadro=None, ts_ns=None, delta_t_s=None, periodo_s=None,
                     exposicion_us=None):
    """Arma las columnas de tiempos.csv y dice de dónde salen los tiempos.

    archivo [n]; cuadro [n] (contador de hardware o None por cuadro); ts_ns [n] (timestamp de
    cámara [ns] o None); delta_t_s [n] (OME DeltaT [s], NaN si falta); periodo_s (nominal, del
    usuario o del contenedor de video); exposicion_us: escalar o [n].
    Prioridad: timestamps TSI → t_camara [ns] ('camara'); si no, OME DeltaT → t_camara [ns]
    ('ome_deltat'); si no, t_nominal_s = i·periodo_s ('nominal'); si no hay nada, 'ninguno'.
    contador_camara: el número de cuadro si TODOS lo tienen, si no vacío (NaN: no hay contador).
    t_host_ns: siempre vacío (no hay reloj del host). Devuelve (columnas, tiempo_fuente).
    """
    n = len(archivo)

    def completa(v):
        return v is not None and len(v) == n and all(x is not None and np.isfinite(x) for x in v)
    nan = np.full(n, np.nan)
    cont = np.array(cuadro, dtype=np.int64) if completa(cuadro) else nan
    tc, fuente = nan, 'ninguno'
    if ts_ns is not None and any(x is not None for x in ts_ns):
        tc = np.array([np.nan if x is None else float(x) for x in ts_ns])
        fuente = 'camara'
    elif delta_t_s is not None and len(delta_t_s) == n and np.isfinite(delta_t_s).any():
        tc, fuente = np.round(np.asarray(delta_t_s, float) * 1e9), 'ome_deltat'
    tn = nan
    if periodo_s and np.isfinite(periodo_s):
        tn = np.arange(n) * float(periodo_s)
        fuente = 'nominal' if fuente == 'ninguno' else fuente
    ex = np.broadcast_to(np.asarray(np.nan if exposicion_us is None else exposicion_us, float), (n,))
    return dict(archivo=list(archivo), contador_camara=cont, t_camara=tc, t_host_ns=nan,
                exposicion_us=np.array(ex), t_nominal_s=tn), fuente


# ---------------------------------------------------------------- metadatos
def fusionar(archivo, usuario, ruta='', conflictos=None):
    """Mezcla dicts anidados: lo leído del ARCHIVO manda; el usuario completa lo que falta (None).

    Si los dos tienen un valor distinto, queda el del archivo y se anota en `conflictos`
    (lista de str). Devuelve el dict mezclado (nuevo)."""
    conflictos = [] if conflictos is None else conflictos
    out = dict(archivo)
    for k, vu in (usuario or {}).items():
        va = out.get(k)
        if isinstance(va, dict) and isinstance(vu, dict):
            out[k] = fusionar(va, vu, f'{ruta}{k}.', conflictos)
        elif va is None or va == '' or va == []:
            out[k] = vu
        elif vu is not None and vu != va:
            conflictos.append(f'{ruta}{k}: archivo={va!r}, usuario={vu!r} (queda el del archivo)')
    return out


def _base_meta():
    return dict(version=1, nombre=None, tipo=None, hora_inicio_iso=None,
                muestra=dict(material=None, diametro_m=None, glicerina_pct=None, preparada=None),
                condiciones=dict(temperatura_K=None, laser=None, profundidad_m=None,
                                 velocidad_platina_m_s=None),
                camara=dict(modelo=None, serie=None, exposicion_us=None, periodo_us=None, roi=None,
                            binning=None, ganancia=None, bits=None, formato_pixel=None,
                            bayer_fase=None, tick_hz=None),
                escala_m_por_px=None, notas='')


def _formato_pixel(info, camara, usuario_cam):
    """'rgb' si hay 3-4 muestras por píxel; si es 2D: lo del usuario, 'bayer' si el modelo es
    color, si no 'mono'. Devuelve (formato, supuesto: bool)."""
    if info['muestras'] >= 3:
        return 'rgb', False
    if usuario_cam.get('formato_pixel'):
        return str(usuario_cam['formato_pixel']).lower(), False
    if camara.get('tipo_polarizacion') is not None:
        return 'mono', True
    return ('bayer', False) if _es_color(usuario_cam.get('modelo')) else ('mono', True)


def meta_tiff(lt, usuario_cam):
    """Resultado de leer_tiff → (camara dict, avisos, encontrados) desde los tags."""
    info = lt['info']
    tp = info['tags_privados']
    cam = camara_tsi(tp)
    enc = [f'{TAGS_TSI[k]} ({k})' for k in sorted(tp) if k in TAGS_TSI]
    avisos = []
    if 'bits' not in cam and info['dtype'] in ('uint8', 'uint16'):
        cam['bits'] = int(np.atleast_1d(info['bits_por_muestra'])[0] or 16)
        avisos.append(f"bits = {cam['bits']} (BitsPerSample del contenedor: no hay TSI_SIGNIFICANT_BITS)")
    fp, supuesto = _formato_pixel(info, cam, usuario_cam)
    cam['formato_pixel'] = fp
    if supuesto:
        avisos.append('formato_pixel = mono SUPUESTO (2D, modelo no dado o no color): si la '
                      'cámara es color (CS165CU), poné formato_pixel = bayer')
    if fp == 'rgb':
        avisos.append(AVISO_SRGB)
    desconocidos = sorted(k for k in tp if k not in TAGS_TSI and k < 33000)
    if desconocidos:
        avisos.append(f'tags privados no documentados (guardados en fuente.tags_privados, sin '
                      f'interpretar): {desconocidos}')
    if info['modelo'] and not usuario_cam.get('modelo'):
        cam['modelo'] = info['modelo']
    return cam, avisos, enc


def fase_bayer(cam, usuario_cam):
    """bayer_fase del ARCHIVO: la del usuario (camara.bayer_fase) o la del sensor
    (camara.bayer_fase_sensor, si no 'RGGB' SUPUESTA) corrida por el origen del ROI.
    Devuelve (fase, supuesta: bool)."""
    if usuario_cam.get('bayer_fase'):
        return str(usuario_cam['bayer_fase']).upper(), False
    sensor = usuario_cam.get('bayer_fase_sensor')
    x0, y0 = (cam.get('roi') or [0, 0])[:2]
    return mover_fase(sensor or 'RGGB', int(x0), int(y0)), not sensor


# ---------------------------------------------------------------- importar
def importar(ruta, meta_usuario=None, raiz=None, serie=True):
    """Importa un archivo o carpeta de ThorCam / ThorImageCAM (o cualquier TIFF, secuencia de
    TIFF o video) como grabación nueva SIN copiar los cuadros. Devuelve la carpeta creada.

    ruta : archivo (.tif/.tiff/.avi/.mp4...) o carpeta. serie: ver listar_archivos.
    meta_usuario : dict con la forma de grabacion.json (SI): tipo, muestra{diametro_m,
           glicerina_pct, material}, condiciones{temperatura_K, profundidad_m}, camara{modelo,
           exposicion_us, periodo_us, ganancia, formato_pixel, bayer_fase | bayer_fase_sensor},
           escala_m_por_px, notas, y además 'etiqueta' (nombre de la carpeta). Completa lo que
           los archivos no traen; si los dos traen algo distinto manda el archivo (se anota en
           importacion.conflictos).
    raiz : carpeta de grabaciones (None = rutas.grabaciones()).
    Los originales nunca se abren para escribir.
    """
    mu = dict(meta_usuario or {})
    etiqueta = mu.pop('etiqueta', None)
    ucam = dict(mu.get('camara') or {})
    formato, carpeta_src, nombres = listar_archivos(ruta, serie)
    periodo_s = float(ucam['periodo_us']) * 1e-6 if ucam.get('periodo_us') else None
    expo = ucam.get('exposicion_us')
    avisos, enc = [], []
    fuente = dict(formato=formato, programa=None, carpeta=carpeta_src, archivos=nombres)

    if formato == 'tiff':
        lt = leer_tiff(carpeta_src, nombres)
        cam, av, enc = meta_tiff(lt, ucam)
        avisos += av
        ome = leer_ome(lt['ome_xml']) if lt['ome_xml'] else None
        dts = None
        if ome is not None:
            enc.append('OME-XML')
            if len(ome['delta_t_s']) == len(lt['archivo']):
                dts = ome['delta_t_s']
                ex = ome['exposicion_s'] * 1e6
                if np.isfinite(ex).all():
                    expo = ex if expo is None else expo
                    enc.append('OME ExposureTime')
            elif len(ome['delta_t_s']):
                avisos.append(f"OME: {len(ome['delta_t_s'])} planos y {len(lt['archivo'])} páginas: no uso DeltaT")
            if ome['tamano_px_m']:
                fuente['ome_tamano_px_m'] = ome['tamano_px_m']
                avisos.append('OME PhysicalSizeX guardado en fuente.ome_tamano_px_m: NO se usa '
                              'como escala (depende de la configuración del programa)')
            cam.setdefault('modelo', ome['modelo'])
        cols, tfuente = columnas_tiempos(lt['archivo'], lt['cuadro'], lt['ts_ns'], dts,
                                         periodo_s, expo)
        info = lt['info']
        fuente.update(paginas=lt['paginas'], dtype=info['dtype'], forma=info['forma'],
                      muestras=info['muestras'], fotometrica=info['fotometrica'],
                      compresion=info['compresion'], software=info['software'],
                      descripcion=info['descripcion'], fecha_tiff=info['fecha_tiff'],
                      tags_privados=info['tags_privados'], paginas_saltadas=lt['saltadas'],
                      programa='ThorCam (tags TSI)' if any(k in TAGS_TSI for k in info['tags_privados'])
                      else (info['software'] or None))
        if info['compresion'] not in ('NONE', '1', 'None'):
            avisos.append(f"TIFF comprimido ({info['compresion']}): verificar que sea sin pérdida")
        forma = info['forma'][:2]
        hora = hora_de_nombre(nombres[0]) or _fecha_tiff(info['fecha_tiff'])
    else:
        vi = leer_video(os.path.join(carpeta_src, nombres[0]))
        if periodo_s is None and np.isfinite(vi['fps_contenedor']):
            periodo_s = 1.0 / vi['fps_contenedor']
            avisos.append(f"fps nominal = {vi['fps_contenedor']:g} (del contenedor del video)")
        cam = dict(bits=8, formato_pixel='mono' if vi['mono'] else 'rgb')
        avisos.append(AVISO_VIDEO)
        enc.append(f"video {vi['fourcc']} {vi['n']} cuadros")
        nuevos = np.flatnonzero(~vi['repetido'])
        if len(nuevos) < vi['n']:
            avisos.append(f"{vi['n'] - len(nuevos)} de {vi['n']} cuadros son REPETIDOS (copia del anterior: "
                          f"la cámara no entregó uno nuevo) y se descartan; cada cuadro nuevo conserva el "
                          f"tiempo en que aparece en el video")
        cols, tfuente = columnas_tiempos([f'{nombres[0]}#{i}' for i in nuevos],
                                         periodo_s=periodo_s, exposicion_us=expo)
        if periodo_s and np.isfinite(periodo_s):
            cols['t_nominal_s'] = nuevos * float(periodo_s)
        fuente['cuadros_repetidos'] = int(vi['n'] - len(nuevos))
        fuente.update(fourcc=vi['fourcc'], cv2_api=vi['cv2_api'], fps_contenedor=vi['fps_contenedor'], forma=vi['forma'],
                      formato_pixel=cam['formato_pixel'], aviso=AVISO_VIDEO)
        forma = vi['forma']
        hora = hora_de_nombre(nombres[0])

    cam.setdefault('roi', [0, 0, forma[1], forma[0]])
    cam['tick_hz'] = 1e9 if tfuente in ('camara', 'ome_deltat') else None
    if np.ndim(expo) == 0 and expo is not None:
        cam['exposicion_us'] = float(expo)
    elif np.ndim(expo):
        cam['exposicion_us'] = float(np.median(expo))
    if cam['formato_pixel'] == 'bayer':
        cam['bayer_fase'], supuesta = fase_bayer(cam, ucam)
        cam['bayer_fase_supuesta'] = supuesta
        avisos += _chequeo_bayer(carpeta_src, cols['archivo'][0], fuente, cam['bayer_fase'])
    ucam.pop('bayer_fase_sensor', None)
    fuente['formato_pixel'] = cam['formato_pixel']

    conflictos = []
    meta = _base_meta()
    meta.update(camara=fusionar(dict(meta['camara'], **{k: v for k, v in cam.items() if v is not None}),
                                ucam, 'camara.', conflictos))
    meta = fusionar(meta, {k: v for k, v in mu.items() if k != 'camara'}, '', conflictos)
    meta['hora_inicio_iso'] = meta['hora_inicio_iso'] or hora
    meta['fuente'] = fuente
    meta['tiempo_fuente'] = tfuente

    carpeta = io.nueva_carpeta_grabacion(etiqueta or 'importada_' + _prefijo(nombres[0], nombres).rstrip('_-. '),
                                         raiz=raiz)
    os.rmdir(os.path.join(carpeta, 'cuadros'))       # vacía: los cuadros quedan en la fuente
    meta['nombre'] = os.path.basename(carpeta)
    io.escribir_tiempos_csv(carpeta, cols)
    rt = tiempo.reconstruir(io.leer_csv(os.path.join(carpeta, 'tiempos.csv')), tick_hz=cam['tick_hz'])
    meta['informe_tiempo'] = dict(rt['informe'], tiempo_fuente=tfuente,
                                  explicacion=EXPLICACION_TIEMPO[tfuente])
    encontrados, faltan = resumen(meta)
    meta['importacion'] = dict(fecha_iso=datetime.now().isoformat(timespec='seconds'),
                               desde=os.path.abspath(ruta), version_importador=VERSION_IMPORTADOR,
                               tags_encontrados=enc, encontrados=encontrados, faltan=faltan,
                               conflictos=conflictos, avisos=avisos)
    io.escribir_grabacion_json(carpeta, meta)
    return carpeta


EXPLICACION_TIEMPO = {
    'camara': 'timestamp de la cámara por cuadro (TSI_TIME_STAMP_RELATIVE_NS, ns) y contador '
              '(TSI_FRAME_NUMBER) si está: tiempos reales y cuadros perdidos por contador. No hay '
              'reloj del host: no se mide la deriva del reloj (ppm_reloj = NaN).',
    'ome_deltat': 'OME Plane DeltaT por cuadro (origen no confirmado: cámara o PC). Sin contador '
                  'de hardware: los cuadros perdidos NO se detectan.',
    'nominal': 'NO HAY timestamps reales: t = índice / fps nominal. No se detectan cuadros '
               'perdidos ni se mide el período, el jitter o la deriva del reloj.',
    'ninguno': 'NO HAY tiempos: ni timestamps ni fps. Cargá camara.periodo_us (--fps) y reimportá.',
}


def _fecha_tiff(s):
    """DateTime de TIFF 'AAAA:MM:DD HH:MM:SS' → ISO, o None."""
    m = re.fullmatch(r'(\d{4}):(\d{2}):(\d{2}) (\d{2}:\d{2}:\d{2})', str(s or '').strip())
    return f'{m.group(1)}-{m.group(2)}-{m.group(3)}T{m.group(4)}' if m else None


def _chequeo_bayer(carpeta, archivo0, fuente, fase):
    """Aviso si el primer cuadro no parece Bayer o si su diagonal de verde no cierra con la fase."""
    lec = io.LectorFuente(dict(fuente, carpeta=carpeta))
    try:
        img = lec.leer(*io.partir_archivo(archivo0, 0, fuente))
    finally:
        lec.cerrar()
    medias, diag = diagonal_verde(img)
    esperada = '01-10' if fase in ('RGGB', 'BGGR') else '00-11'
    if diag == 'parece_mono':
        return [f'el cuadro 0 no parece un mosaico Bayer (medias de las 4 subredes {np.round(medias, 1).tolist()})']
    if diag != esperada:
        return [f'la diagonal del verde del cuadro 0 es {diag} y la fase {fase} dice {esperada}: '
                f'revisá bayer_fase (medias {np.round(medias, 1).tolist()})']
    return []


def resumen(meta):
    """grabacion.json → (encontrados, faltan): listas de str para el informe de la importación."""
    cam, mu, co = meta.get('camara') or {}, meta.get('muestra') or {}, meta.get('condiciones') or {}
    enc, falt = [], []

    def ver(cond, si, no):
        (enc if cond else falt).append(si if cond else no)
    tf = meta.get('tiempo_fuente')
    ver(tf in ('camara', 'ome_deltat'), f'tiempos por cuadro ({tf})',
        'timestamps por cuadro: NO hay (t nominal o ninguno; sin detección de cuadros perdidos)')
    ver(meta.get('informe_tiempo', {}).get('contador') == 'hardware', 'contador de cuadros',
        'contador de cuadros: NO hay (no se detectan cuadros perdidos)')
    ver(cam.get('exposicion_us') is not None, 'exposición', 'camara.exposicion_us (--exposicion_us): LIBRE y TRAMPA la necesitan')
    ver(cam.get('periodo_us') is not None, 'período nominal', 'camara.periodo_us (--fps): solo informativo si hay timestamps')
    ver(meta.get('escala_m_por_px') is not None, 'escala', 'escala_m_por_px (--escala_nm_px) o una retícula')
    ver(bool(meta.get('tipo')), 'tipo', 'tipo (--tipo libre|atrapada|stokes|pegada|reticula)')
    ver(mu.get('diametro_m') is not None, 'diámetro', 'muestra.diametro_m (--diametro_um)')
    ver(mu.get('glicerina_pct') is not None, '% glicerina', 'muestra.glicerina_pct (--glicerina_pct) o eta a mano')
    ver(co.get('temperatura_K') is not None, 'temperatura', 'condiciones.temperatura_K (si no, ambiente por defecto)')
    ver(cam.get('ganancia') is not None, 'ganancia', 'camara.ganancia (informativo)')
    for k in ('bits', 'roi', 'binning'):
        if cam.get(k) is not None:
            enc.append(k)
    if cam.get('formato_pixel') == 'bayer':
        ver(not cam.get('bayer_fase_supuesta'), 'fase Bayer',
            f"bayer_fase = {cam.get('bayer_fase')} SUPUESTA (RGGB en el sensor corrida por el ROI): "
            f"no afecta al canal verde; para rojo/azul confirmala")
    return enc, falt
