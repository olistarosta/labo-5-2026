"""Leer y escribir grabaciones y análisis (formato de CONTRATO.md).

Una grabación es una carpeta con `grabacion.json`, `tiempos.csv` y `cuadros/`. Nada de este
módulo modifica una grabación existente: los análisis van a `analisis/<AAAA-MM-DD_HHMMSS>/`.

Convenciones de lectura de `tiempos.csv`:
  - `archivo` → array de str.
  - `t_camara` → SIEMPRE float64; celda vacía → NaN (sin timestamp de hardware). Los ticks son
    enteros: float64 los representa exactos mientras |t| < 2^53 ≈ 9.0e15 (104 días a 1 GHz).
  - otras columnas → int64 si todas las celdas son enteras y no vacías; si no, float64 con NaN
    en las vacías.
En JSON, numpy pasa a tipos nativos y los no finitos (NaN, ±inf) se escriben como null.
"""
import csv
import glob
import json
import math
import os
import re
from collections import OrderedDict
from datetime import datetime

import numpy as np

from . import rutas

COLUMNAS_TIEMPOS = ['archivo', 'contador_camara', 't_camara', 't_host_ns', 'exposicion_us']
_SIEMPRE_FLOAT = {'t_camara'}


# ---------------------------------------------------------------- JSON
def a_nativo(obj):
    """Convierte (recursivamente) numpy y no finitos a tipos que JSON acepta.

    ndarray → lista; np.integer/np.floating/np.bool_ → int/float/bool; NaN/inf → None;
    tuplas → listas. Lo demás queda igual.
    """
    if isinstance(obj, dict):
        return {str(k): a_nativo(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [a_nativo(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return a_nativo(obj.tolist())
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        f = float(obj)
        return f if math.isfinite(f) else None
    return obj


def escribir_json(ruta, datos):
    """Escribe `datos` en `ruta` como JSON UTF-8 legible (ensure_ascii=False, indent=2)."""
    with open(ruta, 'w', encoding='utf-8') as f:
        json.dump(a_nativo(datos), f, ensure_ascii=False, indent=2, allow_nan=False)


def leer_json(ruta):
    with open(ruta, 'r', encoding='utf-8') as f:
        return json.load(f)


def escribir_grabacion_json(carpeta, meta):
    """Escribe `<carpeta>/grabacion.json` con el dict `meta` (ver CONTRATO.md). Devuelve la ruta."""
    ruta = os.path.join(carpeta, 'grabacion.json')
    escribir_json(ruta, meta)
    return ruta


escribir_grabacion = escribir_grabacion_json      # nombre que usa PLAN.md


# ---------------------------------------------------------------- CSV
def _celda(v):
    """Valor → texto de CSV. NaN/None → vacío; enteros sin '.0'; floats con repr (ida y vuelta exacta)."""
    if v is None:
        return ''
    if isinstance(v, (np.bool_, bool)):
        return str(int(v))
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        f = float(v)
        if not math.isfinite(f):
            return ''
        if f.is_integer() and abs(f) < 2 ** 53:
            return str(int(f))
        return repr(f)
    return str(v)


def escribir_csv(ruta, columnas, orden=None):
    """Escribe un CSV por columnas. columnas: dict nombre → array/lista (mismo largo).

    orden: lista de nombres que van primero (el resto sigue en el orden del dict).
    NaN y None se escriben como celda vacía.
    """
    nombres = [c for c in (orden or []) if c in columnas]
    nombres += [c for c in columnas if c not in nombres]
    datos = [np.asarray(columnas[c]) if not isinstance(columnas[c], list) else columnas[c]
             for c in nombres]
    largos = {len(d) for d in datos}
    if len(largos) > 1:
        raise ValueError(f'columnas de distinto largo: {dict(zip(nombres, map(len, datos)))}')
    n = largos.pop() if largos else 0
    with open(ruta, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(nombres)
        for i in range(n):
            w.writerow([_celda(d[i]) for d in datos])
    return ruta


def _columna(nombre, textos):
    """Lista de str de una columna de CSV → array numpy con las convenciones del módulo."""
    if nombre == 'archivo':
        return np.array(textos, dtype=str)
    vacias = [t.strip() == '' for t in textos]
    if nombre not in _SIEMPRE_FLOAT and not any(vacias):
        try:
            return np.array([int(t) for t in textos], dtype=np.int64)
        except ValueError:
            pass
    try:
        return np.array([float(t) if not v else np.nan for t, v in zip(textos, vacias)])
    except ValueError:
        return np.array(textos, dtype=str)


def leer_csv(ruta):
    """Lee un CSV con encabezado → dict nombre → array (ver convenciones arriba)."""
    with open(ruta, 'r', newline='', encoding='utf-8') as f:
        filas = list(csv.reader(f))
    if not filas:
        return {}
    nombres, cuerpo = filas[0], [r for r in filas[1:] if r]
    return {c: _columna(c, [r[k] if k < len(r) else '' for r in cuerpo])
            for k, c in enumerate(nombres)}


def escribir_tiempos_csv(carpeta, columnas):
    """Escribe `<carpeta>/tiempos.csv` (una fila por cuadro GUARDADO; columnas de CONTRATO.md primero).

    columnas: dict con al menos archivo, contador_camara, t_camara, t_host_ns, exposicion_us.
    t_camara NaN → celda vacía (no hay timestamp de hardware). Devuelve la ruta.
    """
    return escribir_csv(os.path.join(carpeta, 'tiempos.csv'), columnas, orden=COLUMNAS_TIEMPOS)


def escribir_cuadro(carpeta, indice, img):
    """Guarda un cuadro crudo uint16 como `<carpeta>/cuadros/<indice:06d>.tif` (TIFF sin compresión).

    Devuelve el nombre de archivo (lo que va en la columna `archivo`).
    """
    import tifffile
    d = os.path.join(carpeta, 'cuadros')
    os.makedirs(d, exist_ok=True)
    nombre = f'{int(indice):06d}.tif'
    tifffile.imwrite(os.path.join(d, nombre), np.asarray(img, dtype=np.uint16))
    return nombre


# ---------------------------------------------------------------- carpetas
def _sello(ahora=None):
    return (ahora or datetime.now()).strftime('%Y-%m-%d_%H%M%S')


def _crear_sin_pisar(base):
    """Crea la carpeta `base`; si existe, prueba base_2, base_3, ... Devuelve la ruta creada."""
    k = 1
    ruta = base
    while True:
        try:
            os.makedirs(ruta, exist_ok=False)
            return ruta
        except FileExistsError:
            k += 1
            ruta = f'{base}_{k}'


def nueva_carpeta_grabacion(etiqueta, raiz=None):
    """Crea `<raiz>/<AAAA-MM-DD_HHMMSS>_<etiqueta>/cuadros/` y devuelve la carpeta de la grabación.

    raiz: None = rutas.grabaciones() (PINZAS_DATOS o C:/datos_pinzas, fuera de OneDrive).
    Nunca pisa una carpeta existente: si el nombre ya existe agrega _2, _3, ...
    Los caracteres que Windows no admite en nombres (<>:"/\\|?*) y los espacios pasan a '_'.
    """
    raiz = raiz or rutas.grabaciones()
    os.makedirs(raiz, exist_ok=True)
    etiqueta = re.sub(r'[<>:"/\\|?*\s]+', '_', str(etiqueta)).strip('_')
    nombre = _sello() + (f'_{etiqueta}' if etiqueta else '')
    carpeta = _crear_sin_pisar(os.path.join(raiz, nombre))
    os.makedirs(os.path.join(carpeta, 'cuadros'))
    return carpeta


def nueva_carpeta_analisis(grabacion):
    """Crea `<grabacion>/analisis/<AAAA-MM-DD_HHMMSS>/figuras/` y devuelve la carpeta del análisis.

    grabacion: objeto Grabacion o ruta de la carpeta. Nunca pisa un análisis anterior.
    """
    carpeta = getattr(grabacion, 'carpeta', grabacion)
    base = os.path.join(carpeta, 'analisis', _sello())
    os.makedirs(os.path.dirname(base), exist_ok=True)
    ruta = _crear_sin_pisar(base)
    os.makedirs(os.path.join(ruta, 'figuras'))
    return ruta


def _a_columnas(resultados):
    """dict de escalares → 1 fila; dict de arrays → columnas; lista de dicts → filas."""
    if isinstance(resultados, dict):
        vals = list(resultados.values())
        if all(np.ndim(v) == 0 for v in vals):
            return {k: [v] for k, v in resultados.items()}
        return {k: np.atleast_1d(v) for k, v in resultados.items()}
    nombres = []
    for fila in resultados:
        nombres += [k for k in fila if k not in nombres]
    return {k: [fila.get(k) for fila in resultados] for k in nombres}


def guardar_analisis(carpeta, parametros, resultados):
    """Escribe `parametros.json` y `resultados.csv` en una carpeta de análisis (nueva_carpeta_analisis).

    parametros: dict (todos los parámetros del análisis, en un solo lugar).
    resultados: dict de escalares (una fila), dict de arrays del mismo largo (columnas) o
                lista de dicts (una fila cada uno). NaN → celda vacía.
    Si la carpeta ya tiene parametros.json o resultados.csv, no pisa: FileExistsError.
    Devuelve (ruta_parametros, ruta_resultados).
    """
    rp = os.path.join(carpeta, 'parametros.json')
    rr = os.path.join(carpeta, 'resultados.csv')
    for r in (rp, rr):
        if os.path.exists(r):
            raise FileExistsError(f'{r} ya existe: los análisis no se pisan')
    escribir_json(rp, parametros)
    escribir_csv(rr, _a_columnas(resultados))
    return rp, rr


# ---------------------------------------------------------------- lectura
def _indice_bloques(dir_cuadros):
    """Lista [(ruta_bloque, página)] en orden global, para cuadros/bloque_*.tif multipágina."""
    import tifffile
    indice = []
    for ruta in sorted(glob.glob(os.path.join(dir_cuadros, 'bloque_*.tif'))):
        with tifffile.TiffFile(ruta) as tf:
            indice += [(ruta, p) for p in range(len(tf.pages))]
    return indice


class Grabacion:
    """Una grabación en disco, con lectura perezosa de cuadros.

    .carpeta  ruta de la carpeta          .nombre  nombre de la carpeta
    .meta     dict de grabacion.json      .tiempos dict de arrays de tiempos.csv (ver módulo)
    len(g)    cuadros guardados (filas de tiempos.csv; si no hay csv, TIFFs encontrados)
    g.cuadro(i) → np.ndarray uint16 [alto, ancho] ([alto, ancho, 3] si es RGB), de solo lectura
                (vive en un caché LRU).

    Dónde está el cuadro i (en este orden):
      1. `archivo` = 'bloque_000000.tif#17' → página 17 de ese bloque multipágina.
      2. `archivo` = '000123.tif' y existe cuadros/000123.tif → ese TIFF (página 0).
      3. si hay cuadros/bloque_*.tif → la página global i (bloques ordenados por nombre).
    Un bloque se escribe con tifffile.imwrite(ruta, pila [k, alto, ancho] uint16,
    photometric='minisblack'): sin eso, una pila de 3 o 4 cuadros se toma como RGB(A).

    Grabación IMPORTADA (pinzas.importar; grabacion.json tiene `fuente`): los cuadros NO están
    en cuadros/, se leen de los archivos originales (que nunca se modifican), ver LectorFuente.
    `archivo` = 'nombre.tif#17' (página 17 de un TIFF de fuente.archivos), 'nombre.tif' (página
    0) o 'video.avi#17' (cuadro 17). Si los cuadros son RGB (camara.formato_pixel = 'rgb'),
    g.cuadro(i) es [alto, ancho, 3].
    """

    def __init__(self, carpeta, cache=16):
        self.carpeta = os.path.abspath(carpeta)
        self.nombre = os.path.basename(self.carpeta.rstrip('/\\'))
        rj = os.path.join(self.carpeta, 'grabacion.json')
        rt = os.path.join(self.carpeta, 'tiempos.csv')
        self.meta = leer_json(rj) if os.path.exists(rj) else {}
        self.tiempos = leer_csv(rt) if os.path.exists(rt) else {}
        self._dir = os.path.join(self.carpeta, 'cuadros')
        self._bloques = None
        self._cache = OrderedDict()
        self._max_cache = int(cache)
        self.fuente = self.meta.get('fuente') or None
        self._lector = None

    def __len__(self):
        if 'archivo' in self.tiempos:
            return len(self.tiempos['archivo'])
        if self.tiempos:
            return len(next(iter(self.tiempos.values())))
        b = self._indice()
        return len(b) if b else len(glob.glob(os.path.join(self._dir, '[0-9]*.tif')))

    def __repr__(self):
        return f'Grabacion({self.nombre!r}, {len(self)} cuadros)'

    def _indice(self):
        if self._bloques is None:
            self._bloques = _indice_bloques(self._dir)
        return self._bloques

    def _ubicar(self, i):
        """Índice de cuadro → (ruta, página)."""
        arch = str(self.tiempos['archivo'][i]) if 'archivo' in self.tiempos else ''
        if '#' in arch:
            nombre, pag = arch.rsplit('#', 1)
            return os.path.join(self._dir, nombre), int(pag)
        if arch and os.path.exists(os.path.join(self._dir, arch)):
            return os.path.join(self._dir, arch), 0
        if self._indice():
            return self._indice()[i]
        if not arch:
            return os.path.join(self._dir, f'{i:06d}.tif'), 0
        raise FileNotFoundError(f'no encuentro el cuadro {i} ({arch}) en {self._dir}')

    def cuadro(self, i):
        """Cuadro crudo i (0 = primero guardado) como uint16 [alto, ancho], de solo lectura."""
        import tifffile
        n = len(self)
        i = int(i)
        if i < 0:
            i += n
        if not 0 <= i < n:
            raise IndexError(f'cuadro {i} fuera de rango (hay {n})')
        if i in self._cache:
            self._cache.move_to_end(i)
            return self._cache[i]
        if self.fuente:
            img = self._leer_fuente(i)
        else:
            ruta, pag = self._ubicar(i)
            img = np.asarray(tifffile.imread(ruta, key=pag)).astype(np.uint16, copy=False)
        img.flags.writeable = False
        self._cache[i] = img
        if len(self._cache) > self._max_cache:
            self._cache.popitem(last=False)
        return img

    def _leer_fuente(self, i):
        """Cuadro i de la fuente externa (archivos originales), como uint16."""
        if self._lector is None:
            self._lector = LectorFuente(self.fuente)
        arch = str(self.tiempos['archivo'][i]) if 'archivo' in self.tiempos else ''
        img = self._lector.leer(*partir_archivo(arch, i, self.fuente))
        return np.asarray(img).astype(np.uint16, copy=False)

    def cerrar(self):
        """Cierra los archivos de la fuente externa que hayan quedado abiertos (si hay)."""
        if self._lector is not None:
            self._lector.cerrar()
            self._lector = None


# ---------------------------------------------------------------- fuente externa (importadas)
def partir_archivo(arch, i, fuente):
    """Columna `archivo` de una grabación importada → (nombre, página).

    'x.tif#17' → ('x.tif', 17); 'x.tif' → ('x.tif', 0); vacío → (único archivo de la fuente, i).
    """
    m = re.fullmatch(r'(.*)#(\d+)', arch or '')
    if m:
        return m.group(1), int(m.group(2))
    if arch:
        return arch, 0
    archivos = fuente.get('archivos') or []
    if len(archivos) != 1:
        raise ValueError(f'cuadro {i}: sin columna archivo y la fuente tiene {len(archivos)} archivos')
    return archivos[0], int(i)


def ruta_fuente(fuente, nombre):
    """Ruta absoluta de un archivo de la fuente: fuente['carpeta'] + nombre (o nombre si es absoluto).

    Si moviste los originales, cambiá `fuente.carpeta` en grabacion.json (es lo único que hay que tocar).
    """
    if os.path.isabs(nombre):
        return nombre
    return os.path.join(fuente.get('carpeta') or '', nombre)


def ruta_ascii(ruta):
    """cv2.VideoCapture no abre rutas con acentos en Windows ('pinzas ópticas'): devuelve una ruta
    ASCII al mismo archivo. Primero el nombre corto 8.3 de Windows; si el disco no los tiene, una
    copia en <datos>/_cache/videos (se hace una vez; el original no se toca)."""
    ruta = os.path.abspath(ruta)
    if ruta.isascii() or os.name != 'nt':
        return ruta
    import ctypes
    buf = ctypes.create_unicode_buffer(32768)
    if ctypes.windll.kernel32.GetShortPathNameW(ruta, buf, len(buf)) and buf.value.isascii():
        return buf.value
    import hashlib
    import shutil
    from . import rutas
    st = os.stat(ruta)
    clave = hashlib.sha1(('%s|%d|%d' % (ruta, st.st_size, st.st_mtime_ns)).encode()).hexdigest()[:16]
    copia = os.path.join(rutas.datos(), '_cache', 'videos', clave + os.path.splitext(ruta)[1].lower())
    if not os.path.exists(copia):
        os.makedirs(os.path.dirname(copia), exist_ok=True)
        shutil.copyfile(ruta, copia + '.tmp')
        os.replace(copia + '.tmp', copia)
    return copia


def abrir_video(ruta, api=None):
    """cv2.VideoCapture que de verdad lee: prueba `api` y si no CAP_ANY, CAP_OPENCV_MJPEG
    (el Python de TD no trae FFMPEG y MSMF no lee los AVI Motion-JPEG), CAP_MSMF, CAP_FFMPEG.
    Devuelve (captura posicionada en el cuadro 0, api usada). IOError si ninguna lee."""
    import cv2
    ruta = ruta_ascii(ruta)
    apis = ([int(api)] if api is not None else []) + [cv2.CAP_ANY, cv2.CAP_OPENCV_MJPEG,
                                                       cv2.CAP_MSMF, cv2.CAP_FFMPEG]
    for a in apis:
        cap = cv2.VideoCapture(ruta, a)
        if cap.isOpened() and cap.read()[0]:
            cap.set(cv2.CAP_PROP_POS_FRAMES, 0)
            return cap, a
        cap.release()
    raise IOError(f'cv2 no puede leer {ruta} con ningún backend')


class LectorFuente:
    """Lee cuadros de los archivos ORIGINALES de una grabación importada, sin modificarlos.

    fuente (de grabacion.json): formato 'tiff' (un TIFF por cuadro, multipágina o varios
    multipágina) o 'video' (AVI/MP4 vía cv2); carpeta; archivos; formato_pixel.
    TIFF: deja abiertos hasta `max_abiertos` archivos (tifffile.TiffFile); la página k se lee
    sola, sin cargar el resto: tifffile recorre la cadena de IFDs una vez y guarda solo los
    offsets. RGB con planos separados ([3, alto, ancho]) se pasa a [alto, ancho, 3].
    Video: cv2.VideoCapture, BGR → RGB; si se pide un cuadro que no es el siguiente, busca con
    CAP_PROP_POS_FRAMES (exacto en MJPEG; en MP4 depende del códec).
    """

    def __init__(self, fuente, max_abiertos=8):
        self.fuente = fuente
        self.formato = str(fuente.get('formato', 'tiff')).lower()
        self._abiertos = OrderedDict()
        self._max = int(max_abiertos)
        self._siguiente = {}

    def _abrir(self, nombre):
        if nombre in self._abiertos:
            self._abiertos.move_to_end(nombre)
            return self._abiertos[nombre]
        ruta = ruta_fuente(self.fuente, nombre)
        if not os.path.exists(ruta):
            raise FileNotFoundError(f'no encuentro {ruta}: si moviste los originales, cambiá '
                                    f'fuente.carpeta en grabacion.json')
        if self.formato == 'video':
            h, _ = abrir_video(ruta, self.fuente.get('cv2_api'))
            self._siguiente[nombre] = 0
        else:
            import tifffile
            h = tifffile.TiffFile(ruta)
        self._abiertos[nombre] = h
        if len(self._abiertos) > self._max:
            viejo, hv = self._abiertos.popitem(last=False)
            self._cerrar_uno(viejo, hv)
        return h

    def leer(self, nombre, pagina):
        """Cuadro (nombre, página) → ndarray [alto, ancho] o [alto, ancho, 3] con el dtype del archivo."""
        h = self._abrir(nombre)
        if self.formato == 'video':
            return self._leer_video(nombre, h, int(pagina))
        pag = h.pages[int(pagina)]
        img = pag.asarray()
        if img.ndim == 3 and img.shape[0] in (3, 4) and img.shape[-1] not in (3, 4):
            img = np.moveaxis(img, 0, -1)                 # planarconfig separado
        return img[..., :3] if img.ndim == 3 else img

    def _leer_video(self, nombre, cap, k):
        """Cuadro k EXACTO. En MJPEG (cada cuadro es independiente) se busca directo. En MP4/H.264
        la búsqueda de cv2 puede caer en el cuadro de al lado (medido: pedir 600 da el 599), así
        que se avanza leyendo en orden, y para ir hacia atrás se reabre el archivo. Leer todo el
        video en orden (PROCESAR) es igual de rápido; mover el cuadro de TD hacia atrás tarda
        ~8 ms por cuadro desde el principio."""
        import cv2
        actual = self._siguiente.get(nombre)
        if actual != k:
            if str(self.fuente.get('fourcc', '')).upper() in ('MJPG', 'MJPEG'):
                cap.set(cv2.CAP_PROP_POS_FRAMES, k)
            else:
                if actual is None or k < actual:
                    cap.release()
                    cap, _ = abrir_video(ruta_fuente(self.fuente, nombre), self.fuente.get('cv2_api'))
                    self._abiertos[nombre] = cap
                    actual = 0
                for _ in range(k - actual):
                    if not cap.grab():
                        raise IndexError(f'{nombre}: no pude llegar al cuadro {k}')
        ok, bgr = cap.read()
        if not ok:
            raise IndexError(f'{nombre}: no pude leer el cuadro {k}')
        self._siguiente[nombre] = k + 1
        if bgr.ndim == 2 or str(self.fuente.get('formato_pixel')) == 'mono':
            return bgr if bgr.ndim == 2 else bgr[..., 1]
        return bgr[..., ::-1]

    def _cerrar_uno(self, nombre, h):
        try:
            h.release() if self.formato == 'video' else h.close()
        except Exception:
            pass
        self._siguiente.pop(nombre, None)

    def cerrar(self):
        while self._abiertos:
            self._cerrar_uno(*self._abiertos.popitem())


def leer_grabacion(carpeta, cache=16):
    """Abre una grabación (no lee los cuadros hasta que se piden). Devuelve un objeto Grabacion."""
    return Grabacion(carpeta, cache=cache)
