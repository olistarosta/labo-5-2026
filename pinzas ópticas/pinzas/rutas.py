"""Dónde vive cada cosa en disco."""
import os
import sys

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # pinzas ópticas/
TD = os.path.join(RAIZ, 'TouchDesigner')                               # el .toe y su código (no va a git)
LIB = os.path.join(TD, 'lib')                                          # scipy, tifffile, ... para el Python de TD


def preparar_sys_path():
    """Agrega TouchDesigner/lib/ (scipy, tifffile, ...) al FINAL de sys.path: así nunca tapa el numpy de TD."""
    # solo en el Python de TD (3.11): lib/ trae binarios compilados para esa versión y en otro
    # Python (VS Code) romperían los imports
    if sys.version_info[:2] == (3, 11) and os.path.isdir(LIB) and LIB not in sys.path:
        sys.path.append(LIB)


def datos():
    """Carpeta raíz de los datos (fuera de OneDrive). Variable de entorno PINZAS_DATOS o C:/datos_pinzas."""
    return os.environ.get('PINZAS_DATOS', 'C:/datos_pinzas')


def grabaciones():
    return os.path.join(datos(), 'GRABACIONES')


def cache(nombre_grabacion):
    return os.path.join(datos(), '_cache', nombre_grabacion)
