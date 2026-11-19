"""Dónde va cada cosa. El análisis está ordenado por objetivo: una carpeta por objetivo, y adentro
los notebooks, `figuras/` (SVG) y `resultados/` (CSV).

    from carpetas import figuras, resultados
    FIGURAS = figuras("patadas")                       # .../2 - potencial con patadas/figuras
    tabla.to_csv(os.path.join(resultados("patadas"), "patadas.csv"))

Lo que sale de cada grabación (trayectorias, platina, láser, imágenes de control) no es de un
objetivo: va en datos/<grabación>/.
"""
import os

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))          # pinzas ópticas/
DATOS = os.path.join(RAIZ, 'datos')

OBJETIVOS = {
    'comun': '0 - calibración y precisión',                 # lo que usan los tres: escala, catálogo, resolución
    'browniano': '1 - movimiento browniano',
    'patadas': '2 - potencial con patadas',
    'potencia': '3 - potencial con barrido de potencias',
}


def carpeta(objetivo, *partes):
    """Carpeta del objetivo ('comun', 'browniano', 'patadas' o 'potencia'), o una subcarpeta. La crea si falta."""
    ruta = os.path.join(RAIZ, OBJETIVOS[objetivo], *partes)
    os.makedirs(ruta, exist_ok=True)
    return ruta


def figuras(objetivo, *sub):
    """Carpeta de figuras del objetivo (solo SVG)."""
    return carpeta(objetivo, 'figuras', *sub)


def resultados(objetivo):
    """Carpeta de tablas (CSV) con los resultados del objetivo."""
    return carpeta(objetivo, 'resultados')
