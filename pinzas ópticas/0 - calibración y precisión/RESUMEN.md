# Calibración y precisión — lo que usan los tres objetivos

No es un objetivo en sí: es lo que hace falta para confiar en los otros tres. Cuánto mide un píxel, qué es cada grabación, y con qué precisión se ubica el centro de una partícula.

## Qué hay en esta carpeta

| Archivo | Qué es |
|---|---|
| `escala y catálogo.ipynb` | la escala con la retícula de 0.01 mm y la lista de grabaciones (qué es cada una, cuadros por segundo) |
| `resolución y precisión.ipynb` | resolución del video, error del centro de la partícula y forma de las partículas |
| `REPORTE resolución y precisión.md` | lo mismo que el segundo notebook, para leer sin correr nada |
| `figuras/` | las figuras de los dos notebooks, en SVG |
| `resultados/` | tablas (`resolucion_*.csv`); `LEEME.txt` explica las columnas |

La escala y el catálogo que salen del primer notebook quedan en `../datos/` (`calibracion.csv`, `catalogo.csv`), porque los usan todos. Las funciones están en `../codigo/pinzas_datos.py` y `../codigo/resolucion.py`.

## Resultados

**Escala y catálogo**

- **Escala:** 94.61 ± 0.04 nm/px con la retícula de 0.01 mm (error solo estadístico). Campo de 121 × 97 µm.
- **Cuadros por segundo:** 10.02, 13.13 y 12.50 según la grabación, estimados con las fechas de los archivos. Los TIFF no traen el tiempo de cada cuadro ni la exposición.

**Resolución y precisión**

| | Valor |
|---|---|
| Píxel | 94.6 nm |
| Nitidez de la óptica | 0.5 µm (cota) |
| Error del centro por cuadro, partícula pegada al vidrio | 16 a 47 nm, mediana 33 |
| Error del centro por cuadro, partícula atrapada | 19 a 39 nm |
| Error del centro por cuadro, partícula libre | menos de 20 a 95 nm, mediana 55 |
| Error sistemático (linealidad por debajo del píxel) | menos de 1 nm |

- El centro se mide por debajo del píxel porque la partícula ocupa más de mil píxeles. No hay salto de píxel en píxel.
- Los resultados de fuerza, rigidez y difusión no cambian con esto.
- Las partículas pegadas se ven redondas dentro de un 5 %; no se puede decir si ese 5 % es de la partícula.
- Comparando con la imagen de la propia partícula el error bajaría a 5–7 nm. No se hizo.

## Pendiente

- **¿Qué tamaños nominales de partículas se mezclaron?** Se ven dos tamaños.
- **Temperatura** del laboratorio (se usó 23 °C).
