# Estilo de figuras (gráficos de datos)

Guía para que los gráficos de los informes se vean iguales. Es el estilo de las
figuras del informe de láser (`laser/figuras/`), con los ajustes que usamos en el
barrido angular de muones (`muones/figuras/fig_barrido_angular.py`, el ejemplo de
referencia).

Los esquemas del montaje (`muones/figuras/fig_modulo.py`, `fig_disposiciones.py`)
tienen su propio estilo (`iso.py`); esta guía es para gráficos de datos.

## Configuración

```python
import numpy as np
import matplotlib.pyplot as plt

plt.style.use("seaborn-v0_8")          # fondo gris claro, grilla blanca
fs = 14
plt.rcParams.update({
    "font.family": "serif",
    "font.size": fs, "axes.labelsize": fs,
    "xtick.labelsize": fs, "ytick.labelsize": fs, "legend.fontsize": fs,
    "mathtext.default": "regular",     # símbolos en el mismo tipo que el texto
    "mathtext.fontset": "cm",
    "legend.frameon": True,
})

fig, ax = plt.subplots(figsize=(8, 4.5))   # apaisado; (8, 4) también vale
...
fig.tight_layout()
fig.savefig("nombre.pdf", bbox_inches="tight")              # la que va al informe
fig.savefig("nombre.png", bbox_inches="tight", dpi=300)     # vista rápida (*.png está en .gitignore)
```

## Colores

| Uso | Color | Hex |
|---|---|---|
| Serie 1 | verde azulado | `#2A9D8F` |
| Serie 2 | coral | `#E76F51` |
| Ajuste / curva teórica | gris pizarra | `#264653` |

- Si hay una sola serie, va en `#2A9D8F` y el ajuste en `#264653`.
- Cada serie tiene su color y su marcador (por ejemplo, rombos `D` y círculos `o`),
  así se distinguen en blanco y negro y con daltonismo.
- Los datos excluidos del ajuste van con **el color de su serie y el marcador
  hueco** (`mfc="white"`), no en otro color.
- Si hacen falta más colores, la paleta de seaborn de las figuras de láser:
  azul `#4C72B0`, verde `#55A868`, rojo `#C44E52`, violeta `#8172B2`,
  amarillo `#CCB974`, celeste `#64B5CD`.

## Marcas

```python
ax.errorbar(x, y, yerr=dy, fmt="o", ms=7, capsize=3,
            elinewidth=1.6, capthick=1.6, mew=1.6,
            color=c, mfc=c, mec=c, zorder=3, label="...")     # datos arriba
ax.plot(x_fit, y_fit, color="#264653", lw=2, zorder=2, label="Ajuste ...")
```

- Datos con barras de error, nunca unidos con líneas.
- Las curvas de ajuste son finas y van detrás de los puntos.
- Si dos mediciones caen en el mismo $x$, se corren un poco (~1 % del rango)
  para que no se tapen.

## Ejes

- Etiqueta con magnitud y unidad entre corchetes: `Tasa [eventos/min]`,
  `Corriente de bombeo [A]`.
- Las variables con su símbolo en math: `r"Ángulo polar $\theta$"`.
- Ticks en valores redondos. Si el eje tiene un sentido físico, se pone en los
  ticks (`45° S … 0° … 45° N`) en vez de explicarlo en otro lado.
- El eje $y$ empieza en 0 cuando la magnitud es una tasa, un conteo o una potencia.

## Leyenda

- Pocas entradas y descriptivas: qué es cada serie ("Hacia el sur"), no detalles
  técnicos.
- El ajuste se nombra por su modelo: `r"Ajuste $R_0\,\cos^n\theta$"`, `"Ajuste lineal"`.
- **No van en la leyenda ni en la figura**: valores de parámetros, $\chi^2$,
  p-valores ni el método de ajuste. Todo eso va en el pie de figura o en el texto.
- Ubicarla en la zona vacía del gráfico (`loc=` explícito) sin tapar datos; si
  hace falta, `ncol=2` y `fontsize=12`. Fondo blanco: `facecolor="white", framealpha=0.9`.

## Qué no poner

- Títulos dentro de la figura (el título es el pie de figura).
- Paneles de residuos, bandas de confianza o cuadros de texto con resultados,
  salvo que el análisis los necesite para el argumento.
- Más de un eje $y$.

## Pie de figura

Una oración que dice qué se grafica, después cómo se obtuvo cada punto y sus barras
de error, qué es la curva y los resultados del ajuste con incerteza. El análisis va
en el texto, que presenta la figura en una oración y pasa directo a los resultados.
