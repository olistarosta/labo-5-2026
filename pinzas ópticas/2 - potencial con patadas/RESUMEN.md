# Objetivo 2 — Caracterizar el potencial con patadas

Con una partícula atrapada, se mueve la platina de golpe (una «patada»). El líquido arrastra a la partícula fuera del centro de la trampa y se mide con qué fuerza la trampa la trae de vuelta. De ahí sale la fuerza en función de la distancia, es decir, la forma del pozo.

## Qué hay en esta carpeta

| Archivo | Qué es |
|---|---|
| `platina y patadas.ipynb` | **primero**: cuánto se movió de verdad la platina en cada patada, contra lo anotado en el cuaderno |
| `fuerza y pozo.ipynb` | **después**: fuerza de la trampa contra la distancia, forma del pozo, cuándo se escapa, y el paso a pN |
| `figuras/` | las figuras de los dos notebooks, en SVG |
| `resultados/` | tablas (`patadas.csv`, `pozo_*.csv`, `fuerza_*.csv`, `difusion_libres_*.csv`); `LEEME.txt` explica las columnas |
| `DUDAS.md` | preguntas para responder y decisiones que se tomaron |

Las grabaciones que se usan: `BARRIDO 1`, `BARRIDO 2`, `BARRIDO STEPSIZE`, `BARRIDO STEPSIZE 2`, `BARRIDO STEPSIZE 3`, `LASER NE`, `MUEVO A MANO`. Sus tablas por cuadro están en `../datos/<grabación>/`. Las funciones están en `../codigo/pinzas_datos.py` y `../codigo/pinzas_fuerza.py`.

## Cómo se hizo

1. **La platina se mide en el video**, no con los valores del motor: se mira cuánto se corre toda la imagen de un cuadro al siguiente.
2. La partícula no tiene inercia. En cada instante la fuerza de la trampa iguala al roce con el líquido: $F = \gamma\,(v_{partícula} - v_{líquido})$. Midiendo la velocidad de la partícula respecto del líquido se tiene $F/\gamma$ en cada posición.
3. Se ajusta un pozo gaussiano a todas las patadas juntas.
4. Para pasar a unidades de fuerza hace falta $\gamma$. Sale de la difusión de las partículas libres de la misma muestra: $\gamma = k_BT/D$.

## Resultados (a 30 mW)

**La platina**

- Las horas del cuaderno coinciden con el video (0.15 s de diferencia mediana en 31 patadas).
- Los valores del motor no sirven: el desplazamiento real es en mediana el 5 % del paso pedido y no crece con él.
- Una patada típica corre la imagen 5.7 µm en 0.32 s: son 2 a 6 cuadros.

**Fuerza y forma del pozo**

| Magnitud | Valor |
|---|---|
| $\kappa/\gamma$ | 38 ± 4 s⁻¹ |
| Tiempo de relajación | 26 ms (frecuencia de esquina 6.1 Hz) |
| Ancho del pozo (máximo de la fuerza) | 1.14 ± 0.09 µm |
| Alcance | la fuerza se anula pasados 3.5 µm; la partícula volvió sola desde 5.6 µm |
| Arrastre $\gamma$ | (5.0 ± 0.2) × 10⁻⁸ N·s/m, de la difusión de las partículas libres de la misma muestra |
| Rigidez | 1.9 ± 0.2 pN/µm |
| Fuerza máxima | 1.3 ± 0.1 pN |
| Profundidad del pozo | 600 ± 80 $k_BT$ |

- El pozo es gaussiano dentro de lo que se puede ver, e igual en vertical y en horizontal dentro de un error de 25 a 30 %.
- La partícula se escapa con patadas que corren la imagen más de unos 8 µm (pico de unos 35 µm/s por cuadro).
- Errores estadísticos. Sistemáticos conocidos: hasta 8 % por el método y hasta 9 % por el punto oscuro del láser.
- La rigidez coincide con la del barrido de potencias (2.0 pN/µm a 30 mW), que sale por otro camino.

## Lo que no se pudo hacer

- **Ajustar el pozo** en `BARRIDO 2` (partículas chicas), `BARRIDO STEPSIZE` (todas las patadas con escape) y `LASER NE` (solo una deriva lenta).
- **Explicar tres escapes** con patadas chicas (`BARRIDO 1` n = 7, `BARRIDO STEPSIZE` n = 8 y n = 15).
- **Medir la fuerza a lo largo del eje del láser** y la altura de la partícula sobre el vidrio.

## Pendiente

- **¿Qué mezcla tenía la muestra?** La receta del cuaderno (60 ml de agua + 65 ml de glicerina) daría una difusión 4.7 veces más lenta que la medida; los datos corresponden a alrededor de 24 % de glicerina. No cambia los resultados, sí lo que se informe como muestra.
- **¿Qué potencia había en las grabaciones con patadas?** Se supuso 30 mW.
- El detalle está en `DUDAS.md`.
