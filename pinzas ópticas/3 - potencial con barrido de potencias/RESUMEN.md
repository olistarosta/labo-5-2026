# Objetivo 3 — Caracterizar el potencial con el barrido de potencias

Con la partícula atrapada y la platina quieta, se baja la potencia del láser por escalones. La partícula tiembla alrededor del centro de la trampa por agitación térmica: cuanto más floja la trampa, más se aleja. De cuánto se aleja sale la rigidez del pozo a cada potencia.

## Qué hay en esta carpeta

| Archivo | Qué es |
|---|---|
| `los ocho barridos con la potencia calibrada.ipynb` | **el resultado**: la rigidez contra la potencia real, con los ocho barridos (12/11 y 19/11) |
| `barrido de potencias.ipynb` | **el método**, con los dos primeros barridos: cómo se mide la posición con el láser encima, el error de medición, las correcciones, la forma del pozo. Sus «mW» son los del display |
| `figuras/ocho barridos/` y `figuras/método/` | las figuras de cada notebook, en SVG |
| `resultados/` | tablas: `ocho_*.csv` (resultado) y `potencia_*.csv` (método); `LEEME.txt` explica las columnas |
| `DUDAS.md` | preguntas para responder y decisiones que se tomaron |

Las grabaciones: `LASER 5` (barrido 1) y `BARRIDO POTENCIAS 2` a `8`. Las funciones están en `../codigo/potencia.py`. La tabla de potencias de cada escalón está en `../datos/cuaderno/potencias_cuaderno.csv` y la calibración en `../datos/cuaderno/calibracion_potencia.csv`.

## La potencia

El 19/11 se midió la potencia real y el display no la da: «30 mW» son 4.80 mW, «25» son 4.27, «20» son 3.15, «15» son 1.83 y «5» son 0.004 mW. A los barridos 1 a 6 (anotados con el display) se les asigna la potencia real con ± 0.5 mW. En los barridos 7 y 8 se anotó la potencia ya medida.

## Cómo se hizo

1. Se mide la posición de la partícula atrapada en cada cuadro, corrigiendo el punto oscuro que deja el láser.
2. Para cada escalón, la varianza de la posición da la rigidez: $\kappa = k_BT/\mathrm{Var}(x)$, después de restar el error de medición y de corregir el promedio que hace la cámara durante la exposición (30.5 ms el 19/11; no anotada el 12/11, se supone medio período).
3. Se miró el video de cada grabación para ver si la partícula está sola y en foco. Se usan 26 de los 36 escalones.
4. Se ajusta una recta por el origen, $\kappa = b\,P$, teniendo en cuenta el error de la potencia.

Figuras para empezar: `figuras/ocho barridos/rigidez_potencia_x.svg` (el resultado), `calibracion_potencia.svg` (display contra potencia real), `recortes_barrido_7.svg` (la partícula que se le pega a la atrapada) y `potencial_BARRIDO_POTENCIAS_8_x.svg` (las parábolas). Del método: `figuras/método/histograma_posiciones_LASER_5_x_30mW.svg` (el histograma del que sale una parábola).

## Resultados

- **La rigidez es proporcional a la potencia real.** El 19/11, con 16 escalones limpios: $\kappa_x = (0.228 \pm 0.011)$ y $\kappa_y = (0.276 \pm 0.013)$ pN/µm por mW, entre 1.56 y 4.85 mW. La recta por el origen describe bien los datos.
- **La calibración del display sirve.** Los barridos 4 a 6 (display convertido) dan $0.235 \pm 0.017$ en $x$; los 7 y 8 (potencia medida), $0.222 \pm 0.015$.
- **A potencia máxima** (4.8 mW reales), $\kappa_x$ es de 1.0 a 1.5 pN/µm el 19/11 y de 2.0 a 2.1 pN/µm el 12/11.
- **El 12/11 la trampa era 1.3 a 1.7 veces más rígida a igual potencia calibrada** ($0.40 \pm 0.04$ pN/µm por mW en $x$). No es del método: ese día la rigidez medida con las patadas (objetivo 2) dio lo mismo, 1.9 ± 0.2 pN/µm. La causa queda abierta.
- **La medición se repite dentro de un 18 %** (la misma partícula en los barridos 4, 5 y 6). El error calculado en cada grabación es 12 %; la diferencia se sumó al error de cada punto.
- **El arrastre $\gamma$ es el mismo los dos días:** 5.9 ± 0.5 y 5.5 ± 0.4 ×10⁻⁸ N·s/m.
- **Forma del pozo cerca del centro:** parabólica hasta donde llega el movimiento térmico (2 a 5 $k_BT$), los dos días.
- **Con «5» en el display (0.004 mW) la partícula se suelta** en los barridos 4, 5 y 6.

Rigidez $\kappa_x$ de los barridos limpios, en pN/µm (el error incluye el 13 % de repetición):

| Barrido | 4.8 mW | 4.1–4.3 | 3.2–3.6 | 2.3 | 1.6–1.8 | 0.4 |
|---|---|---|---|---|---|---|
| 1 (`LASER 5`, 12/11) | 2.03 ± 0.33 | 1.46 ± 0.25 | 1.17 ± 0.19 | | 0.61 ± 0.11 | 0.13 ± 0.03 |
| 2 (12/11) | 2.11 ± 0.73 | 1.68 ± 0.43 | 1.46 ± 0.39 | | 1.15 ± 0.26 | 0.64 ± 0.11 |
| 4 (19/11) | 1.50 ± 0.22 | 1.09 ± 0.20 | 0.76 ± 0.12 | | 0.32 ± 0.06 | |
| 5 (19/11) | 1.03 ± 0.17 | 0.70 ± 0.13 | | | | |
| 6 (19/11) | 1.19 ± 0.19 | 1.36 ± 0.21 | 0.69 ± 0.13 | | 0.43 ± 0.09 | |
| 7 (19/11) | 1.01 ± 0.18 | | | | | |
| 8 (19/11) | 1.05 ± 0.17 | 1.07 ± 0.21 | 1.00 ± 0.15 | 0.54 ± 0.09 | 0.28 ± 0.05 | |

(Los valores exactos, con los dos ejes, están en `resultados/ocho_escalones.csv`.)

## Lo que no se usa y por qué

- **Barrido 3, entero:** tiene motas pegadas a la partícula toda la grabación, y entre 150 y 180 s sigue retenida cuando con «5» en el display tendría que soltarse.
- **Barrido 5, escalones de «20» y «15»:** la imagen se desenfoca por momentos.
- **Barrido 7, desde el segundo escalón:** a los 76 s llega otra partícula más chica y queda pegada a la atrapada hasta los 170 s. De los dos barridos con potencia medida, el que vale entero es el 8.

## Pendiente

- **¿Por qué el 12/11 la trampa era más rígida?** ¿Se tocó la alineación o el foco entre los dos días? ¿La calibración del 19/11 vale para el 12/11?
- **La exposición del 12/11** sigue sin conocerse.
- **El error de la potencia medida** (barridos 7 y 8): se supuso ± 0.05 mW.
- **El «7» del display** (barridos 1 y 2) sale de interpolar la calibración entre «5» y «15»: 0.37 ± 0.5 mW.
- El detalle está en `DUDAS.md`.
