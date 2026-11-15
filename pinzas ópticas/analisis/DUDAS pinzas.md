# Dudas y decisiones — pinzas ópticas (platina, fuerza y pozo)

Acompaña a `PINZAS 1`, `PINZAS 2` y `PINZAS 3` (12/11/2026). Las dudas del barrido de potencia están en `DUDAS potencia.md` y las del browniano en `DUDAS browniano.md`. Primero van las **preguntas para Lenny**, por orden de cuánto cambian los resultados; después, las decisiones que se tomaron sin preguntar.

Resultados a los que se refieren los impactos (30 mW, partícula de tamaño nominal):

| Magnitud | Valor |
|---|---|
| $\kappa/\gamma$ | 38 ± 4 s⁻¹ ($\tau_c$ = 26 ms, $f_c$ = 6.1 Hz) |
| Ancho del pozo | 1.14 ± 0.09 µm |
| Fuerza máxima / $\gamma$ | 26 ± 2 µm/s |
| Rigidez | 1.9 ± 0.2 pN/µm |
| Fuerza máxima | 1.3 ± 0.1 pN |
| Profundidad | 600 ± 80 $k_BT$ |

## Preguntas para Lenny

### 1. ¿Qué mezcla tenía la muestra de la pinza?

El cuaderno dice «60 ml agua y 65 ml glicerina; usamos 27 micro l». Eso es 58 % de glicerina en masa (8.4 mPa·s a 23 °C), y una partícula de 3 µm difundiría a 0.017 µm²/s.

Lo que se mide en las partículas libres **dentro de las mismas grabaciones** de la pinza (`LASER 5` y `BARRIDO POTENCIAS 2`, 22 partículas de tamaño nominal) es 0.082 ± 0.004 µm²/s, y 0.079 en `BROWNIANO 2`. Es 4.7 veces más rápido: corresponde a una mezcla con alrededor de 24 % de glicerina.

- ¿Los volúmenes eran otros (por ejemplo 0.65 ml de glicerina, o µl en vez de ml)?
- ¿O se usó la muestra de «0.25 ml glicerina, 1 ml agua» anotada para BROWNIANO 1?

**Impacto:** ninguno en la forma del pozo, $\tau_c$ ni $f_c$. En pN, **un factor 4.7**: se usó el arrastre medido con las partículas libres. Si la receta del cuaderno fuera la correcta y la difusión medida estuviera mal, la rigidez sería 9 pN/µm y la fuerza máxima 6 pN.

**Lo que dicen los datos:** tres medidas independientes del arrastre coinciden entre sí y no con la receta (en 10⁻⁸ N·s/m): partículas libres 5.0; la partícula atrapada ya casi suelta 5.5 a 5.9; rigidez térmica dividida por $\kappa/\gamma$ del arrastre 5.3. La receta daría 23.6. Así que la respuesta no cambia los resultados: hace falta para saber qué composición de la muestra se anota en el informe.

### 2. ¿Qué potencia tenía el láser en las grabaciones con patadas?

Se supuso la máxima (30 mW) en `LASER NE`, `BARRIDO 1`, `BARRIDO 2`, `MUEVO A MANO` y los tres `BARRIDO STEPSIZE`, porque el cuaderno no dice otra cosa.

**Impacto:** todos los valores de fuerza son «a la potencia que había». Si no era 30 mW, cambia la comparación con el barrido de potencia.

### 3. ¿Los 30 mW son a la salida del láser o en la muestra?

**Impacto:** no cambia ningún resultado; sí la eficiencia de la trampa que se informe (fuerza por unidad de potencia).

### 4. ¿Qué tamaños de partículas había?

En `BARRIDO 2` tres de las cuatro partículas que pasaron por la trampa tienen un radio aparente de 1.05 µm, contra 1.8 a 2.0 µm del resto. Con esas, la partícula se escapa con patadas más chicas.

**Impacto:** `BARRIDO 2` se dejó fuera del ajuste del pozo. Con el tamaño nominal de cada tipo se podría comparar la fuerza entre tamaños.

### 5. ¿Las patadas de `LASER NE` son el viaje de 89 µm?

El cuaderno anota dos patadas sin hora (a = 49 y 100 µm/s², paso 200 µm). En el video hay un solo movimiento grande: a los 55 s la platina va 89 µm en vertical y vuelve 86 µm tres segundos después. El resto es una deriva lenta de toda la muestra.

**Impacto:** `LASER NE` no se usó para la fuerza.

### 6. ¿`BARRIDO STEPSIZE 3` es el barrido vertical?

No tiene renglones en el cuaderno. Tiene 24 movimientos verticales y 16 horizontales, con la misma partícula.

**Impacto:** se usó para comparar el pozo en vertical y en horizontal. Si los movimientos tenían otro propósito, la comparación vale igual.

### 7. ¿Qué exposición tenía la cámara?

No está en los archivos. Todo indica que era larga, cercana al tiempo entre cuadros (80 a 100 ms).

**Impacto:** en el arrastre, poco (ya está incluido en el sesgo de hasta 8 %). En el análisis térmico, mucho: ver `DUDAS potencia.md`.

## Decisiones que se tomaron

### A. El movimiento de la platina se mide en el video, no con el motor

Por correlación de cuadros consecutivos (`platina.py`), validada con las partículas seguidas antes y después (diferencia mediana 0.5 µm). El desplazamiento medido es en mediana el 5 % del paso pedido y no guarda relación con él.

### B. Cuadros por segundo por fechas de archivo

10.02, 13.13 y 12.50 según la grabación (`datos/catalogo.csv`). Los TIFF no traen tiempo por cuadro.

**Impacto:** un error del 1 % en los cuadros por segundo es 1 % en todas las velocidades y en $\kappa/\gamma$.

### C. La partícula atrapada se sigue con una ventana propia

El seguimiento general pierde a la partícula cuando se mueve más de 10 px por cuadro, que es justo cuando interesa. `trampa_roi.py` detecta en cada cuadro, sin enlazar, en una ventana de ±128 px alrededor del láser, y `pinzas_fuerza.seguir` arma la trayectoria: una detección a menos de 8 px del equilibrio con mancha roja es la atrapada; si no, se sigue a la que venía.

**Riesgo:** después de un escape puede seguir a otra partícula que quedó cerca. Por eso `BARRIDO STEPSIZE` (todas las patadas con escape) da un ajuste no confiable y se excluyó.

### D. El fondo no borra a la partícula atrapada

Con la mediana temporal, una partícula quieta en la imagen pasa a ser fondo. Se agregó al bloque de fondo (`b2_fondo.py`) la opción de rellenar un disco de 45 px alrededor del láser con una interpolación suave. El fondo espacial no alcanza: la iluminación tiene manchas del tamaño de las partículas.

### E. Pozo gaussiano

Es el modelo más simple con un máximo de fuerza y alcance finito. Describe bien los datos hasta unos 3 µm; más lejos la trampa todavía atrae un poco (la partícula volvió desde 5.6 µm), y eso el modelo no lo tiene.

### F. Qué patadas entran en el ajuste

Horizontales de `BARRIDO 1` y `BARRIDO STEPSIZE 2` (23 patadas), dos grabaciones independientes con partículas de tamaño nominal. `BARRIDO STEPSIZE 3` se ajusta aparte (vertical y horizontal) y coincide dentro de su error. `MUEVO A MANO` da valores parecidos pero con movimientos irregulares.

### G. Arrastre de las partículas libres para pasar a pN

$\gamma = k_BT/D$ con las 22 partículas libres de tamaño nominal de los barridos de potencia. La partícula atrapada, cuando casi no la retiene la trampa, da 11 % más. No se corrigió por la diferencia de altura sobre el vidrio entre libres y atrapada, que no se conoce.

### H. Temperatura

23 °C para $k_BT$. ± 3 °C es ± 1 % en $k_BT$; no pesa (la viscosidad no se usa).

### I. Radio de la partícula atrapada

`medir_radios.py` toma el punto oscuro del centro de la partícula atrapada como si fuera el borde y da 0.4 a 0.7 µm. Releyendo el perfil radial y buscando el anillo oscuro fuera de los primeros 8 px sale 1.8 a 2.0 µm, igual que las libres. Es un radio aparente: depende del foco.

### J. El punto oscuro del láser

La cámara le resta rojo al verde (lo encontró el análisis del barrido de potencia): la partícula atrapada tiene un punto oscuro en el láser, y la localización lo usa. Cuando la partícula se aparta poco, la posición medida queda corrida hacia el láser entre 0.1 y 0.25 µm.

**Control:** se repitió la detección de `BARRIDO 1` y `BARRIDO STEPSIZE 2` rellenando un disco de 11 px en el láser (`trampa_roi.py --tapar=11`, salida en `trampa_detecciones_tapado.csv`). El ajuste del pozo cambia −9 % en $\kappa/\gamma$, +6 % en el ancho y −4 % en la fuerza máxima, menos que el error estadístico. Se dejó el ajuste sin tapar como resultado (la posición quieta es menos ruidosa) y esa diferencia como sistemático.

**Ojo:** para el movimiento térmico de la partícula quieta el efecto sí es grande (hasta 2 a 3 veces en la rigidez). Ahí hay que usar `xc_um`, `yc_um` de `datos/pinzas/potencia_cuadros.csv`, no `atrapada.csv`.

### K. Escapes que el pozo no explica

Tres escapes con patadas chicas: `BARRIDO 1` n = 7 (5.4 µm; el cuaderno también anota escape), y `BARRIDO STEPSIZE` n = 8 y n = 15. No se investigó la causa (otra partícula que choca, una partícula distinta, un error de seguimiento). Conviene mirarlos en el MP4.
