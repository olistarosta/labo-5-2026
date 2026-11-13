# Dudas y decisiones — barrido de potencia

Acompaña a `PINZAS 4 - barrido de potencia.ipynb` (grabaciones `LASER 5` y `BARRIDO POTENCIAS 2`, 12/11/2026). Cada punto dice qué se supuso o se decidió, cuánto cambia el resultado y qué habría que medir o preguntar para cerrarlo. Las **preguntas para Lenny** van primero.

Resultados a los que se refieren los impactos ($W = dt/2$, error estadístico + ruido):

| | 30 mW | 25 mW | 20 mW | 15 mW | 7 mW |
|---|---|---|---|---|---|
| `LASER 5`, $\kappa_x$ [pN/µm] | 2.03 ± 0.20 | 1.46 ± 0.17 | 1.17 ± 0.11 | 0.61 ± 0.07 | 0.13 ± 0.02 |
| `LASER 5`, $\kappa_y$ | 1.71 ± 0.12 | 1.19 ± 0.14 | 0.89 ± 0.08 | 0.50 ± 0.07 | 0.24 ± 0.02 |
| `BARRIDO POTENCIAS 2`, $\kappa_x$ | 2.1 ± 0.7 | 1.7 ± 0.4 | 1.5 ± 0.3 (21.5 mW) | 1.15 ± 0.21 | 0.64 ± 0.08 |
| `BARRIDO POTENCIAS 2`, $\kappa_y$ | 1.59 ± 0.25 | 1.80 ± 0.31 | 1.47 ± 0.22 | 1.14 ± 0.13 | 0.58 ± 0.07 |

$\gamma = (5.9 \pm 0.5)\times10^{-8}$ N·s/m (difusión de la partícula suelta) y $(5.4 \pm 0.5)\times10^{-8}$ ($\kappa\,\tau_c$ a 15 y 7 mW).

## Preguntas para Lenny

### 1. ¿Qué exposición tenía la cámara?

Es **la mayor incerteza de $\kappa$ a potencia alta**. Los TIFF no la traen. La cámara promedia la posición durante la exposición $W$ y eso achica la varianza medida cuando $W$ no es mucho menor que $\tau_c$ (≈ 30 ms a 30 mW).

- **Impacto:** $\kappa_x$ a 30 mW va de 3.0 ($W = 0$) a 1.7 pN/µm ($W = dt$ = 76–80 ms); a 15 mW, de 0.69 a 0.55; a 7 mW, menos de 6 %.
- **Lo que se sabe:** (a) en `BARRIDO 2`, por el largo de los trazos de las partículas durante una patada, $W/dt = 0.40 \pm 0.10$ (unos 40 ms a 10 cuadros/s; `po.exposicion_por_estirado`), pero esa grabación tiene otra configuración de cámara; (b) con $W = dt/2$, $\kappa_x/\gamma$ a 30 mW (35 s⁻¹) coincide con el del arrastre (38 ± 4 s⁻¹, PINZAS 3); (c) con $W = dt/2$, $\kappa\,\tau_c$ coincide con la difusión. Los tres apuntan a $W \approx dt/2$. La forma del espectro prefiere $W$ corta en `LASER 5`, pero ese test confunde exposición con ruido.
- **Decisión:** $W = dt/2$ central, $[0, dt]$ como sistemático.
- **Para cerrarlo:** ¿se acuerdan del valor de *Exposure* en ThorCam, o quedó en alguna captura de pantalla? ¿Cambiaron la exposición entre las grabaciones de las 17 h, la de las 18:22 (`LASER 5`) y las de las 19 h? (Los cuadros por segundo y el balance de color son distintos en los tres grupos.) Anotarla la próxima vez.

### 2. ¿Cómo se leyó la potencia del láser? ¿Es potencia medida o un valor del controlador?

La luz roja que llega a la cámara **no cambia entre 0 y 150 s** (dentro de 5–12 % en el núcleo), mientras el cuaderno dice 30 → 7 mW. No es saturación de los píxeles (máximo 110–230 cuentas de 255; el halo, que vale pocas cuentas, tampoco cambia). O el rojo que ve la cámara no es proporcional a la potencia (hipótesis: el filtro bloquea la línea del láser y pasa la emisión espontánea del diodo, que por encima del umbral no crece), o la potencia en la muestra no es la del cuaderno.

- **Impacto:** no se pudo verificar ningún escalón de 30 a 7 mW ni estimar los "< 7 mW". Todo el eje de potencia de las figuras es el del cuaderno.
- **Indicios de que la potencia sí cambia:** la dispersión de la posición crece escalón a escalón y, en `BARRIDO POTENCIAS 2`, $\kappa$ sale proporcional a la potencia del cuaderno.
- **Para cerrarlo:** ¿los mW son de un medidor de potencia, de una tabla corriente → potencia o del display del controlador? ¿Qué filtro hay delante de la cámara? Si se repite: medir la potencia con el medidor en cada escalón y anotar la corriente.

### 3. ¿Por qué `LASER 5` no es lineal y `BARRIDO POTENCIAS 2` sí?

A 30 mW las dos dan lo mismo. A 15 mW `LASER 5` da la mitad de $\kappa$ que `BARRIDO POTENCIAS 2`, y a 7 mW, entre 2.5 y 5 veces menos. Los dos estimadores independientes (varianza y tiempo de relajación) coinciden: es la trampa.

- **Posibles causas, sin poder elegir:** los valores de potencia de una de las dos no son los reales; otra altura de la partícula sobre el vidrio; en `LASER 5` el flujo era algo más rápido (0.26 contra 0.23 µm/s) y a 7 mW corre el equilibrio 0.1 µm; otra alineación o foco (entre las dos se cambió la configuración de la cámara, ¿se tocó algo más?).
- **Impacto:** la recta $\kappa(P)$ de `LASER 5` cruza el cero en 3–9 mW según qué puntos entren; la de `BARRIDO POTENCIAS 2`, en −4 ± 3 mW.
- **Para cerrarlo:** ¿se acuerdan de algo distinto entre las dos (enfoque, altura, otra partícula)? ¿En `LASER 5` los escalones bajos fueron realmente 15 y 7 mW?

### 4. ¿En `LASER 5` el láser se apagó a los 180 s?

El cuaderno dice "< 7 mW" de 150 a 210 s y "0" desde 210 s. En el video no queda nada de rojo desde los 181 s.

- **Impacto:** ninguno en $\kappa$. El escalón 7 de `LASER 5` es láser apagado (o invisible), no "< 7 mW".

### 5. ¿Qué temperatura había?

No se anotó. Se supuso 23 ± 3 °C.

- **Impacto:** ± 1 % en $k_BT$ y por lo tanto en $\kappa$. La viscosidad no se usa (se usa el $\gamma$ medido); solo cambia ± 9 % la comparación con Stokes.

### 6. ¿Qué muestra era?

La receta del cuaderno (60 ml de agua + 65 ml de glicerina, 58 % en masa, 8.4 mPa·s) da un arrastre **4 veces mayor** que el medido de dos maneras con la misma partícula. Con 24 % de glicerina (1.83 mPa·s) y $R$ = 1.5 µm, medido/Stokes = 1.05–1.15. Es la misma conclusión del análisis del browniano (`DUDAS browniano.md`, punto 1).

- **Impacto:** ninguno en $\kappa$ ni en $\tau_c$. Solo en la comparación con Stokes y en el factor de Faxén (1.05–1.15 con 24 %: centro a ~8 µm del vidrio, o radio 10 % mayor; no se separan).
- **Para cerrarlo:** qué mezcla se usó de verdad y a qué altura se atrapó la partícula.

## Decisiones que se tomaron

### 7. Las posiciones de `atrapada.csv` no sirven para la varianza: se midieron de nuevo

La cámara entrega RGB procesado y el rojo del láser se **resta** del verde (0.30 cuentas de verde por cuenta de rojo en `LASER 5`, 0.14 en `BARRIDO POTENCIAS 2`; regresión temporal píxel a píxel). En el verde queda un punto oscuro en el láser, que la simetría radial usa para localizar. Se volvió a medir la posición sobre verde + $c_G$·(rojo − rojo del fondo), sin los píxeles a menos de 8 px del láser (`xc_um`, `yc_um` de `potencia_cuadros.csv`).

- **Impacto:** el desvío estándar cambia de −20 % (`LASER 5`, $y$) a +60 % (`BARRIDO POTENCIAS 2`, $x$); $\kappa$ hasta un factor 2–3. Entre tapar o no el núcleo después de corregir: menos de 15 % en $\sigma$. Fuera del láser los estimadores coinciden (pendiente 0.99–1.01 entre corregido con y sin tapar; 0.95–1.03 contra la conversión).
- **Lo que queda abierto:** en `LASER 5` el centro del perfil corregido queda más claro que el de la partícula suelta. O la corrección no es exacta o la partícula atrapada enfoca distinto. El núcleo tapado lo cubre, pero no se demostró que el estimador `c` tenga ganancia exactamente 1 con el láser prendido. El canal azul, que podría decirlo, es 5–10 veces más ruidoso.
- **Para el análisis por arrastre (PINZAS 3):** usa `atrapada.csv`. Para desplazamientos grandes (la partícula sale del láser) no hay problema; cerca del centro, la misma contaminación.
- **Para cerrarlo de raíz:** grabar *Unprocessed* (Bayer crudo), o poner un filtro que saque el láser, o grabar en monocromo.

### 8. Ruido de localización: el de píxel, medido con una plantilla

Imagen promedio de la partícula + ruido real de la cámara (diferencia de cuadros seguidos con la trampa vacía) → se localiza → desvío. Da 19–22 nm (`LASER 5`) y 32–39 nm (`BARRIDO POTENCIAS 2`). El método predice dentro de ± 15 % el desvío medido de dos manchas de polvo quietas. Se usa el promedio de los escalones 1 y 2 de cada grabación (en los escalones bajos la plantilla está emborronada por el movimiento y el ruido sale 10 % mayor).

- **Es una cota inferior**: no incluye lo que agregue la mancha del láser, que cambia de forma cuadro a cuadro. Si el ruido real es mayor, $\kappa$ es mayor.
- **Impacto:** el ruido es 20–26 % de la varianza a 30 mW en `LASER 5` y 36–55 % en `BARRIDO POTENCIAS 2`. Con ± 20 % en la varianza del ruido, $\kappa$ a 30 mW cambia ± 4–7 % en `LASER 5` y ± 11–24 % en `BARRIDO POTENCIAS 2` (columna `err_kappa_ruido`).
- **Chequeo independiente:** en `BARRIDO POTENCIAS 2` la ordenada de Var contra $1/P$ (1680 ± 300 nm² en $x$, 1250 ± 280 en $y$) coincide con el ruido de la plantilla (1540 y 1020 nm²).
- **No se usaron partículas pegadas** de `BARRIDO 1/2` (como se había pensado): son otra configuración de cámara y otro estimador. En las grabaciones de potencia no hay partículas pegadas a la vista (lo que queda quieto en la mediana es polvo del sensor).

### 9. Tramo usable de cada escalón

El rojo lejano cae a los 151–152 s, 181 s y 211 s: el cambio real es ~1.5 s después de cada múltiplo de 30 s. Se supuso el mismo retraso en los cambios de 30, 60, 90 y 120 s, que no se ven. Se usa de $30k + 5$ s a $30k + 30.5$ s.

- **Impacto:** si algún cambio fue 3 s más tarde que lo supuesto, ese tramo tiene hasta 6 % de cuadros del escalón anterior (más rígido): $\kappa$ sobreestimada unos pocos %.

### 10. Solo los escalones 1 a 5 son de partícula atrapada

Desde los 151 s la partícula se va a más de 1 µm. `en_trampa` de `atrapada.csv` sigue en True (pide una partícula a menos de 45 px), pero ya no está en equilibrio en el pozo. De esos tramos se sacó solo el coeficiente de difusión y la descripción de la sección 10.

### 11. Deriva lenta: se resta una parábola por tramo

El equilibrio se mueve despacio (hasta 50 nm en 25 s). Sin restarlo, la varianza a 30 mW sube hasta 25 %.

- **Impacto (sintéticos):** la parábola se lleva varianza verdadera cuando $\tau_c$ es largo: $\kappa$ sobreestimada 4 % a 15 mW y 10 % a 7 mW. No se corrige (la dispersión de una medición es 9–17 %). El espectro se ajusta desde 0.2 Hz y no tiene ese sesgo.
- Queda una componente lenta que la parábola no saca: la autocorrelación no es una sola exponencial, y su $\tau_c$ sale más largo que el del espectro.

### 12. Qué quiere decir "$\tau_c$ resuelto"

$\tau_c$ del espectro mayor que un cuadro, error menor que 40 % y ruido por debajo del 30 % de la varianza. Con eso quedan 7 de 20 tramos: `LASER 5` a 15 y 7 mW (y 20 mW en $y$) y `BARRIDO POTENCIAS 2` a 7 mW. El tercer criterio saca el eje $x$ de `BARRIDO POTENCIAS 2` a 30, 21.5 y 15 mW, donde el espectro da $\tau_c$ de 0.1–0.2 s que no son creíbles ($\kappa\,\tau_c$ = 14–26 ×10⁻⁸, incompatible con el eje $y$ y con la difusión).

- **Impacto:** con otro criterio cambia qué puntos aparecen llenos en las figuras de $\tau_c$ y el promedio de $\kappa\,\tau_c$ (5.4 → 5.6 ×10⁻⁸ si entran esos tres).

### 13. $\gamma$ de referencia

Se usa el de la difusión de la misma partícula después de soltarse, $(5.9 \pm 0.5)\times10^{-8}$ N·s/m (varianza de desplazamientos a 1–4 cuadros, dentro de trozos de 6 s). Solo entra en la corrección por exposición y en la comparación de $\tau_c$.

- **Impacto:** ± 10 % en $\gamma$ cambia $\kappa$ a 30 mW ± 3 % (con $W = dt/2$).
- En `BARRIDO POTENCIAS 2` la partícula todavía está algo retenida entre 152 y 210 s; un confinamiento débil casi no cambia la pendiente a 1–4 cuadros, pero $D$ podría estar algo subestimado. $D_x$ y $D_y$ difieren (0.08–0.09 contra 0.05–0.06 µm²/s): con 230 y 650 cuadros, es 1.5–2 desvíos.

### 14. La platina "quieta" y el flujo

Las partículas libres derivan a 0.23–0.26 µm/s hacia −x, +y. No se distingue un flujo del líquido de una deriva de la platina (no hay partículas pegadas visibles). Para la partícula atrapada es una fuerza de ~0.015 pN.

- **Impacto:** corre el equilibrio ~6 nm a 30 mW y ~85 nm a 7 mW en `LASER 5`; si el flujo fluctúa, agrega varianza lenta en los escalones bajos (no cuantificado). Por debajo de 7 mW se lleva a la partícula.
- El corrimiento del centro medido en `BARRIDO POTENCIAS 2` (+50 nm en $x$, +60 nm en $y$ de 30 a 7 mW) no se explica con esto y no se distingue de una deriva.

### 15. Potencia de "20 a 23 mW"

Se tomó 21.5 mW (así está en `potencias_cuaderno.csv`). **Impacto:** ± 1.5 mW en ese punto; no cambia las rectas.

### 16. Figuras: lo que la guía no cubre

- Con tres potencias en un gráfico (potencial, autocorrelación, espectro, perfil radial) se usaron los dos colores de la guía y el azul de la paleta de seaborn, cada serie con su marcador.
- Las líneas blancas verticales de `rojo_tiempo` y `desvio_tiempo` marcan los múltiplos de 30 s del cuaderno.
- En las figuras de $\tau_c$ el marcador hueco es "no resuelto" (el criterio de "datos excluidos").
- Ejes logarítmicos en `desvio_tiempo`, `tau_potencia_*` y `espectro`, con los ticks en valores redondos.
- En `validacion_sinteticos` la línea horizontal es el valor 1 (estimado = verdadero).
- En `perfil_radial`, `exposicion_verosimilitud` y `potencia_baja_distancia` los puntos no llevan barras de error (perfil promedio, verosimilitud y posiciones cuadro a cuadro).
- En `varianza_inversa_potencia` el punto de 7 mW de `LASER 5` (30 000 nm²) queda fuera de la escala.

## Para el otro Claude (no son dudas de física)

- `platina.csv` de `BARRIDO 2`: los cuadros 1510 y 1511 figuran quietos (dx ≈ 0, corr 0.04–0.06) pero en la imagen las partículas siguen estiradas: la correlación falla cuando el estirado es grande y la patada dura más de lo que dice la tabla.
- `correr_tests.py` con el Python global: 134 pasan y 4 fallan por `ModuleNotFoundError` (falta cv2: son tests para el Python de TouchDesigner). No se tocó `pinzas/`.
- `potencia.tabla_cuadros(recalcular=True)` relee los TIFF de `C:\datos_pinzas\GRABACIONES` (menos de un minuto con los archivos en caché del disco, varios si no).
