# Dudas y decisiones — movimiento browniano

Acompaña a `BROWNIANO comparación.ipynb` (12/11/2026). Cada punto dice qué se supuso o se decidió, por qué, y cuánto cambia el resultado. Las **preguntas para Lenny** van primero; después las decisiones que se tomaron sin preguntar.

Resultados a los que se refieren los impactos ($D$ de las beads solas):

| Medición | $D$ [µm²/s] |
|---|---|
| BROWNIANO 1 | 0.056 ± 0.004 |
| BROWNIANO 2 | 0.079 ± 0.007 |
| GROS BROWNIANO | 0.169 ± 0.007 |
| CHICAS BROWNIANO | 0.102 ± 0.010 |

## Preguntas para Lenny

### 1. ¿Qué muestra se grabó en BROWNIANO 2?

El cuaderno dice 60 ml de agua + 65 ml de glicerina (58 % de glicerina en masa, 8.4 mPa·s a 23 °C). Con esa mezcla y beads de 3 µm, Stokes–Einstein da $D$ = 0.017 µm²/s. Se mide 0.079: **4.6 veces más**. La pared solo puede bajar $D$, y para llegar a ese valor con 58 % haría falta una muestra a unos 70 °C.

El valor medido es el de una mezcla de **24 % de glicerina** (0.079 µm²/s), que es justo la receta anotada para BROWNIANO 1 (0.25 ml de glicerina + 1 ml de agua).

- ¿Puede ser que BROWNIANO 2 se haya grabado con la muestra de 0.25 ml + 1 ml?
- ¿Los 60 ml + 65 ml son de otra muestra (por ejemplo la de la pinza), o son otras unidades?

**Impacto:** la comparación con Stokes–Einstein de BROWNIANO 2 no se puede cerrar hasta saber la mezcla. El valor de $D$ no cambia.

### 2. ¿La receta de BROWNIANO 1 es la del cuaderno?

"0.25 ml glicerina, 1 ml agua?" tiene un signo de pregunta. Con esa receta (24 % en masa, 1.83 mPa·s) se espera 0.079 µm²/s y se mide 0.056: **0.71 del teórico**. Tres explicaciones que los datos no distinguen:

- las beads están cerca del vidrio: con la superficie de la bead a 1.4 µm de la pared, Faxén da ese factor;
- la mezcla tenía 34 % de glicerina en vez de 24 %;
- las beads miden más de 3 µm (el radio aparente es 1.8 µm; con ese radio el teórico baja a 0.066 y el cociente sube a 0.85).

**Impacto:** hasta 30 % en la comparación con la teoría. Para decidir haría falta la altura de las beads sobre el vidrio (o grabar a varias alturas) y la receta confirmada.

### 3. ¿Qué temperatura había?

No se anotó. Se supuso 23 ± 3 °C.

**Impacto:** en el $D$ teórico, ± 9 % con 24 % de glicerina y ± 8 % en agua. No toca el $D$ medido.

### 4. ¿Qué tamaños de beads había en la muestra?

En BROWNIANO 2 hay partículas de 0.7 a 2.1 µm de radio aparente (7 del tamaño nominal y 12 más chicas, entre las de más de 20 s). En BROWNIANO 1 hay una chica.

- ¿Se mezclaron beads de dos tamaños nominales (por ejemplo 1 y 3 µm)? ¿Cuáles?

**Impacto:** con los radios nominales se podría calibrar el radio aparente y comparar $D \propto 1/R$ en serio. Hoy sale $D \propto R_{ap}^{-1.9 \pm 0.3}$, más empinado que Stokes–Einstein, y no se sabe si es por el radio aparente o por la pared.

### 5. ¿Se puede conseguir la escala, el fps y la exposición de los videos de los otros grupos?

- **Escala.** Para GROS y CHICAS se usó la nuestra (94.61 nm/px; en CHICAS × 1024/1080 porque el MP4 está agrandado a 1440 × 1080 con bandas negras). Es una suposición. Chequeo: el radio aparente de sus beads da 1.91 y 1.90 µm contra 1.80 y 1.82 de las nuestras, 5 % más. Si las beads son las mismas, la escala está bien dentro de 5 %. **Impacto:** cada 1 % de error en la escala son 2 % en $D$; con 5 %, 10 %.
- **Tiempo de CHICAS.** Ver la decisión 7. **Impacto:** hasta −12 % en $D$.
- **Exposición de GROS.** La ordenada negativa del MSD y la correlación positiva entre pasos (+0.13) equivalen a una exposición de 0.6 veces el período (unos 40 ms). Si se supiera la exposición real, se podría confirmar. **Impacto:** ninguno en $D$ (sale de la pendiente).
- **El video original de CHICAS.** El MP4 pesa 1 MB y está muy comprimido. Si existe el original (AVI o TIFF), conviene reconvertir.

### 6. ¿A qué altura estaban las beads?

En BROWNIANO 1 hay tres partículas desenfocadas (24, 26, 27) que derivan cuatro veces más rápido que las beads en foco (0.7 contra 0.15 µm/s) y difunden 1.7 veces más (0.096 contra 0.056 µm²/s). Es lo que se espera de partículas más lejos del vidrio (el flujo es más rápido y la pared frena menos), pero también podrían ser partículas más chicas vistas fuera de foco.

**Impacto:** si son beads del mismo tamaño lejos de la pared, su $D$ (0.096) sería el "libre", y el de las beads en foco estaría frenado por la pared a 0.58. No está demostrado.

## Decisiones que se tomaron

### 7. CHICAS: el tiempo de cada cuadro

El MP4 está a 30 cuadros/s pero repite cuadros. El importador deja 229 cuadros distintos, con el tiempo del contenedor: intervalos de 1 o 2 cuadros del contenedor (patrón 1, 1, 2) y 17 saltos más largos. La varianza del paso es casi la misma con intervalo 1 que con intervalo 2 (1.09 y 1.24 px²) y el doble en los saltos largos (2.5 px²). Se concluyó que la cámara iba pareja a **22.2 cuadros/s** y que en los saltos largos falta un cuadro. `tiempo_regular` pone un período constante entre saltos largos.

**Impacto:** si no faltara ningún cuadro serían 19.5 cuadros/s y $D$ sería 12 % menor. Con el tiempo del contenedor tal cual, el MSD a los primeros retardos mezcla intervalos distintos.

### 8. CHICAS: parámetros de conversión

`convertir.py ... --nombre "CHICAS BROWNIANO" --set b5_realce.Radio=22 b6_deteccion.Umbral=6`, sin escala (se pasa al cargar).

- Con los parámetros por defecto no se detectaba casi nada (0 posiciones en el primer cuadro). Las beads miden unos 44 px de diámetro en el video: Radio 22.
- El z de las beads después de sacar el fondo llega a 10–19, no a 20–40: Umbral 6 (60 a 80 detecciones por cuadro, de las cuales quedan 20 a la vez en trayectorias largas).
- Se probó el fondo `espacial`: las bandas negras del video y el poco contraste dejan z < 5, no sirve. Quedó la **mediana**, aunque las beads se corren poco (una o dos veces su diámetro en todo el video) y la mediana guarda una sombra de cada una. Esa sombra puede correr un poco la posición.
- Las trayectorias pegadas a las bandas negras (x < 115 px o x > 1325 px) se descartan en el análisis.

**Impacto:** el video da un $D$ razonable porque el MSD promedia, pero los pasos de un cuadro no son gaussianos (curtosis 7.9, anisotropía 1.17): la compresión deja quieta la imagen de cada bead y después la actualiza de golpe. $D$ cambia ± 8 % según cuántos retardos se ajustan. De este video no se pueden sacar conclusiones sobre gaussianidad ni sobre la correlación entre pasos.

### 9. BROWNIANO 2: parámetros de conversión

`convertir.py "...\BROWNIANO 2.tif" --nombre "BROWNIANO 2" --fps 10.021 --escala_nm 94.61`, todo por defecto. Los pasos intermedios (`figuras/BROWNIANO 2/pasos/`) muestran que la detección es sensata: se detectan las beads que se mueven y las pegadas quedan en el fondo. Tardó 17 minutos.

El Radio por defecto (17 px) es el de las beads nominales. Las partículas chicas (7 a 12 px) se detectan igual, pero con menos contraste: algunas trayectorias de chicas pueden estar más cortadas.

### 10. BROWNIANO 1: la pausa de escritura

Las partes 26 y 27 del TIFF (cuadros 4342 a 4675) tardaron 3.69 s y 1.10 s más que las demás en escribirse. Entre los cuadros 4414 y 4415 se corta una de las dos trayectorias que había (la 36) y empiezan dos nuevas (37 y 38); la otra (35) termina diez cuadros después. Es lo que pasa si ahí faltan varios segundos: las partículas aparecen corridas más que el salto máximo del enlace. La segunda pausa cae cerca del cuadro 4540. Se descartan los 334 cuadros de esas dos partes y al tiempo de los cuadros siguientes se le suman 4.79 s.

**Impacto:** se pierde 6.5 % de los cuadros. Las trayectorias 40 y 42 empiezan después y no se tocan. En el resultado de las beads solas (trayectorias 0, 1 y 2, todas anteriores), ninguno.

### 11. El fps de las grabaciones nuestras

Los TIFF no traen el tiempo de cada cuadro. Con las fechas de modificación de los archivos, cada parte de 167 cuadros da 10.0215 cuadros/s, con un desvío de 0.004 entre partes, en las dos grabaciones. BROWNIANO 1 se convirtió con 10.02 y BROWNIANO 2 con 10.021.

**Impacto:** 0.04 % en $D$. Lo que no se puede saber es si se perdió algún cuadro suelto dentro de una parte.

### 12. Qué trayectorias entran

- **Largo mínimo:** 200 cuadros medidos (20 s en las nuestras, 13 s en GROS) y 100 en CHICAS (4.5 s; el video dura 11.7 s). Los notebooks de cada grabación usan 300. Con 200 entran las mismas 11 trayectorias en BROWNIANO 1 y 22 en vez de 15 en BROWNIANO 2 (las que se suman son partículas chicas, que salen antes del campo).
- **Colisión:** se descartan los puntos marcados (17 % en BROWNIANO 1, 3 % en BROWNIANO 2 y GROS, 11 % en CHICAS).
- **Confinadas:** pendiente log-log del MSD menor que 0.5 (la 3 de BROWNIANO 1 y dos de CHICAS).

**Impacto:** el sesgo de las trayectorias cortas está medido con la simulación gemela (1 % con la recta; 2.6 % en CHICAS) y corregido.

### 13. Las clases de partícula se pusieron mirando las imágenes

`figuras/<medición>/radios/recortes.png` tiene seis recortes de cada partícula. Lo que se decidió a ojo está en `MEDICIONES` (primera celda del notebook):

- BROWNIANO 1: agregados 6 (dos beads pegadas), 40 y 42 (racimos); desenfocadas 24, 26, 27.
- BROWNIANO 2: agregados 3 y 34 (dos beads pegadas); desenfocada 15.
- GROS: la 2 tiene una vecina que la acompaña; la 15 es una mancha grande.
- CHICAS: nada a ojo.

Después, por radio aparente: menos de 1.6 µm es "chica" y más de 2.35 µm es "agregado".

**Impacto:** el $D$ de cada medición sale solo de las beads "solas". En BROWNIANO 1 quedan **3** (0, 1 y 2). Si se promediaran todas las libres daría 0.068 en vez de 0.056 (BROWNIANO 1) y 0.181 en vez de 0.079 (BROWNIANO 2).

### 14. Qué es el radio aparente

Es el radio donde el perfil radial (mediana en anillos, relativa al fondo local), yendo hacia afuera desde su mínimo, recupera la mitad de su profundidad. Se mide sobre el canal verde del video original, en 30 cuadros por partícula; la incerteza es el desvío entre cuadros / $\sqrt{n}$ (0.01 a 0.1 µm).

**Qué tan confiable es:** para comparar partículas de una misma grabación, bien (las beads solas de BROWNIANO 1 dan 1.79, 1.80 y 1.82 µm). Como radio absoluto, no: da 1.8 µm para beads de 1.5 µm nominal, porque el borde exterior incluye la difracción y el desenfoque. Con otra definición (el mínimo del perfil, el anillo oscuro, cuando se ve) sale 1.1 a 1.4 µm. El error sistemático es de unos 0.3 µm y depende del foco: una partícula que sube o baja cambia su radio aparente (la columna `disp_r_borde_px` de `radios.csv` dice cuánto).

### 15. El método de deriva que se informa

Recta por partícula en las cuatro, con el factor de la simulación gemela. La incerteza "método" es la mitad del rango entre los métodos que valen para cada medición (1 a 4 %).

- En las simulaciones gemelas el ruido de localización de cada partícula sale de la ordenada de su MSD; si la ordenada es negativa (exposición), se pone cero. La exposición no se simula.
- En la gemela "con deriva", el vaivén común se estima promediando las partículas y trae el browniano de ellas dividido por $N$. Con 2.3 partículas a la vez (BROWNIANO 1) exagera la deriva: el 8 % de más que da ahí la recta es una cota por arriba. El chequeo que vale es que la media móvil corregida no dependa de la ventana.

### 16. Ajuste de $D$

Primeros 4 retardos del MSD de cada partícula, sin pesos, como en los notebooks de cada grabación. Retardo máximo 30 cuadros.

**Impacto:** ajustando 2, 3, 6 u 8 retardos $D$ cambia ± 2 % (BROWNIANO 1), ± 4 % (BROWNIANO 2), ± 0.6 % (GROS) y ± 8 % (CHICAS). Está sumado en la incerteza.

### 17. Viscosidad de la mezcla

`fisica.viscosidad_glicerina_agua` (Cheng 2008) pide la fracción **en masa**. Los volúmenes del cuaderno se pasaron a masa con 1.261 g/ml (glicerina) y 0.9975 g/ml (agua a 23 °C): 0.25 ml + 1 ml → 24.0 %; 65 ml + 60 ml → 57.8 %. Se supone glicerina pura (si era glicerina de farmacia, puede traer agua).

### 18. Figuras: lo que la guía no cubre

- Con cuatro mediciones en un gráfico (`D_por_metodo`) se usaron los dos colores de la guía y el azul y el violeta de la paleta de seaborn, cada serie con su marcador.
- Las partículas que no son beads solas van con el color de su serie y el marcador hueco (el criterio de "datos excluidos del ajuste").
- Las trayectorias se dibujan sobre el primer cuadro en gris, con ejes en px; no llevan barras de error.
- El perfil radial va con marcadores sin barras de error (el CSV guarda el perfil promedio).
- El eje de la ventana es logarítmico, con los ticks en los valores usados.
- En `validacion_sinteticos` el eje y no empieza en cero: es un cociente que vale 1.

## Para el otro Claude (no son dudas de física)

- `tiempo.csv` y `trayectorias.csv` guardan `t_s` con 6 cifras significativas: pasados los 100 s la resolución es de 1 ms. A 10 cuadros/s no molesta (1 % del período, sin sesgo), pero en una grabación a 100 cuadros/s o más sería 10 % del período. Conviene escribir más cifras.
- `guardar_pasos.py` deja unos 100 MB por grabación en `figuras/<grabación>/pasos/` (los PNG a resolución completa). Los PNG no van a git, pero sí a OneDrive.
- El Python de TouchDesigner no trae pandas: `medir_radios.py` lee y escribe los CSV con `csv`.
