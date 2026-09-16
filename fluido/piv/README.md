# PIV del vórtice con glitter

Todo está en `TDS GLITTER/TDS GLITTER.toe`, en dos componentes con su propia interfaz
(click derecho sobre el componente → *Open Viewer*):

| componente | qué hace |
|---|---|
| `/project1/grabar` | graba **una o dos cámaras en color**, cuadro a cuadro (TIFF sin pérdida), con la hora de cada cuadro |
| `/project1/analizar` | PIV + seguimiento de la **partícula central**; guarda CSV por cámara en `data glitter` e imágenes de cada etapa para el póster |

Cámara **A** = vista superior (da x, y). Cámara **B** = vista lateral (da la altura z).
Adentro de cada componente, la red está ordenada de izquierda a derecha con placas de
comentario por etapa, y `agents_md` tiene la guía completa.

## 1. Grabar

1. Página **Cámara A** / **Cámara B**: `Fuente` (Apagada / Cámara / Simulación),
   `Dispositivo` y `Espejar horizontal` (si la cámara entrega la imagen espejada; si no,
   el giro sale invertido).
2. Página **Grabar**: `Nombre`, `Agitador (rpm)` y **GRABAR**. Otra vez GRABAR para parar.

Se graba en **color** (hace falta para elegir después el color del glitter y de la partícula).
Cada grabación es una carpeta dentro de `TDS GLITTER/GRABACIONES/` (~2.6 MB por cuadro a
1280×720):

    GRABACIONES/<fecha>_<nombre>/
      grabacion.json                qué grabó cada cámara, rpm
      camara_A/cuadros/*.tif        imágenes
      camara_A/tiempos.csv          indice, t_s  (+ xc_px, yc_px: centro real, solo simulación)
      camara_B/...

La **simulación** (página Simulación) dibuja glitter en un vórtice de Rankine, una partícula
en el centro y un centro que se desplaza (`Deriva del centro`): sirve para practicar y para
comprobar el análisis.

⚠️ No minimices TouchDesigner mientras graba o analiza: minimizado deja de procesar.

## 2. Analizar

1. **Grabación**: elegí esa carpeta. Se carga sola. **Cámara**: A o B.
2. **Marcar** (click sobre la imagen), para cada cámara:
   - `Punto A de la regla` → click; pasa solo a B → click. Cargá la `Distancia real` en la página de la cámara.
   - `Centro del vórtice` → click; pasa solo al `Borde del ROI` → click.
   - `Color del fondo` → click en el agua; pasa solo a `Color del glitter` → click sobre un
     destello; pasa solo a `Color de la partícula` → click sobre la partícula.
     El estado dice si cada imagen quedó separada del otro color.
3. Página **Partícula**: `Seguir partícula central`, `Umbral` y `Área mínima`. Con
   `Mostrar = 3 partícula` se ve qué pasa el umbral: la partícula entera en rojo y nada más.
4. **ANALIZAR TODO**. Cada análisis va a su propia carpeta (nunca pisa uno anterior):

       data glitter/<grabación>/<fecha y hora del análisis>/
           parametros.json             calibración, colores, seguimiento, dt, método
           camaraA_campo.csv           el campo de velocidades de cada par
           camaraA_centro.csv          el centro del vórtice en cada par
           camaraB_...                 (si hay cámara B)
           imagenes/camaraA_cuadros<k>_<k+s>/    imágenes de un par para el póster (ver abajo)

   Si el par que estabas mirando ya estaba calculado, sus imágenes se guardan en `imagenes/`.
   **Detener** corta y guarda lo hecho.

### Cómo se procesa cada par (y qué muestra cada vista)

| `Mostrar` | etapa |
|---|---|
| 1 cuadro original | el cuadro tal cual se grabó |
| 2 proyección del glitter | color → una intensidad: 1 en el glitter, 0 en el fondo **y en la partícula** |
| 3 partícula | proyección de la partícula (0 en fondo y glitter), umbral y centroide, con ampliación |
| 4 CLAHE | ecualización adaptativa: lo que entra al PIV |
| 5 ventanas de la iteración | ventanas de interrogación de la iteración elegida, de su tamaño real, donde se leen en el cuadro B, sobre los dos cuadros sin deformar |
| 6 A y B superpuestas | cuadros A (rojo) y B (azul) después de deformarlos con el predictor: negro = coinciden |
| 7 correlación | una ventana A, la B y su correlación cruzada FFT, con la zona de búsqueda y el pico |
| 8 vectores | campo validado de la iteración; color = rapidez, rojo = interpolado |

Las vistas 5 a 8 guardan **todas las iteraciones**: se elige cuál con `Iteración que se ve`.

**Separación de cuadros**. Por defecto 2: cada par compara el cuadro k con el k + 2 (se analizan
todos los k, así que sigue habiendo un campo por cuadro). Con cuadros seguidos el desplazamiento
es de unos 3–4 px y el error del PIV (~0.05 px) pesa más en relación. Medido en la simulación
(mediana de 3 pares, contra el desplazamiento exacto):

| separación | desplazamiento | error mediano | error relativo | error p95 | interpolados |
|---|---|---|---|---|---|
| 1 | 3.4 px | 0.048 px | 1.45 % | 0.20 px | ~2 |
| 2 | 6.9 px | 0.058 px | **0.90 %** | 0.38 px | ~3 |
| 3 | 11.3 px | 0.073 px | 0.71 % | 0.66 px | ~7 |
| 4 | 14.6 px | 0.089 px | 0.66 % | 0.95 px | ~11 |

Conviene que el desplazamiento no pase de ~1/4 de la ventana final (8 px con 32 px): más
separación casi no mejora y crecen los vectores malos (sobre todo cerca del núcleo, donde el
giro dentro de la ventana es mayor). Con cámaras de 30 fps, separación 2 = pares de 1/15 s.

**Iteraciones**. Por defecto 3, como PIVlab: ventanas de 128, 64 y 32 px, cada una con el campo
de la anterior como predictor. Se puede poner de 3 a 8: desde la cuarta se repite la ventana
final (como "repetir la última pasada" de PIVlab) y se refina el predictor. Medido en la simulación (un par, contra el
desplazamiento exacto): con 3 pasadas el error es 0.049 px (mediana) y 0.20 px (p95); con 5,
0.053 y 0.15 px. La mediana no cambia; mejora algo donde el campo varía rápido.

**Qué hace la deformación de ventana** (vista 6). Con el campo de la iteración anterior d(p),
el cuadro A se lee en p − d(p)/2 y el B en p + d(p)/2 (Lanczos 8×8): cada destello queda en la
mitad de su recorrido en las dos imágenes. Si el predictor es bueno, A y B deformadas coinciden
y la correlación solo mide lo que falta corregir. En la superposición se ve: en la iteración 1
(sin predictor) los destellos aparecen dobles, rojo y azul; desde la 2 quedan negros.

**Proyección de color**: valor = (color − fondo) · w, con w perpendicular a (otro color − fondo)
y escalado para que el color buscado dé 1. Es la dirección de máximo contraste que además
borra el otro color.

**Partícula**: manchas del umbral dentro del ROI, la más cercana al centro del par anterior,
centroide pesado por la intensidad (subpíxel). El centro de un par es el promedio de sus dos cuadros.

**PIV**, como PIVlab: CLAHE → FFT en 3 pasadas con deformación de ventana (Lanczos) → pico
gaussiano de 3 puntos → validación (mediana local 3, desvío estándar 8) e interpolación.

### Imágenes para el póster

Página **Poster**: `Zoom` (px), `Ventana de ejemplo` (o `Marcar = Ventana de ejemplo`, que
elige dónde se hace el zoom) y **GUARDAR IMÁGENES**. Guarda en `imagenes/camaraA_cuadros<k>_<k+s>/` del
**último análisis** terminado de esa grabación (si no hay ninguno, crea la carpeta), una
imagen de cada tipo por iteración (3 por defecto, o las que diga `Iteraciones`):

| archivo | |
|---|---|
| `ventanas_iterN.png` | **lo que considera la iteración N**. Zoom de 384 × 384 px con los cuadros A (rojo claro) y B (azul claro) sin deformar: cada destello aparece dos veces, separado por lo que se movió. Encima: una grilla negra con las ventanas de interrogación de esa iteración, de su tamaño real (128, 64, 32 px), en el lugar donde se leen en el cuadro B, p + d(p)/2, con d = el campo de la iteración anterior (en el cuadro A se leen en p − d(p)/2); en la 1 la grilla es recta; en las siguientes, debajo, la misma grilla sin mover en gris claro como referencia. En cada ventana, una flecha violeta con el **desplazamiento promedio que calculó la iteración** (ampliado, mismo factor en todas, valor en `leeme.txt`). |
| `correlacion_iterN_ventanaA.png`, `_ventanaB.png`, `_plano.png` | la ventana de ejemplo (`Marcar = Ventana de ejemplo`; su posición está en `leeme.txt`) tal como la compara la iteración N: leída en A y en B (ya corridas por el predictor), y su correlación cruzada. En el plano, el centro es desplazamiento cero, el cuadrado blanco la zona de búsqueda y el anillo el pico: lo que falta corregir. En la 1 el pico está lejos del centro; en las siguientes, cerca |
| `desplazamiento_iterN.png` | el desplazamiento medido en la iteración N en todo el ROI (flechas azules; rojas = interpoladas), con la misma escala de flechas en todas las iteraciones |

Con 50 % de solapamiento hay además una ventana centrada en cada esquina de las dibujadas:
se dibuja una sí y una no para que se lea el tamaño.

Son **solo imágenes, sin texto**, de 2048 px, con fondo blanco. `leeme.txt` tiene los números
para rotularlas: región y ampliación, escala en mm por px, escala de las flechas, y por
iteración ventana, paso, vectores, interpolados y desplazamiento mediano.

### Los CSV

`camaraA_campo.csv`

| columna | |
|---|---|
| `frame` | primer cuadro del par |
| `t_s` | tiempo del centro del par, desde el primer cuadro de esa cámara |
| `x_px`, `y_px` | centro de la ventana; origen abajo-izquierda, y hacia arriba |
| `x_mm`, `y_mm` | ídem en mm (en la cámara B, `y` es la altura z) |
| `u_mm_s`, `v_mm_s` | velocidad (en px/s si esa cámara no está calibrada) |
| `interpolado` | 1 = la validación lo descartó y se reemplazó por sus vecinos |

`camaraA_centro.csv`: `frame, t_s, xc_px, yc_px, xc_mm, yc_mm, seguido, area_px`
(`seguido` = 0: no se encontró la partícula y se repite la última posición).

Para graficar: `fluido/analisis/graficar_glitter.ipynb` (centro en el tiempo, campo medio,
perfil v_θ(r), evolución de v(r), calidad).

**Paso temporal**: con cámara, dt fijo por cuadro (duración / intervalos), porque la hora
viene redondeada a ~4 ms y la cámara tiene un ritmo fijo; un cuadro perdido se detecta y ese
par usa 2·dt. Con la simulación, la diferencia de horas de cada par (su reloj es exacto, pero
el ritmo de TD varía).

## Validación

Simulación en color, centro que deriva 40 px, partícula central:

- **TD contra `piv_core.py`** sobre los mismos cuadros: diferencia mediana 0.003 px.
- **Centro seguido contra el real**: mediana 0.05 px, máximo 0.14 px, partícula encontrada en todos los pares.
- **ω(r)** medida desde el centro seguido contra el vórtice inyectado: ver `regresion_td.py`.
- **Grabación con dos webcams reales** (640×480): ~29 fps efectivos por cámara.

    "C:/Program Files/Derivative/TouchDesigner.2025.33070/bin/python.exe" regresion_td.py <carpeta_grabación> <carpeta_del_análisis>/parametros.json [pares]

## Archivos de esta carpeta

`piv_core.py` motor de referencia (lo copia el DAT `piv_numerico` de TD) ·
`piv_sim.py` simulador con campo conocido · `regresion_td.py` prueba de TD contra la referencia.
