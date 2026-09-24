# PIV del vórtice con glitter

Todo está en `TDS GLITTER/TDS GLITTER.toe`, en dos componentes con su propia interfaz
(click derecho sobre el componente → *Open Viewer*):

| componente | qué hace |
|---|---|
| `/project1/grabar` | graba **una o dos cámaras en color**, cuadro a cuadro (TIFF sin pérdida), con la hora de cada cuadro |
| `/project1/analizar` | PIV + seguimiento de la **partícula central**; guarda CSV por cámara en `data glitter` e imágenes de cada etapa para el póster |

Para entender qué hace el análisis con cada imagen (proyección de color, CLAHE, PIV con
deformación de ventana y detección de la partícula), explicado simple y con figuras:
[como_funciona_el_analisis.md](como_funciona_el_analisis.md).

Cámara **A** = vista superior (da x, y). Cámara **B** = vista lateral (da la altura z).
Adentro de cada componente, la red está ordenada de izquierda a derecha con placas de
comentario por etapa, y `agents_md` tiene la guía completa.

## 1. Grabar

1. Página **Cámara A** / **Cámara B**: `Fuente` (Apagada / Cámara / Simulación),
   `Dispositivo` y `Espejar horizontal` (si la cámara entrega la imagen espejada; si no,
   el giro sale invertido).
2. Página **Grabar**: `Nombre`, `Agitador (rpm)`, `Duración` y **GRABAR**. Otra vez GRABAR
   para parar antes.

`Duración (s, 0 = a mano)`: con 0 graba hasta que se aprieta GRABAR de nuevo; con un número
corta sola al llegar (el estado va diciendo cuánto falta). Se mide con la **hora de los
cuadros** (la de la cámara), no con el reloj de TouchDesigner, así que el corte cae donde dice
el archivo: pedir 3 s dio 3.008 s, el primer cuadro que pasa el límite. La duración pedida
queda anotada en `grabacion.json`.

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

   ⚠️ Las marcas son **píxeles de la imagen**, así que dependen de con qué resolución se
   marcó. TouchDesigner sin licencia comercial limita los TOP a **1280 px**: una grabación de
   1920×1080 se analiza a 1280×720. Cada cámara guarda el ancho con el que se marcó
   (`Ancho de la imagen al marcar`) y al cargar la grabación las marcas se reescalan solas si
   hace falta, avisando en el estado. Sin eso el ROI queda corrido y la escala en mm/px 1.5
   veces mal, sin ningún error a la vista.

   - `Punto A de la regla` → click; pasa solo a B → click. Cargá la `Distancia real` en la página de la cámara.
   - `Centro del vórtice` → click; pasa solo al `Borde del ROI` → click.
   - `Color del fondo` → click en el agua; pasa solo a `Color del glitter` → click sobre un
     destello; pasa solo a `Color de la partícula` → click sobre la partícula.
     El estado dice si cada imagen quedó separada del otro color.
3. Página **Partícula**: `Seguir partícula central`, `Umbral (fracción del pico)` y
   `Área mínima`. El umbral es **relativo**: 0.5 = la mitad del brillo de la partícula en ese
   cuadro, así que no depende de haber marcado el color exacto ni de la luz. Con
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
| 4 imagen que entra al PIV | lo que se correlaciona: CLAHE o la proyección lineal (ver abajo) |
| 5 ventanas de la iteración | ventanas de interrogación de la iteración elegida, de su tamaño real, donde se leen en el cuadro B, sobre los dos cuadros sin deformar |
| 6 A y B superpuestas | cuadros A (rojo) y B (azul) después de deformarlos con el predictor: negro = coinciden |
| 7 correlación | una ventana A, la B y su correlación cruzada FFT, con la zona de búsqueda y el pico |
| 8 vectores | campo validado de la iteración; color = rapidez. Lo descartado **no se dibuja** |

Las vistas 5 a 8 guardan **todas las iteraciones**: se elige cuál con `Iteración que se ve`.

**Separación de cuadros**. Por defecto 2: cada par compara el cuadro k con el k + 2 (se analizan
todos los k, así que sigue habiendo un campo por cuadro). Con cuadros seguidos el desplazamiento
es de unos 3–4 px y el error del PIV (~0.05 px) pesa más en relación. Medido en la simulación
(mediana de 3 pares, contra el desplazamiento exacto):

| separación | desplazamiento | error mediano | error relativo | error p95 | descartados |
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

**Partícula**: en cada cuadro se busca el pico de la proyección de la partícula adentro del ROI y
se umbraliza a una fracción de ese pico (`Umbral`, 0.5 por defecto); de las manchas que quedan,
la más cercana al centro del par anterior, con el centroide pesado por la intensidad (subpíxel).
El centro de un par es el promedio de sus dos cuadros. Si el pico no se destaca del resto del ROI
(al menos 4 veces el percentil 99.9), la partícula no está a la vista y ese par queda sin seguir.

Antes el umbral era absoluto y **no se encontraba la partícula en ninguna medición**: el color
marcado era amarillo puro (0.91, 1.00, 0.00) pero la partícula adentro del agua se ve amarillo
pálido (0.61, 0.71, 0.47), así que su proyección llegaba a 0.27 y el umbral estaba en 0.49.
Con el umbral relativo, en la grabación DINAMICA queda una sola mancha de 90 a 200 px en todos
los cuadros, entre 10 y miles de veces más brillante que el fondo.

**PIV**, como PIVlab: CLAHE → FFT en 3 pasadas con deformación de ventana (Lanczos) → pico
gaussiano de 3 puntos → validación (mediana local 3, desvío estándar 8). Lo que no pasa la
validación se descarta.

### Cuánto tarda

Con los cuadros de la palangana (1280×720, ROI de 293 px, 3 iteraciones, separación 2):
**0.48 s por par**, unos 58 pares en 28 s. Antes eran ~2.2 s por par. De dónde salió:

| cambio | antes | ahora |
|---|---|---|
| dibujar las vistas del proceso en cada iteración (ANALIZAR TODO no las mira) | 0.61 s por iteración | no se dibujan; la pantalla muestra el cuadro original y así se ve avanzar |
| correlación: solo las ventanas que tocan el ROI (con 2 ventanas de margen) | 3476 ventanas de 32 px | 1563 |
| correlación: las FFT repartidas en 8 hilos (numpy suelta el GIL mientras calcula) | 402 ms por par | 174 ms |

Los vectores de adentro del ROI dan **idénticos** recortando o no (diferencia máxima medida
0.001 px): afuera del ROI no hay flujo que medir. Al terminar el análisis, el par que se esté
mirando se recalcula con todas sus vistas.

`numba` no entra acá: no está instalado en el Python de TouchDesigner y, sobre todo, no
acelera `numpy.fft`, que es el 80 % de la cuenta. Repartir las FFT en hilos da lo mismo que
buscaba el `prange` y no agrega dependencias. Si alguna vez hace falta más: correlacionar en
la GPU (GLSL) o bajar la resolución de análisis.

### Imágenes para el póster

Página **Poster**: `Zoom` (px de la imagen, por defecto 512) y **GUARDAR IMÁGENES**. Guarda
las imágenes del par que se está mirando (`Par a mirar`) en
`imagenes/camaraA_cuadros<k>_<j>/` del **último análisis** terminado de esa grabación (si no
hay ninguno, crea la carpeta). ANALIZAR TODO también las guarda, del par que se estaba mirando.

Una imagen por iteración, `grilla_iter1.png`, `grilla_iter2.png`, `grilla_iter3.png`: un
recorte cuadrado **centrado en el vórtice** (la partícula del centro si se la encontró), con el
cuadro de fondo tal como lo ve el PIV y encima las **ventanas de interrogación de esa
iteración**, de su tamaño real (128, 64 y 32 px). En la primera la grilla es recta; en las
siguientes cada ventana está corrida la mitad del desplazamiento que midió la iteración
anterior, así que la grilla se curva siguiendo el giro. El borde del recorte cae en múltiplos
de la ventana más grande, así que las ventanas entran enteras en todas.

Encima de la grilla, una **flecha blanca por ventana** con el desplazamiento que midió esa
iteración, centrada en la ventana y ampliada con el mismo factor en las tres imágenes (para que
se puedan comparar; el factor está en `leeme.txt`). Donde la validación descartó el vector no
hay flecha.

Son **solo imágenes**, de 2048 × 2048 px: sin texto ni grillas de referencia. El glitter va en
tonos fríos sobre casi negro, la grilla en un único color cálido y las flechas en blanco con un
contorno oscuro fino, para que se lean sobre zonas claras y oscuras. Con 50 % de
solapamiento hay otra ventana centrada en cada esquina de las dibujadas: se dibuja una sí y
una no para que se lea el tamaño. `leeme.txt` tiene lo necesario para rotularlas: el par, el
intervalo entre los cuadros, la región recortada, la escala en px por cm y el tamaño de las
ventanas de cada iteración.

**Lo que no se midió no se inventa.** La validación descarta vectores (ventana sin textura,
pico fuera del rango de búsqueda, test de mediana local, desvío estándar). Antes esos huecos se
rellenaban con el promedio de los vecinos y salían al CSV marcados `interpolado = 1`; ahora
salen como **NaN** con `valido = 0`, y tampoco se dibujan. El relleno sigue existiendo *adentro*
del multipasada, porque la iteración siguiente necesita un predictor en toda la grilla para
deformar, pero eso no sale al resultado. `Rellenar huecos` (página Análisis) vuelve al
comportamiento viejo: escribe el promedio de los vecinos y los marca igual con `valido = 0`.

En los datos de la palangana la validación descarta un **13 %** de los vectores, casi todos por
el test de mediana local y concentrados cerca del núcleo, donde el giro adentro de una ventana
de 32 px es demasiado grande. Eso no es ruido del método: es que ahí, con esa ventana y esa
separación de cuadros, no hay medición. Cuánto costaba inventarlos: en la simulación, donde se
conoce el campo, los vectores rellenados con los vecinos erraban **~1.4 px** mientras los
medidos erraban 0.15 px.

Cuánto se descarta depende de la siembra: run1 (poco glitter) **29 %**, run2 15 %, run3 13 %,
la simulación 0.5 %. Es una medida directa de cuán bien sembrada quedó la superficie.

`Relación de picos mínima` (página Análisis, por defecto 1.00 = apagada) descarta además los
vectores cuyo pico de correlación no le gana al segundo pico por ese factor. Con glitter denso
el pico es ancho y el criterio resulta romo: medido sobre 10 pares reales, comparando
TouchDesigner contra el motor de referencia (las mismas cuentas sobre imágenes que difieren en
1 LSB), el 98.8 % de los vectores coincide a mejor que 0.5 px, y exigir 1.10 corta la cola
(máxima diferencia de 15.8 a 3.6 px) pero a costa de perder el 35 % de los vectores. Queda
disponible, apagada por defecto.

**Imagen para correlacionar** (página Análisis):

| opción | qué hace |
|---|---|
| CLAHE (como PIVlab) | escala al percentil 99.9, cuantiza a 8 bits y ecualiza el contraste por zonas de ~64 px |
| lineal (sin cuantizar) | la proyección de color tal cual, en 32 bits: sin recortar los destellos más brillantes y sin perder los decimales |

En el camino del glitter **no hay ningún umbral**: el único umbral del proyecto es el de la
partícula central (página Partícula). Lo más parecido a binarizar que había era el recorte al
percentil 99.9 y la cuantización a 8 bits que pide CLAHE; `lineal` los saca. Medido: en la
simulación (campo conocido) CLAHE sigue siendo mejor — error mediano 0.148 px contra 0.160 px
de lineal — y sobre los cuadros reales descarta menos vectores (12.8 % contra 14.1 %). Ecualizar
por zonas rescata el glitter de las partes oscuras, y eso pesa más que los decimales perdidos.

### Los CSV

`camaraA_campo.csv`

| columna | |
|---|---|
| `frame` | primer cuadro del par |
| `t_s` | tiempo del centro del par, desde el primer cuadro de esa cámara |
| `x_px`, `y_px` | centro de la ventana; origen abajo-izquierda, y hacia arriba |
| `x_mm`, `y_mm` | ídem en mm (en la cámara B, `y` es la altura z) |
| `u_mm_s`, `v_mm_s` | velocidad (en px/s si esa cámara no está calibrada) |
| `valido` | 1 = medido; 0 = la validación lo descartó (`u`, `v` vienen NaN) |

`camaraA_centro.csv`: `frame, t_s, xc_px, yc_px, xc_mm, yc_mm, seguido, area_px`
(`seguido` = 0: no se encontró la partícula y se repite la última posición).

Dos notebooks leen esto:

- `fluido/analisis/graficar_glitter.ipynb`: mira **un** análisis (centro en el tiempo, campo
  medio, perfil v_θ(r), evolución de v(r), calidad).
- `fluido/analisis/modelos_vortice.ipynb`: compara **todas** las mediciones contra los modelos
  de vórtice (Rankine, Burgers, Lamb–Oseen, confinado, Vatistas), y saca circulación, radio del
  núcleo, vorticidad, velocidad radial, viscosidad efectiva y cómo escalan con las rpm, el
  fluido y el recipiente. Usa que los dos puntos de la calibración se marcaron diametralmente
  opuestos, así que de ahí salen el radio y el centro del recipiente.

**Paso temporal**: con cámara, dt fijo por cuadro, porque la cámara tiene un ritmo fijo y la
hora de cada cuadro no es confiable salto a salto; un cuadro perdido se detecta y ese par usa
2·dt. Con la simulación, la diferencia de horas de cada par (su reloj es exacto, pero el ritmo
de TD varía).

Cómo se saca ese dt, y por qué no es obvio. La hora que da el driver de la webcam viene
**cuantizada** (a 1/128 s = 7.8 ms) y además con **jitter**. Medido sobre los 3785 intervalos
de todas las grabaciones, tomando como unidad el período real T:

| salto entre cuadros | qué es |
|---|---|
| 0.89 – 1.65 T (99.3 % de los casos, promedio 1.000 T, desvío 0.16 T) | un solo cuadro: cuantización y jitter |
| nada entre 1.65 y 1.80 T | — |
| 1.80 – 2.12 T | ahí sí falta un cuadro |

Así que el período se estima **promediando los saltos de un solo cuadro** (promediar muchos
borra la cuantización) y un salto cuenta como cuadro perdido recién a partir de **1.7 T**, que
cae justo en el valle. Los dos números están medidos, no elegidos: comparando con el PIV el
desplazamiento de cada par contra el de sus vecinos, los saltos de 1.55 – 1.58 T dan el mismo
desplazamiento (no falta nada) y los de 1.80 – 2.26 T dan el doble (falta un cuadro).

Tomar la mediana de los saltos como período —lo que se hacía antes— elige el valor más
repetido, que con la cuantización es 31.25 ms y no 33.4 ms. Con ese período, cada salto de
1.5 T parecía un cuadro perdido: ese par se dividía por un dt 1.5 veces más grande y su
velocidad salía **2/3** de la real (rachas regulares, bien visibles en la evolución temporal),
y de paso **todos** los demás pares quedaban 7 % altos. Las tres tomas de la palangana pasan
de "32 fps con 4 a 7 cuadros perdidos" a **30.0 fps sin ningún cuadro perdido**.

El método aguanta hasta ~30 % de cuadros realmente perdidos. Si una grabación da un fps muy
distinto del que declara la cámara, conviene mirarla.

**Pares que cruzan un cuadro perdido de verdad.** El segundo cuadro de cada par es el que está a
`Separación` intervalos **de cámara** del primero, no a `Separación` posiciones en la lista de
cuadros guardados. Si justo ese cuadro se perdió, se usa el anterior que exista (nunca uno más
lejos). Antes, en DINAMICA (5 cuadros perdidos de verdad, tres seguidos), los pares que los
cruzaban abarcaban 3 o 4 intervalos: el desplazamiento pasaba de 12 a 21 px, el PIV se quedaba
sin rango de búsqueda, descartaba los vectores rápidos y la velocidad del par caía un 15 %.

## Validación

**Simulación en color** (vórtice de Rankine, centro que deriva 40 px, partícula central,
340 pares a 44 fps, separación 2), con el motor actual:

- **TD contra `piv_core.py`** sobre los mismos cuadros: diferencia mediana **0.003 px**, p95 0.008 px.
- **Centro seguido contra el real**: 340 de 340 pares, error mediano **0.054 px**, máximo 0.13 px.
- **ω(r) contra el vórtice inyectado** (desde el centro seguido, restando su velocidad):
  +0.87 % (r 20–60 px), +0.30 % (60–110), −0.19 % (110–160), +0.13 % (160–260).
- Descartados por la validación: **0.5 %**.

Para repetirlo: grabar la simulación (página Simulación, `Duración` ~8 s), poner
`Borde del ROI` en 0 para que el análisis tome el centro y el radio de la simulación,
ANALIZAR TODO y después

    "C:/Program Files/Derivative/TouchDesigner.2025.33070/bin/python.exe" regresion_td.py <carpeta_grabación> <carpeta_del_análisis>/parametros.json [pares]

**Cuadros reales** (palangana, glitter denso, 1280×720): ahí no hay verdad conocida, pero la
misma prueba contra `piv_core` da mediana 0.007 px y p95 0.10 px, con un 1.2 % de vectores que
difieren más de 0.5 px. Son ventanas donde el pico de correlación es ancho y una diferencia de
1 LSB en la imagen alcanza para que se elija otro pico: el glitter real forma filamentos, no
destellos separados como en la simulación. La validación descarta un 13 % de los vectores
(contra 0.5 % en la simulación), casi todo cerca del núcleo.

**Grabación con dos webcams reales** (640×480): ~29 fps efectivos por cámara.

## Archivos de esta carpeta

`piv_core.py` motor de referencia (lo copia el DAT `piv_numerico` de TD) ·
`piv_sim.py` simulador con campo conocido · `regresion_td.py` prueba de TD contra la referencia.
