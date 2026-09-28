# PIV del vórtice con glitter

Todo está en `TDS GLITTER/TDS GLITTER.toe`, en tres componentes con su propia interfaz
(click derecho sobre el componente → *Open Viewer*):

| componente | qué hace |
|---|---|
| `/project1/grabar` | graba **una o dos cámaras en color**, cuadro a cuadro (TIFF sin pérdida), con la hora de cada cuadro |
| `/project1/analizar` | PIV + seguimiento de la **partícula central**; guarda CSV por cámara en `data glitter` e imágenes de cada etapa para el póster |
| `/project1/figuras` | **figuras** del análisis de imagen de un par, una a la vez y cada una con su configuración; independiente de ANALIZAR (ver *3. Figuras*) |

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
      camara_A/tiempos.csv          cámara:     indice, t_s, cuadro_camara, t_driver_s, t_llegada_s
                                    simulación: indice, t_s, xc_px, yc_px (centro real)
      camara_B/...

**La hora de cada cuadro de cámara** se arma con tres relojes que se guardan juntos:
la hora del driver (`t_driver_s`: viene redondeada a 1/64 s = 15.6 ms y con jitter), el
**contador de cuadros** de la cámara (`cuadro_camara`: exacto; un salto de 2 es un cuadro
perdido de verdad) y la hora en que llegó a TouchDesigner (`t_llegada_s`, de respaldo). Al
parar, `t_s` = una recta ajustada a la hora del driver contra el número de cuadro, en una
ventana de ±15 cuadros (así se borra el redondeo y el ritmo puede cambiar: con poca luz la
webcam baja de 30 a 20 fps). Medido con la webcam: intervalos entre 49.9 y 50.6 ms, contra
un error de 6.8 ms de la hora cruda. `grabacion.json` anota `fps_medido` y `cuadros_perdidos`.
El análisis usa `cuadro_camara` para saber cuántos intervalos hay entre los dos cuadros de
cada par; con grabaciones viejas (sin esa columna) lo estima de la hora, como antes.

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
4. Página **Recipiente**: `Volumen de agua (ml)`, `Altura del agua (mm)` y `Circunferencia en
   la línea de agua (mm)`. Van a `parametros.json` → `recipiente` (con `radio_linea_agua_mm` =
   circunferencia / 2π; 0 = sin dato → `null`). Al elegir una grabación se recuperan del último
   análisis de esa grabación; si no tiene, quedan los de antes y el estado avisa que los revises.
5. **ANALIZAR TODO**. Cada análisis va a su propia carpeta (nunca pisa uno anterior):

       data glitter/<grabación>/<fecha y hora del análisis>/
           parametros.json             calibración, colores, seguimiento, dt, método, recipiente
           camaraA_campo.csv           el campo de velocidades de cada par
           camaraA_centro.csv          el centro del vórtice en cada par
           camaraB_...                 (si hay cámara B)
           imagenes/camaraA_cuadros<k>_<k+s>/    todas las etapas de un par, para el póster
           imagenes/camaraA_trayectoria.png       recorrido de la partícula central

   Con `Guardar imágenes al analizar` prendido (página Poster) se guardan las imágenes de un
   par mientras se lo analiza (ver *Imágenes para el póster*).
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

**Partícula**: en cada cuadro se busca el pico de la proyección de la partícula **cerca del
centro del par anterior** (un disco de 0.35 radios del ROI; al mirar un par suelto, cerca del
centro marcado) y se umbraliza a una fracción de ese pico (`Umbral`, 0.5 por defecto); se toma
la mancha que contiene el pico, con el centroide pesado por la intensidad (subpíxel). El centro
de un par es el promedio de sus dos cuadros. La partícula no está a la vista (y ese par queda
sin seguir) si el pico es menor que 0.012 o si no es al menos 3 veces más brillante que su
entorno con la mancha tapada.

Antes el pico se comparaba con el percentil 99.9 de todo el ROI, y en los tarros (filmados a
1920 px) la partícula ocupaba más que ese 0.1 % del ROI: se comparaba consigo misma y se
descartaba siempre. Con la búsqueda local los tarros pasaron de 0 % a 100 % de pares seguidos;
de paso, un reflejo en el borde del recipiente ya no le puede ganar a la partícula. La
dirección de color de la partícula, además de borrar el glitter, **ignora el blanco** cuando se
puede (los reflejos de la luz en la superficie son blancos). Lo que sigue sin funcionar es real:
en GLICERINA_barrido_5 la partícula no se ve (hundida en el hoyo del vórtice), y en
TARRO_ALTO_AGUA_2 y 3 la partícula y el glitter quedaron marcados casi del mismo color: hay
que volver a marcarlos (click justo en el centro de la partícula).

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

Página **Poster**:

| parámetro | qué hace |
|---|---|
| `Guardar imágenes al analizar` | ANALIZAR TODO guarda las imágenes de un par mientras lo analiza (con el centro que da el seguimiento) |
| `Par de las imágenes` | cuál; −1 = el de `Par a mirar` |
| `Región de las imágenes` | centrada en el vórtice, centrada en un punto marcado (`Marcar = Centro de las imágenes`) o todo el ROI |
| `Lado de la región (px)` | el lado del recorte, en px de la imagen (512 por defecto) |
| `Ventana de la correlación (px)` | la ventana de ejemplo de las imágenes 09 (`Marcar = Ventana de ejemplo`) |
| **GUARDAR IMÁGENES AHORA** | las guarda ya, en el **último análisis** de la grabación; si el par no es el que se ve, lo calcula primero |

La región se ve como un **recuadro ámbar** sobre la imagen (en `original`, `glitter` y `clahe`).
Todas son **solo imágenes**, sin texto, con el mismo estilo: casi negro, el glitter en tonos
fríos y un solo color fuerte para lo que cada una muestra. Las de detalle son el recorte de la
región ampliado a 2048 px; `leeme.txt` dice qué es cada una y trae los números para rotularlas
(par, intervalo, escala en px por cm, rangos de color, ventanas).

    imagenes/camaraA_cuadros<k>_<j>/
      01_cuadro_completo.png         el cuadro entero: ROI (gris) y la región (recuadro naranja)
      02_cuadro.png                  detalle: el cuadro original
      03_proyeccion_glitter.png      detalle: la proyección de color del glitter
      04_imagen_piv.png              detalle: lo que entra al PIV (CLAHE o lineal)
      05_particula.png               primer plano de la partícula: mancha, centroide en A y en B
      06_superposicion_sin_deformar  A en naranja + B en celeste: cada destello aparece dos veces
      07_grilla_iter1..3.png         ventanas de cada iteración donde se leen + desplazamiento medido
      08_superposicion_deformada     A y B deformados con el campo: los destellos se juntan (blanco)
      09_correlacion_iterN_*.png     ventana A, ventana B y plano de correlación, por iteración
      10_vectores.png                todo el ROI: el campo final, color = rapidez
      11_vorticidad.png              todo el ROI: vorticidad (naranja antihoraria, celeste horaria)
      escala_rapidez / escala_vorticidad.png   las escalas de color, sin números (van en leeme.txt)
    imagenes/camaraA_trayectoria.png el recorrido de la partícula en toda la grabación, color = tiempo

En las grillas, cada ventana tiene su tamaño real (128, 64 y 32 px) y está corrida la mitad
del desplazamiento que midió la iteración anterior, así que la grilla se curva siguiendo el
giro; se dibuja una sí y una no (con 50 % de solapamiento hay otra centrada en cada esquina).
Las flechas van ampliadas con el mismo factor en las tres. Donde la validación descartó el
vector no hay flecha, y la vorticidad queda sin color donde falta algún vecino.

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

## 3. Figuras

`/project1/figuras` hace figuras del análisis de imagen de **un par de cuadros**, una a la vez,
**independiente de ANALIZAR**: tiene sus propios datos, sus propios parámetros de análisis y su
propio motor (numpy, a la **resolución nativa** de la grabación: sin el límite de 1280 px de TD).

| página | qué se elige |
|---|---|
| **Datos** | grabación, cámara, cuadro A y cuadro B (B = A + 2 salvo que escribas otro) |
| **Calibración** | escala (mm/px), centro y radio del ROI; `Marcar` + click sobre la figura (ROI, colores, centro del encuadre, ventana de ejemplo) |
| **Proyección** | colores de fondo, glitter y partícula, ignorar el blanco, CLAHE o lineal, zonas de CLAHE, umbral y área de la partícula |
| **PIV** | ventana final, iteraciones, pico contra segundo pico |
| **Figura** | cuál (cuadro, proyección, imagen del PIV, partícula, superposición A+B, grilla, correlación, vectores, vorticidad) y sus opciones: encuadre, cuadro A o B, iteración, flechas, densidad… |
| **Aspecto** | tema oscuro o claro, colores, paleta, grosor de líneas y grosor de flechas (por separado), brillo, afuera del ROI, escalas |

En la **grilla de una iteración** se elige además el **fondo** (cuadro A, B o los dos
superpuestos con los colores A y B), **qué se desplaza** (la grilla: imágenes sin deformar y
ventanas corridas con el predictor, en A a −d/2 y en B a +d/2; o la imagen: grilla recta sobre
las imágenes deformadas como las correlaciona esa iteración) y si las **flechas** son el
desplazamiento total o solo el residuo que midió la correlación después de deformar. La
grilla se muestra **redonda**: solo el ROI, con la grilla y las flechas recortadas al círculo.
| **Salida** | tamaño en px, nombre, carpeta (vacía = `fluido/figuras/…`); exportar esta figura o todas |

Calibración, proyección y PIV arrancan, al elegir la grabación, con los valores de su último
análisis (reescalados a la resolución nativa); después se cambian sin tocar ANALIZAR.
**Figura, Aspecto y Salida son de cada figura**: al elegir otra figura se carga su
configuración (se guarda con el proyecto). Solo se recalcula lo necesario: cambiar el aspecto
redibuja; cambiar un color de la proyección rehace la proyección y el PIV del par (~1.5 s a
1920 px); las figuras que no usan el PIV no lo calculan.

Exporta a la carpeta **`fluido/figuras/`**, todas juntas: `figuras/<grabación>/cam<L>_cuadros<a>_<b>/<nombre>.png`
(una subcarpeta por par, porque las figuras se llaman igual en todos), cada una con un
`<nombre>.txt`: qué es, los números para rotularla (px por cm, rangos de color, ampliación de
las flechas, iteración, ventana) y los parámetros con que se calculó el par.

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
