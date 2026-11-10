// Informe del análisis de decaimientos de muones.
// Para compilar (desde la carpeta muones/, hace falta uv):
//   uv run --no-project --with typst python -c "import typst; typst.compile('informe/informe_decaimientos.typ', output='informe/informe_decaimientos.pdf', root='.')"
// Las figuras salen de analisis/reporte_filtros.py, analisis/dos_casos.py, analisis/galeria_rapidos.py,
// analisis/decaimientos.ipynb, analisis/energia_muones.ipynb e informe/fig_ejemplo_decaimiento.py.

#set document(title: "Decaimiento de muones: selección de eventos y vida media")
#set page(paper: "a4", margin: (x: 2.3cm, top: 2.3cm, bottom: 2.5cm), numbering: "1")
#set text(lang: "es", size: 10.5pt)
#set par(justify: true, leading: 0.62em)
#set heading(numbering: "1.1")
#show heading.where(level: 1): it => block(above: 1.5em, below: 0.8em, it)
#show heading.where(level: 2): it => block(above: 1.2em, below: 0.7em, it)
#show figure.caption: set text(size: 9pt)
#show figure: set block(above: 1.2em, below: 1.2em)
#set figure(gap: 0.7em)

#let fig(archivo, ancho: 72%) = image("../salida_decaimientos/" + archivo + ".pdf", width: ancho)
#let par2(a, b) = grid(columns: 2, gutter: 6pt, fig(a, ancho: 100%), fig(b, ancho: 100%))
#let tabla(columnas, alineacion, encabezado, ..celdas) = table(
  columns: columnas, align: alineacion, stroke: none, inset: (x: 7pt, y: 4.5pt),
  table.hline(), table.header(..encabezado.map(c => strong(c))), table.hline(stroke: 0.5pt),
  ..celdas, table.hline(),
)

#align(center)[
  #text(size: 17pt, weight: "bold")[Decaimiento de muones: selección de eventos y vida media]
  #v(0.3em)
  #text(size: 11pt)[Laboratorio 5 · análisis de las corridas del 11/9 al 29/10 de 2026]
  #v(0.2em)
  #text(size: 10pt, fill: luma(90))[Noviembre de 2026]
]
#v(0.8em)

#block(fill: luma(244), inset: 11pt, radius: 3pt, width: 100%)[
  *Resumen.* De 978.340~eventos registrados, 2421 son decaimientos: un pulso del muón que
  entra y, entre 100 y 1400~ns después, un pulso del electrón. Su distribución de tiempos da
  una vida media $tau = 2070^(+264)_(-211)$~ns, compatible con los 2197~ns del muón libre, y el
  valor no depende de los filtros ni de sus umbrales. Por debajo de 200~ns hay un exceso de
  unos 47~eventos (2,4~$sigma$) que no son pospulsos del fotomultiplicador. Separando los
  muones que frenan en el hierro de los que lo cruzan y frenan en el centellador~3, una
  segunda exponencial mejora el ajuste solo en el primer grupo (p-valor 0,07) y los datos son
  compatibles con la componente de unos 206~ns que se espera de los $mu^-$ capturados en el
  hierro, sin que alcancen para medirla. El límite principal es la ventana de medición
  (1,3~µs), más corta que una vida media. El informe incluye también lo que se pudo sacar
  sobre la energía y la velocidad de los muones, una lista de cosas para acordarse y
  propuestas a futuro.
]

= Montaje y datos

Tres barras centelladoras de 4,5 × 4,5 × 46~cm, cada una con su fotomultiplicador, apiladas
una sobre otra con una placa de hierro entre la segunda y la tercera (@fig-montaje). Las
tres señales se digitalizan con una placa DRS4 a 0,6~GS/s: cada evento son 1024 muestras
por canal, una cada 1,66~ns, o sea una ventana de 1,7~µs. La placa dispara cuando el
centellador~1 baja de −100~mV, y se guardaron los eventos que además tienen pulso en el
centellador~2. En las figuras, `u1`, `u2` y `u3` son los centelladores 1, 2 y 3.

#figure(
  image("../figuras/disp_apilados.svg", width: 50%),
  caption: [Disposición de los módulos para medir decaimientos. Los números indican cada
    centellador; la placa de hierro está entre el 2 y el 3. Esquema fuera de escala.],
) <fig-montaje>

Son nueve corridas, del 11 de septiembre al 29 de octubre de 2026: 10.216 archivos con
978.340~eventos (31~GB de texto) en unas 816~horas de medición. Las trazas están publicadas
en el _release_ `muones-decaimientos-v1` del repositorio; el archivo `muones/data/README.md`
explica cómo bajarlas.

= Cómo se reconoce un decaimiento

Casi todos los eventos son muones que atraviesan la pila: un solo pulso por centellador,
siempre cerca de los 207~ns de la traza (el instante del disparo). Un muón que frena y
decae deja ese mismo pulso de entrada y, más tarde, un segundo pulso: el del electrón
(@fig-ejemplo). El tiempo entre los dos, $Delta t$, es lo que se usa para medir la vida
media.

#figure(
  fig("reporte_ejemplo_decaimiento", ancho: 80%),
  caption: [Trazas de los tres centelladores en un decaimiento. La curva clara es la señal
    cruda y la oscura, la suavizada; los triángulos marcan los pulsos encontrados, con su
    amplitud. El muón deja pulso en los tres centelladores (línea punteada) y 943~ns después
    aparece el electrón en el centellador~3 (línea llena).],
) <fig-ejemplo>

Para buscar los pulsos, cada traza se suaviza con un promedio móvil de 10~ns, que elimina
el ruido periódico de unos 100~MHz de la placa, y se le resta la línea de base. Un pulso es
un tramo donde la señal baja más de 20~mV. De cada pulso se guarda cuándo empieza, su
amplitud (medida sobre la señal suavizada) y su ancho. El pulso tardío de un evento es el
primero que aparece más de 40~ns después de la entrada, en cualquiera de los tres
centelladores.

= Filtros

Los filtros se aplican uno sobre otro, en el orden de la @tab-filtros. Tres de ellos hacen
todo el trabajo: pedir que haya un segundo pulso, que tenga al menos 60~mV y que $Delta t$
esté entre 100 y 1400~ns. Los de ancho, pulsos extra y pulsos previos no sacan ningún
evento una vez aplicados los anteriores (@fig-filtros).

#figure(
  tabla(4, (left, left, right, right),
    ([Filtro], [Qué pide], [Quedan], [Saca]),
    [Todos los eventos], [], [978.340], [],
    [Traza limpia], [sin ráfagas de ruido], [977.131], [0,1~%],
    [Entrada en el centellador 1], [pulso de entrada de al menos 60~mV], [976.608], [0,1~%],
    [Hay pulso tardío], [un segundo pulso de al menos 30~mV], [4.664], [99,5~%],
    [Amplitud del pulso tardío], [al menos 60~mV], [2.690], [42,3~%],
    [$Delta t$ en rango], [entre 100 y 1400~ns], [2.421], [10,0~%],
    [Ancho del pulso tardío], [entre 10 y 30~ns], [2.421], [0],
    [Sin pulsos extra], [ningún tercer pulso], [2.421], [0],
    [Sin pulsos previos], [nada antes del disparo], [2.421], [0],
  ),
  caption: [Eventos que quedan después de cada filtro.],
) <tab-filtros>

El filtro de entrada pide 60~mV aunque el disparo era de 100~mV porque miden cosas
distintas: el disparo actúa sobre la señal cruda y el filtro sobre la suavizada, donde el
pico queda 1,36 veces más bajo. Un pulso que apenas dispara mide unos 74~mV suavizado y
pasa. Los 736~eventos que saca son picos de ruido angostos (12~ns contra 18~ns de un pulso
normal) o pulsos fuera del instante de disparo.

#figure(
  fig("reporte_filtros_eventos", ancho: 66%),
  caption: [Cantidad de eventos que quedan al aplicar los filtros en orden. Los puntos
    huecos son dos filtros opcionales de topología que no forman parte de la selección.],
) <fig-filtros>

= Vida media

El tiempo entre pulsos se ajusta por máxima verosimilitud, sin agrupar los eventos en
bins, con una exponencial normalizada en la ventana $[a, b]$ de 100 a 1400~ns:

$ p(Delta t) = e^(-Delta t slash tau) / (tau thin (e^(-a slash tau) - e^(-b slash tau))) $

El error de $tau$ es el intervalo donde el logaritmo de la verosimilitud baja 1/2 respecto
de su máximo. El resultado, con los 2421~eventos de la selección (@fig-ajuste), es

$ tau = 2070^(+264)_(-211) "ns". $

El valor de referencia no es exactamente 2197~ns. Los $mu^+$ decaen con la vida media
libre, pero los $mu^-$ que frenan en el plástico pueden ser capturados por un núcleo de
carbono, y su vida media efectiva ahí es de unos 2026~ns. Con la proporción de cargas de los
rayos cósmicos, una sola exponencial ajustada debería dar cerca de 2100~ns.

#figure(
  fig("ajuste_vida_media", ancho: 68%),
  caption: [Distribución del tiempo entre el pulso de entrada y el pulso tardío para los
    2421~eventos seleccionados, en bins de 50~ns con barras $sqrt(N)$. La curva es el ajuste
    de una exponencial, hecho sobre los eventos sin agrupar: $tau = 2070^(+264)_(-211)$~ns.],
) <fig-ajuste>

== Estabilidad frente a los filtros

Si la muestra son decaimientos, $tau$ no debería cambiar al endurecer los cortes. Al
aplicar los filtros de a uno, $tau$ pasa de 2183~ns con todos los pulsos tardíos de más de
30~mV (4210~eventos) a 2070~ns con la selección completa. Los dos cortes que más eventos
sacan se pueden mover (@fig-barridos). Con la amplitud mínima entre 30 y 60~mV, $tau$ queda
entre 2030 y 2096~ns. Por encima de 70~mV sube, pero el error crece igual de rápido y ningún
valor queda a más de 1,4~$sigma$ de 2197~ns. Con el $Delta t$ mínimo pasa algo distinto: si
la ventana empieza antes de 100~ns, $tau$ baja a unos 1950~ns, y si empieza en 200~ns sube a
2305~ns. Eso indica que hay eventos de más a $Delta t$ corto.

#figure(
  par2("reporte_tau_vs_amplitud", "reporte_tau_vs_dtmin"),
  caption: [Vida media ajustada al cambiar la amplitud mínima del pulso tardío (izquierda)
    y el $Delta t$ mínimo de la ventana (derecha), con el resto de los filtros sin cambios.
    Las barras son el error estadístico de cada ajuste; los puntos comparten eventos.],
) <fig-barridos>

= El exceso por debajo de 200 ns

Tomando como referencia una exponencial de 2197~ns normalizada entre 300 y 1400~ns, entre
60 y 200~ns se esperan 341~eventos y hay 388: un exceso de $47 plus.minus 20$ eventos
(2,4~$sigma$). Se analizaron dos explicaciones.

== ¿Pospulsos del fotomultiplicador?

Un pospulso es un pulso espurio que el fotomultiplicador emite un tiempo después de uno
real. Si fuera eso, tendría que aparecer en el mismo centellador que el pulso de entrada y
ser más frecuente cuanto mayor es ese pulso.

*Hay un efecto real, pero chico.* Promediando 12.000~trazas de muones que atraviesan,
alineadas en el pulso de entrada, se ve una cola que tarda entre 200 y 300~ns en volver a
cero (@fig-cola). Es de unos 4~mV a los 40~ns y de 1~mV a los 100~ns, para pulsos de unos
110~mV. Solo con pulsos de entrada de más de 250~mV la señal llega a bajar de −20~mV entre
60 y 200~ns (en el 3 al 10~% de las trazas, según el canal). Ninguno de los dos efectos
alcanza los 60~mV que se le piden al pulso tardío.

#figure(
  fig("reporte_cola_pulso_entrada", ancho: 62%),
  caption: [Promedio de 12.000~trazas de muones que atraviesan los tres centelladores,
    alineadas en el pulso de entrada (que queda fuera de escala, a la izquierda). La cola
    del pulso vuelve a la línea de base entre 200 y 300~ns después.],
) <fig-cola>

*Comparación con una muestra de control.* Cuando el electrón aparece en un centellador
que no tuvo pulso de entrada, no puede ser un pospulso. La proporción de eventos cortos
(60 a 200~ns) sobre largos (300 a 1400~ns) de esa muestra indica cuántos cortos reales
esperar en el resto; lo que sobre son pospulsos (@tab-control, @fig-control).

#figure(
  tabla(4, (left, right, right, right),
    ([Comparación], [Control], [Prueba], [Pospulsos estimados]),
    [Entrada en 1 y 2: electrón solo en 3 (control) o solo en 2], [58 / 251], [154 / 793], [$-29 plus.minus 28$],
    [Todos: electrón en centellador sin entrada (control) o con entrada], [67 / 296], [321 / 1653], [$-53 plus.minus 50$],
  ),
  caption: [Estimación de pospulsos. En "control" y "prueba" se indican eventos cortos
    sobre eventos largos. Una exponencial pura daría una proporción de 0,175; el control da
    0,23 ± 0,03 y la prueba, 0,19 ± 0,01.],
) <tab-control>

La estimación es compatible con cero: con 95~% de confianza hay menos de 82 pospulsos
entre los 321~eventos cortos (26~%). El exceso es igual o mayor justamente donde no puede
haber pospulsos.

#figure(
  fig("reporte_pospulsos_control", ancho: 66%),
  caption: [Distribución de $Delta t$ según si el centellador donde aparece el electrón
    tuvo (círculos) o no (rombos) pulso de entrada, normalizadas a la cantidad de eventos
    entre 300 y 1400~ns. Las barras son $sqrt(N)$ de cada bin. La curva es una exponencial de
    2197~ns con la misma normalización.],
) <fig-control>

*Chequeos del análisis.* Se miraron una por una las trazas de los eventos con $Delta t$
menor a 200~ns (58 con el electrón solo en el centellador~3 y 154 con el electrón solo en
el~2): en todas el segundo pulso está bien separado y tiene forma normal. Además se
agregaron pulsos artificiales a 1500~trazas reales para medir la eficiencia del buscador:
es de 96 a 97~% para cualquier $Delta t$, y solo baja a 93~% entre 60 y 70~ns cuando el
pulso está en el mismo canal que la entrada. Esa diferencia equivale a uno o dos eventos.

== Captura de $mu^-$ en el hierro

Un $mu^-$ frenado queda ligado a un átomo y tiene dos formas de desaparecer, que compiten:
decaer (da el electrón que se detecta) o ser capturado por el núcleo (no da electrón). El
muón desaparece con la suma de los dos ritmos. Para ver el electrón, el muón tiene que
haber sobrevivido a la captura hasta ese momento, así que los electrones que sí aparecen lo
hacen temprano. El resultado sigue siendo una exponencial, con una vida media más corta:

$ 1 / tau = 1 / tau_"decaimiento" + 1 / tau_"captura" $

La captura crece muy rápido con la carga del núcleo (@tab-captura). En el plástico casi no
se nota; en el hierro domina. Los $mu^+$ no son capturados.

#figure(
  tabla(3, (left, right, right),
    ([Caso], [Vida media observada], [Fracción que da electrón]),
    [$mu^+$ en cualquier material], [2197~ns], [100~%],
    [$mu^-$ en plástico (carbono)], [2026~ns], [92~%],
    [$mu^-$ en hierro], [unos 206~ns], [9~%],
  ),
  caption: [Vida media y fracción de decaimientos según la carga del muón y el material
    donde frena. Valores de bibliografía (Suzuki, Measday y Roalsvig, Phys. Rev. C 35,
    2212, 1987); el del hierro no se pudo verificar contra el texto original.],
) <tab-captura>

Esto afecta a la medición de dos maneras. En la cantidad: casi todos los $mu^-$ que frenan
en el hierro se pierden, lo que baja la estadística pero no cambia la forma. Y en la forma:
los pocos que decaen agregan una componente rápida a $Delta t$ corto, *solo en los eventos
donde el muón frenó en el hierro*. Su tamaño se puede estimar sin ajustar nada: con 44~% de
$mu^-$ y 9~% que decae, alrededor del 4~% de los decaimientos del grupo que frena en el
hierro debería estar en la componente de 206~ns, y ninguno en el grupo que frena en el
centellador~3. Esa es la prueba de la sección siguiente.

= Hierro contra plástico: una exponencial o dos

Los decaimientos se separan en dos grupos según hasta dónde llegó el muón:

- *Grupo 1, frena en el hierro:* pulso de entrada en los centelladores 1 y 2 y no en el 3
  (1888~eventos).
- *Grupo 2, cruza el hierro y frena en el centellador 3:* pulso de entrada en los tres
  (654~eventos).

Cada grupo se ajusta con una exponencial y con la suma de dos,
$f e^(-Delta t slash tau_"rápida") + (1 - f) e^(-Delta t slash tau_"lenta")$, con las dos
vidas medias y la fracción libres. La ventana empieza en 60~ns, porque la componente rápida
vive a $Delta t$ corto. La bondad se mide con el $chi^2$ sobre un histograma de 20~bins y
con la mejora de la verosimilitud (@tab-grupos, @fig-grupos).

#figure(
  tabla(3, (left, right, right),
    ([], [Grupo 1: frena en el hierro], [Grupo 2: frena en el centellador 3]),
    [Eventos], [1.888], [654],
    table.hline(stroke: 0.5pt),
    [Una exponencial: $tau$], [1817~ns (1638 a 2041)], [2454~ns (1963 a 3270)],
    [Una exponencial: $chi^2$ / grados de libertad], [17,3 / 18 (p = 0,50)], [14,6 / 18 (p = 0,69)],
    table.hline(stroke: 0.5pt),
    [Dos exponenciales: $tau$ rápida], [42~ns (hasta 105)], [203~ns (sin cota)],
    [Dos exponenciales: $tau$ lenta], [2025~ns (1789 a 2277)], [2875~ns (desde 2101)],
    [Dos exponenciales: fracción rápida], [1,5~% (0,8 a 3,0)], [2,2~% (0 a 37,5)],
    [Dos exponenciales: $chi^2$ / grados de libertad], [15,0 / 16 (p = 0,52)], [14,3 / 16 (p = 0,57)],
    table.hline(stroke: 0.5pt),
    [Mejora de la verosimilitud, $2 Delta ln L$], [3,20], [0,18],
    [p-valor de la mejora (simulado)], [0,07], [0,53],
    table.hline(stroke: 0.5pt),
    [Con $tau$ fijas en 206 y 2197~ns: fracción rápida], [3,7~%], [0~%],
    [Con $tau$ fijas: $-ln L$ respecto del mejor], [+0,86], [+0,20],
  ),
  caption: [Ajustes de una y de dos exponenciales en cada grupo, con $Delta t$ entre 60 y
    1400~ns. Entre paréntesis, el intervalo de 68~%; "sin cota" indica que los datos no
    limitan ese parámetro dentro del rango permitido (20 a 800~ns para la rápida, 800 a
    20.000~ns para la lenta).],
) <tab-grupos>

#figure(
  par2("casos_frena_en_hierro", "casos_cruza_el_hierro"),
  caption: [Distribución de $Delta t$ del grupo que frena en el hierro (izquierda) y del que
    frena en el centellador~3 (derecha), en bins de 67~ns con barras $sqrt(N)$, con el ajuste
    de una exponencial y el de dos exponenciales libres. Los valores están en la
    @tab-grupos.],
) <fig-grupos>

El resultado va en el sentido esperado: la segunda exponencial mejora el ajuste en el grupo
del hierro y no aporta nada en el del centellador~3. Pero la mejora en el hierro tiene un
p-valor de 0,07, que no llega al 0,05 habitual: es un indicio, no una detección.

== Cómo leer esos números

*El $chi^2$ reducido no sirve para comparar los dos ajustes.* Al agregar dos parámetros se
divide por 16 en lugar de 18, así que el cociente puede subir aunque el ajuste no empeore.
En el grupo~2 pasa de 0,81 a 0,90, y eso es exactamente lo que da no mejorar nada
(14,6 / 16 = 0,91). En simulaciones de una exponencial pura, el $chi^2$ reducido sube al
agregar la segunda en el 88~% de los casos. Además, cualquier valor entre 0,7 y 1,3 es
normal con 18 grados de libertad. Lo que se compara es cuánto baja el $chi^2$ o cuánto
mejora la verosimilitud.

*El p-valor de la mejora hay que simularlo.* La fórmula habitual (un $chi^2$ con dos grados
de libertad) da 0,20 para el grupo~1, pero no vale acá porque la fracción no puede ser
negativa. Simulando 200 muestras de una sola exponencial con la misma cantidad de eventos,
una mejora de 3,20 o más aparece en el 7~% de los casos (entre 4 y 7~% según la vida media
que se simule).

*Los 203~ns del grupo~2 no significan nada.* En muestras simuladas de 654~eventos con una
sola exponencial, sin componente rápida, el ajuste de dos exponenciales encuentra igual una
"rápida" de unos 208~ns en mediana, con alrededor del 1~% de los eventos. Es lo que
devuelve el ajuste cuando no hay nada que ajustar, y por eso su intervalo queda sin cota.

*Qué se esperaría si la componente del hierro está.* Simulando 2197~ns más un 3,9~% de
206~ns con 1888~eventos, una mejora como la observada o mayor aparece en el 31~% de los
casos, y una significativa (p menor a 0,05) solo en el 15~%. La vida media de una sola
exponencial baja a unos 1790~ns; en los datos del grupo~1 da 1817~ns. Es decir: aunque la
componente exista, con esta estadística lo más probable era no verla con claridad.

== El punto de partida y las vidas medias esperadas

Dejando todo libre, la vida media rápida del grupo~1 da 42~ns, lejos de los 206~ns
esperados. La verosimilitud es casi plana a lo largo de un valle, y el minimizador se
detiene en lugares distintos según de dónde arranque (@tab-partida).

#figure(
  tabla(2, (left, right),
    ([Dónde termina el ajuste (rápida, lenta, fracción)], [$-ln L$ respecto del mejor]),
    [42~ns, 2027~ns, 1,5~%], [0],
    [60~ns, 2000~ns, 1,7~%], [+0,11],
    [100~ns, 2200~ns, 2,5~%], [+0,43],
    [198~ns, 2266~ns, 3,9~%], [+0,83],
  ),
  caption: [Algunos de los puntos a los que llega el ajuste de dos exponenciales del
    grupo~1 desde 120 puntos de partida distintos.],
) <tab-partida>

Arrancando cerca de los valores esperados el ajuste se queda cerca de ellos, pero ese no es
el mínimo. Elegir el punto de partida para que dé lo esperado sería elegir la respuesta;
por eso se arranca desde muchos puntos y se toma el de menor $-ln L$. La forma correcta de
preguntar si los valores esperados son compatibles es fijarlos y medir cuánto empeora el
ajuste (últimas dos filas de la @tab-grupos). Para dos parámetros, hasta +1,15 es
compatible al 68~% y hasta +3,00 al 95~%:

- en el grupo del hierro, con 206 y 2197~ns el ajuste empeora solo 0,86 y la fracción
  rápida da 3,7~%, contra el 3,9~% predicho;
- en el grupo del centellador~3, la fracción rápida da cero.

Los dos grupos son compatibles con lo esperado. Lo que los datos no pueden hacer es medir
esas vidas medias por su cuenta.

== Lo que no cierra

- *Depende de la ventana.* Si el ajuste empieza en 100~ns en lugar de 60, el patrón se
  invierte: con la rápida fija en 206~ns, el grupo~1 da 0~% (0 a 3) y el grupo~2 da 5,6~%
  (1 a 10), con una mejora de 1,0~$sigma$. Lo que se ve depende de unos pocos eventos entre
  60 y 100~ns.
- *La vida media rápida libre da poco.* En las simulaciones con 206~ns verdaderos, el
  ajuste devuelve entre 93 y 554~ns en el 68~% de los casos; 42~ns queda por debajo.
- *La predicción es incompleta.* Solo cuenta electrones de decaimiento, no los neutrones y
  rayos gamma que emite el núcleo al capturar al muón, que también podrían dar pulsos.

== Por qué los datos no alcanzan

+ *La ventana es más corta que una vida media.* Entre 0,1 y 1,4~µs la exponencial lenta
  cae solo de 1 a 0,55: es casi una recta. Lo que distingue una exponencial de la suma de
  dos es la curvatura, y en un tramo tan corto casi no se ve.
+ *Las vidas medias del plástico son casi iguales.* Para separar $mu^+$ de $mu^-$ en
  plástico habría que distinguir 2197 de 2026~ns, un 8~%, y el error de una sola vida media
  es de 11~%. Fijando esas dos vidas medias y dejando libre la proporción, la verosimilitud
  cambia 0,11 entre "todo 2197" y "todo 2026", y la fracción queda entre 0 y 100~%.
+ *La componente del hierro es chica y está donde hay menos datos.* Es alrededor del 4~%
  de los eventos, toda por debajo de unos 400~ns.

Sobre la muestra completa, sin separar grupos, pasa lo mismo (@tab-dosexp).

#figure(
  tabla(3, (left, right, right),
    ([Ventana], [Una exponencial], [Dos, con la rápida fija en 206~ns]),
    [60 a 1400~ns (2542~eventos)], [1948~ns], [lenta 2432~ns, rápida 3,7~% (1 a 6)],
    [100 a 1400~ns (2421~eventos)], [2070~ns], [lenta 2237~ns, rápida 1,3~% (0 a 4)],
    [200 a 1400~ns (2154~eventos)], [2305~ns], [lenta 2305~ns, rápida 0~% (0 a 1)],
  ),
  caption: [Ajustes con una y dos exponenciales sobre la selección completa, para tres
    ventanas de $Delta t$. Desde 200~ns no queda componente rápida.],
) <tab-dosexp>

= Dónde frena el muón y hacia dónde sale el electrón

En el grupo~1 el electrón aparece mucho más seguido hacia arriba que hacia abajo
(@tab-arriba). La amplitud que deja el muón en el centellador~1 dice qué tan lento venía y,
con la relación de Bethe-Bloch, hasta dónde llegó (sección siguiente).

#figure(
  tabla(4, (left, right, right, right),
    ([Dónde aparece el electrón], [Eventos], [Amplitud en el centellador 1], [Alcance del muón]),
    [Solo arriba (en el 2 o en el 1)], [1.415], [1,30], [unos 14~g/cm²],
    [Solo abajo (en el 3)], [341], [1,11], [unos 24~g/cm²],
    [En el 2 y en el 3], [132], [1,38], [unos 12~g/cm²],
  ),
  caption: [Grupo~1 según hacia dónde sale el electrón. La amplitud es la mediana del pulso
    de entrada en el centellador~1, relativa a la de los muones que atraviesan. El alcance
    se mide en plástico equivalente desde la entrada a la pila; los centelladores 1 y 2
    suman 9,3~g/cm².],
) <tab-arriba>

Los dos alcances caen más allá de los 9,3~g/cm² de las dos barras, o sea adentro del
hierro: los muones cuyo electrón sale hacia arriba frenan en la parte de arriba de la placa
y los que lo mandan hacia abajo, en la parte de abajo. Es lo esperable si el hierro absorbe
a los electrones y solo salen los que nacen cerca de una cara. El alcance tiene un error de
alrededor del 20~%.

Lo que no queda explicado es por qué hay cuatro veces más hacia arriba. La absorción afecta
igual a las dos caras, y que frenen algo más de muones en las primeras capas da una
diferencia del orden del 25~%. Una parte puede venir de muones que frenan al final del
centellador~2, donde el electrón se ve casi siempre.

= Energía y velocidad de los muones

== Qué mide la amplitud

La amplitud del pulso de entrada mide la energía que el muón deja en el centellador, no la
que trae. Un muón rápido (de GeV) es de ionización mínima: deja unos 2~MeV por cm de
plástico, unos 10~MeV en la barra, casi sin importar su energía. La forma de la
distribución de amplitudes (un pico con cola larga) viene de las fluctuaciones de esa
pérdida, del camino recorrido en la barra y de la colección de luz; no es una distribución
de energías. El borde izquierdo en los centelladores 1 y 2 es el umbral de disparo y el
pico en unos 493~mV, la saturación de la placa.

Lo que sí se ve es que los muones que frenan dejan más (@fig-amplitudes): son lentos, y un
muón lento pasa más tiempo cerca de cada electrón del material y le transfiere más energía.
En el último centellador que tocan, los candidatos a decaimiento tienen una amplitud
mediana 1,5 veces mayor en el centellador~2 (161 contra 108~mV) y 1,7 veces mayor en el~3
(148 contra 87~mV). Mirando solo el centellador~1, la fracción de muones que termina
frenando sube 16 veces entre el pico de los rápidos y el doble de esa amplitud, y cae por
encima de unas 2,5 veces el pico, porque un muón que deja más de 30~MeV ahí ya no llega al
centellador~2.

#figure(
  par2("amplitud_entrada_frenados", "fraccion_frenan_vs_amplitud_u1"),
  caption: [Izquierda: amplitud del pulso de entrada en el último centellador que tocó el
    muón, relativa a la más probable al atravesar, para los muones que atraviesan
    (histograma) y para los que frenan y decaen (puntos, barras $sqrt(N)$); las dos
    distribuciones tienen área~1. Derecha: fracción de muones que frenan y decaen según la
    amplitud en el centellador~1.],
) <fig-amplitudes>

Aun en el mejor tramo de amplitud solo 7 de cada 1000 muones frenan: el resto son rápidos
con una fluctuación grande. La amplitud sola no dice, evento por evento, qué energía tenía
el muón; hace falta además saber dónde frenó.

== Energía de los que frenan

Separando por dónde frenaron y leyendo la energía de la curva de Bethe-Bloch a partir de la
amplitud mediana en el centellador~1 (@fig-bethe):

#figure(
  tabla(5, (left, right, right, right, right),
    ([Dónde frenan], [Amplitud relativa], [Energía cinética], [Velocidad], [Alcance]),
    [En el hierro (grupo 1)], [1,27], [unos 60~MeV], [0,77~$c$], [15~g/cm²],
    [En el centellador 3 (grupo 2)], [1,10], [unos 85~MeV], [0,83~$c$], [25~g/cm²],
  ),
  caption: [Energía al entrar a la pila de los muones que frenan. El error estadístico es
    de 1 a 2~MeV, pero comparar medianas con pérdidas medias y no conocer la respuesta
    exacta del centellador agregan fácilmente un 20~%.],
) <tab-energia>

#figure(
  par2("bethe_bloch_barra", "espectro_electron_u3"),
  caption: [Izquierda: energía media que deja un muón en una barra de 4,5~cm según su
    energía (Bethe-Bloch), con los dos grupos de muones que frenan ubicados a la altura de
    su depósito medido. Por encima de unos 300~MeV la curva es casi plana. Derecha: amplitud
    del pulso tardío en el centellador~3 (puntos, barras $sqrt(N)$) y simulación de
    electrones de decaimiento que se escapan de la barra.],
) <fig-bethe>

Por encima de unos 300~MeV la curva es casi plana, y ahí está la gran mayoría de los
muones: por eso de la amplitud no se puede sacar la energía de los que atraviesan.

*Chequeo con los electrones.* Una simulación donde el electrón nace en un punto cualquiera
de la barra, sale en una dirección cualquiera y deja unos 2~MeV por cm hasta escaparse
reproduce el espectro de amplitudes del pulso tardío (@fig-bethe). La escala que ajusta es
8,6~mV/MeV, contra 8,1~mV/MeV del pico de los muones rápidos, que sale de otra física. De
esa simulación sale también que solo el 34~% de los electrones supera los 60~mV.

== Estimación del espectro

Ninguna medición da la energía muón por muón, pero tres datos juntos fijan los dos
parámetros del modelo más simple del espectro a nivel del mar, $d N slash d E prop (E_0 + E)^(-n)$
(Shukla y Sankrith, arXiv:1606.06907):

- *Muones que frenan en el centellador 3:* 531 decaimientos sobre 503.953 muones que
  atraviesan los tres. Corrigiendo por la ventana temporal (0,42), la eficiencia del
  electrón (0,34) y el sesgo del umbral de disparo (1,18), frena ahí el 0,62~% de los
  muones, con energías entre 79 y 90~MeV: 0,57 ± 0,12 de los muones por GeV.
- *Barrido angular:* el exponente $cos^(2,02 plus.minus 0,37) theta$ fija $n$ = 3,02 ± 0,37.
- *Tiempo de vuelo:* 2,41 ± 0,25~ns en 61,2~cm, o sea 0,85 ± 0,09~$c$ sin corregir el
  retardo entre canales. Sirve solo como chequeo: el espectro predice 2,08~ns.

#figure(
  tabla(3, (left, right, right),
    ([], [Con nuestros datos], [Bibliografía]),
    [$n$], [3,02 ± 0,37], [3,01 (Tsukuba), 3,06 (Durham)],
    [$E_0$], [3,4 (+1,1 / −0,9) GeV], [4,29 (Tsukuba), 3,87 (Durham)],
    [Energía mediana], [1,4 ± 0,3~GeV], [],
  ),
  caption: [Parámetros del espectro. $n$ sale del barrido angular y $E_0$, de la fracción
    de muones que frenan.],
) <tab-espectro>

#figure(
  fig("espectro_energia_muones", ancho: 68%),
  caption: [Fracción de muones por GeV en función de la energía cinética. El punto lleno
    es la medición con los muones que frenan en el centellador~3 (barra horizontal: franja
    de energía). El punto hueco es la misma cuenta para los que frenan en el hierro, que
    queda por debajo porque el hierro absorbe electrones y captura a los $mu^-$; no entra
    en el ajuste. La curva llena pasa por el punto lleno por construcción; la punteada es la
    parametrización de referencia.],
) <fig-espectro>

Con ese espectro, la mitad de los muones tiene más de 1,4~GeV y va a más de 0,998~$c$; el
10~% más lento va a menos de 0,93~$c$ y el 5~% más lento, a menos de 0,84~$c$. Lo medido es
el punto a 85~MeV, el exponente angular y las velocidades de los que frenan; la forma de la
curva entre 0,1 y 100~GeV es del modelo. No se tuvo en cuenta el techo y las losas sobre el
detector, que corren $E_0$ hacia arriba.

= Conclusiones

- *Vida media.* $tau = 2070^(+264)_(-211)$~ns con 2421 decaimientos entre 100 y 1400~ns,
  compatible con el muón libre (2197~ns) y con lo esperado en plástico (unos 2100~ns). El
  error, de 11~%, lo limita la ventana más que la cantidad de eventos.
- *Robustez.* El resultado no depende de qué filtros se apliquen ni de la amplitud mínima
  entre 30 y 60~mV.
- *Por debajo de 200~ns.* Hay un exceso de $47 plus.minus 20$ eventos. No son pospulsos (la
  estimación da $-53 plus.minus 50$) ni un artefacto del análisis.
- *Hierro.* Una segunda exponencial mejora el ajuste solo en los muones que frenan en el
  hierro (p-valor 0,07), y con las vidas medias fijas en 206 y 2197~ns la fracción rápida
  da 3,7~% contra 3,9~% predicho. Es compatible con la captura de $mu^-$, sin que se pueda
  afirmar: depende de los eventos entre 60 y 100~ns.
- *Ventana del ajuste.* Empezar en 200~ns deja afuera cualquier efecto rápido y da
  $tau = 2305^(+388)_(-291)$~ns. Es la opción conservadora, a costa de un error de 15~%.
- *Energía.* Los muones que frenan llegan con 60 a 85~MeV (0,77 a 0,83~$c$) y dejan hasta
  1,7 veces más que los rápidos. El espectro estimado da una energía mediana de 1,4~GeV.

= Cosas para acordarse

- *El disparo y los filtros miden cosas distintas.* El disparo actúa sobre la señal cruda
  y los filtros sobre la suavizada, que da picos 1,36 veces más bajos.
- *Mirar los eventos antes de creerle a un corte.* Las galerías de "rechazados por" cada
  filtro y de los eventos rápidos fueron lo que permitió descartar artefactos.
- *Un grupo de control vale más que un filtro.* Los eventos con el electrón en un
  centellador sin pulso de entrada no pueden ser pospulsos; compararlos con el resto
  alcanzó para descartarlos.
- *El $chi^2$ reducido no compara modelos.* Al agregar parámetros puede subir sin que el
  ajuste empeore. Se compara cuánto baja el $chi^2$ o cuánto mejora la verosimilitud.
- *El p-valor de agregar una componente hay que simularlo*, porque su fracción no puede
  ser negativa y la fórmula habitual no vale.
- *Con una verosimilitud plana, el ajuste depende del punto de partida.* Hay que arrancar
  de muchos y quedarse con el mejor. Un valor "razonable" puede ser simplemente el punto de
  partida.
- *Para saber si un valor esperado es compatible, fijarlo* y ver cuánto empeora el ajuste,
  en lugar de buscar que el ajuste libre lo devuelva.
- *Una señal que cambia con la ventana es débil.* El indicio del hierro aparece desde
  60~ns y no desde 100~ns.
- *La amplitud mide energía depositada, no energía del muón.*
- *Lo que quedó sin verificar:* las vidas medias del $mu^-$ en hierro, aluminio y cobre
  (de memoria), la proporción de cargas (se usó 56~% de $mu^+$) y el espesor de la placa de
  hierro, que no está anotado.

= Propuestas a futuro

== Con los datos que ya hay

- *Bajar la amplitud mínima a 40~mV.* Suma un 44~% de eventos (3496 en lugar de 2421) con
  $tau$ estable: $2093^(+220)_(-182)$~ns. Antes conviene revisar a ojo los eventos de 40 a
  60~mV con el visor del notebook.
- *Ajustar con la componente del hierro incluida*, con su vida media y su fracción fijas
  en los valores predichos, para usar toda la ventana desde 60~ns sin sesgar $tau$.
- *Estimar el fondo de coincidencias casuales* con los pulsos anteriores al disparo. Una
  primera cuenta da como mucho un 7~% de la muestra, que bajaría $tau$ a unos 1930~ns si
  fuera todo fondo.
- *Verificar los valores de bibliografía* y medir el espesor de la placa de hierro.

== Si se vuelve a medir

- *Una ventana de varias vidas medias.* Es lo que más ayuda. Con una ventana de 10~µs o
  más, el error relativo de $tau$ tiende a $1 slash sqrt(N)$: con los mismos 2421~eventos
  sería de un 2~% en lugar de 11~%, y además entrarían más del doble de decaimientos, porque
  hoy la ventana solo alcanza al 43~% de ellos. Con la DRS4 habría que ver si admite una
  frecuencia de muestreo más baja; si no, un osciloscopio de memoria larga o un
  discriminador con un contador de tiempo resuelven el problema.
- *Medir con y sin hierro.* Sin la placa casi no se ven decaimientos, pero no debería
  haber componente rápida; con ella sí. Es la forma directa de confirmar la captura.
- *Cambiar el material.* La vida media del $mu^-$ depende mucho del núcleo: según la
  bibliografía, unos 860~ns en aluminio, 206~ns en hierro y 160~ns en cobre (valores a
  verificar). Con una ventana larga, comparar materiales permitiría ver esa dependencia.
- *Variar el espesor de absorbente* entre los módulos y contar cuántos muones frenan en
  cada configuración: es la forma de medir el espectro de energía a bajas energías.
- *Guardar todos los disparos* y no solo los seleccionados en línea, al menos en una
  corrida, para medir las eficiencias de disparo y selección en lugar de estimarlas.
- *Más estadística.* Se obtienen unos 3 decaimientos por hora. Reducir el error a la mitad
  con la ventana actual pediría cuadruplicar el tiempo de medición, a unos 136~días en
  total.

= Datos y código

Todo está en el repositorio `olistarosta/labo-5-2026`, en la carpeta `muones/`.

- `data/`: trazas crudas (se bajan del _release_ `muones-decaimientos-v1`), con su
  `README.md` y el manifiesto.
- `analisis/resumen_decaimientos.ipynb`: notebook corto con lo más importante.
- `analisis/decaimientos.py`: búsqueda de pulsos, variables por evento, filtros y ajustes.
- `analisis/decaimientos.ipynb`: notebook completo, con el visor de eventos.
- `analisis/reporte_filtros.py`: tablas y figuras de las secciones 3 a 5.
- `analisis/dos_casos.py`: los dos grupos, con una y dos exponenciales (sección 6).
- `analisis/galeria_rapidos.py`: trazas de los eventos con $Delta t$ menor a 200~ns.
- `analisis/energia_muones.ipynb`: energía y velocidad (sección 8).
- `salida_decaimientos/`: figuras en PDF y la lista de candidatos (`candidatos.csv`).
- `informe/`: este informe (`informe_decaimientos.typ`) y la figura de ejemplo.
