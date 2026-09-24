# Cómo funciona el análisis de imágenes

Este documento explica qué le pasa a cada par de cuadros desde que sale de la cámara hasta
que se convierte en un campo de velocidades y en la posición del centro del vórtice. Las
figuras son de la grabación `DINAMICA` (palangana con agua, 1600 rpm), cuadros 400 y 402.

## El recorrido completo

Se toman **dos cuadros** separados por 2 intervalos de cámara (1/15 s a 30 fps) y se hacen
dos cosas en paralelo:

```
                      ┌─ proyección del glitter ─► CLAHE ─► PIV (3 pasadas con deformación) ─► vectores
cuadro en color ──────┤
                      └─ proyección de la partícula ─► umbral relativo ─► centroide ─► centro del vórtice
```

- Arriba: **cómo se mueve el agua**, con el glitter como marcador (PIV).
- Abajo: **dónde está el centro del vórtice**, siguiendo la partícula amarilla que flota en él.

---

## 1. Proyección de color: de tres canales a una sola intensidad

La cámara da tres números por píxel (rojo, verde, azul). El PIV y el seguimiento necesitan
uno solo: "cuánto glitter hay acá" o "cuánta partícula hay acá".

La idea es pensar cada color como un **punto en el espacio RGB**. Se marcan tres colores con
un click: el del agua (fondo), el del glitter y el de la partícula. Después, para cada píxel:

$$\text{valor} = (\text{color del píxel} - \text{fondo}) \cdot \vec w$$

es decir, cuánto se aparta el píxel del agua **en la dirección** $\vec w$. Esa dirección se
elige con dos condiciones:

1. **El color buscado da 1**: $(\text{glitter} - \text{fondo}) \cdot \vec w = 1$.
2. **El otro color da 0**: $\vec w$ es perpendicular a $(\text{partícula} - \text{fondo})$.

El agua da 0 por construcción. Los valores negativos (más oscuro que el agua) se llevan a 0.

**Ejemplo con estos datos.** El glitter es gris y la partícula amarilla. Como el amarillo casi
no tiene azul, la dirección que ignora la partícula resulta ser casi el **canal azul**:
$\vec w_{glitter} \approx (0.07,\ -0.05,\ 1.36)$. Para la partícula pasa lo contrario:
$\vec w_{partícula} \approx (0.45,\ 0.59,\ -1.04)$, o sea "rojo + verde − azul", que es
justamente lo que distingue al amarillo.

![proyección de color](figuras_analisis/proyeccion.png)

Así la partícula **desaparece** de la imagen que usa el PIV (si no, el PIV la seguiría como si
fuera glitter), y el glitter desaparece de la imagen que se usa para encontrar la partícula.

---

## 2. CLAHE: igualar el contraste por zonas

El brillo del glitter no es parejo: hay zonas más iluminadas que otras. En la correlación, las
zonas brillantes pesan mucho más que las oscuras y dominan el resultado.

**CLAHE** (*Contrast Limited Adaptive Histogram Equalization*) corrige eso: divide la imagen en
zonas de unos 64 px y en cada una estira el contraste para que ocupe todo el rango. El
"*contrast limited*" pone un tope a ese estiramiento, para no amplificar el ruido donde casi no
hay nada. Es el mismo preprocesamiento que usa PIVlab por defecto.

![CLAHE](figuras_analisis/clahe.png)

Antes de CLAHE, la proyección se escala para que el percentil 99.9 quede en el máximo y se pasa
a 8 bits. Hay una opción para saltear CLAHE (`Imagen del PIV = lineal`), pero medido sobre la
simulación, CLAHE da mejor resultado: error de 0.148 px contra 0.160 px.

---

## 3. PIV: medir cuánto se movió cada pedacito de imagen

### 3.1 La correlación de una ventana

Se toma un cuadrado de la imagen (una **ventana de interrogación**) en el cuadro A y el mismo
cuadrado en el cuadro B. Si se corre la ventana B un poco hacia cada lado y se mira cuánto se
parece a la A, el corrimiento donde más se parecen es **cuánto se movió el glitter** adentro
de esa ventana.

Probar todos los corrimientos uno por uno sería lento. La **transformada de Fourier (FFT)** hace
todos a la vez: la correlación de A con B es

$$C = \mathcal{F}^{-1}\left[\ \overline{\mathcal{F}(A)} \cdot \mathcal{F}(B)\ \right]$$

y el pico de $C$ está en el desplazamiento.

![correlación de una ventana](figuras_analisis/correlacion.png)

En este ejemplo el pico está 9 px hacia arriba del centro: el glitter de esa ventana se movió
9 px entre los dos cuadros. Tres detalles:

- **Precisión menor a un píxel.** Al pico y a sus dos vecinos en cada dirección se les ajusta
  una gaussiana, y su máximo da el desplazamiento con decimales (del orden de 0.05 px).
- **Rango de búsqueda.** Solo se miran corrimientos de hasta **1/4 de la ventana**: más lejos,
  demasiado glitter entra o sale de la ventana y la correlación deja de ser confiable.
- **Corrección de borde.** La FFT "da la vuelta" a la imagen por los bordes. Eso hace que los
  corrimientos grandes pesen menos que los chicos, y se corrige antes de buscar el pico (si no,
  el resultado queda sesgado hacia cero).

### 3.2 Por qué varias pasadas

Hay un compromiso con el tamaño de la ventana:

- una **ventana grande** (128 px) ve desplazamientos grandes, pero promedia mucho: no resuelve
  detalles del flujo;
- una **ventana chica** (32 px) resuelve detalles, pero su rango de búsqueda es chico (8 px).

La solución es hacer **tres pasadas**: 128, 64 y 32 px, con ventanas que se solapan a la mitad.
Cada pasada usa el resultado de la anterior como **predictor**: ya sabe aproximadamente cuánto
se movió cada zona, y solo tiene que medir lo que falta.

### 3.3 La deformación de ventana

Con el predictor $\vec d(\vec p)$ de la pasada anterior, antes de correlacionar se **deforman**
las dos imágenes:

- el cuadro A se lee en $\vec p - \vec d/2$,
- el cuadro B se lee en $\vec p + \vec d/2$.

Cada imagen se mueve la mitad del recorrido, en sentidos opuestos, así que cada destello queda
en el punto medio de su trayectoria en las dos. Si el predictor es bueno, **A y B deformadas
coinciden** y la correlación da un pico muy cerca de cero: lo que mide es solo la corrección.

Se deforma píxel por píxel, no se corre la ventana entera. Eso importa en un vórtice: adentro de
una misma ventana, un lado gira distinto que el otro. Ningún corrimiento único alinea toda la
ventana, pero deformarla sí. La interpolación de los píxeles (Lanczos 8×8) está elegida porque
casi no introduce sesgo: medido, 0.001 px, contra 0.03–0.06 px de las interpolaciones más
simples.

![las tres pasadas](figuras_analisis/iteraciones.png)

En la figura, la grilla naranja son las ventanas de cada pasada, en el lugar donde se leen: en
la primera la grilla es recta (todavía no hay predictor), y en las siguientes se curva siguiendo
el giro. Las flechas son lo que midió cada pasada (ampliadas ×2).

### 3.4 Validación: descartar lo que no se midió bien

Después de cada pasada se revisa cada vector. Se descartan:

- las ventanas **sin textura** (casi sin glitter): la correlación ahí es ruido;
- los vectores cuyo pico quedó **en el borde del rango de búsqueda**: el desplazamiento real
  era más grande de lo que la ventana podía ver;
- los que no se parecen a sus **8 vecinos** (test de la mediana normalizada, el de PIVlab);
- los que se apartan demasiado del promedio de todo el cuadro.

Lo descartado **no se inventa**: en el resultado queda como `NaN`, con `valido = 0`. Entre
pasadas sí se rellena con el promedio de los vecinos, pero solo porque la pasada siguiente
necesita un predictor en toda la grilla; ese relleno no sale al CSV.

### 3.5 De píxeles a mm/s

El desplazamiento se multiplica por la escala (mm por píxel, de la regla marcada en la
calibración) y se divide por el tiempo entre los dos cuadros del par.

---

## 4. La partícula del centro

La partícula amarilla flota en el centro del vórtice, así que su posición da **dónde está el
centro** en cada instante. Ese centro se mueve (en DINAMICA, hasta unos 40 px), y hace falta
conocerlo para pasar las velocidades a coordenadas polares.

En cada uno de los dos cuadros del par:

1. Se toma la **proyección de la partícula** (sección 1): la partícula brilla y todo lo demás
   queda en negro.
2. Se busca el **pico**: el punto más brillante adentro de la zona de análisis.
3. Se umbraliza a **la mitad de ese pico**: quedan los píxeles que brillan al menos la mitad.
4. De las manchas que quedan, se elige la **más cercana al centro del par anterior**.
5. Su **centroide pesado por el brillo** da la posición con precisión menor a un píxel.

![detección de la partícula](figuras_analisis/particula.png)

El centro del par es el promedio de las posiciones en los dos cuadros: corresponde al mismo
instante que la velocidad. Si el pico no se destaca del resto de la imagen (al menos 4 veces
más brillante), la partícula no está a la vista y ese par queda marcado como no seguido.

**Por qué el umbral es relativo** (la mitad del pico, y no un valor fijo). La proyección vale 1
en el color que se marcó, pero en el agua la partícula se ve más pálida que ese color: acá su
pico llega a 0.27. Con un umbral fijo de 0.5 no se la encontraba nunca. Con el umbral relativo,
en todos los cuadros de DINAMICA queda una sola mancha limpia, de 90 a 200 px, sin importar el
color exacto que se haya marcado ni cuánta luz haya.

---

## 5. Qué tan bien funciona

Todo se probó con una **simulación**: glitter en un vórtice de Rankine con velocidad conocida,
una partícula en el centro y un centro que se desplaza. Así se puede comparar lo medido con la
respuesta exacta.

| qué | resultado |
|---|---|
| TouchDesigner contra el código de referencia en Python, sobre los mismos cuadros | diferencia mediana 0.003 px |
| centro seguido contra el centro real | error mediano 0.05 px, partícula encontrada en todos los pares |
| velocidad angular $\omega(r)$ contra la del vórtice simulado | diferencias menores al 1 % |
| vectores descartados | 0.5 % |

En las mediciones reales no hay respuesta exacta. El indicador más útil es el **porcentaje de
vectores descartados**: si pasa del ~40 %, suele ser porque el vórtice gira demasiado rápido
para la ventana. La solución es volver a analizar con `Separación de cuadros = 1`, sin
necesidad de volver a filmar.

Los detalles de cada parámetro están en [README.md](README.md).
