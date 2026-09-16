# Vórtice — seguimiento de partícula (Labo 5, Fluidos)

Proyecto de TouchDesigner (`DETECCION FLUIDO.toe`) para medir la posición de la "pelotita"
trazadora en el vórtice del agitador magnético, cuadro a cuadro, y exportar su trayectoria
en unidades físicas.

**Dos canales de cámara independientes**, cada uno con su propia calibración:

| canal | qué mira | qué da |
|---|---|---|
| **A** | vista superior (desde arriba del recipiente) | `x`, `y` del plano horizontal |
| **B** | vista lateral (a través de la pared) | `z` (altura) |

El programa **no** calcula el centro del vórtice: entrega posiciones absolutas y el centro
se ajusta después, a partir de la trayectoria.

Todo vive en `/project1/vortex_tracker`. Abrí su panel (click derecho → *Open as Panel*) y
trabajás desde ahí: no hace falta tocar la red. El proyecto corre a **30 fps**.

## Flujo de trabajo

La interfaz son dos columnas idénticas, una por cámara (arriba lo que ve, abajo sus
controles), y una columna a la derecha con los controles de la medición.

1. **Origen** — `Simulación` genera un vórtice sintético para practicar sin equipo;
   `Cámara en vivo` para el experimento (si da error, cambiá `Driver de cámara` y pulsá
   *Refrescar cámaras*); `Archivo de video`; `Sin usar` apaga el canal si solo medís con
   una cámara.

   Si la webcam entrega la imagen espejada, activá `Espejar horizontal`: si no, el sentido
   de giro medido sale al revés del real.

2. **Detección** — poné `Vista previa = Máscara binaria` y ajustá `Umbral`,
   `Niveles negro/blanco` e `Invertir` (partícula oscura sobre fondo claro) hasta que quede
   **sólo la partícula** en blanco. `Área mínima` descarta clusters chicos. Si el borde del
   recipiente mete reflejos, poné un `radio` distinto de 0 en `Región de interés` y marcá
   su centro con click.

   El algoritmo se queda con **el cluster conexo más grande** y calcula su **centro de masa
   pesado por luminancia** — los reflejos puntuales y el ruido no lo mueven.

   **Usá el color, no sólo el brillo.** Antes del umbral, la imagen pasa por una proyección
   lineal `R = pesos · rgb + offset`; con los pesos por defecto eso es exactamente
   luminancia. El botón **`Ajustar color a la partícula`** calcula los pesos solos: busca
   la dirección del espacio RGB que **más separa** la partícula del fondo, tomando como
   partícula el cluster más grande y como fondo todo el resto (reflejos incluidos), y
   ajusta ganancia y offset para que el fondo quede en 0 y la partícula en 1.

   Sirve muchísimo cuando la partícula y el fondo tienen brillo parecido pero color
   distinto — por ejemplo el montaje de tinta fluorescente con luz azul que menciona la
   guía. Para volver a luminancia, click derecho sobre `Pesos color` → reset.

3. **Calibración — hay que hacerla en las dos cámaras, por separado.**
   `Marcar con el mouse = Punto A`, clickeá un extremo de una distancia conocida, después el
   otro extremo (avanza solo a B), y cargá esa distancia en `Distancia real A-B (mm)`. La
   escala sale sola.

   > Ojo: cada cámara ve un plano distinto. La escala de la vista superior se mide en el
   > plano horizontal (p. ej. el diámetro interno del recipiente); la de la vista lateral
   > se mide **en el plano vertical**, a la misma distancia de la cámara que la partícula.

4. **Medición** — nombre, `Guardar datos` y/o `Guardar video`, y **GRABAR**.

## El tiempo

- **Cámara en vivo**: `frame_timestamp` de la capa de captura — el instante en que llegó
  ese cuadro, **corriente arriba de todo el procesamiento de TD**. Además la estampa es la
  compuerta: si no cambió, TD está releyendo la misma imagen y esa fila no se guarda, así
  que hay **una fila por cuadro de cámara, sin repetidos**.
- **Modo archivo**: `t_s = cuadro / fps del video`. Exacto, lo define el archivo.
- **Simulación**: reloj real monótono (`time.perf_counter`).

Nunca se usa el contador de cuadros de TD: ése se queda clavado si TD se pausa y
comprimiría el tiempo sin avisar (medido: 78 s reales sin avanzar un cuadro).

Cada cámara lleva su propia estampa (`t_s` la A, `t_b_s` la B, mismo cero) porque son
independientes y no están sincronizadas entre sí.

**Limitación**: `frame_timestamp` viaja por un CHOP en `float32` con base de reloj de
sistema, así que queda cuantizado a ~1/256 s = 3.9 ms, y empeora cuanto más uptime tenga
la máquina. El `dt` alterna entre 31.25 y 35.16 ms alrededor del período real de 33.33; el
promedio sobre varios cuadros da el período correcto.

Hay una **segunda red** contra un problema que la estampa no ve: en poca luz la cámara
alarga la exposición y la capa de captura re-entrega el mismo contenido con timestamp
nuevo. Si la medición sale idéntica bit a bit a la anterior, la imagen era la misma y esa
fila se descarta. Queda contado en `imagenes_repetidas_descartadas`: compararlo con
`cuadros` te dice la tasa **real** de imágenes distintas, que puede ser mucho menor que
`camara.a.fps`. Si esa tasa es baja, necesitás más luz o menos tiempo de exposición.

El JSON reporta `dt_mediano_s`, `dt_maximo_s`, `dt_desvio_s`, `saltos` y
**`cortes_de_reloj`**. Este último cuenta discontinuidades (`dt ≤ 0` o `dt > 1 s`): si el
stream USB de la cámara se reinicia en medio de una medición, su reloj vuelve a cero y la
base de tiempo queda rota. Si sale distinto de cero, el estado te lo avisa al terminar y
esa medición hay que descartarla o cortarla en ese punto.

**La barra de tiempo de TD representa la medición.** La línea de tiempo se estira a
`Duración max. (min)` (10 minutos por defecto, contra los 20 s que trae TD de fábrica) y
GRABAR la rebobina, así que la posición del cabezal es el tiempo transcurrido de *esta*
medición y el número de cuadro que muestra TD es el mismo que la columna `frame` del CSV.
Es sólo presentación: el tiempo que se guarda no sale de ahí, así que rebobinar no puede
corromper los datos. Si la medición se pasa del largo, la barra da la vuelta y el estado
lo avisa.

## ⚠ TouchDesigner tiene que estar corriendo

Si la línea de tiempo está detenida no se graba nada, sin ningún error visible. **TD pausa
solo cuando minimizás la ventana** (preferencia *Stop Playing when Minimized*, activada por
defecto). Apretar GRABAR reanuda la línea de tiempo y avisa en el estado, pero durante una
medición dejá la ventana visible, o apagá esa preferencia.

## Salida

    data/2026-09-10_143012_<nombre>/
      <nombre>_track.csv    una fila por cuadro
      <nombre>_meta.json    calibración y parámetros de detección DE CADA CÁMARA
      <nombre>_A.mov        video limpio de la cámara A (sólo si `Guardar video`)
      <nombre>_B.mov        video limpio de la cámara B

Columnas del CSV:

| columna | unidad | qué es |
|---|---|---|
| `frame` | — | cuadro de la línea de tiempo de TD (del video, en modo archivo) |
| `t_s` | s | tiempo de la cámara A desde el inicio de la medición |
| `t_b_s` | s | ídem para la cámara B (mismo cero; igual a `t_s` si B no es cámara) |
| `a_ok` | 0/1 | la cámara A detectó la partícula en ese cuadro |
| `a_x_px`, `a_y_px` | px | posición en la imagen de A |
| `a_x_mm`, `a_y_mm` | mm | idem, con la escala de A → **plano horizontal** |
| `a_area` | px | tamaño del cluster (en la resolución de análisis) |
| `b_ok` | 0/1 | la cámara B detectó la partícula |
| `b_h_px`, `b_z_px` | px | posición en la imagen de B |
| `b_h_mm` | mm | coordenada **horizontal** en la vista lateral |
| `b_z_mm` | mm | **altura (z)** |
| `b_area` | px | tamaño del cluster |

Convención: origen abajo-izquierda de la imagen **de cada cámara**, x a la derecha, y hacia
arriba. Las posiciones son absolutas respecto del borde de la imagen, así que el centro del
vórtice y las velocidades se calculan en el análisis:

```python
import numpy as np
d = np.genfromtxt('..._track.csv', delimiter=',', names=True)
ok = d['a_ok'] == 1
x, y, t = d['a_x_mm'][ok], d['a_y_mm'][ok], d['t_s'][ok]

# centro del vórtice: ajuste de circunferencia por mínimos cuadrados (Kåsa)
A = np.column_stack([2*x, 2*y, np.ones(len(x))])
cx, cy, c = np.linalg.lstsq(A, x*x + y*y, rcond=None)[0]

r = np.hypot(x - cx, y - cy)
vx, vy = np.gradient(x, t), np.gradient(y, t)
v_tan = (-vx*(y - cy) + vy*(x - cx)) / r     # velocidad tangencial vs radio
v_rad = ( vx*(x - cx) + vy*(y - cy)) / r
```

En modo archivo las dos cámaras usan el mismo índice de cuadro, así que los dos videos
tienen que estar sincronizados y al mismo fps.

## Precisión

Contrastado contra la simulación, donde la posición verdadera se conoce exactamente: el
error del centroide es de **0.03 – 0.15 px** analizando a 480×360 una fuente de 640×480, es
decir sub-píxel incluso respecto de la grilla de análisis reducida. Con dos canales activos
a 30 fps el uso medido es menos del 1 % del presupuesto de cuadro.

La documentación interna (estructura de la red, algoritmo y las trampas de TD que hubo que
resolver) está en el DAT `agents_md`, adentro del componente.
