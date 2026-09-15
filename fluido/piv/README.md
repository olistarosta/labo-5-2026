# PIV del vórtice con glitter

Todo está en `TDS GLITTER/TDS GLITTER.toe`, en dos componentes con su propia interfaz
(click derecho sobre el componente → *Open Viewer*):

| componente | qué hace |
|---|---|
| `/project1/grabar` | graba **una o dos cámaras** como cuadros TIFF sin pérdida + la hora de cada cuadro |
| `/project1/analizar` | analiza una grabación con PIV y guarda **un CSV por cámara** en `data glitter` |

Cámara **A** = vista superior (da x, y). Cámara **B** = vista lateral (da la altura z).

## 1. Grabar

1. Página **Cámara A** / **Cámara B**: `Fuente` (Apagada / Cámara / Simulación),
   `Dispositivo` y `Espejar horizontal` (si la cámara entrega la imagen espejada; si no,
   el giro sale invertido).
2. Página **Grabar**: `Nombre`, `Agitador (rpm)` y **GRABAR**. Otra vez GRABAR para parar.

Mientras graba, el estado muestra cuadros y fps de cada cámara. `Carpeta de grabaciones`
vacía = `Videos\vortice_glitter` (fuera de OneDrive: ~1 MB por cuadro).

    <fecha>_<nombre>/
      grabacion.json                qué grabó cada cámara, rpm
      camara_A/cuadros/*.tif        imágenes
      camara_A/tiempos.csv          indice, t_s
      camara_B/...

## 2. Analizar

1. **Grabación**: elegí esa carpeta. Se carga sola.
2. **Cámara que ves y calibrás**: A o B. Para cada cámara:
   - `Marcar con el mouse = Punto A de la regla` → click; pasa solo a B → click.
     Cargá la `Distancia real A-B` en la página de esa cámara.
   - `Marcar = Centro del vórtice` → click; pasa solo al borde del ROI → click.
   El estado dice qué falta calibrar.
3. **ANALIZAR TODO**. Analiza todos los pares de todas las cámaras y guarda en `data glitter`:

       <grabación>_camaraA_campo.csv
       <grabación>_camaraB_campo.csv
       <grabación>_parametros.json     calibración de cada cámara, dt, método

   Nunca pisa un análisis anterior (agrega `_2`, `_3`). **Detener** corta y guarda lo hecho.

`Par a inspeccionar` recalcula un par al instante; `Mostrar` elige qué etapa se ve
(imagen con vectores, cuadro original, CLAHE, imagen deformada).

Método, como PIVlab: CLAHE → FFT en 3 pasadas (128, 64, 32 px, 50 % de solapamiento)
con deformación de ventana (Lanczos) → pico gaussiano de 3 puntos → validación
(mediana local 3, desvío estándar 8) e interpolación de los descartados.

### El CSV (uno por cámara)

| columna | |
|---|---|
| `frame` | primer cuadro del par |
| `t_s` | tiempo del centro del par, desde el primer cuadro de esa cámara |
| `x_px`, `y_px` | centro de la ventana; origen abajo-izquierda, y hacia arriba |
| `x_mm`, `y_mm` | ídem en mm (en la cámara B, `y` es la altura z) |
| `u_mm_s`, `v_mm_s` | velocidad (en px/s si esa cámara no está calibrada) |
| `interpolado` | 1 = la validación lo descartó y se reemplazó por sus vecinos |

```python
import pandas as pd
a = pd.read_csv('..._camaraA_campo.csv')
un_par = a[a.frame == 100]
medio = a.groupby(['x_mm', 'y_mm'])[['u_mm_s', 'v_mm_s']].mean().reset_index()
```

**Paso temporal**: dt fijo por cuadro (duración / intervalos), porque la hora de cada
cuadro viene redondeada a ~4 ms; un cuadro perdido se detecta y ese par usa 2·dt.

## Validación

- **Análisis**: la red de TD contra `piv_core.py` sobre los mismos cuadros da diferencia
  mediana 0.003 px; ω(r) contra el vórtice inyectado, −0.1 % a +0.4 % entre r = 60 y 260 px.
- **Grabación con dos webcams reales** (640×480): ~29 fps efectivos por cámara, 6–9 cuadros
  perdidos cada 300 (detectados y compensados en el análisis).

    <python de TD> regresion_td.py <carpeta_grabación> <..._parametros.json> [pares]

## Archivos de esta carpeta

`piv_core.py` motor de referencia (lo copia el DAT `piv_numerico` de TD) ·
`piv_sim.py` simulador con campo conocido · `regresion_td.py` prueba de TD contra la referencia.
