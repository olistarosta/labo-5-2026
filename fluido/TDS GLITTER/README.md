# Vórtice con glitter — seguimiento de trazadores (Labo 5, Fluidos)

Proyecto de TouchDesigner (`TDS GLITTER.toe`) para medir el campo de velocidades de un
vórtice sembrando **glitter** en el agua. En cada cuadro detecta **todos** los destellos y
exporta sus posiciones; el perfil ω(r) sale del análisis posterior.

Todo vive en `/project1/vortex_tracker`. Abrí su panel (click derecho → *Open as Panel*).

## Qué cambió respecto del método de una partícula

Con glitter hay cientos de trazadores por cuadro, así que la regla vieja —quedarse con el
cluster conexo más grande— pasa a ser exactamente lo que **no** hay que hacer. Ahora:

- Se aceptan **todas** las componentes conexas cuya área cae en `Área de un destello
  mín/máx`, hasta `Máximo de destellos`. El máximo importa tanto como el mínimo: descarta
  reflejos grandes y grumos de destellos pegados, que no son un trazador solo.
- Los centroides de todas las componentes salen de **una sola pasada** con `np.bincount`
  sobre las etiquetas: O(píxeles), sin bucle por destello.
- **No hay seguimiento entre cuadros.** Asociar quién es quién con glitter denso es frágil,
  así que cada cuadro es una foto independiente. La asociación queda para el análisis, que
  puede usar métodos de conjunto (ver abajo).
- El **centro del vórtice** es ahora un parámetro por cámara que se marca con el mouse
  (`Marcar con el mouse → Centro del vórtice`). No se le resta a los datos: va al JSON como
  calibración, para que el análisis calcule r y θ de cada destello.

## Flujo de trabajo

Dos columnas idénticas, una por cámara (arriba lo que ve, abajo sus controles), y una
columna a la derecha con la medición.

1. **Origen** — `Simulación` genera un vórtice de glitter sintético para practicar sin
   equipo; `Cámara en vivo` para el experimento; `Archivo de video`; `Sin usar` apaga el
   canal.
2. **Detección** — poné `Vista previa = Máscara binaria` y ajustá `Umbral` y `Niveles`
   hasta que queden **sólo los destellos**. `Área de un destello mín/máx` filtra por
   tamaño. Si el borde del recipiente mete reflejos, poné un `radio` distinto de 0 en
   `Región de interés`.

   **Usá el color.** El botón `Ajustar color a la partícula` calcula la dirección del
   espacio RGB que más separa los destellos del fondo, tomando como muestra los destellos
   que pasan el filtro de área y como fondo todo el resto.
3. **Calibración** — `Punto A` y `Punto B` sobre una distancia conocida en el plano del
   fluido, y la distancia real en mm. Después marcá el `Centro del vórtice`.
4. **Medición** — nombre, velocidad del agitador, y **GRABAR**.

## Salida

    data/2026-09-10_155532_<nombre>/
      <nombre>_track.csv    una fila por destello por cuadro
      <nombre>_meta.json    calibración, centro del vórtice y parámetros de cada cámara
      <nombre>_A.mov        video limpio (sólo si `Guardar video`)

| columna | unidad | qué es |
|---|---|---|
| `frame` | — | cuadro de la línea de tiempo de TD |
| `t_s` | s | tiempo desde el inicio de la medición, de **esa** cámara |
| `cam` | 0/1 | 0 = cámara A (vista superior), 1 = cámara B (vista lateral) |
| `x_px`, `y_px` | px | posición del destello en la imagen de esa cámara |
| `area` | px | tamaño del destello, en la resolución de análisis |
| `lum` | 0-1 | luminancia media del destello |

Origen abajo-izquierda de la imagen de cada cámara, x a la derecha, y hacia arriba. Para la
cámara B (vista lateral) `y_px` es la altura z.

Se guardan **píxeles, no mm**: `mm = px × escala`, y la escala de cada cámara está en el
JSON. Guardar el píxel crudo deja recalibrar después sin volver a medir, y con cientos de
destellos por cuadro cada columna de más cuesta megabytes.

⚠️ El archivo crece rápido: ~250 destellos × 30 fps ≈ 15 MB/min por cámara. Si molesta,
bajá `Máximo de destellos`.

## Cómo sacar ω(r) sin seguir destellos individuales

El vórtice es (casi) axisimétrico, así que no hace falta saber qué destello es cuál: se
correlaciona el patrón angular de un anillo de radio entre cuadros consecutivos.

```python
import numpy as np, glob, json
p = sorted(glob.glob('data/*/*_track.csv'))[-1]
d = np.genfromtxt(p, delimiter=',', names=True)
m = json.load(open(p.replace('_track.csv', '_meta.json'), encoding='utf-8'))
cx, cy = m['cam_a_vista_superior']['calibracion']['centro_vortice_px']

sel = d['cam'] == 0
fr, t = d['frame'][sel].astype(int), d['t_s'][sel]
r = np.hypot(d['x_px'][sel] - cx, d['y_px'][sel] - cy)
th = np.arctan2(d['y_px'][sel] - cy, d['x_px'][sel] - cx)

frames = np.unique(fr)
tof = {f: t[fr == f][0] for f in frames}
idx = {f: np.where(fr == f)[0] for f in frames}
NB = 720                                   # 0.5 grados
k = np.exp(-0.5 * (np.arange(-20, 21) / 4.0) ** 2); k /= k.sum()

def hist(sel):                             # histograma angular suavizado
    h = np.bincount(((th[sel] + np.pi) / (2*np.pi) * NB).astype(int) % NB,
                    minlength=NB).astype(float)
    return np.convolve(np.r_[h[-20:], h, h[:20]], k, 'same')[20:-20]

for lo, hi in [(20,50), (50,85), (85,120), (120,150), (150,180)]:
    ws = []
    for a, b in zip(frames[:-1], frames[1:]):
        sa = idx[a][(r[idx[a]] >= lo) & (r[idx[a]] < hi)]
        sb = idx[b][(r[idx[b]] >= lo) & (r[idx[b]] < hi)]
        if len(sa) < 8 or len(sb) < 8:
            continue
        c = np.fft.irfft(np.fft.rfft(hist(sb)) * np.conj(np.fft.rfft(hist(sa))), NB)
        kk = int(np.argmax(c))             # pico, con interpolacion sub-bin
        y0, y1, y2 = c[(kk-1) % NB], c[kk], c[(kk+1) % NB]
        den = y0 - 2*y1 + y2
        sh = ((kk + (0.5*(y0-y2)/den if den else 0) + NB/2) % NB) - NB/2
        dt = tof[b] - tof[a]
        if dt > 0:
            ws.append(sh / NB * 2*np.pi / dt)
    print('r %3d-%3d   omega = %.4f rad/s' % (lo, hi, np.median(ws)))
```

Después, `v_θ(r) = ω(r) · r · escala` te da el perfil de velocidad tangencial en mm/s, y ahí
se ve el núcleo de rotación rígida y la caída ~1/r del vórtice libre.

## Validación

La `Simulación` no es decorativa: es un vórtice de Rankine con parámetros conocidos
(`OMEGA0 = 2.2 rad/s`, `RCORE = 85 px`, centro en 320,240), con destellos que titilan y
derivan lentamente hacia adentro. Corriendo el análisis de arriba sobre 850 cuadros
grabados de esa simulación, el perfil se recupera dentro del **1.4 %**:

| r (px) | ω medido | ω inyectado | error |
|---|---|---|---|
| 20–50 | 2.1986 | 2.2000 | −0.1 % |
| 50–85 | 2.2057 | 2.2000 | +0.3 % |
| 85–120 | 1.5202 | 1.5129 | +0.5 % |
| 120–150 | 0.8599 | 0.8722 | −1.4 % |
| 150–180 | 0.5805 | 0.5838 | −0.6 % |

El anillo más externo (180–210) da +10 %: está truncado por el borde del sembrado, así que
los destellos se amontonan en su lado interno, donde ω es mayor. Es el sesgo esperado del
bin de borde — no lo uses.

## El tiempo

- **Cámara en vivo**: `frame_timestamp` de la capa de captura, el instante en que llegó ese
  cuadro, corriente arriba de todo el procesamiento de TD. La estampa además es la
  compuerta: si no cambió, TD está releyendo la misma imagen y esa fila no se guarda.
- **Archivo**: `t_s = cuadro / fps del video`. Exacto.
- **Simulación**: reloj real monótono (`time.perf_counter`).

Nunca se usa el contador de cuadros de TD. **Ese error es real y lo medimos**: la primera
versión del simulador avanzaba la física con `absTime` mientras el grabador medía con reloj
real, y como TD corría a 25 fps en vez de 30, el ω reconstruido salió un 15 % bajo — el
cociente 25.4/30 exacto. Con las dos puntas en reloj real el sesgo desaparece.

El JSON reporta `dt_mediano_s`, `dt_maximo_s`, `saltos`, `cortes_de_reloj` e
`imagenes_repetidas_descartadas`. Miralos antes de derivar nada.

## ⚠ TouchDesigner tiene que estar corriendo

Si la línea de tiempo está detenida no se graba nada, sin ningún error visible. **TD pausa
solo cuando minimizás la ventana** (preferencia *Stop Playing when Minimized*). Apretar
GRABAR reanuda la línea de tiempo, pero durante una medición dejá la ventana visible.

La documentación interna (estructura, algoritmo y las trampas de TD que hubo que resolver)
está en el DAT `agents_md`, adentro del componente.
