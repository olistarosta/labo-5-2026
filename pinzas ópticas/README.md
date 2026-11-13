# Pinzas ópticas — de video a CSV (Labo 5)

Se graba con el programa de Thorlabs (**ThorImageCAM** o **ThorCam**) y la platina se mueve con **Kinesis**, todo desde tu PC: ver **[docs/configurar_pc.md](docs/configurar_pc.md)**.

TouchDesigner hace una sola cosa: **PROCESAMIENTO** convierte el video en **trayectorias en CSV**.

El análisis físico (MSD, κ, Stokes) se hace después en VS Code, con los CSV (carpeta `analisis/`).

## Carpetas

```
pinzas ópticas/
  analisis/           análisis y visualización en VS Code: trayectorias.py + notebooks (uno por grabación)
  datos/              los CSV de trayectorias de cada grabación (salida de PROCESAR / _lote.py)
  figuras/            las figuras de cada análisis, en PDF (va al informe) y PNG (no se sube)
  pinzas/             paquete Python (lo usan TouchDesigner y el análisis)
  tests/              tests del paquete y de los bloques   →  correr_tests.py
  docs/               configurar_pc (cámara, Kinesis, PC), CONTRATO (unidades y formatos), LECCIONES, PLAN
  GRUPOS ANTERIORES/  videos de referencia
  TouchDesigner/      NO va a git (.gitignore): el .toe, sus Backup, el código de los nodos
                      (bloques/, construir_bloques.py), herramientas/ y lib/
```

`TouchDesigner/` no está versionada: su único respaldo es OneDrive.

Todo lo de TD está en **`TouchDesigner/PINZAS ÓPTICAS.toe`**: en `/project1` está el container **PROCESAMIENTO**.

**Dónde va cada cosa.** Los videos van a `C:\datos_pinzas\`, fuera de OneDrive y de GitHub; los análisis y los gráficos van acá, en el repo.

| Qué | Dónde |
|---|---|
| Grabaciones de ThorCam (TIFF) | `C:\datos_pinzas\GRABACIONES\<grabación>\` |
| Video para mirar (MP4, lo arma `analisis/ver_video.py`) | `C:\datos_pinzas\VISTAS\<grabación>.mp4` |
| Grabación importada (apunta a los TIFF, no los copia) | `C:\datos_pinzas\IMPORTADAS\` |
| CSV de trayectorias | `datos/<grabación>/` |
| Notebook del análisis | `analisis/<grabación>.ipynb` |
| Figuras | `figuras/<grabación>/` |
| Pasos intermedios del análisis de imagen (lo arma `TouchDesigner/herramientas/guardar_pasos.py`) | `figuras/<grabación>/pasos/` |
| Radio aparente de cada partícula (lo arma `analisis/medir_radios.py`) | `datos/<grabación>/radios.csv`, `perfiles.csv` y `figuras/<grabación>/radios/recortes.png` |

## La práctica del 12/11/2026: por dónde empezar

Se grabaron 11 videos (`C:\datos_pinzas\GRABACIONES\`, 133 GB; `datos/catalogo.csv` dice qué es cada uno) y se analizaron además dos de otros grupos. El cuaderno de laboratorio está en `docs/cuaderno_de_laboratorio.md` y pasado a tablas en `datos/cuaderno/`.

**Notebooks, en orden de lectura** (carpeta `analisis/`; cada uno termina con una conclusión y lo que no se puede afirmar):

| Notebook | Qué responde | Dudas |
|---|---|---|
| `PINZAS 1 - calibración y catálogo` | escala (94.61 nm/px con la retícula) y qué es cada grabación | `DUDAS pinzas.md` |
| `BROWNIANO comparación` | movimiento browniano en 4 mediciones, seis formas de sacar la deriva, $D$ y Stokes–Einstein, radios | `DUDAS browniano.md` |
| `PINZAS 2 - platina y patadas` | cuánto se movió de verdad la platina en cada patada, contra lo anotado en el cuaderno | `DUDAS pinzas.md` |
| `PINZAS 3 - fuerza y pozo por arrastre` | fuerza de la trampa contra posición, forma del pozo, tiempo de relajación, cuándo se escapa, y el paso a pN | `DUDAS pinzas.md` |
| `PINZAS 4 - barrido de potencia` | rigidez y frecuencia de la trampa contra la potencia, con el movimiento térmico de la partícula quieta | `DUDAS potencia.md` |

`BROWNIANO 1` y `GROS BROWNIANO` son los primeros análisis, de una grabación cada uno; `BROWNIANO comparación` los incluye y los supera.

**Tablas de cada grabación con láser** (`datos/<grabación>/`; las columnas están documentadas en el archivo que las arma):

| Archivo | Qué es | Lo arma |
|---|---|---|
| `trayectorias.csv`, `particulas.csv`, `tiempo.csv`, `parametros.json` | todas las partículas seguidas | `TouchDesigner/herramientas/convertir.py` |
| `platina.csv` | cuánto se corre toda la imagen de un cuadro al siguiente (movimiento de la platina) | `TouchDesigner/herramientas/platina.py` |
| `laser.csv` | cuánto rojo dispersa lo que hay en la trampa (≈ 0 = vacía) | `analisis/laser_csv.py` |
| `trampa_detecciones.csv` | detecciones en una ventana alrededor del láser, en todos los cuadros y sin enlazar | `analisis/trampa_roi.py` |
| `atrapada.csv` | una fila por cuadro: la partícula atrapada (posición respecto del equilibrio) y la platina | `analisis/pinzas_datos.py` |
| `seguida.csv` | la partícula atrapada seguida también cuando se aparta varios µm | `analisis/pinzas_fuerza.py` |
| `patadas.csv` | una fila por movimiento de la platina: desplazamiento, duración, velocidad, qué hizo la partícula, renglón del cuaderno | `analisis/pinzas_fuerza.py` |
| `fuerza_intervalos.csv` | una fila por par de cuadros alrededor de cada patada: fuerza de la trampa / γ y posición | `analisis/pinzas_fuerza.py` |
| `radios.csv`, `perfiles.csv` | radio aparente y perfil radial de cada partícula | `analisis/medir_radios.py` |

Los resultados de todas las grabaciones juntas están en `datos/pinzas/` (ver `LEEME.txt` para patadas y fuerza, y `LEEME potencia.txt` para el barrido de potencia) y en `datos/BROWNIANO comparación/`. El resumen de resultados y las preguntas abiertas están en `RESUMEN 12-11-2026.md`.

Tres cosas a tener en cuenta con la pinza prendida:

- **El láser ensucia la posición.** La cámara le resta rojo al verde y el láser deja un punto oscuro en el centro de la partícula atrapada, que corre la posición medida. Para el movimiento térmico de la partícula quieta hay que usar `xc_um`, `yc_um` de `datos/pinzas/potencia_cuadros.csv`, no `atrapada.csv`. Para los apartamientos grandes de las patadas el efecto es menor que el error (`DUDAS pinzas.md`, punto J).
- **Los valores del motor no sirven.** El movimiento de la platina se mide en el video (`platina.csv`).
- **Para la próxima vez:** anotar la exposición de la cámara y medir la potencia con el medidor. Son las dos incertezas mayores de la rigidez.

**Para rehacer una grabación con láser desde el video** (en este orden; `Python de TD` = `"C:/Program Files/Derivative/TouchDesigner/bin/python.exe"`, `Python global` = `C:\Dev\Global\.venv\Scripts\python.exe`):

```bash
# 1. trayectorias, sin borrar a la partícula atrapada (x, y del láser en px; ver datos/catalogo.csv)
<Python de TD> TouchDesigner/herramientas/convertir.py "C:/datos_pinzas/GRABACIONES/BARRIDO 1/ACELERACION.tif" --nombre "BARRIDO 1" --fps 10.021 --escala_nm 94.61 --set b2_fondo.Trampax=496.9 b2_fondo.Trampay=477.2 b2_fondo.Tramparadio=45 b8_trayectorias.Largomin=10
# 2. movimiento de la platina
<Python de TD> TouchDesigner/herramientas/platina.py "C:/datos_pinzas/GRABACIONES/BARRIDO 1/ACELERACION.tif" --nombre "BARRIDO 1" --fps 10.021 --trampa 496.9 477.2 60 --escala_nm 94.61
# 3. mancha roja y detecciones alrededor de la trampa (desde analisis/)
<Python global> laser_csv.py "BARRIDO 1"
<Python global> trampa_roi.py "BARRIDO 1"
# 4. las tablas atrapada / seguida / patadas / fuerza_intervalos las escribe el notebook PINZAS 2
```

Una grabación nueva hay que agregarla a `CATALOGO` en `analisis/pinzas_datos.py` (posición del láser, cuadros por segundo, tipo).

## Convertir a CSV (container PROCESAMIENTO)

Un container por función, de izquierda a derecha:

```
B0_CUADRO → B1_CANAL → B2_FONDO → B3_FOURIER → B4_CLAHE → B5_REALCE → B6_DETECCION → B7_LOCALIZACION → B8_TRAYECTORIAS
```

Cada container tiene adentro:
- **`funcion`:** su código, una función numpy autocontenida `aplicar(d, p, ctx)` en `TouchDesigner/bloques/bN_*.py`. Recibe el resultado del container de la izquierda y le pasa el suyo al de la derecha. El comentario de arriba del archivo explica qué hace y cuándo conviene cada opción.
- **`salida`:** un Script TOP que muestra **solo el resultado** de ese paso (el tile del container). La entrada es la salida del container de la izquierda.
- **Pasos intermedios** (algunos containers): Script TOPs debajo de `salida` que muestran lo que se calcula en el medio, sin calcular nada nuevo. Los declara la función en `INTERMEDIOS`.
- **Sus parámetros:** en las páginas del container. Cada uno tiene una etiqueta en castellano, y la ayuda aparece al pasar el mouse.

Los datos viajan en numpy, a resolución completa; nunca pasan por un TOP. El límite de 1280 de la Non-Commercial solo afecta lo que se ve.

**ESTADO** (el panel ancho, arriba de la red) dice qué video se está analizando (nombre, formato, cuadros, tamaño, cuadro que se muestra) y, mientras corre PROCESAR, el avance en vivo: %, cuadro x de N, segundos y cuánto falta, con una barra. Al terminar muestra cuántas partículas salieron y dónde quedaron los CSV.

### Paso a paso

1. **B0_CUADRO → Cómo está guardado:** `thorcam` (TIFF de ThorImageCAM / ThorCam: **lo recomendado**), `video` (AVI/MP4, como los de GRUPOS ANTERIORES), `secuencia` (carpeta con un TIFF por cuadro) o `grabacion` (una carpeta ya importada o sintética, con `grabacion.json`). Elegí **Carpeta** o **Archivo** según el caso. Si el archivo no trae exposición o tiempo por cuadro, completá **Exposición** y **Cuadros por segundo**.
2. Recorré los containers mirando sus imágenes y ajustá lo necesario (sobre todo **Radio** en B5 y **Umbral** en B6).
3. **B8_TRAYECTORIAS → página Procesar:** si querés, elegí **Carpeta de salida**. Después **PROCESAR TODO EL VIDEO**, y mirá ESTADO.

| Container | Qué hace | Parámetros |
|---|---|---|
| B0_CUADRO | de dónde sale el video, y el cuadro crudo (Bayer, 10 bits) | Cómo está guardado, Carpeta, Archivo, Exposición, Cuadros por segundo, Cuadro que se muestra |
| B1_CANAL | Bayer → un canal (verde/rojo/azul/suma), media resolución | Canal |
| B2_FONDO | estima la iluminación y la saca (divide o resta). **mediana/media**: integra N cuadros repartidos en todo el video; **espacial**: suavizado grande del mismo cuadro | Método, Cuadros a integrar, Escala espacial, Cómo se saca |
| B3_FOURIER | pasabanda: atenúa lo más chico que *Corte chico* (ruido de píxel) y lo más grande que *Corte grande* (iluminación). Activo, Corte chico, Corte grande |
| B4_CLAHE | contraste por zonas, **solo para detectar** (la posición se mide sobre la imagen lineal). Viene apagado: ver abajo | Activo, Límite de contraste, Teselas por lado |
| B5_REALCE | diferencia de gaussianas sobre la imagen invertida → cada bead es un pico; z = realce/ruido | Radio de la bead |
| B6_DETECCION | máximos locales con z > Umbral | Umbral, Máximo de beads |
| B7_LOCALIZACION | centro subpíxel por simetría radial (Parthasarathy 2012), en px del sensor | Lado del recorte, Suavizado |
| B8_TRAYECTORIAS | **PROCESAR** corre B0–B7 sobre todo el video (proceso aparte), enlaza y escribe los CSV | ver abajo |

Pasos intermedios, adentro de cada container:

| Container | TOP | Qué muestra |
|---|---|---|
| B2_FONDO | `fondo` | el fondo estimado que se saca, y la *desigualdad del fondo* de entrada y salida (en el sintético libre, 4.2 % → 0.2 %) |
| B3_FOURIER | `espectro` | log\|FFT\|² de la entrada, con los dos cortes como círculos |
| B3_FOURIER | `filtro` | el filtro H(k): blanco pasa, negro se atenúa |
| B3_FOURIER | `espectro_filtrado` | el espectro × H, y el *ruido de píxel* de entrada y salida (0.041 → 0.006) |
| B5_REALCE | `invertida` | la imagen aplanada e invertida, antes de la diferencia de gaussianas |

**CLAHE:** sirve cuando el contraste cambia mucho de una zona a otra (iluminación muy despareja que FONDO no sacó, beads más o menos enfocadas). Con iluminación pareja **empeora** la detección, porque estira también el ruido: en los sintéticos el z máximo baja de ~20 a ~7. Prendelo y mirá en REALCE y DETECCION si aparecen más beads.

Para cambiar un paso, editá su archivo `TouchDesigner/bloques/bN_*.py`: TD lo relee solo. PROCESAR usa esos mismos archivos, así que lo que se ve en cada container es exactamente lo que se procesa.

### Los parámetros de B8_TRAYECTORIAS

El enlace une las detecciones de cuadros consecutivos en trayectorias: a cada bead del cuadro anterior le asigna una del siguiente, eligiendo la asignación que minimiza la distancia total entre todas a la vez.

**Página Procesar**

| Parámetro | Qué es |
|---|---|
| PROCESAR TODO EL VIDEO | corre todo y escribe los CSV |
| Estado | lo escribe el proceso: avance, o la carpeta de salida al terminar |
| Carpeta de salida | vacía = `<grabación>\csv\<fecha>\`; si elegís una, se crea `<carpeta>\<grabación>_<fecha>\` |
| Escala (nm/px) | nm por px del **sensor**; 0 = la de `grabacion.json` (sale de la retícula) |
| Procesar solo los primeros | N cuadros, para probar rápido; 0 = todo |

**Página Enlace** (lo que más se toca)

| Parámetro | Qué es | Cómo elegirlo |
|---|---|---|
| Salto máximo (px) | lo más lejos que puede estar una bead de un cuadro al siguiente para considerarla la misma | ~4 × lo que se mueve por cuadro. Si es chico, las trayectorias se cortan; si es grande, se confunden beads vecinas |
| Cuadros perdidos tolerados | si una bead no se detecta en tantos cuadros seguidos, su trayectoria termina | 3 alcanza; más si parpadea o se desenfoca |
| Largo mínimo (cuadros) | las trayectorias más cortas se descartan | sacan detecciones sueltas de ruido |
| Si dos se acercan | **marcar**: se dejan y quedan con `colision = 1`; **descartar**: esos puntos quedan vacíos; **cortar**: cada una sigue con un número nuevo (nunca se intercambian ids) | ver la tabla de abajo |

**Página Avanzado**

| Parámetro | Qué es |
|---|---|
| Rellenar huecos de hasta (cuadros) | huecos cortos se rellenan interpolando, con `interpolado = 1`. 0 = nunca (lo recomendado: se interpola en el análisis si hace falta) |
| Distancia de colisión (px) | a qué distancia cuentan como "se acercan"; 0 = 2.5 radios |
| Predecir con la velocidad (cuadros) | busca la bead donde *debería* estar según su velocidad en los últimos N cuadros, no donde estaba. Imprescindible en Stokes; 0 = apagado |

### Qué cambia según lo que se grabó

| | Bead atrapada (pinza prendida) | Beads libres (pinza apagada) | Stokes (platina en movimiento) |
|---|---|---|---|
| B1 Canal | verde o azul (el láser se ve en el rojo y corre el centro) | cualquiera | verde |
| B2 Método | **espacial**: la bead casi no se mueve, así que la mediana la contiene y la borra | mediana | mediana (la bead cruza el campo) o espacial |
| B8 Salto máximo | 3 | ~4·√(2D/fps) en px: 5–8 | v·dt/escala + 3 (p. ej. 10–15) |
| B8 Si dos se acercan | marcar | **cortar** (un id intercambiado arruina el MSD) | marcar |
| B8 Largo mínimo | 1–20 | 50–100 | 20 |
| B8 Predecir con la velocidad | 10 | 10 | **10** (imprescindible) |

### Videos MP4 (los defaults están calibrados para estos)

Calibrado con `GRUPOS ANTERIORES/recorte de laser y browniano 10.mp4` (1314 × 1016, color, 15.6 fps, 828 cuadros): beads de ~35 px de diámetro, poco contraste, iluminación despareja, muchas beads pegadas al vidrio y una **deriva** de toda la muestra de ~15–20 px/s que no es igual para todas.

| Bloque | Valor | Por qué |
|---|---|---|
| B0 Cómo está guardado | `video` | elegí el `.mp4` en *Archivo* |
| B1 Canal | verde | el canal con más luz |
| B2 Método | **mediana**, 51 cuadros, dividir | borra las beads pegadas y la iluminación; las que se mueven quedan como discos negros (z ≈ 20–40) |
| B3 Corte chico / grande | 1.5 / **60** px | el corte grande tiene que ser bastante mayor que la bead (35 px): con 20 la borraba |
| B5 Radio | **17** px | radio aparente de la bead |
| B6 Umbral | **12** | las beads que se mueven dan z ≥ 15; lo que queda de las pegadas, z ≤ 10 |
| B8 Salto máximo | **10** px | paso típico 1.2 px por cuadro, 99 % < 8 px; las beads están a ~300 px entre sí |
| B8 Cuadros perdidos tolerados | 5 | |
| B8 Largo mínimo | 50 | |
| B8 Predecir con la velocidad | 10 | la deriva es lo que más mueve a las beads |

Resultado (PROCESAR en TD + `analisis/browniano.ipynb`): 719 cuadros sin repetir, ~13 beads por cuadro, 40 trayectorias (20 de más de 150 cuadros); tarda ~0.4 s por cuadro. **D = 25 ± 2 px²/s** (19 beads libres) y **1 bead atrapada** (la del láser: MSD que se aplana). Con beads de 3 µm, la escala que sale de D (79 nm/px) y la del tamaño aparente (88 nm/px) coinciden en ~10 %: el tracking mide bien el browniano.

Dos cosas de los MP4 que el importador resuelve solo:
- **Cuadros repetidos.** El video repite el cuadro anterior cuando la cámara no entregó uno nuevo (109 de 828). Se descartan al importar (un cuadro es repetido si cambia menos del 0.01 % de los píxeles) y cada cuadro nuevo conserva el tiempo en que aparece. Si no, el tiempo queda mal y el MSD tiene un pico falso en Δ = 0.
- **Búsqueda inexacta.** En MP4, saltar al cuadro *k* puede dar el *k*−1; el lector avanza leyendo en orden, así que el cuadro es siempre el pedido. Ir hacia atrás en *Cuadro que se muestra* tarda un poco más (relee desde el principio).

Y una que no: **el MP4 no guarda la exposición** (para el blur) **ni la escala**. Sin escala, todo sale en px; ver la sección 5 del notebook.

### Los CSV

**PROCESAR** crea `<grabación>\csv\<fecha>\` (o la carpeta de salida elegida) y nunca pisa una carpeta anterior. El avance se ve en ESTADO y en el parámetro **Estado** de B8.

| Archivo | Contenido |
|---|---|
| `trayectorias.csv` | una fila por partícula y cuadro: cuadro, t_s, particula, x_px, y_px, x_um, y_um, z, interpolado, colision |
| `particulas.csv` | una fila por partícula: dónde empieza y termina, cuadros medidos, interpolados, colisiones |
| `tiempo.csv` | cuadro, contador de la cámara, t_s, dt_s |
| `grabacion.json` · `parametros.json` | los metadatos y todos los parámetros de los bloques (reproducible) |
| `LEEME.txt` | qué es cada columna |

- **Posiciones:** px es el píxel del **sensor**; µm = px × escala. **Sin deriva restada**: eso es análisis (`analisis/`).
- **Tiempo:** el tiempo real de cada cuadro, no el del fps nominal. Usalo siempre en vez del número de cuadro.
- **Celda vacía:** no medido (NaN).

```python
import pandas as pd
d = pd.read_csv(r"C:\datos_pinzas\IMPORTADAS\<grabación>\csv\<fecha>\trayectorias.csv")
for pid, tr in d.groupby("particula"):
    print(pid, len(tr), tr.x_um.std())
```

Para repetir una conversión sin TD:
```bash
"C:/Program Files/Derivative/TouchDesigner/bin/python.exe" TouchDesigner/bloques/_lote.py --parametros <csv/fecha/parametros.json>
```

- **`cuadro0.png`:** el primer cuadro, para dibujar las trayectorias encima.

Para guardar lo que muestra cada container (B0 … B7) y sus pasos intermedios en algunos cuadros, con los mismos parámetros de una conversión:
```bash
"C:/Program Files/Derivative/TouchDesigner/bin/python.exe" TouchDesigner/herramientas/guardar_pasos.py "datos/<grabación>/parametros.json" "figuras/<grabación>/pasos" --cuadros 0 300 599
```
Escribe un PNG por paso a resolución completa (`cuadro_NNNN/`) y un panel con todos los pasos de cada cuadro (`pasos_cuadro_NNNN.png` y `.pdf`). Sin `--cuadros`: el primero, el del medio y el último.

## Análisis y visualización (carpeta `analisis/`)

Con los CSV, en VS Code. Necesita numpy, pandas y matplotlib en el Python que uses.

| Archivo | Qué es |
|---|---|
| `trayectorias.py` | funciones: `cargar` / `ultima_salida`, `matrices` (formato ancho, NaN = no medido), `deriva`, `sin_deriva_por_particula`, `msd`, `coeficientes_D`, `D_teorico` y los gráficos (`graficar_trayectorias`, `graficar_velocidades`, `graficar_msd`, `graficar_histograma`, …). Para comparar formas de sacar la deriva: `tiempo_regular`, `sin_deriva_media_movil`, `msd_pares`, `ajustar_msd_deriva`, `correlacion_pasos`, `simular` y `analizar_deriva` (aplica cualquiera de los seis métodos de `METODOS_DERIVA`) |
| `browniano.ipynb` | browniano libre de punta a punta: trayectorias, deriva, MSD, D por partícula, histograma de desplazamientos y chequeo de la escala |
| `BROWNIANO 1.ipynb` | la grabación BROWNIANO 1 (12/11/2026): ¿es browniano? MSD, histograma de desplazamientos contra la gaussiana y correlación entre pasos. Lee `datos/BROWNIANO 1/` y guarda en `figuras/BROWNIANO 1/` |
| `GROS BROWNIANO.ipynb` | el video de otro grupo `OTROS GRUPOS\GROS BROWNIANO` (AVI), el mismo análisis **sin restar la deriva** (`RESTAR_DERIVA = False`) |
| `BROWNIANO comparación.ipynb` | las cuatro mediciones de browniano libre (BROWNIANO 1 y 2, GROS, CHICAS): seis formas de sacar la deriva comparadas y validadas con sintéticos y con simulación gemela, chequeos de browniano, radio aparente de cada partícula, D en µm²/s y Stokes–Einstein. Lee `datos/<medición>/`, guarda en `figuras/BROWNIANO comparación/` y `datos/BROWNIANO comparación/` (`resumen.csv`, `particulas.csv`, `D_por_medicion.csv`) |
| `PINZAS 1` a `PINZAS 4` | la pinza prendida: ver «La práctica del 12/11/2026» más arriba |
| `pinzas_datos.py` | catálogo de las grabaciones con láser y tablas por cuadro: `atrapada` (la partícula en la trampa), `patadas` (movimientos de la platina), `emparejar_cuaderno` |
| `pinzas_fuerza.py` | `seguir` (la partícula atrapada también cuando se aparta), `intervalos` (fuerza de la trampa / γ entre cuadros), `ajustar_pozo` (pozo gaussiano integrando el movimiento), `resumen_patadas` |
| `potencia.py` | funciones del barrido de potencia (`PINZAS 4`): relee los TIFF alrededor del láser (rojo y posición con el láser corregido), ruido de localización por plantilla, Ornstein–Uhlenbeck sintético, estimadores de rigidez por tramo, exposición por estirado |
| `laser_csv.py`, `trampa_roi.py` | leen los TIFF directo (solo la zona del láser): mancha roja por cuadro y detecciones sin enlazar alrededor de la trampa |
| `medir_radios.py` | radio aparente de cada partícula seguida, por perfil radial sobre el video. Se corre con el Python de TD: `python.exe analisis/medir_radios.py "datos/<grabación>" [--largo_min 200]`. Para la partícula ATRAPADA no sirve tal cual (toma el punto oscuro del centro): ver `DUDAS pinzas.md`, punto I |
| `ver_video.py` | arma un MP4 para mirar una grabación de ThorCam: `python ver_video.py "C:\datos_pinzas\GRABACIONES\<grabación>"` |

El MSD y su ajuste son los de `pinzas.libre` (todos los orígenes, tiempos reales, Savin & Doyle). Además el paquete trae `pinzas.trampa` (κ por cinco métodos) y `pinzas.stokes`.

**La deriva hay que restarla antes del MSD**, y en estos videos **no es común a todas las beads**: restando solo el promedio, el MSD queda con pendiente ~1.5. Por eso el notebook le resta a cada partícula su propia recta x(t), y(t). `BROWNIANO comparación.ipynb` compara seis formas de hacerlo; la recta por partícula alcanza en las cuatro mediciones.


## Grabaciones de ThorCam

Con ThorImageCAM / ThorCam:
- Grabá en **TIFF** (multipágina u OME-TIFF). Si el programa ofrece *Image Type*, elegí **Unprocessed** (crudo, 16 bits). **Nunca AVI ni MP4**: comprimen a 8 bits y pierden el tiempo de cada cuadro.
- Anotá la exposición, el fps y la ganancia: el archivo puede no traerlos.

Para importarla:
```bash
"C:/Program Files/Derivative/TouchDesigner/bin/python.exe" TouchDesigner/herramientas/importar_thorcam.py <archivo.tif> --tipo atrapada --diametro_um 3 --exposicion_us 2000 --fps 100 --escala_nm_px 90
```
Crea una grabación nuestra que **apunta** al archivo original, sin copiarlo ni tocarlo. Esa carpeta se pone en B0 → Carpeta. **Más simple:** en B0 poné *Cómo está guardado* = `thorcam`, elegí el `.tif` en *Archivo* y completá exposición y fps; la importación se hace sola, una vez. Si el TIFF trae los tags de Thorlabs con el contador y el tiempo de cada cuadro, se usan; si no, el tiempo es **nominal** y no se pueden detectar cuadros perdidos.

Cómo configurar la cámara para que la grabación sirva (exposición, fps, formato, dónde guardar): [docs/configurar_pc.md](docs/configurar_pc.md).

## Una sola vez por computadora

1. Thorlabs: ThorImageCAM (o ThorCam, según la cámara) y Kinesis; ver [docs/configurar_pc.md](docs/configurar_pc.md).
2. Las dependencias del Python de TD (scipy, tifffile, matplotlib, pytest) están en `TouchDesigner/lib/` y sirven para 33230 y para 32460.
3. **Para rearmar la red**, por ejemplo en la Pro, corré en el Textport:
   ```python
   exec(open(project.folder + '/construir_bloques.py', encoding='utf-8').read()); construir_procesamiento(op('/project1'))
   ```

## Dónde está cada cosa

| | |
|---|---|
| Bloques de procesamiento (un archivo por container) | `TouchDesigner/bloques/b0_cuadro.py` … `b8_trayectorias.py` |
| Lo común de los containers y el lote | `TouchDesigner/bloques/_td_bloque.py`, `_lote.py` |
| Enlace de trayectorias (memoria, interpolación, colisiones) | `pinzas/enlace.py` |
| Lectura, tiempo, importación de ThorCam y videos | `pinzas/io.py`, `pinzas/tiempo.py`, `pinzas/importar.py` |
| Análisis y gráficos de los CSV (VS Code) | `analisis/trayectorias.py`, `analisis/browniano.ipynb` |
| MSD, κ, Stokes | `pinzas/libre.py`, `pinzas/trampa.py`, `pinzas/stokes.py` |
| Versión anterior, sin uso (se puede borrar) | `TouchDesigner/_anterior/` |

Tests:
```bash
"C:/Program Files/Derivative/TouchDesigner/bin/python.exe" correr_tests.py
```

Demo con sintéticos (bead atrapada, pegada y browniano libre, por los mismos bloques):
```bash
"C:/Program Files/Derivative/TouchDesigner/bin/python.exe" TouchDesigner/herramientas/demo_sinteticos.py
```
