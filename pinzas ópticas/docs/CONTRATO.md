# Contrato de datos (interno, para el desarrollo)

Todo lo que comparten los módulos de `pinzas/`. Si cambiás algo acá, cambialo en el código.

## Convenciones

- Unidades SI dentro del código (m, s, N, Pa·s, K). px solo en tracking y en los bordes; la
  conversión px → m la hace DERIVA con la escala de ESCALA.
- Imágenes: numpy, fila 0 = arriba de la imagen, x = columna, y = fila (convención numpy).
  Posición subpíxel: el centro del píxel (i, j) está en (x=j, y=i).
- Lo no medido es NaN. Nada se interpola ni se inventa.
- Funciones puras: entran arrays y escalares, salen arrays, escalares o dicts de eso.
  Docstrings en español: qué entra, qué sale, unidades, de qué ecuación sale.
- Python de TD (3.11, numpy 2.1.2). `lib/` (scipy, tifffile, matplotlib, pytest) se agrega al
  FINAL de `sys.path` con `pinzas.rutas.preparar_sys_path()`.
- Tests: `"C:/Program Files/Derivative/TouchDesigner/bin/python.exe" -m pytest tests` desde
  `pinzas ópticas/` (conftest.py agrega `lib/` y la raíz). Ojo: `-m pytest` necesita encontrar
  pytest ANTES de leer conftest.py; en esta PC hace falta `PYTHONPATH=lib` (lib/ no trae numpy,
  así que no tapa el de TD).
- **px = píxel del SENSOR** (del cuadro guardado) en todo lo que sale de un paso y en
  `escala_m_por_px`. En Bayer, `imagen.preparar(..., canal='verde'|'rojo'|'azul')` da un plano de
  media resolución: las posiciones medidas ahí se pasan al sensor con
  `imagen.a_coordenadas_sensor(x, y, canal, fase, recorte)` (verde: x_s = 2x + 1/2; rojo/azul:
  x_s = 2x + dj) ANTES de multiplicar por la escala. `imagen.escala_reticula` ya devuelve m por
  px del sensor (argumento `px_sensor_por_px`).
- `bayer_fase` es la fase del cuadro GUARDADO (píxel (0,0) del archivo), no la del chip
  completo: un ROI que empieza en fila/columna impar la cambia.
- JSON: numpy → tipos nativos; NaN/inf → `null` (`io.a_nativo`). CSV: NaN → celda vacía.

## Grabación en disco

```
<DATOS>/GRABACIONES/<AAAA-MM-DD_HHMMSS>_<etiqueta>/
    grabacion.json
    tiempos.csv
    cuadros/000000.tif ...        uint16, datos crudos del sensor (Bayer si es color), 1 página
    analisis/<AAAA-MM-DD_HHMMSS>/ parametros.json, resultados.csv, figuras/   (nunca se pisa)
<DATOS>/_cache/<nombre_grabacion>/<PASO>_<clave>.npz
```

`<DATOS>` = variable de entorno `PINZAS_DATOS`, o `C:/datos_pinzas` (fuera de OneDrive).

`tiempos.csv` (una fila por cuadro GUARDADO, en orden):

| columna | tipo | qué es |
|---|---|---|
| archivo | str | `000123.tif` |
| contador_camara | int | contador de hardware; un salto > 1 = cuadros perdidos |
| t_camara | int | timestamp de hardware en ticks crudos; vacío si no hay |
| t_host_ns | int | `time.perf_counter_ns()` al recibir el cuadro |
| exposicion_us | float | exposición de ese cuadro |

Al leer (`io.leer_grabacion`): `t_camara` es SIEMPRE float64 con NaN en las vacías (exacto
para |ticks| < 2^53); las demás columnas numéricas son int64 si no tienen vacías. Bloques
multipágina (opcional): `cuadros/bloque_000000.tif`, ... escritos con
`photometric='minisblack'`; `archivo` = `bloque_000000.tif#17` (página 17) o vacío (página
global i, bloques en orden de nombre).

`grabacion.json`:

```json
{
  "version": 1,
  "nombre": "2026-11-05_143000_trampa_3um",
  "tipo": "libre | atrapada | stokes | pegada | reticula",
  "hora_inicio_iso": "2026-11-05T14:30:00.123",
  "muestra": {"material": "PS", "diametro_m": 3e-6, "glicerina_pct": 20, "preparada": "2026-11-05"},
  "condiciones": {"temperatura_K": 296.15, "laser": "...", "profundidad_m": 5e-6, "velocidad_platina_m_s": null},
  "camara": {"modelo": "CS165CU", "serie": "...", "exposicion_us": 5000, "periodo_us": 10000,
             "roi": [x, y, ancho, alto], "binning": [1, 1], "ganancia": 0, "bits": 10,
             "formato_pixel": "bayer | mono | rgb", "bayer_fase": "RGGB | GRBG | GBRG | BGGR | null",
             "tick_hz": 1e9},
  "escala_m_por_px": 9.0e-8,
  "notas": "",
  "informe_tiempo": {...},          // lo agrega el grabador al terminar (tiempo.informe)
  "verdad": {...}                   // SOLO en sintéticos: parámetros verdaderos
}
```

## Grabación importada (fuente externa: ThorCam, TIFF, AVI)

`pinzas.importar.importar(ruta, meta_usuario)` (o `herramientas/importar_thorcam.py`) crea una
carpeta de grabación NUEVA **sin `cuadros/`**: los cuadros se leen de los archivos originales,
que nunca se abren para escribir. `io.Grabacion` lo detecta por `fuente` en grabacion.json.

```json
"fuente": {"formato": "tiff | video", "programa": "ThorCam (tags TSI) | ...",
           "carpeta": "D:/labo/ThorCam",            // si movés los originales, cambiá SOLO esto
           "archivos": ["image_..._0.tif", "image_..._1.tif"], "paginas": [150, 90],
           "dtype": "uint16", "forma": [1080, 1440], "muestras": 1, "fotometrica": "MINISBLACK",
           "compresion": "NONE", "tags_privados": {"32768": 10, "32779": 1, ...},
           "formato_pixel": "bayer | mono | rgb", "cv2_api": 2200 /* solo video */},
"tiempo_fuente": "camara | ome_deltat | nominal | ninguno",
"importacion": {"desde": "...", "tags_encontrados": [...], "encontrados": [...], "faltan": [...],
                "conflictos": [...], "avisos": [...]}
```

- `tiempos.csv` de una importada: `archivo` = `nombre.tif#página` (o `video.avi#cuadro`);
  `contador_camara` = TSI_FRAME_NUMBER (vacío si no hay: NO se detectan perdidos);
  `t_camara` = TSI_TIME_STAMP_RELATIVE_NS en ns (`camara.tick_hz` = 1e9) u OME DeltaT en ns;
  `t_host_ns` vacío (no hay reloj del host: `ppm_reloj` = NaN); columna extra `t_nominal_s` =
  índice · período nominal (solo se usa si no hay ningún reloj → `tiempo.reconstruir` da
  `metodo = 'nominal'` con `informe.aviso`). `informe.contador` = 'hardware' | 'ninguno'; sin
  contador `informe.cuadros_perdidos` = null (el resultado TIEMPO dice 0 = ninguno detectado).
- `camara.formato_pixel = 'rgb'`: `g.cuadro(i)` es [alto, ancho, 3] (ya demosaicado por el
  programa). `imagen.preparar(..., 'verde')` = plano G a resolución completa,
  `px_sensor_por_px` = 1, `a_coordenadas_sensor` = identidad (+ recorte).
- `camara.bayer_fase_supuesta = true`: la fase no vino del archivo (RGGB del sensor corrida por
  el origen del ROI). No afecta al canal verde.
- Video (AVI/MP4): 8 bits, compresión con pérdida, sin timestamps: solo para mirar.

## Resultados de cada paso (dicts serializables a npz)

| Paso | Claves principales |
|---|---|
| TIEMPO | `t` [n] s desde el primer cuadro (NaN = no confiable) · `dt` [n-1] · `contador` [n] · `huecos` [k] índices i donde falta algo entre i e i+1 · `tramos` [m,2] (i0, i1 inclusive) sin huecos ni anomalías de contador · `cuadros_perdidos` · `informe` dict (ver `tiempo.reconstruir`) |
| ESCALA | `m_por_px` (por px del sensor), `err_m_por_px`, `fuente` |
| DETECCION | `candidatos` [k,3] (x, y, radio) px · `elegidas` [nb,2] (x, y) px |
| TRACKING | `x`, `y` [n, nb] px (NaN = perdido) · `metodo` |
| DERIVA | `x`, `y` [n, nb] m sin deriva · `deriva_x`, `deriva_y` [n] m |
| LIBRE | ver `libre.analizar` |
| TRAMPA | ver `trampa.analizar` |
