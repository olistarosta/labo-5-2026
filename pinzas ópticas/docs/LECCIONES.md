# Lecciones del proyecto de fluidos → pinzas ópticas

Fuentes leídas (5 de noviembre de 2026): `fluido/TDS GLITTER/README.md`, `fluido/piv/README.md`,
`fluido/piv/como_funciona_el_analisis.md`, `piv_core.py`, `regresion_td.py`, el historial de git
y la guía interna `agents_md` de los tres componentes de `TDS GLITTER.toe` (sacada con
`toeexpand`, sin abrir ni tocar el proyecto).

## 1. Qué funcionó y se repite

| Acierto | Cómo se aplica en pinzas |
|---|---|
| Tres relojes por cuadro (driver, contador, llegada a TD) y el tiempo final de un ajuste, nunca del fps nominal | Igual, pero la Zelux da contador y timestamp **de hardware**: el ajuste pasa a ser ticks → segundos (`tiempo.py`) |
| Una carpeta por grabación con `grabacion.json`; una subcarpeta por análisis con `parametros.json`; nunca se pisa un análisis | Igual. Mismos nombres de archivo |
| Simulación con parámetros conocidos + motor de referencia en Python puro (`piv_core.py`) + prueba TD contra referencia (`regresion_td.py`, 0.003 px) | Igual, pero al revés: el paquete `pinzas/` **es** el motor, y TD solo lo llama. La prueba TD = script tiene que dar idéntico, no "parecido" |
| Lo que no se midió queda NaN (`valido = 0`). Rellenar con vecinos erraba 1.4 px contra 0.15 px de lo medido | Igual: cuadros perdidos, beads que salen del ROI y ajustes que no convergen quedan NaN y se cuentan |
| `figuras` (el tercer componente) calcula todo en numpy a resolución nativa, por etapas con clave: solo recalcula la etapa cuyos parámetros cambiaron. El código vive en un `.py` del repo que TD lee con `syncfile` ("una sola fuente") | **Este es el patrón para todo ANALISIS.toe.** Ver PLAN.md |
| Cada medición documentada con números (qué cambió, cuánto mejoró, en qué grabación) | Igual: README con tablas medidas, no adjetivos |

## 2. Lo que costó (y la regla que deja)

**El tiempo**

- La hora del driver venía cuantizada (1/64 s según `grabar`, 1/128 s según `analizar`: las dos
  notas no coinciden) y con jitter. Tomar la **mediana** de los saltos como período eligió el
  valor cuantizado (31.25 ms en vez de 33.4 ms): cada salto de 1.5 T parecía un cuadro perdido,
  esos pares salían a 2/3 de la velocidad y el resto un 7 % alto. "32 fps con 4 a 7 cuadros
  perdidos" era en realidad 30.0 fps sin ninguno.
  → **Regla:** los cuadros perdidos se detectan por el contador de hardware, nunca por umbrales
  sobre dt. Si no hay contador, el período se estima promediando, no con la mediana.
- El contador del driver a veces se leía antes de actualizarse (saltos 0 y 2 seguidos) y hubo
  que repararlo (`numerar()`).
  → **Regla:** cada lectura cruda se guarda tal cual; la reparación va aparte y se cuenta en el
  informe (`lecturas_del_contador_reparadas`).
- El simulador avanzaba con `absTime` de TD y el grabador con reloj real: TD corría a 25 fps en
  vez de 30 y ω salió **15 % bajo**.
  → **Regla:** nada usa el reloj ni el contador de cuadros de TD.

**TouchDesigner**

- `analizar` hizo el PIV *dentro* de la red: GLSL, Movie File In, un bucle con Feedback + Cache
  cocinado en `onFrameStart`. Funcionó, pero acumuló la mayoría de las trampas de la lista de
  abajo, y obligó a analizar a 1280 px por el límite de la Non-Commercial (las marcas en píxeles
  quedaban 1.5× mal sin ningún error a la vista; hubo que reescalarlas).
  → **Regla:** el cálculo nunca pasa por un TOP. TOPs y CHOPs son solo pantallas.
- Trampas medidas que no hay que volver a pisar:
  - `copyNumpyArray` con un canal: `reshape(forma + (1,))`, no `[..., None]` (stride 0).
  - Una expresión que lee un atributo suelto de una extensión **no crea dependencia**: el valor
    queda viejo sin avisar.
  - `allowCooking = False` deja los nodos vacíos: no sirve para congelar un resultado.
  - Parámetros nuevos por Python: poner `.val` además de `.default`.
  - TD **pausa al minimizar** la ventana ("Stop Playing when Minimized") y no graba nada, sin error.
  - Dibujar las vistas del proceso costaba 0.61 s por iteración: durante el análisis completo no
    se dibujan; se dibuja solo lo que se está mirando.
- El código vivía adentro del `.toe`. Resultado: 85 backups `.toe` en `Backup/`, un
  `ext_analizar_antes_del_parche.py` copiado a mano, y git no puede mostrar qué cambió.
  → **Regla:** todo el Python vive en archivos `.py` del repo; cada DAT de código los lee con
  `syncfile`. El `.toe` solo guarda la red y la disposición.

**Rendimiento**

- PNG comprimía en los hilos de TD y frenaba la grabación: se pasó a TIFF sin compresión.
- Las FFT repartidas en 8 hilos (numpy suelta el GIL) dieron lo mismo que buscaba numba, que no
  está en el Python de TD. → Antes de sumar una dependencia, medir si numpy + hilos alcanza.

**Datos**

- Los cuadros se guardaron como TIFF **RGBA de 8 bits** pasando por un TOP (`TOP.save`). Para
  fluidos alcanzaba; acá no: el sensor da 10 bits y el TOP no puede llevar Bayer crudo.
- Las grabaciones se guardaban dentro de OneDrive (`TDS GLITTER/GRABACIONES`, ~2.6 MB por cuadro).
  Hoy esas carpetas tienen `video.mp4` (fechado el 9 de octubre) en lugar de `cuadros/*.tif`:
  los TIFF se reemplazaron por video comprimido después de analizar. Además, git muestra
  borrados `tiempos.csv` y `grabacion.json` de 26 grabaciones (siguen en el historial de git).
  → **Regla:** el presupuesto de disco se calcula antes de medir. Los crudos van **fuera de
  OneDrive** (C: tiene 584 GB libres). Para archivar se usa compresión **sin pérdida** (TIFF
  con zlib/zstd), nunca mp4. En git solo van `grabacion.json`, `tiempos.csv` y los resultados.

## 2b. Lo que se aprendió armando pinzas (5 de noviembre de 2026)

- **`os.replace` en Windows falla si otro proceso tiene el archivo abierto.** Medido: 1 de cada 2 veces con un lector continuo. Eso mataba al grabador cuando TD leía `estado.json` en el mismo instante. Ahora se reintenta, y un estado que no se pudo publicar nunca frena la grabación. La caché también reintenta, y cada proceso usa su propio archivo temporal.
- **El estado no se guarda con `store()`:** queda pegado en el `.toe`. Va en atributos de la extensión, y los `Popen` en globales de módulo.
- **El trabajo largo en TD va en un proceso aparte con el Python de TD.** No se combinan hilos con subprocess (manual del MCP).
- **`-m pytest` no encuentra pytest en `lib/`,** porque se busca antes de `conftest.py`. Está `correr_tests.py`.
- **El `.gitignore` del repo ignora `*.toe`.** Por eso la red se arma con `td/construir.py` y queda en git igual.

## 2c. Para revisar (fórmulas y decisiones que no son del encargo)

**Física**
1. **DERIVA 'media' le resta a cada bead 1/nb de su propio movimiento**, así que el MSD queda multiplicado por (1 − 1/nb). LIBRE lo corrige multiplicando por nb/(nb − 1): `corregir_deriva_media = 1`. La alternativa es usar beads pegadas como referencia.
2. **Blur en la autocorrelación.** Para τ ≥ W, C(τ) = v·B(α)·e^(−τ/τc) con B(α) = 2(cosh α − 1)/α². El blur solo cambia la amplitud, así que τc sale sin sesgo.
3. **MSD atrapado con ruido y blur:** 2σ² + 2v[S(α) − B(α)e^(−τ/τc)].
4. **Espectro.** Es la forma cerrada de la Lorentziana con aliasing de BS&F por el filtro sinc² de la exposición. Coincide con la suma explícita dentro de 2e-4. Se ajusta por máxima verosimilitud de Whittle, no con los bloques y cuadrados mínimos de BS&F.
5. **Boltzmann:** ajuste de Poisson con corrección de Sheppard. Con cuadrados mínimos sobre ln(cuentas) había un sesgo de −0.5 a −1.3 %.
6. **Piso de ruido:** ½⟨Δx²⟩ de cuadros consecutivos, porque no lo infla la deriva. Existe la opción `varianza`.
7. **Pasaaltos contra deriva:** es una línea de base gaussiana en tiempo real, no un Butterworth (el Butterworth metía transitorios). Le saca ≈1.03·f_corte/f_c a la varianza y eso no se corrige.
8. **Viscosidad:** Cheng 2008 con el % de glicerina **en masa**. ¿El 20 % del encargo es en masa o en volumen?
9. **h para Faxén:** con `h_um = auto` se toma igual a `condiciones.profundidad_m`. Hay que confirmar que la profundidad anotada es la del centro de la bead.

**Tracking**

10. **Simetría radial con un suavizado previo de 0.2·R.** Es un agregado al paper: con R = 16 el original tiene 3 veces más ruido.
11. **`contraste_z` mínimo de 8 para aceptar una posición.** Hay que revisarlo con datos reales.
12. **La mancha del láser corre el centro** entre 0.04 y 1.8 px, según su forma (simulado). Con el canal verde o azul se evita.

**Tiempo**

13. **TIEMPO usa el tick ajustado contra el host, no el nominal del SDK.** Con el nominal, 14 ppm en 30 s ya dan 0.4 ms, más que el criterio.
14. **Sin timestamp de hardware, se usa el método de fluidos (recta local ±15).** Da ~1 ms. El ajuste `global` da ≤ 5 µs si el período de la cámara es fijo.

**Datos**

15. **Los sintéticos usan 1 e⁻/ADU.** Si la Zelux tiene más ganancia, el ruido real es menor que el simulado.

## 3. Entorno verificado en esta PC

| Qué | Resultado |
|---|---|
| TouchDesigner instalado | 2025.33230 (`TouchDesigner/`) y 2025.32460 (`TouchDesigner.2025.32460/`): se puede probar en la versión de la Pro sin salir de esta máquina |
| Python de TD | 3.11.15, numpy 2.1.2, OpenCV 4.11.0, pip 24.2 |
| Faltaban en TD | scipy, tifffile, matplotlib, pytest: **instalados en `lib/`** (anda en 33230 y 32460). `thorlabs_tsi_sdk`: después de instalar ThorCam, con `herramientas/instalar_sdk_thorlabs.py` |
| Python del sistema | 3.12.12: **no sirve** para instalar paquetes para TD (versión distinta) |
| ThorCam / SDK de Thorlabs | No instalado en esta PC |
| Disco | C: 952 GB, 584 GB libres; no hay segundo disco |
| MCP de TouchDesigner | twozero, en `127.0.0.1:40404/mcp`. Desde Claude se usa por HTTP |

**Instalar dependencias sin pisar el numpy de TD:**

```bash
"C:/Program Files/Derivative/TouchDesigner/bin/python.exe" -m pip install --target "<proyecto>/lib" --no-deps tifffile scipy==1.14.1
```

`--no-deps` es obligatorio: sin eso pip instala otro numpy en `lib/`, y si `lib/` queda antes en
`sys.path` lo carga en lugar del de TD. `lib/` se agrega al **final** de `sys.path`.
scipy 1.14.1 está compilado para numpy 2.x y Python 3.11 (a confirmar al instalar).

## 4. La cámara — pendiente de verificar en la PC del labo

No se puede verificar acá (ni cámara ni SDK). Lo que sé y lo que hay que medir:

| Pregunta | Lo esperado | Cómo se verifica |
|---|---|---|
| ¿La ve el Video Device In TOP? | No: la Zelux usa el driver TSI, no DirectShow ni Media Foundation | Abrir el TOP con la cámara conectada |
| Acceso por Python | `thorlabs_tsi_sdk` (viene con ThorCam) o `pylablib` `ThorlabsTLCamera` | Script de 20 líneas que abre la cámara, saca 100 cuadros e imprime los metadatos |
| Metadatos por cuadro | `frame_count` (hardware) y `time_stamp_relative_ns_or_null` en el SDK de Thorlabs; índice, framestamp y pixel clock en pylablib | Mismo script: ver si el timestamp viene `None`, y su resolución |
| Color | Bayer crudo de 10 bits; el demosaico lo hace el procesador de color del SDK en software | Mismo script: `bit_depth`, forma y dtype del buffer |
| fps máximo con ROI | ~35 fps con el sensor completo; con ROI sube sobre todo al recortar filas (lectura por filas) | Barrer el alto del ROI (1080, 512, 256, 128) con exposición mínima |
| ¿ThorCam tiene que estar cerrado? | Sí: un solo programa puede abrir la cámara | Probar con ThorCam abierto |
| ¿El hilo de TD sostiene la tasa? | El SDK y la escritura a disco sueltan el GIL; el riesgo es el disco en OneDrive | Prueba de aceptación (60 s, sensor completo y ROI chico, con ANALISIS.toe abierto) |
