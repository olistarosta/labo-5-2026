# Pinzas ópticas — plan y estado

## ¿TouchDesigner sí o no?

**Sí, con una condición: TD es la pantalla con perillas de un paquete de Python, no el motor.**

Fluidos ya probó las dos formas:
- **`analizar`** calculaba adentro de la red (GLSL, Feedback, Cache) y juntó casi todas las trampas de LECCIONES.md.
- **`figuras`** calculaba en numpy y TD solo mostraba, y salió limpio.

Acá se hace todo como `figuras`. Lo que se gana con TD es ver la cadena entera, el resultado de cada paso en su propio nodo y los controles a la vista. Lo que cuesta es su modelo de cocinado; por eso todo cálculo y la grabación corren en **procesos aparte**.

## Cómo quedó

Hay un solo `.toe` (`Análisis TDS/PINZAS ÓPTICAS.toe`) con dos containers, **GRABACION** y **ANALISIS**. Se separaron en procesos, no en archivos:

```
TD  ── UDP (órdenes) ──▶  proceso grabador (pinzas/grabador.py) ── dueño de la cámara ──▶ TIFF + tiempos.csv
    ◀── estado.json + vista.bin (memmap) ──
TD  ── pedido.json ──▶  proceso de cálculo (pinzas/cadena.py) ── caché .npz ──▶ TD dibuja
```

```
pinzas ópticas/
  pinzas/        física y procesamiento, sin TD (cadena.py = el procedimiento; pasos.py = qué llama cada paso)
  tests/         101 tests con pytest sobre sintéticos
  td/            el código de cada nodo de TD (lo leen los DAT con syncfile) + construir.py (arma la red)
  herramientas/  prueba del grabador, demo con sintéticos, barrido de condiciones, instalador del SDK
  lib/           scipy, tifffile, matplotlib, pytest para el Python de TD (fuera de git)
C:\datos_pinzas\ GRABACIONES\, _cache\, SINTETICOS\   (fuera de OneDrive)
```

### Un nodo de ANALISIS

```
TIEMPO (Base COMP, tile = su resultado)
├─ codigo   td/analisis/tiempo.py: QUÉ se muestra (tabla, serie, imagen) y qué función calcula
├─ estado   listo / calculando % / error, y la clave
├─ tabla · serie · imagen   pantallas (Script DAT/CHOP/TOP): solo para mirar
└─ entrada_* / salida       el cable lleva el estado: la red dibuja las dependencias reales
ctrl_TIEMPO (zona CONTROLES)   parametro | valor | unidad | que_cambia  →  parametros.json
calc_TIEMPO                    botón CALCULAR
```

**Qué pasa al calcular un nodo:**
1. TD escribe los parámetros de todas las tablas `ctrl_*` y lanza `pinzas.cadena.calcular(<grabación>, TIEMPO)` con el Python de TD.
2. Ese proceso calcula lo que falte arriba; lo que ya está en la caché no se repite.
3. Deja el resultado en `_cache/<grabación>/TIEMPO_<clave>.npz`.

La clave es un hash de los parámetros del paso, las claves de arriba y la versión del código. Así el resultado de TD es idéntico al de `python -m pinzas.analizar`, y hay un test que lo verifica.

## Hecho

| Fase | Estado |
|---|---|
| 1. Reconocimiento (LECCIONES.md) | hecho |
| 2. Paquete, sintéticos, tests y barrido de condiciones | hecho: 101 tests y los criterios del encargo cumplidos (ver README) |
| 3. Grabador + GRABACION | hecho con la cámara **simulada**: 60 s con 0 perdidos. Falta con la Zelux real |
| 4. ANALISIS en TD | hecho, probado sobre sintéticos. Falta una grabación real |
| 5. STOKES | el nodo está y está testeado con sintéticos. Control de la platina (fase 2): no, falta tu aprobación |
| 6. README | hecho |

## Pendiente o para que decidas

1. **La cámara real.** Instalar ThorCam y correr `herramientas/instalar_sdk_thorlabs.py`. Ahí se ve qué metadatos da la Zelux (`time_stamp_relative_ns_or_null`), el fps con cada ROI y la prueba de aceptación real.
2. **Las fórmulas.** Están marcadas en LECCIONES.md → "Para revisar".
3. **Elegir beads con un clic.** Hoy se eligen escribiendo los índices que muestra DETECCION (`elegir = 0,3`). Clic sobre la imagen y centroide en vivo en GRABACION quedan para después.
4. **Los operadores de ejemplo de `/project1`** (`geo1`, `noise1`, `moviefilein1`…, del proyecto nuevo) siguen ahí. No los borré.
