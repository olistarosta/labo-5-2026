# Medir con tu PC (sin la del laboratorio)

La PC del labo solo maneja dos cosas por USB: la **cámara** y la **platina motorizada**. El LED, el enfoque (micrómetro) y la preparación de la muestra son a mano. Si enchufás la cámara y el controlador de la platina a tu PC, la del labo no hace falta.

Equipo (según el punteo del grupo anterior): cámara CMOS Thorlabs de 1280 × 1024 en color, controlador **T-Cube DC Motor Controller** con actuador **Z812** (platina), LED con su driver y láser de 658 nm, 40 mW.

Revisado en tu PC (5/11/2026): Kinesis instalado, con el driver del T-Cube DC Servo y los drivers USB de Thorlabs; ThorImageCAM 1.3, con el driver de las cámaras CS / Zelux; 571 GB libres, USB 3.

## 1. Fijate qué cámara es (define el programa)

Mirá la etiqueta de la cámara:

| Modelo | Programa | Qué hacer |
|---|---|---|
| **CS135CU / CS135MU** (Zelux, 1280 × 1024) u otra **CS…** | **ThorImageCAM** | nada: ya está instalado, con su driver |
| **DCC1645C / DCC1545M / DCC…** (serie DCx, USB 2) | **ThorCam** clásico | instalar ThorCam (thorlabs.com → Software → ThorCam); trae el driver de las DCx (uc480). ThorImageCAM **no** maneja las DCx |

Prueba rápida: enchufá la cámara a tu PC y abrí ThorImageCAM. Si aparece en la lista de cámaras, listo.

## 2. Conexiones

- **Cámara → USB de tu PC**, directo, sin hub. Las Zelux van en un puerto USB 3 (azul).
- **T-Cube → USB de tu PC.** El cubo necesita **su propia fuente** de 15 V enchufada: por USB solo pasan los datos.
- **Láser:** si se maneja con perillas, no necesita PC. Si está conectado a un cubo Thorlabs (driver de diodo láser T-Cube o K-Cube), también aparece en Kinesis: Kinesis ya trae esos drivers.
- **LED y enfoque:** a mano, como siempre.

## 3. Kinesis (la platina)

1. Con el T-Cube enchufado y prendido, abrí **Kinesis**: tiene que aparecer *T-Cube DC Motor Controller* con su número de serie.
2. **Settings → actuador: Z812.** Si no, las distancias y velocidades en mm salen mal.
3. **Home** al empezar (va al fin de carrera; recorrido 12 mm). Antes, verificá que el objetivo no choque con la muestra.
4. Para la **velocidad de escape** (punteo, §3.3): *Settings → Velocity*: aceleración **0.004 mm/s²** y velocidad máxima de **0.001 a 0.020 mm/s**, de a pasos.

Para practicar en tu casa: **Kinesis Simulator** (está instalado) simula el T-Cube con el Z812. Sirve para aprender la interfaz y armar los perfiles de velocidad sin el equipo.

## 4. La cámara: cómo configurar cada grabación

**Todo fijo y en manual**; si algo se ajusta solo durante la grabación, cambia la imagen a mitad del video:

- **Exposición fija.** Corta (unos pocos ms) para que la bead no se mueva durante la exposición. **Auto exposure apagado.**
- **Ganancia fija.** **Balance de blancos automático apagado** (cámara color).
- **fps fijo**, lo más alto que dé la exposición. El MP4 de los grupos anteriores iba a 15.6 fps y con eso D se mide con pocos puntos. Con un **ROI** más chico la cámara va más rápido.
- **Guardar en TIFF** (multipágina u OME-TIFF), **no en MP4 ni AVI**, y **nunca grabar la pantalla**. El MP4 que analizamos tenía 109 cuadros repetidos de 828, tiempos inexactos y compresión con pérdida. Los TIFF de Thorlabs guardan **el tiempo real de cada cuadro** y el contador de la cámara, así que se detectan los cuadros perdidos. B0 los lee directamente (*Cómo está guardado* = `thorcam`).
- **Guardar en `C:\datos_pinzas\THORCAM\`**, fuera de OneDrive: sincronizar decenas de GB traba la PC y puede bloquear los archivos.

**Anotá en cada grabación** (el archivo puede no traerlo): exposición, fps, ganancia, ROI, potencia del láser, intensidad del LED, diámetro de las beads y la muestra (agua, % de glicerina). ThorImageCAM 1.3 tiene un campo de *notas del experimento*.

**Cuánto ocupa:** a 1280 × 1024, cada cuadro pesa de ~1.3 MB (8 bits, un canal) a ~4 MB (color). Una grabación de 5–10 min (punteo, §3.2: 5000–9000 cuadros) son **10–35 GB**. Tenés 571 GB libres.

**Al menos una vez por jornada:**
- **Retícula** (punteo, §3.1): una imagen de la regla de 10 µm con el mismo objetivo y foco. Da la escala en nm/px, que va en B8 → *Escala (nm/px)*.
- Una **bead pegada** al vidrio con la misma exposición: da el piso de ruido de la localización.

## 5. La PC durante la medición

- **Enchufada**, con el plan de energía en **alto rendimiento**.
- **Suspensión selectiva de USB desactivada**: Panel de control → Opciones de energía → configuración avanzada → USB. Si Windows apaga el puerto, se cortan la cámara o la platina.
- Que **no se suspenda ni apague la pantalla** durante la medición.
- Mientras grabás, **cerrá lo pesado** (TouchDesigner, el navegador, OneDrive pausado). Procesás después.

## 6. El flujo completo en tu PC

1. **Kinesis:** Home y posicionar la platina.
2. **ThorImageCAM / ThorCam:** en vivo, enfocar, fijar exposición, ganancia y fps, y grabar el TIFF en `C:\datos_pinzas\THORCAM\`.
3. **TouchDesigner** (`TouchDesigner/PINZAS ÓPTICAS.toe`): B0 → *Cómo está guardado* = `thorcam`, elegí el `.tif`. Si faltan exposición o fps, completalos. Después **PROCESAR TODO EL VIDEO**.
4. **VS Code:** `analisis/browniano.ipynb` sobre la salida.

**Antes del día de medición:** grabá un TIFF corto de prueba (unos segundos) y pasalo por B0. Así confirmás que el importador lee bien los tiempos de **tu** cámara y de **tu** versión del programa; ESTADO muestra el formato y los cuadros. Los tags de Thorlabs y los OME-TIFF están probados con archivos sintéticos, no con uno real de esta cámara.
