# Cuaderno de análisis — pinzas ópticas

Registro de lo que se va haciendo con los datos, en orden y con fecha: qué se preguntó, qué se hizo, qué salió, qué no funcionó y qué se decidió. Es el equivalente, para el análisis, del cuaderno de laboratorio de la práctica (`docs/cuaderno_de_laboratorio.md`).

**Cómo se usa**

- Una entrada por día. Lo nuevo va **al final**.
- Cada entrada dice dónde quedó cada cosa (notebook, CSV, figura).
- Lo que no funcionó también se anota, con el motivo.
- Todo está ordenado por objetivo (ver `README.md`). Cada carpeta de objetivo tiene su `RESUMEN.md` (qué salió) y su `DUDAS.md` (qué falta responder).
- Las rutas de este cuaderno son las de hoy: cuando se reordenan carpetas se actualizan, para que sigan sirviendo.

---

## 12/11/2026 — día de la práctica y primer análisis completo

### Primera grabación: `BROWNIANO 1`

- **Pedido:** correr el análisis sobre `C:\datos_pinzas\GRABACIONES\BROWNIANO 1` y tener una forma de mirar los videos, que en TIFF no se pueden abrir.
- **Hecho:** conversión de video a CSV y notebook `1 - movimiento browniano/BROWNIANO 1.ipynb`. Visor: `codigo/ver_video.py` arma un MP4 en `C:\datos_pinzas\VISTAS\`.
- **Problema:** el importador no juntaba las 31 partes del TIFF porque el nombre de la grabación termina en número. Corregido (con test).
- **Dato:** los TIFF de ThorCam no traen el tiempo de cada cuadro ni la exposición. Los cuadros por segundo se sacan de las fechas de los archivos.
- **Resultado:** el histograma de desplazamientos es gaussiano; es browniano. Se resta la deriva (una recta por partícula).

### Orden de las carpetas

- **Decisión de Lenny:** los videos van en `C:\datos_pinzas\` (no se suben a OneDrive ni a GitHub). Análisis, CSV y gráficos van en el repo.
- La carpeta del proyecto se movió adentro del repo.

### `GROS BROWNIANO` (video de otro grupo)

- **Pedido:** analizarlo **sin restar la deriva** y guardar los pasos intermedios del análisis de imagen.
- **Hecho:** `1 - movimiento browniano/GROS BROWNIANO.ipynb`; imágenes de cada paso en `datos/GROS BROWNIANO/imágenes de control/pasos/` (las arma `TouchDesigner/herramientas/guardar_pasos.py`).

### Análisis de todo lo grabado

- **Pedido:** caracterizar el browniano con las cuatro mediciones y distintas formas de sacar la deriva; fuerza de la pinza y forma del pozo; no confiar en los valores del motor sino en las partículas no atrapadas; medir radios; usar la retícula; barrido de potencia.
- **Hecho, un notebook por tema:**

| Notebook | Resultado principal |
|---|---|
| `escala y catálogo` | escala 94.61 ± 0.04 nm/px; qué es cada grabación |
| `comparación de mediciones` | las cuatro mediciones son brownianas; restar una recta por partícula alcanza |
| `platina y patadas` | el desplazamiento real de la platina es el 5 % de lo que dice el motor |
| `fuerza y pozo` | $\kappa/\gamma$ = 38 ± 4 s⁻¹, ancho del pozo 1.14 ± 0.09 µm, rigidez 1.9 ± 0.2 pN/µm, fuerza máxima 1.3 ± 0.1 pN |
| `barrido de potencias` | rigidez térmica 2.0 pN/µm a 30 mW, baja con la potencia |

- **Herramientas nuevas:** `convertir.py` (conversión por terminal), `platina.py` (cuánto se corre la imagen), relleno del disco del láser en el fondo para no borrar a la partícula atrapada.
- **Hallazgos que cambian cómo se analiza:**
  - La receta del cuaderno (58 % de glicerina) no coincide con la difusión medida: da un arrastre 4.7 veces mayor. Se usa el arrastre medido en las partículas libres de la misma muestra, $\gamma$ = (5.0 ± 0.2) × 10⁻⁸ N·s/m.
  - La cámara le resta rojo al verde: el láser deja un punto oscuro sobre la partícula atrapada y corre la posición. Para el movimiento térmico hay que usar `xc_um`, `yc_um` de `3 - potencial con barrido de potencias/resultados/potencia_cuadros.csv`.
- **Lo que no funcionó:** estimar la exposición de la cámara por el desenfoque (dos métodos). Se borró. La exposición sigue sin conocerse.
- **Correcciones en el camino:** el emparejamiento de las patadas con el cuaderno pasó de «por orden» a «por cercanía en el tiempo» (tolerancia 4 s).
- **Dónde está todo:** el `RESUMEN.md` de cada objetivo. Dudas en `2 - potencial con patadas/DUDAS.md`, `3 - potencial con barrido de potencias/DUDAS.md`, `1 - movimiento browniano/DUDAS.md`.
- Commit `58b0795` en `main`.

## 14/11/2026 — orden del repo

- El repo quedó en una sola carpeta (`FISICA\labo-5-2026`); ya no hay carpeta anidada.

## 19/11/2026 — repaso de resultados, figuras, y resolución del video

### Repaso con Lenny (preguntas y respuestas)

- **Cómo se saca la forma del pozo.** Dos caminos. *Cerca del centro:* el histograma de posiciones de la partícula quieta; cuanto más tiempo pasa en un lugar, más hondo es el pozo ahí ($U = -k_BT \ln N$). *Lejos del centro:* las patadas de la platina sacan a la partícula y se mide con qué velocidad vuelve.
- **La rigidez sale de la parábola**: es su curvatura. Coincide con la del arrastre (2.0 contra 1.9 ± 0.2 pN/µm).
- **Cómo entra la disipación lejos del centro.** La partícula no tiene inercia: en cada instante la fuerza de la trampa iguala al roce, $F = \gamma\,(v_{partícula} - v_{líquido})$. Medir la velocidad es medir $F/\gamma$.
- **De dónde sale $\gamma$.** De la difusión de las partículas libres de la misma grabación: $\gamma = k_BT/D$ (Einstein). No usa viscosidad ni radio.

### Figuras

- Primero se juntaron todas en una carpeta; quedaban demasiadas sueltas.
- **Decisión de Lenny:** una carpeta por análisis y **solo SVG**. (Más tarde, el mismo día, pasaron a ordenarse por objetivo: ver la última parte de esta entrada.)

### Histogramas pedidos

| Qué | Figura | Qué muestra |
|---|---|---|
| Posiciones de la atrapada, `LASER 5`, $x$, 30 mW | `3 - potencial con barrido de potencias/figuras/método/histograma_posiciones_LASER_5_x_30mW.svg` | campana angosta; de acá sale la parábola |
| Pasos de las partículas libres | `2 - potencial con patadas/figuras/histograma_pasos_libres_LASER_5.svg` (y `BARRIDO_POTENCIAS_2`) | promedio 0, ancho 145 nm; más cola que una gaussiana |
| Posiciones con la partícula casi suelta, `BARRIDO POTENCIAS 2`, $x$ | `3 - potencial con barrido de potencias/figuras/método/histograma_posiciones_BARRIDO_POTENCIAS_2_x_suelta.svg` | repartidas de −3.3 a +2.7 µm, sin campana; la atrapada a 30 mW cabe en un intervalo |

- **Ojo con el histograma de pasos:** su ancho da $D$ = 0.139 µm²/s, contra 0.089 de la pendiente del MSD. La diferencia es el error de medición de cada posición (ver abajo). $\gamma$ se saca del MSD, no de ese ancho.
- **Error mío:** Lenny pidió «el histograma de movimiento libre» y primero hice el de pasos; quería el de posiciones. Quedaron los dos.
- **Tropiezo:** la celda nueva de `barrido de potencias` pisaba variables del notebook y rompía una celda posterior. Se metió en una función.

### Resolución del video y precisión del centro — `resolución y precisión`

- **Pedido:** cómo se resuelven distancias tan cortas, qué resolución tiene el video, con qué precisión se mide el centro y qué tan esféricas son las partículas.
- **Reporte:** `0 - calibración y precisión/REPORTE resolución y precisión.md`. **Notebook:** `0 - calibración y precisión/resolución y precisión.ipynb`. **Funciones:** `codigo/resolucion.py`. **CSV:** `0 - calibración y precisión/resultados/resolucion_*.csv`.

**Qué se hizo, en orden**

1. **Nitidez de la óptica con la retícula.** Perfil de una línea sobremuestreado aprovechando que está inclinada 4°. Un borde ocupa 0.5 µm.
   - *No funcionó la primera versión:* normalizaba con los costados del perfil y daba valores absurdos. Se rehízo con el ancho a media altura de la pendiente.
2. **Ruido de un píxel:** 1.6 a 1.8 % del brillo, comparando cuadros seguidos. Los píxeles vecinos se parecen (0.45 a 0.5).
3. **Partículas pegadas al vidrio de `BARRIDO 1`** (13) como patrón de «no se mueve». Error por pares: 16 a 47 nm, mediana 33.
4. **Prueba con corrimientos conocidos:** la imagen de una partícula real corrida una fracción de píxel y pegada sobre zonas vacías del video. Lineal a menos de 1 nm; error de una medición 16 nm; corrimiento fijo según el lugar 5 a 7 nm.
   - *No sirvió:* pegar la partícula sobre recortes de cuadros al azar daba 45 a 65 nm, porque mezclaba el error con el cambio de fondo de un lugar a otro. Se separó en «cambia en el tiempo» y «fijo por lugar».
5. **Cota teórica** con este ruido: 3 nm. **Comparando con la imagen de la propia partícula:** 5 nm en la prueba, 7 nm en una partícula pegada real (contra 23 con el método actual).
6. **Hipótesis descartada:** que el error viniera de que la ventana de medición salta de a un píxel. Recentrándola dos y tres veces el error queda igual.
7. **Forma:** contorno de las 13 pegadas en 36 direcciones. Elipse de 5 % en mediana, real en la imagen, no común a todas.
8. **Partículas libres** (agregado al armar el histograma de pasos): el error sale de cuánto se compensan dos pasos seguidos. De menos de 20 a 95 nm, mediana 55. Peor que las pegadas.

**Conclusiones**

- No se «ven» 20 nm: se ubica el centro de una mancha de más de mil píxeles.
- Error por cuadro: 16 a 47 nm pegadas, 19 a 39 nm la atrapada, hasta 95 nm las libres.
- Los resultados de `fuerza y pozo` y `barrido de potencias` no cambian: el error ya estaba restado en la rigidez térmica y no afecta ni al arrastre ni a la pendiente del MSD.
- Las partículas se ven redondas dentro de un 5 %; no se puede decir si ese 5 % es de la partícula.

**Pista nueva sobre la exposición.** En cuatro partículas libres bien enfocadas, dos pasos seguidos se parecen en vez de compensarse. Eso es lo que hace una exposición larga, y alcanza para decir que fue de al menos unos 15 ms (la quinta parte del tiempo entre cuadros). Es compatible con la estimación indirecta de `BARRIDO 2` (0.4 del período). No es una medición: una deriva lenta que no sea una recta da el mismo efecto, más chico.

### Todo reordenado por objetivo

- **Pedido de Lenny:** «Sigue caótico. No entiendo qué son las carpetas PINZAS 1, PINZAS 2… Ordenalo todo por objetivo»: caracterizar el movimiento browniano, caracterizar el potencial con patadas, caracterizar el potencial con barrido de potencias.
- **Hecho:** una carpeta por objetivo, y adentro de cada una lo mismo: `RESUMEN.md`, los notebooks, `figuras/`, `resultados/` y `DUDAS.md`. Lo que usan los tres (escala, catálogo, resolución) quedó en una cuarta carpeta, `0 - calibración y precisión/` (decisión mía).

| Antes | Ahora |
|---|---|
| `analisis/PINZAS 1 - calibración y catálogo.ipynb` | `0 - calibración y precisión/escala y catálogo.ipynb` |
| `analisis/PINZAS 5 - resolución y precisión.ipynb` | `0 - calibración y precisión/resolución y precisión.ipynb` |
| `analisis/BROWNIANO comparación.ipynb` | `1 - movimiento browniano/comparación de mediciones.ipynb` |
| `analisis/BROWNIANO 1.ipynb`, `GROS BROWNIANO.ipynb` | `1 - movimiento browniano/` (mismo nombre) |
| `analisis/browniano.ipynb`, `browniano_1_prueba.ipynb`, `figuras/browniano MP4/` | borrados (eran pruebas viejas; los notebooks quedan en el historial de git, los PNG no) |
| `analisis/PINZAS 2 - platina y patadas.ipynb` | `2 - potencial con patadas/platina y patadas.ipynb` |
| `analisis/PINZAS 3 - fuerza y pozo por arrastre.ipynb` | `2 - potencial con patadas/fuerza y pozo.ipynb` |
| `analisis/PINZAS 4 - barrido de potencia.ipynb` | `3 - potencial con barrido de potencias/barrido de potencias.ipynb` |
| `figuras/<análisis>/` | `<objetivo>/figuras/` |
| `datos/pinzas/` y `datos/BROWNIANO comparación/` (tablas de resultados) | `<objetivo>/resultados/` |
| `analisis/DUDAS *.md` | `<objetivo>/DUDAS.md` |
| `RESUMEN 12-11-2026.md` | repartido en el `RESUMEN.md` de cada objetivo; las preguntas abiertas, en el `README.md` |
| `reportes/2026-11-19 resolución y precisión.md` | `0 - calibración y precisión/REPORTE resolución y precisión.md` |
| `figuras/pasos y recortes/<grabación>/` | `datos/<grabación>/imágenes de control/` |
| `analisis/*.py` | `codigo/*.py` |

- **No cambió:** `datos/<grabación>/` (lo que sale de cada video) y los videos en `C:\datos_pinzas\`.
- **Código:** `codigo/carpetas.py` dice dónde escribe cada notebook (`figuras("patadas")`, `resultados("patadas")`, …). Los notebooks agregan `codigo/` al camino de Python en su primera celda.
- **Control:** se volvieron a ejecutar todos los notebooks desde su carpeta nueva.
- Se borraron las carpetas vacías `td/` y `Análisis TDS/`.

**Pendiente / para decidir**

- ¿Reprocesar la partícula atrapada comparando con su propia imagen (5–7 nm en vez de 20–35)? Mejoraría la rigidez a 25 y 30 mW. Antes hay que resolver el punto oscuro del láser.
- Siguen abiertas las preguntas del `README.md`: mezcla real de la muestra, exposición, cómo se leyó la potencia, si en `LASER 5` el láser se apagó a los 180 s.
- Lo de esta parte del día quedó confirmado en git y subido (commit «Pinzas ópticas: todo ordenado por objetivo…»).

## 19/11/2026 (tarde) — grabaciones nuevas, potencia calibrada y browniano con todo

**Pedido de Lenny:** mirar las mediciones nuevas (clase 3: `BARRIDO POTENCIAS 3` a `8`); las calibraciones de potencia anteriores «eran cualquier cosa»; los barridos 7 y 8 son los más confiables; a los barridos 1 a 6 asignarles la calibración nueva con 0.5 de error; y un buen análisis browniano con todas las grabaciones, usando las partículas no atrapadas, juntando todo y distinguiendo por radio.

### Grabaciones nuevas

- Seis videos de 240 s, 13.752 cuadros por segundo, exposición 30.5 ms (anotada en el cuaderno). Convertidos igual que los del 12/11: `datos/BARRIDO POTENCIAS 3` … `8` (trayectorias, platina, láser, detecciones en la trampa, `atrapada.csv`, radios, imágenes de control). MP4 para mirar en `C:\datos_pinzas\VISTAS\`.
- Posición del láser (mancha roja): (475.5, 484.5) px en 3 a 6 y (472.0, 479.2) px en 7 y 8. Escala: se supuso la misma, 94.61 nm/px (el pico de radio aparente de las partículas da 1.82 µm los dos días).
- El cuaderno de laboratorio nuevo quedó copiado en `docs/cuaderno_de_laboratorio.md`. Trae además los renglones de «ahora nos vamos a mover verticalmente» (22 patadas con sus tiempos), que antes no estaban: **no se pasaron todavía** a `datos/cuaderno/patadas_cuaderno.csv`.

### Calibración de la potencia

- `datos/cuaderno/calibracion_potencia.csv`: display 30 → 4.80 mW, 25 → 4.27, 20 → 3.15, 15 → 1.83, 5 → 0.0044. Entre puntos se interpola con una recta.
- `datos/cuaderno/potencias_cuaderno.csv` ahora tiene `display_mW`, `potencia_mW` (real) y `err_potencia_mW`: ± 0.5 mW en los barridos 1 a 6 (indicación de Lenny), ± 0.05 mW en 7 y 8 (**supuesto mío**).
- **Decisión:** el «7» del display (barridos 1 y 2) no tiene punto de calibración; interpolado da 0.37 mW. El «20 a 23» del barrido 2 da 3.49 mW.

### Barrido de potencias con la potencia real — `los ocho barridos con la potencia calibrada`

- **Notebook nuevo** en `3 - potencial con barrido de potencias/`. El anterior (`barrido de potencias`) queda como **el método**, con un aviso arriba de que sus mW son los del display; sus figuras pasaron a `figuras/método/` y las nuevas están en `figuras/ocho barridos/`.
- **Código:** `codigo/potencia.py` generalizado a las ocho grabaciones (`TODOS`, `NUEVOS`, `EXPOSICION_S`, `ATRAPADA_HASTA`, `escalones("display")`, `tabla_cuadros(grabaciones=…, sufijo=…)`).
- **Resultado:** el 19/11 la rigidez es proporcional a la potencia real: $\kappa_x = (0.228 \pm 0.011)$ y $\kappa_y = (0.276 \pm 0.013)$ pN/µm por mW (16 escalones, de 1.56 a 4.85 mW). Los barridos con display calibrado (4 a 6) y los de potencia medida (7 y 8) dan la misma pendiente.
- **El 12/11 da 0.40 ± 0.04** en $x$: 1.3 a 1.7 veces más rígida a igual potencia calibrada. Queda abierto (alineación, foco, o que la calibración no valga para ese día).
- **Repetición:** la misma partícula en los barridos 4, 5 y 6 difiere 18 %, contra 12 % de error calculado. Se sumó 13 % al error de todos los puntos.

**Lo que no funcionó o hubo que cambiar**

1. **El fondo de la imagen.** El 12/11 la partícula se iba de la trampa y quedaban cuadros vacíos para medir el fondo. El 19/11 tarda en irse (en el barrido 3 no se va nunca). Se cambió a: mediana de los cuadros con el láser apagado, sin los píxeles cercanos a las partículas; y el ruido de la cámara del barrido 3, de una ventana corrida 230 px.
2. **Buscar vecinas con el seguimiento no sirve.** La partícula que se le pega a la atrapada en el barrido 7 no aparece en `trayectorias.csv` ni en las detecciones (las toma como una sola mancha). Hubo que **mirar el video**: recortes cada 10 s (`figuras/ocho barridos/recortes_barrido_*.svg`).
3. **Primera pasada con todos los escalones:** los ajustes daban $\chi^2$ de 35 a 230 con 8 a 15 grados de libertad. Mirando los videos: el barrido 3 tiene motas pegadas a la partícula todo el tiempo (y no se suelta a los 150 s), el 5 se desenfoca después de los 95 s, y al 7 se le pega otra partícula a los 76 s. Sin esos escalones (quedan 26 de 36) y con el error de repetición, los $\chi^2$ quedan en 4 a 14.
4. **Tiempos de los escalones del 19/11:** el cambio cae entre 0 y 2.5 s después de la hora del cuaderno; se descartan 5 s al principio y 1 s al final de cada escalón.

### Browniano con todas las grabaciones, por radio — `todas las grabaciones por radio`

- **Notebook nuevo** en `1 - movimiento browniano/` (lo armó un agente en paralelo; revisé sus figuras y números). Funciones en `codigo/browniano_todas.py`; tablas `resultados/todas_*.csv`; figuras en `figuras/todas las grabaciones/`.
- **254 trayectorias libres** de 14 grabaciones: 141 de la muestra del 12/11, 109 del 19/11 y 4 de `BROWNIANO 1`. No aportan `LASER NE` (la imagen entera viaja), `MUEVO A MANO` ni `BARRIDO 2` (sin tramos quietos de 20 s).
- **Por radio aparente:** nominales (1.6 a 2.1 µm) $D$ = 0.076 ± 0.003 µm²/s el 12/11 y 0.091 ± 0.004 el 19/11; intermedias 0.117 y 0.180; chicas (< 1.2 µm) 0.26 los dos días.
- **Los dos días:** comparando barridos de potencia con barridos de potencia, las nominales coinciden (0.086 contra 0.091; juntas 0.089 ± 0.003). Lo que baja el promedio del 12/11 son las grabaciones con patadas (0.073 ± 0.005); no se sabe por qué.
- **Es browniano:** MSD recto (pendiente log-log 0.97 a 1.02), pasos independientes una vez descontado el error de localización, partículas que no se mueven juntas, histograma de pasos gaussiano en trozos de 6 s.
- **$D$ contra el radio:** con el radio aparente tal cual no sigue $1/R$. Sí describe la tendencia si el radio real es el aparente menos unos 0.4 µm (el anillo de difracción): una partícula de 1.5 µm se ve de 1.9 µm.
- **Decisiones del agente:** $D$ por la varianza del desplazamiento en trozos de 6 s (validado con simulaciones, 1 %); criterio de patada más estricto que el de `pinzas_datos`; clases con bordes fijos en 1.2, 1.6 y 2.1 µm; 32 pares o racimos sacados mirando recortes.

**Pendiente / para decidir**

- ¿Por qué el 12/11 la trampa era más rígida a igual potencia calibrada?
- El $\gamma$ = 5.0 ± 0.2 ×10⁻⁸ N·s/m que usa el objetivo 2 (patadas) para pasar a pN queda dentro de lo que da el browniano con todo (4.6 a 5.6), pero su error parece optimista; ± 0.4 sería más realista. **No se cambió nada en el objetivo 2.**
- Pasar al `patadas_cuaderno.csv` los renglones nuevos del cuaderno (movimiento vertical, `BARRIDO STEPSIZE 3`).
- Esta tanda quedó confirmada en git (sin subir a GitHub).
