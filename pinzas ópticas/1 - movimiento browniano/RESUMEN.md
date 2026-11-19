# Objetivo 1 — Caracterizar el movimiento browniano

¿Las partículas libres se mueven como un browniano? ¿Con qué coeficiente de difusión $D$? ¿Cómo conviene sacar la deriva?

## Qué hay en esta carpeta

| Archivo | Qué es |
|---|---|
| `comparación de mediciones.ipynb` | **el análisis principal**: las cuatro mediciones juntas y seis formas de sacar la deriva |
| `todas las grabaciones por radio.ipynb` | las partículas libres de **todas** nuestras grabaciones (también las que tienen el láser prendido), separadas por radio aparente. Tablas `resultados/todas_*.csv` |
| `BROWNIANO 1.ipynb` | primer análisis, solo de nuestra grabación `BROWNIANO 1` (restando la deriva) |
| `GROS BROWNIANO.ipynb` | primer análisis del video de otro grupo, **sin** restar la deriva |
| `figuras/` | una subcarpeta por notebook, con sus figuras en SVG |
| `resultados/` | tablas del análisis principal: `D_por_medicion.csv`, `resumen.csv`, `particulas.csv`; y las de todas las grabaciones: `todas_*.csv` (el `LEEME.txt` explica las columnas) |
| `DUDAS.md` | preguntas para responder y decisiones que se tomaron |

Las trayectorias de cada grabación están en `../datos/<grabación>/trayectorias.csv`. Las funciones (MSD, deriva, histogramas) están en `../codigo/trayectorias.py`.

## Cómo se hizo

1. Cada video se convierte en trayectorias (posición de cada partícula en cada cuadro).
2. Se resta la deriva (el arrastre del líquido), de seis formas distintas para comparar.
3. Se calcula el desplazamiento cuadrático medio (MSD): cuánto se aleja una partícula, en promedio, en función del tiempo. En un browniano crece en línea recta, y la pendiente da $D$.
4. Se mira el histograma de desplazamientos: en un browniano es una campana de Gauss.

## Resultados

| Medición | Partículas solas | $D$ [µm²/s] | ¿Browniano? |
|---|---|---|---|
| BROWNIANO 1 | 3 | 0.056 ± 0.004 | sí |
| BROWNIANO 2 | 7 | 0.079 ± 0.007 | sí |
| GROS (otro grupo) | 11 | 0.169 ± 0.007 | sí |
| CHICAS (otro grupo) | 26 | 0.102 ± 0.010 | solo el MSD es creíble (video muy comprimido) |

- **Deriva:** de seis formas de sacarla, restar una recta por partícula alcanza en las cuatro mediciones y se lleva 1 a 3 % del browniano (corregido). Sin restar nada, $D$ sale entre 2 y 27 % de más.
- **Tamaños:** en BROWNIANO 2 hay partículas de 0.7 a 2.1 µm de radio aparente y $D$ baja con el tamaño.

## Con todas las grabaciones, separando por radio

`todas las grabaciones por radio.ipynb` (19/11/2026, con las 17 grabaciones). En las grabaciones con la pinza prendida también hay partículas libres en el campo. Se juntan con las de las grabaciones sin láser: **254 trayectorias libres**, 141 de la muestra del 12/11 (7 grabaciones), 109 de la del 19/11 (6 grabaciones) y 4 de `BROWNIANO 1`.

Cómo se hizo: con patadas se usan solo los tramos con la platina quieta (nunca un desplazamiento que cruce una patada); $D$ sale de la pendiente de la varianza del desplazamiento contra el retardo, restando el desplazamiento medio en trozos de 6 s para que no entre el arrastre (probado con simulaciones: devuelve el $D$ puesto dentro de 1 %); se descartan la atrapada y lo que pasa a menos de 5.7 µm del láser, las pegadas al vidrio, los cuadros con colisión, las trayectorias de menos de 200 cuadros y los pares o racimos de partículas pegadas que se ven en los recortes (32 en las grabaciones con láser).

| Muestra | Clase (radio aparente) | Trayectorias | $D$ [µm²/s] |
|---|---|---|---|
| 12/11 (la de la pinza y `BROWNIANO 2`) | chica (< 1.2 µm) | 13 | 0.258 ± 0.021 |
| | intermedia (1.2 a 1.6 µm) | 22 | 0.117 ± 0.012 |
| | **nominal (1.6 a 2.1 µm)** | 102 | **0.076 ± 0.003** |
| | grande (> 2.1 µm) | 4 | 0.090 ± 0.009 |
| 19/11 (no anotada) | chica | 36 | 0.264 ± 0.011 |
| | intermedia | 20 | 0.180 ± 0.011 |
| | **nominal** | 52 | **0.091 ± 0.004** |
| `BROWNIANO 1` | nominal | 3 | 0.055 ± 0.002 |

- **¿Los dos días dan lo mismo?** Las chicas sí (0.2 errores de diferencia). Las nominales, con todas las grabaciones, no: 20 % más el 19/11 (2.9 errores). Pero las del 19/11 son todas barridos de potencia, y comparando con los barridos de potencia del 12/11 (0.086 ± 0.004) sí coinciden (0.8 errores): juntos dan **0.089 ± 0.003 µm²/s**. Lo que baja el promedio del 12/11 son las grabaciones con patadas (0.073 ± 0.005). No hay evidencia de que la muestra del 19/11 sea otra; la diferencia que no se entiende es entre barridos de potencia y grabaciones con patadas.
- **Es browniano, los dos días:** el MSD es una recta (pendiente log-log 0.97 a 1.02), los pasos seguidos son independientes una vez descontado el error de localización, y las partículas no se mueven juntas (correlación +0.002). $D_y/D_x$ = 1.00 ± 0.02 el 12/11 y 0.93 ± 0.03 el 19/11 en las nominales (0.96 ± 0.03 con todas las libres del 19/11). El histograma de pasos tiene algo más de cola que una gaussiana (curtosis 0.35 a 0.5) porque el ancho de los pasos cambia despacio a lo largo de cada trayectoria; en trozos de 6 s es gaussiano.
- **$D$ contra el radio:** las que difunden rápido son las de radio aparente chico (el 19/11, 50 de las 51 con $D$ > 0.15 µm²/s tienen $R_{ap}$ < 1.6 µm). Con el radio aparente tal cual no sigue $1/R$ (sale $R_{ap}^{-1.5}$ a $R_{ap}^{-1.8}$). Stokes–Einstein describe la tendencia si el radio real es el aparente menos un borde de 0.39 ± 0.08 µm (el anillo de difracción): con eso una partícula de 1.5 µm se ve de 1.9 µm, que es donde está el pico del histograma. No se puede demostrar sin conocer los tamaños reales.
- **Arrastre y viscosidad** (nominales, 1.5 µm): $\gamma$ = (5.4 ± 0.2) × 10⁻⁸ N·s/m el 12/11 y (4.5 ± 0.2) × 10⁻⁸ el 19/11; en los barridos de potencia de los dos días, (4.6 ± 0.2) × 10⁻⁸. Viscosidad efectiva 1.6 a 1.9 mPa·s, la de una mezcla con 20 a 25 % de glicerina (la receta anotada el 12/11, 58 %, daría una difusión 4 a 5 veces más lenta).
- **La partícula atrapada, una vez apagado el láser** (19/11, cinco trayectorias): $D$ = 0.079 ± 0.004 µm²/s, $\gamma$ = (5.2 ± 0.3) × 10⁻⁸ N·s/m. Difunde 13 % menos que las libres nominales del mismo día.
- **Contra la hora:** el 12/11 el $D$ de las nominales baja unos 10 % por hora (la pendiente es 2.4 veces su error); el 19/11 no. Entre grabaciones el promedio cambia más que su error estadístico los dos días: por eso el error de los promedios sale de la dispersión entre grabaciones.
- **Qué no entró:** `LASER NE` (la imagen entera viaja a 0.7 µm/s), `MUEVO A MANO` y `BARRIDO 2` (no hay 20 s seguidos con la platina quieta). Los videos de otros grupos no se mezclan.

## Pendiente

- Saber qué tamaños nominales de partículas se mezclaron y qué muestra se usó el 19/11.
- Entender por qué las nominales difunden menos en las grabaciones con patadas del 12/11 que en los barridos de potencia (hora, altura del foco sobre el vidrio).
- Los agregados se marcaron mirando seis recortes por trayectoria; si se vuelve a convertir un video cambian los números de partícula y hay que mirarlos de nuevo (lista `AGREGADOS_VISTOS` del notebook).
- De los videos de otros grupos no se conoce la escala, los cuadros por segundo ni la muestra.
- El detalle está en `DUDAS.md`.
