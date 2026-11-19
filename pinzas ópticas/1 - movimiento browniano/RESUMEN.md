# Objetivo 1 — Caracterizar el movimiento browniano

¿Las partículas libres se mueven como un browniano? ¿Con qué coeficiente de difusión $D$? ¿Cómo conviene sacar la deriva?

## Qué hay en esta carpeta

| Archivo | Qué es |
|---|---|
| `comparación de mediciones.ipynb` | **el análisis principal**: las cuatro mediciones juntas y seis formas de sacar la deriva |
| `BROWNIANO 1.ipynb` | primer análisis, solo de nuestra grabación `BROWNIANO 1` (restando la deriva) |
| `GROS BROWNIANO.ipynb` | primer análisis del video de otro grupo, **sin** restar la deriva |
| `figuras/` | una subcarpeta por notebook, con sus figuras en SVG |
| `resultados/` | tablas del análisis principal: `D_por_medicion.csv`, `resumen.csv`, `particulas.csv` |
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

## Pendiente

- De los videos de otros grupos no se conoce la escala, los cuadros por segundo ni la muestra.
- El detalle está en `DUDAS.md`.
