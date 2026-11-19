# Objetivo 3 — Caracterizar el potencial con el barrido de potencias

Con la partícula atrapada y la platina quieta, se baja la potencia del láser por escalones. La partícula tiembla alrededor del centro de la trampa por agitación térmica: cuanto más floja la trampa, más se aleja. De cuánto se aleja sale la rigidez del pozo a cada potencia.

## Qué hay en esta carpeta

| Archivo | Qué es |
|---|---|
| `barrido de potencias.ipynb` | el análisis completo |
| `figuras/` | sus figuras, en SVG |
| `resultados/` | tablas (`potencia_*.csv`); `LEEME.txt` explica las columnas |
| `DUDAS.md` | preguntas para responder y decisiones que se tomaron |

Las grabaciones que se usan: `LASER 5` y `BARRIDO POTENCIAS 2`. Las funciones están en `../codigo/potencia.py`.

## Cómo se hizo

1. Se mide la posición de la partícula atrapada en cada cuadro, corrigiendo el punto oscuro que deja el láser.
2. Para cada escalón de potencia se arma el **histograma de posiciones**. Donde la partícula pasa más tiempo, el pozo es más hondo: $U = -k_BT \ln N$. Cerca del centro eso da una parábola, y su curvatura es la rigidez $\kappa$.
3. A la varianza de la posición se le resta el error de medición antes de calcular la rigidez.
4. Se compara con otros caminos (autocorrelación y espectro) y con la rigidez que sale de las patadas.

Figuras para empezar: `figuras/histograma_posiciones_LASER_5_x_30mW.svg` (el histograma de un escalón), `figuras/histograma_posiciones_BARRIDO_POTENCIAS_2_x_suelta.svg` (la partícula casi suelta, sin pozo), `figuras/potencial_LASER_5_x.svg` (las parábolas) y `figuras/kappa_potencia_x.svg` (rigidez contra potencia).

## Resultados

Rigidez por el movimiento térmico de la partícula quieta, en pN/µm (exposición supuesta de medio período):

| | 30 mW | 25 | 20 | 15 | 7 |
|---|---|---|---|---|---|
| `LASER 5`, $\kappa_x$ | 2.03 ± 0.20 | 1.46 ± 0.17 | 1.17 ± 0.11 | 0.61 ± 0.07 | 0.13 ± 0.02 |
| `BARRIDO POTENCIAS 2`, $\kappa_x$ | 2.1 ± 0.7 | 1.7 ± 0.4 | 1.5 ± 0.3 | 1.15 ± 0.21 | 0.64 ± 0.08 |

- **Coincide con las patadas:** 2.0 pN/µm a 30 mW contra 1.9 ± 0.2. Son dos caminos independientes.
- **Rigidez contra potencia:** `BARRIDO POTENCIAS 2` es compatible con una recta por el origen (0.055 a 0.060 pN/µm por mW). `LASER 5` no: a 7 mW es 2.5 a 5 veces más blanda. La causa queda abierta.
- **Frecuencia de esquina:** se resuelve solo a potencia baja (1.6 Hz a 15 mW y 0.5 a 0.9 Hz a 7 mW en `LASER 5`; 1.6 Hz a 7 mW en `BARRIDO POTENCIAS 2`). A 20–30 mW queda por encima de lo que da la cámara; la extrapolación da 4.3 a 5.7 Hz a 30 mW, compatible con los 6.1 Hz de las patadas.
- **Forma del pozo cerca del centro:** parabólica hasta donde llega el movimiento térmico (3 a 5 $k_BT$). $\kappa_y/\kappa_x$ = 0.81 ± 0.06 en `LASER 5` y 0.93 ± 0.10 en `BARRIDO POTENCIAS 2`.
- **Por debajo de 7 mW** la partícula deja de estar atrapada: en `LASER 5` se la lleva el flujo; en `BARRIDO POTENCIAS 2` queda 60 s a 1–4 µm del centro.
- **La exposición** es la incerteza dominante a potencia alta: la rigidez a 30 mW va de 1.7 a 3.0 pN/µm según se suponga un período entero o cero.

## Lo que no se pudo hacer

- **Medir la exposición** de forma directa en estas grabaciones. Hay una estimación indirecta en `BARRIDO 2` (0.40 ± 0.10 del período) y una pista de que fue de al menos unos 15 ms (ver `../0 - calibración y precisión/REPORTE resolución y precisión.md`).

## Pendiente

- **¿Qué exposición tenía la cámara, y cambió entre grabaciones?** Cambia la rigidez a potencia alta hasta un factor 1.8.
- **¿Cómo se leyó la potencia?** (medidor, tabla o display). El rojo que dispersa la partícula casi no cambia entre 30 y 7 mW nominales, así que el video no puede confirmar esos escalones.
- **¿En `LASER 5` el láser se apagó a los 180 s?** El cuaderno dice «< 7 mW» y no queda rojo.
- El detalle está en `DUDAS.md`.
