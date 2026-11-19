# Reporte — resolución del video y precisión del centro de la partícula

**Fecha:** 19/11/2026 · **Notebook:** `resolución y precisión.ipynb` · **Funciones:** `codigo/resolucion.py` · **Datos:** `resultados/resolucion_*.csv` · **Figuras:** `figuras/`

## Las preguntas

Los análisis de la pinza hablan de movimientos de 20 a 50 nm, y un píxel de la cámara mide 95 nm.

1. ¿Qué resolución tiene el video?
2. ¿Con qué precisión se mide el centro de una partícula, y cómo puede ser menor que un píxel?
3. ¿Qué tan esféricas son las partículas?

## Respuesta corta

- **El video** tiene píxeles de 94.6 nm y una óptica que no distingue detalles de menos de 0.5 µm.
- **El centro de una partícula se mide con un error de 16 a 35 nm por cuadro** si está pegada al vidrio (la quinta a la tercera parte de un píxel), y de **20 a 95 nm si está libre**, según qué tan en foco esté. La medición no tiene error sistemático: un corrimiento de 9.5 nm se mide como 9.5 nm.
- No es que se «vean» 20 nm. Se ubica el **centro** de una mancha de 39 px de diámetro, y en eso participan más de mil píxeles.
- **Las partículas se ven redondas dentro de un 5 %.** No se puede decir si ese 5 % es de la partícula o de cómo se la ve.
- Con otro método de localización el error podría bajar a 5–7 nm sin tocar el experimento. No se hizo.

## 1. Resolución del video

| | Valor | Cómo se midió |
|---|---|---|
| Tamaño del píxel | 94.61 ± 0.04 nm | retícula de 0.01 mm (`escala y catálogo`) |
| Imagen | 1280 × 1024 px = 121 × 97 µm, 8 bits | archivos TIFF |
| Cuadros por segundo | 10 a 13 | fechas de los archivos (`escala y catálogo`) |
| Nitidez de la óptica | 0.5 µm (5 px) | cuánto se ensancha el borde de una línea de la retícula |
| Partícula nominal | 39 px de diámetro, unos 1200 px | contorno de una partícula pegada |
| Ruido de un píxel | 1.6 a 1.8 % del brillo (1 a 2 cuentas de 255) | diferencia entre dos cuadros seguidos, en zonas vacías |
| Contraste de la partícula | 15 % del brillo | imagen promedio de una partícula |

![Una partícula en un cuadro, píxel por píxel](<figuras/particula_en_pixeles.svg>)

*Una partícula pegada al vidrio en un solo cuadro. Cada cuadradito es un píxel (95 nm).*

**Nitidez.** El borde de una línea de la retícula es un escalón abrupto. En la imagen ocupa 0.5 µm. Como las líneas están inclinadas 4°, cada fila de la foto corta el borde en una fracción distinta de píxel; juntando 400 filas se arma el perfil cada cuarto de píxel.

![Perfil de una línea de la retícula](<figuras/borde_reticula.svg>)

*Perfil de brillo a través de una línea de la retícula. Los flancos tardan 0.5 µm en bajar y subir.*

Es una cota: si la retícula estaba algo fuera de foco al sacar la foto, la óptica es mejor que eso.

**Dos cosas distintas.** *Resolver* (distinguir dos puntos cercanos) está limitado por esos 0.5 µm. *Ubicar el centro* de un objeto aislado, no: se puede ubicar el centro de una montaña con mucha más precisión que su ancho.

## 2. Precisión del centro

### Cómo se mide

La conversión de video a CSV usa **simetría radial**: en cada píxel el brillo cambia más rápido en la dirección que apunta al centro de la partícula. Se traza esa recta en cada píxel y el centro es el punto que mejor coincide con todas. El resultado no tiene por qué caer en un número entero de píxeles.

### Prueba 1: corrimientos conocidos

Se toma la imagen promedio de una partícula real, se la corre una fracción conocida de píxel, se la pega sobre zonas reales del video sin partículas (para que tenga el ruido verdadero de la cámara) y se la mide con el mismo código de la conversión.

![Corrimiento medido contra corrimiento impuesto](<figuras/corrimientos.svg>)

*Corrimiento medido contra el impuesto, de un décimo de píxel a un píxel. Promedio de 60 mediciones; la barra es el error de una sola (16 nm).*

- El promedio coincide con lo impuesto con menos de 1 nm de diferencia. **No hay salto de píxel en píxel.**
- **Una medición sola tiene 16 nm de error**, al azar de un cuadro al otro.
- El fondo de la imagen corre el centro 5 a 7 nm, distinto en cada lugar pero fijo en el tiempo. No afecta a una partícula que se queda en su sitio.

### Prueba 2: partículas que no se mueven

En `BARRIDO 1` hay 13 partículas pegadas al vidrio. Con la platina quieta su centro no cambia: lo que cambia el centro medido es error.

Para no confundir el error con una vibración de la platina se usan **pares** de partículas pegadas que están a la vez en el campo. La distancia real entre ellas es fija; lo que varía la distancia medida es el error de las dos mediciones.

![Posición medida de una partícula pegada, durante un minuto](<figuras/quieta_en_el_tiempo.svg>)

*Centro medido de una partícula pegada durante un minuto. Las líneas marcan medio píxel.*

| | Error de una medición |
|---|---|
| Mejor par | 16 nm (0.17 px) |
| Mediana de los 17 pares | 33 nm (0.35 px) |
| Peor par | 47 nm |
| Lo que la platina mueve a las dos juntas | 8 nm |

El mejor par da lo mismo que la prueba 1. O sea: lo que «se mueve» una partícula pegada es error de medición, no movimiento.

### Las partículas libres se miden peor

Las pegadas están siempre en foco. Las libres suben y bajan. Para ellas no hay posición verdadera con la que comparar, pero el error se puede sacar de los propios datos:

> Si en un cuadro la posición salió corrida a la derecha por error, el paso que llega a ese cuadro sale más largo y el que sale de él, más corto. Dos pasos seguidos tienden a compensarse. El movimiento browniano no hace eso. De cuánto se compensan sale el error.

![Error de las partículas libres contra qué tan marcadas se ven](<figuras/error_libres.svg>)

*Error de una medición de cada partícula libre de los barridos de potencia, contra qué tan marcada se ve (señal sobre ruido).*

- De menos de 20 nm a 95 nm según la partícula; mediana 55 nm (22 partículas).
- Las seis que se ven más marcadas (señal sobre ruido mayor que 33) se miden tan bien como las pegadas; las otras dieciséis dan 66 nm de mediana. La relación no es limpia: hay partículas con la misma señal y errores de 25 y de 70 nm.
- Son **mínimos**: el tiempo de exposición de la cámara, que no se conoce, tapa parte del efecto.

### Resumen de precisión

| Partícula | Error de una medición |
|---|---|
| Pegada al vidrio | 16 a 47 nm, mediana 33 nm |
| Atrapada (`barrido de potencias`) | 19 a 22 nm en `LASER 5`; 32 a 39 nm en `BARRIDO POTENCIAS 2` |
| Libre | menos de 20 a 95 nm, mediana 55 nm |

### Qué significa para los otros análisis

- **Rigidez por movimiento térmico (`barrido de potencias`).** A 30 mW la partícula atrapada se mueve unos 45 nm por agitación térmica. El error (20 a 35 nm) es comparable. Ya estaba tenido en cuenta: se lo resta de la varianza antes de calcular la rigidez. Es la razón por la que la rigidez a potencia alta tiene más error.
- **Fuerza por arrastre (`fuerza y pozo`).** Los desplazamientos son de micrones. El error no pesa.
- **Difusión y $\gamma$.** Salen de cómo crece el MSD con el tiempo. Un error de medición constante no cambia esa pendiente. No se ven afectados.
- **Histograma de pasos de las libres.** Su ancho (103 a 175 nm) mezcla la difusión con el error. No sirve para sacar $D$.
- **Promedios.** El error es al azar de cuadro a cuadro: promediando $N$ cuadros baja como $1/\sqrt N$. Por eso el centro del pozo se conoce a pocos nm.

### ¿Se podría medir mejor?

| Forma de medir | Prueba con ruido real | Partícula pegada real |
|---|---|---|
| Simetría radial (la que se usa) | 16 nm | 23 nm |
| Comparar con la imagen promedio de la propia partícula | 5 nm | 7 nm |
| Lo mejor posible con este ruido de cámara | 3 nm | |

Hay margen para medir tres veces mejor. No se cambió el método: hay que conocer cómo se ve cada partícula, y para la atrapada hay que resolver antes el punto oscuro que deja el láser.

**Lo que se descartó.** Se pensó que el error de la simetría radial venía de que la ventana donde se mide salta de a un píxel entero. No es eso: recentrando la ventana varias veces el error queda igual (16 nm).

## 3. Forma de las partículas

Se promedia la imagen de cada partícula pegada (150 cuadros) y se mide el radio del borde en 36 direcciones.

![Las ocho partículas pegadas más grandes con su contorno](<figuras/contornos_imagenes.svg>)

*Imagen promedio de ocho partículas pegadas, con el borde medido (puntos). Cada imagen mide 7.7 µm de lado.*

![Radio del borde contra dirección](<figuras/contorno_radio_direccion.svg>)

*Cuánto se aparta el radio del borde de su promedio, según la dirección, para la partícula más redonda, una intermedia y la menos redonda.*

- **Elipse:** el eje mayor es en mediana 5 % más largo que el menor (de 0.3 a 12 % según la partícula).
- **No es ruido:** con la mitad de los cuadros da lo mismo, con la misma orientación.
- **No es el método:** una partícula perfectamente redonda sobre el fondo real da 0.4 %.
- **No es la óptica:** cada partícula está alargada para un lado distinto. La parte común a todas es 1.3 ± 0.9 %. De paso, eso dice que la escala es la misma en $x$ y en $y$ dentro de 1 a 2 %.
- **Dos tamaños:** radio aparente de 1.7 a 1.85 µm (6 partículas) y de 1.1 a 1.4 µm (7).

## Lo que no se puede afirmar

- **Si el 5 % es de la partícula.** Se mide el borde de una imagen borrosa, no del material. Suciedad pegada al lado o una iluminación inclinada darían lo mismo.
- **La forma de las libres y de la atrapada.** Giran y se mueven; no se puede promediar su imagen. Las pegadas pueden no ser representativas.
- **La resolución óptica exacta.** Los 0.5 µm son una cota superior.
- **El error exacto de las libres.** Lo medido es un mínimo.
- **La precisión en las grabaciones más oscuras.** Solo se probó con una partícula de `BARRIDO STEPSIZE 3`; da parecido.

## Para decidir

1. ¿Vale la pena reprocesar la partícula atrapada comparando con su propia imagen (error de 5–7 nm en vez de 20–35)? Serviría sobre todo para la rigidez a 25 y 30 mW.
2. Si se vuelve a medir: sacar la foto de la retícula bien en foco, y anotar la exposición de la cámara.
