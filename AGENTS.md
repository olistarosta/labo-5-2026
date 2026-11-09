# Instrucciones para agentes

## Figuras

Todo gráfico de datos (matplotlib u otro) que se haga o se modifique en este
repositorio tiene que seguir [ESTILO_FIGURAS.md](ESTILO_FIGURAS.md). Leela antes
de escribir el código de la figura y respetá la configuración, los colores, las
marcas, los ejes, la leyenda y lo que no va en la figura tal como están ahí.

- Ejemplo de referencia: `muones/figuras/fig_barrido_angular.py`.
- Guardar cada figura en `.pdf` (la que va al informe) y `.png` a 300 dpi
  (vista rápida).
- Los resultados del ajuste (parámetros, $\chi^2$, p-valores) no van en la
  figura: van en el pie de figura o en el texto.
- Excepción: los esquemas del montaje de muones (`fig_modulo.py`,
  `fig_disposiciones.py`) usan su propio estilo (`iso.py`), no esta guía.

Si una figura necesita algo que la guía no cubre, seguí el criterio del ejemplo
de referencia y avisá qué decisión tomaste.
