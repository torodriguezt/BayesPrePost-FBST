# Guía del código R

`main.R` llama, en el orden del artículo, a una función `fbst_run_*` por tabla o
grupo de tablas. Cada script de `Simulation/`, `Sensitivity/` y `Real_data/`
indica en su encabezado qué tabla o figura reproduce y qué archivos escribe.
Cargar un script solo define funciones; nada se calcula hasta llamarlas.

## Recorrido sugerido

1. `Method/priors.R`: priors (KL-óptimo, informativo, conflicto), escenarios y
   controles numéricos. Para cambiar un ajuste del artículo, editar aquí.
2. `Method/exact_fbst.R`: carga el núcleo numérico: núcleo posterior
   (`posterior_kernel.R`), e-valor (`evidence.R`), predictivas bajo H y A
   (`predictive.R`), corte adaptativo (`decision_rules.R`), resúmenes de δ
   (`posterior_summary.R`) y las rutinas en C++ (`exact_fbst_core.cpp`).
3. `Method/design_calibration.R`: `fbst_get_design()` calcula y guarda en caché
   los e-valores de todo el espacio muestral, las predictivas y el corte k*.
   Todos los estudios parten de esta función.
4. `Simulation/sampling.R`: distribución exacta de los márgenes, con o sin
   asociación dentro del sujeto.
5. Los estudios: `Simulation/`, `Sensitivity/` y `Real_data/`.

## Trabajar con un estudio

Desde la raíz del repositorio:

```r
source("R/load.R")
fbst_run_power()                 # Tabla 2
d <- fbst_get_design(50, 50, FBST_PRIORS$KL)
d$calibration$kstar
```

En `Real_data/`, `fbst_run_vaping()` debe ejecutarse antes que las tablas, las
figuras y los demás análisis, que leen `results/vaping_results.csv`.

## Caché

`results/cache/` guarda cada cálculo costoso con una firma que incluye sus
argumentos y los archivos de `R/Method/`. Modificar el núcleo numérico crea
cachés nuevas; reorganizar estudios, tablas o figuras no las invalida.

## Convenciones

- Un archivo por tabla o análisis del artículo; las funciones de uso común
  viven en `Method/` y `Utils/`.
- Las tablas se escriben sin redondear con `fbst_write_table()` y las figuras
  con `fbst_save_figure()`; los valores del artículo se redondean al escribirlo.
- Dos espacios de sangría y comentarios que explican decisiones estadísticas.
