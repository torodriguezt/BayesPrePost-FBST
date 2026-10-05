# FBST adaptativo para proporciones pretest–posttest

Código R/C++ que reproduce las tablas y figuras del artículo sobre el FBST con
punto de corte adaptativo para proporciones pretest–posttest sin vínculo entre
respuestas individuales, con la aplicación a los conocimientos sobre vapeo de
McCauley et al. (2023).

## Ejecutar

Instalar R (≥ 4.4) con una herramienta de compilación C++ para Rcpp (Rtools en
Windows) y los paquetes:

```r
install.packages(c("Rcpp", "statmod"))
```

Desde la raíz del repositorio, en una sesión nueva:

```sh
Rscript main.R                  # todo, en el orden del artículo
Rscript main.R real_data        # un bloque: simulation, sensitivity o real_data
```

o, dentro de R, `source("main.R")`. Las tablas se guardan en `results/` (CSV y
TeX) y las figuras en `figures/`. El cálculo completo tarda muchas horas; los
diseños calibrados se guardan en `results/cache/` y una ejecución interrumpida
continúa donde quedó.

## Estructura

```text
main.R                  Ejecuta los tres bloques en el orden del artículo
R/
├── config.R            Rutas de entrada y salida
├── load.R              Carga todas las funciones sin ejecutar estudios
├── Data/               Conteos publicados del estudio de vapeo
├── Data_preparation/   Lectura y comprobación de los datos
├── Method/             Núcleo del FBST: priors, e-valor, predictivas y corte adaptativo (R + C++)
├── Simulation/         Sección 3: características operativas, calculadas exactamente
├── Sensitivity/        Sección 3 y material suplementario: análisis de sensibilidad
├── Real_data/          Sección 4: aplicación a los datos de vapeo
└── Utils/              Caché, escritura de tablas y figuras
results/, figures/      Salidas generadas (fuera de Git)
```

La guía del código está en [R/README.md](R/README.md).

## De cada tabla al código

Numeración de la versión JAS; entre paréntesis, la etiqueta LaTeX.

| Tabla o figura | Script | Salida en `results/` o `figures/` |
|---|---|---|
| Tabla 2 (`tab:design`) | `Simulation/power.R` | `power_exact.csv` |
| Tabla 3 (`tab:coverage`) | `Simulation/association.R` | `delta_coverage.csv`, `paired_vs_mcnemar.csv` |
| Tabla 4 (`tab:sens`), SI `tab:sens_summary`, `tab:S_sens` | `Sensitivity/prior_sensitivity.R` | `sensitivity_average.csv`, `table_S1.csv` |
| Tabla 5 (`tab:boundary`) | `Sensitivity/boundaries.R` | `boundary_summary.csv`, `boundary_null_rejection.pdf` |
| Tabla 6 (`tab:indep_decisions`), SI `tab:indep` | `Sensitivity/independent_priors.R` | `independent_prior_by_n.csv`, `independent_evalue_differences.csv` |
| Tabla 7 (`tab:vaping_data`) | `Real_data/fits.R` | `vaping_transition.csv` |
| Tabla 8 (`tab:vaping_ni`), SI `tab:S_vaping_inf`, `tab:S_vaping_conf`, `tab:vaping_freq` | `Real_data/tables.R` | `app_vaping_*.csv`, `vaping_frequentist.csv` |
| Figura 1 (`fig:vaping_errors`), Figura 2 (`fig:vaping_forest`) | `Real_data/figures.R` | `vaping_error_curves.pdf`, `vaping_delta_intervals.pdf` |
| Sección 4: efecto de diseño | `Real_data/tables.R`, `Real_data/design_effect.R` | `vaping_design_effect.csv`, `design_effect_recalibrated.csv` |
| Sección 4: encuestas no vinculadas (MAR) | `Real_data/missing_data.R` | `vaping_missingness_MAR.csv` |
| SI `tab:arcmarg` y colas del prior restringido | `Sensitivity/null_prior.R` | `arc_vs_marginal.csv`, `restricted_prior.csv` |
| SI `tab:reference` | `Sensitivity/reference_density.R` | `reference_sensitivity.csv` |
| SI `tab:unequal` | `Sensitivity/unequal_samples.R` | `unequal_samples.csv` |
| SI `tab:estimation`, `fig:est_prior_post` | `Simulation/estimation.R` | `estimation_posterior_mean.csv`, `est_prior_posterior.pdf` |
| SI `tab:S_reconstruction` | `Real_data/reconstruction.R` | `reconstruction_summary.csv`, `reconstruction_solutions.csv` |
| SI `tab:S_recoding` | `Real_data/recoding.R` | `recoding.csv` |

La Tabla 1 y las figuras del DAG y de la restricción del prior no se calculan.

## Datos

`R/Data/vaping_mccauley2023.csv` contiene conteos agregados publicados por
McCauley, Baiocchi, Cruse y Halpern-Felsher (2023), no registros individuales.
Véase [R/Data/README.md](R/Data/README.md).
