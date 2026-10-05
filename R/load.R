# Load every function, in dependency order, without running any study.
# For interactive work, from the repository root: source("R/load.R").
if (!file.exists("main.R") || !dir.exists("R/Method")) {
  stop("Run from the repository root, containing main.R and R/.")
}

source("R/config.R")

# FBST numerical core: priors, e-value, predictives and adaptive cutoff (R + C++).
source("R/Method/priors.R")
source("R/Method/exact_fbst.R")
source("R/Method/design_calibration.R")

# Caches, result tables and figure files.
source("R/Utils/cache.R")
source("R/Utils/tables.R")
source("R/Utils/figures.R")

# Section 3: operating characteristics, computed exactly over the sample space.
source("R/Simulation/sampling.R")
source("R/Simulation/frequentist_tests.R")
source("R/Simulation/power.R")
source("R/Simulation/association.R")
source("R/Simulation/estimation.R")

# Section 3 and Supporting Information: sensitivity analyses.
source("R/Sensitivity/prior_sensitivity.R")
source("R/Sensitivity/boundaries.R")
source("R/Sensitivity/independent_priors.R")
source("R/Sensitivity/reference_density.R")
source("R/Sensitivity/unequal_samples.R")
source("R/Sensitivity/null_prior.R")

# Section 4: vaping knowledge application.
source("R/Data_preparation/vaping_data.R")
source("R/Real_data/fits.R")
source("R/Real_data/tables.R")
source("R/Real_data/figures.R")
source("R/Real_data/design_effect.R")
source("R/Real_data/missing_data.R")
source("R/Real_data/recoding.R")
source("R/Real_data/reconstruction.R")
