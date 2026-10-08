# Load functions without running tables or figures.
source("R/config.R")
source("R/Method/study_settings.R")
source("R/Method/load_core.R")
source("R/Method/calibration.R")
source("R/Method/sampling_distribution.R")
source("R/Method/frequentist_tests.R")
source("R/Utils/cache.R")
source("R/Utils/write_tables.R")
source("R/Utils/plot_helpers.R")
source("R/Data_preparation/vaping_counts.R")
source("R/Data_preparation/reconstruct_counts.R")
source("R/Method/marginal_analysis.R")

for (file in list.files(
  c("R/Tables", "R/Figures"),
  pattern = "\\.R$",
  full.names = TRUE
)) {
  source(file)
}

rm(file)
