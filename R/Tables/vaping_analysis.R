# Compute posterior changes and test equality of vaping knowledge proportions.

analyze_vaping <- function() {
  data <- read_vaping_counts()
  labels <- c(
    KL = "vaping_ni",
    informative = "vaping_inf",
    conflict = "vaping_conf"
  )
  results <- list()

  for (name in names(labels)) {
    rows <- lapply(seq_len(nrow(data)), function(i) {
      message("Vaping: ", data$item[i], ", ", data$sample[i], ", ", name)
      analyze_margins(data[i, ], PRIORS[[name]])
    })
    results[[name]] <- do.call(rbind, rows)
    results[[name]]$P_decrease <- NULL

    write_result_table(results[[name]], labels[[name]])
  }

  invisible(results)
}
