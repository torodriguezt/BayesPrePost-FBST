# Summarize evidence across admissible full-sample reconstructions.

summarize_reconstructions <- function() {
  data <- read_vaping_counts()
  solutions <- rows <- list()

  for (item in c("daily_use_addiction", "addiction_definition")) {
    counts <- reconstruct_counts(item, data)
    message("Reconstruction: ", item, ", ", nrow(counts), " solutions")
    fits <- lapply(seq_len(nrow(counts)), function(i) {
      analyze_margins(counts[i, ], PRIORS$KL, calibrate = FALSE)
    })
    fits <- do.call(rbind, fits)
    counts$evalue <- fits$evalue
    counts$P_increase <- fits$P_increase
    solutions[[item]] <- counts
    columns <- c("n1", "n2", "p_cc", "delta_observed", "evalue", "P_increase")
    used <- counts[counts$used_in_article, columns]
    stopifnot(nrow(used) == 1L)
    rows[[item]] <- cbind(
      data.frame(
        item = item,
        summary = c("used", "minimum", "maximum"),
        solutions = nrow(counts)
      ),
      rbind(
        used, as.data.frame(lapply(counts[columns], min)),
        as.data.frame(lapply(counts[columns], max))
      )
    )
  }

  # These are ranges of evidence; no monotonicity of the adaptive cutoff is assumed.
  utils::write.csv(
    do.call(rbind, solutions),
    file.path(results_directory(), "reconstruction_solutions.csv"),
    row.names = FALSE
  )

  write_result_table(do.call(rbind, rows), "S_reconstruction")
}
