# Compare evidence after reversing the outcome coding.

compare_outcome_coding <- function() {
  data <- read_vaping_counts()
  rows <- list()

  for (name in c("KL", "informative")) {
    prior <- PRIORS[[name]]

    for (i in seq_len(nrow(data))) {
      original <- data[i, ]
      reversed <- original
      reversed$x1 <- original$n1 - original$x1
      reversed$x2 <- original$n2 - original$x2
      message("Recoding: ", original$item, ", ", original$sample, ", ", name)
      fit <- analyze_margins(original, prior)
      recoded <- analyze_margins(reversed, prior, calibrate = FALSE)

      # Hyperparameters stay fixed; P(delta' < 0) uses the original outcome scale.
      rows[[length(rows) + 1L]] <- data.frame(
        item = original$item,
        sample = original$sample,
        prior = name,
        evalue_original = fit$evalue,
        evalue_reversed = recoded$evalue,
        P_increase_original = fit$P_increase,
        P_increase_reversed = recoded$P_decrease,
        k_star = fit$k_star,
        reject_reversed = recoded$evalue <= fit$k_star
      )
    }
  }

  write_result_table(do.call(rbind, rows), "S_recoding")
}
