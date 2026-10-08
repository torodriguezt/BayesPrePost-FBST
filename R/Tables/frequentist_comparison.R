# Compute McNemar and two-proportion z-test p-values.

compare_frequentist_tests <- function() {
  data <- read_vaping_counts()
  result <- data[c("item", "sample")]
  result$p_mcnemar <- mapply(mcnemar_pvalue, data$n01, data$n10)
  result$reject_mcnemar <- result$p_mcnemar <= .05
  result$p_z <- mapply(z_test_pvalue, data$x1, data$n1, data$x2, data$n2)
  result$reject_z <- result$p_z <= .05

  write_result_table(result, "vaping_freq")
}
