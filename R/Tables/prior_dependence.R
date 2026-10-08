# Compare evidence and decisions under dependent and independent priors.

compare_prior_dependence <- function() {
  theta <- SCENARIOS$baseline
  decisions <- differences <- list()

  for (name in c("KL", "informative")) {
    prior <- PRIORS[[name]]
    sizes <- if (name == "KL") {
      sort(unique(c(50L, SCENARIOS$independent_n)))
    } else {
      SCENARIOS$independent_n
    }

    for (n in sizes) {
      message("Independent priors: ", name, ", n = ", n)
      dependent_design <- calibrate_test(n, n, prior)
      independent_design <- calibrate_test(n, n, prior, family = "independent")
      difference <- abs(dependent_design$EV - independent_design$EV)

      if (n %in% SCENARIOS$independent_n) {
        mass <- sampling_distribution(n, n, theta, theta + .20)
        differences[[length(differences) + 1L]] <- data.frame(
          prior = name,
          lambda = sum(prior),
          n = n,
          mean_abs_evalue_difference = sum(mass * difference)
        )
      }

      if (name == "KL" && n %in% c(50, 100, 400)) {
        reject_dependent <- dependent_design$EV <= dependent_design$calibration$kstar
        reject_independent <- independent_design$EV <= independent_design$calibration$kstar
        null <- sampling_distribution(n, n, theta, theta)
        alternative <- sampling_distribution(n, n, theta, theta + .10)
        decisions[[length(decisions) + 1L]] <- data.frame(
          n = n,
          k_OL = dependent_design$calibration$kstar,
          k_ind = independent_design$calibration$kstar,
          alpha_OL = dependent_design$calibration$alpha,
          alpha_ind = independent_design$calibration$alpha,
          beta_OL = dependent_design$calibration$beta,
          beta_ind = independent_design$calibration$beta,
          type_I_OL = sum(null * reject_dependent),
          type_I_ind = sum(null * reject_independent),
          power_OL = sum(alternative * reject_dependent),
          power_ind = sum(alternative * reject_independent),
          disagreement_0 = sum(null * (reject_dependent != reject_independent)),
          disagreement_010 = sum(alternative * (reject_dependent != reject_independent)),
          mean_abs_evalue_difference = sum(alternative * difference)
        )
      }
    }
  }

  write_result_table(do.call(rbind, decisions), "indep_decisions")

  write_result_table(do.call(rbind, differences), "indep")

  invisible(decisions)
}
