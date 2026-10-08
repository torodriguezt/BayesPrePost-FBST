# Export the same drawing to PDF, SVG and PNG.
save_figure <- function(name, draw, width = 8, height = 6) {
  directory <- figures_directory()

  for (extension in c("pdf", "svg", "png")) {
    path <- file.path(directory, paste0(name, ".", extension))

    if (extension == "pdf") {
      grDevices::pdf(path, width = width, height = height)
    } else if (extension == "svg") {
      grDevices::svg(path, width = width, height = height)
    } else {
      grDevices::png(path, width = width, height = height, units = "in", res = 180)
    }

    tryCatch(draw(), finally = grDevices::dev.off())
  }

  invisible(NULL)
}

# Averaged errors of the decision rule e-value <= k.
error_curve <- function(design) {
  order <- order(as.numeric(design$EV))
  evidence <- as.numeric(design$EV)[order]
  alpha <- cumsum(as.numeric(design$pH)[order])
  beta <- sum(design$pA) - cumsum(as.numeric(design$pA)[order])
  keep <- !duplicated(evidence, fromLast = TRUE)

  data.frame(
    k = c(-1, evidence[keep]),
    alpha = c(0, alpha[keep]),
    beta = c(sum(design$pA), beta[keep])
  )
}

ITEM_LABELS <- c(
  nicotine_delivery_form = "Nicotine delivery form",
  daily_use_addiction = "Daily use and addiction",
  addiction_definition = "Definition of addiction"
)
SAMPLE_LABELS <- c(linked = "Linked students", full = "Margins as observed")
PRIOR_LABELS <- c(
  KL = "KL-optimal",
  informative = "Informative",
  conflict = "Conflict"
)
PRIOR_COLORS <- c(
  KL = "#0072B2",
  informative = "#D55E00",
  conflict = "#009E73"
)
