# Save figures as PDF, SVG and PNG in figures/, and build error curves.

fbst_save_figure <- function(name, draw, width = 8, height = 6) {
  directory <- fbst_figure_dir()
  paths <- character()
  for (ext in c("pdf", "svg", "png")) {
    path <- file.path(directory, paste0(name, ".", ext))
    if (ext == "pdf") grDevices::pdf(path, width = width, height = height)
    else if (ext == "svg") grDevices::svg(path, width = width, height = height)
    else grDevices::png(path, width = width, height = height, units = "in", res = 180)
    tryCatch(draw(), finally = grDevices::dev.off())
    paths <- c(paths, path)
  }
  paths
}

# Averaged errors alpha(k) and beta(k) of the rule ev <= k, as step functions of k.
fbst_error_curve <- function(design) {
  ev <- as.numeric(design$EV)
  ord <- order(ev)
  ev <- ev[ord]
  alpha <- cumsum(as.numeric(design$pH)[ord])
  beta <- 1 - cumsum(as.numeric(design$pA)[ord])
  keep <- !duplicated(ev, fromLast = TRUE)
  data.frame(k = c(-1, ev[keep]), alpha = c(0, alpha[keep]), beta = c(1, beta[keep]))
}
