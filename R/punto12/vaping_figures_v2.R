# Figures 3 and 4 (v2) for the vaping application.
#  - vaping_delta_intervals_v2: the manuscript's run_vaping_figures() with the single
#    label "conflicting prior" replaced by "conflict prior"; order, colours, values and
#    size are those of the original code, which is run unmodified otherwise.
#  - vaping_error_curves_single: one panel with alpha(k)+beta(k) of the six analyses
#    (3 items x linked/full) under the KL prior, k* marked on each curve.
# Designs come from vaping_fits_v2.R (vaping_KL_designs.rds); nothing is written to
# Figures/ or output/.
#
# Usage (from the repository root):
#   Rscript vaping_figures_v2.R <scratch_output_root> <vaping_csv> <figure_dir>
args <- commandArgs(trailingOnly = TRUE)
Sys.setenv(FBST_PROFILE = "full", FBST_OUTPUT_ROOT = args[1], FBST_VAPING_DATA = args[2])
source("R/pipeline_helpers.R"); source("R/manuscript_applications.R"); source("R/application_vaping.R")
fig_dir <- args[3]
fbst_figure_dir <- function() { dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE); fig_dir }
designs <- readRDS(file.path(fbst_output_dir(), "vaping_KL_designs.rds"))
stopifnot(length(designs) == 6L)
# run_vaping_figures() asks fbst_get_design() for the six KL designs; serve the ones
# already calibrated instead of recomputing them.
fbst_get_design <- function(n1, n2, prior, ...) {
  stopifnot(identical(unname(prior), unname(FBST_PRIORS$KL)))
  hit <- Filter(function(d) d$n1 == n1 && d$n2 == n2, designs)
  if (length(hit) != 1L) stop("No calibrated KL design for ", n1, "/", n2)
  hit[[1]]
}

# ---- Figure 3: identical code, one label changed ------------------------------------
src <- deparse(run_vaping_figures, width.cutoff = 500L)
stopifnot(sum(lengths(regmatches(src, gregexpr("\"conflicting prior\"", src, fixed = TRUE)))) == 1L)
f2 <- eval(parse(text = sub("\"conflicting prior\"", "\"conflict prior\"", src, fixed = TRUE)))
environment(f2) <- globalenv()
f2()

# ---- Figure 4 (single panel) --------------------------------------------------------
dat <- fbst_vaping_data()
res <- fbst_app_read("vaping_results"); kl <- res[res$D == 1 & res$prior_key == "KL", ]
cols <- c("#0072B2", "#D55E00", "#009E73")          # item colours of run_vaping_figures()
item_col <- setNames(cols[seq_len(nrow(dat))], dat$group)
item_labels <- c(nicotine_delivery_form = "Nicotine delivery form",
                 daily_use_addiction = "Daily use and addiction",
                 addiction_definition = "Definition of addiction")
sample_lty <- c(linked = 1, full = 2)
sample_labels <- c(linked = "linked sample", full = "full sample")
curves <- list()
for (s in c("linked", "full")) for (g in dat$group) {
  d <- designs[[paste(g, s)]]; cr <- fbst_vaping_error_curve(d); risk <- cr$alpha + cr$beta
  ks <- d$calibration$kstar
  # Stored k* must be the minimiser of the plotted curve.
  stopifnot(abs(min(risk) - d$calibration$risk) < 1e-10,
            abs(risk[max(which(cr$k <= ks))] - d$calibration$risk) < 1e-10,
            abs(ks - kl$k_star[kl$group == g & kl$sample == s]) < 1e-12)
  curves[[length(curves) + 1]] <- list(group = g, sample = s, k = cr$k, risk = risk, kstar = ks,
    rmin = d$calibration$risk, n1 = d$n1, n2 = d$n2)
}
draw_single <- function() {
  graphics::par(mar = c(4.2, 4.6, 0.8, 0.8), mgp = c(2.6, 0.7, 0), las = 1)
  graphics::plot(NA_real_, NA_real_, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
    xlab = expression("Cutoff " * k),
    ylab = expression("Prior-averaged total error " * alpha[varphi[e]](k) + beta[varphi[e]](k)))
  graphics::abline(h = seq(0, 1, .2), v = seq(0, 1, .2), col = "grey92", lwd = .8)
  graphics::box()
  for (cv in curves) graphics::lines(cv$k, cv$risk, type = "s", col = item_col[[cv$group]],
                                     lty = sample_lty[[cv$sample]], lwd = 1.6)
  for (cv in curves) graphics::points(cv$kstar, cv$rmin, pch = 21, cex = 1.25, lwd = 1.2,
                                      bg = item_col[[cv$group]], col = "white")
  lab <- vapply(curves, function(cv) sprintf("%s, %s (%d/%d): k* = %.3f", item_labels[[cv$group]],
    sample_labels[[cv$sample]], cv$n1, cv$n2, cv$kstar), "")
  graphics::legend("bottomright", inset = c(.02, .04), legend = c(lab, expression("adaptive cutoff " * k^"*")),
    col = c(vapply(curves, function(cv) item_col[[cv$group]], ""), "grey30"),
    lty = c(vapply(curves, function(cv) sample_lty[[cv$sample]], 0), NA),
    pch = c(rep(NA, length(curves)), 21), pt.bg = c(rep(NA, length(curves)), "grey30"),
    lwd = c(rep(1.6, length(curves)), 1), bg = "white", box.col = "grey80", cex = .72, seg.len = 2.6)
}
fbst_app_plot_files("vaping_error_curves_single", draw_single, width = 7, height = 5.2)
saveRDS(lapply(curves, function(cv) cv[c("group", "sample", "n1", "n2", "kstar", "rmin")]),
        file.path(fig_dir, "vaping_error_curves_single_data.rds"))
