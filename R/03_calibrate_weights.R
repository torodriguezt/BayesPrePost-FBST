# Calibrates the loss weight a* from the prior alone, for each sample size
# and prior category, so that the Type I Bayes risk stays near 0.05.
# No observed data enter the calibration, so the table can be reused in any
# application. Writes output/a_star_lookup.csv.
# Runtime: about 45 min.

library(Rcpp)
sourceCpp("src/BivBetaBinom.cpp")
dir.create("output", showWarnings = FALSE)

SEED   <- 7
M_CAL  <- 3000
NGRID  <- 401
TARGET <- 0.05
A_GRID <- c(1, 1.25, 1.5, 1.75, 2, 2.5, 3, 4, 5, 7, 10, 15, 20)

n_grid <- c(30, 50, 75, 100, 150, 200, 300, 380, 421, 450)

categories <- list(
  list(label = "Non-informative", N0 = 2,
       mu0_grid = c(0.1, 0.3, 0.5, 0.7, 0.9)),
  list(label = "Informative", N0 = 50,
       mu0_grid = c(0.30, 0.40, 0.45, 0.50, 0.60)),
  list(label = "Conflict", N0 = 50,
       mu0_grid = c(0.05, 0.10, 0.15, 0.20, 0.25))
)

# Smallest weight a with alpha_phi below the target.
calibrate_a <- function(n, a0, a1, a2, seed = SEED) {
  set.seed(seed)
  evH <- simulate_evs_H(n, n, a0, a1, a2, M_CAL, NGRID)
  evA <- simulate_evs_A(n, n, a0, a1, a2, M_CAL, NGRID)
  alpha_of_a <- vapply(A_GRID, function(a)
    mean(evH <= find_kstar(evH, evA, a, 1)$k_star), 0.0)
  below <- which(alpha_of_a <= TARGET)
  if (!length(below)) return(max(A_GRID))
  if (below[1] == 1)  return(A_GRID[1])
  j <- below[1]
  approx(x = c(alpha_of_a[j - 1], alpha_of_a[j]),
         y = c(A_GRID[j - 1], A_GRID[j]), xout = TARGET)$y
}

cat("Prior-based calibration of a* (target alpha_phi = 0.05)\n\n")
rows <- list()
for (cat_info in categories) {
  cat(sprintf("== %s (N0 = %d) ==\n", cat_info$label, cat_info$N0))
  for (n in n_grid) {
    vals <- sapply(cat_info$mu0_grid, function(mu0) {
      a1 <- mu0 * cat_info$N0
      a2 <- mu0 * cat_info$N0
      a0 <- (1 - mu0) * cat_info$N0
      calibrate_a(n, a0, a1, a2)
    })
    cat(sprintf("  n = %4d  mean a* = %.3f  (range %.3f-%.3f)\n",
                n, mean(vals), min(vals), max(vals)))
    rows[[length(rows) + 1]] <- data.frame(
      category = cat_info$label, N0 = cat_info$N0, n = n,
      mu0 = cat_info$mu0_grid, a_star = vals)
  }
}

res <- do.call(rbind, rows)
write.csv(res, "output/a_star_lookup.csv", row.names = FALSE)
cat("\n-> output/a_star_lookup.csv\n\n")

# Category means, as reported in the paper.
summ <- aggregate(a_star ~ category + n, res, mean)
summ <- reshape(summ, idvar = "n", timevar = "category", direction = "wide")
names(summ) <- sub("a_star\\.", "", names(summ))
summ <- summ[order(summ$n), c("n", "Non-informative", "Informative", "Conflict")]
cat("Lookup table (category means):\n")
print(summ, row.names = FALSE, digits = 3)
