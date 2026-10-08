# Full-sample counts compatible with the published percentages and corrected p-value.
counts_for_percentage <- function(n, percent) {
  lo <- max(0, ceiling(n * (percent - .05) / 100 - 1e-9))
  hi <- min(n, floor(n * (percent + .05) / 100 + 1e-9))

  if (hi < lo) {
    integer()
  } else {
    seq.int(lo, hi)
  }
}

reconstruct_counts <- function(item, data = read_vaping_counts()) {
  linked <- data[data$item == item & data$sample == "linked", ]
  full <- data[data$item == item & data$sample == "full", ]
  p_text <- as.character(full$p_published)
  p_published <- as.numeric(p_text)
  decimals <- nchar(sub("^[^.]*\\.", "", p_text))
  tolerance <- .5 * 10^(-decimals)
  rows <- list()

  for (n1 in seq.int(linked$n1, 600L)) {
    x1s <- counts_for_percentage(n1, full$pct_pre)
    x1s <- x1s[x1s >= linked$x1 & n1 - x1s >= linked$n1 - linked$x1]

    if (length(x1s) == 0) next

    for (n2 in seq.int(linked$n2, 410L)) {
      x2s <- counts_for_percentage(n2, full$pct_post)
      x2s <- x2s[x2s >= linked$x2 & n2 - x2s >= linked$n2 - linked$x2]

      for (x1 in x1s) {
        for (x2 in x2s) {
          p <- suppressWarnings(stats::prop.test(c(x1, x2), c(n1, n2), correct = TRUE)$p.value)

          if (abs(p - p_published) <= tolerance + 1e-12) {
            rows[[length(rows) + 1L]] <- data.frame(
              item = item,
              sample = "full",
              n1 = n1,
              x1 = x1,
              n2 = n2,
              x2 = x2,
              p_cc = p,
              delta_observed = x2 / n2 - x1 / n1
            )
          }
        }
      }
    }
  }

  if (!length(rows)) {
    stop("No admissible reconstruction for ", item)
  }

  result <- do.call(rbind, rows)
  result$used_in_article <- with(result, n1 == full$n1 & n2 == full$n2 &
    x1 == full$x1 & x2 == full$x2)
  result[order(1010 - result$n1 - result$n2, -result$n1), ]
}
