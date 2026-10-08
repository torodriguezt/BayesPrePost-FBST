# Published counts; missing transition cells remain NA.
read_vaping_counts <- function(path = vaping_data_path()) {
  data <- utils::read.csv(path, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  counts <- as.matrix(data[c("n1", "n2", "x1", "x2")])

  if (anyNA(counts) || any(counts < 0) || any(counts != floor(counts)) ||
    any(data$n1 <= 0 | data$n2 <= 0 | data$x1 > data$n1 | data$x2 > data$n2)) {
    stop("Invalid marginal counts in ", path)
  }

  linked <- data$sample == "linked"
  paired <- complete.cases(data[c("n00", "n01", "n10", "n11")])
  cells <- data[paired, c("n00", "n01", "n10", "n11")]

  if (any(data$n1[linked] != data$n2[linked]) ||
    any(rowSums(cells) != data$n1[paired]) ||
    any(cells$n10 + cells$n11 != data$x1[paired]) ||
    any(cells$n01 + cells$n11 != data$x2[paired])) {
    stop("Transition cells do not reproduce the linked margins.")
  }

  data$psi <- with(data, n00 * n11 / (n01 * n10))
  data[order(
    match(data$item, unique(data$item)),
    match(data$sample, c("linked", "full"))
  ), ]
}
