# Read the published vaping summaries for linked and full samples.
# Inputs: A summary-count CSV located by fbst_vaping_path().
# Outputs: One validated analysis row per questionnaire item.
#
# Each item may have a linked row (complete pairs) and a full row (all observed
# students). Supply n1, n2, x1 and x2; full-sample counts may be recovered from
# pct_pre/pct_post when x1/x2 are empty. Optional n00/n01/n10/n11 columns provide
# paired transitions. Items without complete margins in either sample are skipped.

fbst_vaping_data <- function(path = fbst_vaping_path()) {
  if (!file.exists(path)) stop("Vaping data not found: ", path, " (see R/config.R).")
  raw <- utils::read.csv(path,
    stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = c("", "NA")
  )
  miss <- setdiff(c("sample", "item", "x1", "x2"), names(raw))
  if (length(miss)) stop("Vaping data lacks columns: ", paste(miss, collapse = ", "))
  # Sample sizes: n1/n2 per row (current layout) or N / n_pre / n_post (earlier layout).
  if (!"n1" %in% names(raw)) raw$n1 <- if ("n_pre" %in% names(raw)) raw$n_pre else NA
  if (!"n2" %in% names(raw)) raw$n2 <- if ("n_post" %in% names(raw)) raw$n_post else NA
  if ("N" %in% names(raw)) {
    raw$n1 <- ifelse(is.na(raw$n1), raw$N, raw$n1)
    raw$n2 <- ifelse(is.na(raw$n2), raw$N, raw$n2)
  }
  # Published McNemar statistic: its own column, or "chi2=..." inside test_published.
  if (!"chi2_mcnemar_published" %in% names(raw) && "test_published" %in% names(raw)) {
    hit <- regmatches(raw$test_published, regexpr("chi2 *= *[0-9.]+", raw$test_published))
    raw$chi2_mcnemar_published <- NA
    raw$chi2_mcnemar_published[grepl("chi2 *= *[0-9.]+", raw$test_published)] <- sub("chi2 *= *", "", hit)
  }
  optional <- c(
    "description", "pct_pre", "pct_post", "n00", "n01", "n10", "n11",
    "psi", "chi2_mcnemar_published"
  )
  for (v in setdiff(optional, names(raw))) raw[[v]] <- NA
  raw$description <- ifelse(is.na(raw$description), raw$item, raw$description)
  num <- function(x) suppressWarnings(as.numeric(x))
  one <- function(d, v) if (nrow(d)) num(d[[v]][1L]) else NA_real_
  rows <- lapply(unique(raw$item), function(it) {
    L <- raw[raw$sample == "linked" & raw$item == it, , drop = FALSE]
    C <- raw[raw$sample == "full" & raw$item == it, , drop = FALSE]
    if (nrow(L) > 1L || nrow(C) > 1L) stop("Duplicated item/sample rows for item ", it)
    # Linked (complete-pairs) sample.
    n_cc <- one(L, "n1")
    x1_cc <- one(L, "x1")
    x2_cc <- one(L, "x2")
    if (nrow(L) && is.finite(one(L, "n2")) && one(L, "n2") != n_cc)
      stop("Item ", it, ": the linked sample must have n1 = n2.")
    linked <- all(is.finite(c(n_cc, x1_cc, x2_cc)))
    cells <- c(n00 = one(L, "n00"), n01 = one(L, "n01"), n10 = one(L, "n10"), n11 = one(L, "n11"))
    has_cells <- linked && all(is.finite(cells))
    if (has_cells && (sum(cells) != n_cc || cells[["n10"]] + cells[["n11"]] != x1_cc ||
      cells[["n01"]] + cells[["n11"]] != x2_cc))
      stop("Transition cells of item ", it, " do not reproduce N, x1 and x2.")
    # Full sample: margins as observed, n_pre != n_post allowed.
    n1 <- one(C, "n1")
    n2 <- one(C, "n2")
    x1 <- one(C, "x1")
    x2 <- one(C, "x2")
    source_x <- "counts"
    if (all(is.finite(c(n1, n2))) && !all(is.finite(c(x1, x2)))) {
      p1 <- one(C, "pct_pre")
      p2 <- one(C, "pct_post")
      if (all(is.finite(c(p1, p2)))) {
        x1 <- round(p1 * n1 / 100)
        x2 <- round(p2 * n2 / 100)
        source_x <- "rounded_from_percent"
      }
    }
    full <- all(is.finite(c(n1, n2, x1, x2)))
    # The linked students are a subset of those seen at each stage.
    if (full && linked && (n_cc > n1 || n_cc > n2 || x1_cc > x1 || x2_cc > x2 ||
      n_cc - x1_cc > n1 - x1 || n_cc - x2_cc > n2 - x2))
      stop("Item ", it, ": the linked sample is not contained in the full sample.")
    if (!linked && !full) {
      message("Vaping item ", it, ": no complete margins in either sample; skipped.")
      return(NULL)
    }
    if (!full) {
      n1 <- n2 <- x1 <- x2 <- NA_real_
      source_x <- NA_character_
    }
    if (!has_cells) cells[] <- NA_real_
    data.frame(
      study = "vaping", group = it,
      description = if (nrow(L)) L$description[1L] else C$description[1L],
      n1 = n1, x1 = x1, n2 = n2, x2 = x2, full_available = full, full_margins_source = source_x,
      n_cc = n_cc, x1_cc = x1_cc, x2_cc = x2_cc, linked_available = linked,
      n00 = cells[["n00"]], n01 = cells[["n01"]], n10 = cells[["n10"]], n11 = cells[["n11"]],
      cells_available = has_cells,
      baseline_only = if (full && linked) n1 - n_cc else NA_real_,
      final_only = if (full && linked) n2 - n_cc else NA_real_,
      baseline_success_lost = if (full && linked) x1 - x1_cc else NA_real_,
      psi = if (has_cells && cells[["n01"]] * cells[["n10"]] > 0)
        cells[["n00"]] * cells[["n11"]] / (cells[["n01"]] * cells[["n10"]]) else NA_real_,
      psi_published = one(L, "psi"),
      chi2_cc = if (has_cells && cells[["n01"]] + cells[["n10"]] > 0)
        (abs(cells[["n01"]] - cells[["n10"]]) - 1)^2 / (cells[["n01"]] + cells[["n10"]]) else NA_real_,
      chi2_published = one(L, "chi2_mcnemar_published"),
      pct_pre_linked = one(L, "pct_pre"), pct_post_linked = one(L, "pct_post"),
      pct_pre_full = one(C, "pct_pre"), pct_post_full = one(C, "pct_post"),
      data_source = basename(path), stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || !nrow(out)) stop("No analysable item in ", path)
  out
}
