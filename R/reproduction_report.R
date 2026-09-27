# Traceability is generated from the current manuscript and actual output files.
source("R/pipeline_helpers.R")
fbst_manuscript_audit <- function(path) {
  text <- readLines(path,warn=FALSE)
  labels <- grep("\\\\label\\{",text)
  asset_lines <- grep("\\\\includegraphics",text)
  pending <- grep("PENDIENTE-|verificar que el codigo",text)
  mapping <- c("tab:design"="power_vs_n","tab:coverage"="delta_coverage",
    "tab:assoc"="paired_vs_mcnemar","tab:sens_summary"="prior_sensitivity",
    "tab:sens_alpha"="prior_sensitivity","tab:indep_decisions"="independent_prior_by_n",
    "tab:arcmarg"="arc_vs_marginal","tab:indep"="independent_evalue_differences",
    "tab:thks_transition"="tvsfp_transition","tab:thks_ni"="tvsfp_results",
    "tab:thks_inf"="tvsfp_results","tab:thks_conf"="tvsfp_results",
    "tab:thks_mcnemar"="tvsfp_results","tab:between"="tvsfp_contrasts",
    "tab:toenail_ni"="toenail_results","tab:toenail_freq"="toenail_results",
    "tab:toenail_inf"="toenail_results","prior1"="estimation_exact_sampling",
    "prior2"="estimation_exact_sampling","prior3"="estimation_exact_sampling",
    "sim1"="estimation_exact_sampling","sim2"="estimation_exact_sampling","sim3"="estimation_exact_sampling")
  inv <- lapply(labels,function(i) {
    label <- sub(".*\\\\label\\{([^}]+)\\}.*","\\1",text[i])
    out <- if(label %in% names(mapping))unname(mapping[label]) else ""
    file <- if(nzchar(out))file.path(fbst_output_dir(),paste0(out,".csv")) else ""
    data.frame(line=i,label=label,source=text[i],result=file,
      status=if(nzchar(file)&&file.exists(file))if(fbst_profile()=="pilot")"pilot_output_only" else "output_available_review_required" else "pending_or_not_computational",
      stringsAsFactors=FALSE)
  })
  fbst_write_table(do.call(rbind,inv),"manuscript_inventory")
  # Preserve every line containing numeric content as a searchable audit ledger.
  numeric <- grep("[0-9]",text)
  numeric <- numeric[!grepl("^\\s*%",text[numeric])]
  fbst_write_table(data.frame(line=numeric,source=text[numeric],verification="requires_result_link"),
                   "manuscript_numeric_claims")
  if(length(pending)) fbst_write_table(data.frame(line=pending,comment=text[pending]),"manuscript_pending_comments")
  if(length(asset_lines)) fbst_write_table(data.frame(line=asset_lines,source=text[asset_lines]),"manuscript_figures")
}
fbst_reproduction_report <- function() {
  manuscript <- fbst_manuscript_path()
  if(file.exists(manuscript)) fbst_manuscript_audit(manuscript)
  check_path <- file.path(fbst_output_dir(),"validation_checks.csv")
  if(file.exists(check_path)) {
    checks <- read.csv(check_path,stringsAsFactors=FALSE)
    d <- checks[checks$kind=="manuscript_anchor",]
    if(nrow(d)) {
      changes <- data.frame(quantity=d$check,previous_reference=d$expected,computed=d$observed,
         difference=d$observed-d$expected,status=ifelse(d$passed,"agrees_with_supplied_precision","review_difference"),
         reason="Computed by new R core; this comparison does not alter the manuscript",
         source_file=check_path)
      fbst_write_table(changes,"manuscript_numeric_changes")
    }
  }
  status_path <- file.path(fbst_output_dir(),"execution_status.csv")
  status <- if(file.exists(status_path))read.csv(status_path,stringsAsFactors=FALSE) else data.frame()
  report <- c("# Reproduction status",paste("Computational version:",FBST_VERSION),
    paste("Profile:",fbst_profile()),paste("Generated UTC:",format(Sys.time(),tz="UTC",usetz=TRUE)),
    "", if(file.exists(manuscript))"The manuscript was read but not modified." else "No manuscript found (set FBST_MANUSCRIPT); traceability audit skipped.",
    "An output file is not itself evidence of numerical validation. Consult validation_checks.csv, design_diagnostics.csv and the execution status.",
    "Pilot results use reduced experimental sizes and omit application calibration; they must not replace manuscript tables.",
    "", "## Execution ledger")
  if(nrow(status)) {
    latest <- status[!duplicated(paste(status$stage,status$task),fromLast=TRUE),]
    report <- c(report,apply(latest,1,function(x)paste("-",x["stage"],"/",x["task"],":",x["status"],x["detail"])),
                "The full history and code fingerprints are retained in execution_status.csv.")
  }
  else report <- c(report,"No complete pipeline stages have been recorded.")
  report <- c(report,"","## Required text updates",
    "- Priors are now exactly KL=(a,a,a), informative=(25,25,25), conflict=(45,5,5); marginal N0 is lambda0+lambda1.",
    "- McNemar uses asymptotic chi-square without continuity correction throughout, with no rejection for zero discordant pairs.",
    "- Exact refers to summation over the finite sample space. Posterior/evidence integrations retain numerical error.",
    "- The complete sets of equivalent and near-optimal cutoffs are saved as unions of intervals; do not replace them by an interval hull if there are gaps.",
    "- Contrasts between groups report posterior summaries and directional probabilities, not a newly introduced four-dimensional FBST.",
    "- Design effects temper the likelihood while keeping the original cutoffs; no error recalibration is implied.",
    "- The appendix enumerates the sampling distribution of the retained 60x60 grid estimator. See appendix_method_changes.md.",
    "- Error risks averaged under different priors do not rank procedures under a common criterion.",
    "", "## Continuing",
    "Run the stage commands in README.md from the repository root. Complete cached designs and row checkpoints are reused.",
    "Failed stages retain their error and partial state; absence of an executed entry means the stage remains unverified.")
  writeLines(report,file.path(fbst_output_dir(),"reproduction_report.md"))
  invisible(status)
}
if(sys.nframe()==0L)fbst_reproduction_report()
