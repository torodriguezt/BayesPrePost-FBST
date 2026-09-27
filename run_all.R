# Run from the repository root. See README.md for stages and execution status.
args <- commandArgs(trailingOnly=TRUE)
value <- function(key,default) {
  hit <- grep(paste0("^--",key,"="),args,value=TRUE)
  if(length(hit))sub(paste0("^--",key,"="),"",tail(hit,1)) else default
}
if("--help" %in% args) {
  cat("Rscript run_all.R --stage=validate|pilot|main|boundaries|applications|vaping|appendix|export|all --profile=full|pilot [--task=NAME] [--keep-going]\n")
  quit(status=0)
}
stage <- value("stage","all"); profile <- value("profile","full"); task <- value("task","all")
Sys.setenv(FBST_PROFILE=profile)
source("R/pipeline_helpers.R")
source("R/validate_numerics.R")
source("R/manuscript_experiments.R")
source("R/manuscript_applications.R")
source("R/application_vaping.R")
source("R/04_simulation_estimation_figures.R")
source("R/reproduction_report.R")
fbst_record_environment()
writeLines(c(paste("pid",Sys.getpid()),paste("stage",stage),paste("task",task),
             paste("profile",profile),paste("started_utc",format(Sys.time(),tz="UTC",usetz=TRUE))),
           file.path(fbst_output_dir(),"active_run.txt"))
tasks <- list(
  validate=list(priors=fbst_validate_priors,sampling=fbst_validate_sampling,numerics=fbst_validate_core),
  pilot=list(runtime=fbst_run_pilot),
  main=list(wald_benchmark=fbst_run_wald_benchmark,power=fbst_run_power,prior_sensitivity=fbst_run_sensitivity,
            paired=fbst_run_paired,coverage=fbst_run_coverage,
            null_priors=fbst_run_null_priors,independent=fbst_run_independent_priors),
  boundaries=list(boundaries=fbst_run_boundaries,references=fbst_run_references,
                  unequal=fbst_run_unequal,restricted_prior=fbst_run_restricted_prior),
  applications=list(data=fbst_audit_application_data,tvsfp=run_application_tvsfp,
                     toenail=run_application_toenail,tables=run_application_tables,
                     figures=run_application_figures,application_validation=fbst_application_validation),
  vaping=list(data=fbst_audit_vaping_data,results=run_application_vaping,
              tables=run_vaping_tables,figures=run_vaping_figures),
  appendix=list(estimation=fbst_run_estimation),
  export=list(report=fbst_reproduction_report))
if(!stage %in% c(names(tasks),"all"))stop("Unknown stage: ",stage)
selected <- if(stage=="all")names(tasks) else stage
if(task!="all") {
  if(length(selected)!=1L || !task %in% names(tasks[[selected]]))
    stop("--task must name a task in one selected stage: ",paste(names(tasks[[selected[1]]]),collapse=", "))
  tasks[[selected]] <- tasks[[selected]][task]
}
status_path <- file.path(fbst_output_dir(),"execution_status.csv")
status <- if(file.exists(status_path))read.csv(status_path,stringsAsFactors=FALSE) else NULL
failed <- FALSE
for(st in selected) for(nm in names(tasks[[st]])) {
  message("\n",st," / ",nm," [",profile,"]")
  started <- Sys.time(); cpu <- proc.time()[[3]]; err <- ""
  status <- if(file.exists(status_path))read.csv(status_path,stringsAsFactors=FALSE) else NULL
  start_row <- data.frame(stage=st,task=nm,profile=profile,
    started_utc=format(started,tz="UTC",usetz=TRUE),seconds=0,
    status="running",detail="",
    code_hash=fbst_hash(fbst_code_hashes()),stringsAsFactors=FALSE)
  write.csv(rbind(status,start_row),status_path,row.names=FALSE)
  tryCatch(tasks[[st]][[nm]](),error=function(e) {err <<- conditionMessage(e)})
  row <- data.frame(stage=st,task=nm,profile=profile,
    started_utc=format(started,tz="UTC",usetz=TRUE),seconds=proc.time()[[3]]-cpu,
    status=if(nzchar(err))"failed" else "executed",
    detail=if(nzchar(err))err else if(profile=="pilot")"Reduced pilot scope; not a full manuscript reproduction" else "",
    code_hash=fbst_hash(fbst_code_hashes()),stringsAsFactors=FALSE)
  status <- rbind(read.csv(status_path,stringsAsFactors=FALSE),row)
  write.csv(status,status_path,row.names=FALSE)
  if(nzchar(err)) {
    failed <- TRUE; message("FAILED: ",err)
    if(!"--keep-going" %in% args) {
      try(fbst_reproduction_report(),silent=TRUE)
      stop("Stage interrupted; state retained in ",fbst_output_dir(),
           ". Re-run this command to resume cached work.")
    }
  }
}
if(stage!="export")fbst_reproduction_report()
if(failed)quit(status=1)
