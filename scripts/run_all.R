# run from the project root: Rscript scripts/run_all.R
if (!file.exists("data/comm.txt")) stop("run this script from the project root")

scripts <- c("01_ses.R", "02_univariate_models.R", "03_model_selection.R",
             "04_diagnostics.R", "05_acoustic_overlap.R", "06_phylogenetic_signal.R")

# picante is only needed for step 1; the SES results are also stored in results/
if (!requireNamespace("picante", quietly = TRUE)) {
  message("picante not installed: skipping 01_ses.R and using the stored SES results")
  scripts <- scripts[-1]
}

for (s in scripts) {
  message("running ", s)
  source(file.path("scripts", s), local = new.env())
}
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
