# SES of PD, MPD and MNTD for each landscape
library(ape)
library(picante)

seed <- 123
n_runs <- 999
null_model <- "taxa.labels"

comm_raw <- read.table("data/comm.txt", row.names = 1, header = TRUE)
phy <- read.tree("data/phy.tre")
land <- read.csv("data/landscape_variables.csv")

# prune the tree to the recorded species
matched <- match.phylo.comm(phy, comm_raw)
phy_matched <- matched$phy
comm_abund <- matched$comm
comm_pa <- comm_abund
comm_pa[comm_pa > 0] <- 1
phydist <- cophenetic(phy_matched)

run_ses <- function(comm, abundance_weighted) {
  set.seed(seed)
  pd <- ses.pd(comm, phy_matched, null.model = null_model, runs = n_runs, include.root = FALSE)
  set.seed(seed)
  mpd <- ses.mpd(comm, phydist, null.model = null_model,
                 abundance.weighted = abundance_weighted, runs = n_runs)
  set.seed(seed)
  mntd <- ses.mntd(comm, phydist, null.model = null_model,
                   abundance.weighted = abundance_weighted, runs = n_runs)
  data.frame(leec_code = rownames(comm),
             SES_PD = pd$pd.obs.z, P_PD = pd$pd.obs.p,
             SES_MPD = mpd$mpd.obs.z, P_MPD = mpd$mpd.obs.p,
             SES_MNTD = mntd$mntd.obs.z, P_MNTD = mntd$mntd.obs.p)
}

ses_pa <- run_ses(comm_pa, abundance_weighted = FALSE)
ses_abund <- run_ses(comm_abund, abundance_weighted = TRUE)

dir.create("results", showWarnings = FALSE)
write.csv(ses_pa, "results/ses_presence_absence.csv", row.names = FALSE)
write.csv(ses_abund, "results/ses_abundance.csv", row.names = FALSE)
write.csv(merge(ses_abund, land, by = "leec_code"), "results/ses_landscape_data.csv", row.names = FALSE)
