# Blomberg's K of the acoustic parameters, tree pruned to the recorded species
library(ape)

traits <- read.table("data/traits.txt", header = TRUE, sep = "\t", row.names = 1)
comm <- read.table("data/comm.txt", header = TRUE, row.names = 1)
phy <- read.tree("data/phy.tre")
phy <- drop.tip(phy, setdiff(phy$tip.label, colnames(comm)))
traits <- traits[phy$tip.label, ]

reps <- 9999
set.seed(123)

cov_mat <- vcv(phy)
inv <- solve(cov_mat)
n <- nrow(cov_mat)
one <- rep(1, n)
expected <- (sum(diag(cov_mat)) - n / (t(one) %*% inv %*% one)) / (n - 1)

blomberg_k <- function(x) {
  a <- as.numeric((t(one) %*% inv %*% x) / (t(one) %*% inv %*% one))
  r <- x - a
  mse0 <- sum(r^2) / (n - 1)
  mse <- as.numeric(t(r) %*% inv %*% r) / (n - 1)
  as.numeric((mse0 / mse) / expected)
}

res <- do.call(rbind, lapply(names(traits), function(p) {
  x <- setNames(traits[[p]], rownames(traits))
  k <- blomberg_k(x)
  sims <- replicate(reps, blomberg_k(sample(x)))
  data.frame(parameter = p, K = round(k, 3), p_value = (sum(sims >= k) + 1) / (reps + 1))
}))
dir.create("results", showWarnings = FALSE)
write.csv(res, "results/phylogenetic_signal.csv", row.names = FALSE)
print(res)

if (requireNamespace("picante", quietly = TRUE)) {
  k_picante <- sapply(names(traits), function(p) round(picante::Kcalc(setNames(traits[[p]], rownames(traits)), phy), 3))
  cat("K from picante::Kcalc:", k_picante, "\n")
}
