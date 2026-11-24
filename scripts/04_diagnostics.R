# Residual diagnostics of the top-ranked model of each index
library(ape)
library(lmtest)
library(sandwich)

dat <- read.csv("results/ses_landscape_data.csv")
dist_km <- as.matrix(read.csv("data/distance_matrix_km.csv", row.names = 1))
w <- 1 / dist_km
diag(w) <- 0

moran_perm <- function(res, perms = 9999) {
  stat <- function(x) {
    x <- x - mean(x)
    (length(x) / sum(w)) * sum(w * outer(x, x)) / sum(x^2)
  }
  obs <- stat(res)
  expected <- -1 / (length(res) - 1)
  sims <- replicate(perms, stat(sample(res)))
  (sum(abs(sims - expected) >= abs(obs - expected)) + 1) / (perms + 1)
}

models <- list("SES.PD ~ forest cover" = lm(SES_PD ~ Cobertura_pct, dat),
               "SES.MPD ~ 1" = lm(SES_MPD ~ 1, dat),
               "SES.MNTD ~ edge density" = lm(SES_MNTD ~ Densidade_borda, dat))

set.seed(123)
out <- do.call(rbind, lapply(names(models), function(nm) {
  m <- models[[nm]]
  r <- residuals(m)
  has_x <- length(coef(m)) > 1
  mi <- Moran.I(r, w)
  data.frame(model = nm,
             shapiro_p = shapiro.test(r)$p.value,
             bp_p = if (has_x) bptest(m)$p.value else NA,
             hc3_p_t = if (has_x) coeftest(m, vcov = vcovHC(m, type = "HC3"))[2, 4] else NA,
             hc3_p_z = if (has_x) coeftest(m, vcov = vcovHC(m, type = "HC3"), df = Inf)[2, 4] else NA,
             moran_I = mi$observed,
             moran_p_analytical = mi$p.value,
             moran_p_perm = moran_perm(r))
}))
write.csv(out, "results/diagnostics.csv", row.names = FALSE)
print(out, digits = 3, row.names = FALSE)
