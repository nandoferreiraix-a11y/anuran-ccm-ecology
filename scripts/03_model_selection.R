# AICc selection among 11 candidate models per index:
# null, four single-variable and six two-variable models (z-scored predictors)
dat <- read.csv("results/ses_landscape_data.csv")
vars <- c("Forest cover" = "Cobertura_pct", "Isolation" = "Isolamento_m",
          "Edge density" = "Densidade_borda", "Habitat heterogeneity" = "Heterogeneidade_SHDI")
resp <- c("SES.PD" = "SES_PD", "SES.MPD" = "SES_MPD", "SES.MNTD" = "SES_MNTD")

for (v in vars) dat[[paste0(v, "_z")]] <- as.numeric(scale(dat[[v]]))
n <- nrow(dat)

aicc <- function(m) {
  k <- length(coef(m)) + 1   # coefficients + residual variance
  AIC(m) + 2 * k * (k + 1) / (n - k - 1)
}
vif2 <- function(s) 1 / (1 - cor(dat[[paste0(s[1], "_z")]], dat[[paste0(s[2], "_z")]])^2)
label <- function(s) if (length(s) == 0) "Null" else paste(names(vars)[match(s, vars)], collapse = " + ")

sets <- c(list(character(0)), as.list(vars), combn(vars, 2, simplify = FALSE))

res <- do.call(rbind, lapply(names(resp), function(r) {
  tab <- do.call(rbind, lapply(sets, function(s) {
    rhs <- if (length(s) == 0) "1" else paste0(s, "_z", collapse = " + ")
    m <- lm(as.formula(paste(resp[[r]], "~", rhs)), dat)
    sm <- summary(m)
    one <- length(s) == 1
    data.frame(index = r, model = label(s), K = length(coef(m)) + 1,
               R2adj = sm$adj.r.squared, AICc = aicc(m),
               VIF_max = if (length(s) == 2) vif2(s) else NA,
               beta_std = if (one) coef(m)[2] else NA,
               CI_low = if (one) confint(m)[2, 1] else NA,
               CI_high = if (one) confint(m)[2, 2] else NA,
               p = if (one) coef(sm)[2, 4] else NA)
  }))
  tab <- tab[order(tab$AICc), ]
  tab$dAICc <- tab$AICc - min(tab$AICc)
  tab$weight <- exp(-0.5 * tab$dAICc) / sum(exp(-0.5 * tab$dAICc))
  tab$evidence_ratio_vs_null <- tab$weight / tab$weight[tab$model == "Null"]
  tab
}))

stopifnot(max(res$VIF_max, na.rm = TRUE) < 3)
write.csv(res, "results/model_selection.csv", row.names = FALSE)
print(res[res$dAICc <= 2 | res$model == "Null",
          c("index", "model", "R2adj", "AICc", "dAICc", "weight", "evidence_ratio_vs_null")],
      digits = 3, row.names = FALSE)
