# Acoustic niche overlap (Horn-Morisita) per landscape and regressions for hypothesis 3
traits <- read.table("data/traits.txt", header = TRUE, sep = "\t", row.names = 1)
comm <- read.table("data/comm.txt", header = TRUE, row.names = 1)
dat <- read.csv("results/ses_landscape_data.csv")

eps <- 0.001
traits <- traits[colnames(comm), ]

# rescale each parameter to 0-1 across species
minmax <- function(x) (x - min(x)) / (max(x) - min(x))
traits01 <- sapply(traits, minmax) + eps
rownames(traits01) <- rownames(traits)

horn_morisita <- function(x, y) {
  X <- sum(x); Y <- sum(y)
  2 * sum(x * y) / ((sum(x^2) / X^2 + sum(y^2) / Y^2) * X * Y)
}
sp <- rownames(traits01)
overlap <- matrix(NA_real_, length(sp), length(sp), dimnames = list(sp, sp))
for (i in seq_along(sp)) for (j in seq_along(sp)) overlap[i, j] <- horn_morisita(traits01[i, ], traits01[j, ])

# landscape value: mean overlap among the species recorded there
by_landscape <- sapply(rownames(comm), function(l) {
  m <- overlap[comm[l, ] > 0, comm[l, ] > 0, drop = FALSE]
  mean(m[lower.tri(m)])
})
out <- data.frame(landscape = seq_len(nrow(comm)), leec_code = rownames(comm),
                  richness = rowSums(comm > 0), Horn_Morisita = round(unname(by_landscape), 5))

dir.create("results", showWarnings = FALSE)
write.csv(overlap, "results/overlap_species_matrix.csv")
write.csv(out, "results/overlap_by_landscape.csv", row.names = FALSE)
cat(sprintf("Horn-Morisita range: %.3f to %.3f\n", min(out$Horn_Morisita), max(out$Horn_Morisita)))

# hypothesis 3
d <- merge(out, dat, by = "leec_code")
d <- d[match(out$leec_code, d$leec_code), ]
preds <- c(forest_cover = "Cobertura_pct", isolation = "Isolamento_m", edge_density = "Densidade_borda",
           SHDI = "Heterogeneidade_SHDI", SES.MPD = "SES_MPD", SES.MNTD = "SES_MNTD")
reg <- do.call(rbind, lapply(names(preds), function(p) {
  s <- summary(lm(reformulate(preds[[p]], "Horn_Morisita"), data = d))
  data.frame(predictor = p, slope = coef(s)[2, 1], R2adj = s$adj.r.squared, p = coef(s)[2, 4])
}))
write.csv(reg, "results/overlap_regressions.csv", row.names = FALSE)
print(reg, digits = 3, row.names = FALSE)
