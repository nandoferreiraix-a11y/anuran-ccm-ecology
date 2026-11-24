# Univariate models, multiple-testing corrections and Figure 2
library(ggplot2)
library(patchwork)

dat <- read.csv("results/ses_landscape_data.csv")
responses <- c("SES_PD", "SES_MPD", "SES_MNTD")
predictors <- c("Cobertura_pct", "Isolamento_m", "Densidade_borda", "Heterogeneidade_SHDI")

uni <- do.call(rbind, lapply(responses, function(y) {
  do.call(rbind, lapply(predictors, function(x) {
    s <- summary(lm(reformulate(x, y), data = dat))
    data.frame(response = y, predictor = x, r2 = s$r.squared, r2_adj = s$adj.r.squared,
               estimate = coef(s)[x, "Estimate"], p_value = coef(s)[x, "Pr(>|t|)"])
  }))
}))
uni$p_bonferroni <- p.adjust(uni$p_value, "bonferroni")
uni$q_fdr_all <- p.adjust(uni$p_value, "BH")
uni$q_fdr_within_response <- ave(uni$p_value, uni$response, FUN = function(p) p.adjust(p, "BH"))
write.csv(uni, "results/univariate_models.csv", row.names = FALSE)

# one test per landscape, corrected across the 18 landscapes
land <- dat[order(dat$landscape), c("landscape", "leec_code", "SES_MPD", "P_MPD", "SES_MNTD", "P_MNTD")]
land$q_fdr_mpd <- p.adjust(land$P_MPD, "BH")
land$q_fdr_mntd <- p.adjust(land$P_MNTD, "BH")
land$p_bonferroni_mntd <- p.adjust(land$P_MNTD, "bonferroni")
write.csv(land, "results/fdr_by_landscape.csv", row.names = FALSE)

# Figure 2
panel <- function(xvar, yvar, xlab, ylab) {
  ggplot(dat, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
    geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey40", alpha = 0.2, linewidth = 0.7) +
    geom_point(size = 2.2, color = "grey15") +
    labs(x = xlab, y = ylab) +
    theme_bw(base_size = 11) +
    theme(panel.grid.major = element_line(color = "grey85", linetype = "dashed", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          axis.title = element_text(face = "bold"))
}
fig2 <- panel("Cobertura_pct", "SES_PD", "Forest cover (%)", "Phylogenetic diversity (SES.PD)") +
  panel("Isolamento_m", "SES_PD", "Isolation (m)", "Phylogenetic diversity (SES.PD)") +
  panel("Densidade_borda", "SES_MNTD", "Edge density", "Mean nearest taxon distance (SES.MNTD)") +
  plot_layout(ncol = 3) + plot_annotation(tag_levels = "A")
dir.create("figures", showWarnings = FALSE)
ggsave("figures/figure2.png", fig2, width = 15, height = 4.8, dpi = 300)
