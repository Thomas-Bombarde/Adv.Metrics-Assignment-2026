# make_figures.R
# Three diagnostic figures from the Monte Carlo run.
#
#   Fig A. Density of the raw CCE slope t-statistic under H0 vs N(0,1) and
#          BCH's fixed-G limit, sqrt(G / (G - 1)) * t_{G-1}.
#   Fig B. Empirical rejection rates at 5% nominal level across gamma and G.
#   Fig C. Sampling density of the asymptotic slope-variance object
#          (Q^{-1} Vhat Q^{-1})_{slope,slope}. The actual variance of beta_hat
#          divides this object by N and therefore always shrinks mechanically.
#
# Usage from the report project root:
#   Rscript code_and_data/make_figures.R

source("code_and_data/paths.R")
#
# Audience guide:
# BCH do not claim that the cluster covariance estimator converges to the true
# variance when G is fixed. They use an asymptotic sequence where each group
# grows but the number of groups stays small. Group scores become approximately
# independent normals, but there are still only G draws from the group-score
# distribution. The CCE therefore has an irreducibly random fixed-G limit.
# The fixed-G t reference distribution accounts for that randomness.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

df <- readRDS(file.path(results_dir, "sim_results.rds"))
required_cols <- c("t_hac", "t_cce", "t_cce_rescaled",
                   "V_slope_asymp", "V_hac_asymp")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("sim_results.rds is missing columns created by the current helpers.R: ",
       paste(missing_cols, collapse = ", "),
       "\nRerun simulation/run_sim.R before making figures.")
}

# Figures are embedded in Report.Rmd with in-document captions. Do not add
# plot-level titles/subtitles here; put that content in the Rmd captions.
fig_dir <- figures_dir
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ----- Figure A: t-statistic density under H0 -----
# Focus on K = 36 (largest group size), gamma = 0.6 (strong spatial dependence),
# and G = 4. This isolates BCH's main message: the raw CCE statistic is not
# standard normal under fixed G. Its limit is sqrt(G / (G - 1)) times t_{G-1}.
fa_data <- df %>% filter(K == 36, gamma == 0.6, G == 4)

G_ref <- 4
xx <- seq(-6, 6, length.out = 401)
scale_bch <- sqrt(G_ref / (G_ref - 1))
ref <- bind_rows(
  data.frame(t = xx, dens = dnorm(xx), ref = "N(0,1)"),
  data.frame(
    t = xx,
    dens = dt(xx / scale_bch, df = G_ref - 1) / scale_bch,
    ref = sprintf("sqrt(%d/%d) * t_{%d}", G_ref, G_ref - 1, G_ref - 1)
  )
)

fig_A <- ggplot(fa_data, aes(x = t_cce)) +
  geom_density(colour = "#1b6ca8", linewidth = 0.8, fill = NA) +
  geom_line(data = ref, aes(x = t, y = dens, linetype = ref), colour = "black",
            linewidth = 0.5, inherit.aes = FALSE) +
  scale_linetype_manual(name = "Reference", values = c("solid", "dashed")) +
  coord_cartesian(xlim = c(-6, 6)) +
  labs(
    x = "Raw CCE t-statistic for H0: beta = 1",
    y = "Density"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "figA_t_density.png"), fig_A,
       width = 8, height = 4.5, dpi = 160)

# ----- Figure B: rejection rates at 5% across gamma -----
# A normal reference distribution treats Vhat as if it had converged to a
# constant. BCH's correction rescales the CCE statistic by sqrt((G - 1) / G)
# and compares it with t_{G-1}, explicitly pricing the uncertainty from having
# only G approximately independent group scores.
crit_n <- qnorm(0.975)
fb <- df %>%
  group_by(gamma, K, G) %>%
  summarise(
    rej_hc0       = mean(abs(t_hc0) > crit_n),
    rej_hac       = mean(abs(t_hac) > crit_n),
    rej_cce_norm  = mean(abs(t_cce) > crit_n),
    rej_cce_t     = mean(abs(t_cce_rescaled) > qt(0.975, df = unique(G) - 1)),
    .groups = "drop"
  ) %>%
  pivot_longer(starts_with("rej_"), names_to = "method", values_to = "rate") %>%
  mutate(method = recode(method,
    rej_hc0      = "HC0 + N(0,1)",
    rej_hac      = "HAC (Bartlett, H=5) + N(0,1)",
    rej_cce_norm = "CCE + N(0,1)",
    rej_cce_t    = "CCE + fixed-G t_{G-1}"
  ))

# Plot for K = 36 and G = 9. Keeping one group count makes the size
# comparison easier to read in the report figure.
fb_plot <- fb %>% filter(K == 36, G == 9)

fig_B <- ggplot(fb_plot, aes(x = factor(gamma), y = rate,
                             fill = method, group = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey30") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = expression(gamma~"(spatial MA decay)"),
       y = "Empirical rejection rate at 5% nominal level",
       fill = "Inference procedure") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave(file.path(fig_dir, "figB_rejection_rates.png"), fig_B,
       width = 9, height = 5, dpi = 160)

# ----- Figure C: sampling density of Vhat on asymptotic-variance scale -----
# Object: the slope element of Q^{-1} Vhat Q^{-1}. This is the object appearing
# in the denominator of sqrt(N)(beta_hat - beta). It should not be divided by N
# for this diagnostic. If we plotted V_slope = this object / N, every estimator
# would collapse to zero by construction.
#
# Hold G = 4 fixed for CCE and vary K. Each K gives a larger group size L_N.
# Under BCH fixed-G asymptotics, the CCE density should remain visibly spread
# out as L_N grows. The HAC panel is included as a non-fixed-G local covariance
# contrast: its sampling variation concentrates as the lattice grows. In this
# Monte Carlo file H is fixed at 5, so read it as a concentration comparison;
# a textbook consistent Bartlett HAC sequence would let H grow slowly with K.
K_sorted <- sort(unique(df$K))
K_labels <- sprintf("L_N = %d  (K = %d)", (K_sorted / 2)^2, K_sorted)

fc_cce <- df %>%
  filter(G == 4, gamma == 0.6) %>%
  mutate(estimator = "CCE  (G = 4, fixed)",
         L_label = factor(sprintf("L_N = %d  (K = %d)", (K / 2)^2, K),
                          levels = K_labels),
         V = V_slope_asymp) %>%
  select(K, L_label, estimator, V)

# HAC: collapse to one row per (rep, K) since V_hac_asymp does not depend on G.
fc_hac <- df %>%
  filter(gamma == 0.6, G == min(G)) %>%
  mutate(estimator = "HAC  (Bartlett, H = 5)",
         L_label = factor(sprintf("L_N = %d  (K = %d)", (K / 2)^2, K),
                          levels = K_labels),
         V = V_hac_asymp) %>%
  select(K, L_label, estimator, V)

fc <- bind_rows(fc_cce, fc_hac)

fig_C <- ggplot(fc, aes(x = V, colour = L_label)) +
  geom_density(linewidth = 0.7, fill = NA) +
  facet_wrap(~ estimator, scales = "free_y", nrow = 1) +
  labs(x = expression("Asymptotic slope variance,  "*Q^{-1}*hat(V)[N]*Q^{-1}*" (slope element)"),
       y = "Density",
       colour = "Group size") +
  coord_cartesian(xlim = c(0, quantile(fc$V, 0.99))) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "figC_Vhat_density.png"), fig_C,
       width = 10, height = 4.5, dpi = 160)

# ----- Print a compact rejection-rate table to console -----
cat("\nEmpirical rejection rates at 5% nominal level (K = 36):\n")
fb_table <- fb %>%
  filter(K == 36) %>%
  pivot_wider(names_from = method, values_from = rate) %>%
  arrange(gamma, G)
print(as.data.frame(fb_table), digits = 3)

notes <- c(
  "---",
  "title: \"Fixed-G Cluster Inference in Bester, Conley, and Hansen\"",
  "subtitle: \"Monte Carlo Figures for a Spatial Moving-Average Design\"",
  "author: \"Graduate Econometrics Simulation Note\"",
  "date: \"`r format(Sys.Date(), '%B %d, %Y')`\"",
  "output:",
  "  pdf_document:",
  "    latex_engine: xelatex",
  "    number_sections: true",
  "    toc: false",
  "fontsize: 11pt",
  "geometry: margin=1in",
  "header-includes:",
  "  - \\usepackage{booktabs}",
  "  - \\usepackage{caption}",
  "  - \\captionsetup{font=small,labelfont=bf}",
  "  - \\usepackage{float}",
  "---",
  "",
  "```{r setup, include=FALSE}",
  "knitr::opts_chunk$set(",
  "  echo = FALSE,",
  "  message = FALSE,",
  "  warning = FALSE,",
  "  fig.align = \"center\",",
  "  out.width = \"0.92\\\\linewidth\"",
  ")",
  "```",
  "",
  "## Overview",
  "",
  "These figures illustrate Bester, Conley, and Hansen's fixed-\\(G\\) argument in the spatial moving-average Monte Carlo. The key point is not that cluster scores converge to constants. Within each large group, the score obeys a central limit theorem; across groups, the scores become asymptotically independent. With \\(G\\) fixed, the cluster covariance estimator (CCE) is therefore built from only \\(G\\) approximately independent random scores, so its limit remains random.",
  "",
  "The simulation uses the paper's spatial moving-average design on a \\(K \\times K\\) lattice:",
  "\\[",
  "x_s = \\sum_{\\lVert h \\rVert \\leq 2}\\gamma^{\\lVert h \\rVert}v_{s+h},",
  "\\qquad",
  "\\varepsilon_s = \\sum_{\\lVert h \\rVert \\leq 2}\\gamma^{\\lVert h \\rVert}u_{s+h},",
  "\\]",
  "with \\(y_s = x_s\\beta + \\varepsilon_s\\) and \\(\\beta = 1\\). The null hypothesis is true in every Monte Carlo replication.",
  "",
  "## Figure A: CCE t-statistic",
  "",
  "```{r fig-a, fig.cap=\"Raw CCE t-statistic under the null. The blue density is the Monte Carlo distribution for K = 36, gamma = 0.6, and G = 4. The dashed reference is the BCH fixed-G limit.\"}",
  "knitr::include_graphics(\"figA_t_density.png\")",
  "```",
  "",
  "The null is true, so a conventional plug-in argument would suggest a \\(N(0,1)\\) reference distribution. BCH instead keep the number of groups fixed and let each group grow. The raw CCE statistic has limiting reference distribution",
  "\\[",
  "\\sqrt{\\frac{G}{G-1}}\\,t_{G-1}.",
  "\\]",
  "For \\(G=4\\), this is much wider than \\(N(0,1)\\), because the denominator is estimated from only four approximately independent group scores.",
  "",
  "## Figure B: Rejection Rates",
  "",
  "```{r fig-b, fig.cap=\"Empirical rejection rates for nominal 5 percent tests under the true null. The dashed horizontal line is the target size.\"}",
  "knitr::include_graphics(\"figB_rejection_rates.png\")",
  "```",
  "",
  "The dashed line is the nominal 5% rejection rate. Procedures using \\(N(0,1)\\) critical values treat the variance estimate as essentially known. The BCH procedure rescales the CCE t-statistic by \\(\\sqrt{(G-1)/G}\\) and uses \\(t_{G-1}\\) critical values, which accounts for fixed-\\(G\\) variance-estimation uncertainty.",
  "",
  "## Figure C: Variance Estimator",
  "",
  "```{r fig-c, fig.cap=\"Sampling distribution of the asymptotic-scale slope variance estimate. The CCE panel keeps G fixed while group size grows.\"}",
  "knitr::include_graphics(\"figC_Vhat_density.png\")",
  "```",
  "",
  "The plotted object is",
  "\\[",
  "\\left(Q^{-1}\\hat V_N Q^{-1}\\right)_{\\text{slope},\\text{slope}},",
  "\\]",
  "the asymptotic variance scale for \\(\\sqrt{N}(\\hat\\beta-\\beta)\\). This is deliberately not divided by \\(N\\). The actual variance of \\(\\hat\\beta\\) shrinks like \\(1/N\\) for any sensible estimator, so plotting that object would miss the BCH point.",
  "",
  "The CCE density remains dispersed as \\(L_N\\) grows because fixed \\(G\\) leaves only \\(G\\) pieces of information about the group-score variance. The HAC panel is a local-covariance contrast whose sampling variation concentrates as the lattice grows. In this simulation \\(H\\) is fixed at 5, so it should be read as a concentration comparison rather than an exact implementation of a growing-bandwidth HAC sequence.",
  "",
  "## Takeaway",
  "",
  "The figures separate two ideas that are easy to conflate. First, large groups justify approximately normal and independent group scores. Second, fixed \\(G\\) means the variance estimator is still random, because it is based on only \\(G\\) group-level pieces of information. BCH's fixed-\\(G\\) critical values target that second source of uncertainty."
)
notes_rmd <- file.path(fig_dir, "figure_notes.Rmd")
writeLines(notes, notes_rmd)

if (requireNamespace("rmarkdown", quietly = TRUE)) {
  rmarkdown::render(notes_rmd, quiet = TRUE)
}

cat("\nFigures written to:", fig_dir, "\n")
cat("Figure notes written to:", notes_rmd, "\n")
