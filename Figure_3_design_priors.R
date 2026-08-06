#### Figure 3: Design priors for computing the Ave TIE and power ####

## Source the data of calibrated analysis priors
load("Figure_1_analysis_priors.RData")

# --- Very vague prior: N(-49, 8800^2) ----
dens_vague <- dnorm(
  mu_grid,
  mean = -49,
  sd   = 8800
)

# --- Skeptical prior: N(-90, 25^2) ----
dens_skeptical <- dnorm(
  mu_grid,
  mean = -90,
  sd   = 25
)

# --- Realistic prior: ----
## 0.51 N(-51.0, 19.9^2)
## + 0.44 N(-46.8, 7.6^2)
## + 0.05 N(-54.1, 51.7^2)

dens_realistic <-
  0.51 * dnorm(mu_grid, mean = -51.0, sd = 19.9) +
  0.44 * dnorm(mu_grid, mean = -46.8, sd = 7.6) +
  0.05 * dnorm(mu_grid, mean = -54.1, sd = 51.7)


## Figure 3

plot(
  mu_grid, dens_fpp,
  type = "l",
  lwd = 2,
  ylim = c(0,max(dens_vague,dens_skeptical,dens_realistic)),
  xlab = expression(theta[C]),
  ylab = "Density"#,
  #main = "Prior Density"
)

abline(v = xbar_ch, lty = 2)

lines(mu_grid, dens_vague, lwd = 2, lty = 2)
lines(mu_grid, dens_skeptical, lwd = 2, col=4, lty = 3)
lines(mu_grid, dens_realistic, lwd = 2, col=5, lty = 4)

lines(mu_grid, dens_map, lwd = 2, lty = 1, col=2)

lines(mu_grid, dens_rmap_1$rMAP,lwd = 2, col=3)
lines(mu_grid, dens_rmap_2$rMAP,lwd = 2, col=3, lty=2)

legend(
  "topright",
  legend = c(
    expression(paste("FPP (", alpha[0], " = 0.0465)")),
    expression(paste("MAP (", psi, " = 85.01)")),
    expression(
      paste(
        "rMAP"[1], " (",
        w, " = 0.80, ",
        psi, " = 61.20)"
      )
    ),
    expression(
      paste(
        "rMAP"[2], " (",
        w, " = 0.4445, ",
        psi, " = 30.80)"
      )
    ),
    expression(paste("Very vague: ", N(-49, 8800^2))),
    expression(paste("Skeptical: ", N(-90, 25^2))),
    "Realistic mixture prior"
  ),
  lwd = rep(2, 7),
  lty = c(1, 1, 1, 2, 2, 3, 4),
  col = c(1, 2, 3, 3, 1, 4, 5),
  bty = "n",
  cex = 0.8
)

