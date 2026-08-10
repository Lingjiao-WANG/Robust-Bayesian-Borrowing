#### Figure 1: Calibrated analysis priors

## Source functions
source("case_study_funs.R")

## Load the calibrated analysis-prior results
load("Calibrated_analysis_priors.RData")

## Construct a grid of true control-arm means
mu_grid <- seq(
  xbar_ch - 200,
  xbar_ch + 200,
  length.out = 1e5
)

#---- FPP ----

## Compute the calibrated FPP prior density
sd_fpp <- se_ch / sqrt(alpha_0)

dens_fpp <- dnorm(
  mu_grid,
  mean = xbar_ch,
  sd = sd_fpp
)

#---- Vague prior ----

## Compute the vague prior density
dens_vague <- dnorm(
  mu_grid,
  mean = xbar_ch,
  sd = 88
)

#---- MAP ----

## Compute the calibrated MAP prior density
dens_map <- dmap_prior(
  theta_C = mu_grid,
  xbar_ch = xbar_ch,
  se_ch = se_ch,
  psi = psi_cal_map,
  n_tau = 401
)

#---- rMAP ----

## Compute the calibrated rMAP_1 prior density
dens_rmap_1 <- drmap_prior(
  mu_grid = mu_grid,
  xbar_ch = xbar_ch,
  se_ch = se_ch,
  psi = psi_cal_rmap,
  w_rob = w_rob,
  vague_sd = 88,
  n_tau = 401
)

## Compute the calibrated rMAP_2 prior density
dens_rmap_2 <- drmap_prior(
  mu_grid = mu_grid,
  xbar_ch = xbar_ch,
  se_ch = se_ch,
  psi = 0.35 * sigma,
  w_rob = w_cal,
  vague_sd = 88,
  n_tau = 401
)

#### Plot Figure 1

plot(
  mu_grid,
  dens_fpp,
  type = "l",
  lwd = 2,
  xlab = expression(theta[C]),
  ylab = "Density"
)

## Add a vertical reference line at the historical control mean
abline(v = xbar_ch, lty = 2)

## Add the remaining analysis-prior densities
lines(mu_grid, dens_vague, lwd = 2, lty = 2)
lines(mu_grid, dens_map, lwd = 2, lty = 1, col = 2)
lines(mu_grid, dens_rmap_1$rMAP, lwd = 2, col = 3)
lines(mu_grid, dens_rmap_2$rMAP, lwd = 2, lty = 2, col = 3)

## Add the legend
legend(
  "topright",
  legend = c(
    bquote("FPP (" * alpha[0] * " = " * .(round(alpha_0, 4)) * ")"),
    bquote("MAP (" * psi * " = " * .(round(psi_cal_map, 2)) * ")"),
    bquote("rMAP"[1] * " (" * w * " = " * .(round(w_rob, 2)) *
             ", " * psi * " = " * .(round(psi_cal_rmap, 2)) * ")"),
    bquote("rMAP"[2] * " (" * w * " = " * .(round(w_cal, 2)) *
             ", " * psi * " = " * .(round(0.35 * sigma, 2)) * ")"),
    expression(paste("Vague prior: ", N(-49, 88^2)))
  ),
  lwd = 2,
  lty = c(1, 1, 1, 2, 2),
  col = c(1, 2, 3, 3, 1),
  bty = "n"
)