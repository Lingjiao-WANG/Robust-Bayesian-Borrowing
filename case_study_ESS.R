######################################
#### Real-data example
#### Borrowing from the control arm
#### Compute the EHSS
######################################

library(quadprog)
library(MCMCpack)
library(doParallel)
library(parallel)

## Source functions
source("case_study_funs.R")

## Load the calibrated analysis-prior results
load("Calibrated_analysis_priors.RData")

## Construct a grid of current or true control-arm means
mu_grid <- seq(
  xbar_ch - 200,
  xbar_ch + 200,
  length.out = 1e3
)


#---- FPP ----

## Compute the EHSS under the FPP
ESS_fpp <- alpha_0 * n_ch


#---- APP ----

## Compute the EHSS as a function of the observed control-arm sample mean
out_ESS_app_1 <- ESS_adaptive_power_prior_xbar(
  xbar_C_grid = mu_grid,
  xbar_ch = xbar_ch,
  n_C = n_C,
  sigma2 = sigma2,
  se_ch = se_ch,
  delta0 = delta0,
  weights = NULL
)

## Compute the expected EHSS given the true control-arm mean
out_ESS_app_2 <- ESS_adaptive_power_prior_muC(
  mu_grid = mu_grid,
  xbar_ch = xbar_ch,
  n_C = n_C,
  sigma2 = sigma2,
  se_ch = se_ch,
  delta0 = delta0,
  weights = NULL
)


#---- MAP ----

## Compute the EHSS as a function of the observed control-arm sample mean
out_ESS_map_1 <- ESS_MAP_xbarC(
  xbar_C_grid = mu_grid,
  xbar_ch = xbar_ch,
  n_C = n_C,
  sigma2 = sigma2,
  se_ch = se_ch,
  psi = psi_cal_map,
  n_tau = 401
)

## Compute the expected EHSS given the true control-arm mean
out_ESS_map_2 <- ESS_MAP_xbarC_given_mu(
  mu_grid = mu_grid,
  xbar_ch = xbar_ch,
  n_C = n_C,
  sigma2 = sigma2,
  se_ch = se_ch,
  psi = psi_cal_map,
  n_tau = 401
)


#---- rMAP ----

## rMAP_1

## Compute the EHSS as a function of the observed control-arm sample mean
out_ESS_rmap_1_1 <- ESS_rMAP_xbarC(
  xbar_C_grid = mu_grid,
  xbar_ch = xbar_ch,
  n_C = n_C,
  se_ch = se_ch,
  sigma2 = sigma2,
  psi = psi_cal_rmap,
  w_rob = w_rob,
  vague_sd = 88,
  weights = NULL,
  n_tau = 401
)

## Compute the expected EHSS given the true control-arm mean
out_ESS_rmap_1_2 <- ESS_rMAP_xbarC_given_mu(
  mu_grid = mu_grid,
  xbar_ch = xbar_ch,
  n_C = n_C,
  sigma2 = sigma2,
  se_ch = se_ch,
  psi = psi_cal_rmap,
  w_rob = w_rob,
  vague_sd = 88,
  weights = NULL,
  n_tau = 401,
  z_grid = seq(-5, 5, length.out = 1001)
)


## rMAP_2

## Compute the EHSS as a function of the observed control-arm sample mean
out_ESS_rmap_2_1 <- ESS_rMAP_xbarC(
  xbar_C_grid = mu_grid,
  xbar_ch = xbar_ch,
  n_C = n_C,
  se_ch = se_ch,
  sigma2 = sigma2,
  psi = 0.35 * sigma,
  w_rob = w_cal,
  vague_sd = 88,
  weights = NULL,
  n_tau = 401
)

## Compute the expected EHSS given the true control-arm mean
out_ESS_rmap_2_2 <- ESS_rMAP_xbarC_given_mu(
  mu_grid = mu_grid,
  xbar_ch = xbar_ch,
  n_C = n_C,
  sigma2 = sigma2,
  se_ch = se_ch,
  psi = 0.35 * sigma,
  w_rob = w_cal,
  vague_sd = 88,
  weights = NULL,
  n_tau = 401,
  z_grid = seq(-5, 5, length.out = 1001)
)


#### Plot Figure 3

## Define a common y-axis range for both panels
ylim_EHSS <- c(
  min(
    out_ESS_rmap_2_2$ESS_grid,
    out_ESS_rmap_1_2$ESS_grid
  ) - 2,
  max(out_ESS_app_2$ESS_grid) + 2
)

## Save current graphical settings
old_par <- par(no.readonly = TRUE)

## Change graphical settings
par(mfrow = c(1, 2))


#---- Panel A: EHSS as a function of observed xbar_C ----

plot(
  mu_grid,
  out_ESS_map_1$ESS_grid,
  type = "l",
  col = 2,
  lwd = 2,
  lty = 2,
  ylim = ylim_EHSS,
  xlab = expression(bar(y)[C]),
  ylab = "EHSS"
)

lines(
  mu_grid,
  out_ESS_rmap_1_1$ESS_grid,
  col = 3,
  lwd = 2,
  lty = 3
)

lines(
  mu_grid,
  out_ESS_rmap_2_1$ESS_grid,
  col = 4,
  lwd = 2,
  lty = 4
)

lines(
  mu_grid,
  out_ESS_app_1$ESS_grid,
  lwd = 2
)

abline(h = ESS_fpp, lty = 2)
abline(h = 0, lty = 3)

legend(
  "topright",
  legend = c(
    "APP",
    "MAP",
    expression(rMAP[1]),
    expression(rMAP[2]),
    "FPP",
    "EHSS = 0"
  ),
  col = c(1, 2, 3, 4, 1, 1),
  lty = c(1, 2, 3, 4, 2, 3),
  lwd = c(2, 2, 2, 2, 1, 1),
  bty = "n",
  cex = 0.9
)


#---- Panel B: Expected EHSS given the true control-arm mean ----

plot(
  mu_grid,
  out_ESS_app_2$ESS_grid,
  type = "l",
  lwd = 2,
  ylim = ylim_EHSS,
  xlab = expression(theta[C]),
  ylab = "EHSS"
)

lines(
  mu_grid,
  out_ESS_map_2$ESS_grid,
  col = 2,
  lwd = 2,
  lty = 2
)

lines(
  mu_grid,
  out_ESS_rmap_1_2$ESS_grid,
  col = 3,
  lwd = 2,
  lty = 3
)

lines(
  mu_grid,
  out_ESS_rmap_2_2$ESS_grid,
  col = 4,
  lwd = 2,
  lty = 4
)

abline(h = ESS_fpp, lty = 2)
abline(h = 0, lty = 3)

legend(
  "topright",
  legend = c(
    "APP",
    "MAP",
    expression(rMAP[1]),
    expression(rMAP[2]),
    "FPP",
    "EHSS = 0"
  ),
  col = c(1, 2, 3, 4, 1, 1),
  lty = c(1, 2, 3, 4, 2, 3),
  lwd = c(2, 2, 2, 2, 1, 1),
  bty = "n",
  cex = 0.9
)

## Restore the original graphical settings
par(old_par)


