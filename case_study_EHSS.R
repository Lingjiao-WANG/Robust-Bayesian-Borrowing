######################################
#### Real-data example
#### Borrowing from the control arm
#### Expected historical sample size
######################################

library(quadprog)
library(MCMCpack)
library(doParallel)
library(parallel)

## Source functions used in the real-data example
source("case_study_funs.R")

## Load the calibrated analysis-prior results
load("Rdatas/Calibrated_analysis_priors.RData")

## Construct a grid used to evaluate EHSS as a function of either
## the observed current control sample mean or the true current control mean
mu_grid <- seq(
  xbar_ch - 200,
  xbar_ch + 200,
  length.out = 1e3
)

#---- FPP ----

## Compute the EHSS under the FPP
ESS_fpp <- alpha_0 * n_ch


#---- APP ----

## Compute the EHSS as a function of the observed current control sample mean
out_ESS_app_1 <- ESS_adaptive_power_prior_xbar(
  xbar_C_grid = mu_grid,
  xbar_ch = xbar_ch,
  n_C = n_C,
  sigma2 = sigma2,
  se_ch = se_ch,
  delta0 = delta0,
  weights = NULL
)

## Compute the expected EHSS as a function of the true current control mean
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

## Compute the EHSS as a function of the observed current control sample mean
out_ESS_map_1 <- ESS_MAP_xbarC(
  xbar_C_grid = mu_grid,
  xbar_ch = xbar_ch,
  n_C = n_C,
  sigma2 = sigma2,
  se_ch = se_ch,
  psi = psi_cal_map,
  n_tau = 401
)

## Compute the expected EHSS as a function of the true current control mean
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

## Compute the EHSS as a function of the observed current control sample mean
out_ESS_rmap_1_1 <- ESS_rMAP_xbarC(
  xbar_C_grid = mu_grid,
  xbar_ch = xbar_ch,
  n_C = n_C,
  se_ch = se_ch,
  sigma2 = sigma2,
  psi = psi_cal_rmap,
  w_rob = w_rob,
  vague_sd = vague_sd,
  weights = NULL,
  n_tau = 401
)

## Compute the expected EHSS as a function of the true current control mean
out_ESS_rmap_1_2 <- ESS_rMAP_xbarC_given_mu(
  mu_grid = mu_grid,
  xbar_ch = xbar_ch,
  n_C = n_C,
  sigma2 = sigma2,
  se_ch = se_ch,
  psi = psi_cal_rmap,
  w_rob = w_rob,
  vague_sd = vague_sd,
  weights = NULL,
  n_tau = 401,
  z_grid = seq(-5, 5, length.out = 1001)
)


## rMAP_2

## Compute the EHSS as a function of the observed current control sample mean
out_ESS_rmap_2_1 <- ESS_rMAP_xbarC(
  xbar_C_grid = mu_grid,
  xbar_ch = xbar_ch,
  n_C = n_C,
  se_ch = se_ch,
  sigma2 = sigma2,
  psi = psi_rmap_2,
  w_rob = w_cal,
  vague_sd = vague_sd,
  weights = NULL,
  n_tau = 401
)

## Compute the expected EHSS as a function of the true current control mean
out_ESS_rmap_2_2 <- ESS_rMAP_xbarC_given_mu(
  mu_grid = mu_grid,
  xbar_ch = xbar_ch,
  n_C = n_C,
  sigma2 = sigma2,
  se_ch = se_ch,
  psi = psi_rmap_2,
  w_rob = w_cal,
  vague_sd = vague_sd,
  weights = NULL,
  n_tau = 401,
  z_grid = seq(-5, 5, length.out = 1001)
)

#save.image(file = "Rdatas/Case_study_EHSS.RData")

