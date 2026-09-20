######################################
#### Real-data example
#### Borrowing from the control arm
#### Type I error, power, and calibrated power gain
######################################

library(quadprog)
library(MCMCpack)
library(doParallel)
library(parallel)

## Source functions used in the real-data example
source("case_study_funs.R")

## Load the calibrated analysis-prior results
load("Rdatas/Calibrated_analysis_priors.RData")

## Compute the standard error of the treatment-effect estimator
## under the no-borrowing analysis
SE_nob <- sqrt(
  sigma2 / n_T +
    sigma2 / n_C
)

## Construct a grid of true current control-arm means
mu_grid <- seq(
  xbar_ch - 200,
  xbar_ch + 200,
  length.out = 1e3
)

## Generate common Monte Carlo samples used across all methods
set.seed(123)

R <- 1e6

z_C <- rnorm(R)
z_T <- rnorm(R)

## Identify values within the prespecified discrepancy region
## Delta = theta_C - xbar_ch in [DEL_l, DEL_u]
idx_ps <- mu_grid >= xbar_ch + DEL_l &
  mu_grid <= xbar_ch + DEL_u


#---- FPP ----

## Compute the Type I error rate over the grid of true current control-arm means
TIE_fpp <- sapply(
  mu_grid,
  calc_TIE_power_pp,
  theta_true = theta0,
  xbar_ch = xbar_ch,
  se_ch = se_ch,
  n_T = n_T,
  n_C = n_C,
  sigma2 = sigma2,
  B_thr = p_0_fpp,
  alpha_0 = alpha_0
)

## Compute power over the grid of true current control-arm means
Power_fpp <- sapply(
  mu_grid,
  calc_TIE_power_pp,
  theta_true = -theta1,
  xbar_ch = xbar_ch,
  se_ch = se_ch,
  n_T = n_T,
  n_C = n_C,
  sigma2 = sigma2,
  B_thr = p_0_fpp,
  alpha_0 = alpha_0
)

## Summarise the operating characteristics over the full range
max_TIE_fpp <- max(TIE_fpp)
max_pow_fpp <- max(Power_fpp)
min_pow_fpp <- min(Power_fpp)

## Summarise the operating characteristics over the prespecified discrepancy region
max_TIE_fpp_ps <- max(TIE_fpp[idx_ps])
max_pow_fpp_ps <- max(Power_fpp[idx_ps])
min_pow_fpp_ps <- min(Power_fpp[idx_ps])

## Compute the power gain relative to the no-borrowing analysis
## over the prespecified discrepancy region
power_gain_fpp_ps <- Power_fpp - (1-beta)

#---- APP ----

## Compute the Type I error rate over the grid of true current control-arm means
TIE_app <- sapply(
  mu_grid,
  function(muC) {
    
    estimate_OC_cpp_box_p_fast_1(
      muC0 = muC,
      n_T = n_T,
      n_C = n_C,
      theta = theta0,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      B_thr = p_0_app,
      delta0 = delta0,
      z_C = z_C,
      z_T = z_T,
      theta0 = theta0
    )
  }
)

## Compute power over the grid of true current control-arm means
Power_app <- sapply(
  mu_grid,
  function(muC) {
    
    estimate_OC_cpp_box_p_fast_1(
      muC0 = muC,
      n_T = n_T,
      n_C = n_C,
      theta = -theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      B_thr = p_0_app,
      delta0 = delta0,
      z_C = z_C,
      z_T = z_T,
      theta0 = theta0
    )
  }
)

## Summarise the operating characteristics over the full range
max_TIE_app <- max(TIE_app)
max_pow_app <- max(Power_app)
min_pow_app <- min(Power_app)

## Summarise the operating characteristics over the prespecified discrepancy region
max_TIE_app_ps <- max(TIE_app[idx_ps])
max_pow_app_ps <- max(Power_app[idx_ps])
min_pow_app_ps <- min(Power_app[idx_ps])

## Compute the power gain relative to the no-borrowing analysis
## over the prespecified discrepancy region
power_gain_app_ps <- Power_app - (1-beta)

#---- MAP ----

## Set up parallel computation
cl <- parallel::makeCluster(ncore_set)

## Export the required objects and functions to the parallel workers
parallel::clusterExport(
  cl,
  varlist = c(
    "mu_grid",
    "n_T",
    "n_C",
    "theta0",
    "theta1",
    "sigma2",
    "xbar_ch",
    "se_ch",
    "p_0_map",
    "make_tau_grid",
    "psi_cal_map",
    "ppos_map_vec",
    "z_C",
    "z_T",
    "estimate_OC_map_p_fast_1",
    "ppos_map_vec_batch"
  ),
  envir = .GlobalEnv
)

## Compute the Type I error rate and power over the grid
## of true current control-arm means
out_map <- parallel::parLapply(
  cl,
  X = seq_along(mu_grid),
  fun = function(imu) {
    
    ## Compute the Type I error rate
    tie <- estimate_OC_map_p_fast_1(
      muC0 = mu_grid[imu],
      n_T = n_T,
      n_C = n_C,
      theta = theta0,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      B_thr = p_0_map,
      psi = psi_cal_map,
      z_C = z_C,
      z_T = z_T,
      theta0 = theta0,
      n_tau = 401
    )
    
    ## Compute power
    power <- estimate_OC_map_p_fast_1(
      muC0 = mu_grid[imu],
      n_T = n_T,
      n_C = n_C,
      theta = -theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      B_thr = p_0_map,
      psi = psi_cal_map,
      z_C = z_C,
      z_T = z_T,
      theta0 = theta0,
      n_tau = 401
    )
    
    c(
      TIE = tie,
      Power = power
    )
  }
)

## Stop the parallel cluster
parallel::stopCluster(cl)

## Combine the parallel results
out_map <- do.call(
  rbind,
  out_map
)

## Extract the Type I error rate and power
TIE_map <- out_map[, "TIE"]
Power_map <- out_map[, "Power"]

## Summarise the operating characteristics over the full range
max_TIE_map <- max(TIE_map)
max_pow_map <- max(Power_map)
min_pow_map <- min(Power_map)

## Summarise the operating characteristics over the prespecified discrepancy region
max_TIE_map_ps <- max(TIE_map[idx_ps])
max_pow_map_ps <- max(Power_map[idx_ps])
min_pow_map_ps <- min(Power_map[idx_ps])

## Compute the calibrated power gain relative to the no-borrowing analysis
## over the prespecified discrepancy region
power_gain_map_ps <- Power_map - (1-beta)

#---- rMAP ----

## Define the fixed rMAP quantities
vague_sd <- 88
psi_rmap_2 <- 0.35 * sigma

## Set up parallel computation
cl <- parallel::makeCluster(ncore_set)

## Export the required objects and functions to the parallel workers
parallel::clusterExport(
  cl,
  varlist = c(
    "estimate_OC_rmap_p_fast_1",
    "make_tau_grid",
    "mu_grid",
    "n_T",
    "n_C",
    "theta0",
    "theta1",
    "sigma2",
    "xbar_ch",
    "se_ch",
    "p_0_rmap_1",
    "p_0_rmap_2",
    "psi_cal_rmap",
    "psi_rmap_2",
    "ppos_rmap_vec_batch",
    "ppos_rmap_vec",
    "w_rob",
    "w_cal",
    "vague_sd",
    "z_C",
    "z_T"
  ),
  envir = .GlobalEnv
)

## Compute the Type I error rate and power for rMAP_1 and rMAP_2
## over the grid of true current control-arm means
out_rmap <- parallel::parLapply(
  cl,
  X = seq_along(mu_grid),
  fun = function(i) {
    
    ## Compute the Type I error rate for rMAP_1
    tie_rmap_1 <- estimate_OC_rmap_p_fast_1(
      muC0 = mu_grid[i],
      n_T = n_T,
      n_C = n_C,
      theta = theta0,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      B_thr = p_0_rmap_1,
      psi = psi_cal_rmap,
      z_C = z_C,
      z_T = z_T,
      theta0 = theta0,
      w_rob = w_rob,
      vague_sd = vague_sd,
      n_tau = 401
    )
    
    ## Compute power for rMAP_1
    power_rmap_1 <- estimate_OC_rmap_p_fast_1(
      muC0 = mu_grid[i],
      n_T = n_T,
      n_C = n_C,
      theta = -theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      B_thr = p_0_rmap_1,
      psi = psi_cal_rmap,
      z_C = z_C,
      z_T = z_T,
      theta0 = theta0,
      w_rob = w_rob,
      vague_sd = vague_sd,
      n_tau = 401
    )
    
    ## Compute the Type I error rate for rMAP_2
    tie_rmap_2 <- estimate_OC_rmap_p_fast_1(
      muC0 = mu_grid[i],
      n_T = n_T,
      n_C = n_C,
      theta = theta0,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      B_thr = p_0_rmap_2,
      psi = psi_rmap_2,
      z_C = z_C,
      z_T = z_T,
      theta0 = theta0,
      w_rob = w_cal,
      vague_sd = vague_sd,
      n_tau = 401
    )
    
    ## Compute power for rMAP_2
    power_rmap_2 <- estimate_OC_rmap_p_fast_1(
      muC0 = mu_grid[i],
      n_T = n_T,
      n_C = n_C,
      theta = -theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      B_thr = p_0_rmap_2,
      psi = psi_rmap_2,
      z_C = z_C,
      z_T = z_T,
      theta0 = theta0,
      w_rob = w_cal,
      vague_sd = vague_sd,
      n_tau = 401
    )
    
    ## Return the operating characteristics for both rMAP specifications
    c(
      TIE_rmap_1 = tie_rmap_1,
      Power_rmap_1 = power_rmap_1,
      TIE_rmap_2 = tie_rmap_2,
      Power_rmap_2 = power_rmap_2
    )
  }
)

## Stop the parallel cluster
parallel::stopCluster(cl)

## Combine the parallel results
out_rmap <- do.call(
  rbind,
  out_rmap
)

#---- rMAP_1 ----

## Extract the Type I error rate and power
TIE_rmap_1 <- out_rmap[, "TIE_rmap_1"]
Power_rmap_1 <- out_rmap[, "Power_rmap_1"]

## Summarise the operating characteristics over the full range
max_TIE_rmap_1 <- max(TIE_rmap_1)
max_pow_rmap_1 <- max(Power_rmap_1)
min_pow_rmap_1 <- min(Power_rmap_1)

## Summarise the operating characteristics over the prespecified discrepancy region
max_TIE_rmap_1_ps <- max(TIE_rmap_1[idx_ps])
max_pow_rmap_1_ps <- max(Power_rmap_1[idx_ps])
min_pow_rmap_1_ps <- min(Power_rmap_1[idx_ps])

## Compute the calibrated power gain relative to the no-borrowing analysis
## using the maximum Type I error rate over the prespecified discrepancy region
power_gain_rmap_1_ps <- Power_rmap_1 - (1-beta)

#---- rMAP_2 ----

## Extract the Type I error rate and power
TIE_rmap_2 <- out_rmap[, "TIE_rmap_2"]
Power_rmap_2 <- out_rmap[, "Power_rmap_2"]

## Summarise the operating characteristics over the full range
max_TIE_rmap_2 <- max(TIE_rmap_2)
max_pow_rmap_2 <- max(Power_rmap_2)
min_pow_rmap_2 <- min(Power_rmap_2)

## Summarise the operating characteristics over the prespecified discrepancy region
max_TIE_rmap_2_ps <- max(TIE_rmap_2[idx_ps])
max_pow_rmap_2_ps <- max(Power_rmap_2[idx_ps])
min_pow_rmap_2_ps <- min(Power_rmap_2[idx_ps])

## Compute the calibrated power gain relative to the no-borrowing analysis
## using the maximum Type I error rate over the prespecified discrepancy region
power_gain_rmap_2_ps <- Power_rmap_2 - (1-beta)

#save.image(file = "Rdatas/case_study_TIE_pow_new.RData")
