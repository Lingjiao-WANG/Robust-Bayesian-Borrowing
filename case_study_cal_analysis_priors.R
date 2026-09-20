######################################
#### Real-data example
#### Borrowing from the control arm
#### Calibrated analysis priors
######################################

library(quadprog)
library(MCMCpack)
library(doParallel)
library(parallel)
library(mgcv)

## Source functions used in the real-data example
source("case_study_funs.R")

## Number of CPU cores used for parallel computation
ncore_set <- 8

#### Historical control data

## Six historical control groups
hist_dat <- data.frame(
  study = c(
    "Gastr06",
    "AIMed07",
    "NEJM07",
    "Gastr01a",
    "APhTh04",
    "Gastr01b"
  ),
  n = c(74, 166, 328, 20, 25, 58),
  y = c(-51, -49, -36, -47, -90, -54),
  se = c(10.2, 6.8, 4.9, 19.7, 17.6, 11.6)
)

## Known outcome standard deviation and variance
sigma <- 88
sigma2 <- sigma^2

## Calculate standard errors based on the known outcome standard deviation
hist_dat$se_calc <- sigma / sqrt(hist_dat$n)


#### Historical information used for borrowing

## Consider the MAP prior reported in Best (2025).
## Borrowing is based on a single historical control group.
## AIMed07 is selected because its mean is close to the centre
## of the reported MAP prior.

xbar_ch <- hist_dat$y[2]
n_ch <- hist_dat$n[2]
se_ch <- hist_dat$se[2]

#### Design of the new study

## Treatment- and control-arm sample sizes
n_C <- 20
n_T <- 40

## Nominal Type I error rate
alpha <- 0.025

## Treatment-effect value under the null hypothesis
theta0 <- 0

## Treatment-effect value under the alternative hypothesis
theta1 <- 70

## Standard error of the treatment-effect estimator without borrowing
SE <- sigma * sqrt(
  1 / n_T + 1 / n_C
)

## Critical value for the one-sided test
z_alp <- qnorm(1 - alpha)

## Power of the design without historical borrowing
power_no_borrowing <- pnorm(
  (theta0 - z_alp * SE + theta1) / SE
)

## Set the target power equal to the power of the no-borrowing design
beta <- 1 - power_no_borrowing

#### Discrepancy region

## Define Delta = theta_C - xbar_ch and consider
## Delta in [Delta_l, Delta_u] = [-sigma/4, sigma/4] = [-22, 22]
k <- 0.25

DEL_u <- k * sqrt(sigma2)
DEL_l <- -k * sqrt(sigma2)

#### Monte Carlo simulation settings

## Generate common standard normal random variables for
## simulating the current control and treatment sample means
set.seed(123)

R <- 1e6
z_C <- rnorm(R)
z_T <- rnorm(R)

## Independent seeds used for repeated Monte Carlo evaluations
seed.num <- seq(
  from = 10,
  to = 100,
  by = 10
)

## Number of repeated Monte Carlo evaluations
S <- length(seed.num)

#### Calibrated analysis priors

#---- FPP ----

## Calibrate the FPP borrowing parameter alpha_0
alpha_0 <- uniroot(
  Calibration_f,
  interval = c(0, 1),
  DEL_l    = DEL_l,
  DEL_u    = DEL_u,
  alpha    = alpha,
  beta     = beta,
  n_T      = n_T,
  n_C      = n_C,
  se_ch    = se_ch,
  xbar_ch  = xbar_ch,
  sigma2   = sigma2,
  theta0   = theta0,
  theta1   = theta1
)$root

## Compute the calibrated posterior probability threshold for the FPP
p_0_fpp <- B_thr_f(
  DEL     = DEL_u,
  alpha   = alpha,
  alpha_0 = alpha_0,
  n_T     = n_T,
  n_C     = n_C,
  se_ch   = se_ch,
  sigma2  = sigma2,
  theta_0 = theta0
)


#---- APP ----

## Calibrate the APP borrowing parameter delta0
delta0_list <- find_delta0_roots_cal_f_ad(
  DEL_l = DEL_l,
  DEL_u = DEL_u,
  n_T = n_T,
  n_C = n_C,
  alpha = alpha,
  beta = beta,
  theta0 = theta0,
  theta1 = theta1,
  sigma2 = sigma2,
  xbar_ch = xbar_ch,
  se_ch = se_ch,
  z_C = z_C,
  z_T = z_T,
  delta0_min = 0,
  delta0_max = 3,
  n_grid = 201,
  ncore0 = ncore_set,
  tol = 1e-5
)

## Use the default value delta0 = 1 if no calibrated value is found
if (is.na(delta0_list$delta0)) {
  
  delta0 <- 1
  
} else {
  
  delta0 <- delta0_list$delta0
}

## Compute the calibrated posterior probability threshold for the APP
p_0_app <- find_global_min_p_quantile(
  DEL_l = DEL_l,
  DEL_u = DEL_u,
  z_C = z_C,
  z_T = z_T,
  xbar_ch = xbar_ch,
  n_T = n_T,
  n_C = n_C,
  sigma2 = sigma2,
  se_ch = se_ch,
  delta0 = delta0,
  theta0 = theta0,
  alpha = alpha,
  grid_len = 11
)$p_0

#---- MAP ----

## Candidate values of the MAP borrowing parameter psi
psi_grid <- 50:90

## Initialise the matrix for storing the calibration residuals
## Rows correspond to Monte Carlo seeds and columns to candidate psi values
res_mat <- matrix(
  NA_real_,
  nrow = S,
  ncol = length(psi_grid)
)

rownames(res_mat) <- paste0(
  "seed_", seed.num
)

colnames(res_mat) <- paste0(
  "psi_", psi_grid
)


#====================================================================
# Evaluate the calibration function
#====================================================================
#
# For each Monte Carlo seed, the same simulated standard normal
# variables are reused across all candidate psi values, providing
# common random numbers within each seed.
#
# For each psi and seed, compute the calibration residual
#
#   r_s(psi) = estimated worst-case power - target power.
#
# The calibrated value of psi is subsequently determined from a
# smooth GAM fitted to the calibration residuals across all seeds.
#====================================================================

for (i in seq_along(seed.num)) {
  
  ## Set the Monte Carlo seed
  set.seed(seed.num[i])
  
  ## Generate Monte Carlo samples once for the current seed
  z_C <- rnorm(R)
  z_T <- rnorm(R)
  
  for (j in seq_along(psi_grid)) {
    
    psi_test <- psi_grid[j]
    
    ## Compute the worst-case power residual
    tmp <- cal_f_map(
      DEL_l = DEL_l,
      DEL_u = DEL_u,
      n_T = n_T,
      n_C = n_C,
      alpha = alpha,
      beta = beta,
      theta0 = theta0,
      theta1 = theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      psi = psi_test,
      z_C = z_C,
      z_T = z_T,
      n_tau = 401,
      batch_size = 5000,
      ncore0 = ncore_set
    )
    
    ## Store the calibration residual
    res_mat[i, j] <- tmp$res
  }
}


#====================================================================
# Convert the calibration results to long format
#====================================================================
#
# Each row contains one Monte Carlo estimate for a given combination
# of psi and Monte Carlo seed.
#====================================================================

gam_dat <- data.frame(
  psi = rep(
    psi_grid,
    each = S
  ),
  seed = rep(
    seed.num,
    times = length(psi_grid)
  ),
  res = as.vector(res_mat)
)


#====================================================================
# Fit the GAM
#====================================================================
#
# Fit the relationship between psi and the calibration residual
#
#   r(psi) = worst-case power - target power
#
# using the Monte Carlo estimates from all seeds.
#====================================================================

fit_gam <- gam(
  res ~ s(psi, k = 6),
  data = gam_dat,
  method = "REML"
)


#====================================================================
# Evaluate the fitted GAM over a fine psi grid
#====================================================================

psi_fine <- seq(
  min(psi_grid),
  max(psi_grid),
  length.out = 5000
)

pred_gam <- predict(
  fit_gam,
  newdata = data.frame(
    psi = psi_fine
  ),
  type = "response"
)


#====================================================================
# Identify zero crossings of the fitted calibration curve
#====================================================================

## Identify intervals over which the fitted GAM changes sign
cross_index <- which(
  pred_gam[-length(pred_gam)] *
    pred_gam[-1] <= 0
)

## Stop if no zero crossing is found within the candidate range
if (length(cross_index) == 0) {
  stop(
    "No zero crossing was found within the specified psi range."
  )
}

## Refine all zero crossings using uniroot
psi_roots <- sapply(
  cross_index,
  function(k) {
    
    uniroot(
      function(x) {
        
        predict(
          fit_gam,
          newdata = data.frame(
            psi = x
          ),
          type = "response"
        )
      },
      interval = c(
        psi_fine[k],
        psi_fine[k + 1L]
      )
    )$root
  }
)

## Remove numerically duplicated roots, if any
psi_roots <- unique(
  round(
    psi_roots,
    digits = 8
  )
)

## Select the largest zero crossing as the calibrated psi
psi_selected <- min(psi_roots)

## Calibrate the borrowing parameter psi
## Extract the calibrated borrowing parameter
psi_cal_map <- psi_selected
#psi_cal_map <- 73.87874

set.seed(123)
R <- 1e6
z_C <- rnorm(R)
z_T <- rnorm(R)

## Compute the calibrated posterior probability threshold
p_0_map <- find_global_min_p_quantile_map(
  DEL_l = DEL_l,
  DEL_u = DEL_u,
  z_C = z_C,
  z_T = z_T,
  xbar_ch = xbar_ch,
  n_T = n_T,
  n_C = n_C,
  sigma2 = sigma2,
  se_ch = se_ch,
  psi = psi_cal_map,
  theta0 = theta0,
  alpha = alpha,
  grid_len = 21,
  n_tau = 401,
  ncore = ncore_set
)$p_0

#---- rMAP ----

## Standard deviation of the vague component
vague_sd <- 88

## rMAP_1: fix w = 0.8 and calibrate psi
w_rob <- 0.8

#====================================================================
# rMAP_1: Calibration settings
#====================================================================

## Candidate values of the rMAP borrowing parameter psi
psi_grid <- 40:90

## Initialise the matrix for storing the calibration residuals
## Rows correspond to Monte Carlo seeds and columns to candidate psi values
res_mat <- matrix(
  NA_real_,
  nrow = S,
  ncol = length(psi_grid)
)

rownames(res_mat) <- paste0(
  "seed_", seed.num
)

colnames(res_mat) <- paste0(
  "psi_", psi_grid
)


#====================================================================
# Evaluate the calibration function
#====================================================================
#
# For each Monte Carlo seed, the same simulated standard normal
# variables are reused across all candidate psi values, providing
# common random numbers within each seed.
#
# For each psi and seed, compute the calibration residual
#
#   r_s(psi) = estimated worst-case power - target power.
#
# The calibrated value of psi is subsequently determined from a
# smooth GAM fitted to the calibration residuals across all seeds.
#====================================================================

for (i in seq_along(seed.num)) {
  
  ## Set the Monte Carlo seed
  set.seed(seed.num[i])
  
  ## Generate Monte Carlo samples once for the current seed
  z_C <- rnorm(R)
  z_T <- rnorm(R)
  
  for (j in seq_along(psi_grid)) {
    
    psi_test <- psi_grid[j]
    
    ## Compute the worst-case power residual
    tmp <- cal_f_rmap(
      DEL_l = DEL_l,
      DEL_u = DEL_u,
      n_T = n_T,
      n_C = n_C,
      alpha = alpha,
      beta = beta,
      theta0 = theta0,
      theta1 = theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      psi = psi_test,
      z_C = z_C,
      z_T = z_T,
      w_rob = w_rob,
      vague_sd = vague_sd,
      n_tau = 401,
      batch_size = 5000,
      ncore0 = ncore_set
    )
    
    ## Store the calibration residual
    res_mat[i, j] <- tmp$res
  }
}


#====================================================================
# Convert the calibration results to long format
#====================================================================
#
# Each row contains one Monte Carlo estimate for a given combination
# of psi and Monte Carlo seed.
#====================================================================

gam_dat_rmap <- data.frame(
  psi = rep(
    psi_grid,
    each = S
  ),
  seed = rep(
    seed.num,
    times = length(psi_grid)
  ),
  res = as.vector(res_mat)
)


#====================================================================
# Fit the GAM
#====================================================================
#
# Fit the relationship between psi and the calibration residual
#
#   r(psi) = worst-case power - target power
#
# using the Monte Carlo estimates from all seeds.
#====================================================================

fit_gam_rmap <- gam(
  res ~ s(psi, k = 6),
  data = gam_dat_rmap,
  method = "REML"
)


#====================================================================
# Evaluate the fitted GAM over a fine psi grid
#====================================================================

psi_fine <- seq(
  min(psi_grid),
  max(psi_grid),
  length.out = 5000
)

pred_gam <- predict(
  fit_gam_rmap,
  newdata = data.frame(
    psi = psi_fine
  ),
  type = "response"
)


#====================================================================
# Identify zero crossings of the fitted calibration curve
#====================================================================

## Identify intervals over which the fitted GAM changes sign
cross_index <- which(
  pred_gam[-length(pred_gam)] *
    pred_gam[-1] <= 0
)

## Stop if no zero crossing is found within the candidate range
if (length(cross_index) == 0) {
  stop(
    "No zero crossing was found within the specified psi range."
  )
}

## Refine all zero crossings using uniroot
psi_roots <- sapply(
  cross_index,
  function(k) {
    
    uniroot(
      function(x) {
        
        predict(
          fit_gam_rmap,
          newdata = data.frame(
            psi = x
          ),
          type = "response"
        )
      },
      interval = c(
        psi_fine[k],
        psi_fine[k + 1L]
      )
    )$root
  }
)

## Remove numerically duplicated roots, if any
psi_roots <- unique(
  round(
    psi_roots,
    digits = 8
  )
)

## Select the largest zero crossing as the calibrated psi
psi_selected_rmap <- min(psi_roots)

## Store calibrated value
psi_cal_rmap <- as.numeric(
  psi_selected_rmap
)
#psi_cal_rmap <- 54.63299

#====================================================================
# Compute the calibrated posterior probability threshold
#====================================================================
set.seed(123)

z_C <- rnorm(R)
z_T <- rnorm(R)

p_0_rmap_1 <- find_global_min_p_quantile_rmap(
  DEL_l = DEL_l,
  DEL_u = DEL_u,
  z_C = z_C,
  z_T = z_T,
  xbar_ch = xbar_ch,
  n_T = n_T,
  n_C = n_C,
  sigma2 = sigma2,
  se_ch = se_ch,
  psi = psi_cal_rmap,
  theta0 = theta0,
  alpha = alpha,
  w_rob = w_rob,
  vague_sd = vague_sd,
  grid_len = 21,
  n_tau = 401,
  ncore = ncore_set
)$p_0

#====================================================================
# rMAP_2: Fix psi = 0.35 * sigma and calibrate w
#====================================================================

## Fix the scale parameter of the half-normal prior for tau
psi_rmap_2 <- 0.35 * sigma

## Candidate values of the MAP mixture weight w
w_grid <- seq(
  from = 0.40,
  to = 0.60,
  by = 0.01
)

## Initialise the matrix for storing the calibration residuals
## Rows correspond to Monte Carlo seeds and columns to candidate w values
res_mat_w <- matrix(
  NA_real_,
  nrow = S,
  ncol = length(w_grid)
)

rownames(res_mat_w) <- paste0(
  "seed_",
  seed.num
)

colnames(res_mat_w) <- paste0(
  "w_",
  w_grid
)


#====================================================================
# Evaluate the calibration function
#====================================================================
#
# For each Monte Carlo seed, the same simulated standard normal
# variables are reused across all candidate w values, providing
# common random numbers within each seed.
#
# For each w and seed, compute the calibration residual
#
#   r_s(w) = estimated worst-case power - target power.
#
# The calibrated value of w is subsequently determined from a
# smooth GAM fitted to the calibration residuals across all seeds.
#====================================================================

for (i in seq_along(seed.num)) {
  
  ## Set the Monte Carlo seed
  set.seed(seed.num[i])
  
  ## Generate Monte Carlo samples once for the current seed
  z_C <- rnorm(R)
  z_T <- rnorm(R)
  
  for (j in seq_along(w_grid)) {
    
    w_test <- w_grid[j]
    
    ## Compute the worst-case power residual
    tmp <- cal_f_rmap(
      DEL_l = DEL_l,
      DEL_u = DEL_u,
      n_T = n_T,
      n_C = n_C,
      alpha = alpha,
      beta = beta,
      theta0 = theta0,
      theta1 = theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      psi = psi_rmap_2,
      z_C = z_C,
      z_T = z_T,
      w_rob = w_test,
      vague_sd = vague_sd,
      n_tau = 401,
      batch_size = 5000,
      ncore0 = ncore_set
    )
    
    ## Store the calibration residual
    res_mat_w[i, j] <- tmp$res
  }
}


#====================================================================
# Convert the calibration results to long format
#====================================================================
#
# Each row contains one Monte Carlo estimate for a given combination
# of w and Monte Carlo seed.
#====================================================================

gam_dat_rmap_2 <- data.frame(
  w = rep(
    w_grid,
    each = S
  ),
  seed = rep(
    seed.num,
    times = length(w_grid)
  ),
  res = as.vector(
    res_mat_w
  )
)


#====================================================================
# Fit the GAM
#====================================================================
#
# Fit the relationship between w and the calibration residual
#
#   r(w) = worst-case power - target power
#
# using the Monte Carlo estimates from all seeds.
#====================================================================

fit_gam_rmap_2 <- gam(
  res ~ s(w, k = 6),
  data = gam_dat_rmap_2,
  method = "REML"
)


#====================================================================
# Evaluate the fitted GAM over a fine w grid
#====================================================================

w_fine <- seq(
  min(w_grid),
  max(w_grid),
  length.out = 5000
)

pred_gam_w <- predict(
  fit_gam_rmap_2,
  newdata = data.frame(
    w = w_fine
  ),
  type = "response"
)


#====================================================================
# Identify zero crossings of the fitted calibration curve
#====================================================================

## Identify intervals over which the fitted GAM changes sign
cross_index_w <- which(
  pred_gam_w[-length(pred_gam_w)] *
    pred_gam_w[-1] <= 0
)

## Stop if no zero crossing is found within the candidate range
if (length(cross_index_w) == 0) {
  stop(
    "No zero crossing was found within the specified w range."
  )
}

## Refine all zero crossings using uniroot
w_roots <- sapply(
  cross_index_w,
  function(k) {
    
    uniroot(
      function(x) {
        
        predict(
          fit_gam_rmap_2,
          newdata = data.frame(
            w = x
          ),
          type = "response"
        )
      },
      interval = c(
        w_fine[k],
        w_fine[k + 1L]
      )
    )$root
  }
)

## Remove numerically duplicated roots, if any
w_roots <- unique(
  round(
    w_roots,
    digits = 8
  )
)

## Select the smallest zero crossing as the calibrated MAP mixture weight
w_selected_rmap <- max(w_roots)

## Store calibrated value
w_cal <- as.numeric(
  w_selected_rmap
)
#w_cal <- 0.4841975 

set.seed(123)

z_C <- rnorm(R)
z_T <- rnorm(R)

## Compute the calibrated posterior probability threshold
p_0_rmap_2 <- find_global_min_p_quantile_rmap(
  DEL_l = DEL_l,
  DEL_u = DEL_u,
  z_C = z_C,
  z_T = z_T,
  xbar_ch = xbar_ch,
  n_T = n_T,
  n_C = n_C,
  sigma2 = sigma2,
  se_ch = se_ch,
  psi = psi_rmap_2,
  theta0 = theta0,
  alpha = alpha,
  w_rob = w_cal,
  vague_sd = 88,
  grid_len = 21,
  n_tau = 401,
  ncore = ncore_set
)$p_0

#save.image(file = "Rdatas/Calibrated_analysis_priors.RData")