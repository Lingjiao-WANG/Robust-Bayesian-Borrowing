######################################
#### Real-data example
#### Borrowing from the control arm
#### Calibrated analysis priors
######################################

library(quadprog)
library(MCMCpack)
library(doParallel)
library(parallel)

## Source functions
source("case_study_funs.R")

#### Six historical control groups

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

## Known outcome standard deviation
sigma <- 88
sigma2 <- sigma^2

## Calculate the standard errors based on the known standard deviation
hist_dat$se_calc <- sigma / sqrt(hist_dat$n)

#### Consider the MAP prior reported in Best (2025)
#### Borrowing is based on a single historical control group.
#### AIMed07 is selected because its mean is close to the centre
#### of the MAP prior.

xbar_ch <- hist_dat$y[2]
n_ch <- hist_dat$n[2]
se_ch <- hist_dat$se[2]

#### Design of the new study

## Treatment- and control-arm sample sizes
n_C <- 20
n_T <- 40

## Nominal Type I error rate
alpha <- 0.025

## Treatment effects under the null and alternative hypotheses
theta0 <- 0

## Under the convention H_0: theta <= theta0, theta1 = 70 corresponds
## to theta1 = -70 under the equivalent convention H_0: theta >= theta0
theta1 <- 70

## Standard error of the treatment-effect estimator without borrowing
SE <- sigma * sqrt(1 / n_T + 1 / n_C)

## Success rule: theta_T - theta_C < theta0
## Reject H_0 if the observed treatment difference is sufficiently small

z_alp <- qnorm(1 - alpha)

## Power of the design without historical borrowing
power_no_borrowing <- pnorm(
  (theta0 - z_alp * SE + theta1) / SE
)

## Target Type II error rate
beta <- 1 - power_no_borrowing

#### Specify the discrepancy region between the current and historical controls

## Delta in [Delta_l, Delta_u] = [-sigma/4, sigma/4] = [-22, 22]
k <- 0.25
DEL_u <- k * sqrt(sigma2)
DEL_l <- -k * sqrt(sigma2)

#### Generate standard normal values for simulating xbar_C and xbar_T
set.seed(123)
R <- 1e4
z_C <- rnorm(R)
z_T <- rnorm(R)

#### Calibrated analysis priors

#---- FPP ----

## Calibrate the borrowing parameter alpha_0
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

## Compute the calibrated posterior probability threshold p
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

## Calibrate the borrowing parameter delta0
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
  ncore0 = 8,
  tol = 1e-5
)

## Use the default value delta0 = 1 if no calibrated value is found
if (is.na(delta0_list$delta0)) {
  delta0 <- 1
} else {
  delta0 <- delta0_list$delta0
}

## APP_ps: compute the calibrated posterior probability threshold
p_0_app_ps <- find_global_min_p_quantile(
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

## APP_ss: determine the widest feasible discrepancy interval 
##         within [Delta_l, Delta_u] = [-22, 22]
ss_asy_list_within_ps <- find_widest_fast(
  DEL_l_min = DEL_l,
  DEL_l_max = DEL_u - 0.001,
  DEL_u_min = DEL_l + 0.001,
  DEL_u_max = DEL_u,
  delta0 = delta0,
  n_l = 15,
  n_u = 15,
  ncore = 8
)

## Extract the calibrated threshold and sweet-spot interval
p_0_app_ss <- ss_asy_list_within_ps$p_0
D_ss <- -c(
  ss_asy_list_within_ps$DEL_u,
  ss_asy_list_within_ps$DEL_l
)


#---- MAP ----

## Calibrate the borrowing parameter psi
psi_res <- find_all_psi(
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
  psi_grid = seq(
    0.001,
    3 * sigma,
    length.out = 30
  ),
  n_tau = 401,
  ncore0 = 8
)

## Extract the calibrated borrowing parameter
psi_cal_map <- psi_res$psi_use

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
  grid_len = 201,
  n_tau = 401,
  ncore = 8
)$p_0

#---- rMAP ----

## Standard deviation of the vague component
vague_sd <- 88

## rMAP_1: fix w = 0.8 and calibrate psi
w_rob <- 0.8

## Calibrate the borrowing parameter psi
psi_rmap_res <- find_all_psi_rmap(
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
  w_rob = w_rob,
  z_C = z_C,
  z_T = z_T,
  psi_grid = seq(
    0.001,
    3 * sigma,
    length.out = 30
  ),
  vague_sd = vague_sd,
  n_tau = 401,
  ncore0 = 8,
  tol = 1e-6
)

## Extract the calibrated value of psi
psi_cal_rmap <- psi_rmap_res$psi_use

## Compute the calibrated posterior probability threshold
p_0_rmap_1 <- find_global_min_p_quantile_rmap(
  DEL_l = DEL_l,
  DEL_u = DEL_u,
  z_C = z_C,
  z_T = z_T,
  xbar_ch = xbar_ch,
  n_T = n_T,
  n_C = n_C,
  se_ch = se_ch,
  sigma2 = sigma2,
  psi = psi_cal_rmap,
  theta0 = theta0,
  alpha = alpha,
  w_rob = w_rob,
  vague_sd = vague_sd,
  grid_len = 201,
  n_tau = 401,
  ncore = 8
)$p_0


## rMAP_2: fix psi = 0.35 * sigma and calibrate w
psi_rmap_2 <- 0.35 * sigma

## Calibrate the borrowing parameter w
w_res <- find_all_w_rmap(
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
  w_grid = seq(
    0.001,
    0.99,
    length.out = 30
  ),
  vague_sd = vague_sd,
  n_tau = 401,
  ncore0 = 8,
  tol = 1e-6
)

## Extract the calibrated value of w
w_cal <- w_res$w_use

## Compute the calibrated posterior probability threshold
p_0_rmap_2 <- find_global_min_p_quantile_rmap(
  DEL_l = DEL_l,
  DEL_u = DEL_u,
  z_C = z_C,
  z_T = z_T,
  xbar_ch = xbar_ch,
  n_T = n_T,
  n_C = n_C,
  se_ch = se_ch,
  sigma2 = sigma2,
  psi = psi_rmap_2,
  theta0 = theta0,
  alpha = alpha,
  w_rob = w_cal,
  vague_sd = vague_sd,
  grid_len = 201,
  n_tau = 401,
  ncore = 8
)$p_0

## Print the calibrated parameters in Table 1

Table_1 <- rbind(
  
  data.frame(
    Method = "FPP",
    Borrowing_parameter = paste0(
      "a0 = ", round(alpha_0, 4)
    ),
    Threshold = paste0(
      "p = ", round(p_0_fpp, 4)
    )
  ),
  
  data.frame(
    Method = "APP_ps",
    Borrowing_parameter = paste0(
      "delta0 = ", round(delta0, 4)
    ),
    Threshold = paste0(
      "p = ", round(p_0_app_ps, 4)
    )
  ),
  
  data.frame(
    Method = "APP_ss",
    Borrowing_parameter = paste0(
      "delta0 = ", round(delta0, 4),
      ", sweet spot = [",
      round(D_ss[1], 4), ", ",
      round(D_ss[2], 4), "]"
    ),
    Threshold = paste0(
      "p = ", round(p_0_app_ss, 4)
    )
  ),
  
  data.frame(
    Method = "MAP",
    Borrowing_parameter = paste0(
      "psi = ", round(psi_cal_map, 4)
    ),
    Threshold = paste0(
      "p = ", round(p_0_map, 4)
    )
  ),
  
  data.frame(
    Method = "rMAP_1",
    Borrowing_parameter = paste0(
      "w = ", round(w_rob, 4),
      ", psi = ", round(psi_cal_rmap, 4)
    ),
    Threshold = paste0(
      "p = ", round(p_0_rmap_1, 4)
    )
  ),
  
  data.frame(
    Method = "rMAP_2",
    Borrowing_parameter = paste0(
      "w = ", round(w_cal, 4),
      ", psi = ", round(0.35 * sigma, 4)
    ),
    Threshold = paste0(
      "p = ", round(p_0_rmap_2, 4)
    )
  )
)

print(Table_1, row.names = FALSE)


