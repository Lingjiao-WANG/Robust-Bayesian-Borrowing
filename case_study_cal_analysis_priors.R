######################################
#### Real data example:           ####
#### Borrow from the control arm  ####
#### Calibrated analysis priors   ####
######################################
library(quadprog)
library(MCMCpack)
library(doParallel)
library(parallel)
library(RBesT)

## Source functions
source("case_study_funs.R")

#### Six sets of historical controls ####

hist_dat <- data.frame(
  study = c("Gastr06",
            "AIMed07",
            "NEJM07",
            "Gastr01a",
            "APhTh04",
            "Gastr01b"),
  n = c(74, 166, 328, 20, 25, 58),
  y = c(-51, -49, -36, -47, -90, -54),
  se = c(10.2, 6.8, 4.9, 19.7, 17.6, 11.6)
)

# Known standard deviation
sigma <- 88
sigma2 <- sigma^2

hist_dat$se_calc <- sigma / sqrt(hist_dat$n)

#### MAP prior in Best (2025) for six historical datasets ####

map_data <- data.frame(
  study = hist_dat$study,
  mean  = hist_dat$y,
  se    = hist_dat$se
)


gMAP_fit <- gMAP(
  cbind(mean, se) ~ 1 | study,
  data = map_data,
  family = gaussian,
  tau.dist = "HalfNormal",
  tau.prior = sqrt(88^2/4)
)

fit <- gMAP(
  cbind(mean, se) ~ 1 | study,
  data      = map_data,
  family    = gaussian,
  tau.dist  = "HalfNormal",
  tau.prior = sigma / 2,
  beta.prior = cbind(0, sigma)
)

map_prior <- automixfit(gMAP_fit, Nc = 3)
# map_prior

#### Consider the MAP prior given in Best (2025) ####
#### we only borrow from one historical control  ####
#### then we select the AIMed07, which has the   ####
#### mean is close to the center of MAP.         ####

xbar_ch=hist_dat$y[2]
n_ch=hist_dat$n[2]
se_ch=hist_dat$se[2]

#### Design for the new study ####

# Grid for true control mean
mu_grid <- seq(
  xbar_ch - 200,
  xbar_ch + 200,
  length.out = 1e5
)

n_C=20
n_T=40

alpha <- 0.025

theta0 <- 0      # null treatment difference
theta1 <- 70    # true treatment difference: theta_T - theta_C

SE <- sigma * sqrt(1 / n_T + 1 / n_C)

## Success rule: theta_T - theta_C < 0
## Reject if observed difference is less than theta0 - z * SE
z_alp <- qnorm(1 - alpha)

power_no_borrowing <- pnorm((theta0 - z_alp * SE + theta1) / SE)

beta <- 1-power_no_borrowing

# Set the space of drift Delta between mu_C and ybar_Ch
## Delta \in [Delta_l, Delta_u] = [-sigma/4, sigma/4] = [-22, 22]

k=0.25 
DEL_u=k*sqrt(sigma2)
DEL_l=-k*sqrt(sigma2)

# Set the simulated values of ybar_C and ybar_T
set.seed(123)
R <- 1e4
z_C <- rnorm(R)
z_T <- rnorm(R)

#### Calibrated analysis priors ####
# --- FPP ----

## Calibrated parameters
# Tuning parameter -- a_0
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

# Posterior probability threshold -- p 
p_0_fpp <- B_thr_f(
  DEL     = DEL_u,
  alpha   = alpha,
  alpha_0 = alpha_0,
  n_T     = n_T,
  n_C     = n_C,
  se_ch   = se_ch,
  sigma2  = sigma2,
  theta_d = theta0
)

# --- APP ----

delta0_list=find_delta0_roots_cal_f_ad(DEL_l, DEL_u,
                                       n_T, n_C,
                                       alpha, beta,
                                       theta0, theta1,
                                       sigma2,
                                       xbar_ch, se_ch,
                                       z_C, z_T,
                                       delta0_min = 0,
                                       delta0_max = 3,
                                       n_grid = 201,
                                       ncore0 = 8,
                                       tol = 1e-5)

if(is.na(delta0_list$delta0)==TRUE){
  delta0=1
}else{
  delta0=delta0_list$delta0
}

## APP_ps

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

## APP_ss

# Compute sweet spot 

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

p_0_app_ss <- ss_asy_list_within_ps$p_0
D_ss <- -c(ss_asy_list_within_ps$DEL_u,ss_asy_list_within_ps$DEL_l)

# --- MAP ----

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
  psi_grid = seq(0.001, 3*sigma, length.out = 30),
  n_tau = 401,
  ncore0 = 8
)

if(is.na(psi_res$psi_use) == TRUE){
  psi_cal_map=30.8
}else{
  psi_cal_map=psi_res$psi_use
}

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
#0.0186886

# --- rMAP ----

# rMAP_1 with w=0.8
w_rob=0.8

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
  psi_grid = seq(0.001, 3 * sqrt(sigma2), length.out = 30),
  vague_sd = 88,
  n_tau = 401,
  ncore0 = 8,
  tol = 1e-6
)

if (is.na(psi_rmap_res$psi_use)) {
  psi_cal_rmap <- psi
} else {
  psi_cal_rmap <- psi_rmap_res$psi_use
}

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
  vague_sd = 88,
  grid_len = 201,
  n_tau = 401,
  ncore = 8
)$p_0
#0.01869885

# rMAP_2 with psi=0.35*sigma

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
  psi = 0.35*sigma,
  z_C = z_C,
  z_T = z_T,
  w_grid = seq(0.001, 0.99, length.out = 30),
  vague_sd = 88,
  n_tau = 401,
  ncore0 = 8,
  tol = 1e-6
)

if (is.na(w_res$w_use)) {
  w_cal <- w_rob
} else {
  w_cal <- w_res$w_use
}



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
  psi = 0.35*sigma,
  theta0 = theta0,
  alpha = alpha,
  w_rob = w_cal,
  vague_sd = 88,
  grid_len = 201,
  n_tau = 401,
  ncore = 8
)$p_0
#0.01864751

## Print the outcomes in Table 1

Table_1 <- rbind(
  data.frame(
    Method = "FPP",
    Tuning_parameter = paste0("a0 = ", round(alpha_0, 4)),
    Threshold = paste0("p = ", round(p_0_fpp, 4))
  ),
  data.frame(
    Method = "APP_ps",
    Tuning_parameter = paste0("delta0 = ", round(delta0, 4)),
    Threshold = paste0("p = ", round(p_0_app_ps, 4))
  ),
  data.frame(
    Method = "APP_ss",
    Tuning_parameter = paste0(
      "delta0 = ", round(delta0, 4),
      ", sweet spot = [",
      round(D_ss[1], 4), ", ",
      round(D_ss[2], 4), "]"
    ),
    Threshold = paste0("p = ", round(p_0_app_ss, 4))
  ),
  data.frame(
    Method = "MAP",
    Tuning_parameter = paste0("psi = ", round(psi_cal_map, 4)),
    Threshold = paste0("p = ", round(p_0_map, 4))
  ),
  data.frame(
    Method = "rMAP_1",
    Tuning_parameter = paste0(
      "w = ", round(w_rob, 4),
      ", psi = ", round(psi_cal_rmap, 4)
    ),
    Threshold = paste0("p = ", round(p_0_rmap_1, 4))
  ),
  data.frame(
    Method = "rMAP_2",
    Tuning_parameter = paste0(
      "w = ", round(w_cal, 4),
      ", psi = ", round(0.35 * sigma, 4)
    ),
    Threshold = paste0("p = ", round(p_0_rmap_2, 4))
  )
)

print(Table_1, row.names = FALSE)



