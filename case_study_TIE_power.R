######################################
#### Real data example:           ####
#### Borrow from the control arm  ####
#### TIE, power and calibrated    #### 
#### power gain                   ####
######################################
library(quadprog)
library(MCMCpack)
library(doParallel)
library(parallel)

## Source functions
source("case_study_funs.R")

## Source the data of calibrated analysis priors
load("Calibrated_analysis_priors.RData")

# Standard error for no-borrowing case
SE_nob <- sqrt(sigma2 / n_T + sigma2 / n_C)

mu_grid <- seq(
  xbar_ch - 200,
  xbar_ch + 200,
  length.out = 1e5
)

# --- FPP ----
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

max_TIE_fpp=max(TIE_fpp)
max_pow_fpp=max(Power_fpp)
min_pow_fpp=min(Power_fpp)

max_TIE_fpp_ps=max(TIE_fpp[which(mu_grid>= xbar_ch + DEL_l & mu_grid<= xbar_ch + DEL_u)])
max_pow_fpp_ps=max(Power_fpp[which(mu_grid>= xbar_ch + DEL_l & mu_grid<= xbar_ch + DEL_u)])
min_pow_fpp_ps=min(Power_fpp[which(mu_grid>= xbar_ch + DEL_l & mu_grid<= xbar_ch + DEL_u)])

pow_nob_cal_fpp <- pnorm(
  (theta1 / SE_nob) - qnorm(1 - max_TIE_fpp)
)

pow_nob_cal_fpp_ps <- pnorm(
  (theta1 / SE_nob) - qnorm(1 - max_TIE_fpp_ps)
)

power_gain_fpp=Power_fpp-pow_nob_cal_fpp
power_gain_fpp_ps=Power_fpp-pow_nob_cal_fpp_ps

# --- APP ----
# APP_ps

TIE_app_ps <- sapply(mu_grid, function(muC) {
  estimate_OC_cpp_box_p_fast_1(
    muC0   = muC,
    n_T    = n_T,
    n_C    = n_C,
    theta  = theta0,
    sigma2 = sigma2,
    xbar_ch = xbar_ch,
    se_ch  = se_ch,
    B_thr  = p_0_app_ps,
    delta0 = delta0,
    z_C    = z_C,
    z_T    = z_T,
    theta0 = theta0
  )
})

Power_app_ps <- sapply(mu_grid, function(muC) {
  estimate_OC_cpp_box_p_fast_1(
    muC0   = muC,
    n_T    = n_T,
    n_C    = n_C,
    theta  = -theta1,
    sigma2 = sigma2,
    xbar_ch = xbar_ch,
    se_ch  = se_ch,
    B_thr  = p_0_app_ps,
    delta0 = delta0,
    z_C    = z_C,
    z_T    = z_T,
    theta0 = theta0
  )
})

max_TIE_app_ps=max(TIE_app_ps)
max_pow_app_ps=max(Power_app_ps)
min_pow_app_ps=min(Power_app_ps)

max_TIE_app_ps_ps=max(TIE_app_ps[which(mu_grid>= xbar_ch + DEL_l & mu_grid<= xbar_ch + DEL_u)])
max_pow_app_ps_ps=max(Power_app_ps[which(mu_grid>= xbar_ch + DEL_l & mu_grid<= xbar_ch + DEL_u)])
min_pow_app_ps_ps=min(Power_app_ps[which(mu_grid>= xbar_ch + DEL_l & mu_grid<= xbar_ch + DEL_u)])

pow_nob_cal_app_ps <- pnorm(
  (theta1 / SE_nob) - qnorm(1 - max_TIE_app_ps)
)

pow_nob_cal_app_ps_ps <- pnorm(
  (theta1 / SE_nob) - qnorm(1 - max_TIE_app_ps_ps)
)

power_gain_app_ps=Power_app_ps-pow_nob_cal_app_ps
power_gain_app_ps_ps=Power_app_ps-pow_nob_cal_app_ps_ps

# APP_ss
TIE_app_ss <- sapply(mu_grid, function(muC) {
  estimate_OC_cpp_box_p_fast_1(
    muC0   = muC,
    n_T    = n_T,
    n_C    = n_C,
    theta  = theta0,
    sigma2 = sigma2,
    xbar_ch = xbar_ch,
    se_ch  = se_ch,
    B_thr  = p_0_app_ss,
    delta0 = delta0,
    z_C    = z_C,
    z_T    = z_T,
    theta0 = theta0
  )
})

Power_app_ss <- sapply(mu_grid, function(muC) {
  estimate_OC_cpp_box_p_fast_1(
    muC0   = muC,
    n_T    = n_T,
    n_C    = n_C,
    theta  = -theta1,
    sigma2 = sigma2,
    xbar_ch = xbar_ch,
    se_ch  = se_ch,
    B_thr  = p_0_app_ss,
    delta0 = delta0,
    z_C    = z_C,
    z_T    = z_T,
    theta0 = theta0
  )
})

max_TIE_app_ss=max(TIE_app_ss)
max_pow_app_ss=max(Power_app_ss)
min_pow_app_ss=min(Power_app_ss)

max_TIE_app_ss_ss=max(TIE_app_ss[which(mu_grid>= xbar_ch - ss_asy_list_within_ps$DEL_u & mu_grid<= xbar_ch - ss_asy_list_within_ps$DEL_l)])
max_pow_app_ss_ss=max(Power_app_ss[which(mu_grid>= xbar_ch - ss_asy_list_within_ps$DEL_u & mu_grid<= xbar_ch - ss_asy_list_within_ps$DEL_l)])
min_pow_app_ss_ss=min(Power_app_ss[which(mu_grid>= xbar_ch - ss_asy_list_within_ps$DEL_u & mu_grid<= xbar_ch - ss_asy_list_within_ps$DEL_l)])

pow_nob_cal_app_ss <- pnorm(
  (theta1 / SE_nob) - qnorm(1 - max_TIE_app_ss)
)

pow_nob_cal_app_ss_ss <- pnorm(
  (theta1 / SE_nob) - qnorm(1 - max_TIE_app_ss_ss)
)

power_gain_app_ss=Power_app_ss-pow_nob_cal_app_ss
power_gain_app_ss_ss=Power_app_ss-pow_nob_cal_app_ss_ss

mu_grid <- seq(
  xbar_ch - 200,
  xbar_ch + 200,
  length.out = 1e3
)

# --- MAP ----

cl <- makeCluster(8)

clusterExport(
  cl,
  varlist = c(
    "estimate_OC_map_p_fast_1",
    "make_tau_grid",
    "mu_grid",
    "n_T",
    "n_C",
    "theta0",
    "theta1",
    "sigma2",
    "xbar_ch",
    "se_ch",
    "p_0_map",
    "psi_cal_map",
    "z_C",
    "z_T"
  ),
  envir = environment()
)

out_map <- parLapply(
  cl,
  X = seq_along(mu_grid),
  fun = function(imu) {
    
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
    
    c(TIE = tie, Power = power)
  }
)

stopCluster(cl)

out_map <- do.call(rbind, out_map)

TIE_map <- out_map[, "TIE"]
Power_map <- out_map[, "Power"]

max_TIE_map=max(TIE_map)
max_pow_map=max(Power_map)
min_pow_map=min(Power_map)

max_TIE_map_ps=max(TIE_map[which(mu_grid>= xbar_ch + DEL_l & mu_grid<= xbar_ch + DEL_u)])
max_pow_map_ps=max(Power_map[which(mu_grid>= xbar_ch + DEL_l & mu_grid<= xbar_ch + DEL_u)])
min_pow_map_ps=min(Power_map[which(mu_grid>= xbar_ch + DEL_l & mu_grid<= xbar_ch + DEL_u)])

pow_nob_cal_map <- pnorm(
  (theta1 / SE_nob) - qnorm(1 - max_TIE_map)
)

pow_nob_cal_map_ps <- pnorm(
  (theta1 / SE_nob) - qnorm(1 - max_TIE_map_ps)
)

power_gain_map=Power_map-pow_nob_cal_map
power_gain_map_ps=Power_map-pow_nob_cal_map_ps

# --- rMAP ----

cl <- makeCluster(8)

clusterExport(
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
    "w_rob",
    "w_cal",
    "z_C",
    "z_T"
  ),
  envir = environment()
)

out_rmap <- parLapply(
  cl,
  X = seq_along(mu_grid),
  fun = function(i) {
    
    ## rMAP_1: Type I error
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
      vague_sd = 88,
      n_tau = 401
    )
    
    ## rMAP_1: Power
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
      vague_sd = 88,
      n_tau = 401
    )
    
    ## rMAP_2: Type I error
    tie_rmap_2 <- estimate_OC_rmap_p_fast_1(
      muC0 = mu_grid[i],
      n_T = n_T,
      n_C = n_C,
      theta = theta0,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      B_thr = p_0_rmap_2,
      psi = 0.35 * sigma,
      z_C = z_C,
      z_T = z_T,
      theta0 = theta0,
      w_rob = w_cal,
      vague_sd = 88,
      n_tau = 401
    )
    
    ## rMAP_2: Power
    power_rmap_2 <- estimate_OC_rmap_p_fast_1(
      muC0 = mu_grid[i],
      n_T = n_T,
      n_C = n_C,
      theta = -theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      B_thr = p_0_rmap_2,
      psi = 0.35 * sigma,
      z_C = z_C,
      z_T = z_T,
      theta0 = theta0,
      w_rob = w_cal,
      vague_sd = 88,
      n_tau = 401
    )
    
    c(
      TIE_rmap_1 = tie_rmap_1,
      Power_rmap_1 = power_rmap_1,
      TIE_rmap_2 = tie_rmap_2,
      Power_rmap_2 = power_rmap_2
    )
  }
)

stopCluster(cl)

## Convert the output list to a matrix
out_rmap <- do.call(rbind, out_rmap)

# rMAP_1
TIE_rmap_1   <- out_rmap[, "TIE_rmap_1"]
Power_rmap_1 <- out_rmap[, "Power_rmap_1"]

max_TIE_rmap_1=max(TIE_rmap_1)
max_pow_rmap_1=max(Power_rmap_1)
min_pow_rmap_1=min(Power_rmap_1)

max_TIE_rmap_ps_1=max(TIE_rmap_1[which(mu_grid>= xbar_ch + DEL_l & mu_grid<= xbar_ch + DEL_u)])
max_pow_rmap_ps_1=max(Power_rmap_1[which(mu_grid>= xbar_ch + DEL_l & mu_grid<= xbar_ch + DEL_u)])
min_pow_rmap_ps_1=min(Power_rmap_1[which(mu_grid>= xbar_ch + DEL_l & mu_grid<= xbar_ch + DEL_u)])

pow_nob_cal_rmap_1 <- pnorm(
  (theta1 / SE_nob) - qnorm(1 - max_TIE_rmap_1)
)

pow_nob_cal_rmap_1_ps <- pnorm(
  (theta1 / SE_nob) - qnorm(1 - max_TIE_rmap_ps_1)
)

power_gain_rmap_1=Power_rmap_1-pow_nob_cal_rmap_1
power_gain_rmap_1_ps=Power_rmap_1-pow_nob_cal_rmap_1_ps


# rMAP_2
TIE_rmap_2   <- out_rmap[, "TIE_rmap_2"]
Power_rmap_2 <- out_rmap[, "Power_rmap_2"]

max_TIE_rmap_2=max(TIE_rmap_2)
max_pow_rmap_2=max(Power_rmap_2)
min_pow_rmap_2=min(Power_rmap_2)

max_TIE_rmap_ps_2=max(TIE_rmap_2[which(mu_grid>= xbar_ch + DEL_l & mu_grid<= xbar_ch + DEL_u)])
max_pow_rmap_ps_2=max(Power_rmap_2[which(mu_grid>= xbar_ch + DEL_l & mu_grid<= xbar_ch + DEL_u)])
min_pow_rmap_ps_2=min(Power_rmap_2[which(mu_grid>= xbar_ch + DEL_l & mu_grid<= xbar_ch + DEL_u)])

pow_nob_cal_rmap_2 <- pnorm(
  (theta1 / SE_nob) - qnorm(1 - max_TIE_rmap_2)
)

pow_nob_cal_rmap_2_ps <- pnorm(
  (theta1 / SE_nob) - qnorm(1 - max_TIE_rmap_ps_2)
)

power_gain_rmap_2=Power_rmap_2-pow_nob_cal_rmap_2
power_gain_rmap_2_ps=Power_rmap_2-pow_nob_cal_rmap_2_ps

#### Print Table 2 ####

tab_fpp <- data.frame(
  Prior        = "FPP",
  max_TIE      = round(max_TIE_fpp,3),
  max_power    = round(max(power_gain_fpp),3),
  min_power    = round(min(power_gain_fpp),3),
  max_power_ps = round(max(power_gain_fpp_ps),3)
)

tab_app_ps <- data.frame(
  Prior        = "APP_ps",
  max_TIE      = round(max_TIE_app_ps, 3),
  max_power    = round(max(power_gain_app_ps), 3),
  min_power    = round(min(power_gain_app_ps), 3),
  max_power_ps = round(max(power_gain_app_ps_ps), 3)
)

tab_app_ss <- data.frame(
  Prior        = "APP_ss",
  max_TIE      = round(max_TIE_app_ss, 3),
  max_power    = round(max(power_gain_app_ss), 3),
  min_power    = round(min(power_gain_app_ss), 3),
  max_power_ps = round(max(power_gain_app_ss_ss), 3)
)

tab_map <- data.frame(
  Prior        = "MAP",
  max_TIE      = round(max_TIE_map, 3),
  max_power    = round(max(power_gain_map), 3),
  min_power    = round(min(power_gain_map), 3),
  max_power_ps = round(max(power_gain_map_ps), 3)
)

tab_rmap_1 <- data.frame(
  Prior        = "rMAP_1",
  max_TIE      = round(max_TIE_rmap_1, 3),
  max_power    = round(max(power_gain_rmap_1), 3),
  min_power    = round(min(power_gain_rmap_1), 3),
  max_power_ps = round(max(power_gain_rmap_1_ps), 3)
)

tab_rmap_2 <- data.frame(
  Prior        = "rMAP_2",
  max_TIE      = round(max_TIE_rmap_2, 3),
  max_power    = round(max(power_gain_rmap_2), 3),
  min_power    = round(min(power_gain_rmap_2), 3),
  max_power_ps = round(max(power_gain_rmap_2_ps), 3)
)

## Combine all methods
Table_2 <- rbind(
  tab_fpp,
  tab_app_ps,
  tab_app_ss,
  tab_map,
  tab_rmap_1,
  tab_rmap_2
)

print(Table_2)





