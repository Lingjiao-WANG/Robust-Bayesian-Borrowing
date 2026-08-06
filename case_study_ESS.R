######################################
#### Real data example:           ####
#### Borrow from the control arm  ####
#### Compute the EHSS             #### 
######################################

library(quadprog)
library(MCMCpack)
library(doParallel)
library(parallel)

## Source functions
source("case_study_funs.R")

## Source the data of calibrated analysis priors
load("Calibrated_analysis_priors.RData")

mu_grid <- seq(
  xbar_ch - 200,
  xbar_ch + 200,
  length.out = 1e3
)

# --- FPP ----
ESS_fpp=alpha_0*n_ch

# --- APP ----
out_ESS_app_1=ESS_adaptive_power_prior_xbar(xbar_C_grid=mu_grid,
                                            xbar_ch,
                                            n_C,
                                            sigma2,
                                            se_ch,
                                            delta0,
                                            weights = NULL)

out_ESS_app_2=ESS_adaptive_power_prior_muC(mu_grid=mu_grid,
                                           xbar_ch,
                                           n_C,
                                           sigma2,
                                           se_ch,
                                           delta0,
                                           weights = NULL)


# --- MAP ----

out_ESS_map_2 <- ESS_MAP_xbarC_given_mu(
  mu_grid = mu_grid,
  xbar_ch = xbar_ch,
  n_C = n_C,
  sigma2 = sigma2,
  se_ch = se_ch,
  psi = psi_cal_map,
  n_tau = 401
)

out_ESS_map_1 <- ESS_MAP_xbarC(
  xbar_C_grid = mu_grid,
  xbar_ch = xbar_ch,
  n_C = n_C,
  sigma2 = sigma2,
  se_ch = se_ch,
  psi = psi_cal_map,
  n_tau = 401
)


# --- rMAP ---
# rMAP_1

out_ESS_rmap_1_1=ESS_rMAP_xbarC(xbar_C_grid=mu_grid,
                                xbar_ch,
                                n_C,
                                se_ch,
                                sigma2,
                                psi=psi_cal_rmap,
                                w_rob = w_rob,
                                vague_sd = 88,
                                weights = NULL,
                                n_tau = 401)

out_ESS_rmap_1_2=ESS_rMAP_xbarC_given_mu(mu_grid = mu_grid,
                                         xbar_ch = xbar_ch,
                                         n_C = n_C,
                                         sigma2 = sigma2,
                                         se_ch = se_ch,
                                         psi=psi_cal_rmap,
                                         w_rob = w_rob,
                                         vague_sd = 88,
                                         weights = NULL,
                                         n_tau = 401,
                                         z_grid = seq(-5, 5, length.out = 1001))


# rMAP_2
out_ESS_rmap_2_1=ESS_rMAP_xbarC(xbar_C_grid=mu_grid,
                                xbar_ch,
                                n_C,
                                se_ch,
                                sigma2,
                                psi=0.35*sigma,
                                w_rob = w_cal,
                                vague_sd = 88,
                                weights = NULL,
                                n_tau = 401)

out_ESS_rmap_2_2=ESS_rMAP_xbarC_given_mu(mu_grid = mu_grid,
                                         xbar_ch = xbar_ch,
                                         n_C = n_C,
                                         sigma2 = sigma2,
                                         se_ch = se_ch,
                                         psi=0.35*sigma,
                                         w_rob = w_cal,
                                         vague_sd = 88,
                                         weights = NULL,
                                         n_tau = 401,
                                         z_grid = seq(-5, 5, length.out = 1001))

#### Plot Figure 4 ####

par(mfrow = c(1, 2))


plot(mu_grid,out_ESS_map_1$ESS_grid,type = 'l',col=2,lwd=2,lty=2,,ylim = c(min(out_ESS_rmap_2_2$ESS_grid,out_ESS_rmap_1_2$ESS_grid)-2,max(out_ESS_app_2$ESS_grid)+2),xlab = expression(bar(y)[C]),ylab = "EHSS")
lines(mu_grid,out_ESS_rmap_1_1$ESS_grid,col=3,lwd=2,lty=3)
lines(mu_grid,out_ESS_rmap_2_1$ESS_grid,col=4,lwd=2,lty=4)
abline(h=ESS_fpp,lty=2)
abline(h=0,lty=3)
lines(mu_grid,out_ESS_app_1$ESS_grid,lwd=2)

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

plot(mu_grid,out_ESS_app_2$ESS_grid,ylim = c(min(out_ESS_rmap_2_2$ESS_grid,out_ESS_rmap_1_2$ESS_grid)-2,max(out_ESS_app_2$ESS_grid)+2),type = 'l',lwd=2,xlab = expression(theta[C]),ylab = "EHSS")
lines(mu_grid,out_ESS_map_2$ESS_grid,col=2,lwd=2,lty=2)
lines(mu_grid,out_ESS_rmap_1_2$ESS_grid,col=3,lwd=2,lty=3)
lines(mu_grid,out_ESS_rmap_2_2$ESS_grid,col=4,lwd=2,lty=4)
abline(h=ESS_fpp,lty=2)
abline(h=0,lty=3)

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


par(mfrow = c(1, 1))



