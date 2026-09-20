######################################
#### Real-data example
#### Borrowing from the control arm
#### Compute design-prior weights
#### for average Type I error and power
######################################

library(quadprog)
library(MCMCpack)
library(doParallel)
library(parallel)

## Source functions used in the real-data example
source("case_study_funs.R")

## Load the Type I error and power results
load("Case_study_TIE_power.RData")

## Construct a grid of true current control-arm means
mu_grid <- seq(
  xbar_ch - 200,
  xbar_ch + 200,
  length.out = 1e3
)

## Identify values within the prespecified discrepancy region
## Delta = theta_C - xbar_ch in [DEL_l, DEL_u]
idx_ps <- mu_grid >= xbar_ch + DEL_l &
  mu_grid <= xbar_ch + DEL_u

#### Design-prior density weights ####

#---- 1. Vague design prior ----

## Compute the design-prior weights over the full range
w_vague <- mean_OC_normal_weight(
  mu_grid,
  mean_prior = xbar_ch,
  sd_prior = 8800
)

## Restrict and renormalise the weights over the prespecified discrepancy region
w_vague_ps <- w_vague[idx_ps]
w_vague_ps <- w_vague_ps / sum(w_vague_ps)


#---- 2. Skeptical design prior ----

## Compute the design-prior weights over the full range
w_skep <- mean_OC_normal_weight(
  mu_grid,
  mean_prior = -90,
  sd_prior = 25
)

## Restrict and renormalise the weights over the prespecified discrepancy region
w_skep_ps <- w_skep[idx_ps]
w_skep_ps <- w_skep_ps / sum(w_skep_ps)


#---- 3. Realistic design prior ----

## Compute the design-prior weights over the full range
w_real <- real_prior_muC_weights(
  mu_grid,
  mix_w = c(0.51, 0.44, 0.05),
  mix_mean = c(-51, -46.8, -54.1),
  mix_sd = c(19.9, 7.6, 51.7),
  normalise = TRUE
)$w

## Restrict and renormalise the weights over the prespecified discrepancy region
w_real_ps <- w_real[idx_ps]
w_real_ps <- w_real_ps / sum(w_real_ps)


#---- 4. Calibrated FPP design prior ----

## Compute the design-prior weights over the full range
w_fpp <- mean_OC_normal_weight(
  mu_grid,
  mean_prior = xbar_ch,
  sd_prior = se_ch / sqrt(alpha_0)
)

## Restrict and renormalise the weights over the prespecified discrepancy region
w_fpp_ps <- w_fpp[idx_ps]
w_fpp_ps <- w_fpp_ps / sum(w_fpp_ps)


#---- 5. Calibrated MAP design prior ----

## Compute the design-prior weights over the full range
w_map <- map_prior_muC_weights(
  mu_grid,
  xbar_ch,
  se_ch,
  psi_cal_map,
  n_tau = 401,
  normalise = TRUE
)$w

## Restrict and renormalise the weights over the prespecified discrepancy region
w_map_ps <- w_map[idx_ps]
w_map_ps <- w_map_ps / sum(w_map_ps)


#---- 6. Calibrated rMAP design priors ----

## rMAP_1

## Compute the design-prior weights over the full range
w_rmap_1 <- rmap_prior_muC_weights(
  mu_grid,
  xbar_ch,
  se_ch,
  psi = psi_cal_rmap,
  w_rob = w_rob,
  vague_sd = vague_sd,
  n_tau = 401,
  normalise = TRUE
)$w

## Restrict and renormalise the weights over the prespecified discrepancy region
w_rmap_1_ps <- w_rmap_1[idx_ps]
w_rmap_1_ps <- w_rmap_1_ps / sum(w_rmap_1_ps)


## rMAP_2

## Compute the design-prior weights over the full range
w_rmap_2 <- rmap_prior_muC_weights(
  mu_grid,
  xbar_ch,
  se_ch,
  psi = psi_rmap_2,
  w_rob = w_cal,
  vague_sd = vague_sd,
  n_tau = 401,
  normalise = TRUE
)$w

## Restrict and renormalise the weights over the prespecified discrepancy region
w_rmap_2_ps <- w_rmap_2[idx_ps]
w_rmap_2_ps <- w_rmap_2_ps / sum(w_rmap_2_ps)

#### Compute the average Type I error and power over the full range
#### and over the prespecified discrepancy region

#---- FPP: Full range ----

## Compute the average Type I error under each design prior
mean_TIE_fpp_vague <- sum(w_vague * TIE_fpp)
mean_TIE_fpp_skep <- sum(w_skep * TIE_fpp)
mean_TIE_fpp_real <- sum(w_real * TIE_fpp)
mean_TIE_fpp_fpp <- sum(w_fpp * TIE_fpp)
mean_TIE_fpp_map <- sum(w_map * TIE_fpp)
mean_TIE_fpp_rmap_1 <- sum(w_rmap_1 * TIE_fpp)
mean_TIE_fpp_rmap_2 <- sum(w_rmap_2 * TIE_fpp)

## Compute the average power under each design prior
mean_pow_fpp_vague <- sum(w_vague * Power_fpp)
mean_pow_fpp_skep <- sum(w_skep * Power_fpp)
mean_pow_fpp_real <- sum(w_real * Power_fpp)
mean_pow_fpp_fpp <- sum(w_fpp * Power_fpp)
mean_pow_fpp_map <- sum(w_map * Power_fpp)
mean_pow_fpp_rmap_1 <- sum(w_rmap_1 * Power_fpp)
mean_pow_fpp_rmap_2 <- sum(w_rmap_2 * Power_fpp)

## Summarise the average operating characteristics over the full range
tab_fpp <- as.data.frame(
  rbind(
    TIE = c(
      Prior = "FPP Ave TIE",
      Vague = mean_TIE_fpp_vague,
      Skeptical = mean_TIE_fpp_skep,
      Realistic = mean_TIE_fpp_real,
      FPP = mean_TIE_fpp_fpp,
      MAP = mean_TIE_fpp_map,
      rMAP_1 = mean_TIE_fpp_rmap_1,
      rMAP_2 = mean_TIE_fpp_rmap_2
    ),
    Power = c(
      Prior = "FPP Ave Pow",
      Vague = mean_pow_fpp_vague,
      Skeptical = mean_pow_fpp_skep,
      Realistic = mean_pow_fpp_real,
      FPP = mean_pow_fpp_fpp,
      MAP = mean_pow_fpp_map,
      rMAP_1 = mean_pow_fpp_rmap_1,
      rMAP_2 = mean_pow_fpp_rmap_2
    )
  ),
  stringsAsFactors = FALSE
)

## Convert the numerical columns and round to three decimal places
tab_fpp[-1] <- lapply(
  tab_fpp[-1],
  function(x) round(as.numeric(x), 3)
)


#---- FPP: Prespecified discrepancy region ----

## Restrict the operating characteristics to the prespecified discrepancy region
TIE_fpp_ps <- TIE_fpp[idx_ps]
Power_fpp_ps_region <- Power_fpp[idx_ps]

## Compute the average Type I error under each design prior
mean_TIE_fpp_vague_ps <- sum(w_vague_ps * TIE_fpp_ps)
mean_TIE_fpp_skep_ps <- sum(w_skep_ps * TIE_fpp_ps)
mean_TIE_fpp_real_ps <- sum(w_real_ps * TIE_fpp_ps)
mean_TIE_fpp_fpp_ps <- sum(w_fpp_ps * TIE_fpp_ps)
mean_TIE_fpp_map_ps <- sum(w_map_ps * TIE_fpp_ps)
mean_TIE_fpp_rmap_1_ps <- sum(w_rmap_1_ps * TIE_fpp_ps)
mean_TIE_fpp_rmap_2_ps <- sum(w_rmap_2_ps * TIE_fpp_ps)

## Compute the average power under each design prior
mean_pow_fpp_vague_ps <- sum(w_vague_ps * Power_fpp_ps_region)
mean_pow_fpp_skep_ps <- sum(w_skep_ps * Power_fpp_ps_region)
mean_pow_fpp_real_ps <- sum(w_real_ps * Power_fpp_ps_region)
mean_pow_fpp_fpp_ps <- sum(w_fpp_ps * Power_fpp_ps_region)
mean_pow_fpp_map_ps <- sum(w_map_ps * Power_fpp_ps_region)
mean_pow_fpp_rmap_1_ps <- sum(w_rmap_1_ps * Power_fpp_ps_region)
mean_pow_fpp_rmap_2_ps <- sum(w_rmap_2_ps * Power_fpp_ps_region)

## Compute the average power gain relative to the target power
power_gain_fpp_vague_ps <- mean_pow_fpp_vague_ps - (1 - beta)
power_gain_fpp_skep_ps <- mean_pow_fpp_skep_ps - (1 - beta)
power_gain_fpp_real_ps <- mean_pow_fpp_real_ps - (1 - beta)
power_gain_fpp_fpp_ps <- mean_pow_fpp_fpp_ps - (1 - beta)
power_gain_fpp_map_ps <- mean_pow_fpp_map_ps - (1 - beta)
power_gain_fpp_rmap_1_ps <- mean_pow_fpp_rmap_1_ps - (1 - beta)
power_gain_fpp_rmap_2_ps <- mean_pow_fpp_rmap_2_ps - (1 - beta)

## Summarise the average operating characteristics over the
## prespecified discrepancy region
tab_fpp_ps <- as.data.frame(
  rbind(
    TIE = c(
      Prior = "FPP Ave TIE",
      Vague = mean_TIE_fpp_vague_ps,
      Skeptical = mean_TIE_fpp_skep_ps,
      Realistic = mean_TIE_fpp_real_ps,
      FPP = mean_TIE_fpp_fpp_ps,
      MAP = mean_TIE_fpp_map_ps,
      rMAP_1 = mean_TIE_fpp_rmap_1_ps,
      rMAP_2 = mean_TIE_fpp_rmap_2_ps
    ),
    Power = c(
      Prior = "FPP Ave Pow",
      Vague = mean_pow_fpp_vague_ps,
      Skeptical = mean_pow_fpp_skep_ps,
      Realistic = mean_pow_fpp_real_ps,
      FPP = mean_pow_fpp_fpp_ps,
      MAP = mean_pow_fpp_map_ps,
      rMAP_1 = mean_pow_fpp_rmap_1_ps,
      rMAP_2 = mean_pow_fpp_rmap_2_ps
    ),
    `Power gain` = c(
      Prior = "FPP Pow Gain",
      Vague = power_gain_fpp_vague_ps,
      Skeptical = power_gain_fpp_skep_ps,
      Realistic = power_gain_fpp_real_ps,
      FPP = power_gain_fpp_fpp_ps,
      MAP = power_gain_fpp_map_ps,
      rMAP_1 = power_gain_fpp_rmap_1_ps,
      rMAP_2 = power_gain_fpp_rmap_2_ps
    )
  ),
  stringsAsFactors = FALSE
)

## Convert the numerical columns and round to three decimal places
tab_fpp_ps[-1] <- lapply(
  tab_fpp_ps[-1],
  function(x) round(as.numeric(x), 3)
)

#---- APP: Full range ----

## Compute the average Type I error under each design prior
mean_TIE_app_vague <- sum(w_vague * TIE_app_ps)
mean_TIE_app_skep <- sum(w_skep * TIE_app_ps)
mean_TIE_app_real <- sum(w_real * TIE_app_ps)
mean_TIE_app_fpp <- sum(w_fpp * TIE_app_ps)
mean_TIE_app_map <- sum(w_map * TIE_app_ps)
mean_TIE_app_rmap_1 <- sum(w_rmap_1 * TIE_app_ps)
mean_TIE_app_rmap_2 <- sum(w_rmap_2 * TIE_app_ps)

## Compute the average power under each design prior
mean_pow_app_vague <- sum(w_vague * Power_app_ps)
mean_pow_app_skep <- sum(w_skep * Power_app_ps)
mean_pow_app_real <- sum(w_real * Power_app_ps)
mean_pow_app_fpp <- sum(w_fpp * Power_app_ps)
mean_pow_app_map <- sum(w_map * Power_app_ps)
mean_pow_app_rmap_1 <- sum(w_rmap_1 * Power_app_ps)
mean_pow_app_rmap_2 <- sum(w_rmap_2 * Power_app_ps)

## Summarise the average operating characteristics over the full range
tab_app <- as.data.frame(
  rbind(
    TIE = c(
      Prior = "APP Ave TIE",
      Vague = mean_TIE_app_vague,
      Skeptical = mean_TIE_app_skep,
      Realistic = mean_TIE_app_real,
      FPP = mean_TIE_app_fpp,
      MAP = mean_TIE_app_map,
      rMAP_1 = mean_TIE_app_rmap_1,
      rMAP_2 = mean_TIE_app_rmap_2
    ),
    Power = c(
      Prior = "APP Ave Pow",
      Vague = mean_pow_app_vague,
      Skeptical = mean_pow_app_skep,
      Realistic = mean_pow_app_real,
      FPP = mean_pow_app_fpp,
      MAP = mean_pow_app_map,
      rMAP_1 = mean_pow_app_rmap_1,
      rMAP_2 = mean_pow_app_rmap_2
    )
  ),
  stringsAsFactors = FALSE
)

## Convert the numerical columns and round to three decimal places
tab_app[-1] <- lapply(
  tab_app[-1],
  function(x) round(as.numeric(x), 3)
)


#---- APP: Prespecified discrepancy region ----

## Restrict the operating characteristics to the prespecified discrepancy region
TIE_app_ps_region <- TIE_app_ps[idx_ps]
Power_app_ps_region <- Power_app_ps[idx_ps]

## Compute the average Type I error under each design prior
mean_TIE_app_vague_ps <- sum(w_vague_ps * TIE_app_ps_region)
mean_TIE_app_skep_ps <- sum(w_skep_ps * TIE_app_ps_region)
mean_TIE_app_real_ps <- sum(w_real_ps * TIE_app_ps_region)
mean_TIE_app_fpp_ps <- sum(w_fpp_ps * TIE_app_ps_region)
mean_TIE_app_map_ps <- sum(w_map_ps * TIE_app_ps_region)
mean_TIE_app_rmap_1_ps <- sum(w_rmap_1_ps * TIE_app_ps_region)
mean_TIE_app_rmap_2_ps <- sum(w_rmap_2_ps * TIE_app_ps_region)

## Compute the average power under each design prior
mean_pow_app_vague_ps <- sum(w_vague_ps * Power_app_ps_region)
mean_pow_app_skep_ps <- sum(w_skep_ps * Power_app_ps_region)
mean_pow_app_real_ps <- sum(w_real_ps * Power_app_ps_region)
mean_pow_app_fpp_ps <- sum(w_fpp_ps * Power_app_ps_region)
mean_pow_app_map_ps <- sum(w_map_ps * Power_app_ps_region)
mean_pow_app_rmap_1_ps <- sum(w_rmap_1_ps * Power_app_ps_region)
mean_pow_app_rmap_2_ps <- sum(w_rmap_2_ps * Power_app_ps_region)

## Compute the average power gain relative to the target power
power_gain_app_vague_ps <- mean_pow_app_vague_ps - (1 - beta)
power_gain_app_skep_ps <- mean_pow_app_skep_ps - (1 - beta)
power_gain_app_real_ps <- mean_pow_app_real_ps - (1 - beta)
power_gain_app_fpp_ps <- mean_pow_app_fpp_ps - (1 - beta)
power_gain_app_map_ps <- mean_pow_app_map_ps - (1 - beta)
power_gain_app_rmap_1_ps <- mean_pow_app_rmap_1_ps - (1 - beta)
power_gain_app_rmap_2_ps <- mean_pow_app_rmap_2_ps - (1 - beta)

## Summarise the average operating characteristics over the
## prespecified discrepancy region
tab_app_ps <- as.data.frame(
  rbind(
    TIE = c(
      Prior = "APP Ave TIE",
      Vague = mean_TIE_app_vague_ps,
      Skeptical = mean_TIE_app_skep_ps,
      Realistic = mean_TIE_app_real_ps,
      FPP = mean_TIE_app_fpp_ps,
      MAP = mean_TIE_app_map_ps,
      rMAP_1 = mean_TIE_app_rmap_1_ps,
      rMAP_2 = mean_TIE_app_rmap_2_ps
    ),
    Power = c(
      Prior = "APP Ave Pow",
      Vague = mean_pow_app_vague_ps,
      Skeptical = mean_pow_app_skep_ps,
      Realistic = mean_pow_app_real_ps,
      FPP = mean_pow_app_fpp_ps,
      MAP = mean_pow_app_map_ps,
      rMAP_1 = mean_pow_app_rmap_1_ps,
      rMAP_2 = mean_pow_app_rmap_2_ps
    ),
    `Power gain` = c(
      Prior = "APP Pow Gain",
      Vague = power_gain_app_vague_ps,
      Skeptical = power_gain_app_skep_ps,
      Realistic = power_gain_app_real_ps,
      FPP = power_gain_app_fpp_ps,
      MAP = power_gain_app_map_ps,
      rMAP_1 = power_gain_app_rmap_1_ps,
      rMAP_2 = power_gain_app_rmap_2_ps
    )
  ),
  stringsAsFactors = FALSE
)

## Convert the numerical columns and round to three decimal places
tab_app_ps[-1] <- lapply(
  tab_app_ps[-1],
  function(x) round(as.numeric(x), 3)
)

#---- MAP: Full range ----

## Compute the average Type I error under each design prior
mean_TIE_map_vague <- sum(w_vague * TIE_map)
mean_TIE_map_skep <- sum(w_skep * TIE_map)
mean_TIE_map_real <- sum(w_real * TIE_map)
mean_TIE_map_fpp <- sum(w_fpp * TIE_map)
mean_TIE_map_map <- sum(w_map * TIE_map)
mean_TIE_map_rmap_1 <- sum(w_rmap_1 * TIE_map)
mean_TIE_map_rmap_2 <- sum(w_rmap_2 * TIE_map)

## Compute the average power under each design prior
mean_pow_map_vague <- sum(w_vague * Power_map)
mean_pow_map_skep <- sum(w_skep * Power_map)
mean_pow_map_real <- sum(w_real * Power_map)
mean_pow_map_fpp <- sum(w_fpp * Power_map)
mean_pow_map_map <- sum(w_map * Power_map)
mean_pow_map_rmap_1 <- sum(w_rmap_1 * Power_map)
mean_pow_map_rmap_2 <- sum(w_rmap_2 * Power_map)

## Summarise the average operating characteristics over the full range
tab_map <- as.data.frame(
  rbind(
    TIE = c(
      Prior = "MAP Ave TIE",
      Vague = mean_TIE_map_vague,
      Skeptical = mean_TIE_map_skep,
      Realistic = mean_TIE_map_real,
      FPP = mean_TIE_map_fpp,
      MAP = mean_TIE_map_map,
      rMAP_1 = mean_TIE_map_rmap_1,
      rMAP_2 = mean_TIE_map_rmap_2
    ),
    Power = c(
      Prior = "MAP Ave Pow",
      Vague = mean_pow_map_vague,
      Skeptical = mean_pow_map_skep,
      Realistic = mean_pow_map_real,
      FPP = mean_pow_map_fpp,
      MAP = mean_pow_map_map,
      rMAP_1 = mean_pow_map_rmap_1,
      rMAP_2 = mean_pow_map_rmap_2
    )
  ),
  stringsAsFactors = FALSE
)

## Convert the numerical columns and round to three decimal places
tab_map[-1] <- lapply(
  tab_map[-1],
  function(x) round(as.numeric(x), 3)
)


#---- MAP: Prespecified discrepancy region ----

## Restrict the operating characteristics to the prespecified discrepancy region
TIE_map_ps_region <- TIE_map[idx_ps]
Power_map_ps_region <- Power_map[idx_ps]

## Compute the average Type I error under each design prior
mean_TIE_map_vague_ps <- sum(w_vague_ps * TIE_map_ps_region)
mean_TIE_map_skep_ps <- sum(w_skep_ps * TIE_map_ps_region)
mean_TIE_map_real_ps <- sum(w_real_ps * TIE_map_ps_region)
mean_TIE_map_fpp_ps <- sum(w_fpp_ps * TIE_map_ps_region)
mean_TIE_map_map_ps <- sum(w_map_ps * TIE_map_ps_region)
mean_TIE_map_rmap_1_ps <- sum(w_rmap_1_ps * TIE_map_ps_region)
mean_TIE_map_rmap_2_ps <- sum(w_rmap_2_ps * TIE_map_ps_region)

## Compute the average power under each design prior
mean_pow_map_vague_ps <- sum(w_vague_ps * Power_map_ps_region)
mean_pow_map_skep_ps <- sum(w_skep_ps * Power_map_ps_region)
mean_pow_map_real_ps <- sum(w_real_ps * Power_map_ps_region)
mean_pow_map_fpp_ps <- sum(w_fpp_ps * Power_map_ps_region)
mean_pow_map_map_ps <- sum(w_map_ps * Power_map_ps_region)
mean_pow_map_rmap_1_ps <- sum(w_rmap_1_ps * Power_map_ps_region)
mean_pow_map_rmap_2_ps <- sum(w_rmap_2_ps * Power_map_ps_region)

## Compute the average power gain relative to the target power
power_gain_map_vague_ps <- mean_pow_map_vague_ps - (1 - beta)
power_gain_map_skep_ps <- mean_pow_map_skep_ps - (1 - beta)
power_gain_map_real_ps <- mean_pow_map_real_ps - (1 - beta)
power_gain_map_fpp_ps <- mean_pow_map_fpp_ps - (1 - beta)
power_gain_map_map_ps <- mean_pow_map_map_ps - (1 - beta)
power_gain_map_rmap_1_ps <- mean_pow_map_rmap_1_ps - (1 - beta)
power_gain_map_rmap_2_ps <- mean_pow_map_rmap_2_ps - (1 - beta)

## Summarise the average operating characteristics over the
## prespecified discrepancy region
tab_map_ps <- as.data.frame(
  rbind(
    TIE = c(
      Prior = "MAP Ave TIE",
      Vague = mean_TIE_map_vague_ps,
      Skeptical = mean_TIE_map_skep_ps,
      Realistic = mean_TIE_map_real_ps,
      FPP = mean_TIE_map_fpp_ps,
      MAP = mean_TIE_map_map_ps,
      rMAP_1 = mean_TIE_map_rmap_1_ps,
      rMAP_2 = mean_TIE_map_rmap_2_ps
    ),
    Power = c(
      Prior = "MAP Ave Pow",
      Vague = mean_pow_map_vague_ps,
      Skeptical = mean_pow_map_skep_ps,
      Realistic = mean_pow_map_real_ps,
      FPP = mean_pow_map_fpp_ps,
      MAP = mean_pow_map_map_ps,
      rMAP_1 = mean_pow_map_rmap_1_ps,
      rMAP_2 = mean_pow_map_rmap_2_ps
    ),
    `Power gain` = c(
      Prior = "MAP Pow Gain",
      Vague = power_gain_map_vague_ps,
      Skeptical = power_gain_map_skep_ps,
      Realistic = power_gain_map_real_ps,
      FPP = power_gain_map_fpp_ps,
      MAP = power_gain_map_map_ps,
      rMAP_1 = power_gain_map_rmap_1_ps,
      rMAP_2 = power_gain_map_rmap_2_ps
    )
  ),
  stringsAsFactors = FALSE
)

## Convert the numerical columns and round to three decimal places
tab_map_ps[-1] <- lapply(
  tab_map_ps[-1],
  function(x) round(as.numeric(x), 3)
)

#---- rMAP_1: Full range ----

## Compute the average Type I error under each design prior
mean_TIE_rmap_1_vague <- sum(w_vague * TIE_rmap_1)
mean_TIE_rmap_1_skep <- sum(w_skep * TIE_rmap_1)
mean_TIE_rmap_1_real <- sum(w_real * TIE_rmap_1)
mean_TIE_rmap_1_fpp <- sum(w_fpp * TIE_rmap_1)
mean_TIE_rmap_1_map <- sum(w_map * TIE_rmap_1)
mean_TIE_rmap_1_rmap_1 <- sum(w_rmap_1 * TIE_rmap_1)
mean_TIE_rmap_1_rmap_2 <- sum(w_rmap_2 * TIE_rmap_1)

## Compute the average power under each design prior
mean_pow_rmap_1_vague <- sum(w_vague * Power_rmap_1)
mean_pow_rmap_1_skep <- sum(w_skep * Power_rmap_1)
mean_pow_rmap_1_real <- sum(w_real * Power_rmap_1)
mean_pow_rmap_1_fpp <- sum(w_fpp * Power_rmap_1)
mean_pow_rmap_1_map <- sum(w_map * Power_rmap_1)
mean_pow_rmap_1_rmap_1 <- sum(w_rmap_1 * Power_rmap_1)
mean_pow_rmap_1_rmap_2 <- sum(w_rmap_2 * Power_rmap_1)

## Summarise the average operating characteristics over the full range
tab_rmap_1 <- as.data.frame(
  rbind(
    TIE = c(
      Prior = "rMAP_1 Ave TIE",
      Vague = mean_TIE_rmap_1_vague,
      Skeptical = mean_TIE_rmap_1_skep,
      Realistic = mean_TIE_rmap_1_real,
      FPP = mean_TIE_rmap_1_fpp,
      MAP = mean_TIE_rmap_1_map,
      rMAP_1 = mean_TIE_rmap_1_rmap_1,
      rMAP_2 = mean_TIE_rmap_1_rmap_2
    ),
    Power = c(
      Prior = "rMAP_1 Ave Pow",
      Vague = mean_pow_rmap_1_vague,
      Skeptical = mean_pow_rmap_1_skep,
      Realistic = mean_pow_rmap_1_real,
      FPP = mean_pow_rmap_1_fpp,
      MAP = mean_pow_rmap_1_map,
      rMAP_1 = mean_pow_rmap_1_rmap_1,
      rMAP_2 = mean_pow_rmap_1_rmap_2
    )
  ),
  stringsAsFactors = FALSE
)

## Convert the numerical columns and round to three decimal places
tab_rmap_1[-1] <- lapply(
  tab_rmap_1[-1],
  function(x) round(as.numeric(x), 3)
)


#---- rMAP_1: Prespecified discrepancy region ----

## Restrict the operating characteristics to the prespecified discrepancy region
TIE_rmap_1_ps_region <- TIE_rmap_1[idx_ps]
Power_rmap_1_ps_region <- Power_rmap_1[idx_ps]

## Compute the average Type I error under each design prior
mean_TIE_rmap_1_vague_ps <- sum(w_vague_ps * TIE_rmap_1_ps_region)
mean_TIE_rmap_1_skep_ps <- sum(w_skep_ps * TIE_rmap_1_ps_region)
mean_TIE_rmap_1_real_ps <- sum(w_real_ps * TIE_rmap_1_ps_region)
mean_TIE_rmap_1_fpp_ps <- sum(w_fpp_ps * TIE_rmap_1_ps_region)
mean_TIE_rmap_1_map_ps <- sum(w_map_ps * TIE_rmap_1_ps_region)
mean_TIE_rmap_1_rmap_1_ps <- sum(w_rmap_1_ps * TIE_rmap_1_ps_region)
mean_TIE_rmap_1_rmap_2_ps <- sum(w_rmap_2_ps * TIE_rmap_1_ps_region)

## Compute the average power under each design prior
mean_pow_rmap_1_vague_ps <- sum(w_vague_ps * Power_rmap_1_ps_region)
mean_pow_rmap_1_skep_ps <- sum(w_skep_ps * Power_rmap_1_ps_region)
mean_pow_rmap_1_real_ps <- sum(w_real_ps * Power_rmap_1_ps_region)
mean_pow_rmap_1_fpp_ps <- sum(w_fpp_ps * Power_rmap_1_ps_region)
mean_pow_rmap_1_map_ps <- sum(w_map_ps * Power_rmap_1_ps_region)
mean_pow_rmap_1_rmap_1_ps <- sum(w_rmap_1_ps * Power_rmap_1_ps_region)
mean_pow_rmap_1_rmap_2_ps <- sum(w_rmap_2_ps * Power_rmap_1_ps_region)

## Compute the average power gain relative to the target power
power_gain_rmap_1_vague_ps <- mean_pow_rmap_1_vague_ps - (1 - beta)
power_gain_rmap_1_skep_ps <- mean_pow_rmap_1_skep_ps - (1 - beta)
power_gain_rmap_1_real_ps <- mean_pow_rmap_1_real_ps - (1 - beta)
power_gain_rmap_1_fpp_ps <- mean_pow_rmap_1_fpp_ps - (1 - beta)
power_gain_rmap_1_map_ps <- mean_pow_rmap_1_map_ps - (1 - beta)
power_gain_rmap_1_rmap_1_ps <- mean_pow_rmap_1_rmap_1_ps - (1 - beta)
power_gain_rmap_1_rmap_2_ps <- mean_pow_rmap_1_rmap_2_ps - (1 - beta)

## Summarise the average operating characteristics over the
## prespecified discrepancy region
tab_rmap_1_ps <- as.data.frame(
  rbind(
    TIE = c(
      Prior = "rMAP_1 Ave TIE",
      Vague = mean_TIE_rmap_1_vague_ps,
      Skeptical = mean_TIE_rmap_1_skep_ps,
      Realistic = mean_TIE_rmap_1_real_ps,
      FPP = mean_TIE_rmap_1_fpp_ps,
      MAP = mean_TIE_rmap_1_map_ps,
      rMAP_1 = mean_TIE_rmap_1_rmap_1_ps,
      rMAP_2 = mean_TIE_rmap_1_rmap_2_ps
    ),
    Power = c(
      Prior = "rMAP_1 Ave Pow",
      Vague = mean_pow_rmap_1_vague_ps,
      Skeptical = mean_pow_rmap_1_skep_ps,
      Realistic = mean_pow_rmap_1_real_ps,
      FPP = mean_pow_rmap_1_fpp_ps,
      MAP = mean_pow_rmap_1_map_ps,
      rMAP_1 = mean_pow_rmap_1_rmap_1_ps,
      rMAP_2 = mean_pow_rmap_1_rmap_2_ps
    ),
    `Power gain` = c(
      Prior = "rMAP_1 Pow Gain",
      Vague = power_gain_rmap_1_vague_ps,
      Skeptical = power_gain_rmap_1_skep_ps,
      Realistic = power_gain_rmap_1_real_ps,
      FPP = power_gain_rmap_1_fpp_ps,
      MAP = power_gain_rmap_1_map_ps,
      rMAP_1 = power_gain_rmap_1_rmap_1_ps,
      rMAP_2 = power_gain_rmap_1_rmap_2_ps
    )
  ),
  stringsAsFactors = FALSE
)

## Convert the numerical columns and round to three decimal places
tab_rmap_1_ps[-1] <- lapply(
  tab_rmap_1_ps[-1],
  function(x) round(as.numeric(x), 3)
)


#---- rMAP_2: Full range ----

## Compute the average Type I error under each design prior
mean_TIE_rmap_2_vague <- sum(w_vague * TIE_rmap_2)
mean_TIE_rmap_2_skep <- sum(w_skep * TIE_rmap_2)
mean_TIE_rmap_2_real <- sum(w_real * TIE_rmap_2)
mean_TIE_rmap_2_fpp <- sum(w_fpp * TIE_rmap_2)
mean_TIE_rmap_2_map <- sum(w_map * TIE_rmap_2)
mean_TIE_rmap_2_rmap_1 <- sum(w_rmap_1 * TIE_rmap_2)
mean_TIE_rmap_2_rmap_2 <- sum(w_rmap_2 * TIE_rmap_2)

## Compute the average power under each design prior
mean_pow_rmap_2_vague <- sum(w_vague * Power_rmap_2)
mean_pow_rmap_2_skep <- sum(w_skep * Power_rmap_2)
mean_pow_rmap_2_real <- sum(w_real * Power_rmap_2)
mean_pow_rmap_2_fpp <- sum(w_fpp * Power_rmap_2)
mean_pow_rmap_2_map <- sum(w_map * Power_rmap_2)
mean_pow_rmap_2_rmap_1 <- sum(w_rmap_1 * Power_rmap_2)
mean_pow_rmap_2_rmap_2 <- sum(w_rmap_2 * Power_rmap_2)

## Summarise the average operating characteristics over the full range
tab_rmap_2 <- as.data.frame(
  rbind(
    TIE = c(
      Prior = "rMAP_2 Ave TIE",
      Vague = mean_TIE_rmap_2_vague,
      Skeptical = mean_TIE_rmap_2_skep,
      Realistic = mean_TIE_rmap_2_real,
      FPP = mean_TIE_rmap_2_fpp,
      MAP = mean_TIE_rmap_2_map,
      rMAP_1 = mean_TIE_rmap_2_rmap_1,
      rMAP_2 = mean_TIE_rmap_2_rmap_2
    ),
    Power = c(
      Prior = "rMAP_2 Ave Pow",
      Vague = mean_pow_rmap_2_vague,
      Skeptical = mean_pow_rmap_2_skep,
      Realistic = mean_pow_rmap_2_real,
      FPP = mean_pow_rmap_2_fpp,
      MAP = mean_pow_rmap_2_map,
      rMAP_1 = mean_pow_rmap_2_rmap_1,
      rMAP_2 = mean_pow_rmap_2_rmap_2
    )
  ),
  stringsAsFactors = FALSE
)

## Convert the numerical columns and round to three decimal places
tab_rmap_2[-1] <- lapply(
  tab_rmap_2[-1],
  function(x) round(as.numeric(x), 3)
)


#---- rMAP_2: Prespecified discrepancy region ----

## Restrict the operating characteristics to the prespecified discrepancy region
TIE_rmap_2_ps_region <- TIE_rmap_2[idx_ps]
Power_rmap_2_ps_region <- Power_rmap_2[idx_ps]

## Compute the average Type I error under each design prior
mean_TIE_rmap_2_vague_ps <- sum(w_vague_ps * TIE_rmap_2_ps_region)
mean_TIE_rmap_2_skep_ps <- sum(w_skep_ps * TIE_rmap_2_ps_region)
mean_TIE_rmap_2_real_ps <- sum(w_real_ps * TIE_rmap_2_ps_region)
mean_TIE_rmap_2_fpp_ps <- sum(w_fpp_ps * TIE_rmap_2_ps_region)
mean_TIE_rmap_2_map_ps <- sum(w_map_ps * TIE_rmap_2_ps_region)
mean_TIE_rmap_2_rmap_1_ps <- sum(w_rmap_1_ps * TIE_rmap_2_ps_region)
mean_TIE_rmap_2_rmap_2_ps <- sum(w_rmap_2_ps * TIE_rmap_2_ps_region)

## Compute the average power under each design prior
mean_pow_rmap_2_vague_ps <- sum(w_vague_ps * Power_rmap_2_ps_region)
mean_pow_rmap_2_skep_ps <- sum(w_skep_ps * Power_rmap_2_ps_region)
mean_pow_rmap_2_real_ps <- sum(w_real_ps * Power_rmap_2_ps_region)
mean_pow_rmap_2_fpp_ps <- sum(w_fpp_ps * Power_rmap_2_ps_region)
mean_pow_rmap_2_map_ps <- sum(w_map_ps * Power_rmap_2_ps_region)
mean_pow_rmap_2_rmap_1_ps <- sum(w_rmap_1_ps * Power_rmap_2_ps_region)
mean_pow_rmap_2_rmap_2_ps <- sum(w_rmap_2_ps * Power_rmap_2_ps_region)

## Compute the average power gain relative to the target power
power_gain_rmap_2_vague_ps <- mean_pow_rmap_2_vague_ps - (1 - beta)
power_gain_rmap_2_skep_ps <- mean_pow_rmap_2_skep_ps - (1 - beta)
power_gain_rmap_2_real_ps <- mean_pow_rmap_2_real_ps - (1 - beta)
power_gain_rmap_2_fpp_ps <- mean_pow_rmap_2_fpp_ps - (1 - beta)
power_gain_rmap_2_map_ps <- mean_pow_rmap_2_map_ps - (1 - beta)
power_gain_rmap_2_rmap_1_ps <- mean_pow_rmap_2_rmap_1_ps - (1 - beta)
power_gain_rmap_2_rmap_2_ps <- mean_pow_rmap_2_rmap_2_ps - (1 - beta)

## Summarise the average operating characteristics over the
## prespecified discrepancy region
tab_rmap_2_ps <- as.data.frame(
  rbind(
    TIE = c(
      Prior = "rMAP_2 Ave TIE",
      Vague = mean_TIE_rmap_2_vague_ps,
      Skeptical = mean_TIE_rmap_2_skep_ps,
      Realistic = mean_TIE_rmap_2_real_ps,
      FPP = mean_TIE_rmap_2_fpp_ps,
      MAP = mean_TIE_rmap_2_map_ps,
      rMAP_1 = mean_TIE_rmap_2_rmap_1_ps,
      rMAP_2 = mean_TIE_rmap_2_rmap_2_ps
    ),
    Power = c(
      Prior = "rMAP_2 Ave Pow",
      Vague = mean_pow_rmap_2_vague_ps,
      Skeptical = mean_pow_rmap_2_skep_ps,
      Realistic = mean_pow_rmap_2_real_ps,
      FPP = mean_pow_rmap_2_fpp_ps,
      MAP = mean_pow_rmap_2_map_ps,
      rMAP_1 = mean_pow_rmap_2_rmap_1_ps,
      rMAP_2 = mean_pow_rmap_2_rmap_2_ps
    ),
    `Power gain` = c(
      Prior = "rMAP_2 Pow Gain",
      Vague = power_gain_rmap_2_vague_ps,
      Skeptical = power_gain_rmap_2_skep_ps,
      Realistic = power_gain_rmap_2_real_ps,
      FPP = power_gain_rmap_2_fpp_ps,
      MAP = power_gain_rmap_2_map_ps,
      rMAP_1 = power_gain_rmap_2_rmap_1_ps,
      rMAP_2 = power_gain_rmap_2_rmap_2_ps
    )
  ),
  stringsAsFactors = FALSE
)

## Convert the numerical columns and round to three decimal places
tab_rmap_2_ps[-1] <- lapply(
  tab_rmap_2_ps[-1],
  function(x) round(as.numeric(x), 3)
)

#save.image(file = "Rdatas/Case_study_ave_TIE_power.RData")
