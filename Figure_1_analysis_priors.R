#### Figure 1: Calibrated analysis priors ####

## Source functions
source("case_study_funs.R")

## Source the data of calibrated analysis priors
load("Calibrated_analysis_priors.RData")

# --- FPP ----
sd_fpp <- se_ch / sqrt(alpha_0)

dens_fpp <- dnorm(
  mu_grid,
  mean = xbar_ch,
  sd   = sd_fpp
)

# --- Vague ----

dens_vague <- dnorm(mu_grid,
                    mean = xbar_ch,
                    sd   = 88)

# --- MAP ----

dens_map <- dmap_prior(
  theta_C = mu_grid,
  xbar_ch = xbar_ch,
  se_ch = se_ch,
  psi = psi_cal_map,
  n_tau = 401
)

# --- rMAP ----
# rMAP_1
dens_rmap_1 <- drmap_prior(
  mu_grid = mu_grid,
  xbar_ch = xbar_ch,
  se_ch = se_ch,
  psi = psi_cal_rmap,
  w_rob = w_rob,
  vague_sd = 88,
  n_tau = 401
)

# rMAP_2
dens_rmap_2 <- drmap_prior(
  mu_grid = mu_grid,
  xbar_ch = xbar_ch,
  se_ch = se_ch,
  psi = 0.35*sigma,
  w_rob = w_cal,
  vague_sd = 88,
  n_tau = 401
)

## Figure 1

plot(
  mu_grid, dens_fpp,
  type = "l",
  lwd = 2,
  xlab = expression(theta[C]),
  ylab = "Density"#,
  #main = "Prior Density"
)

abline(v = xbar_ch, lty = 2)

lines(mu_grid, dens_vague, lwd = 2, lty = 2)

lines(mu_grid, dens_map, lwd = 2, lty = 1, col=2)

lines(mu_grid, dens_rmap_1$rMAP,lwd = 2,col=3)
lines(mu_grid, dens_rmap_2$rMAP,lwd = 2,col=3,lty=2)

legend(
  "topright",
  legend = c(
    expression(paste("FPP (", alpha[0], " = 0.0465)")),
    expression(paste("MAP (", psi, " = 85.01)")),
    expression(paste("rMAP"[1], " (", w, " = 0.80, ", psi, " = 61.20)")),
    expression(paste("rMAP"[2], " (", w, " = 0.44, ", psi, " = 30.80)")),
    expression(paste("Vague prior: ", N(-49, 88^2)))
  ),
  lwd = 2,
  lty = c(1, 1, 1, 2, 2),
  col = c(1, 2, 3, 3, 1),
  bty = "n"
)
