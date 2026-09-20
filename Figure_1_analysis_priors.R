#### Figure 1: Calibrated analysis priors

## Source functions used in the real-data example
source("case_study_funs.R")

## Load the calibrated analysis-prior results
load("Rdatas/Calibrated_analysis_priors.RData")

## Construct a grid of current control means for evaluating the prior densities
mu_grid <- seq(
  xbar_ch - 200,
  xbar_ch + 200,
  length.out = 1e5
)

## Standard deviation of the vague component
vague_sd <- 88


#---- FPP ----

## Compute the calibrated FPP prior density
sd_fpp <- se_ch / sqrt(alpha_0)

dens_fpp <- dnorm(
  mu_grid,
  mean = xbar_ch,
  sd = sd_fpp
)


#---- Vague prior ----

## Compute the vague prior density
dens_vague <- dnorm(
  mu_grid,
  mean = xbar_ch,
  sd = vague_sd
)


#---- MAP ----

## Compute the calibrated MAP prior density
dens_map <- dmap_prior(
  theta_C = mu_grid,
  xbar_ch = xbar_ch,
  se_ch = se_ch,
  psi = psi_cal_map,
  n_tau = 401
)


#---- rMAP ----

## Compute the calibrated rMAP_1 prior density
dens_rmap_1 <- drmap_prior(
  mu_grid = mu_grid,
  xbar_ch = xbar_ch,
  se_ch = se_ch,
  psi = psi_cal_rmap,
  w_rob = w_rob,
  vague_sd = vague_sd,
  n_tau = 401
)

## Compute the calibrated rMAP_2 prior density
dens_rmap_2 <- drmap_prior(
  mu_grid = mu_grid,
  xbar_ch = xbar_ch,
  se_ch = se_ch,
  psi = psi_rmap_2,
  w_rob = w_cal,
  vague_sd = vague_sd,
  n_tau = 401
)


#### Plot Figure 1

##------------------------------------------------------------
## Graphical settings
##------------------------------------------------------------

## Line types
lty_fpp   <- 1   # solid
lty_map   <- 2   # dashed
lty_rmap1 <- 3   # dotted
lty_rmap2 <- 4   # dot-dashed
lty_vague <- 5   # long-dashed
lty_hist  <- 2   # historical mean

## Greyscale levels
## gray(0) = black; gray(1) = white
col_fpp   <- gray(0.00)   # black
col_map   <- gray(0.20)   # very dark grey
col_rmap1 <- gray(0.40)   # medium-dark grey
col_rmap2 <- gray(0.55)   # medium grey
col_vague <- gray(0.70)   # light grey
col_hist  <- gray(0.50)   # reference line

## Line widths
lwd_fpp   <- 2.0
lwd_map   <- 2.5
lwd_rmap1 <- 2.2
lwd_rmap2 <- 2.2
lwd_vague <- 1.5
lwd_hist  <- 1


##------------------------------------------------------------
## Plot the calibrated prior densities
##------------------------------------------------------------

plot(
  mu_grid,
  dens_fpp,
  type = "l",
  lwd = lwd_fpp,
  lty = lty_fpp,
  col = col_fpp,
  ylim = c(
    0,
    max(
      dens_fpp,
      dens_map,
      dens_rmap_1$rMAP,
      dens_rmap_2$rMAP,
      dens_vague
    ) + 0.001
  ),
  xlab = expression(theta[C]),
  ylab = "Density"
)

## Add the historical control sample mean
abline(
  v = xbar_ch,
  lty = lty_hist,
  lwd = lwd_hist,
  col = col_hist
)

## Add the vague prior
lines(
  mu_grid,
  dens_vague,
  lwd = lwd_vague,
  lty = lty_vague,
  col = col_vague
)

## Add the MAP prior
lines(
  mu_grid,
  dens_map,
  lwd = lwd_map,
  lty = lty_map,
  col = col_map
)

## Add the rMAP_1 prior
lines(
  mu_grid,
  dens_rmap_1$rMAP,
  lwd = lwd_rmap1,
  lty = lty_rmap1,
  col = col_rmap1
)

## Add the rMAP_2 prior
lines(
  mu_grid,
  dens_rmap_2$rMAP,
  lwd = lwd_rmap2,
  lty = lty_rmap2,
  col = col_rmap2
)


##------------------------------------------------------------
## Add the legend
##------------------------------------------------------------

legend(
  "topright",
  legend = c(
    bquote(
      "FPP (" * alpha[0] * " = " *
        .(round(alpha_0, 3)) * ")"
    ),
    bquote(
      "MAP (" * psi * " = " *
        .(round(psi_cal_map, 3)) * ")"
    ),
    bquote(
      "rMAP"[1] * " (" *
        w * " = " * .(round(w_rob, 2)) * ", " *
        psi * " = " * .(round(psi_cal_rmap, 3)) * ")"
    ),
    bquote(
      "rMAP"[2] * " (" *
        w * " = " * .(round(w_cal, 3)) * ", " *
        psi * " = " * .(round(psi_rmap_2, 2)) * ")"
    ),
    bquote(
      "Vague prior: " * N(-49, 88^2)
    )
  ),
  lty = c(
    lty_fpp,
    lty_map,
    lty_rmap1,
    lty_rmap2,
    lty_vague
  ),
  lwd = c(
    lwd_fpp,
    lwd_map,
    lwd_rmap1,
    lwd_rmap2,
    lwd_vague
  ),
  col = c(
    col_fpp,
    col_map,
    col_rmap1,
    col_rmap2,
    col_vague
  ),
  bty = "n"
)
