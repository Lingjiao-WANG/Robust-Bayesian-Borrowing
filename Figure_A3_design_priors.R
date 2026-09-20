#### Figure A3: Design priors used to compute the average
#### Type I error and power

## Load the calibrated analysis-prior densities used here as design priors
load("Rdatas/Figure_1_analysis_priors.RData")

## Define the fixed scale parameter for rMAP_2
psi_rmap_2 <- 0.35 * sigma


#---- Design-prior densities ----

## Very vague prior: N(xbar_ch, 8800^2)
dens_vague <- dnorm(
  mu_grid,
  mean = xbar_ch,
  sd = 8800
)

## Skeptical prior: N(-90, 25^2)
dens_skeptical <- dnorm(
  mu_grid,
  mean = -90,
  sd = 25
)

## Realistic mixture prior:
## 0.51 N(-51.0, 19.9^2)
## + 0.44 N(-46.8, 7.6^2)
## + 0.05 N(-54.1, 51.7^2)
dens_realistic <-
  0.51 * dnorm(mu_grid, mean = -51.0, sd = 19.9) +
  0.44 * dnorm(mu_grid, mean = -46.8, sd = 7.6) +
  0.05 * dnorm(mu_grid, mean = -54.1, sd = 51.7)


#---- Plot Figure A3 ----

plot(
  mu_grid,
  dens_fpp,
  type = "l",
  lwd = 2,
  ylim = c(
    0,
    max(
      dens_fpp,
      dens_map,
      dens_rmap_1$rMAP,
      dens_rmap_2$rMAP,
      dens_vague,
      dens_skeptical,
      dens_realistic
    )
  ),
  xlab = expression(theta[C]),
  ylab = "Density"
)

## Historical control sample mean
abline(v = xbar_ch, lty = 2)

## Calibrated analysis priors used as design priors
lines(mu_grid, dens_map, lwd = 2, lty = 1, col = 2)
lines(mu_grid, dens_rmap_1$rMAP, lwd = 2, lty = 1, col = 3)
lines(mu_grid, dens_rmap_2$rMAP, lwd = 2, lty = 2, col = 3)

## Additional design priors
lines(mu_grid, dens_vague, lwd = 2, lty = 2)
lines(mu_grid, dens_skeptical, lwd = 2, lty = 3, col = 4)
lines(mu_grid, dens_realistic, lwd = 2, lty = 4, col = 5)

legend(
  "topright",
  legend = c(
    bquote(
      "FPP (" * alpha[0] * " = " * .(round(alpha_0, 3)) * ")"
    ),
    bquote(
      "MAP (" * psi * " = " * .(round(psi_cal_map, 3)) * ")"
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
      "Very vague: " * N(.(xbar_ch), 8800^2)
    ),
    expression(
      "Skeptical: " * N(-90, 25^2)
    ),
    "Realistic mixture prior"
  ),
  lwd = rep(2, 7),
  lty = c(1, 1, 1, 2, 2, 3, 4),
  col = c(1, 2, 3, 3, 1, 4, 5),
  bty = "n",
  cex = 0.8
)