#### Figure A2: Type I error and power for the calibrated APP

## Load the saved APP operating-characteristic results
source("Rdatas/APP_ss_results.txt")


#---- Determine the sweet-spot boundaries ----

## Lower boundary of the sweet-spot region
theta_C_ss_l <- xbar_ch - ss_asy_list_within_ps$DEL_u

## Upper boundary of the sweet-spot region
theta_C_ss_u <- mu_grid[
  max(which(Power_app_ss >= 1 - beta))
]


#---- Plot the operating characteristics ----

## Plot the Type I error
plot(
  mu_grid,
  TIE_app_ss,
  type = "l",
  lwd = 2,
  ylim = c(
    0,
    max(1, Power_app_ss)
  ),
  xlab = expression(theta[C]),
  ylab = "TIE or Power"
)

## Add the power curve
lines(
  mu_grid,
  Power_app_ss,
  lwd = 2,
  col = 2
)


#---- Add reference lines ----

## Nominal Type I error and target power
abline(h = alpha, lty = 2)
abline(h = 1 - beta, lty = 2, col = 2)

## Historical control sample mean
abline(v = xbar_ch, lty = 2)

## Prespecified discrepancy region
abline(v = xbar_ch + DEL_l, lty = 3)
abline(v = xbar_ch + DEL_u, lty = 3)

## Sweet-spot region
abline(
  v = theta_C_ss_l,
  lty = 3,
  col = 2
)

abline(
  v = theta_C_ss_u,
  lty = 3,
  col = 2
)


#---- Label the discrepancy and sweet-spot regions ----

## Mark the boundaries of the sweet-spot region
axis(
  side = 1,
  at = c(
    theta_C_ss_l,
    theta_C_ss_u
  ),
  labels = c(
    expression(theta[C * ";" * l]^ss),
    expression(theta[C * ";" * u]^ss)
  ),
  tick = TRUE,
  line = 0.3,
  col.axis = "red",
  col.ticks = "red",
  cex.axis = 0.9
)

## Mark the boundaries of the prespecified discrepancy region
axis(
  side = 3,
  at = c(
    xbar_ch + DEL_l,
    xbar_ch + DEL_u
  ),
  labels = c(
    expression(bar(y)[Ch] + Delta[l]),
    expression(bar(y)[Ch] + Delta[u])
  ),
  tick = TRUE,
  line = 0.3,
  cex.axis = 0.9
)
