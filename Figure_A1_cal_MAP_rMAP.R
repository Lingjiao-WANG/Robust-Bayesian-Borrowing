#### Figure A1: Numerical calibration of the borrowing parameters
#### for the MAP, rMAP_1, and rMAP_2 priors

library(mgcv)


#---- Load calibration results ----

## Load the simulation results used for GAM calibration
dat_map <- read.table(
  "Rdatas/MAP_cal_psi.txt",
  header = TRUE,
  sep = "\t"
)

dat_rmap1 <- read.table(
  "Rdatas/rMAP1_cal_psi.txt",
  header = TRUE,
  sep = "\t"
)

dat_rmap2 <- read.table(
  "Rdatas/rMAP2_cal_w.txt",
  header = TRUE,
  sep = "\t"
)


#---- Set graphical parameters ----

## Save the current graphical settings
old_par <- par(no.readonly = TRUE)

## Arrange the three calibration plots in one row
par(
  mfrow = c(1, 3),
  mar = c(4.5, 4.5, 3, 1),
  oma = c(0, 0, 0, 0)
)

#---- (A) MAP prior: Calibration of psi ----

## Fit the GAM to the Monte Carlo calibration results
fit_map <- gam(
  res ~ s(psi, k = 6),
  data = dat_map,
  method = "REML"
)

## Plot the Monte Carlo estimates
plot(
  dat_map$psi,
  dat_map$res,
  xlab = expression(psi),
  ylab = expression(
    underline(Pow)(psi, p) - (1 - beta)
  ),
  main = expression("(A) MAP prior"),
  pch = 16,
  cex = 0.8
)

## Construct a fine grid for evaluating the fitted GAM
psi_grid_map <- seq(
  min(dat_map$psi),
  max(dat_map$psi),
  length.out = 500
)

## Obtain predictions from the fitted GAM
pred_map <- predict(
  fit_map,
  newdata = data.frame(psi = psi_grid_map),
  se.fit = TRUE
)

## Add the fitted GAM curve
lines(
  psi_grid_map,
  pred_map$fit,
  lwd = 2
)

## Add the target power constraint
abline(
  h = 0,
  lty = 2,
  lwd = 1
)

## Identify intervals containing zero crossings of the fitted GAM
y_grid_map <- predict(
  fit_map,
  newdata = data.frame(psi = psi_grid_map)
)

idx_map <- which(
  y_grid_map[-length(y_grid_map)] *
    y_grid_map[-1] <= 0
)

## Refine the zero crossings
roots_map <- sapply(idx_map, function(i) {
  
  uniroot(
    function(x) {
      predict(
        fit_map,
        newdata = data.frame(psi = x)
      )
    },
    interval = c(
      psi_grid_map[i],
      psi_grid_map[i + 1]
    )
  )$root
  
})

## Check that at least one zero crossing is identified
if (length(roots_map) == 0) {
  stop("No zero crossing found for the MAP calibration.")
}

## Select the largest zero crossing as the calibrated value of psi
psi_cal_map <- max(roots_map)
# psi_cal_map = 73.87874

## Add the calibrated value of psi
abline(
  v = psi_cal_map,
  lty = 3,
  lwd = 2,
  col = 2
)

points(
  psi_cal_map,
  0,
  pch = 16,
  col = 2,
  cex = 1.1
)

axis(
  side = 1,
  at = psi_cal_map,
  labels = round(psi_cal_map, 3),
  col.axis = 2,
  col.ticks = 2
)

## Add the legend
legend(
  "bottomright",
  legend = as.expression(list(
    "Monte Carlo estimates",
    "Fitted GAM",
    bquote(
      "Calibrated " ~ psi ==
        .(round(psi_cal_map, 3))
    )
  )),
  pch = c(16, NA, 16),
  lty = c(NA, 1, NA),
  lwd = c(NA, 2, NA),
  col = c(1, 1, 2),
  bty = "n",
  cex = 0.9
)


#---- (B) rMAP_1 prior: Calibration of psi ----

## Fit the GAM to the Monte Carlo calibration results
fit_rmap1 <- gam(
  res ~ s(psi, k = 6),
  data = dat_rmap1,
  method = "REML"
)

## Plot the Monte Carlo estimates
plot(
  dat_rmap1$psi,
  dat_rmap1$res,
  xlab = expression(psi),
  ylab = expression(
    underline(Pow)(psi, p) - (1 - beta)
  ),
  main = expression("(B) " * rMAP[1] * " prior"),
  pch = 16,
  cex = 0.8
)

## Construct a fine grid for evaluating the fitted GAM
psi_grid_rmap1 <- seq(
  min(dat_rmap1$psi),
  max(dat_rmap1$psi),
  length.out = 500
)

## Obtain predictions from the fitted GAM
pred_rmap1 <- predict(
  fit_rmap1,
  newdata = data.frame(
    psi = psi_grid_rmap1
  )
)

## Add the fitted GAM curve
lines(
  psi_grid_rmap1,
  pred_rmap1,
  lwd = 2
)

## Add the target power constraint
abline(
  h = 0,
  lty = 2,
  lwd = 1
)

## Identify intervals containing zero crossings of the fitted GAM
idx_rmap1 <- which(
  pred_rmap1[-length(pred_rmap1)] *
    pred_rmap1[-1] <= 0
)

## Refine the zero crossings
roots_rmap1 <- sapply(idx_rmap1, function(i) {
  
  uniroot(
    function(x) {
      predict(
        fit_rmap1,
        newdata = data.frame(psi = x)
      )
    },
    interval = c(
      psi_grid_rmap1[i],
      psi_grid_rmap1[i + 1]
    )
  )$root
  
})

## Check that at least one zero crossing is identified
if (length(roots_rmap1) == 0) {
  stop("No zero crossing found for the rMAP_1 calibration.")
}

## Select the largest zero crossing as the calibrated value of psi
psi_cal_rmap1 <- max(roots_rmap1)
# psi_cal_rmap1 = 54.63299

## Add the calibrated value of psi
abline(
  v = psi_cal_rmap1,
  lty = 3,
  lwd = 2,
  col = 2
)

points(
  psi_cal_rmap1,
  0,
  pch = 16,
  col = 2,
  cex = 1.1
)

axis(
  side = 1,
  at = psi_cal_rmap1,
  labels = round(psi_cal_rmap1, 3),
  col.axis = 2,
  col.ticks = 2
)

## Add the legend
legend(
  "bottomright",
  legend = as.expression(list(
    "Monte Carlo estimates",
    "Fitted GAM",
    bquote(
      "Calibrated " ~ psi ==
        .(round(psi_cal_rmap1, 3))
    )
  )),
  pch = c(16, NA, 16),
  lty = c(NA, 1, NA),
  lwd = c(NA, 2, NA),
  col = c(1, 1, 2),
  bty = "n",
  cex = 0.9
)

#---- (C) rMAP_2 prior: Calibration of w ----

## Fit the GAM to the Monte Carlo calibration results
fit_rmap2 <- gam(
  res ~ s(w, k = 6),
  data = dat_rmap2,
  method = "REML"
)

## Plot the Monte Carlo estimates
plot(
  dat_rmap2$w,
  dat_rmap2$res,
  xlab = expression(w),
  ylab = expression(
    underline(Pow)(w, p) - (1 - beta)
  ),
  main = expression("(C) " * rMAP[2] * " prior"),
  pch = 16,
  cex = 0.8
)

## Construct a fine grid for evaluating the fitted GAM
w_grid_rmap2 <- seq(
  min(dat_rmap2$w),
  max(dat_rmap2$w),
  length.out = 500
)

## Obtain predictions from the fitted GAM
pred_rmap2 <- predict(
  fit_rmap2,
  newdata = data.frame(
    w = w_grid_rmap2
  )
)

## Add the fitted GAM curve
lines(
  w_grid_rmap2,
  pred_rmap2,
  lwd = 2
)

## Add the target power constraint
abline(
  h = 0,
  lty = 2,
  lwd = 1
)

## Identify intervals containing zero crossings of the fitted GAM
idx_rmap2 <- which(
  pred_rmap2[-length(pred_rmap2)] *
    pred_rmap2[-1] <= 0
)

## Refine the zero crossings
roots_rmap2 <- sapply(idx_rmap2, function(i) {
  
  uniroot(
    function(x) {
      predict(
        fit_rmap2,
        newdata = data.frame(w = x)
      )
    },
    interval = c(
      w_grid_rmap2[i],
      w_grid_rmap2[i + 1]
    )
  )$root
  
})

## Check that at least one zero crossing is identified
if (length(roots_rmap2) == 0) {
  stop("No zero crossing found for the rMAP_2 calibration.")
}

## Select the smallest zero crossing as the calibrated value of w
w_cal_rmap2 <- min(roots_rmap2)
# w_cal_rmap2 = 0.4855757

## Add the calibrated value of w
abline(
  v = w_cal_rmap2,
  lty = 3,
  lwd = 2,
  col = 2
)

points(
  w_cal_rmap2,
  0,
  pch = 16,
  col = 2,
  cex = 1.1
)

axis(
  side = 1,
  at = w_cal_rmap2,
  labels = round(w_cal_rmap2, 3),
  col.axis = 2,
  col.ticks = 2
)

## Add the legend
legend(
  "topright",
  legend = as.expression(list(
    "Monte Carlo estimates",
    "Fitted GAM",
    bquote(
      "Calibrated " ~ w ==
        .(round(w_cal_rmap2, 3))
    )
  )),
  pch = c(16, NA, 16),
  lty = c(NA, 1, NA),
  lwd = c(NA, 2, NA),
  col = c(1, 1, 2),
  bty = "n",
  cex = 0.9
)


#---- Restore graphical settings ----

## Restore the original graphical settings
par(old_par)

