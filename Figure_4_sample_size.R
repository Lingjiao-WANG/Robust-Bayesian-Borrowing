#### Figure 4: Minimum required sample size for FPP, APP, and MAP

## Load the minimum required sample-size results
load("Rdatas/Case_study_sample_size.RData")

## Define greyscale levels
cols <- c(
  FPP = gray(0.65),  # light grey
  APP = gray(0.35),  # dark grey
  MAP = gray(0.00)   # black
)

## Plot the minimum required control-arm sample size for the FPP
plot(
  Del_grid,
  nC_min_FPP,
  type = "s",
  lwd = 2,
  lty = 1,
  col = cols["FPP"],
  xlab = expression(Delta[u]),
  ylab = expression(n[C]),
  ylim = c(7, 21),
  las = 1
)

## Add the minimum required control-arm sample size for the APP
lines(
  Del_grid,
  nC_min_APP,
  type = "s",
  lwd = 2,
  lty = 2,
  col = cols["APP"]
)

## Add the minimum required control-arm sample size for the MAP
lines(
  Del_grid,
  nC_min_MAP,
  type = "s",
  lwd = 2,
  lty = 3,
  col = cols["MAP"]
)

## Add the legend
legend(
  "topleft",
  legend = c("FPP", "APP", "MAP"),
  col = cols,
  lty = c(1, 2, 3),
  lwd = 2,
  bty = "n"
)


