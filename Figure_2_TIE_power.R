#### Figure 2: Type I error and power for five calibrated priors

## Source functions
source("case_study_funs.R")

## Load the Type I error and power results
load("Rdatas/Case_study_TIE_power.RData")

## Save current graphical settings
old_par <- par(no.readonly = TRUE)

## Set up a 2 x 3 panel layout
par(
  mfrow = c(2, 3),
  oma = c(4, 0, 0, 0),
  mar = c(4, 4, 2, 1)
)


#------------------------------------------------------------
# Graphical settings
#------------------------------------------------------------

## Line types
lty_TIE_borrow <- 1   # solid
lty_Pow_borrow <- 5   # long-dashed
lty_TIE_NB     <- 2   # dashed
lty_Pow_NB     <- 4   # dot-dashed
lty_hist       <- 2   # dashed
lty_region     <- 3   # dotted

## Greyscale levels
## gray(0) = black; gray(1) = white
col_TIE_borrow <- gray(0.00)   # black
col_Pow_borrow <- gray(0.30)   # dark grey
col_TIE_NB     <- gray(0.45)   # medium grey
col_Pow_NB     <- gray(0.60)   # lighter grey
col_hist       <- gray(0.35)   # dark-medium grey
col_region     <- gray(0.65)   # light grey

## Line widths
lwd_borrow <- 2
lwd_NB     <- 1.2
lwd_ref    <- 1


#------------------------------------------------------------
# Panel A: FPP
#------------------------------------------------------------

plot(
  mu_grid,
  TIE_fpp,
  type = "l",
  lwd = lwd_borrow,
  lty = lty_TIE_borrow,
  col = col_TIE_borrow,
  ylim = c(0, max(1, Power_fpp)),
  xlab = expression(theta[C]),
  ylab = "TIE or Power",
  main = "(A) FPP",
  font.main = 1
)

## No-borrowing TIE
abline(
  h = alpha,
  lty = lty_TIE_NB,
  lwd = lwd_NB,
  col = col_TIE_NB
)

## No-borrowing power
abline(
  h = 1 - beta,
  lty = lty_Pow_NB,
  lwd = lwd_NB,
  col = col_Pow_NB
)

## Historical mean
abline(
  v = xbar_ch,
  lty = lty_hist,
  lwd = lwd_ref,
  col = col_hist
)

## Prespecified discrepancy-region boundaries
abline(
  v = xbar_ch + DEL_l,
  lty = lty_region,
  lwd = lwd_ref,
  col = col_region
)

abline(
  v = xbar_ch + DEL_u,
  lty = lty_region,
  lwd = lwd_ref,
  col = col_region
)

## Borrowing power
lines(
  mu_grid,
  Power_fpp,
  lwd = lwd_borrow,
  lty = lty_Pow_borrow,
  col = col_Pow_borrow
)

axis(
  side = 1,
  at = c(xbar_ch + DEL_l, xbar_ch + DEL_u),
  labels = c(
    expression(bar(y)[Ch] + Delta[l]),
    expression(bar(y)[Ch] + Delta[u])
  ),
  tick = TRUE,
  line = 1
)


#------------------------------------------------------------
# Panel B: APP
#------------------------------------------------------------

plot(
  mu_grid,
  TIE_app_ps,
  type = "l",
  lwd = lwd_borrow,
  lty = lty_TIE_borrow,
  col = col_TIE_borrow,
  ylim = c(0, max(1, Power_app_ps)),
  xlab = expression(theta[C]),
  ylab = "TIE or Power",
  main = expression("(B) " * APP)
)

## No-borrowing TIE
abline(
  h = alpha,
  lty = lty_TIE_NB,
  lwd = lwd_NB,
  col = col_TIE_NB
)

## No-borrowing power
abline(
  h = 1 - beta,
  lty = lty_Pow_NB,
  lwd = lwd_NB,
  col = col_Pow_NB
)

## Historical mean
abline(
  v = xbar_ch,
  lty = lty_hist,
  lwd = lwd_ref,
  col = col_hist
)

## Prespecified discrepancy-region boundaries
abline(
  v = xbar_ch + DEL_l,
  lty = lty_region,
  lwd = lwd_ref,
  col = col_region
)

abline(
  v = xbar_ch + DEL_u,
  lty = lty_region,
  lwd = lwd_ref,
  col = col_region
)

## Borrowing power
lines(
  mu_grid,
  Power_app_ps,
  lwd = lwd_borrow,
  lty = lty_Pow_borrow,
  col = col_Pow_borrow
)

axis(
  side = 1,
  at = c(xbar_ch + DEL_l, xbar_ch + DEL_u),
  labels = c(
    expression(bar(y)[Ch] + Delta[l]),
    expression(bar(y)[Ch] + Delta[u])
  ),
  tick = TRUE,
  line = 1
)


#------------------------------------------------------------
# Third position: common legend
#------------------------------------------------------------

plot.new()

legend(
  "center",
  cex = 1.4,
  legend = c(
    "TIE: Bayesian Prior",
    "Power: Bayesian Prior",
    "TIE: No Borrowing",
    "Power: No Borrowing"
  ),
  col = c(
    col_TIE_borrow,
    col_Pow_borrow,
    col_TIE_NB,
    col_Pow_NB
  ),
  lty = c(
    lty_TIE_borrow,
    lty_Pow_borrow,
    lty_TIE_NB,
    lty_Pow_NB
  ),
  lwd = c(
    lwd_borrow,
    lwd_borrow,
    lwd_NB,
    lwd_NB
  ),
  bty = "n"
)


#------------------------------------------------------------
# Panel C: MAP
#------------------------------------------------------------

plot(
  mu_grid,
  TIE_map,
  type = "l",
  lwd = lwd_borrow,
  lty = lty_TIE_borrow,
  col = col_TIE_borrow,
  ylim = c(0, max(1, Power_map)),
  xlab = expression(theta[C]),
  ylab = "TIE or Power",
  main = "(C) MAP",
  font.main = 1
)

## No-borrowing TIE
abline(
  h = alpha,
  lty = lty_TIE_NB,
  lwd = lwd_NB,
  col = col_TIE_NB
)

## No-borrowing power
abline(
  h = 1 - beta,
  lty = lty_Pow_NB,
  lwd = lwd_NB,
  col = col_Pow_NB
)

## Historical mean
abline(
  v = xbar_ch,
  lty = lty_hist,
  lwd = lwd_ref,
  col = col_hist
)

## Prespecified discrepancy-region boundaries
abline(
  v = xbar_ch + DEL_l,
  lty = lty_region,
  lwd = lwd_ref,
  col = col_region
)

abline(
  v = xbar_ch + DEL_u,
  lty = lty_region,
  lwd = lwd_ref,
  col = col_region
)

## Borrowing power
lines(
  mu_grid,
  Power_map,
  lwd = lwd_borrow,
  lty = lty_Pow_borrow,
  col = col_Pow_borrow
)

axis(
  side = 1,
  at = c(xbar_ch + DEL_l, xbar_ch + DEL_u),
  labels = c(
    expression(bar(y)[Ch] + Delta[l]),
    expression(bar(y)[Ch] + Delta[u])
  ),
  tick = TRUE,
  line = 1
)


#------------------------------------------------------------
# Panel D: rMAP_1
#------------------------------------------------------------

plot(
  mu_grid,
  TIE_rmap_1,
  type = "l",
  lwd = lwd_borrow,
  lty = lty_TIE_borrow,
  col = col_TIE_borrow,
  ylim = c(0, max(1, Power_rmap_1)),
  xlab = expression(theta[C]),
  ylab = "TIE or Power",
  main = expression("(D) " * rMAP[1])
)

## No-borrowing TIE
abline(
  h = alpha,
  lty = lty_TIE_NB,
  lwd = lwd_NB,
  col = col_TIE_NB
)

## No-borrowing power
abline(
  h = 1 - beta,
  lty = lty_Pow_NB,
  lwd = lwd_NB,
  col = col_Pow_NB
)

## Historical mean
abline(
  v = xbar_ch,
  lty = lty_hist,
  lwd = lwd_ref,
  col = col_hist
)

## Prespecified discrepancy-region boundaries
abline(
  v = xbar_ch + DEL_l,
  lty = lty_region,
  lwd = lwd_ref,
  col = col_region
)

abline(
  v = xbar_ch + DEL_u,
  lty = lty_region,
  lwd = lwd_ref,
  col = col_region
)

## Borrowing power
lines(
  mu_grid,
  Power_rmap_1,
  lwd = lwd_borrow,
  lty = lty_Pow_borrow,
  col = col_Pow_borrow
)

axis(
  side = 1,
  at = c(xbar_ch + DEL_l, xbar_ch + DEL_u),
  labels = c(
    expression(bar(y)[Ch] + Delta[l]),
    expression(bar(y)[Ch] + Delta[u])
  ),
  tick = TRUE,
  line = 1
)


#------------------------------------------------------------
# Panel E: rMAP_2
#------------------------------------------------------------

plot(
  mu_grid,
  TIE_rmap_2,
  type = "l",
  lwd = lwd_borrow,
  lty = lty_TIE_borrow,
  col = col_TIE_borrow,
  ylim = c(0, max(1, Power_rmap_2)),
  xlab = expression(theta[C]),
  ylab = "TIE or Power",
  main = expression("(E) " * rMAP[2])
)

## No-borrowing TIE
abline(
  h = alpha,
  lty = lty_TIE_NB,
  lwd = lwd_NB,
  col = col_TIE_NB
)

## No-borrowing power
abline(
  h = 1 - beta,
  lty = lty_Pow_NB,
  lwd = lwd_NB,
  col = col_Pow_NB
)

## Historical mean
abline(
  v = xbar_ch,
  lty = lty_hist,
  lwd = lwd_ref,
  col = col_hist
)

## Prespecified discrepancy-region boundaries
abline(
  v = xbar_ch + DEL_l,
  lty = lty_region,
  lwd = lwd_ref,
  col = col_region
)

abline(
  v = xbar_ch + DEL_u,
  lty = lty_region,
  lwd = lwd_ref,
  col = col_region
)

## Borrowing power
lines(
  mu_grid,
  Power_rmap_2,
  lwd = lwd_borrow,
  lty = lty_Pow_borrow,
  col = col_Pow_borrow
)

axis(
  side = 1,
  at = c(xbar_ch + DEL_l, xbar_ch + DEL_u),
  labels = c(
    expression(bar(y)[Ch] + Delta[l]),
    expression(bar(y)[Ch] + Delta[u])
  ),
  tick = TRUE,
  line = 1
)


## Restore the original graphical settings
par(old_par)
