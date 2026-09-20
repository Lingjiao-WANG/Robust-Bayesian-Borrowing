#### Figure 3: EHSS and expected EHSS for the five calibrated priors

## Load the EHSS results if not already available in the workspace
load("Rdatas/Case_study_EHSS.RData")

## Define a common y-axis range for both panels
ylim_EHSS <- c(
  min(
    out_ESS_rmap_2_2$ESS_grid,
    out_ESS_rmap_1_2$ESS_grid
  ) - 2,
  max(out_ESS_app_2$ESS_grid) + 2
)

## Save the current graphical settings
old_par <- par(no.readonly = TRUE)

## Set up a two-panel figure
par(mfrow = c(1, 2))

## Define line types
lty_app   <- 1   # solid
lty_map   <- 2   # dashed
lty_rmap1 <- 3   # dotted
lty_rmap2 <- 4   # dot-dashed
lty_fpp   <- 5   # long-dashed
lty_zero  <- 3   # dotted

## Define greyscale levels
## 0 = black, 1 = white
col_app   <- gray(0.00)   # black
col_map   <- gray(0.25)   # dark grey
col_rmap1 <- gray(0.45)   # medium grey
col_rmap2 <- gray(0.60)   # lighter grey
col_fpp   <- gray(0.25)   # dark grey reference
col_zero  <- gray(0.65)   # light grey reference

## Define line widths
lwd_method <- 2
lwd_ref    <- 1


#---- Panel A: EHSS as a function of observed xbar_C ----

plot(
  mu_grid,
  out_ESS_map_1$ESS_grid,
  type = "l",
  col = col_map,
  lwd = lwd_method,
  lty = lty_map,
  ylim = ylim_EHSS,
  xlab = expression(bar(y)[C]),
  ylab = "EHSS"
)

lines(
  mu_grid,
  out_ESS_rmap_1_1$ESS_grid,
  col = col_rmap1,
  lwd = lwd_method,
  lty = lty_rmap1
)

lines(
  mu_grid,
  out_ESS_rmap_2_1$ESS_grid,
  col = col_rmap2,
  lwd = lwd_method,
  lty = lty_rmap2
)

lines(
  mu_grid,
  out_ESS_app_1$ESS_grid,
  col = col_app,
  lwd = lwd_method,
  lty = lty_app
)

## FPP reference
abline(
  h = ESS_fpp,
  col = col_fpp,
  lty = lty_fpp,
  lwd = lwd_ref
)

## EHSS = 0 reference
abline(
  h = 0,
  col = col_zero,
  lty = lty_zero,
  lwd = lwd_ref
)

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
  col = c(
    col_app,
    col_map,
    col_rmap1,
    col_rmap2,
    col_fpp,
    col_zero
  ),
  lty = c(
    lty_app,
    lty_map,
    lty_rmap1,
    lty_rmap2,
    lty_fpp,
    lty_zero
  ),
  lwd = c(
    lwd_method,
    lwd_method,
    lwd_method,
    lwd_method,
    lwd_ref,
    lwd_ref
  ),
  bty = "n",
  cex = 0.9
)



#---- Panel B: Expected EHSS given the true control-arm mean ----

plot(
  mu_grid,
  out_ESS_app_2$ESS_grid,
  type = "l",
  col = col_app,
  lwd = lwd_method,
  lty = lty_app,
  ylim = ylim_EHSS,
  xlab = expression(theta[C]),
  ylab = "EHSS"
)

lines(
  mu_grid,
  out_ESS_map_2$ESS_grid,
  col = col_map,
  lwd = lwd_method,
  lty = lty_map
)

lines(
  mu_grid,
  out_ESS_rmap_1_2$ESS_grid,
  col = col_rmap1,
  lwd = lwd_method,
  lty = lty_rmap1
)

lines(
  mu_grid,
  out_ESS_rmap_2_2$ESS_grid,
  col = col_rmap2,
  lwd = lwd_method,
  lty = lty_rmap2
)

## FPP reference
abline(
  h = ESS_fpp,
  col = col_fpp,
  lty = lty_fpp,
  lwd = lwd_ref
)

## EHSS = 0 reference
abline(
  h = 0,
  col = col_zero,
  lty = lty_zero,
  lwd = lwd_ref
)

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
  col = c(
    col_app,
    col_map,
    col_rmap1,
    col_rmap2,
    col_fpp,
    col_zero
  ),
  lty = c(
    lty_app,
    lty_map,
    lty_rmap1,
    lty_rmap2,
    lty_fpp,
    lty_zero
  ),
  lwd = c(
    lwd_method,
    lwd_method,
    lwd_method,
    lwd_method,
    lwd_ref,
    lwd_ref
  ),
  bty = "n",
  cex = 0.9
)

## Restore the original graphical settings
par(old_par)
