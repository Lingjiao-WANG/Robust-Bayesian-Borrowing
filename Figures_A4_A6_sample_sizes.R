#### Figures A4-A6: Minimum required sample size, borrowing parameter,
#### and posterior probability threshold for FPP, APP, and MAP

## Save the current graphical settings
old_par <- par(no.readonly = TRUE)


#---- Figure A4: FPP ----

## Load the FPP sample-size results
FPP_results <- read.table(
  "Rdatas/FPP_sample_size_results.txt",
  header = TRUE,
  sep = "\t"
)

## Reserve additional space on the right for two right-hand axes
par(
  mar = c(5, 5, 4, 8) + 0.1,
  xpd = FALSE
)


## Minimum required control-arm sample size
plot(
  FPP_results$Del_FPP,
  FPP_results$nC_FPP,
  type = "s",
  lwd = 2,
  xlab = expression(Delta[u]),
  ylab = expression(n[C]),
  axes = FALSE
)

axis(1)
axis(2)
box()


## Calibrated FPP borrowing parameter alpha_0
par(new = TRUE)

plot(
  FPP_results$Del_FPP,
  FPP_results$alpha0,
  type = "l",
  lwd = 2,
  lty = 2,
  col = 3,
  axes = FALSE,
  xlab = "",
  ylab = "",
  xlim = range(FPP_results$Del_FPP, na.rm = TRUE),
  ylim = range(FPP_results$alpha0, na.rm = TRUE)
)

axis(
  side = 4
)

mtext(
  expression(alpha[0]),
  side = 4,
  line = 3
)


## Calibrated posterior probability threshold parameter p
par(new = TRUE)

plot(
  FPP_results$Del_FPP,
  FPP_results$p_FPP,
  type = "l",
  lwd = 2,
  lty = 3,
  col = 2,
  axes = FALSE,
  xlab = "",
  ylab = "",
  xlim = range(FPP_results$Del_FPP, na.rm = TRUE),
  ylim = range(FPP_results$p_FPP, na.rm = TRUE)
)

axis(
  side = 4,
  at = pretty(
    range(FPP_results$p_FPP, na.rm = TRUE)
  ),
  line = 4
)

mtext(
  expression(p),
  side = 4,
  line = 7
)


## Add the legend
legend(
  "topright",
  inset = c(0, 0.05),
  legend = c(
    expression(n[C]),
    expression(alpha[0]),
    expression(p)
  ),
  lty = c(1, 2, 3),
  col = c(1, 3, 2),
  lwd = 2,
  bty = "n"
)

#---- Figure A5: APP ----

## Load the APP sample-size results
APP_results <- read.table(
  "Rdatas/APP_sample_size_results.txt",
  header = TRUE,
  sep = "\t"
)

## Reserve additional space on the right for two right-hand axes
par(
  mar = c(5, 5, 4, 8) + 0.1,
  xpd = FALSE
)


## Minimum required control-arm sample size
plot(
  APP_results$Del_APP,
  APP_results$nC_APP,
  type = "s",
  lwd = 2,
  xlab = expression(Delta[u]),
  ylab = expression(n[C]),
  axes = FALSE
)

axis(1)
axis(2)
box()


## Calibrated APP borrowing parameter delta_0
par(new = TRUE)

plot(
  APP_results$Del_APP,
  APP_results$delta0,
  type = "l",
  lwd = 2,
  lty = 2,
  col = 3,
  axes = FALSE,
  xlab = "",
  ylab = "",
  xlim = range(APP_results$Del_APP, na.rm = TRUE),
  ylim = range(APP_results$delta0, na.rm = TRUE)
)

axis(
  side = 4
)

mtext(
  expression(delta[0]),
  side = 4,
  line = 3
)


## Calibrated posterior probability threshold parameter p
par(new = TRUE)

plot(
  APP_results$Del_APP,
  APP_results$p_APP,
  type = "l",
  lwd = 2,
  lty = 3,
  col = 2,
  axes = FALSE,
  xlab = "",
  ylab = "",
  xlim = range(APP_results$Del_APP, na.rm = TRUE),
  ylim = range(APP_results$p_APP, na.rm = TRUE)
)

axis(
  side = 4,
  at = pretty(
    range(APP_results$p_APP, na.rm = TRUE)
  ),
  line = 4
)

mtext(
  expression(p),
  side = 4,
  line = 7
)


## Add the legend
legend(
  "top",
  legend = c(
    expression(n[C]),
    expression(delta[0]),
    expression(p)
  ),
  lty = c(1, 2, 3),
  col = c(1, 3, 2),
  lwd = 2,
  bty = "n"
)

#---- Figure A6: MAP ----

## Load the MAP sample-size results
MAP_results <- read.table(
  "Rdatas/MAP_sample_size_results.txt",
  header = TRUE,
  sep = "\t"
)

## Reserve additional space on the right for two right-hand axes
par(
  mar = c(5, 5, 4, 8) + 0.1,
  xpd = FALSE
)


## Minimum required control-arm sample size
plot(
  MAP_results$Del_MAP,
  MAP_results$nC_MAP,
  type = "s",
  lwd = 2,
  xlab = expression(Delta[u]),
  ylab = expression(n[C]),
  axes = FALSE
)

axis(1)
axis(2)
box()


## Calibrated MAP borrowing parameter psi
par(new = TRUE)

plot(
  MAP_results$Del_MAP,
  MAP_results$psi,
  type = "l",
  lwd = 2,
  lty = 2,
  col = 3,
  axes = FALSE,
  xlab = "",
  ylab = "",
  xlim = range(MAP_results$Del_MAP, na.rm = TRUE),
  ylim = range(MAP_results$psi, na.rm = TRUE)
)

axis(
  side = 4
)

mtext(
  expression(psi),
  side = 4,
  line = 3
)


## Calibrated posterior probability threshold parameter p
par(new = TRUE)

plot(
  MAP_results$Del_MAP,
  MAP_results$p_MAP,
  type = "l",
  lwd = 2,
  lty = 3,
  col = 2,
  axes = FALSE,
  xlab = "",
  ylab = "",
  xlim = range(MAP_results$Del_MAP, na.rm = TRUE),
  ylim = range(MAP_results$p_MAP, na.rm = TRUE)
)

axis(
  side = 4,
  at = pretty(
    range(MAP_results$p_MAP, na.rm = TRUE)
  ),
  line = 4
)

mtext(
  expression(p),
  side = 4,
  line = 7
)


## Add the legend
legend(
  "top",
  legend = c(
    expression(n[C]),
    expression(psi),
    expression(p)
  ),
  lty = c(1, 2, 3),
  col = c(1, 3, 2),
  lwd = 2,
  bty = "n"
)


#---- Restore graphical settings ----

## Restore the original graphical settings
par(old_par)

