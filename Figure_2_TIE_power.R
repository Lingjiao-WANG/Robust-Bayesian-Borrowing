#### Figure 2: TIE, power, overall calibrated power gain, #### 
####           and restricted calibrated power gain       ####

## Source functions
source("case_study_funs.R")

## Source the data of calibrated analysis priors
load("TIE_power.RData")

par(
  mfrow = c(2,3),
  oma = c(4,0,0,0)   # reserve bottom outer margin
)

# --- FPP ----

plot(
  mu_grid,
  TIE_fpp,
  lwd = 2,
  type = "l",
  ylim = c(min(power_gain_fpp,power_gain_fpp_ps),max(1,Power_fpp)),
  xlab = expression(theta[C]),
  ylab = "TIE or Power",
  main = "(A) FPP",
  font.main = 1
)
abline(h=alpha,lty=2)
abline(v=xbar_ch,lty=2)
abline(h=1-beta,lty=2,col=2)
abline(v=xbar_ch+DEL_l,lty=3)
abline(v=xbar_ch+DEL_u,lty=3)
abline(h=0)
lines(mu_grid,Power_fpp,lwd=2,col=2)
lines(mu_grid,power_gain_fpp,lwd=2,col=3,lty=1)
lines(mu_grid,power_gain_fpp_ps,lwd=2,col=3,lty=2)

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

# --- APP ----
# APP_ps
plot(mu_grid,TIE_app_ps,
     lwd = 2,
     type = "l",
     ylim = c(min(power_gain_app_ps,power_gain_app_ps_ps),max(1,Power_app_ps)),
     xlab = expression(theta[C]),
     ylab = "TIE or Power",
     main = expression("(B) " * APP[ps])
)
abline(h=alpha,lty=2)
abline(v=xbar_ch,lty=2)
abline(h=1-beta,lty=2,col=2)
abline(v=xbar_ch+DEL_l,lty=3)
abline(v=xbar_ch+DEL_u,lty=3)
abline(h=0)
lines(mu_grid,Power_app_ps,lwd=2,col=2)
lines(mu_grid,power_gain_app_ps,lwd=2,col=3,lty=1)
lines(mu_grid,power_gain_app_ps_ps,lwd=2,col=3,lty=2)

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

# APP_ss
plot(mu_grid,TIE_app_ss,
     lwd = 2,
     type = "l",
     ylim = c(min(power_gain_app_ss,power_gain_app_ps),max(1,Power_app_ss)),
     xlab = expression(theta[C]),
     ylab = "TIE or Power",
     main = expression("(C) " * APP[ss])
)
abline(h=alpha,lty=2)
abline(v=xbar_ch,lty=2)
abline(h=1-beta,lty=2,col=2)
abline(v=xbar_ch+DEL_l,lty=3)
abline(v=xbar_ch+DEL_u,lty=3)
abline(v=xbar_ch-ss_asy_list_within_ps$DEL_u,lty=3,col=2)
abline(v=mu_grid[max(which(Power_app_ss>=1-beta))],lty=3,col=2)
abline(h=0)
lines(mu_grid,Power_app_ss,lwd=2,col=2)
lines(mu_grid,power_gain_app_ss,lwd=2,col=3,lty=1)
lines(mu_grid,power_gain_app_ss_ss,lwd=2,col=3,lty=2)

axis(
  side = 1,
  at = c(xbar_ch-ss_asy_list_within_ps$DEL_u, mu_grid[max(which(Power_app_ss>=1-beta))]),
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

mu_grid <- seq(
  xbar_ch - 200,
  xbar_ch + 200,
  length.out = 1e3
)

plot(
  mu_grid,TIE_map,
  lwd = 2,
  type = "l",
  ylim = c(min(power_gain_map,power_gain_map_ps),max(1,Power_map)),
  xlab = expression(theta[C]),
  ylab = "TIE or Power",
  main = "(D) MAP Prior",
  font.main = 1
)
abline(h=alpha,lty=2)
abline(v=xbar_ch,lty=2)
abline(h=1-beta,lty=2,col=2)
abline(v=xbar_ch+DEL_l,lty=3)
abline(v=xbar_ch+DEL_u,lty=3)
abline(h=0)
lines(mu_grid,Power_map,lwd=2,col=2)
lines(mu_grid,power_gain_map,lwd=2,col=3,lty=1)
lines(mu_grid,power_gain_map_ps,lwd=2,col=3,lty=2)

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

plot(
  mu_grid,TIE_rmap_1,
  lwd = 2,
  type = "l",
  ylim = c(min(power_gain_rmap_1,power_gain_rmap_1_ps),max(1,Power_rmap_1)),
  xlab = expression(theta[C]),
  ylab = "TIE or Power",
  main = expression(
    paste("(E) ", rMAP[1])
  )
)
abline(h=alpha,lty=2)
abline(v=xbar_ch,lty=2)
abline(h=1-beta,lty=2,col=2)
abline(v=xbar_ch+DEL_l,lty=3)
abline(v=xbar_ch+DEL_u,lty=3)
abline(h=0)
lines(mu_grid,Power_rmap_1,lwd=2,col=2)
lines(mu_grid,power_gain_rmap_1,lwd=2,col=3,lty=1)
lines(mu_grid,power_gain_rmap_1_ps,lwd=2,col=3,lty=2)

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

plot(
  mu_grid,TIE_rmap_2,
  lwd = 2,
  type = "l",
  ylim = c(min(power_gain_rmap_2,power_gain_rmap_2_ps),max(1,Power_rmap_2)),
  xlab = expression(theta[C]),
  ylab = "TIE or Power",
  main = expression(
    paste("(F) ", rMAP[2])
  )
)
abline(h=alpha,lty=2)
abline(v=xbar_ch,lty=2)
abline(h=1-beta,lty=2,col=2)
abline(v=xbar_ch+DEL_l,lty=3)
abline(v=xbar_ch+DEL_u,lty=3)
abline(h=0)
lines(mu_grid,Power_rmap_2,lwd=2,col=2)
lines(mu_grid,power_gain_rmap_2,lwd=2,col=3,lty=1)
lines(mu_grid,power_gain_rmap_2_ps,lwd=2,col=3,lty=2)

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


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend(
  "bottom",
  #inset = c(0.12, -0.35),   # positive moves left for "bottom"
  horiz = TRUE,
  xpd = TRUE,
  cex = 1.1,
  legend = c(
    "TIE: Bayesian Prior",
    "Power: Bayesian Prior",
    "TIE: No Borrowing",
    "Power: No Borrowing",
    "Power Gain",
    "Restricted Power Gain"
  ),
  lwd = c(2, 2, 1, 1, 2, 2),
  lty = c(1, 1, 2, 2, 1, 2),
  col = c(1, 2, 1, 2, 3, 3),
  bty = "n"
)

