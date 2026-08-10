######################################
#### Case study:                  ####
#### Borrow from the control arm  ####
#### Functions                    ####
######################################

#### Fixed power prior ####

#---- Compute the posterior probability threshold p ----
B_thr_f <- function(DEL, alpha, alpha_0,
                    n_T, n_C, se_ch,
                    sigma2, theta_0) {
  # DEL     -- prior-data discrepancy, Delta
  # alpha   -- nominal Type I error rate
  # alpha_0 -- borrowing parameter for the FPP
  # n_T     -- treatment-arm sample size
  # n_C     -- control-arm sample size
  # se_ch   -- standard error of the historical control mean;
  #            se_ch = sigma / sqrt(n_ch)
  # sigma2  -- known outcome variance
  # theta_0 -- superiority margin
  
  prec_ch <- alpha_0 / se_ch^2
  prec_C  <- n_C / sigma2
  prec_post <- prec_ch + prec_C
  
  w_ch <- prec_ch / prec_post
  
  V  <- sigma2 / n_T + 1 / prec_post
  V0 <- sigma2 / n_T + (prec_C / prec_post)^2 * sigma2 / n_C
  
  z_alp <- qnorm(1 - alpha)
  
  pnorm((theta_0 - z_alp * sqrt(V0) - w_ch * DEL) / sqrt(V))
}

#---- Compute the worst case power over Delta in [DEL_l, DEL_u] ----
min_Pow_f <- function(DEL_l, DEL_u,
                      alpha, alpha_0,
                      n_T, n_C, se_ch,
                      xbar_ch,
                      theta0, theta1, beta,
                      sigma2) {
  # DEL_l   -- lower bound of the discrepancy range
  # DEL_u   -- upper bound of the discrepancy range
  # alpha   -- nominal Type I error rate
  # alpha_0 -- borrowing parameter for the FPP
  # n_T     -- treatment-arm sample size
  # n_C     -- control-arm sample size
  # se_ch   -- standard error of the historical control mean;
  #            se_ch = sigma / sqrt(n_ch)
  # xbar_ch -- historical control mean
  # theta0  -- superiority margin
  # theta1  -- treatment effect under the alternative hypothesis
  # beta    -- target Type II error rate, i.e., 1 - power
  # sigma2  -- known outcome variance
  
  B_thr <- B_thr_f(
    DEL     = DEL_u,
    alpha   = alpha,
    alpha_0 = alpha_0,
    n_T     = n_T,
    n_C     = n_C,
    se_ch   = se_ch,
    sigma2  = sigma2,
    theta_0 = theta0
  )
  
  z_B_thr <- qnorm(B_thr)
  
  prec_ch <- alpha_0 / se_ch^2
  prec_C  <- n_C / sigma2
  prec_post <- prec_ch + prec_C
  
  w_ch <- prec_ch / prec_post
  
  V  <- sigma2 / n_T + 1 / prec_post
  V0 <- sigma2 / n_T + (prec_C / prec_post)^2 * sigma2 / n_C
  
  1 - pnorm(((theta0 - theta1) -
               z_B_thr * sqrt(V) -
               w_ch * DEL_l) / sqrt(V0))
}

#---- Compute the difference between the worst-case power over Delta in [DEL_l, DEL_u] and target power ----
Calibration_f <- function(DEL_l, DEL_u,
                          alpha, beta, alpha_0,
                          n_T, n_C, se_ch,
                          xbar_ch,
                          sigma2, theta0, theta1) {
  # DEL_l   -- lower bound of the discrepancy range
  # DEL_u   -- upper bound of the discrepancy range
  # alpha   -- nominal Type I error rate
  # beta    -- target Type II error rate, i.e., 1 - power
  # alpha_0 -- borrowing parameter for the FPP
  # n_T     -- treatment-arm sample size
  # n_C     -- control-arm sample size
  # se_ch   -- standard error of the historical control mean;
  #            se_ch = sigma / sqrt(n_ch)
  # xbar_ch -- historical control mean
  # sigma2  -- known outcome variance
  # theta0  -- superiority margin
  # theta1  -- treatment effect under the alternative hypothesis
  
  min_Pow_f(
    DEL_l   = DEL_l,
    DEL_u   = DEL_u,
    alpha   = alpha,
    alpha_0 = alpha_0,
    n_T     = n_T,
    n_C     = n_C,
    se_ch   = se_ch,
    xbar_ch = xbar_ch,
    theta0  = theta0,
    theta1  = theta1,
    beta    = beta,
    sigma2  = sigma2
  ) - (1 - beta)
}

#---- Compute the probability of rejecting H_0 ----
calc_TIE_power_pp <- function(mu_C,
                              theta_true,
                              xbar_ch,
                              se_ch,
                              n_T,
                              n_C,
                              sigma2,
                              B_thr,
                              alpha_0) {
  # mu_C       -- true control-arm mean
  # theta_true -- true treatment effect
  # xbar_ch    -- historical control mean
  # se_ch      -- standard error of the historical control mean;
  #               se_ch = sigma / sqrt(n_ch)
  # n_T        -- treatment-arm sample size
  # n_C        -- control-arm sample size
  # sigma2     -- known outcome variance
  # B_thr      -- calibrated posterior probability threshold
  # alpha_0    -- borrowing parameter for the FPP
  # theta0     -- superiority margin
  
  prior_var <- se_ch^2 / alpha_0
  
  ## True current treatment mean
  mu_T_true <- mu_C + theta_true
  
  ## Posterior variance for mu_C after observing current control
  post_var_C <- 1 / (1 / prior_var + n_C / sigma2)
  
  ## Posterior mean weight on current control data
  w_C <- (n_C / sigma2) / (1 / prior_var + n_C / sigma2)
  
  ## Posterior mean of mu_C:
  ## E[post_mean_C] = (1-w_C)*xbar_ch + w_C*mu_C
  E_post_mean_C <- (1 - w_C) * xbar_ch + w_C * mu_C
  
  ## Decision statistic:
  ## posterior for theta = mu_T - mu_C
  ## posterior mean = xbar_T - post_mean_C
  ## posterior variance = sigma2/n_T + post_var_C
  post_var_theta <- sigma2 / n_T + post_var_C
  post_sd_theta  <- sqrt(post_var_theta)
  
  ## Rejection rule:
  ## Pr(theta < theta0 | data) > 1 - alpha
  ## equivalent to posterior mean < theta0 - z_{1-alpha} * post_sd
  crit <- theta0 - qnorm(1 - B_thr) * post_sd_theta
  
  ## Sampling variance of posterior mean statistic
  var_stat <- sigma2 / n_T + w_C^2 * sigma2 / n_C
  
  ## Expected posterior mean of theta
  mean_stat <- mu_T_true - E_post_mean_C
  
  ## Probability of rejection
  pnorm((crit - mean_stat) / sqrt(var_stat))
}


#### Adaptive power prior ####

#---- Compute Pr(theta < theta0 | data) ----
ppos_box_vec <- function(xbar_T, xbar_C, xbar_ch,
                         n_T, n_C, sigma2, se_ch,
                         delta0, theta0) {
  
  # xbar_T  -- treatment-arm sample mean
  # xbar_C  -- control-arm sample mean
  # xbar_ch -- historical control mean
  # n_T     -- treatment-arm sample size
  # n_C     -- control-arm sample size
  # sigma2  -- known outcome variance
  # se_ch   -- standard error of the historical control mean;
  #            se_ch = sigma / sqrt(n_ch)
  # delta0  -- borrowing parameter for the APP
  # theta0  -- superiority margin
  
  v_ch <- se_ch^2
  
  se_pred <- sqrt(sigma2 / n_C + v_ch)
  z_obs <- (xbar_C - xbar_ch) / se_pred
  omega0 <- 2 * (1 - pnorm(abs(z_obs)))
  
  n_ch_eff <- sigma2 / v_ch
  a0 <- omega0 ^ ((n_ch_eff / n_C) ^ delta0)
  a0 <- pmin(1, pmax(0, a0))
  
  prec_C <- n_C / sigma2 + a0 / v_ch
  v_C <- 1 / prec_C
  m_C <- v_C * (n_C * xbar_C / sigma2 + a0 * xbar_ch / v_ch)
  
  v_T <- sigma2 / n_T
  
  m_theta <- xbar_T - m_C
  v_theta <- v_T + v_C
  
  pnorm((theta0 - m_theta) / sqrt(v_theta))
}

#---- Compute the posterior probability threshold ----
p_quantile_fast <- function(muC, z_C, z_T, xbar_ch,
                            n_T, n_C, sigma2, se_ch, delta0,
                            theta0, alpha) {
  
  # muC     -- true control-arm mean
  # z_C     -- simulated standard normal values used to generate
  #            the control-arm sample mean xbar_C
  # z_T     -- simulated standard normal values used to generate
  #            the treatment-arm sample mean xbar_T
  # xbar_ch -- historical control mean
  # n_T     -- treatment-arm sample size
  # n_C     -- control-arm sample size
  # sigma2  -- known outcome variance
  # se_ch   -- standard error of the historical control mean;
  #            se_ch = sigma / sqrt(n_ch)
  # delta0  -- borrowing parameter for the APP
  # theta0  -- superiority margin
  # alpha   -- nominal Type I error rate
  
  xbar_C <- muC + sqrt(sigma2 / n_C) * z_C
  xbar_T <- muC + theta0 + sqrt(sigma2 / n_T) * z_T
  
  ppos <- ppos_box_vec(
    xbar_T = xbar_T,
    xbar_C = xbar_C,
    xbar_ch = xbar_ch,
    n_T = n_T,
    n_C = n_C,
    sigma2 = sigma2,
    se_ch = se_ch,
    delta0 = delta0,
    theta0 = theta0
  )
  
  as.numeric(quantile(ppos, alpha, names = FALSE))
}

#---- Compute the probability of rejecting H_0: theta<=theta_0 ----
estimate_OC_cpp_box_p_fast <- function(muC0, n_T, n_C, theta, sigma2,
                                       xbar_ch, se_ch, B_thr, delta0,
                                       z_C, z_T, theta0 = 0) {
  
  # muC0    -- true control-arm mean
  # n_T     -- treatment-arm sample size
  # n_C     -- control-arm sample size
  # theta   -- true treatment effect
  # sigma2  -- known outcome variance
  # xbar_ch -- historical control mean
  # se_ch   -- standard error of the historical control mean;
  #            se_ch = sigma / sqrt(n_ch)
  # B_thr   -- calibrated posterior probability threshold
  # delta0  -- borrowing parameter for the APP
  # z_C     -- simulated standard normal values used to generate
  #            the control-arm sample mean xbar_C
  # z_T     -- simulated standard normal values used to generate
  #            the treatment-arm sample mean xbar_T
  # theta0  -- superiority margin
  
  v_ch <- se_ch^2
  
  ## Generate the current treatment- and control-arm sample means
  xbar_C <- muC0 + sqrt(sigma2 / n_C) * z_C
  xbar_T <- muC0 + theta + sqrt(sigma2 / n_T) * z_T
  
  ## Compute the prior-data compatibility measure
  se_pred <- sqrt(sigma2 / n_C + v_ch)
  z_obs <- (xbar_C - xbar_ch) / se_pred
  omega0 <- 2 * (1 - pnorm(abs(z_obs)))
  
  ## Compute the adaptive borrowing weight
  n_ch_eff <- sigma2 / v_ch
  a0 <- omega0^((n_ch_eff / n_C)^delta0)
  a0 <- pmin(1, pmax(0, a0))
  
  ## Posterior mean and variance of the control-arm mean
  prec_C <- n_C / sigma2 + a0 / v_ch
  v_C <- 1 / prec_C
  m_C <- v_C * (n_C * xbar_C / sigma2 +
                  a0 * xbar_ch / v_ch)
  
  ## Posterior mean and variance of the treatment effect
  v_T <- sigma2 / n_T
  m_theta <- xbar_T - m_C
  v_theta <- v_T + v_C
  
  ## Posterior probability used in the decision rule
  ppos <- pnorm((m_theta - theta0) / sqrt(v_theta))
  
  ## Estimated probability of rejecting H_0
  mean(ppos >= 1 - B_thr)
}

#---- Compute the calibrated posterior probability threshold ----
find_global_min_p_quantile <- function(DEL_l, DEL_u,
                                       z_C, z_T, xbar_ch,
                                       n_T, n_C, sigma2, se_ch,
                                       delta0, theta0,
                                       alpha,
                                       grid_len = 11,
                                       ncore = 8) {
  
  # DEL_l   -- lower bound of the discrepancy range
  # DEL_u   -- upper bound of the discrepancy range
  # z_C     -- simulated standard normal values used to generate
  #            the control-arm sample mean xbar_C
  # z_T     -- simulated standard normal values used to generate
  #            the treatment-arm sample mean xbar_T
  # xbar_ch -- historical control mean
  # n_T     -- treatment-arm sample size
  # n_C     -- control-arm sample size
  # sigma2  -- known outcome variance
  # se_ch   -- standard error of the historical control mean;
  #            se_ch = sigma / sqrt(n_ch)
  # delta0  -- borrowing parameter for the APP
  # theta0  -- superiority margin
  # alpha   -- nominal Type I error rate
  # grid_len -- number of grid points used in the initial search
  # ncore   -- number of CPU cores used for parallel computation
  
  ## Define the range of true control-arm means
  lower <- xbar_ch + DEL_l
  upper <- xbar_ch + DEL_u
  
  ## Construct the initial search grid
  mu_grid <- seq(lower, upper, length.out = grid_len)
  
  ## Evaluate the posterior probability threshold over the grid
  if (ncore > 1) {
    
    cl <- parallel::makeCluster(ncore)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterExport(
      cl,
      c(
        "p_quantile_fast",
        "ppos_box_vec",
        "z_C", "z_T",
        "xbar_ch",
        "n_T", "n_C",
        "sigma2", "se_ch",
        "delta0",
        "theta0",
        "alpha"
      ),
      envir = environment()
    )
    
    q_grid <- unlist(
      parallel::parLapply(cl, mu_grid, function(mu) {
        p_quantile_fast(
          muC = mu,
          z_C = z_C,
          z_T = z_T,
          xbar_ch = xbar_ch,
          n_T = n_T,
          n_C = n_C,
          sigma2 = sigma2,
          se_ch = se_ch,
          delta0 = delta0,
          theta0 = theta0,
          alpha = alpha
        )
      })
    )
    
  } else {
    
    q_grid <- sapply(mu_grid, function(mu) {
      p_quantile_fast(
        muC = mu,
        z_C = z_C,
        z_T = z_T,
        xbar_ch = xbar_ch,
        n_T = n_T,
        n_C = n_C,
        sigma2 = sigma2,
        se_ch = se_ch,
        delta0 = delta0,
        theta0 = theta0,
        alpha = alpha
      )
    })
  }
  
  ## Identify local minima on the initial grid
  idx_local <- which(
    q_grid[2:(length(q_grid) - 1)] <
      q_grid[1:(length(q_grid) - 2)] &
      q_grid[2:(length(q_grid) - 1)] <
      q_grid[3:length(q_grid)]
  ) + 1
  
  ## Include the two boundary points as candidate minima
  idx_candidate <- unique(
    c(1, length(mu_grid), idx_local)
  )
  
  ## Refine each candidate minimum using one-dimensional optimisation
  local_res <- lapply(idx_candidate, function(j) {
    
    if (j == 1 || j == length(mu_grid)) {
      return(
        list(
          minimum = mu_grid[j],
          objective = q_grid[j]
        )
      )
    }
    
    optimise(
      p_quantile_fast,
      interval = c(mu_grid[j - 1], mu_grid[j + 1]),
      z_C = z_C,
      z_T = z_T,
      xbar_ch = xbar_ch,
      n_T = n_T,
      n_C = n_C,
      sigma2 = sigma2,
      se_ch = se_ch,
      delta0 = delta0,
      theta0 = theta0,
      alpha = alpha,
      maximum = FALSE,
      tol = 1e-3
    )
  })
  
  ## Select the global minimum among all candidate minima
  local_vals <- sapply(local_res, function(x) x$objective)
  local_mu   <- sapply(local_res, function(x) x$minimum)
  
  k <- which.min(local_vals)
  
  list(
    p_0 = local_vals[k],
    muC_at_min = local_mu[k],
    mu_grid = mu_grid,
    q_grid = q_grid,
    all_candidate_mu = local_mu,
    all_candidate_val = local_vals
  )
}

#---- Compute the difference between the worst-case power over Delta in [DEL_l, DEL_u] and target power ----
cal_f_ad <- function(DEL_l, DEL_u,
                     n_T, n_C,
                     alpha, beta,
                     theta0, theta1,
                     sigma2,
                     xbar_ch, se_ch, delta0,
                     z_C, z_T,
                     ncore0 = 8) {
  
  # DEL_l   -- lower bound of the discrepancy range
  # DEL_u   -- upper bound of the discrepancy range
  # n_T     -- treatment-arm sample size
  # n_C     -- control-arm sample size
  # alpha   -- nominal Type I error rate
  # beta    -- target Type II error rate, i.e., 1 - power
  # theta0  -- superiority margin
  # theta1  -- treatment effect under the alternative hypothesis
  # sigma2  -- known outcome variance
  # xbar_ch -- historical control mean
  # se_ch   -- standard error of the historical control mean;
  #            se_ch = sigma / sqrt(n_ch)
  # delta0  -- borrowing parameter for the APP
  # z_C     -- simulated standard normal values used to generate
  #            the control-arm sample mean xbar_C
  # z_T     -- simulated standard normal values used to generate
  #            the treatment-arm sample mean xbar_T
  # ncore0  -- number of CPU cores used for parallel computation
  
  ## Compute the calibrated posterior probability threshold
  p_0 <- find_global_min_p_quantile(
    DEL_l = DEL_l,
    DEL_u = DEL_u,
    z_C = z_C,
    z_T = z_T,
    xbar_ch = xbar_ch,
    n_T = n_T,
    n_C = n_C,
    sigma2 = sigma2,
    se_ch = se_ch,
    delta0 = delta0,
    theta0 = theta0,
    alpha = alpha,
    grid_len = 11,
    ncore = ncore0
  )$p_0
  
  ## Construct a grid of true control-arm means over the discrepancy range
  muC0_v <- seq(
    from = DEL_l + xbar_ch,
    to = DEL_u + xbar_ch,
    length.out = 10
  )
  
  ## Define the power function for a given true control-arm mean
  obj <- function(muC0) {
    estimate_OC_cpp_box_p_fast(
      muC0 = muC0,
      n_T = n_T,
      n_C = n_C,
      theta = theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      B_thr = p_0,
      delta0 = delta0,
      z_C = z_C,
      z_T = z_T,
      theta0 = theta0
    )
  }
  
  ## Evaluate power over the initial grid
  min_pow_v <- sapply(muC0_v, obj)
  
  ## Identify local minima on the initial grid
  idx_local <- which(
    min_pow_v[2:(length(min_pow_v) - 1)] <
      min_pow_v[1:(length(min_pow_v) - 2)] &
      min_pow_v[2:(length(min_pow_v) - 1)] <
      min_pow_v[3:length(min_pow_v)]
  ) + 1
  
  ## Use the minimum grid value as an initial candidate
  candidate_vals <- min(min_pow_v)
  candidate_mu   <- muC0_v[which.min(min_pow_v)]
  
  ## Refine each local minimum using one-dimensional optimisation
  if (length(idx_local) > 0) {
    
    local_res <- lapply(idx_local, function(j) {
      optimise(
        obj,
        interval = c(muC0_v[j - 1], muC0_v[j + 1]),
        maximum = FALSE
      )
    })
    
    local_vals <- sapply(local_res, function(x) x$objective)
    local_mu   <- sapply(local_res, function(x) x$minimum)
    
    all_vals <- c(candidate_vals, local_vals)
    all_mu   <- c(candidate_mu, local_mu)
    
    k <- which.min(all_vals)
    
    min_pow_o <- all_vals[k]
    muC0_opt  <- all_mu[k]
    
  } else {
    
    min_pow_o <- candidate_vals
    muC0_opt  <- candidate_mu
  }
  
  ## Return the power difference and calibrated quantities
  list(
    res = min_pow_o - (1 - beta),
    p_0 = p_0,
    min_pow = min_pow_o,
    muC0_opt = muC0_opt
  )
}

#---- Compute the calibrated borrowing parameter delta0 for the APP ----
find_delta0_roots_cal_f_ad <- function(DEL_l, DEL_u,
                                       n_T, n_C,
                                       alpha, beta,
                                       theta0, theta1,
                                       sigma2,
                                       xbar_ch, se_ch,
                                       z_C, z_T,
                                       delta0_min = 0,
                                       delta0_max = 3,
                                       n_grid = 201,
                                       ncore0 = 8,
                                       tol = 1e-6) {
  
  # DEL_l      -- lower bound of the discrepancy range
  # DEL_u      -- upper bound of the discrepancy range
  # n_T        -- treatment-arm sample size
  # n_C        -- control-arm sample size
  # alpha      -- nominal Type I error rate
  # beta       -- target Type II error rate, i.e., 1 - power
  # theta0     -- superiority margin
  # theta1     -- treatment effect under the alternative hypothesis
  # sigma2     -- known outcome variance
  # xbar_ch    -- historical control mean
  # se_ch      -- standard error of the historical control mean;
  #               se_ch = sigma / sqrt(n_ch)
  # z_C        -- simulated standard normal values used to generate
  #               the control-arm sample mean xbar_C
  # z_T        -- simulated standard normal values used to generate
  #               the treatment-arm sample mean xbar_T
  # delta0_min -- lower bound of the search range for delta0
  # delta0_max -- upper bound of the search range for delta0
  # n_grid     -- number of grid points used in the search
  # ncore0     -- number of CPU cores used for parallel computation
  # tol        -- numerical tolerance for identifying a root
  
  ## Construct the search grid for delta0
  delta_grid <- seq(
    delta0_min,
    delta0_max,
    length.out = n_grid
  )
  
  delta_grid <- sort(unique(delta_grid))
  n_delta <- length(delta_grid)
  
  if (n_delta < 2L) {
    stop("delta_grid must contain at least two distinct values.")
  }
  
  ## Evaluate the calibration function for a given delta0
  evaluate_delta0 <- function(delta0) {
    cal_f_ad(
      DEL_l = DEL_l,
      DEL_u = DEL_u,
      n_T = n_T,
      n_C = n_C,
      alpha = alpha,
      beta = beta,
      theta0 = theta0,
      theta1 = theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      delta0 = delta0,
      z_C = z_C,
      z_T = z_T,
      ncore0 = ncore0
    )
  }
  
  ## Objective function: worst-case power minus target power
  f_obj <- function(delta0) {
    evaluate_delta0(delta0)$res
  }
  
  ## Store only the grid points that are actually evaluated
  delta_checked <- numeric(0)
  f_checked <- numeric(0)
  
  # ===============================================================
  # Evaluate the smallest value of delta0
  # ===============================================================
  
  delta_prev <- delta_grid[1L]
  selected_info <- evaluate_delta0(delta_prev)
  f_prev <- selected_info$res
  
  if (!is.finite(f_prev)) {
    stop(
      "A non-finite objective value is obtained at delta0 = ",
      delta_prev,
      "."
    )
  }
  
  delta_checked <- c(delta_checked, delta_prev)
  f_checked <- c(f_checked, f_prev)
  
  # ===============================================================
  # Case 1:
  # The objective at the smallest delta0 is already non-negative
  # ===============================================================
  
  if (f_prev >= -tol) {
    
    delta0_use <- delta_prev
    
    root_case <- if (abs(f_prev) <= tol) {
      "minimum_grid_value_is_zero"
    } else {
      "minimum_grid_value_is_nonnegative"
    }
    
    return(list(
      has_root = abs(f_prev) <= tol,
      root_case = root_case,
      delta0 = delta0_use,
      root = if (abs(f_prev) <= tol) {
        delta0_use
      } else {
        numeric(0)
      },
      root_type = if (abs(f_prev) <= tol) {
        "grid_zero"
      } else {
        "boundary_nonnegative"
      },
      root_info = selected_info,
      f_at_root = selected_info$res,
      crossing_interval = NULL,
      transition_index = 1L,
      delta_checked = delta_checked,
      f_checked = f_checked,
      delta_grid = delta_grid,
      res = paste0(
        "The objective at delta0_min is non-negative. ",
        "Therefore, delta0 = ",
        round(delta0_use, 6),
        " is selected."
      )
    ))
  }
  
  ## At this point, the objective is negative at delta0_min
  
  # ===============================================================
  # Sequentially evaluate the remaining grid points
  # ===============================================================
  
  for (i in 2L:n_delta) {
    
    delta_curr <- delta_grid[i]
    selected_info <- evaluate_delta0(delta_curr)
    f_curr <- selected_info$res
    
    if (!is.finite(f_curr)) {
      stop(
        "A non-finite objective value is obtained at delta0 = ",
        delta_curr,
        "."
      )
    }
    
    delta_checked <- c(delta_checked, delta_curr)
    f_checked <- c(f_checked, f_curr)
    
    # =============================================================
    # Case 2:
    # The current grid value is numerically equal to zero
    # =============================================================
    
    if (abs(f_curr) <= tol) {
      
      delta0_use <- delta_curr
      
      return(list(
        has_root = TRUE,
        root_case = "negative_to_numerical_zero",
        delta0 = delta0_use,
        root = delta0_use,
        root_type = "grid_zero",
        root_info = selected_info,
        f_at_root = selected_info$res,
        crossing_interval = c(delta_prev, delta_curr),
        transition_index = i,
        delta_checked = delta_checked,
        f_checked = f_checked,
        delta_grid = delta_grid,
        res = paste0(
          "The objective starts negative. The first grid value ",
          "at which the objective is numerically zero is selected: ",
          "delta0 = ",
          round(delta0_use, 6),
          "."
        )
      ))
    }
    
    # =============================================================
    # Case 3:
    # Strict negative-to-positive crossing
    # =============================================================
    
    if (f_prev < -tol && f_curr > tol) {
      
      delta0_use <- uniroot(
        f = f_obj,
        interval = c(delta_prev, delta_curr),
        tol = tol
      )$root
      
      selected_info <- evaluate_delta0(delta0_use)
      
      return(list(
        has_root = TRUE,
        root_case = "negative_to_positive_crossing",
        delta0 = delta0_use,
        root = delta0_use,
        root_type = "uniroot",
        root_info = selected_info,
        f_at_root = selected_info$res,
        crossing_interval = c(delta_prev, delta_curr),
        transition_index = i - 1L,
        delta_checked = delta_checked,
        f_checked = f_checked,
        delta_grid = delta_grid,
        res = paste0(
          "The first interval over which the objective changes ",
          "from negative to positive is identified. ",
          "The corresponding root is selected: delta0 = ",
          round(delta0_use, 6),
          "."
        )
      ))
    }
    
    ## Move to the next interval
    delta_prev <- delta_curr
    f_prev <- f_curr
  }
  
  # ===============================================================
  # Case 4:
  # No negative-to-non-negative transition over the search range
  # ===============================================================
  
  list(
    has_root = FALSE,
    root_case = "no_negative_to_nonnegative_transition",
    delta0 = NA_real_,
    root = NA_real_,
    root_type = NA_character_,
    root_info = NULL,
    f_at_root = NA_real_,
    crossing_interval = NULL,
    transition_index = NA_integer_,
    delta_checked = delta_checked,
    f_checked = f_checked,
    delta_grid = delta_grid,
    res = paste0(
      "The objective is negative at delta0_min, and no transition ",
      "from negative to non-negative is found over the specified ",
      "delta0 range."
    )
  )
}

#---- Cache evaluations of the calibration function ----
cal_cache <- new.env(parent = emptyenv())

#---- Evaluate the calibration function with caching ----
cal_f_cached <- function(dl, du) {
  
  # dl -- lower bound of the discrepancy range
  # du -- upper bound of the discrepancy range
  
  ## Create a unique cache key for the discrepancy range
  key <- paste(round(dl, 6), round(du, 6), sep = "_")
  
  ## Return the cached result if it has already been evaluated
  if (exists(key, envir = cal_cache)) {
    return(get(key, envir = cal_cache))
  }
  
  ## Evaluate the calibration function
  val <- cal_f_ad(
    DEL_l = dl,
    DEL_u = du,
    n_T = n_T,
    n_C = n_C,
    alpha = alpha,
    beta = beta,
    theta0 = theta0,
    theta1 = theta1,
    sigma2 = sigma2,
    xbar_ch = xbar_ch,
    se_ch = se_ch,
    delta0 = delta0,
    z_C = z_C,
    z_T = z_T
  )
  
  ## Store the result in the cache
  assign(key, val, envir = cal_cache)
  
  val
}

#---- Find the widest feasible discrepancy interval ----
find_widest_fast <- function(DEL_l_min, DEL_l_max,
                             DEL_u_min, DEL_u_max,
                             delta0,
                             n_l = 15,
                             n_u = 15,
                             ncore = 8) {
  
  # DEL_l_min -- lower bound of the search range for DEL_l
  # DEL_l_max -- upper bound of the search range for DEL_l
  # DEL_u_min -- lower bound of the search range for DEL_u
  # DEL_u_max -- upper bound of the search range for DEL_u
  # delta0    -- borrowing parameter for the APP
  # n_l       -- number of grid points used to search over DEL_l
  # n_u       -- number of grid points used to search over DEL_u
  # ncore     -- number of CPU cores used for parallel computation
  
  ## Construct the search grid for the lower discrepancy bound
  dl_grid <- seq(
    DEL_l_min,
    DEL_l_max,
    length.out = n_l
  )
  
  ## Set up parallel computation
  cl <- parallel::makeCluster(ncore)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterExport(
    cl,
    c(
      "cal_f_ad",
      "find_global_min_p_quantile",
      "p_quantile_fast",
      "estimate_OC_cpp_box_p_fast",
      "ppos_box_vec",
      "n_T", "n_C",
      "alpha", "beta",
      "theta0", "theta1",
      "sigma2", "xbar_ch", "se_ch",
      "delta0", "z_C", "z_T",
      "DEL_u_min", "DEL_u_max", "n_u"
    ),
    envir = environment()
  )
  
  ## For each lower bound, find feasible upper bounds
  res_list <- parallel::parLapply(cl, dl_grid, function(dl) {
    
    ## Construct the search grid for the upper discrepancy bound
    du_grid <- seq(
      max(dl, DEL_u_min),
      DEL_u_max,
      length.out = n_u
    )
    
    ## Evaluate the calibration objective over the upper-bound grid
    f_grid <- sapply(du_grid, function(du) {
      cal_f_ad(
        DEL_l = dl,
        DEL_u = du,
        n_T = n_T,
        n_C = n_C,
        alpha = alpha,
        beta = beta,
        theta0 = theta0,
        theta1 = theta1,
        sigma2 = sigma2,
        xbar_ch = xbar_ch,
        se_ch = se_ch,
        delta0 = delta0,
        z_C = z_C,
        z_T = z_T,
        ncore0 = 1
      )$res
    })
    
    ## Identify intervals containing a root
    idx <- which(
      f_grid[-length(f_grid)] * f_grid[-1] <= 0
    )
    
    if (length(idx) == 0) {
      return(NULL)
    }
    
    ## Refine all identified roots
    roots <- sapply(idx, function(j) {
      uniroot(
        function(du) {
          cal_f_ad(
            DEL_l = dl,
            DEL_u = du,
            n_T = n_T,
            n_C = n_C,
            alpha = alpha,
            beta = beta,
            theta0 = theta0,
            theta1 = theta1,
            sigma2 = sigma2,
            xbar_ch = xbar_ch,
            se_ch = se_ch,
            delta0 = delta0,
            z_C = z_C,
            z_T = z_T,
            ncore0 = 1
          )$res
        },
        interval = c(du_grid[j], du_grid[j + 1]),
        tol = 1e-3
      )$root
    })
    
    ## Select the largest feasible upper bound
    du_star <- max(roots)
    
    ## Evaluate the operating characteristics at the selected interval
    info <- cal_f_ad(
      DEL_l = dl,
      DEL_u = du_star,
      n_T = n_T,
      n_C = n_C,
      alpha = alpha,
      beta = beta,
      theta0 = theta0,
      theta1 = theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      delta0 = delta0,
      z_C = z_C,
      z_T = z_T,
      ncore0 = 1
    )
    
    data.frame(
      DEL_l = dl,
      DEL_u = du_star,
      width = du_star - dl,
      p_0 = info$p_0,
      min_pow = info$min_pow,
      muC0_opt = info$muC0_opt,
      res = info$res
    )
  })
  
  ## Combine all feasible intervals
  res_df <- do.call(rbind, res_list)
  
  if (is.null(res_df) || nrow(res_df) == 0) {
    stop("No feasible interval found.")
  }
  
  ## Return the feasible interval with the largest width
  res_df[which.max(res_df$width), ]
}

#---- Compute the probability of rejecting H_0: theta>=theta_0 ----
estimate_OC_cpp_box_p_fast_1 <- function(muC0, n_T, n_C, theta, sigma2,
                                         xbar_ch, se_ch, B_thr, delta0,
                                         z_C, z_T, theta0 = 0) {
  
  # muC0    -- true control-arm mean
  # n_T     -- treatment-arm sample size
  # n_C     -- control-arm sample size
  # theta   -- true treatment effect
  # sigma2  -- known outcome variance
  # xbar_ch -- historical control mean
  # se_ch   -- standard error of the historical control mean;
  #            se_ch = sigma / sqrt(n_ch)
  # B_thr   -- calibrated posterior probability threshold
  # delta0  -- borrowing parameter for the APP
  # z_C     -- simulated standard normal values used to generate
  #            the control-arm sample mean xbar_C
  # z_T     -- simulated standard normal values used to generate
  #            the treatment-arm sample mean xbar_T
  # theta0  -- superiority margin
  
  v_ch <- se_ch^2
  
  ## Generate the current control- and treatment-arm sample means
  xbar_C <- muC0 + sqrt(sigma2 / n_C) * z_C
  xbar_T <- muC0 + theta + sqrt(sigma2 / n_T) * z_T
  
  ## Compute the prior-data compatibility measure
  se_pred <- sqrt(sigma2 / n_C + v_ch)
  z_obs <- (xbar_C - xbar_ch) / se_pred
  omega0 <- 2 * (1 - pnorm(abs(z_obs)))
  
  ## Compute the adaptive borrowing weight
  n_ch_eff <- sigma2 / v_ch
  a0 <- omega0^((n_ch_eff / n_C)^delta0)
  a0 <- pmin(1, pmax(0, a0))
  
  ## Posterior mean and variance of the control-arm mean
  prec_C <- n_C / sigma2 + a0 / v_ch
  v_C <- 1 / prec_C
  
  m_C <- v_C * (
    n_C * xbar_C / sigma2 +
      a0 * xbar_ch / v_ch
  )
  
  ## Posterior mean and variance of the treatment effect
  v_T <- sigma2 / n_T
  m_theta <- xbar_T - m_C
  v_theta <- v_T + v_C
  
  ## Posterior probability Pr(theta < theta0 | data)
  pneg <- pnorm(
    (theta0 - m_theta) / sqrt(v_theta)
  )
  
  ## Estimated probability of rejecting H_0
  mean(pneg >= 1 - B_thr)
}

#---- Compute the EHSS as a function of xbar_C under the APP ----
ESS_adaptive_power_prior_xbar <- function(xbar_C_grid,
                                          xbar_ch,
                                          n_C,
                                          sigma2,
                                          se_ch,
                                          delta0,
                                          weights = NULL) {
  
  # xbar_C_grid -- grid of current control-arm sample means
  # xbar_ch     -- historical control mean
  # n_C         -- control-arm sample size
  # sigma2      -- known outcome variance
  # se_ch       -- standard error of the historical control mean;
  #                se_ch = sigma / sqrt(n_ch)
  # delta0      -- borrowing parameter for the APP
  # weights     -- optional weights used to compute the mean EHSS
  
  ## Historical control variance and corresponding effective sample size
  v_ch <- se_ch^2
  n_ch_eff <- sigma2 / v_ch
  
  ## Compute the prior-data compatibility measure
  se_pred <- sqrt(sigma2 / n_C + v_ch)
  z_xbar_C <- (xbar_C_grid - xbar_ch) / se_pred
  omega0 <- 2 * (1 - pnorm(abs(z_xbar_C)))
  
  ## Compute the adaptive borrowing weight
  a0 <- omega0^((n_ch_eff / n_C)^delta0)
  a0 <- pmin(1, pmax(0, a0))
  
  ## Compute the EHSS over the grid of current control-arm means
  ESS_grid <- a0 * n_ch_eff
  
  ## Specify weights for computing the mean EHSS
  if (is.null(weights)) {
    weights <- rep(1 / length(xbar_C_grid), length(xbar_C_grid))
  } else {
    weights <- weights / sum(weights)
  }
  
  ## Return the EHSS summaries
  list(
    mean_ESS = sum(weights * ESS_grid),
    max_ESS = max(ESS_grid),
    xbar_C_at_max_ESS = xbar_C_grid[which.max(ESS_grid)],
    ESS_grid = ESS_grid,
    a0_grid = a0,
    xbar_C_grid = xbar_C_grid,
    n_ch_eff = n_ch_eff
  )
}

#---- Compute the expected EHSS given theta_C under the APP ----
ESS_adaptive_power_prior_muC <- function(mu_grid,
                                         xbar_ch,
                                         n_C,
                                         sigma2,
                                         se_ch,
                                         delta0,
                                         weights = NULL,
                                         z_grid = seq(-5, 5, length.out = 1001)) {
  
  # mu_grid -- grid of true control-arm means
  # xbar_ch -- historical control mean
  # n_C     -- control-arm sample size
  # sigma2  -- known outcome variance
  # se_ch   -- standard error of the historical control mean;
  #            se_ch = sigma / sqrt(n_ch)
  # delta0  -- borrowing parameter for the APP
  # weights -- optional weights used to compute the mean expected EHSS
  # z_grid  -- standard normal grid used for numerical integration
  
  ## Historical control variance and corresponding effective sample size
  v_ch <- se_ch^2
  n_ch_eff <- sigma2 / v_ch
  
  se_pred <- sqrt(sigma2 / n_C + v_ch)
  
  ## Compute the EHSS for a given observed control-arm sample mean
  ESS_one_xbarC <- function(xbar_C) {
    
    z_obs <- (xbar_C - xbar_ch) / se_pred
    
    omega0 <- 2 * (1 - pnorm(abs(z_obs)))
    
    a0 <- omega0^((n_ch_eff / n_C)^delta0)
    a0 <- pmin(1, pmax(0, a0))
    
    a0 * n_ch_eff
  }
  
  ## Construct weights for numerical integration over xbar_C | mu_C
  dz <- z_grid[2] - z_grid[1]
  z_w <- dnorm(z_grid) * dz
  z_w <- z_w / sum(z_w)
  
  ## Compute E_{xbar_C | mu_C}[EHSS(xbar_C)] for each true control-arm mean
  ESS_grid <- sapply(mu_grid, function(mu_C) {
    
    xbar_C_grid <- mu_C + sqrt(sigma2 / n_C) * z_grid
    
    ESS_xbarC <- sapply(xbar_C_grid, ESS_one_xbarC)
    
    sum(z_w * ESS_xbarC)
  })
  
  ## Specify weights for computing the mean expected EHSS
  if (is.null(weights)) {
    weights <- rep(1 / length(mu_grid), length(mu_grid))
  } else {
    weights <- weights / sum(weights)
  }
  
  ## Return the expected EHSS summaries
  list(
    mean_ESS = sum(weights * ESS_grid),
    max_ESS = max(ESS_grid),
    muC_at_max_ESS = mu_grid[which.max(ESS_grid)],
    ESS_grid = ESS_grid,
    n_ch_eff = n_ch_eff,
    z_grid = z_grid,
    z_weights = z_w
  )
}


#### MAP prior ####

#---- Construct the tau grid and integration weights under a half-normal prior ----
make_tau_grid <- function(psi, n_tau = 1001, prob = 0.999) {
  
  # psi   -- scale parameter of the half-normal prior for tau
  # n_tau -- number of grid points used for numerical integration
  # prob  -- cumulative probability used to determine the upper bound of tau
  
  ## Determine the upper bound of the tau grid
  tau_max <- psi * qnorm((prob + 1) / 2)
  
  ## Construct the tau grid
  tau <- seq(0, tau_max, length.out = n_tau)
  
  ## Compute integration interval widths
  d_tau <- rep(NA_real_, n_tau)
  d_tau[1] <- tau[2] - tau[1]
  d_tau[n_tau] <- tau[n_tau] - tau[n_tau - 1]
  d_tau[2:(n_tau - 1)] <-
    (tau[3:n_tau] - tau[1:(n_tau - 2)]) / 2
  
  ## Evaluate the half-normal density over the tau grid
  dens <- sqrt(2 / pi) / psi *
    exp(-tau^2 / (2 * psi^2))
  
  ## Compute and normalise the numerical integration weights
  w <- dens * d_tau
  w <- w / sum(w)
  
  list(
    tau = tau,
    w = w
  )
}


#---- Compute Pr(theta < theta0 | data) under the MAP prior ----
ppos_map_vec <- function(xbar_T, xbar_C, xbar_ch,
                         n_T, n_C, sigma2, se_ch, psi, theta0,
                         n_tau = 1001) {
  
  # xbar_T  -- treatment-arm sample mean
  # xbar_C  -- control-arm sample mean
  # xbar_ch -- historical control mean
  # n_T     -- treatment-arm sample size
  # n_C     -- control-arm sample size
  # sigma2  -- known outcome variance
  # se_ch   -- standard error of the historical control mean;
  #            se_ch = sigma / sqrt(n_ch)
  # psi     -- scale parameter of the half-normal prior for tau
  # theta0  -- superiority margin
  # n_tau   -- number of grid points used for numerical integration over tau
  
  ## Construct the tau grid and prior integration weights
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w0 <- tg$w
  
  v_ch <- se_ch^2
  
  ## Prior variance of the current control mean conditional on tau
  v0_tau <- v_ch + 2 * tau^2
  
  ## Posterior variance contribution from the treatment arm
  vT <- sigma2 / n_T
  
  ## Predictive distribution of the current control mean conditional on tau
  pred_var <- outer(
    rep(1, length(xbar_C)),
    sigma2 / n_C + v0_tau
  )
  
  pred_mean <- matrix(
    xbar_ch,
    nrow = length(xbar_C),
    ncol = length(tau)
  )
  
  ## Compute posterior weights for tau
  log_w <- sweep(
    dnorm(
      xbar_C,
      mean = pred_mean,
      sd = sqrt(pred_var),
      log = TRUE
    ),
    2,
    log(w0),
    "+"
  )
  
  ## Stabilise and normalise the posterior weights
  log_w <- log_w - apply(log_w, 1, max)
  w_post <- exp(log_w)
  w_post <- w_post / rowSums(w_post)
  
  ## Posterior mean and variance of the control-arm mean conditional on tau
  prec_C <- n_C / sigma2
  prec_0 <- 1 / v0_tau
  
  v_C_tau <- 1 / outer(
    rep(1, length(xbar_C)),
    prec_C + prec_0
  )
  
  m_C_tau <- v_C_tau * (
    xbar_C * prec_C +
      outer(
        rep(1, length(xbar_C)),
        xbar_ch * prec_0
      )
  )
  
  ## Posterior mean and variance of the treatment effect conditional on tau
  m_theta_tau <- xbar_T - m_C_tau
  v_theta_tau <- vT + v_C_tau
  
  ## Integrate Pr(theta < theta0 | data, tau) over the posterior distribution of tau
  rowSums(
    w_post *
      pnorm((theta0 - m_theta_tau) / sqrt(v_theta_tau))
  )
}

#---- Compute the calibrated posterior probability threshold given theta_C under the MAP prior ----
p_quantile_fast_map <- function(muC, z_C, z_T, xbar_ch,
                                n_T, n_C, sigma2, se_ch, psi,
                                theta0, alpha,
                                n_tau = 1001) {
  
  # muC     -- true control-arm mean
  # z_C     -- simulated standard normal values used to generate
  #            the control-arm sample mean xbar_C
  # z_T     -- simulated standard normal values used to generate
  #            the treatment-arm sample mean xbar_T
  # xbar_ch -- historical control mean
  # n_T     -- treatment-arm sample size
  # n_C     -- control-arm sample size
  # sigma2  -- known outcome variance
  # se_ch   -- standard error of the historical control mean;
  #            se_ch = sigma / sqrt(n_ch)
  # psi     -- scale parameter of the half-normal prior for tau
  # theta0  -- superiority margin
  # alpha   -- nominal Type I error rate
  # n_tau   -- number of grid points used for numerical integration over tau
  
  ## Generate the current control- and treatment-arm sample means under H_0
  xbar_C <- muC + sqrt(sigma2 / n_C) * z_C
  xbar_T <- muC + theta0 + sqrt(sigma2 / n_T) * z_T
  
  ## Compute the posterior probability Pr(theta < theta0 | data)
  ppos <- ppos_map_vec(
    xbar_T = xbar_T,
    xbar_C = xbar_C,
    xbar_ch = xbar_ch,
    n_T = n_T,
    n_C = n_C,
    sigma2 = sigma2,
    se_ch = se_ch,
    psi = psi,
    theta0 = theta0,
    n_tau = n_tau
  )
  
  ## Calibrate the posterior probability threshold using the alpha quantile
  as.numeric(
    quantile(ppos, alpha, names = FALSE)
  )
}


#---- Compute the probability of rejecting H_0: theta<=theta_0 under the MAP prior ----
estimate_OC_map_p_fast <- function(muC0, n_T, n_C, theta, sigma2,
                                   xbar_ch, se_ch, B_thr, psi,
                                   z_C, z_T,
                                   theta0 = 0,
                                   n_tau = 1001) {
  
  # muC0    -- true control-arm mean
  # n_T     -- treatment-arm sample size
  # n_C     -- control-arm sample size
  # theta   -- true treatment effect
  # sigma2  -- known outcome variance
  # xbar_ch -- historical control mean
  # se_ch   -- standard error of the historical control mean;
  #            se_ch = sigma / sqrt(n_ch)
  # B_thr   -- calibrated posterior probability threshold
  # psi     -- scale parameter of the half-normal prior for tau
  # z_C     -- simulated standard normal values used to generate
  #            the control-arm sample mean xbar_C
  # z_T     -- simulated standard normal values used to generate
  #            the treatment-arm sample mean xbar_T
  # theta0  -- superiority margin
  # n_tau   -- number of grid points used for numerical integration over tau
  
  ## Generate the current control- and treatment-arm sample means
  xbar_C <- muC0 + sqrt(sigma2 / n_C) * z_C
  xbar_T <- muC0 + theta + sqrt(sigma2 / n_T) * z_T
  
  ## Construct the tau grid and prior integration weights
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w0 <- tg$w
  
  v_ch <- se_ch^2
  
  ## Prior variance of the current control mean conditional on tau
  v0_tau <- v_ch + 2 * tau^2
  vT <- sigma2 / n_T
  
  ## Predictive distribution of the current control mean conditional on tau
  pred_var <- outer(
    rep(1, length(xbar_C)),
    sigma2 / n_C + v0_tau
  )
  
  pred_mean <- matrix(
    xbar_ch,
    nrow = length(xbar_C),
    ncol = length(tau)
  )
  
  ## Compute posterior weights for tau
  log_w <- sweep(
    dnorm(
      xbar_C,
      mean = pred_mean,
      sd = sqrt(pred_var),
      log = TRUE
    ),
    2,
    log(w0),
    "+"
  )
  
  ## Stabilise and normalise the posterior weights
  log_w <- log_w - apply(log_w, 1, max)
  w_post <- exp(log_w)
  w_post <- w_post / rowSums(w_post)
  
  ## Posterior mean and variance of the control-arm mean conditional on tau
  prec_C <- n_C / sigma2
  prec_0 <- 1 / v0_tau
  
  v_C_tau <- 1 / outer(
    rep(1, length(xbar_C)),
    prec_C + prec_0
  )
  
  m_C_tau <- v_C_tau * (
    xbar_C * prec_C +
      outer(
        rep(1, length(xbar_C)),
        xbar_ch * prec_0
      )
  )
  
  ## Posterior mean and variance of the treatment effect conditional on tau
  m_theta_tau <- xbar_T - m_C_tau
  v_theta_tau <- vT + v_C_tau
  
  ## Integrate the posterior probability over the posterior distribution of tau
  ppos <- rowSums(
    w_post *
      pnorm((m_theta_tau - theta0) / sqrt(v_theta_tau))
  )
  
  ## Estimated probability of rejecting H_0
  mean(ppos >= 1 - B_thr)
}

#---- Compute the probability of rejecting H_0: theta>=theta_0 under the MAP prior ----
estimate_OC_map_p_fast_1 <- function(muC0, n_T, n_C, theta, sigma2,
                                     xbar_ch, se_ch, B_thr, psi,
                                     z_C, z_T,
                                     theta0 = 0,
                                     n_tau = 1001) {
  
  # muC0    -- true control-arm mean
  # n_T     -- treatment-arm sample size
  # n_C     -- control-arm sample size
  # theta   -- true treatment effect
  # sigma2  -- known outcome variance
  # xbar_ch -- historical control mean
  # se_ch   -- standard error of the historical control mean;
  #            se_ch = sigma / sqrt(n_ch)
  # B_thr   -- calibrated posterior probability threshold
  # psi     -- scale parameter of the half-normal prior for tau
  # z_C     -- simulated standard normal values used to generate
  #            the control-arm sample mean xbar_C
  # z_T     -- simulated standard normal values used to generate
  #            the treatment-arm sample mean xbar_T
  # theta0  -- superiority margin
  # n_tau   -- number of grid points used for numerical integration over tau
  
  ## Generate the current control- and treatment-arm sample means
  xbar_C <- muC0 + sqrt(sigma2 / n_C) * z_C
  xbar_T <- muC0 + theta + sqrt(sigma2 / n_T) * z_T
  
  ## Construct the tau grid and prior integration weights
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w0 <- tg$w
  
  v_ch <- se_ch^2
  
  ## Prior variance of the current control mean conditional on tau
  v0_tau <- v_ch + 2 * tau^2
  vT <- sigma2 / n_T
  
  ## Predictive distribution of the current control mean conditional on tau
  pred_var <- outer(
    rep(1, length(xbar_C)),
    sigma2 / n_C + v0_tau
  )
  
  pred_mean <- matrix(
    xbar_ch,
    nrow = length(xbar_C),
    ncol = length(tau)
  )
  
  ## Compute posterior weights for tau
  log_w <- sweep(
    dnorm(
      xbar_C,
      mean = pred_mean,
      sd = sqrt(pred_var),
      log = TRUE
    ),
    2,
    log(w0),
    "+"
  )
  
  ## Stabilise and normalise the posterior weights
  log_w <- log_w - apply(log_w, 1, max)
  w_post <- exp(log_w)
  w_post <- w_post / rowSums(w_post)
  
  ## Posterior mean and variance of the control-arm mean conditional on tau
  prec_C <- n_C / sigma2
  prec_0 <- 1 / v0_tau
  
  v_C_tau <- 1 / outer(
    rep(1, length(xbar_C)),
    prec_C + prec_0
  )
  
  m_C_tau <- v_C_tau * (
    xbar_C * prec_C +
      outer(
        rep(1, length(xbar_C)),
        xbar_ch * prec_0
      )
  )
  
  ## Posterior mean and variance of the treatment effect conditional on tau
  m_theta_tau <- xbar_T - m_C_tau
  v_theta_tau <- vT + v_C_tau
  
  ## Integrate Pr(theta < theta0 | data, tau) over the posterior distribution of tau
  pneg <- rowSums(
    w_post *
      pnorm(
        (theta0 - m_theta_tau) /
          sqrt(v_theta_tau)
      )
  )
  
  ## Estimated probability of rejecting H_0
  mean(pneg >= 1 - B_thr)
}

#---- Compute the calibrated posterior probability threshold under the MAP prior over Delta in [DEL_l, DEL_u] ----
find_global_min_p_quantile_map <- function(DEL_l, DEL_u,
                                           z_C, z_T, xbar_ch,
                                           n_T, n_C, sigma2, se_ch,
                                           psi, theta0,
                                           alpha,
                                           grid_len = 11,
                                           n_tau = 1001,
                                           ncore = 8) {
  
  # DEL_l    -- lower bound of the discrepancy range
  # DEL_u    -- upper bound of the discrepancy range
  # z_C      -- simulated standard normal values used to generate
  #             the control-arm sample mean xbar_C
  # z_T      -- simulated standard normal values used to generate
  #             the treatment-arm sample mean xbar_T
  # xbar_ch  -- historical control mean
  # n_T      -- treatment-arm sample size
  # n_C      -- control-arm sample size
  # sigma2   -- known outcome variance
  # se_ch    -- standard error of the historical control mean;
  #             se_ch = sigma / sqrt(n_ch)
  # psi      -- scale parameter of the half-normal prior for tau
  # theta0   -- superiority margin
  # alpha    -- nominal Type I error rate
  # grid_len -- number of grid points used in the initial search
  # n_tau    -- number of grid points used for numerical integration over tau
  # ncore    -- number of CPU cores used for parallel computation
  
  ## Define the range of true control-arm means
  lower <- xbar_ch + DEL_l
  upper <- xbar_ch + DEL_u
  
  ## Construct the initial search grid
  mu_grid <- seq(lower, upper, length.out = grid_len)
  
  ## Evaluate the posterior probability threshold over the grid
  if (ncore > 1) {
    
    cl <- parallel::makeCluster(ncore)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterExport(
      cl,
      c(
        "make_tau_grid",
        "p_quantile_fast_map",
        "ppos_map_vec",
        "z_C", "z_T",
        "xbar_ch",
        "n_T", "n_C",
        "sigma2", "se_ch",
        "psi",
        "theta0",
        "alpha",
        "n_tau"
      ),
      envir = environment()
    )
    
    q_grid <- unlist(
      parallel::parLapply(cl, mu_grid, function(mu) {
        p_quantile_fast_map(
          muC = mu,
          z_C = z_C,
          z_T = z_T,
          xbar_ch = xbar_ch,
          n_T = n_T,
          n_C = n_C,
          sigma2 = sigma2,
          se_ch = se_ch,
          psi = psi,
          theta0 = theta0,
          alpha = alpha,
          n_tau = n_tau
        )
      })
    )
    
  } else {
    
    q_grid <- sapply(mu_grid, function(mu) {
      p_quantile_fast_map(
        muC = mu,
        z_C = z_C,
        z_T = z_T,
        xbar_ch = xbar_ch,
        n_T = n_T,
        n_C = n_C,
        sigma2 = sigma2,
        se_ch = se_ch,
        psi = psi,
        theta0 = theta0,
        alpha = alpha,
        n_tau = n_tau
      )
    })
  }
  
  ## Identify local minima on the initial grid
  idx_local <- which(
    q_grid[2:(length(q_grid) - 1)] <
      q_grid[1:(length(q_grid) - 2)] &
      q_grid[2:(length(q_grid) - 1)] <
      q_grid[3:length(q_grid)]
  ) + 1
  
  ## Include the two boundary points as candidate minima
  idx_candidate <- unique(
    c(1, length(mu_grid), idx_local)
  )
  
  ## Refine each candidate minimum using one-dimensional optimisation
  local_res <- lapply(idx_candidate, function(j) {
    
    if (j == 1 || j == length(mu_grid)) {
      return(
        list(
          minimum = mu_grid[j],
          objective = q_grid[j]
        )
      )
    }
    
    optimise(
      p_quantile_fast_map,
      interval = c(mu_grid[j - 1], mu_grid[j + 1]),
      z_C = z_C,
      z_T = z_T,
      xbar_ch = xbar_ch,
      n_T = n_T,
      n_C = n_C,
      sigma2 = sigma2,
      se_ch = se_ch,
      psi = psi,
      theta0 = theta0,
      alpha = alpha,
      n_tau = n_tau,
      maximum = FALSE,
      tol = 1e-3
    )
  })
  
  ## Select the global minimum among all candidate minima
  local_vals <- sapply(local_res, function(x) x$objective)
  local_mu <- sapply(local_res, function(x) x$minimum)
  
  k <- which.min(local_vals)
  
  list(
    p_0 = local_vals[k],
    muC_at_min = local_mu[k],
    mu_grid = mu_grid,
    q_grid = q_grid,
    all_candidate_mu = local_mu,
    all_candidate_val = local_vals
  )
}

#---- Compute the difference between the worst-case and target power over Delta in [DEL_l, DEL_u] ----
cal_f_map <- function(DEL_l, DEL_u,
                      n_T, n_C,
                      alpha, beta, theta0, theta1,
                      sigma2,
                      xbar_ch, se_ch, psi,
                      z_C, z_T,
                      n_tau = 1001,
                      ncore0 = 8) {
  
  # DEL_l   -- lower bound of the discrepancy range
  # DEL_u   -- upper bound of the discrepancy range
  # n_T     -- treatment-arm sample size
  # n_C     -- control-arm sample size
  # alpha   -- nominal Type I error rate
  # beta    -- target Type II error rate, i.e., 1 - power
  # theta0  -- superiority margin
  # theta1  -- treatment effect under the alternative hypothesis
  # sigma2  -- known outcome variance
  # xbar_ch -- historical control mean
  # se_ch   -- standard error of the historical control mean;
  #            se_ch = sigma / sqrt(n_ch)
  # psi     -- scale parameter of the half-normal prior for tau
  # z_C     -- simulated standard normal values used to generate
  #            the control-arm sample mean xbar_C
  # z_T     -- simulated standard normal values used to generate
  #            the treatment-arm sample mean xbar_T
  # n_tau   -- number of grid points used for numerical integration over tau
  # ncore0  -- number of CPU cores used for parallel computation
  
  ## Compute the calibrated posterior probability threshold
  p_0 <- find_global_min_p_quantile_map(
    DEL_l = DEL_l,
    DEL_u = DEL_u,
    z_C = z_C,
    z_T = z_T,
    xbar_ch = xbar_ch,
    n_T = n_T,
    n_C = n_C,
    sigma2 = sigma2,
    se_ch = se_ch,
    psi = psi,
    theta0 = theta0,
    alpha = alpha,
    grid_len = 11,
    n_tau = n_tau,
    ncore = ncore0
  )$p_0
  
  ## Construct a grid of true control-arm means over the discrepancy range
  muC0_v <- seq(
    from = DEL_l + xbar_ch,
    to = DEL_u + xbar_ch,
    length.out = 10
  )
  
  ## Define power as a function of the true control-arm mean
  obj <- function(muC0) {
    estimate_OC_map_p_fast(
      muC0 = muC0,
      n_T = n_T,
      n_C = n_C,
      theta = theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      B_thr = p_0,
      psi = psi,
      z_C = z_C,
      z_T = z_T,
      theta0 = theta0,
      n_tau = n_tau
    )
  }
  
  ## Evaluate power over the initial grid
  min_pow_v <- sapply(muC0_v, obj)
  
  ## Identify local minima on the initial grid
  idx_local <- which(
    min_pow_v[2:(length(min_pow_v) - 1)] <
      min_pow_v[1:(length(min_pow_v) - 2)] &
      min_pow_v[2:(length(min_pow_v) - 1)] <
      min_pow_v[3:length(min_pow_v)]
  ) + 1
  
  ## Use the minimum grid value as an initial candidate
  candidate_vals <- min(min_pow_v)
  candidate_mu <- muC0_v[which.min(min_pow_v)]
  
  ## Refine each local minimum using one-dimensional optimisation
  if (length(idx_local) > 0) {
    
    local_res <- lapply(idx_local, function(j) {
      optimise(
        obj,
        interval = c(muC0_v[j - 1], muC0_v[j + 1]),
        maximum = FALSE,
        tol = 1e-3
      )
    })
    
    local_vals <- sapply(local_res, function(x) x$objective)
    local_mu <- sapply(local_res, function(x) x$minimum)
    
    all_vals <- c(candidate_vals, local_vals)
    all_mu <- c(candidate_mu, local_mu)
    
    k <- which.min(all_vals)
    
    min_pow_o <- all_vals[k]
    muC0_opt <- all_mu[k]
    
  } else {
    
    min_pow_o <- candidate_vals
    muC0_opt <- candidate_mu
  }
  
  ## Return the power difference and calibrated quantities
  list(
    res = min_pow_o - (1 - beta),
    p_0 = p_0,
    min_pow = min_pow_o,
    muC0_opt = muC0_opt
  )
}

#---- Compute the calibrated borrowing parameter psi for the MAP prior ----
find_all_psi <- function(DEL_l, DEL_u,
                         n_T, n_C,
                         alpha, beta,
                         theta0, theta1,
                         sigma2,
                         xbar_ch, se_ch,
                         z_C, z_T,
                         psi_grid = seq(
                           0.001,
                           5 * sqrt(sigma2),
                           length.out = 30
                         ),
                         n_tau = 401,
                         ncore0 = 8,
                         tol = 1e-6) {
  
  # DEL_l    -- lower bound of the discrepancy range
  # DEL_u    -- upper bound of the discrepancy range
  # n_T      -- treatment-arm sample size
  # n_C      -- control-arm sample size
  # alpha    -- nominal Type I error rate
  # beta     -- target Type II error rate, i.e., 1 - power
  # theta0   -- superiority margin
  # theta1   -- treatment effect under the alternative hypothesis
  # sigma2   -- known outcome variance
  # xbar_ch  -- historical control mean
  # se_ch    -- standard error of the historical control mean;
  #             se_ch = sigma / sqrt(n_ch)
  # z_C      -- simulated standard normal values used to generate
  #             the control-arm sample mean xbar_C
  # z_T      -- simulated standard normal values used to generate
  #             the treatment-arm sample mean xbar_T
  # psi_grid -- grid of candidate values for the MAP borrowing parameter psi
  # n_tau    -- number of grid points used for numerical integration over tau
  # ncore0   -- number of CPU cores used for parallel computation
  # tol      -- numerical tolerance for identifying a root
  
  ## Prepare the search grid for psi
  psi_grid <- sort(unique(psi_grid))
  n_psi <- length(psi_grid)
  
  if (n_psi < 2L) {
    stop("psi_grid must contain at least two distinct values.")
  }
  
  ## Evaluate the calibration function for a given psi
  evaluate_psi <- function(psi) {
    cal_f_map(
      DEL_l = DEL_l,
      DEL_u = DEL_u,
      n_T = n_T,
      n_C = n_C,
      alpha = alpha,
      beta = beta,
      theta0 = theta0,
      theta1 = theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      psi = psi,
      z_C = z_C,
      z_T = z_T,
      n_tau = n_tau,
      ncore0 = ncore0
    )
  }
  
  ## Objective function: worst-case power minus target power
  f <- function(psi) {
    evaluate_psi(psi)$res
  }
  
  ## Store only the grid points evaluated before stopping
  psi_checked <- numeric(0)
  res_checked <- numeric(0)
  
  ## Evaluate the smallest value of psi
  psi_prev <- psi_grid[1L]
  out <- evaluate_psi(psi_prev)
  res_prev <- out$res
  
  if (!is.finite(res_prev)) {
    stop(
      "A non-finite residual is obtained at psi = ",
      psi_prev,
      "."
    )
  }
  
  psi_checked <- c(psi_checked, psi_prev)
  res_checked <- c(res_checked, res_prev)
  
  # ---------------------------------------------------------------
  # Case 1:
  # The residual at the smallest psi is already non-negative
  # ---------------------------------------------------------------
  
  if (res_prev >= -tol) {
    
    psi_use <- psi_prev
    
    root_case <- if (abs(res_prev) <= tol) {
      "minimum_grid_value_is_zero"
    } else {
      "minimum_grid_value_is_nonnegative"
    }
    
    return(list(
      has_root = abs(res_prev) <= tol,
      root_case = root_case,
      psi_root = if (abs(res_prev) <= tol) {
        psi_use
      } else {
        numeric(0)
      },
      root_type = if (abs(res_prev) <= tol) {
        "grid_zero"
      } else {
        "boundary_nonnegative"
      },
      psi_use = psi_use,
      res_use = out$res,
      p_0 = out$p_0,
      min_pow = out$min_pow,
      muC0_opt = out$muC0_opt,
      psi_grid = psi_checked,
      res_grid = res_checked,
      message = paste(
        "The residual at the smallest psi-grid value is non-negative.",
        "Therefore, the smallest value in psi_grid is selected."
      )
    ))
  }
  
  # ---------------------------------------------------------------
  # The residual starts negative.
  # Sequentially evaluate the remaining grid points until the first
  # numerical zero or negative-to-positive crossing is identified.
  # ---------------------------------------------------------------
  
  for (i in 2L:n_psi) {
    
    psi_curr <- psi_grid[i]
    out <- evaluate_psi(psi_curr)
    res_curr <- out$res
    
    if (!is.finite(res_curr)) {
      stop(
        "A non-finite residual is obtained at psi = ",
        psi_curr,
        "."
      )
    }
    
    psi_checked <- c(psi_checked, psi_curr)
    res_checked <- c(res_checked, res_curr)
    
    # -------------------------------------------------------------
    # Case 2:
    # The current grid value is numerically equal to zero
    # -------------------------------------------------------------
    
    if (abs(res_curr) <= tol) {
      
      psi_use <- psi_curr
      
      return(list(
        has_root = TRUE,
        root_case = "negative_to_numerical_zero",
        psi_root = psi_use,
        root_type = "grid_zero",
        psi_use = psi_use,
        res_use = out$res,
        p_0 = out$p_0,
        min_pow = out$min_pow,
        muC0_opt = out$muC0_opt,
        psi_grid = psi_checked,
        res_grid = res_checked,
        transition_index = i,
        message = paste(
          "The residual starts negative.",
          "The first psi-grid value at which the residual is",
          "numerically equal to zero is selected."
        )
      ))
    }
    
    # -------------------------------------------------------------
    # Case 3:
    # Strict negative-to-positive crossing
    # -------------------------------------------------------------
    
    if (res_prev < -tol && res_curr > tol) {
      
      psi_use <- uniroot(
        f,
        interval = c(psi_prev, psi_curr),
        tol = tol
      )$root
      
      out <- evaluate_psi(psi_use)
      
      return(list(
        has_root = TRUE,
        root_case = "negative_to_positive_crossing",
        psi_root = psi_use,
        root_type = "uniroot",
        psi_use = psi_use,
        res_use = out$res,
        p_0 = out$p_0,
        min_pow = out$min_pow,
        muC0_opt = out$muC0_opt,
        psi_grid = psi_checked,
        res_grid = res_checked,
        transition_interval = c(psi_prev, psi_curr),
        transition_index = i - 1L,
        message = paste(
          "The residual starts negative.",
          "The first interval over which the residual changes from",
          "negative to positive is identified, and its root is selected."
        )
      ))
    }
    
    ## Move to the next interval
    psi_prev <- psi_curr
    res_prev <- res_curr
  }
  
  # ---------------------------------------------------------------
  # Case 4:
  # No negative-to-non-negative transition over the specified range
  # ---------------------------------------------------------------
  
  list(
    has_root = FALSE,
    root_case = "no_negative_to_nonnegative_transition",
    psi_root = numeric(0),
    root_type = NA_character_,
    psi_use = NA_real_,
    res_use = max(res_checked),
    p_0 = NA_real_,
    min_pow = NA_real_,
    muC0_opt = NA_real_,
    psi_grid = psi_checked,
    res_grid = res_checked,
    message = paste(
      "The residual is negative at the smallest psi-grid value,",
      "and no transition from negative to non-negative is found",
      "over the specified psi_grid."
    )
  )
}

#---- Cache evaluations of the MAP calibration function ----
cal_cache_map <- new.env(parent = emptyenv())

#---- Evaluate the MAP calibration function with caching ----
cal_f_cached_map <- function(dl, du) {
  
  # dl -- lower bound of the discrepancy range
  # du -- upper bound of the discrepancy range
  
  ## Create a unique cache key for the discrepancy range and psi
  key <- paste(
    round(dl, 6),
    round(du, 6),
    round(psi, 6),
    sep = "_"
  )
  
  ## Return the cached result if it has already been evaluated
  if (exists(key, envir = cal_cache_map)) {
    return(get(key, envir = cal_cache_map))
  }
  
  ## Evaluate the MAP calibration function
  val <- cal_f_map(
    DEL_l = dl,
    DEL_u = du,
    n_T = n_T,
    n_C = n_C,
    alpha = alpha,
    beta = beta,
    theta0 = theta0,
    theta1 = theta1,
    sigma2 = sigma2,
    xbar_ch = xbar_ch,
    se_ch = se_ch,
    psi = psi,
    z_C = z_C,
    z_T = z_T,
    n_tau = n_tau,
    ncore0 = 1
  )
  
  ## Store the result in the cache
  assign(key, val, envir = cal_cache_map)
  
  val
}

#---- Compute the EHSS given xbar_C under the MAP prior ----
ESS_MAP_xbarC <- function(xbar_C_grid,
                          xbar_ch,
                          n_C,
                          sigma2,
                          se_ch,
                          psi,
                          weights = NULL,
                          n_tau = 1001) {
  
  # xbar_C_grid -- grid of current control-arm sample means
  # xbar_ch     -- historical control mean
  # n_C         -- control-arm sample size
  # sigma2      -- known outcome variance
  # se_ch       -- standard error of the historical control mean;
  #                se_ch = sigma / sqrt(n_ch)
  # psi         -- scale parameter of the half-normal prior for tau
  # weights     -- optional weights used to compute the mean EHSS
  # n_tau       -- number of grid points used for numerical integration over tau
  
  ## Construct the tau grid and prior integration weights
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau_prior <- tg$w
  
  ## Prior variance of the current control mean conditional on tau
  v0_tau <- se_ch^2 + 2 * tau^2
  
  ## Predictive variance of xbar_C conditional on tau
  pred_var_tau <- sigma2 / n_C + v0_tau
  
  ## Evaluate the predictive density of xbar_C conditional on tau
  dens_mat <- sapply(seq_along(tau), function(j) {
    dnorm(
      xbar_C_grid,
      mean = xbar_ch,
      sd = sqrt(pred_var_tau[j])
    )
  })
  
  ## Compute posterior weights for tau given xbar_C
  w_tau_given_xbarC <- dens_mat * matrix(
    w_tau_prior,
    nrow = length(xbar_C_grid),
    ncol = length(tau),
    byrow = TRUE
  )
  
  w_tau_given_xbarC <-
    w_tau_given_xbarC / rowSums(w_tau_given_xbarC)
  
  ## Posterior mean and variance of the control-arm mean conditional on tau
  prec_C <- n_C / sigma2
  prec_0 <- 1 / v0_tau
  
  v_C_tau <- 1 / (prec_C + prec_0)
  
  m_C_tau <- sapply(seq_along(tau), function(j) {
    v_C_tau[j] * (
      xbar_C_grid * prec_C +
        xbar_ch * prec_0[j]
    )
  })
  
  ## Compute the posterior mixture mean of the control-arm mean
  m_mix <- rowSums(
    w_tau_given_xbarC * m_C_tau
  )
  
  ## Compute the posterior mixture variance of the control-arm mean
  V_mix <- rowSums(
    w_tau_given_xbarC * (
      m_C_tau^2 +
        matrix(
          v_C_tau,
          nrow = length(xbar_C_grid),
          ncol = length(tau),
          byrow = TRUE
        )
    )
  ) - m_mix^2
  
  ## Compute the EHSS due to historical borrowing
  ESS_grid <- sigma2 / V_mix - n_C
  
  ## Specify weights for computing the mean EHSS
  if (is.null(weights)) {
    weights <- rep(
      1 / length(xbar_C_grid),
      length(xbar_C_grid)
    )
  } else {
    weights <- weights / sum(weights)
  }
  
  ## Return the EHSS and posterior summaries
  list(
    mean_ESS = sum(weights * ESS_grid),
    max_ESS = max(ESS_grid),
    xbarC_at_max_ESS = xbar_C_grid[which.max(ESS_grid)],
    ESS_grid = ESS_grid,
    V_post = V_mix,
    mean_post = m_mix,
    tau = tau,
    w_tau_prior = w_tau_prior,
    w_tau_given_xbarC = w_tau_given_xbarC,
    v0_tau = v0_tau,
    v_C_tau = v_C_tau
  )
}

#---- Compute the expected EHSS given theta_C under the MAP prior ----
ESS_MAP_xbarC_given_mu <- function(mu_grid,
                                   xbar_ch,
                                   n_C,
                                   sigma2,
                                   se_ch,
                                   psi,
                                   weights = NULL,
                                   n_tau = 41,
                                   z_grid = seq(-5, 5, length.out = 1001)) {
  
  # mu_grid -- grid of true control-arm means
  # xbar_ch -- historical control mean
  # n_C     -- control-arm sample size
  # sigma2  -- known outcome variance
  # se_ch   -- standard error of the historical control mean;
  #            se_ch = sigma / sqrt(n_ch)
  # psi     -- scale parameter of the half-normal prior for tau
  # weights -- optional weights used to compute the mean expected EHSS
  # n_tau   -- number of grid points used for numerical integration over tau
  # z_grid  -- standard normal grid used for numerical integration over xbar_C
  
  ## Construct the tau grid and prior integration weights
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau_prior <- tg$w
  
  ## Prior variance of the current control mean conditional on tau
  v0_tau <- se_ch^2 + 2 * tau^2
  
  ## Compute the EHSS for a given observed control-arm sample mean
  ESS_one_xbarC <- function(xbar_C) {
    
    ## Predictive variance of xbar_C conditional on tau
    pred_var_tau <- sigma2 / n_C + v0_tau
    
    ## Compute posterior weights for tau given xbar_C
    log_w <- dnorm(
      xbar_C,
      mean = xbar_ch,
      sd = sqrt(pred_var_tau),
      log = TRUE
    ) + log(w_tau_prior)
    
    ## Stabilise and normalise the posterior weights
    log_w <- log_w - max(log_w)
    w_post <- exp(log_w)
    w_post <- w_post / sum(w_post)
    
    ## Posterior mean and variance of the control-arm mean conditional on tau
    prec_C <- n_C / sigma2
    prec_0 <- 1 / v0_tau
    
    v_C_tau <- 1 / (prec_C + prec_0)
    
    m_C_tau <- v_C_tau * (
      xbar_C * prec_C +
        xbar_ch * prec_0
    )
    
    ## Compute the posterior mixture mean and variance
    m_mix <- sum(w_post * m_C_tau)
    
    V_mix <- sum(
      w_post * (v_C_tau + m_C_tau^2)
    ) - m_mix^2
    
    ## Compute the EHSS due to historical borrowing
    sigma2 / V_mix - n_C
  }
  
  ## Construct weights for numerical integration over xbar_C | mu_C
  dz <- z_grid[2] - z_grid[1]
  z_w <- dnorm(z_grid) * dz
  z_w <- z_w / sum(z_w)
  
  ## Compute E_{xbar_C | mu_C}[EHSS(xbar_C)] for each true control-arm mean
  ESS_grid <- sapply(mu_grid, function(mu_C) {
    
    xbar_C_grid <- mu_C +
      sqrt(sigma2 / n_C) * z_grid
    
    ESS_xbarC <- sapply(
      xbar_C_grid,
      ESS_one_xbarC
    )
    
    sum(z_w * ESS_xbarC)
  })
  
  ## Specify weights for computing the mean expected EHSS
  if (is.null(weights)) {
    weights <- rep(
      1 / length(mu_grid),
      length(mu_grid)
    )
  } else {
    weights <- weights / sum(weights)
  }
  
  ## Return the expected EHSS summaries
  list(
    mean_ESS = sum(weights * ESS_grid),
    max_ESS = max(ESS_grid),
    muC_at_max_ESS = mu_grid[which.max(ESS_grid)],
    ESS_grid = ESS_grid,
    tau = tau,
    w_tau_prior = w_tau_prior,
    v0_tau = v0_tau,
    z_grid = z_grid,
    z_weights = z_w
  )
}

#---- Compute the MAP prior density of theta_C ----
dmap_prior <- function(theta_C, xbar_ch, se_ch, psi, n_tau = 1001) {
  
  # theta_C -- value or vector of values of the true control-arm mean
  # xbar_ch -- historical control mean
  # se_ch   -- standard error of the historical control mean;
  #            se_ch = sigma / sqrt(n_ch)
  # psi     -- scale parameter of the half-normal prior for tau
  # n_tau   -- number of grid points used for numerical integration over tau
  
  ## Construct the tau grid and prior integration weights
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w <- tg$w
  
  ## Prior variance of theta_C conditional on tau
  v_ch <- se_ch^2
  v0_tau <- v_ch + 2 * tau^2
  
  ## Evaluate the weighted conditional prior density for each tau
  dens_mat <- sapply(seq_along(tau), function(j) {
    w[j] * dnorm(
      theta_C,
      mean = xbar_ch,
      sd = sqrt(v0_tau[j])
    )
  })
  
  ## Integrate over the prior distribution of tau
  rowSums(dens_mat)
}


#### rMAP prior ####

#---- Construct the tau grid and integration weights under a half-normal prior ----
make_tau_grid <- function(psi, n_tau = 1001, prob = 0.999) {
  
  # psi   -- scale parameter of the half-normal prior for tau
  # n_tau -- number of grid points used for numerical integration
  # prob  -- cumulative probability used to determine the upper bound of tau
  
  ## Determine the upper bound of the tau grid
  tau_max <- psi * qnorm((prob + 1) / 2)
  
  ## Construct the tau grid
  tau <- seq(0, tau_max, length.out = n_tau)
  
  ## Compute integration interval widths
  d_tau <- rep(NA_real_, n_tau)
  d_tau[1] <- tau[2] - tau[1]
  d_tau[n_tau] <- tau[n_tau] - tau[n_tau - 1]
  d_tau[2:(n_tau - 1)] <-
    (tau[3:n_tau] - tau[1:(n_tau - 2)]) / 2
  
  ## Evaluate the half-normal density over the tau grid
  dens <- sqrt(2 / pi) / psi *
    exp(-tau^2 / (2 * psi^2))
  
  ## Compute and normalise the numerical integration weights
  w <- dens * d_tau
  w <- w / sum(w)
  
  ## Return the tau grid and integration weights
  list(
    tau = tau,
    w = w
  )
}


#---- Compute Pr(theta < theta0 | data) under the rMAP prior ----
ppos_rmap_vec <- function(xbar_T, xbar_C, xbar_ch,
                          n_T, n_C, sigma2, se_ch, psi, theta0,
                          w_rob = 0.8,
                          vague_sd = 88,
                          n_tau = 1001) {
  
  # xbar_T  -- treatment-arm sample mean
  # xbar_C  -- control-arm sample mean
  # xbar_ch -- historical control mean
  # n_T     -- treatment-arm sample size
  # n_C     -- control-arm sample size
  # sigma2  -- known outcome variance
  # se_ch   -- standard error of the historical control mean;
  #            se_ch = sigma / sqrt(n_ch)
  # psi     -- scale parameter of the half-normal prior for tau
  # theta0  -- superiority margin
  # w_rob   -- mixture weight assigned to the MAP component
  # vague_sd -- standard deviation of the vague component
  # n_tau   -- number of grid points used for numerical integration over tau
  
  ## Construct the tau grid and prior integration weights
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau0 <- tg$w
  
  ## Prior variances for the MAP and vague components
  v_map <- se_ch^2 + 2 * tau^2
  v_vag <- vague_sd^2
  
  ## Combine the MAP and vague components of the rMAP prior
  v0 <- c(v_map, v_vag)
  w0 <- c(w_rob * w_tau0, 1 - w_rob)
  
  ## Posterior variance contribution from the treatment arm
  vT <- sigma2 / n_T
  
  ## Predictive distribution of the current control mean
  pred_var <- outer(
    rep(1, length(xbar_C)),
    sigma2 / n_C + v0
  )
  
  pred_mean <- matrix(
    xbar_ch,
    nrow = length(xbar_C),
    ncol = length(v0)
  )
  
  ## Compute posterior mixture weights
  log_w <- sweep(
    dnorm(
      xbar_C,
      mean = pred_mean,
      sd = sqrt(pred_var),
      log = TRUE
    ),
    2,
    log(w0),
    "+"
  )
  
  ## Stabilise and normalise the posterior weights
  log_w <- log_w - apply(log_w, 1, max)
  w_post <- exp(log_w)
  w_post <- w_post / rowSums(w_post)
  
  ## Posterior mean and variance of the control-arm mean
  ## conditional on each mixture component
  prec_C <- n_C / sigma2
  prec_0 <- 1 / v0
  
  v_C <- 1 / outer(
    rep(1, length(xbar_C)),
    prec_C + prec_0
  )
  
  m_C <- v_C * (
    xbar_C * prec_C +
      outer(
        rep(1, length(xbar_C)),
        xbar_ch * prec_0
      )
  )
  
  ## Posterior mean and variance of the treatment effect
  ## conditional on each mixture component
  m_theta <- xbar_T - m_C
  v_theta <- vT + v_C
  
  ## Integrate Pr(theta < theta0 | data) over the posterior mixture weights
  rowSums(
    w_post *
      pnorm((theta0 - m_theta) / sqrt(v_theta))
  )
}


#---- Compute the calibrated posterior probability threshold under the rMAP prior ----
p_quantile_fast_rmap <- function(muC, z_C, z_T, xbar_ch,
                                 n_T, n_C, sigma2, se_ch, psi,
                                 theta0, alpha,
                                 w_rob = 0.8,
                                 vague_sd = 88,
                                 n_tau = 1001) {
  
  # muC      -- true control-arm mean
  # z_C      -- simulated standard normal values used to generate
  #             the control-arm sample mean xbar_C
  # z_T      -- simulated standard normal values used to generate
  #             the treatment-arm sample mean xbar_T
  # xbar_ch  -- historical control mean
  # n_T      -- treatment-arm sample size
  # n_C      -- control-arm sample size
  # sigma2   -- known outcome variance
  # se_ch    -- standard error of the historical control mean;
  #             se_ch = sigma / sqrt(n_ch)
  # psi      -- scale parameter of the half-normal prior for tau
  # theta0   -- superiority margin
  # alpha    -- nominal Type I error rate
  # w_rob    -- mixture weight assigned to the MAP component
  # vague_sd -- standard deviation of the vague component
  # n_tau    -- number of grid points used for numerical integration over tau
  
  ## Generate the current control- and treatment-arm sample means under H_0
  xbar_C <- muC + sqrt(sigma2 / n_C) * z_C
  xbar_T <- muC + theta0 + sqrt(sigma2 / n_T) * z_T
  
  ## Compute the posterior probability Pr(theta < theta0 | data)
  ppos <- ppos_rmap_vec(
    xbar_T = xbar_T,
    xbar_C = xbar_C,
    xbar_ch = xbar_ch,
    n_T = n_T,
    n_C = n_C,
    sigma2 = sigma2,
    se_ch = se_ch,
    psi = psi,
    theta0 = theta0,
    w_rob = w_rob,
    vague_sd = vague_sd,
    n_tau = n_tau
  )
  
  ## Calibrate the posterior probability threshold using the alpha quantile
  as.numeric(
    quantile(ppos, alpha, names = FALSE)
  )
}

#---- Compute the probability of rejecting H_0: theta<=theta_0 under the rMAP prior ----
estimate_OC_rmap_p_fast <- function(muC0, n_T, n_C, theta, sigma2,
                                    xbar_ch, se_ch, B_thr, psi,
                                    z_C, z_T,
                                    theta0 = 0,
                                    w_rob = 0.8,
                                    vague_sd = 88,
                                    n_tau = 1001) {
  
  # muC0     -- true control-arm mean
  # n_T      -- treatment-arm sample size
  # n_C      -- control-arm sample size
  # theta    -- true treatment effect
  # sigma2   -- known outcome variance
  # xbar_ch  -- historical control mean
  # se_ch    -- standard error of the historical control mean;
  #             se_ch = sigma / sqrt(n_ch)
  # B_thr    -- calibrated posterior probability threshold
  # psi      -- scale parameter of the half-normal prior for tau
  # z_C      -- simulated standard normal values used to generate
  #             the control-arm sample mean xbar_C
  # z_T      -- simulated standard normal values used to generate
  #             the treatment-arm sample mean xbar_T
  # theta0   -- superiority margin
  # w_rob    -- mixture weight assigned to the MAP component
  # vague_sd -- standard deviation of the vague component
  # n_tau    -- number of grid points used for numerical integration over tau
  
  ## Generate the current control- and treatment-arm sample means
  xbar_C <- muC0 + sqrt(sigma2 / n_C) * z_C
  xbar_T <- muC0 + theta + sqrt(sigma2 / n_T) * z_T
  
  ## Construct the tau grid and prior integration weights
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau0 <- tg$w
  
  ## Prior variances for the MAP and vague components
  v_map <- se_ch^2 + 2 * tau^2
  v_vag <- vague_sd^2
  
  ## Combine the MAP and vague components of the rMAP prior
  v0 <- c(v_map, v_vag)
  w0 <- c(w_rob * w_tau0, 1 - w_rob)
  
  ## Posterior variance contribution from the treatment arm
  vT <- sigma2 / n_T
  
  ## Predictive distribution of the current control mean
  ## conditional on each mixture component
  pred_var <- outer(
    rep(1, length(xbar_C)),
    sigma2 / n_C + v0
  )
  
  pred_mean <- matrix(
    xbar_ch,
    nrow = length(xbar_C),
    ncol = length(v0)
  )
  
  ## Compute posterior mixture weights
  log_w <- sweep(
    dnorm(
      xbar_C,
      mean = pred_mean,
      sd = sqrt(pred_var),
      log = TRUE
    ),
    2,
    log(w0),
    "+"
  )
  
  ## Stabilise and normalise the posterior weights
  log_w <- log_w - apply(log_w, 1, max)
  w_post <- exp(log_w)
  w_post <- w_post / rowSums(w_post)
  
  ## Posterior mean and variance of the control-arm mean
  ## conditional on each mixture component
  prec_C <- n_C / sigma2
  prec_0 <- 1 / v0
  
  v_C <- 1 / outer(
    rep(1, length(xbar_C)),
    prec_C + prec_0
  )
  
  m_C <- v_C * (
    xbar_C * prec_C +
      outer(
        rep(1, length(xbar_C)),
        xbar_ch * prec_0
      )
  )
  
  ## Posterior mean and variance of the treatment effect
  ## conditional on each mixture component
  m_theta <- xbar_T - m_C
  v_theta <- vT + v_C
  
  ## Integrate the posterior probability over the posterior mixture weights
  ppos <- rowSums(
    w_post *
      pnorm((m_theta - theta0) / sqrt(v_theta))
  )
  
  ## Estimated probability of rejecting H_0
  mean(ppos >= 1 - B_thr)
}


#---- Compute the probability of rejecting H_0: theta>=theta_0 under the rMAP prior ----
estimate_OC_rmap_p_fast_1 <- function(muC0, n_T, n_C, theta, sigma2,
                                      xbar_ch, se_ch, B_thr, psi,
                                      z_C, z_T,
                                      theta0 = 0,
                                      w_rob = 0.8,
                                      vague_sd = 88,
                                      n_tau = 1001) {
  
  # muC0     -- true control-arm mean
  # n_T      -- treatment-arm sample size
  # n_C      -- control-arm sample size
  # theta    -- true treatment effect
  # sigma2   -- known outcome variance
  # xbar_ch  -- historical control mean
  # se_ch    -- standard error of the historical control mean;
  #             se_ch = sigma / sqrt(n_ch)
  # B_thr    -- calibrated posterior probability threshold
  # psi      -- scale parameter of the half-normal prior for tau
  # z_C      -- simulated standard normal values used to generate
  #             the control-arm sample mean xbar_C
  # z_T      -- simulated standard normal values used to generate
  #             the treatment-arm sample mean xbar_T
  # theta0   -- superiority margin
  # w_rob    -- mixture weight assigned to the MAP component
  # vague_sd -- standard deviation of the vague component
  # n_tau    -- number of grid points used for numerical integration over tau
  
  ## Generate the current control- and treatment-arm sample means
  xbar_C <- muC0 + sqrt(sigma2 / n_C) * z_C
  xbar_T <- muC0 + theta + sqrt(sigma2 / n_T) * z_T
  
  ## Construct the tau grid and prior integration weights
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau0 <- tg$w
  
  ## Prior variances for the MAP and vague components
  v_map <- se_ch^2 + 2 * tau^2
  v_vag <- vague_sd^2
  
  ## Combine the MAP and vague components of the rMAP prior
  v0 <- c(v_map, v_vag)
  w0 <- c(w_rob * w_tau0, 1 - w_rob)
  
  ## Posterior variance contribution from the treatment arm
  vT <- sigma2 / n_T
  
  ## Predictive distribution of the current control mean
  ## conditional on each mixture component
  pred_var <- outer(
    rep(1, length(xbar_C)),
    sigma2 / n_C + v0
  )
  
  pred_mean <- matrix(
    xbar_ch,
    nrow = length(xbar_C),
    ncol = length(v0)
  )
  
  ## Compute posterior mixture weights
  log_w <- sweep(
    dnorm(
      xbar_C,
      mean = pred_mean,
      sd = sqrt(pred_var),
      log = TRUE
    ),
    2,
    log(w0),
    "+"
  )
  
  ## Stabilise and normalise the posterior weights
  log_w <- log_w - apply(log_w, 1, max)
  w_post <- exp(log_w)
  w_post <- w_post / rowSums(w_post)
  
  ## Posterior mean and variance of the control-arm mean
  ## conditional on each mixture component
  prec_C <- n_C / sigma2
  prec_0 <- 1 / v0
  
  v_C <- 1 / outer(
    rep(1, length(xbar_C)),
    prec_C + prec_0
  )
  
  m_C <- v_C * (
    xbar_C * prec_C +
      outer(
        rep(1, length(xbar_C)),
        xbar_ch * prec_0
      )
  )
  
  ## Posterior mean and variance of the treatment effect
  ## conditional on each mixture component
  m_theta <- xbar_T - m_C
  v_theta <- vT + v_C
  
  ## Integrate Pr(theta < theta0 | data) over the posterior mixture weights
  pneg <- rowSums(
    w_post *
      pnorm(
        (theta0 - m_theta) /
          sqrt(v_theta)
      )
  )
  
  ## Estimated probability of rejecting H_0
  mean(pneg >= 1 - B_thr)
}

#---- Compute the calibrated posterior probability threshold under the rMAP prior ----
find_global_min_p_quantile_rmap <- function(DEL_l, DEL_u,
                                            z_C, z_T, xbar_ch,
                                            n_T, n_C, se_ch, sigma2,
                                            psi, theta0,
                                            alpha,
                                            w_rob = 0.8,
                                            vague_sd = 88,
                                            grid_len = 21,
                                            n_tau = 1001,
                                            ncore = 8) {
  
  # DEL_l    -- lower bound of the discrepancy range
  # DEL_u    -- upper bound of the discrepancy range
  # z_C      -- simulated standard normal values used to generate
  #             the control-arm sample mean xbar_C
  # z_T      -- simulated standard normal values used to generate
  #             the treatment-arm sample mean xbar_T
  # xbar_ch  -- historical control mean
  # n_T      -- treatment-arm sample size
  # n_C      -- control-arm sample size
  # se_ch    -- standard error of the historical control mean;
  #             se_ch = sigma / sqrt(n_ch)
  # sigma2   -- known outcome variance
  # psi      -- scale parameter of the half-normal prior for tau
  # theta0   -- superiority margin
  # alpha    -- nominal Type I error rate
  # w_rob    -- mixture weight assigned to the MAP component
  # vague_sd -- standard deviation of the vague component
  # grid_len -- number of grid points used in the initial search
  # n_tau    -- number of grid points used for numerical integration over tau
  # ncore    -- number of CPU cores used for parallel computation
  
  ## Define the range of true control-arm means
  lower <- xbar_ch + DEL_l
  upper <- xbar_ch + DEL_u
  
  ## Construct the initial search grid
  mu_grid <- seq(
    lower,
    upper,
    length.out = grid_len
  )
  
  ## Define the posterior probability threshold as a function of mu_C
  eval_q <- function(mu) {
    p_quantile_fast_rmap(
      muC = mu,
      z_C = z_C,
      z_T = z_T,
      xbar_ch = xbar_ch,
      n_T = n_T,
      n_C = n_C,
      se_ch = se_ch,
      sigma2 = sigma2,
      psi = psi,
      theta0 = theta0,
      alpha = alpha,
      w_rob = w_rob,
      vague_sd = vague_sd,
      n_tau = n_tau
    )
  }
  
  ## Evaluate the posterior probability threshold over the grid
  if (ncore > 1) {
    
    cl <- parallel::makeCluster(ncore)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterExport(
      cl,
      c(
        "make_tau_grid",
        "p_quantile_fast_rmap",
        "ppos_rmap_vec",
        "z_C", "z_T", "xbar_ch",
        "n_T", "n_C", "se_ch", "sigma2",
        "psi", "theta0", "alpha",
        "w_rob", "vague_sd", "n_tau"
      ),
      envir = environment()
    )
    
    q_grid <- unlist(
      parallel::parLapply(cl, mu_grid, eval_q)
    )
    
  } else {
    
    q_grid <- sapply(mu_grid, eval_q)
  }
  
  ## Identify local minima on the initial grid
  idx_local <- which(
    q_grid[2:(length(q_grid) - 1)] <
      q_grid[1:(length(q_grid) - 2)] &
      q_grid[2:(length(q_grid) - 1)] <
      q_grid[3:length(q_grid)]
  ) + 1
  
  ## Include the two boundary points as candidate minima
  idx_candidate <- unique(
    c(1, length(mu_grid), idx_local)
  )
  
  ## Refine each candidate minimum using one-dimensional optimisation
  local_res <- lapply(idx_candidate, function(j) {
    
    if (j == 1 || j == length(mu_grid)) {
      return(
        list(
          minimum = mu_grid[j],
          objective = q_grid[j]
        )
      )
    }
    
    optimise(
      eval_q,
      interval = c(mu_grid[j - 1], mu_grid[j + 1]),
      maximum = FALSE,
      tol = 1e-3
    )
  })
  
  ## Select the global minimum among all candidate minima
  local_vals <- sapply(local_res, function(x) x$objective)
  local_mu <- sapply(local_res, function(x) x$minimum)
  
  k <- which.min(local_vals)
  
  list(
    p_0 = local_vals[k],
    muC_at_min = local_mu[k],
    mu_grid = mu_grid,
    q_grid = q_grid,
    all_candidate_mu = local_mu,
    all_candidate_val = local_vals
  )
}


#---- Compute the difference between the worst-case and target power over Delta in [DEL_l, DEL_u] under the rMAP prior ----
cal_f_rmap <- function(DEL_l, DEL_u,
                       n_T, n_C,
                       alpha, beta, theta0, theta1,
                       sigma2,
                       xbar_ch, se_ch, psi,
                       z_C, z_T,
                       w_rob = 0.8,
                       vague_sd = 88,
                       n_tau = 1001,
                       ncore0 = 8) {
  
  # DEL_l    -- lower bound of the discrepancy range
  # DEL_u    -- upper bound of the discrepancy range
  # n_T      -- treatment-arm sample size
  # n_C      -- control-arm sample size
  # alpha    -- nominal Type I error rate
  # beta     -- target Type II error rate, i.e., 1 - power
  # theta0   -- superiority margin
  # theta1   -- treatment effect under the alternative hypothesis
  # sigma2   -- known outcome variance
  # xbar_ch  -- historical control mean
  # se_ch    -- standard error of the historical control mean;
  #             se_ch = sigma / sqrt(n_ch)
  # psi      -- scale parameter of the half-normal prior for tau
  # z_C      -- simulated standard normal values used to generate
  #             the control-arm sample mean xbar_C
  # z_T      -- simulated standard normal values used to generate
  #             the treatment-arm sample mean xbar_T
  # w_rob    -- mixture weight assigned to the MAP component
  # vague_sd -- standard deviation of the vague component
  # n_tau    -- number of grid points used for numerical integration over tau
  # ncore0   -- number of CPU cores used for parallel computation
  
  ## Compute the calibrated posterior probability threshold
  p_0 <- find_global_min_p_quantile_rmap(
    DEL_l = DEL_l,
    DEL_u = DEL_u,
    z_C = z_C,
    z_T = z_T,
    xbar_ch = xbar_ch,
    n_T = n_T,
    n_C = n_C,
    se_ch = se_ch,
    sigma2 = sigma2,
    psi = psi,
    theta0 = theta0,
    alpha = alpha,
    w_rob = w_rob,
    vague_sd = vague_sd,
    grid_len = 11,
    n_tau = n_tau,
    ncore = ncore0
  )$p_0
  
  ## Construct a grid of true control-arm means over the discrepancy range
  muC0_v <- seq(
    xbar_ch + DEL_l,
    xbar_ch + DEL_u,
    length.out = 10
  )
  
  ## Define power as a function of the true control-arm mean
  obj <- function(muC0) {
    estimate_OC_rmap_p_fast(
      muC0 = muC0,
      n_T = n_T,
      n_C = n_C,
      theta = theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      B_thr = p_0,
      psi = psi,
      z_C = z_C,
      z_T = z_T,
      w_rob = w_rob,
      vague_sd = vague_sd,
      n_tau = n_tau
    )
  }
  
  ## Evaluate power over the initial grid
  min_pow_v <- sapply(muC0_v, obj)
  
  ## Identify local minima on the initial grid
  idx_local <- which(
    min_pow_v[2:(length(min_pow_v) - 1)] <
      min_pow_v[1:(length(min_pow_v) - 2)] &
      min_pow_v[2:(length(min_pow_v) - 1)] <
      min_pow_v[3:length(min_pow_v)]
  ) + 1
  
  ## Use the minimum grid value as an initial candidate
  candidate_vals <- min(min_pow_v)
  candidate_mu <- muC0_v[which.min(min_pow_v)]
  
  ## Refine each local minimum using one-dimensional optimisation
  if (length(idx_local) > 0) {
    
    local_res <- lapply(idx_local, function(j) {
      optimise(
        obj,
        interval = c(muC0_v[j - 1], muC0_v[j + 1]),
        maximum = FALSE,
        tol = 1e-3
      )
    })
    
    local_vals <- sapply(local_res, function(x) x$objective)
    local_mu <- sapply(local_res, function(x) x$minimum)
    
    all_vals <- c(candidate_vals, local_vals)
    all_mu <- c(candidate_mu, local_mu)
    
    k <- which.min(all_vals)
    
    min_pow_o <- all_vals[k]
    muC0_opt <- all_mu[k]
    
  } else {
    
    min_pow_o <- candidate_vals
    muC0_opt <- candidate_mu
  }
  
  ## Return the power difference and calibrated quantities
  list(
    res = min_pow_o - (1 - beta),
    p_0 = p_0,
    min_pow = min_pow_o,
    muC0_opt = muC0_opt
  )
}


#---- Compute the calibrated mixture weight w for the rMAP prior ----
find_all_w_rmap <- function(DEL_l, DEL_u,
                            n_T, n_C,
                            alpha, beta,
                            theta0, theta1,
                            sigma2,
                            xbar_ch, se_ch,
                            psi,
                            z_C, z_T,
                            w_grid = seq(
                              0.001,
                              0.99,
                              length.out = 30
                            ),
                            vague_sd = 88,
                            n_tau = 1001,
                            ncore0 = 8,
                            tol = 1e-3) {
  
  # DEL_l    -- lower bound of the discrepancy range
  # DEL_u    -- upper bound of the discrepancy range
  # n_T      -- treatment-arm sample size
  # n_C      -- control-arm sample size
  # alpha    -- nominal Type I error rate
  # beta     -- target Type II error rate, i.e., 1 - power
  # theta0   -- superiority margin
  # theta1   -- treatment effect under the alternative hypothesis
  # sigma2   -- known outcome variance
  # xbar_ch  -- historical control mean
  # se_ch    -- standard error of the historical control mean;
  #             se_ch = sigma / sqrt(n_ch)
  # psi      -- scale parameter of the half-normal prior for tau
  # z_C      -- simulated standard normal values used to generate
  #             the control-arm sample mean xbar_C
  # z_T      -- simulated standard normal values used to generate
  #             the treatment-arm sample mean xbar_T
  # w_grid   -- grid of candidate mixture weights for the MAP component
  # vague_sd -- standard deviation of the vague component
  # n_tau    -- number of grid points used for numerical integration over tau
  # ncore0   -- number of CPU cores used for parallel computation
  # tol      -- numerical tolerance for identifying a root
  
  ## Search from the largest to the smallest value of w
  w_grid <- sort(unique(w_grid), decreasing = TRUE)
  n_w <- length(w_grid)
  
  if (n_w < 2L) {
    stop("w_grid must contain at least two distinct values.")
  }
  
  ## Evaluate the calibration function for a given w
  evaluate_w <- function(w_rob) {
    cal_f_rmap(
      DEL_l = DEL_l,
      DEL_u = DEL_u,
      n_T = n_T,
      n_C = n_C,
      alpha = alpha,
      beta = beta,
      theta0 = theta0,
      theta1 = theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      psi = psi,
      z_C = z_C,
      z_T = z_T,
      w_rob = w_rob,
      vague_sd = vague_sd,
      n_tau = n_tau,
      ncore0 = ncore0
    )
  }
  
  ## Objective function: worst-case power minus target power
  f <- function(w_rob) {
    evaluate_w(w_rob)$res
  }
  
  ## Store only the values that are actually evaluated
  w_checked <- numeric(0)
  res_checked <- numeric(0)
  
  # ===============================================================
  # Evaluate the largest value of w
  # ===============================================================
  
  w_prev <- w_grid[1L]
  previous_info <- evaluate_w(w_prev)
  res_prev <- previous_info$res
  
  if (!is.finite(res_prev)) {
    stop(
      "A non-finite residual is obtained at w = ",
      w_prev,
      "."
    )
  }
  
  w_checked <- c(w_checked, w_prev)
  res_checked <- c(res_checked, res_prev)
  
  # ===============================================================
  # Case 1:
  # The residual at the largest w is already non-negative
  # ===============================================================
  
  if (res_prev >= -tol) {
    
    w_use <- w_prev
    
    root_case <- if (abs(res_prev) <= tol) {
      "maximum_grid_value_is_zero"
    } else {
      "maximum_grid_value_is_nonnegative"
    }
    
    return(list(
      has_root = abs(res_prev) <= tol,
      root_case = root_case,
      w_root = if (abs(res_prev) <= tol) {
        w_use
      } else {
        numeric(0)
      },
      root_type = if (abs(res_prev) <= tol) {
        "grid_zero"
      } else {
        "boundary_nonnegative"
      },
      w_use = w_use,
      psi = psi,
      res_use = res_prev,
      p_0 = previous_info$p_0,
      min_pow = previous_info$min_pow,
      muC0_opt = previous_info$muC0_opt,
      crossing_interval = NULL,
      transition_index = 1L,
      w_checked = w_checked,
      res_checked = res_checked,
      w_grid = w_grid,
      message = paste(
        "The residual at the largest w-grid value is non-negative.",
        "Therefore, the largest value in w_grid is selected."
      )
    ))
  }
  
  ## At this point, the residual at the largest w is negative
  
  # ===============================================================
  # Search from large to small w for the first transition from
  # negative to non-negative
  # ===============================================================
  
  for (i in 2L:n_w) {
    
    w_curr <- w_grid[i]
    current_info <- evaluate_w(w_curr)
    res_curr <- current_info$res
    
    if (!is.finite(res_curr)) {
      stop(
        "A non-finite residual is obtained at w = ",
        w_curr,
        "."
      )
    }
    
    w_checked <- c(w_checked, w_curr)
    res_checked <- c(res_checked, res_curr)
    
    # =============================================================
    # Case 2:
    # The current grid value is numerically equal to zero
    # =============================================================
    
    if (abs(res_curr) <= tol) {
      
      w_use <- w_curr
      
      return(list(
        has_root = TRUE,
        root_case = "negative_to_numerical_zero",
        w_root = w_use,
        root_type = "grid_zero",
        w_use = w_use,
        psi = psi,
        res_use = res_curr,
        p_0 = current_info$p_0,
        min_pow = current_info$min_pow,
        muC0_opt = current_info$muC0_opt,
        crossing_interval = c(w_curr, w_prev),
        transition_index = i,
        w_checked = w_checked,
        res_checked = res_checked,
        w_grid = w_grid,
        message = paste0(
          "Searching from large to small w, the first grid value ",
          "at which the residual is numerically zero is selected: ",
          "w = ",
          round(w_use, 6),
          "."
        )
      ))
    }
    
    # =============================================================
    # Case 3:
    # The residual changes from negative to positive as w decreases
    # =============================================================
    
    if (res_prev < -tol && res_curr > tol) {
      
      ## uniroot() requires an increasing interval
      w_use <- uniroot(
        f = f,
        interval = c(w_curr, w_prev),
        tol = tol
      )$root
      
      ## Evaluate the full calibration output at the root
      root_info <- evaluate_w(w_use)
      
      return(list(
        has_root = TRUE,
        root_case = "negative_to_positive_crossing",
        w_root = w_use,
        root_type = "uniroot",
        w_use = w_use,
        psi = psi,
        res_use = root_info$res,
        p_0 = root_info$p_0,
        min_pow = root_info$min_pow,
        muC0_opt = root_info$muC0_opt,
        crossing_interval = c(w_curr, w_prev),
        transition_index = i - 1L,
        w_checked = w_checked,
        res_checked = res_checked,
        w_grid = w_grid,
        message = paste0(
          "Searching from large to small w, the first interval ",
          "over which the residual changes from negative to positive ",
          "is identified. The largest qualifying root is selected: ",
          "w = ",
          round(w_use, 6),
          "."
        )
      ))
    }
    
    ## Move to the next, smaller value of w
    w_prev <- w_curr
    res_prev <- res_curr
    previous_info <- current_info
  }
  
  # ===============================================================
  # Case 4:
  # No negative-to-non-negative transition is found
  # ===============================================================
  
  list(
    has_root = FALSE,
    root_case = "no_negative_to_nonnegative_transition",
    w_root = numeric(0),
    root_type = NA_character_,
    w_use = NA_real_,
    psi = psi,
    res_use = max(res_checked),
    p_0 = NA_real_,
    min_pow = NA_real_,
    muC0_opt = NA_real_,
    crossing_interval = NULL,
    transition_index = NA_integer_,
    w_checked = w_checked,
    res_checked = res_checked,
    w_grid = w_grid,
    message = paste(
      "The residual is negative at the largest w-grid value,",
      "and no transition from negative to non-negative is found",
      "while searching from large to small w."
    )
  )
}


#---- Compute the calibrated borrowing parameter psi for the rMAP prior ----
find_all_psi_rmap <- function(DEL_l, DEL_u,
                              n_T, n_C,
                              alpha, beta,
                              theta0, theta1,
                              sigma2,
                              xbar_ch, se_ch,
                              w_rob,
                              z_C, z_T,
                              psi_grid = seq(
                                0.001,
                                5 * sqrt(sigma2),
                                length.out = 30
                              ),
                              vague_sd = 88,
                              n_tau = 1001,
                              ncore0 = 8,
                              tol = 1e-3) {
  
  # DEL_l    -- lower bound of the discrepancy range
  # DEL_u    -- upper bound of the discrepancy range
  # n_T      -- treatment-arm sample size
  # n_C      -- control-arm sample size
  # alpha    -- nominal Type I error rate
  # beta     -- target Type II error rate, i.e., 1 - power
  # theta0   -- superiority margin
  # theta1   -- treatment effect under the alternative hypothesis
  # sigma2   -- known outcome variance
  # xbar_ch  -- historical control mean
  # se_ch    -- standard error of the historical control mean;
  #             se_ch = sigma / sqrt(n_ch)
  # w_rob    -- mixture weight assigned to the MAP component
  # z_C      -- simulated standard normal values used to generate
  #             the control-arm sample mean xbar_C
  # z_T      -- simulated standard normal values used to generate
  #             the treatment-arm sample mean xbar_T
  # psi_grid -- grid of candidate values for the rMAP borrowing parameter psi
  # vague_sd -- standard deviation of the vague component
  # n_tau    -- number of grid points used for numerical integration over tau
  # ncore0   -- number of CPU cores used for parallel computation
  # tol      -- numerical tolerance for identifying a root
  
  ## Prepare the search grid for psi
  psi_grid <- sort(unique(psi_grid))
  n_psi <- length(psi_grid)
  
  if (n_psi < 2L) {
    stop("psi_grid must contain at least two distinct values.")
  }
  
  ## Evaluate the calibration function for a given psi
  evaluate_psi <- function(psi) {
    cal_f_rmap(
      DEL_l = DEL_l,
      DEL_u = DEL_u,
      n_T = n_T,
      n_C = n_C,
      alpha = alpha,
      beta = beta,
      theta0 = theta0,
      theta1 = theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      psi = psi,
      z_C = z_C,
      z_T = z_T,
      w_rob = w_rob,
      vague_sd = vague_sd,
      n_tau = n_tau,
      ncore0 = ncore0
    )
  }
  
  ## Objective function: worst-case power minus target power
  f <- function(psi) {
    evaluate_psi(psi)$res
  }
  
  ## Store only the grid points that are actually evaluated
  psi_checked <- numeric(0)
  res_checked <- numeric(0)
  
  # ===============================================================
  # Evaluate the smallest value of psi
  # ===============================================================
  
  psi_prev <- psi_grid[1L]
  previous_info <- evaluate_psi(psi_prev)
  res_prev <- previous_info$res
  
  if (!is.finite(res_prev)) {
    stop(
      "A non-finite residual is obtained at psi = ",
      psi_prev,
      "."
    )
  }
  
  psi_checked <- c(psi_checked, psi_prev)
  res_checked <- c(res_checked, res_prev)
  
  # ===============================================================
  # Case 1:
  # The residual at the smallest psi is already non-negative
  # ===============================================================
  
  if (res_prev >= -tol) {
    
    psi_use <- psi_prev
    
    root_case <- if (abs(res_prev) <= tol) {
      "minimum_grid_value_is_zero"
    } else {
      "minimum_grid_value_is_nonnegative"
    }
    
    return(list(
      has_root = abs(res_prev) <= tol,
      root_case = root_case,
      psi_root = if (abs(res_prev) <= tol) {
        psi_use
      } else {
        numeric(0)
      },
      root_type = if (abs(res_prev) <= tol) {
        "grid_zero"
      } else {
        "boundary_nonnegative"
      },
      psi_use = psi_use,
      w_rob = w_rob,
      res_use = res_prev,
      p_0 = previous_info$p_0,
      min_pow = previous_info$min_pow,
      muC0_opt = previous_info$muC0_opt,
      crossing_interval = NULL,
      transition_index = 1L,
      psi_checked = psi_checked,
      res_checked = res_checked,
      psi_grid = psi_grid,
      message = paste(
        "The residual at the smallest psi-grid value is non-negative.",
        "Therefore, the smallest psi-grid value is selected."
      )
    ))
  }
  
  ## At this point, the residual at the smallest psi is negative
  
  # ===============================================================
  # Sequentially search for the first negative-to-non-negative
  # transition
  # ===============================================================
  
  for (i in 2L:n_psi) {
    
    psi_curr <- psi_grid[i]
    current_info <- evaluate_psi(psi_curr)
    res_curr <- current_info$res
    
    if (!is.finite(res_curr)) {
      stop(
        "A non-finite residual is obtained at psi = ",
        psi_curr,
        "."
      )
    }
    
    psi_checked <- c(psi_checked, psi_curr)
    res_checked <- c(res_checked, res_curr)
    
    # =============================================================
    # Case 2:
    # The current grid value is numerically equal to zero
    # =============================================================
    
    if (abs(res_curr) <= tol) {
      
      psi_use <- psi_curr
      
      return(list(
        has_root = TRUE,
        root_case = "negative_to_numerical_zero",
        psi_root = psi_use,
        root_type = "grid_zero",
        psi_use = psi_use,
        w_rob = w_rob,
        res_use = res_curr,
        p_0 = current_info$p_0,
        min_pow = current_info$min_pow,
        muC0_opt = current_info$muC0_opt,
        crossing_interval = c(psi_prev, psi_curr),
        transition_index = i,
        psi_checked = psi_checked,
        res_checked = res_checked,
        psi_grid = psi_grid,
        message = paste(
          "The residual starts negative.",
          "The first psi-grid value at which the residual is",
          "numerically equal to zero is selected."
        )
      ))
    }
    
    # =============================================================
    # Case 3:
    # Strict negative-to-positive crossing
    # =============================================================
    
    if (res_prev < -tol && res_curr > tol) {
      
      psi_use <- uniroot(
        f = f,
        interval = c(psi_prev, psi_curr),
        tol = tol
      )$root
      
      ## Evaluate the full calibration output at the root
      root_info <- evaluate_psi(psi_use)
      
      return(list(
        has_root = TRUE,
        root_case = "negative_to_positive_crossing",
        psi_root = psi_use,
        root_type = "uniroot",
        psi_use = psi_use,
        w_rob = w_rob,
        res_use = root_info$res,
        p_0 = root_info$p_0,
        min_pow = root_info$min_pow,
        muC0_opt = root_info$muC0_opt,
        crossing_interval = c(psi_prev, psi_curr),
        transition_index = i - 1L,
        psi_checked = psi_checked,
        res_checked = res_checked,
        psi_grid = psi_grid,
        message = paste0(
          "The first interval over which the residual changes ",
          "from negative to positive is identified. ",
          "The corresponding root is selected: psi = ",
          round(psi_use, 6),
          "."
        )
      ))
    }
    
    ## Move to the next interval
    psi_prev <- psi_curr
    res_prev <- res_curr
    previous_info <- current_info
  }
  
  # ===============================================================
  # Case 4:
  # No negative-to-non-negative transition is found
  # ===============================================================
  
  list(
    has_root = FALSE,
    root_case = "no_negative_to_nonnegative_transition",
    psi_root = numeric(0),
    root_type = NA_character_,
    psi_use = NA_real_,
    w_rob = w_rob,
    res_use = max(res_checked),
    p_0 = NA_real_,
    min_pow = NA_real_,
    muC0_opt = NA_real_,
    crossing_interval = NULL,
    transition_index = NA_integer_,
    psi_checked = psi_checked,
    res_checked = res_checked,
    psi_grid = psi_grid,
    message = paste(
      "The residual is negative at the smallest psi-grid value,",
      "and no transition from negative to non-negative is found",
      "over the specified psi range."
    )
  )
}

#---- Cache calibration results for the rMAP prior ----

## Create an environment for storing previously evaluated calibration results
cal_cache_rmap <- new.env(parent = emptyenv())

## Evaluate the rMAP calibration function using cached results when available
cal_f_cached_rmap <- function(dl, du) {
  
  # dl -- lower bound of the discrepancy range
  # du -- upper bound of the discrepancy range
  
  ## Construct a unique key for the current calibration setting
  key <- paste(
    round(dl, 6),
    round(du, 6),
    round(psi, 6),
    round(w_rob, 6),
    round(vague_sd, 6),
    sep = "_"
  )
  
  ## Return the cached result if this setting has already been evaluated
  if (exists(key, envir = cal_cache_rmap)) {
    return(get(key, envir = cal_cache_rmap))
  }
  
  ## Evaluate the calibration function
  val <- cal_f_rmap(
    DEL_l = dl,
    DEL_u = du,
    n_T = n_T,
    n_C = n_C,
    alpha = alpha,
    beta = beta,
    theta0 = theta0,
    theta1 = theta1,
    sigma2 = sigma2,
    xbar_ch = xbar_ch,
    se_ch = se_ch,
    psi = psi,
    z_C = z_C,
    z_T = z_T,
    w_rob = w_rob,
    vague_sd = vague_sd,
    n_tau = n_tau,
    ncore0 = 1
  )
  
  ## Store the result in the cache
  assign(key, val, envir = cal_cache_rmap)
  
  ## Return the calibration result
  val
}


#---- Compute the EHSS as a function of xbar_C under the rMAP prior ----
ESS_rMAP_xbarC <- function(xbar_C_grid,
                           xbar_ch,
                           n_C,
                           se_ch,
                           sigma2,
                           psi,
                           w_rob = 0.8,
                           vague_sd = 88,
                           weights = NULL,
                           n_tau = 1001) {
  
  # xbar_C_grid -- grid of current control-arm sample means
  # xbar_ch     -- historical control mean
  # n_C         -- control-arm sample size
  # se_ch       -- standard error of the historical control mean;
  #                se_ch = sigma / sqrt(n_ch)
  # sigma2      -- known outcome variance
  # psi         -- scale parameter of the half-normal prior for tau
  # w_rob       -- mixture weight assigned to the MAP component
  # vague_sd    -- standard deviation of the vague component
  # weights     -- optional weights used to compute the mean EHSS
  # n_tau       -- number of grid points used for numerical integration over tau
  
  ## Construct the tau grid and prior integration weights
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau_prior <- tg$w
  
  ## Prior variances of the MAP components
  v0_tau <- se_ch^2 + 2 * tau^2
  
  ## Prior variance of the vague component
  v0_vague <- vague_sd^2
  
  ## Combine the MAP and vague components
  v0_all <- c(v0_tau, v0_vague)
  
  comp_prior_w <- c(
    w_rob * w_tau_prior,
    1 - w_rob
  )
  
  ## Predictive variance of xbar_C conditional on each component
  pred_var_all <- sigma2 / n_C + v0_all
  
  ## Evaluate the predictive density of xbar_C conditional on each component
  dens_all <- sapply(seq_along(v0_all), function(j) {
    dnorm(
      xbar_C_grid,
      mean = xbar_ch,
      sd = sqrt(pred_var_all[j])
    )
  })
  
  ## Compute posterior component weights given xbar_C
  w_comp_given_xbarC <- dens_all * matrix(
    comp_prior_w,
    nrow = length(xbar_C_grid),
    ncol = length(comp_prior_w),
    byrow = TRUE
  )
  
  w_comp_given_xbarC <-
    w_comp_given_xbarC / rowSums(w_comp_given_xbarC)
  
  ## Posterior mean and variance of the control-arm mean
  ## conditional on each mixture component
  prec_C <- n_C / sigma2
  prec_0_all <- 1 / v0_all
  
  v_C_all <- 1 / (prec_C + prec_0_all)
  
  m_C_all <- sapply(seq_along(v0_all), function(j) {
    v_C_all[j] * (
      xbar_C_grid * prec_C +
        xbar_ch * prec_0_all[j]
    )
  })
  
  ## Compute the posterior mixture mean
  m_mix <- rowSums(
    w_comp_given_xbarC * m_C_all
  )
  
  ## Compute the posterior mixture variance
  V_mix <- rowSums(
    w_comp_given_xbarC * (
      m_C_all^2 +
        matrix(
          v_C_all,
          nrow = length(xbar_C_grid),
          ncol = length(v0_all),
          byrow = TRUE
        )
    )
  ) - m_mix^2
  
  ## Compute the EHSS due to historical borrowing
  ESS_grid <- sigma2 / V_mix - n_C
  
  ## Specify weights for computing the mean EHSS
  if (is.null(weights)) {
    weights <- rep(
      1 / length(xbar_C_grid),
      length(xbar_C_grid)
    )
  } else {
    weights <- weights / sum(weights)
  }
  
  ## Return the EHSS and posterior mixture summaries
  list(
    mean_ESS = sum(weights * ESS_grid),
    max_ESS = max(ESS_grid),
    xbarC_at_max_ESS = xbar_C_grid[which.max(ESS_grid)],
    ESS_grid = ESS_grid,
    V_post = V_mix,
    mean_post = m_mix,
    
    tau = tau,
    w_tau_prior = w_tau_prior,
    w_rob = w_rob,
    vague_sd = vague_sd,
    
    v0_tau = v0_tau,
    v0_vague = v0_vague,
    v0_all = v0_all,
    
    v_C_all = v_C_all,
    comp_prior_w = comp_prior_w,
    w_comp_given_xbarC = w_comp_given_xbarC,
    
    map_density = as.vector(
      dens_all[, seq_along(tau), drop = FALSE] %*%
        w_tau_prior
    ),
    vague_density = dens_all[, length(v0_all)],
    prior_density = as.vector(
      dens_all %*% comp_prior_w
    )
  )
}

#---- Compute the expected EHSS given theta_C under the rMAP prior ----
ESS_rMAP_xbarC_given_mu <- function(mu_grid,
                                    xbar_ch,
                                    n_C,
                                    sigma2,
                                    se_ch,
                                    psi,
                                    w_rob = 0.8,
                                    vague_sd = 88,
                                    weights = NULL,
                                    n_tau = 41,
                                    z_grid = seq(-5, 5, length.out = 1001)) {
  
  # mu_grid   -- grid of true control-arm means
  # xbar_ch   -- historical control mean
  # n_C       -- control-arm sample size
  # sigma2    -- known outcome variance
  # se_ch     -- standard error of the historical control mean;
  #              se_ch = sigma / sqrt(n_ch)
  # psi       -- scale parameter of the half-normal prior for tau
  # w_rob     -- mixture weight assigned to the MAP component
  # vague_sd  -- standard deviation of the vague component
  # weights   -- optional weights used to compute the mean expected EHSS
  # n_tau     -- number of grid points used for numerical integration over tau
  # z_grid    -- standard normal grid used for numerical integration over xbar_C
  
  ## Construct the tau grid and prior integration weights
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau_prior <- tg$w
  
  ## Prior variances of the MAP components
  v0_tau <- se_ch^2 + 2 * tau^2
  
  ## Prior variance of the vague component
  v0_vague <- vague_sd^2
  
  ## Combine the MAP and vague components
  v0_all <- c(v0_tau, v0_vague)
  comp_prior_w <- c(
    w_rob * w_tau_prior,
    1 - w_rob
  )
  
  ## Compute the EHSS for a given observed control-arm sample mean
  ESS_one_xbarC <- function(xbar_C) {
    
    ## Predictive variance of xbar_C conditional on each mixture component
    pred_var_all <- sigma2 / n_C + v0_all
    
    ## Compute posterior mixture weights given xbar_C
    log_w <- dnorm(
      xbar_C,
      mean = xbar_ch,
      sd = sqrt(pred_var_all),
      log = TRUE
    ) + log(comp_prior_w)
    
    ## Stabilise and normalise the posterior weights
    log_w <- log_w - max(log_w)
    w_post <- exp(log_w)
    w_post <- w_post / sum(w_post)
    
    ## Posterior mean and variance of the control-arm mean
    ## conditional on each mixture component
    prec_C <- n_C / sigma2
    prec_0_all <- 1 / v0_all
    
    v_C_all <- 1 / (prec_C + prec_0_all)
    
    m_C_all <- v_C_all * (
      xbar_C * prec_C +
        xbar_ch * prec_0_all
    )
    
    ## Compute the posterior mixture mean and variance
    m_mix <- sum(w_post * m_C_all)
    
    V_mix <- sum(
      w_post * (v_C_all + m_C_all^2)
    ) - m_mix^2
    
    ## Compute the EHSS due to historical borrowing
    sigma2 / V_mix - n_C
  }
  
  ## Construct weights for numerical integration over xbar_C | mu_C
  dz <- z_grid[2] - z_grid[1]
  z_w <- dnorm(z_grid) * dz
  z_w <- z_w / sum(z_w)
  
  ## Compute E_{xbar_C | mu_C}[EHSS(xbar_C)] for each true control-arm mean
  ESS_grid <- sapply(mu_grid, function(mu_C) {
    
    xbar_C_grid <- mu_C +
      sqrt(sigma2 / n_C) * z_grid
    
    ESS_xbarC <- sapply(
      xbar_C_grid,
      ESS_one_xbarC
    )
    
    sum(z_w * ESS_xbarC)
  })
  
  ## Specify weights for computing the mean expected EHSS
  if (is.null(weights)) {
    weights <- rep(
      1 / length(mu_grid),
      length(mu_grid)
    )
  } else {
    weights <- weights / sum(weights)
  }
  
  ## Return the expected EHSS summaries
  list(
    mean_ESS = sum(weights * ESS_grid),
    max_ESS = max(ESS_grid),
    muC_at_max_ESS = mu_grid[which.max(ESS_grid)],
    ESS_grid = ESS_grid,
    
    tau = tau,
    w_tau_prior = w_tau_prior,
    w_rob = w_rob,
    vague_sd = vague_sd,
    
    v0_tau = v0_tau,
    v0_vague = v0_vague,
    v0_all = v0_all,
    comp_prior_w = comp_prior_w,
    
    z_grid = z_grid,
    z_weights = z_w
  )
}

#---- Compute the rMAP prior density over a grid of true control-arm means ----
drmap_prior <- function(mu_grid, xbar_ch, se_ch, psi,
                        w_rob = 0.8,
                        vague_sd = 88,
                        n_tau = 1001) {
  
  # mu_grid   -- grid of true control-arm means
  # xbar_ch   -- historical control mean
  # se_ch     -- standard error of the historical control mean;
  #              se_ch = sigma / sqrt(n_ch)
  # psi       -- scale parameter of the half-normal prior for tau
  # w_rob     -- mixture weight assigned to the MAP component
  # vague_sd  -- standard deviation of the vague component
  # n_tau     -- number of grid points used for numerical integration over tau
  
  ## Construct the tau grid and prior integration weights
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau <- tg$w
  
  ## Prior variance of the MAP component conditional on tau
  v_map <- se_ch^2 + 2 * tau^2
  
  ## Evaluate and integrate the MAP prior density over tau
  d_map <- sapply(seq_along(tau), function(j) {
    w_tau[j] * dnorm(
      mu_grid,
      mean = xbar_ch,
      sd = sqrt(v_map[j])
    )
  })
  
  d_map <- rowSums(d_map)
  
  ## Evaluate the vague-component density
  d_vague <- dnorm(
    mu_grid,
    mean = xbar_ch,
    sd = vague_sd
  )
  
  ## Combine the MAP and vague components of the rMAP prior
  d_rmap <- w_rob * d_map +
    (1 - w_rob) * d_vague
  
  ## Return the rMAP density and its weighted components
  data.frame(
    mu = mu_grid,
    rMAP = d_rmap,
    MAP = w_rob * d_map,
    Vague = (1 - w_rob) * d_vague
  )
}


#### Design priors

#---- Compute normal design-prior weights over a specified grid ----
mean_OC_normal_weight <- function(mu_grid, mean_prior, sd_prior) {
  
  # Evaluate the normal design-prior density at each grid point
  w <- dnorm(
    mu_grid,
    mean = mean_prior,
    sd = sd_prior
  )
  
  # Normalise the density values to obtain discrete weights
  w <- w / sum(w)
  
  w
}

#---- Compute real-data prior weights over a grid of true control-arm means ----
real_prior_muC_weights <- function(mu_grid,
                                   mix_w = c(0.51, 0.44, 0.05),
                                   mix_mean = c(-51, -46.8, -54.1),
                                   mix_sd = c(19.9, 7.6, 51.7),
                                   normalise = TRUE) {
  
  # mu_grid   -- grid of true control-arm means
  # mix_w     -- mixture weights of the normal components
  # mix_mean  -- means of the normal mixture components
  # mix_sd    -- standard deviations of the normal mixture components
  # normalise -- whether to normalise the prior density values to sum to one
  
  ## Normalise the mixture weights
  mix_w <- mix_w / sum(mix_w)
  
  ## Evaluate the weighted density of each normal mixture component
  dens_mat <- sapply(seq_along(mix_w), function(j) {
    mix_w[j] * dnorm(
      mu_grid,
      mean = mix_mean[j],
      sd = mix_sd[j]
    )
  })
  
  ## Compute the mixture prior density
  prior_density <- rowSums(dens_mat)
  
  ## Convert the prior density values to discrete weights if requested
  if (normalise) {
    w <- prior_density / sum(prior_density)
  } else {
    w <- prior_density
  }
  
  ## Return the prior weights and mixture component information
  list(
    w = w,
    prior_density = prior_density,
    mix_w = mix_w,
    mix_mean = mix_mean,
    mix_sd = mix_sd,
    component_density = dens_mat
  )
}



#---- Compute the MAP prior weights over a grid of true control-arm means ----
map_prior_muC_weights <- function(mu_grid,
                                  xbar_ch,
                                  se_ch,
                                  psi,
                                  n_tau = 1001,
                                  normalise = TRUE) {
  
  # mu_grid   -- grid of true control-arm means
  # xbar_ch   -- historical control mean
  # se_ch     -- standard error of the historical control mean;
  #              se_ch = sigma / sqrt(n_ch)
  # psi       -- scale parameter of the half-normal prior for tau
  # n_tau     -- number of grid points used for numerical integration over tau
  # normalise -- whether to normalise the prior density values to sum to one
  
  ## Construct the tau grid and prior integration weights
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau <- tg$w
  
  v_ch <- se_ch^2
  
  ## MAP prior conditional on tau:
  ## mu_C | D_Ch, tau ~ N(xbar_ch, se_ch^2 + 2 * tau^2)
  v0_tau <- v_ch + 2 * tau^2
  sd0_tau <- sqrt(v0_tau)
  
  ## Evaluate the conditional prior density over mu_grid for each tau
  dens_mat <- sapply(seq_along(tau), function(j) {
    dnorm(
      mu_grid,
      mean = xbar_ch,
      sd = sd0_tau[j]
    )
  })
  
  ## Integrate over the prior distribution of tau
  prior_density <- as.vector(
    dens_mat %*% w_tau
  )
  
  ## Convert the prior density values to discrete weights if requested
  if (normalise) {
    w <- prior_density / sum(prior_density)
  } else {
    w <- prior_density
  }
  
  list(
    w = w,
    prior_density = prior_density,
    tau = tau,
    w_tau = w_tau,
    prior_sd_tau = sd0_tau,
    prior_var_tau = v0_tau
  )
}



#---- Compute rMAP prior weights over a grid of true control-arm means ----
rmap_prior_muC_weights <- function(mu_grid,
                                   xbar_ch,
                                   se_ch,
                                   sigma2,
                                   psi,
                                   w_rob = 0.8,
                                   vague_sd = 88,
                                   n_tau = 1001,
                                   normalise = TRUE) {
  
  # mu_grid   -- grid of true control-arm means
  # xbar_ch   -- historical control mean
  # se_ch     -- standard error of the historical control mean;
  #              se_ch = sigma / sqrt(n_ch)
  # sigma2    -- known outcome variance
  # psi       -- scale parameter of the half-normal prior for tau
  # w_rob     -- mixture weight assigned to the MAP component
  # vague_sd  -- standard deviation of the vague component
  # n_tau     -- number of grid points used for numerical integration over tau
  # normalise -- whether to normalise the prior density values to sum to one
  
  ## Construct the tau grid and prior integration weights
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau <- tg$w
  
  ## MAP component:
  ## mu_C | D_Ch, tau ~ N(xbar_ch, se_ch^2 + 2 * tau^2)
  v0_tau <- se_ch^2 + 2 * tau^2
  sd0_tau <- sqrt(v0_tau)
  
  ## Evaluate the conditional MAP prior density for each tau
  dens_map_mat <- sapply(seq_along(tau), function(j) {
    dnorm(
      mu_grid,
      mean = xbar_ch,
      sd = sd0_tau[j]
    )
  })
  
  ## Integrate over the prior distribution of tau
  map_density <- as.vector(
    dens_map_mat %*% w_tau
  )
  
  ## Evaluate the vague-component density
  vague_density <- dnorm(
    mu_grid,
    mean = xbar_ch,
    sd = vague_sd
  )
  
  ## Combine the MAP and vague components of the rMAP prior
  prior_density <-
    w_rob * map_density +
    (1 - w_rob) * vague_density
  
  ## Convert the prior density values to discrete weights if requested
  if (normalise) {
    w <- prior_density / sum(prior_density)
  } else {
    w <- prior_density
  }
  
  ## Return the rMAP prior weights and component densities
  list(
    w = w,
    prior_density = prior_density,
    map_density = map_density,
    vague_density = vague_density,
    tau = tau,
    w_tau = w_tau,
    prior_sd_tau = sd0_tau,
    prior_var_tau = v0_tau,
    w_rob = w_rob,
    vague_sd = vague_sd
  )
}