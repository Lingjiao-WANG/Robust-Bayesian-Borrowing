######################################
#### Real data example:           ####
#### Borrow from the control arm  ####
#### Functions                    ####
######################################

#### Fixed power prior ####

B_thr_f <- function(DEL, alpha, alpha_0,
                    n_T, n_C, se_ch,
                    sigma2, theta_d) {
  
  prec_ch <- alpha_0 / se_ch^2
  prec_C  <- n_C / sigma2
  prec_post <- prec_ch + prec_C
  
  w_ch <- prec_ch / prec_post
  
  V  <- sigma2 / n_T + 1 / prec_post
  V0 <- sigma2 / n_T + (prec_C / prec_post)^2 * sigma2 / n_C
  
  z_alp <- qnorm(1 - alpha)
  
  pnorm((theta_d - z_alp * sqrt(V0) - w_ch * DEL) / sqrt(V))
}

min_Pow_f <- function(DEL_l, DEL_u,
                      alpha, alpha_0,
                      n_T, n_C, se_ch,
                      xbar_ch,
                      theta0, theta1, beta,
                      sigma2) {
  
  B_thr <- B_thr_f(
    DEL     = DEL_u,
    alpha   = alpha,
    alpha_0 = alpha_0,
    n_T     = n_T,
    n_C     = n_C,
    se_ch   = se_ch,
    sigma2  = sigma2,
    theta_d = theta0
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

Calibration_f <- function(DEL_l, DEL_u,
                          alpha, beta, alpha_0,
                          n_T, n_C, se_ch,
                          xbar_ch,
                          sigma2, theta0, theta1) {
  
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

calc_TIE_power_pp <- function(mu_C,
                              theta_true,
                              xbar_ch,
                              se_ch,
                              n_T,
                              n_C,
                              sigma2,
                              B_thr,
                              alpha_0) {
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

ppos_box_vec <- function(xbar_T, xbar_C, xbar_ch,
                         n_T, n_C, sigma2, se_ch, delta0, theta0) {
  
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

p_quantile_fast <- function(muC, z_C, z_T, xbar_ch,
                            n_T, n_C, sigma2, se_ch, delta0,
                            theta0, alpha) {
  
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

cpp_theta_box_p <- function(xbar_T, xbar_C, xbar_ch,
                            n_T, n_C, sigma2, se_ch,
                            delta0 = 1) {
  
  v_ch <- se_ch^2
  
  se_pred <- sqrt(sigma2 / n_C + v_ch)
  z_obs <- (xbar_C - xbar_ch) / se_pred
  omega0 <- 2 * (1 - pnorm(abs(z_obs)))
  
  n_ch_eff <- sigma2 / v_ch
  a0 <- omega0 ^ ((n_ch_eff / n_C) ^ delta0)
  a0 <- max(0, min(1, a0))
  
  prec_C <- n_C / sigma2 + a0 / v_ch
  v_C <- 1 / prec_C
  m_C <- v_C * (n_C * xbar_C / sigma2 + a0 * xbar_ch / v_ch)
  
  v_T <- sigma2 / n_T
  m_theta <- xbar_T - m_C
  v_theta <- v_T + v_C
  
  list(
    a0 = a0,
    m = m_theta,
    V = v_theta
  )
}

estimate_OC_cpp_box_p_fast <- function(muC0, n_T, n_C, theta, sigma2,
                                       xbar_ch, se_ch, B_thr, delta0,
                                       z_C, z_T, theta0 = 0) {
  
  v_ch <- se_ch^2
  
  xbar_C <- muC0 + sqrt(sigma2 / n_C) * z_C
  xbar_T <- muC0 + theta + sqrt(sigma2 / n_T) * z_T
  
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
  
  ppos <- pnorm((m_theta - theta0) / sqrt(v_theta))
  
  mean(ppos >= 1 - B_thr)
}

find_global_min_p_quantile <- function(DEL_l, DEL_u,
                                       z_C, z_T, xbar_ch,
                                       n_T, n_C, sigma2, se_ch,
                                       delta0, theta0,
                                       alpha,
                                       grid_len = 11,
                                       ncore = 8) {
  
  lower <- xbar_ch + DEL_l
  upper <- xbar_ch + DEL_u
  
  mu_grid <- seq(lower, upper, length.out = grid_len)
  
  if (ncore > 1) {
    cl <- parallel::makeCluster(ncore)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterExport(cl, c(
      "p_quantile_fast",
      "ppos_box_vec",
      "z_C", "z_T",
      "xbar_ch",
      "n_T", "n_C",
      "sigma2", "se_ch",
      "delta0",
      "theta0",
      "alpha"
    ), envir = environment())
    
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
  
  idx_local <- which(
    q_grid[2:(length(q_grid) - 1)] < q_grid[1:(length(q_grid) - 2)] &
      q_grid[2:(length(q_grid) - 1)] < q_grid[3:length(q_grid)]
  ) + 1
  
  idx_candidate <- unique(c(1, length(mu_grid), idx_local))
  
  local_res <- lapply(idx_candidate, function(j) {
    
    if (j == 1 || j == length(mu_grid)) {
      return(list(
        minimum = mu_grid[j],
        objective = q_grid[j]
      ))
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


cal_f_ad <- function(DEL_l, DEL_u,
                     n_T, n_C,
                     alpha, beta,
                     theta0, theta1,
                     sigma2,
                     xbar_ch, se_ch, delta0,
                     z_C, z_T,
                     ncore0 = 8) {
  
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
  
  muC0_v <- seq(
    from = DEL_l + xbar_ch,
    to = DEL_u + xbar_ch,
    length.out = 10
  )
  
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
  
  min_pow_v <- sapply(muC0_v, obj)
  
  idx_local <- which(
    min_pow_v[2:(length(min_pow_v) - 1)] < min_pow_v[1:(length(min_pow_v) - 2)] &
      min_pow_v[2:(length(min_pow_v) - 1)] < min_pow_v[3:length(min_pow_v)]
  ) + 1
  
  candidate_vals <- min(min_pow_v)
  candidate_mu   <- muC0_v[which.min(min_pow_v)]
  
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
  
  list(
    res = min_pow_o - (1 - beta),
    p_0 = p_0,
    min_pow = min_pow_o,
    muC0_opt = muC0_opt
  )
}


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
  
  f_obj <- function(delta0) {
    evaluate_delta0(delta0)$res
  }
  
  # Store only the grid points that are actually evaluated
  delta_checked <- numeric(0)
  f_checked <- numeric(0)
  
  # ===============================================================
  # Evaluate the smallest delta0
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
  # The smallest delta0 is already non-negative
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
  
  # At this point, the objective starts below zero.
  
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
    
    # Move to the next interval
    delta_prev <- delta_curr
    f_prev <- f_curr
  }
  
  # ===============================================================
  # Case 4:
  # No negative-to-non-negative transition over delta_grid
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

cal_cache <- new.env(parent = emptyenv())

cal_f_cached <- function(dl, du) {
  key <- paste(round(dl, 6), round(du, 6), sep = "_")
  
  if (exists(key, envir = cal_cache)) {
    return(get(key, envir = cal_cache))
  }
  
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
  
  assign(key, val, envir = cal_cache)
  val
}


find_widest_fast <- function(DEL_l_min, DEL_l_max,
                             DEL_u_min, DEL_u_max,
                             delta0,
                             n_l = 15,
                             n_u = 15,
                             ncore = 8) {
  
  dl_grid <- seq(DEL_l_min, DEL_l_max, length.out = n_l)
  
  cl <- parallel::makeCluster(ncore)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterExport(cl, c(
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
  ), envir = environment())
  
  res_list <- parallel::parLapply(cl, dl_grid, function(dl) {
    
    du_grid <- seq(max(dl, DEL_u_min), DEL_u_max, length.out = n_u)
    
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
    
    idx <- which(f_grid[-length(f_grid)] * f_grid[-1] <= 0)
    
    if (length(idx) == 0) return(NULL)
    
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
    
    du_star <- max(roots)
    
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
  
  res_df <- do.call(rbind, res_list)
  
  if (is.null(res_df) || nrow(res_df) == 0) {
    stop("No feasible interval found.")
  }
  
  res_df[which.max(res_df$width), ]
}


find_widest_symmetric_fast <- function(d_min = 0,
                                       d_max,
                                       delta0,
                                       n_d = 15) {
  
  d_grid <- seq(d_min, d_max, length.out = n_d)
  
  f_grid <- sapply(d_grid, function(d) {
    cal_f_ad(
      DEL_l = -d,
      DEL_u =  d,
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
  
  idx <- which(f_grid[-length(f_grid)] * f_grid[-1] <= 0)
  
  if (length(idx) == 0) {
    stop("No feasible symmetric interval found. Increase d_max.")
  }
  
  j <- idx[length(idx)]
  
  d_star <- uniroot(
    function(d) {
      cal_f_ad(
        DEL_l = -d,
        DEL_u =  d,
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
    interval = c(d_grid[j], d_grid[j + 1]),
    tol = 1e-3
  )$root
  
  info <- cal_f_ad(
    DEL_l = -d_star,
    DEL_u =  d_star,
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
    DEL_l = -d_star,
    DEL_u =  d_star,
    width = 2 * d_star,
    p_0 = info$p_0,
    min_pow = info$min_pow,
    muC0_opt = info$muC0_opt,
    res = info$res
  )
}

max_TIE_given_interval <- function(lower, upper, delta0,
                                   n_T, n_C, alpha, theta0,
                                   sigma2, xbar_ch, se_ch,
                                   z_C, z_T, p_0,
                                   grid_len = 11) {
  
  obj_TIE <- function(muC0) {
    estimate_OC_cpp_box_p_fast(
      muC0 = muC0,
      n_T = n_T,
      n_C = n_C,
      theta = theta0,
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
  
  mu_grid <- seq(xbar_ch + lower, xbar_ch + upper, length.out = grid_len)
  
  tie_grid <- sapply(mu_grid, obj_TIE)
  
  idx_local_max <- which(
    tie_grid[2:(length(tie_grid) - 1)] > tie_grid[1:(length(tie_grid) - 2)] &
      tie_grid[2:(length(tie_grid) - 1)] > tie_grid[3:length(tie_grid)]
  ) + 1
  
  candidate_idx <- c(1, length(tie_grid), idx_local_max)
  candidate_vals <- tie_grid[candidate_idx]
  candidate_mu   <- mu_grid[candidate_idx]
  
  if (length(idx_local_max) > 0) {
    
    local_res <- lapply(idx_local_max, function(j) {
      optimise(
        obj_TIE,
        interval = c(mu_grid[j - 1], mu_grid[j + 1]),
        maximum = TRUE,
        tol = 1e-3
      )
    })
    
    local_vals <- sapply(local_res, function(x) x$objective)
    local_mu   <- sapply(local_res, function(x) x$maximum)
    
    candidate_vals <- c(candidate_vals, local_vals)
    candidate_mu   <- c(candidate_mu, local_mu)
  }
  
  k <- which.max(candidate_vals)
  
  list(
    max_TIE = candidate_vals[k],
    muC0_at_max = candidate_mu[k],
    p_0 = p_0,
    mu_grid = mu_grid,
    tie_grid = tie_grid
  )
}


power_range_given_p0 <- function(lower, upper, p_0, delta0,
                                 n_T, n_C, theta1, sigma2,
                                 xbar_ch, se_ch,
                                 z_C, z_T,
                                 theta0 = 0,
                                 grid_len = 11) {
  
  obj_pow <- function(muC0) {
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
  
  lower_abs <- xbar_ch + lower
  upper_abs <- xbar_ch + upper
  
  mu_grid <- seq(lower_abs, upper_abs, length.out = grid_len)
  pow_grid <- sapply(mu_grid, obj_pow)
  
  # ---- min power ----
  j_min <- which.min(pow_grid)
  
  if (j_min == 1 || j_min == length(mu_grid)) {
    min_power <- pow_grid[j_min]
    mu_min <- mu_grid[j_min]
  } else {
    res_min <- optimise(
      obj_pow,
      interval = c(mu_grid[j_min - 1], mu_grid[j_min + 1]),
      maximum = FALSE,
      tol = 1e-3
    )
    min_power <- res_min$objective
    mu_min <- res_min$minimum
  }
  
  # ---- max power ----
  j_max <- which.max(pow_grid)
  
  if (j_max == 1 || j_max == length(mu_grid)) {
    max_power <- pow_grid[j_max]
    mu_max <- mu_grid[j_max]
  } else {
    res_max <- optimise(
      obj_pow,
      interval = c(mu_grid[j_max - 1], mu_grid[j_max + 1]),
      maximum = TRUE,
      tol = 1e-3
    )
    max_power <- res_max$objective
    mu_max <- res_max$maximum
  }
  
  list(
    min_power = min_power,
    muC0_min_power = mu_min,
    max_power = max_power,
    muC0_max_power = mu_max,
    mu_grid = mu_grid,
    pow_grid = pow_grid
  )
}

estimate_OC_cpp_box_p_fast_1 <- function(muC0, n_T, n_C, theta, sigma2,
                                       xbar_ch, se_ch, B_thr, delta0,
                                       z_C, z_T, theta0 = 0) {
  
  v_ch <- se_ch^2
  
  xbar_C <- muC0 + sqrt(sigma2 / n_C) * z_C
  xbar_T <- muC0 + theta + sqrt(sigma2 / n_T) * z_T
  
  se_pred <- sqrt(sigma2 / n_C + v_ch)
  z_obs <- (xbar_C - xbar_ch) / se_pred
  
  omega0 <- 2 * (1 - pnorm(abs(z_obs)))
  
  n_ch_eff <- sigma2 / v_ch
  a0 <- omega0 ^ ((n_ch_eff / n_C) ^ delta0)
  a0 <- pmin(1, pmax(0, a0))
  
  prec_C <- n_C / sigma2 + a0 / v_ch
  v_C <- 1 / prec_C
  
  m_C <- v_C * (
    n_C * xbar_C / sigma2 +
      a0 * xbar_ch / v_ch
  )
  
  v_T <- sigma2 / n_T
  
  m_theta <- xbar_T - m_C
  v_theta <- v_T + v_C
  
  # Posterior probability P(theta < theta0 | data)
  pneg <- pnorm(
    (theta0 - m_theta) / sqrt(v_theta)
  )
  
  # Reject H0 when posterior support for H1 is sufficiently large
  mean(pneg >= 1 - B_thr)
}

mean_OC_normal <- function(mu_grid, TIE_grid, Power_grid, mean_prior, sd_prior) {
  
  # prior weights: muC ~ N(xbar_ch, 100^2), restricted to the region
  w <- dnorm(mu_grid, mean = mean_prior, sd = sd_prior)
  w <- w / sum(w)
  
  data.frame(
    Mean_TIE = sum(w * TIE_grid),
    Mean_Power = sum(w * Power_grid)
  )
}

adaptive_power_prior_weights <- function(mu_grid,
                                         xbar_ch,
                                         n_C,
                                         sigma2,
                                         se_ch,
                                         delta0,
                                         a_min = 1e-4) {
  
  v_ch <- se_ch^2
  n_ch_eff <- sigma2 / v_ch
  
  se_pred <- sqrt(sigma2 / n_C + v_ch)
  
  z_mu <- (mu_grid - xbar_ch) / se_pred
  
  omega0 <- 2 * (1 - pnorm(abs(z_mu)))
  
  a0 <- omega0 ^ ((n_ch_eff / n_C)^delta0)
  a0 <- pmin(1, pmax(a_min, a0))
  
  prior_sd <- sqrt(v_ch / a0)
  
  w <- dnorm(mu_grid, mean = xbar_ch, sd = prior_sd)
  w <- w / sum(w)
  
  list(
    w = w,
    a0_design = a0,
    prior_sd = prior_sd,
    n_ch_eff = n_ch_eff
  )
}


ESS_adaptive_power_prior_xbar <- function(xbar_C_grid,
                                          xbar_ch,
                                          n_C,
                                          sigma2,
                                          se_ch,
                                          delta0,
                                          weights = NULL) {
  
  v_ch <- se_ch^2
  n_ch_eff <- sigma2 / v_ch
  
  se_pred <- sqrt(sigma2 / n_C + v_ch)
  
  z_xbar_C <- (xbar_C_grid - xbar_ch) / se_pred
  
  omega0 <- 2 * (1 - pnorm(abs(z_xbar_C)))
  
  a0 <- omega0 ^ ((n_ch_eff / n_C)^delta0)
  a0 <- pmin(1, pmax(0, a0))
  
  ESS_grid <- a0 * n_ch_eff
  
  if (is.null(weights)) {
    weights <- rep(1 / length(xbar_C_grid), length(xbar_C_grid))
  } else {
    weights <- weights / sum(weights)
  }
  
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

ESS_adaptive_power_prior_muC <- function(mu_grid,
                                         xbar_ch,
                                         n_C,
                                         sigma2,
                                         se_ch,
                                         delta0,
                                         weights = NULL,
                                         z_grid = seq(-5, 5, length.out = 1001)) {
  
  v_ch <- se_ch^2
  n_ch_eff <- sigma2 / v_ch
  
  se_pred <- sqrt(sigma2 / n_C + v_ch)
  
  ## local ESS for one observed xbar_C
  ESS_one_xbarC <- function(xbar_C) {
    
    z_obs <- (xbar_C - xbar_ch) / se_pred
    
    omega0 <- 2 * (1 - pnorm(abs(z_obs)))
    
    a0 <- omega0 ^ ((n_ch_eff / n_C)^delta0)
    a0 <- pmin(1, pmax(0, a0))
    
    a0 * n_ch_eff
  }
  
  ## weights for numerical integration over xbar_C | mu_C
  dz <- z_grid[2] - z_grid[1]
  z_w <- dnorm(z_grid) * dz
  z_w <- z_w / sum(z_w)
  
  ## E_{xbar_C | mu_C}[ESS(xbar_C)]
  ESS_grid <- sapply(mu_grid, function(mu_C) {
    
    xbar_C_grid <- mu_C + sqrt(sigma2 / n_C) * z_grid
    
    ESS_xbarC <- sapply(xbar_C_grid, ESS_one_xbarC)
    
    sum(z_w * ESS_xbarC)
  })
  
  if (is.null(weights)) {
    weights <- rep(1 / length(mu_grid), length(mu_grid))
  } else {
    weights <- weights / sum(weights)
  }
  
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

make_tau_grid <- function(psi, n_tau = 1001, prob = 0.999) {
  tau_max <- psi * qnorm((prob + 1) / 2)
  tau <- seq(0, tau_max, length.out = n_tau)
  
  d_tau <- rep(NA_real_, n_tau)
  d_tau[1] <- tau[2] - tau[1]
  d_tau[n_tau] <- tau[n_tau] - tau[n_tau - 1]
  d_tau[2:(n_tau - 1)] <- (tau[3:n_tau] - tau[1:(n_tau - 2)]) / 2
  
  dens <- sqrt(2 / pi) / psi * exp(-tau^2 / (2 * psi^2))
  w <- dens * d_tau
  w <- w / sum(w)
  
  list(tau = tau, w = w)
}


ppos_map_vec <- function(xbar_T, xbar_C, xbar_ch,
                         n_T, n_C, sigma2, se_ch, psi, theta0,
                         n_tau = 1001) {
  
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w0 <- tg$w
  
  v_ch <- se_ch^2
  
  v0_tau <- v_ch + 2 * tau^2
  vT <- sigma2 / n_T
  
  pred_var <- outer(rep(1, length(xbar_C)), sigma2 / n_C + v0_tau)
  pred_mean <- matrix(xbar_ch, nrow = length(xbar_C), ncol = length(tau))
  
  log_w <- sweep(
    dnorm(xbar_C, mean = pred_mean, sd = sqrt(pred_var), log = TRUE),
    2, log(w0),
    "+"
  )
  
  log_w <- log_w - apply(log_w, 1, max)
  w_post <- exp(log_w)
  w_post <- w_post / rowSums(w_post)
  
  prec_C <- n_C / sigma2
  prec_0 <- 1 / v0_tau
  
  v_C_tau <- 1 / outer(rep(1, length(xbar_C)), prec_C + prec_0)
  m_C_tau <- v_C_tau * (
    xbar_C * prec_C +
      outer(rep(1, length(xbar_C)), xbar_ch * prec_0)
  )
  
  m_theta_tau <- xbar_T - m_C_tau
  v_theta_tau <- vT + v_C_tau
  
  rowSums(
    w_post * pnorm((theta0 - m_theta_tau) / sqrt(v_theta_tau))
  )
}


p_quantile_fast_map <- function(muC, z_C, z_T, xbar_ch,
                                n_T, n_C, sigma2, se_ch, psi,
                                theta0, alpha,
                                n_tau = 1001) {
  
  xbar_C <- muC + sqrt(sigma2 / n_C) * z_C
  xbar_T <- muC + theta0 + sqrt(sigma2 / n_T) * z_T
  
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
  
  as.numeric(quantile(ppos, alpha, names = FALSE))
}


map_theta_p <- function(xbar_T, xbar_C, xbar_ch,
                        n_T, n_C, sigma2, se_ch, psi,
                        n_tau = 1001) {
  
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w0 <- tg$w
  
  v_ch <- se_ch^2
  
  v0_tau <- v_ch + 2 * tau^2
  
  pred_var <- sigma2 / n_C + v0_tau
  log_w <- dnorm(
    xbar_C,
    mean = xbar_ch,
    sd = sqrt(pred_var),
    log = TRUE
  ) + log(w0)
  
  log_w <- log_w - max(log_w)
  w_post <- exp(log_w)
  w_post <- w_post / sum(w_post)
  
  prec_C <- n_C / sigma2
  prec_0 <- 1 / v0_tau
  
  v_C_tau <- 1 / (prec_C + prec_0)
  m_C_tau <- v_C_tau * (xbar_C * prec_C + xbar_ch * prec_0)
  
  v_T <- sigma2 / n_T
  
  m_theta_tau <- xbar_T - m_C_tau
  v_theta_tau <- v_T + v_C_tau
  
  m_mix <- sum(w_post * m_theta_tau)
  V_mix <- sum(w_post * (v_theta_tau + m_theta_tau^2)) - m_mix^2
  
  list(
    tau = tau,
    w_tau = w_post,
    m = m_mix,
    V = V_mix
  )
}


estimate_OC_map_p_fast <- function(muC0, n_T, n_C, theta, sigma2,
                                   xbar_ch, se_ch, B_thr, psi,
                                   z_C, z_T,
                                   theta0 = 0,
                                   n_tau = 1001) {
  
  xbar_C <- muC0 + sqrt(sigma2 / n_C) * z_C
  xbar_T <- muC0 + theta + sqrt(sigma2 / n_T) * z_T
  
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w0 <- tg$w
  
  v_ch <- se_ch^2
  
  v0_tau <- v_ch + 2 * tau^2
  vT <- sigma2 / n_T
  
  pred_var <- outer(rep(1, length(xbar_C)), sigma2 / n_C + v0_tau)
  pred_mean <- matrix(xbar_ch, nrow = length(xbar_C), ncol = length(tau))
  
  log_w <- sweep(
    dnorm(xbar_C, mean = pred_mean, sd = sqrt(pred_var), log = TRUE),
    2, log(w0),
    "+"
  )
  
  log_w <- log_w - apply(log_w, 1, max)
  w_post <- exp(log_w)
  w_post <- w_post / rowSums(w_post)
  
  prec_C <- n_C / sigma2
  prec_0 <- 1 / v0_tau
  
  v_C_tau <- 1 / outer(rep(1, length(xbar_C)), prec_C + prec_0)
  m_C_tau <- v_C_tau * (
    xbar_C * prec_C +
      outer(rep(1, length(xbar_C)), xbar_ch * prec_0)
  )
  
  m_theta_tau <- xbar_T - m_C_tau
  v_theta_tau <- vT + v_C_tau
  
  ppos <- rowSums(
    w_post * pnorm((m_theta_tau - theta0) / sqrt(v_theta_tau))
  )
  
  mean(ppos >= 1 - B_thr)
}

estimate_OC_map_p_fast_1 <- function(muC0, n_T, n_C, theta, sigma2,
                                   xbar_ch, se_ch, B_thr, psi,
                                   z_C, z_T,
                                   theta0 = 0,
                                   n_tau = 1001) {
  
  xbar_C <- muC0 + sqrt(sigma2 / n_C) * z_C
  xbar_T <- muC0 + theta + sqrt(sigma2 / n_T) * z_T
  
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w0 <- tg$w
  
  v_ch <- se_ch^2
  
  v0_tau <- v_ch + 2 * tau^2
  vT <- sigma2 / n_T
  
  pred_var <- outer(
    rep(1, length(xbar_C)),
    sigma2 / n_C + v0_tau
  )
  
  pred_mean <- matrix(
    xbar_ch,
    nrow = length(xbar_C),
    ncol = length(tau)
  )
  
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
  
  log_w <- log_w - apply(log_w, 1, max)
  
  w_post <- exp(log_w)
  w_post <- w_post / rowSums(w_post)
  
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
  
  m_theta_tau <- xbar_T - m_C_tau
  v_theta_tau <- vT + v_C_tau
  
  # Posterior probability P(theta < theta0 | data)
  pneg <- rowSums(
    w_post *
      pnorm(
        (theta0 - m_theta_tau) /
          sqrt(v_theta_tau)
      )
  )
  
  # Reject H0 when posterior support for H1 is sufficiently large
  mean(pneg >= 1 - B_thr)
}

find_global_min_p_quantile_map <- function(DEL_l, DEL_u,
                                           z_C, z_T, xbar_ch,
                                           n_T, n_C, sigma2, se_ch,
                                           psi, theta0,
                                           alpha,
                                           grid_len = 11,
                                           n_tau = 1001,
                                           ncore = 8) {
  
  lower <- xbar_ch + DEL_l
  upper <- xbar_ch + DEL_u
  
  mu_grid <- seq(lower, upper, length.out = grid_len)
  
  if (ncore > 1) {
    cl <- parallel::makeCluster(ncore)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterExport(cl, c(
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
    ), envir = environment())
    
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
  
  idx_local <- which(
    q_grid[2:(length(q_grid) - 1)] < q_grid[1:(length(q_grid) - 2)] &
      q_grid[2:(length(q_grid) - 1)] < q_grid[3:length(q_grid)]
  ) + 1
  
  idx_candidate <- unique(c(1, length(mu_grid), idx_local))
  
  local_res <- lapply(idx_candidate, function(j) {
    
    if (j == 1 || j == length(mu_grid)) {
      return(list(
        minimum = mu_grid[j],
        objective = q_grid[j]
      ))
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


cal_f_map <- function(DEL_l, DEL_u,
                      n_T, n_C,
                      alpha, beta, theta0, theta1,
                      sigma2,
                      xbar_ch, se_ch, psi,
                      z_C, z_T,
                      n_tau = 1001,
                      ncore0 = 8) {
  
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
  
  muC0_v <- seq(
    from = DEL_l + xbar_ch,
    to = DEL_u + xbar_ch,
    length.out = 10
  )
  
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
  
  min_pow_v <- sapply(muC0_v, obj)
  
  idx_local <- which(
    min_pow_v[2:(length(min_pow_v) - 1)] < min_pow_v[1:(length(min_pow_v) - 2)] &
      min_pow_v[2:(length(min_pow_v) - 1)] < min_pow_v[3:length(min_pow_v)]
  ) + 1
  
  candidate_vals <- min(min_pow_v)
  candidate_mu <- muC0_v[which.min(min_pow_v)]
  
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
  
  list(
    res = min_pow_o - (1 - beta),
    p_0 = p_0,
    min_pow = min_pow_o,
    muC0_opt = muC0_opt
  )
}


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
  
  psi_grid <- sort(unique(psi_grid))
  n_psi <- length(psi_grid)
  
  if (n_psi < 2L) {
    stop("psi_grid must contain at least two distinct values.")
  }
  
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
  
  f <- function(psi) {
    evaluate_psi(psi)$res
  }
  
  # Store only the grid points evaluated before stopping
  psi_checked <- numeric(0)
  res_checked <- numeric(0)
  
  # Evaluate the first grid point
  psi_prev <- psi_grid[1L]
  out <- evaluate_psi(psi_prev)
  res_prev <- out$res
  
  if (!is.finite(res_prev)) {
    stop("A non-finite residual is obtained at psi = ", psi_prev, ".")
  }
  
  psi_checked <- c(psi_checked, psi_prev)
  res_checked <- c(res_checked, res_prev)
  
  # ---------------------------------------------------------------
  # Case 1:
  # The smallest psi already gives a non-negative residual.
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
  # Evaluate subsequent grid points until the first numerical zero
  # or the first negative-to-positive crossing is found.
  # ---------------------------------------------------------------
  
  for (i in 2L:n_psi) {
    
    psi_curr <- psi_grid[i]
    out <- evaluate_psi(psi_curr)
    res_curr <- out$res
    
    if (!is.finite(res_curr)) {
      stop("A non-finite residual is obtained at psi = ", psi_curr, ".")
    }
    
    psi_checked <- c(psi_checked, psi_curr)
    res_checked <- c(res_checked, res_curr)
    
    # Case 2:
    # The current grid point is numerically equal to zero.
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
    
    # Case 3:
    # The residual changes strictly from negative to positive.
    if (res_prev < -tol && res_curr > tol) {
      
      psi_use <- uniroot(
        f,
        interval = c(psi_prev, psi_curr),
        tol = tol
      )$root
      
      #out <- evaluate_psi(psi_use)
      
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
    
    # Move to the next interval
    psi_prev <- psi_curr
    res_prev <- res_curr
  }
  
  # ---------------------------------------------------------------
  # Case 4:
  # No negative-to-non-negative transition is found.
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

cal_cache_map <- new.env(parent = emptyenv())

cal_f_cached_map <- function(dl, du) {
  
  key <- paste(
    round(dl, 6),
    round(du, 6),
    round(psi, 6),
    sep = "_"
  )
  
  if (exists(key, envir = cal_cache_map)) {
    return(get(key, envir = cal_cache_map))
  }
  
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
  
  assign(key, val, envir = cal_cache_map)
  val
}


find_widest_fast_map <- function(DEL_l_min, DEL_l_max,
                                 DEL_u_min, DEL_u_max,
                                 psi,
                                 n_l = 15,
                                 n_u = 15,
                                 n_tau = 1001,
                                 ncore = 8) {
  
  dl_grid <- seq(DEL_l_min, DEL_l_max, length.out = n_l)
  
  cl <- parallel::makeCluster(ncore)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterExport(cl, c(
    "make_tau_grid",
    "cal_f_map",
    "find_global_min_p_quantile_map",
    "p_quantile_fast_map",
    "estimate_OC_map_p_fast",
    "ppos_map_vec",
    "n_T", "n_C",
    "alpha", "beta",
    "theta0", "theta1",
    "sigma2", "xbar_ch", "se_ch",
    "psi", "z_C", "z_T",
    "DEL_u_min", "DEL_u_max", "n_u", "n_tau"
  ), envir = environment())
  
  res_list <- parallel::parLapply(cl, dl_grid, function(dl) {
    
    cal_cache_map <- new.env(parent = emptyenv())
    
    cal_f_cached_map <- function(dl, du) {
      
      key <- paste(
        round(dl, 6),
        round(du, 6),
        round(psi, 6),
        sep = "_"
      )
      
      if (exists(key, envir = cal_cache_map)) {
        return(get(key, envir = cal_cache_map))
      }
      
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
      
      assign(key, val, envir = cal_cache_map)
      val
    }
    
    du_grid <- seq(max(dl, DEL_u_min), DEL_u_max, length.out = n_u)
    
    f_grid <- sapply(du_grid, function(du) {
      cal_f_cached_map(dl, du)$res
    })
    
    idx <- which(
      f_grid[-length(f_grid)] * f_grid[-1] <= 0
    )
    
    if (length(idx) == 0) return(NULL)
    
    roots <- sapply(idx, function(j) {
      uniroot(
        function(du) {
          cal_f_cached_map(dl, du)$res
        },
        interval = c(du_grid[j], du_grid[j + 1]),
        tol = 1e-3
      )$root
    })
    
    du_star <- max(roots)
    
    info <- cal_f_cached_map(dl, du_star)
    
    data.frame(
      DEL_l = dl,
      DEL_u = du_star,
      width = du_star - dl,
      psi = psi,
      p_0 = info$p_0,
      min_pow = info$min_pow,
      muC0_opt = info$muC0_opt,
      res = info$res
    )
  })
  
  res_df <- do.call(rbind, res_list)
  
  if (is.null(res_df) || nrow(res_df) == 0) {
    stop("No feasible interval found.")
  }
  
  res_df[which.max(res_df$width), ]
}


find_widest_symmetric_fast_map <- function(d_min = 0,
                                           d_max,
                                           psi,
                                           n_d = 15,
                                           n_tau = 1001) {
  
  d_grid <- seq(d_min, d_max, length.out = n_d)
  
  cal_cache_map <- new.env(parent = emptyenv())
  
  cal_f_cached_map <- function(d) {
    
    key <- paste(round(d, 6), round(psi, 6), sep = "_")
    
    if (exists(key, envir = cal_cache_map)) {
      return(get(key, envir = cal_cache_map))
    }
    
    val <- cal_f_map(
      DEL_l = -d,
      DEL_u =  d,
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
    
    assign(key, val, envir = cal_cache_map)
    val
  }
  
  f_grid <- sapply(d_grid, function(d) {
    cal_f_cached_map(d)$res
  })
  
  idx <- which(
    f_grid[-length(f_grid)] * f_grid[-1] <= 0
  )
  
  if (length(idx) == 0) {
    stop("No feasible symmetric interval found. Increase d_max.")
  }
  
  j <- idx[length(idx)]
  
  d_star <- uniroot(
    function(d) {
      cal_f_cached_map(d)$res
    },
    interval = c(d_grid[j], d_grid[j + 1]),
    tol = 1e-3
  )$root
  
  info <- cal_f_cached_map(d_star)
  
  data.frame(
    DEL_l = -d_star,
    DEL_u =  d_star,
    width = 2 * d_star,
    psi = psi,
    p_0 = info$p_0,
    min_pow = info$min_pow,
    muC0_opt = info$muC0_opt,
    res = info$res
  )
}

max_TIE_given_interval_map <- function(lower, upper, psi,
                                       n_T, n_C, alpha, theta0,
                                       sigma2, xbar_ch, se_ch,
                                       z_C, z_T, p_0,
                                       n_tau = 1001,
                                       grid_len = 11) {
  
  obj_TIE <- function(muC0) {
    estimate_OC_map_p_fast(
      muC0 = muC0,
      n_T = n_T,
      n_C = n_C,
      theta = theta0,
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
  
  mu_grid <- seq(xbar_ch + lower, xbar_ch + upper, length.out = grid_len)
  
  tie_grid <- sapply(mu_grid, obj_TIE)
  
  idx_local_max <- which(
    tie_grid[2:(length(tie_grid) - 1)] > tie_grid[1:(length(tie_grid) - 2)] &
      tie_grid[2:(length(tie_grid) - 1)] > tie_grid[3:length(tie_grid)]
  ) + 1
  
  candidate_idx <- unique(c(1, length(tie_grid), idx_local_max))
  candidate_vals <- tie_grid[candidate_idx]
  candidate_mu <- mu_grid[candidate_idx]
  
  if (length(idx_local_max) > 0) {
    
    local_res <- lapply(idx_local_max, function(j) {
      optimise(
        obj_TIE,
        interval = c(mu_grid[j - 1], mu_grid[j + 1]),
        maximum = TRUE,
        tol = 1e-3
      )
    })
    
    local_vals <- sapply(local_res, function(x) x$objective)
    local_mu <- sapply(local_res, function(x) x$maximum)
    
    candidate_vals <- c(candidate_vals, local_vals)
    candidate_mu <- c(candidate_mu, local_mu)
  }
  
  k <- which.max(candidate_vals)
  
  list(
    max_TIE = candidate_vals[k],
    muC0_at_max = candidate_mu[k],
    p_0 = p_0,
    mu_grid = mu_grid,
    tie_grid = tie_grid
  )
}


power_range_given_p0_map <- function(lower, upper, p_0, psi,
                                     n_T, n_C, theta1, sigma2,
                                     xbar_ch, se_ch,
                                     z_C, z_T,
                                     theta0 = 0,
                                     n_tau = 1001,
                                     grid_len = 11) {
  
  obj_pow <- function(muC0) {
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
  
  lower_abs <- xbar_ch + lower
  upper_abs <- xbar_ch + upper
  
  mu_grid <- seq(lower_abs, upper_abs, length.out = grid_len)
  pow_grid <- sapply(mu_grid, obj_pow)
  
  j_min <- which.min(pow_grid)
  
  if (j_min == 1 || j_min == length(mu_grid)) {
    min_power <- pow_grid[j_min]
    mu_min <- mu_grid[j_min]
  } else {
    res_min <- optimise(
      obj_pow,
      interval = c(mu_grid[j_min - 1], mu_grid[j_min + 1]),
      maximum = FALSE,
      tol = 1e-3
    )
    min_power <- res_min$objective
    mu_min <- res_min$minimum
  }
  
  j_max <- which.max(pow_grid)
  
  if (j_max == 1 || j_max == length(mu_grid)) {
    max_power <- pow_grid[j_max]
    mu_max <- mu_grid[j_max]
  } else {
    res_max <- optimise(
      obj_pow,
      interval = c(mu_grid[j_max - 1], mu_grid[j_max + 1]),
      maximum = TRUE,
      tol = 1e-3
    )
    max_power <- res_max$objective
    mu_max <- res_max$maximum
  }
  
  list(
    min_power = min_power,
    muC0_min_power = mu_min,
    max_power = max_power,
    muC0_max_power = mu_max,
    mu_grid = mu_grid,
    pow_grid = pow_grid
  )
}

mean_OC_normal <- function(mu_grid, TIE_grid, Power_grid, mean_prior, sd_prior) {
  
  # prior weights: muC ~ N(xbar_ch, 100^2), restricted to the region
  w <- dnorm(mu_grid, mean = mean_prior, sd = sd_prior)
  w <- w / sum(w)
  
  data.frame(
    Mean_TIE = sum(w * TIE_grid),
    Mean_Power = sum(w * Power_grid)
  )
}

map_prior_muC_weights <- function(mu_grid,
                                  xbar_ch,
                                  se_ch,
                                  psi,
                                  n_tau = 1001,
                                  normalise = TRUE) {
  
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau <- tg$w
  
  v_ch <- se_ch^2
  
  # MAP prior:
  # mu_C | D_Ch, tau ~ N(xbar_ch, se_ch^2 + 2 * tau^2)
  v0_tau <- v_ch + 2 * tau^2
  sd0_tau <- sqrt(v0_tau)
  
  dens_mat <- sapply(seq_along(tau), function(j) {
    dnorm(mu_grid, mean = xbar_ch, sd = sd0_tau[j])
  })
  
  prior_density <- as.vector(dens_mat %*% w_tau)
  
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


ESS_MAP_xbarC <- function(xbar_C_grid,
                          xbar_ch,
                          n_C,
                          sigma2,
                          se_ch,
                          psi,
                          weights = NULL,
                          n_tau = 1001) {
  
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau_prior <- tg$w
  
  v0_tau <- se_ch^2 + 2 * tau^2
  
  ## predictive variance of xbar_C given tau
  pred_var_tau <- sigma2 / n_C + v0_tau
  
  ## likelihood p(xbar_C | tau)
  dens_mat <- sapply(seq_along(tau), function(j) {
    dnorm(
      xbar_C_grid,
      mean = xbar_ch,
      sd = sqrt(pred_var_tau[j])
    )
  })
  
  ## posterior weights p(tau | xbar_C)
  w_tau_given_xbarC <- dens_mat * matrix(
    w_tau_prior,
    nrow = length(xbar_C_grid),
    ncol = length(tau),
    byrow = TRUE
  )
  
  w_tau_given_xbarC <- w_tau_given_xbarC / rowSums(w_tau_given_xbarC)
  
  ## posterior component distribution for mu_C | xbar_C, tau
  prec_C <- n_C / sigma2
  prec_0 <- 1 / v0_tau
  
  v_C_tau <- 1 / (prec_C + prec_0)
  
  m_C_tau <- sapply(seq_along(tau), function(j) {
    v_C_tau[j] * (
      xbar_C_grid * prec_C +
        xbar_ch * prec_0[j]
    )
  })
  
  ## posterior mixture mean
  m_mix <- rowSums(w_tau_given_xbarC * m_C_tau)
  
  ## posterior mixture variance
  V_mix <- rowSums(
    w_tau_given_xbarC * (m_C_tau^2 + 
                           matrix(v_C_tau,
                                  nrow = length(xbar_C_grid),
                                  ncol = length(tau),
                                  byrow = TRUE))
  ) - m_mix^2
  
  ## local ESS due to historical borrowing
  ESS_grid <- sigma2 / V_mix - n_C
  
  if (is.null(weights)) {
    weights <- rep(1 / length(xbar_C_grid), length(xbar_C_grid))
  } else {
    weights <- weights / sum(weights)
  }
  
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

ESS_MAP_xbarC_given_mu <- function(mu_grid,
                                   xbar_ch,
                                   n_C,
                                   sigma2,
                                   se_ch,
                                   psi,
                                   weights = NULL,
                                   n_tau = 41,
                                   z_grid = seq(-5, 5, length.out = 1001)) {
  
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau_prior <- tg$w
  
  v0_tau <- se_ch^2 + 2 * tau^2
  
  ## function: local ESS for one observed xbar_C
  ESS_one_xbarC <- function(xbar_C) {
    
    pred_var_tau <- sigma2 / n_C + v0_tau
    
    log_w <- dnorm(
      xbar_C,
      mean = xbar_ch,
      sd = sqrt(pred_var_tau),
      log = TRUE
    ) + log(w_tau_prior)
    
    log_w <- log_w - max(log_w)
    w_post <- exp(log_w)
    w_post <- w_post / sum(w_post)
    
    prec_C <- n_C / sigma2
    prec_0 <- 1 / v0_tau
    
    v_C_tau <- 1 / (prec_C + prec_0)
    
    m_C_tau <- v_C_tau * (
      xbar_C * prec_C +
        xbar_ch * prec_0
    )
    
    m_mix <- sum(w_post * m_C_tau)
    
    V_mix <- sum(w_post * (v_C_tau + m_C_tau^2)) - m_mix^2
    
    sigma2 / V_mix - n_C
  }
  
  ## weights for numerical integration over xbar_C | mu_C
  dz <- z_grid[2] - z_grid[1]
  z_w <- dnorm(z_grid) * dz
  z_w <- z_w / sum(z_w)
  
  ## E_{xbar_C | mu_C}[ESS(xbar_C)]
  ESS_grid <- sapply(mu_grid, function(mu_C) {
    
    xbar_C_grid <- mu_C + sqrt(sigma2 / n_C) * z_grid
    
    ESS_xbarC <- sapply(xbar_C_grid, ESS_one_xbarC)
    
    sum(z_w * ESS_xbarC)
  })
  
  if (is.null(weights)) {
    weights <- rep(1 / length(mu_grid), length(mu_grid))
  } else {
    weights <- weights / sum(weights)
  }
  
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

dmap_prior <- function(theta_C, xbar_ch, se_ch, psi, n_tau = 1001) {
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w   <- tg$w
  
  v_ch <- se_ch^2
  v0_tau <- v_ch + 2 * tau^2
  
  dens_mat <- sapply(seq_along(tau), function(j) {
    w[j] * dnorm(theta_C, mean = xbar_ch, sd = sqrt(v0_tau[j]))
  })
  
  rowSums(dens_mat)
}


#### rMAP prior ####

make_tau_grid <- function(psi, n_tau = 1001, prob = 0.999) {
  tau_max <- psi * qnorm((prob + 1) / 2)
  tau <- seq(0, tau_max, length.out = n_tau)
  
  d_tau <- rep(NA_real_, n_tau)
  d_tau[1] <- tau[2] - tau[1]
  d_tau[n_tau] <- tau[n_tau] - tau[n_tau - 1]
  d_tau[2:(n_tau - 1)] <- (tau[3:n_tau] - tau[1:(n_tau - 2)]) / 2
  
  dens <- sqrt(2 / pi) / psi * exp(-tau^2 / (2 * psi^2))
  w <- dens * d_tau
  w <- w / sum(w)
  
  list(tau = tau, w = w)
}


ppos_rmap_vec <- function(xbar_T, xbar_C, xbar_ch,
                          n_T, n_C, sigma2, se_ch, psi, theta0,
                          w_rob = 0.8,
                          vague_sd = 88,
                          n_tau = 1001) {
  
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau0 <- tg$w
  
  v_map <- se_ch^2 + 2 * tau^2
  v_vag <- vague_sd^2
  
  v0 <- c(v_map, v_vag)
  w0 <- c(w_rob * w_tau0, 1 - w_rob)
  
  vT <- sigma2 / n_T
  
  pred_var <- outer(rep(1, length(xbar_C)), sigma2 / n_C + v0)
  pred_mean <- matrix(xbar_ch, nrow = length(xbar_C), ncol = length(v0))
  
  log_w <- sweep(
    dnorm(xbar_C, mean = pred_mean, sd = sqrt(pred_var), log = TRUE),
    2, log(w0),
    "+"
  )
  
  log_w <- log_w - apply(log_w, 1, max)
  w_post <- exp(log_w)
  w_post <- w_post / rowSums(w_post)
  
  prec_C <- n_C / sigma2
  prec_0 <- 1 / v0
  
  v_C <- 1 / outer(rep(1, length(xbar_C)), prec_C + prec_0)
  m_C <- v_C * (
    xbar_C * prec_C +
      outer(rep(1, length(xbar_C)), xbar_ch * prec_0)
  )
  
  m_theta <- xbar_T - m_C
  v_theta <- vT + v_C
  
  rowSums(
    w_post * pnorm((theta0 - m_theta) / sqrt(v_theta))
  )
}


p_quantile_fast_rmap <- function(muC, z_C, z_T, xbar_ch,
                                 n_T, n_C, sigma2, se_ch, psi,
                                 theta0, alpha,
                                 w_rob = 0.8,
                                 vague_sd = 88,
                                 n_tau = 1001) {
  
  xbar_C <- muC + sqrt(sigma2 / n_C) * z_C
  xbar_T <- muC + theta0 + sqrt(sigma2 / n_T) * z_T
  
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
  
  as.numeric(quantile(ppos, alpha, names = FALSE))
}


rmap_theta_p <- function(xbar_T, xbar_C, xbar_ch,
                         n_T, n_C, sigma2, se_ch, psi,
                         w_rob = 0.8,
                         vague_sd = 88,
                         n_tau = 1001) {
  
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau0 <- tg$w
  
  v_map <- se_ch^2 + 2 * tau^2
  v_vag <- vague_sd^2
  
  v0 <- c(v_map, v_vag)
  w0 <- c(w_rob * w_tau0, 1 - w_rob)
  
  pred_var <- sigma2 / n_C + v0
  
  log_w <- dnorm(
    xbar_C,
    mean = xbar_ch,
    sd = sqrt(pred_var),
    log = TRUE
  ) + log(w0)
  
  log_w <- log_w - max(log_w)
  w_post <- exp(log_w)
  w_post <- w_post / sum(w_post)
  
  prec_C <- n_C / sigma2
  prec_0 <- 1 / v0
  
  v_C <- 1 / (prec_C + prec_0)
  m_C <- v_C * (xbar_C * prec_C + xbar_ch * prec_0)
  
  v_T <- sigma2 / n_T
  m_theta <- xbar_T - m_C
  v_theta <- v_T + v_C
  
  m_mix <- sum(w_post * m_theta)
  V_mix <- sum(w_post * (v_theta + m_theta^2)) - m_mix^2
  
  list(
    tau = c(tau, NA_real_),
    component = c(rep("MAP", length(tau)), "vague"),
    w_component = w_post,
    m = m_mix,
    V = V_mix
  )
}


estimate_OC_rmap_p_fast <- function(muC0, n_T, n_C, theta, sigma2,
                                    xbar_ch, se_ch, B_thr, psi,
                                    z_C, z_T,
                                    theta0 = 0,
                                    w_rob = 0.8,
                                    vague_sd = 88,
                                    n_tau = 1001) {
  
  xbar_C <- muC0 + sqrt(sigma2 / n_C) * z_C
  xbar_T <- muC0 + theta + sqrt(sigma2 / n_T) * z_T
  
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau0 <- tg$w
  
  v_map <- se_ch^2 + 2 * tau^2
  v_vag <- vague_sd^2
  
  v0 <- c(v_map, v_vag)
  w0 <- c(w_rob * w_tau0, 1 - w_rob)
  
  vT <- sigma2 / n_T
  
  pred_var <- outer(rep(1, length(xbar_C)), sigma2 / n_C + v0)
  pred_mean <- matrix(xbar_ch, nrow = length(xbar_C), ncol = length(v0))
  
  log_w <- sweep(
    dnorm(xbar_C, mean = pred_mean, sd = sqrt(pred_var), log = TRUE),
    2, log(w0),
    "+"
  )
  
  log_w <- log_w - apply(log_w, 1, max)
  w_post <- exp(log_w)
  w_post <- w_post / rowSums(w_post)
  
  prec_C <- n_C / sigma2
  prec_0 <- 1 / v0
  
  v_C <- 1 / outer(rep(1, length(xbar_C)), prec_C + prec_0)
  m_C <- v_C * (
    xbar_C * prec_C +
      outer(rep(1, length(xbar_C)), xbar_ch * prec_0)
  )
  
  m_theta <- xbar_T - m_C
  v_theta <- vT + v_C
  
  ppos <- rowSums(
    w_post * pnorm((m_theta - theta0) / sqrt(v_theta))
  )
  
  mean(ppos >= 1 - B_thr)
}

estimate_OC_rmap_p_fast_1 <- function(muC0, n_T, n_C, theta, sigma2,
                                    xbar_ch, se_ch, B_thr, psi,
                                    z_C, z_T,
                                    theta0 = 0,
                                    w_rob = 0.8,
                                    vague_sd = 88,
                                    n_tau = 1001) {
  
  xbar_C <- muC0 + sqrt(sigma2 / n_C) * z_C
  xbar_T <- muC0 + theta + sqrt(sigma2 / n_T) * z_T
  
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau0 <- tg$w
  
  v_map <- se_ch^2 + 2 * tau^2
  v_vag <- vague_sd^2
  
  v0 <- c(v_map, v_vag)
  w0 <- c(w_rob * w_tau0, 1 - w_rob)
  
  vT <- sigma2 / n_T
  
  pred_var <- outer(
    rep(1, length(xbar_C)),
    sigma2 / n_C + v0
  )
  
  pred_mean <- matrix(
    xbar_ch,
    nrow = length(xbar_C),
    ncol = length(v0)
  )
  
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
  
  log_w <- log_w - apply(log_w, 1, max)
  w_post <- exp(log_w)
  w_post <- w_post / rowSums(w_post)
  
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
  
  m_theta <- xbar_T - m_C
  v_theta <- vT + v_C
  
  # Posterior probability P(theta < theta0 | data)
  pneg <- rowSums(
    w_post *
      pnorm(
        (theta0 - m_theta) /
          sqrt(v_theta)
      )
  )
  
  # Reject H0 when posterior support for H1 is sufficiently large
  mean(pneg >= 1 - B_thr)
}

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
  
  lower <- xbar_ch + DEL_l
  upper <- xbar_ch + DEL_u
  mu_grid <- seq(lower, upper, length.out = grid_len)
  
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
  
  if (ncore > 1) {
    cl <- parallel::makeCluster(ncore)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterExport(cl, c(
      "make_tau_grid", "p_quantile_fast_rmap", "ppos_rmap_vec",
      "z_C", "z_T", "xbar_ch",
      "n_T", "n_C", "se_ch", "sigma2",
      "psi", "theta0", "alpha",
      "w_rob", "vague_sd", "n_tau"
    ), envir = environment())
    
    q_grid <- unlist(parallel::parLapply(cl, mu_grid, eval_q))
  } else {
    q_grid <- sapply(mu_grid, eval_q)
  }
  
  idx_local <- which(
    q_grid[2:(length(q_grid) - 1)] < q_grid[1:(length(q_grid) - 2)] &
      q_grid[2:(length(q_grid) - 1)] < q_grid[3:length(q_grid)]
  ) + 1
  
  idx_candidate <- unique(c(1, length(mu_grid), idx_local))
  
  local_res <- lapply(idx_candidate, function(j) {
    if (j == 1 || j == length(mu_grid)) {
      return(list(minimum = mu_grid[j], objective = q_grid[j]))
    }
    
    optimise(
      eval_q,
      interval = c(mu_grid[j - 1], mu_grid[j + 1]),
      maximum = FALSE,
      tol = 1e-3
    )
  })
  
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
  
  muC0_v <- seq(xbar_ch + DEL_l, xbar_ch + DEL_u, length.out = 10)
  
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
  
  min_pow_v <- sapply(muC0_v, obj)
  
  idx_local <- which(
    min_pow_v[2:(length(min_pow_v) - 1)] < min_pow_v[1:(length(min_pow_v) - 2)] &
      min_pow_v[2:(length(min_pow_v) - 1)] < min_pow_v[3:length(min_pow_v)]
  ) + 1
  
  candidate_vals <- min(min_pow_v)
  candidate_mu <- muC0_v[which.min(min_pow_v)]
  
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
  
  list(
    res = min_pow_o - (1 - beta),
    p_0 = p_0,
    min_pow = min_pow_o,
    muC0_opt = muC0_opt
  )
}


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
  
  # Search from the largest to the smallest w
  w_grid <- sort(unique(w_grid), decreasing = TRUE)
  n_w <- length(w_grid)
  
  if (n_w < 2L) {
    stop("w_grid must contain at least two distinct values.")
  }
  
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
  
  f <- function(w_rob) {
    evaluate_w(w_rob)$res
  }
  
  # Store only values actually evaluated
  w_checked <- numeric(0)
  res_checked <- numeric(0)
  
  # ===============================================================
  # Evaluate the largest w
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
  # The largest w is already non-negative
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
  
  # At this point, the residual at the largest w is negative.
  
  # ===============================================================
  # Search from large to small for the first transition from
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
      
      # uniroot requires an increasing interval
      w_use <- uniroot(
        f = f,
        interval = c(w_curr, w_prev),
        tol = tol
      )$root
      
      # The root is not generally an evaluated grid point,
      # so evaluate the full output at the root.
      #root_info <- evaluate_w(w_use)
      
      return(list(
        has_root = TRUE,
        root_case = "negative_to_positive_crossing",
        w_root = w_use,
        root_type = "uniroot",
        w_use = w_use,
        psi = psi,
        res_use = current_info$res,
        p_0 = current_info$p_0,
        min_pow = current_info$min_pow,
        muC0_opt = current_info$muC0_opt,
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
    
    # Move to the next, smaller w
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
  
  psi_grid <- sort(unique(psi_grid))
  n_psi <- length(psi_grid)
  
  if (n_psi < 2L) {
    stop("psi_grid must contain at least two distinct values.")
  }
  
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
  
  f <- function(psi) {
    evaluate_psi(psi)$res
  }
  
  # Store only the grid points evaluated before stopping
  psi_checked <- numeric(0)
  res_checked <- numeric(0)
  
  # ===============================================================
  # Evaluate the smallest psi
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
  # The smallest psi is already non-negative
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
  
  # At this point, res_prev < -tol.
  
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
    # The current grid point is numerically equal to zero
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
      
      # This additional evaluation is needed because psi_use is
      # generally not one of the previously evaluated grid points.
      #root_info <- evaluate_psi(psi_use)
      
      return(list(
        has_root = TRUE,
        root_case = "negative_to_positive_crossing",
        psi_root = psi_use,
        root_type = "uniroot",
        psi_use = psi_use,
        w_rob = w_rob,
        res_use = current_info$res,
        p_0 = current_info$p_0,
        min_pow = current_info$min_pow,
        muC0_opt = current_info$muC0_opt,
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
    
    # Move to the next interval
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

cal_cache_rmap <- new.env(parent = emptyenv())

cal_f_cached_rmap <- function(dl, du) {
  
  key <- paste(
    round(dl, 6),
    round(du, 6),
    round(psi, 6),
    round(w_rob, 6),
    round(vague_sd, 6),
    sep = "_"
  )
  
  if (exists(key, envir = cal_cache_rmap)) {
    return(get(key, envir = cal_cache_rmap))
  }
  
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
  
  assign(key, val, envir = cal_cache_rmap)
  val
}

find_widest_fast_rmap <- function(DEL_l_min, DEL_l_max,
                                  DEL_u_min, DEL_u_max,
                                  psi,
                                  w_rob = 0.8,
                                  vague_sd = 88,
                                  n_l = 15,
                                  n_u = 15,
                                  n_tau = 1001,
                                  ncore = 8) {
  
  dl_grid <- seq(DEL_l_min, DEL_l_max, length.out = n_l)
  
  cl <- parallel::makeCluster(ncore)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterExport(cl, c(
    "make_tau_grid",
    "cal_f_rmap",
    "find_global_min_p_quantile_rmap",
    "p_quantile_fast_rmap",
    "estimate_OC_rmap_p_fast",
    "ppos_rmap_vec",
    "n_T", "n_C",
    "alpha", "beta", "theta0", "theta1",
    "sigma2", "xbar_ch", "se_ch",
    "psi", "w_rob", "vague_sd",
    "z_C", "z_T",
    "DEL_u_min", "DEL_u_max", "n_u", "n_tau"
  ), envir = environment())
  
  res_list <- parallel::parLapply(cl, dl_grid, function(dl) {
    
    cal_cache_rmap <- new.env(parent = emptyenv())
    
    cal_f_cached_rmap <- function(dl, du) {
      
      key <- paste(
        round(dl, 6),
        round(du, 6),
        round(psi, 6),
        round(w_rob, 6),
        round(vague_sd, 6),
        sep = "_"
      )
      
      if (exists(key, envir = cal_cache_rmap)) {
        return(get(key, envir = cal_cache_rmap))
      }
      
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
      
      assign(key, val, envir = cal_cache_rmap)
      val
    }
    
    du_grid <- seq(max(dl, DEL_u_min), DEL_u_max, length.out = n_u)
    
    f_grid <- sapply(du_grid, function(du) {
      cal_f_cached_rmap(dl, du)$res
    })
    
    idx <- which(f_grid[-length(f_grid)] * f_grid[-1] <= 0)
    
    if (length(idx) == 0) return(NULL)
    
    roots <- sapply(idx, function(j) {
      uniroot(
        function(du) cal_f_cached_rmap(dl, du)$res,
        interval = c(du_grid[j], du_grid[j + 1]),
        tol = 1e-3
      )$root
    })
    
    du_star <- max(roots)
    info <- cal_f_cached_rmap(dl, du_star)
    
    data.frame(
      DEL_l = dl,
      DEL_u = du_star,
      width = du_star - dl,
      psi = psi,
      w_rob = w_rob,
      vague_sd = vague_sd,
      p_0 = info$p_0,
      min_pow = info$min_pow,
      muC0_opt = info$muC0_opt,
      res = info$res
    )
  })
  
  res_df <- do.call(rbind, res_list)
  
  if (is.null(res_df) || nrow(res_df) == 0) {
    stop("No feasible interval found.")
  }
  
  res_df[which.max(res_df$width), ]
}

find_widest_symmetric_fast_rmap <- function(d_min = 0,
                                            d_max,
                                            psi,
                                            w_rob = 0.8,
                                            vague_sd = 88,
                                            n_d = 15,
                                            n_tau = 1001) {
  
  d_grid <- seq(d_min, d_max, length.out = n_d)
  
  cal_cache_rmap <- new.env(parent = emptyenv())
  
  cal_f_cached_rmap <- function(d) {
    
    key <- paste(
      round(d, 6),
      round(psi, 6),
      round(w_rob, 6),
      round(vague_sd, 6),
      sep = "_"
    )
    
    if (exists(key, envir = cal_cache_rmap)) {
      return(get(key, envir = cal_cache_rmap))
    }
    
    val <- cal_f_rmap(
      DEL_l = -d,
      DEL_u =  d,
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
    
    assign(key, val, envir = cal_cache_rmap)
    val
  }
  
  f_grid <- sapply(d_grid, function(d) {
    cal_f_cached_rmap(d)$res
  })
  
  idx <- which(f_grid[-length(f_grid)] * f_grid[-1] <= 0)
  
  if (length(idx) == 0) {
    stop("No feasible symmetric interval found. Increase d_max.")
  }
  
  j <- idx[length(idx)]
  
  d_star <- uniroot(
    function(d) cal_f_cached_rmap(d)$res,
    interval = c(d_grid[j], d_grid[j + 1]),
    tol = 1e-3
  )$root
  
  info <- cal_f_cached_rmap(d_star)
  
  data.frame(
    DEL_l = -d_star,
    DEL_u =  d_star,
    width = 2 * d_star,
    psi = psi,
    w_rob = w_rob,
    vague_sd = vague_sd,
    p_0 = info$p_0,
    min_pow = info$min_pow,
    muC0_opt = info$muC0_opt,
    res = info$res
  )
}

max_TIE_given_interval_rmap <- function(lower, upper, psi,
                                        n_T, n_C,
                                        alpha, theta0,
                                        sigma2, xbar_ch, se_ch,
                                        z_C, z_T, p_0,
                                        w_rob = 0.8,
                                        vague_sd = 88,
                                        n_tau = 1001,
                                        grid_len = 11) {
  
  obj_TIE <- function(muC0) {
    estimate_OC_rmap_p_fast(
      muC0 = muC0,
      n_T = n_T,
      n_C = n_C,
      theta = theta0,
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
  
  mu_grid <- seq(xbar_ch + lower, xbar_ch + upper, length.out = grid_len)
  tie_grid <- sapply(mu_grid, obj_TIE)
  
  idx_local_max <- which(
    tie_grid[2:(length(tie_grid) - 1)] > tie_grid[1:(length(tie_grid) - 2)] &
      tie_grid[2:(length(tie_grid) - 1)] > tie_grid[3:length(tie_grid)]
  ) + 1
  
  candidate_idx <- unique(c(1, length(tie_grid), idx_local_max))
  candidate_vals <- tie_grid[candidate_idx]
  candidate_mu <- mu_grid[candidate_idx]
  
  if (length(idx_local_max) > 0) {
    
    local_res <- lapply(idx_local_max, function(j) {
      optimise(
        obj_TIE,
        interval = c(mu_grid[j - 1], mu_grid[j + 1]),
        maximum = TRUE,
        tol = 1e-3
      )
    })
    
    local_vals <- sapply(local_res, function(x) x$objective)
    local_mu <- sapply(local_res, function(x) x$maximum)
    
    candidate_vals <- c(candidate_vals, local_vals)
    candidate_mu <- c(candidate_mu, local_mu)
  }
  
  k <- which.max(candidate_vals)
  
  list(
    max_TIE = candidate_vals[k],
    muC0_at_max = candidate_mu[k],
    p_0 = p_0,
    psi = psi,
    w_rob = w_rob,
    vague_sd = vague_sd,
    mu_grid = mu_grid,
    tie_grid = tie_grid
  )
}

power_range_given_p0_rmap <- function(lower, upper, p_0, psi,
                                      n_T, n_C,
                                      theta1, sigma2,
                                      xbar_ch, se_ch,
                                      z_C, z_T,
                                      w_rob = 0.8,
                                      vague_sd = 88,
                                      n_tau = 1001,
                                      grid_len = 11) {
  
  obj_pow <- function(muC0) {
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
  
  lower_abs <- xbar_ch + lower
  upper_abs <- xbar_ch + upper
  
  mu_grid <- seq(lower_abs, upper_abs, length.out = grid_len)
  pow_grid <- sapply(mu_grid, obj_pow)
  
  j_min <- which.min(pow_grid)
  
  if (j_min == 1 || j_min == length(mu_grid)) {
    min_power <- pow_grid[j_min]
    mu_min <- mu_grid[j_min]
  } else {
    res_min <- optimise(
      obj_pow,
      interval = c(mu_grid[j_min - 1], mu_grid[j_min + 1]),
      maximum = FALSE,
      tol = 1e-3
    )
    min_power <- res_min$objective
    mu_min <- res_min$minimum
  }
  
  j_max <- which.max(pow_grid)
  
  if (j_max == 1 || j_max == length(mu_grid)) {
    max_power <- pow_grid[j_max]
    mu_max <- mu_grid[j_max]
  } else {
    res_max <- optimise(
      obj_pow,
      interval = c(mu_grid[j_max - 1], mu_grid[j_max + 1]),
      maximum = TRUE,
      tol = 1e-3
    )
    max_power <- res_max$objective
    mu_max <- res_max$maximum
  }
  
  list(
    min_power = min_power,
    muC0_min_power = mu_min,
    max_power = max_power,
    muC0_max_power = mu_max,
    p_0 = p_0,
    psi = psi,
    w_rob = w_rob,
    vague_sd = vague_sd,
    mu_grid = mu_grid,
    pow_grid = pow_grid
  )
}

mean_OC_normal <- function(mu_grid, TIE_grid, Power_grid, mean_prior, sd_prior) {
  
  # prior weights: muC ~ N(xbar_ch, 100^2), restricted to the region
  w <- dnorm(mu_grid, mean = mean_prior, sd = sd_prior)
  w <- w / sum(w)
  
  data.frame(
    Mean_TIE = sum(w * TIE_grid),
    Mean_Power = sum(w * Power_grid)
  )
}


rmap_prior_muC_weights <- function(mu_grid,
                                   xbar_ch,
                                   se_ch,
                                   sigma2,
                                   psi,
                                   w_rob = 0.8,
                                   vague_sd = 88,
                                   n_tau = 1001,
                                   normalise = TRUE) {
  
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau <- tg$w
  
  # MAP component:
  # mu_C | D_Ch, tau ~ N(xbar_ch, se_ch^2 + 2 * tau^2)
  v0_tau <- se_ch^2 + 2 * tau^2
  sd0_tau <- sqrt(v0_tau)
  
  dens_map_mat <- sapply(seq_along(tau), function(j) {
    dnorm(mu_grid, mean = xbar_ch, sd = sd0_tau[j])
  })
  
  map_density <- as.vector(dens_map_mat %*% w_tau)
  
  # Vague component:
  vague_density <- dnorm(
    mu_grid,
    mean = xbar_ch,
    sd = vague_sd
  )
  
  prior_density <- w_rob * map_density + (1 - w_rob) * vague_density
  
  if (normalise) {
    w <- prior_density / sum(prior_density)
  } else {
    w <- prior_density
  }
  
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

ESS_rMAP_prior <- function(mu_grid,
                           xbar_ch,
                           se_ch,
                           sigma2,
                           psi,
                           w_rob = 0.8,
                           vague_sd = 88,
                           weights = NULL,
                           n_tau = 1001) {
  
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau_prior <- tg$w
  
  # MAP components:
  # historical-control variance is se_ch^2
  v0_tau <- se_ch^2 + 2 * tau^2
  
  # ESS is relative to one current-control observation variance sigma2
  ESS_tau <- sigma2 / v0_tau
  
  dens_map_mat <- sapply(seq_along(tau), function(j) {
    dnorm(mu_grid, mean = xbar_ch, sd = sqrt(v0_tau[j]))
  })
  
  # Vague component
  v0_vague <- vague_sd^2
  ESS_vague <- sigma2 / v0_vague
  
  dens_vague <- dnorm(
    mu_grid,
    mean = xbar_ch,
    sd = vague_sd
  )
  
  comp_prior_w <- c(w_rob * w_tau_prior, 1 - w_rob)
  
  dens_all <- cbind(dens_map_mat, dens_vague)
  
  prior_density <- as.vector(dens_all %*% comp_prior_w)
  
  w_comp_given_mu <- dens_all * matrix(
    comp_prior_w,
    nrow = length(mu_grid),
    ncol = length(comp_prior_w),
    byrow = TRUE
  )
  
  w_comp_given_mu <- w_comp_given_mu / rowSums(w_comp_given_mu)
  
  ESS_comp <- c(ESS_tau, ESS_vague)
  ESS_grid <- as.vector(w_comp_given_mu %*% ESS_comp)
  
  if (is.null(weights)) {
    weights <- rep(1 / length(mu_grid), length(mu_grid))
  } else {
    weights <- weights / sum(weights)
  }
  
  list(
    mean_ESS = sum(weights * ESS_grid),
    max_ESS = max(ESS_grid),
    muC_at_max_ESS = mu_grid[which.max(ESS_grid)],
    ESS_grid = ESS_grid,
    prior_density = prior_density,
    map_density = as.vector(dens_map_mat %*% w_tau_prior),
    vague_density = dens_vague,
    tau = tau,
    w_tau_prior = w_tau_prior,
    w_rob = w_rob,
    vague_sd = vague_sd,
    ESS_tau = ESS_tau,
    ESS_vague = ESS_vague,
    ESS_comp = ESS_comp,
    v0_tau = v0_tau,
    v0_vague = v0_vague,
    w_comp_given_mu = w_comp_given_mu
  )
}

drmap_prior <- function(mu_grid, xbar_ch, se_ch, psi,
                        w_rob = 0.8,
                        vague_sd = 88,
                        n_tau = 1001) {
  
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau <- tg$w
  
  v_map <- se_ch^2 + 2 * tau^2
  
  d_map <- sapply(seq_along(tau), function(j) {
    w_tau[j] * dnorm(mu_grid, mean = xbar_ch, sd = sqrt(v_map[j]))
  })
  
  d_map <- rowSums(d_map)
  
  d_vague <- dnorm(mu_grid, mean = xbar_ch, sd = vague_sd)
  
  d_rmap <- w_rob * d_map + (1 - w_rob) * d_vague
  
  data.frame(
    mu = mu_grid,
    rMAP = d_rmap,
    MAP = w_rob * d_map,
    Vague = (1 - w_rob) * d_vague
  )
}


mean_OC_normal_weight <- function(mu_grid, mean_prior, sd_prior) {
  
  # prior weights: muC ~ N(xbar_ch, 100^2), restricted to the region
  w <- dnorm(mu_grid, mean = mean_prior, sd = sd_prior)
  w <- w / sum(w)
  
  w
}


real_prior_muC_weights <- function(mu_grid,
                                   mix_w = c(0.51, 0.44, 0.05),
                                   mix_mean = c(-51, -46.8, -54.1),
                                   mix_sd = c(19.9, 7.6, 51.7),
                                   normalise = TRUE) {
  
  mix_w <- mix_w / sum(mix_w)
  
  dens_mat <- sapply(seq_along(mix_w), function(j) {
    mix_w[j] * dnorm(mu_grid, mean = mix_mean[j], sd = mix_sd[j])
  })
  
  prior_density <- rowSums(dens_mat)
  
  if (normalise) {
    w <- prior_density / sum(prior_density)
  } else {
    w <- prior_density
  }
  
  list(
    w = w,
    prior_density = prior_density,
    mix_w = mix_w,
    mix_mean = mix_mean,
    mix_sd = mix_sd,
    component_density = dens_mat
  )
}


map_prior_muC_weights <- function(mu_grid,
                                  xbar_ch,
                                  se_ch,
                                  psi,
                                  n_tau = 1001,
                                  normalise = TRUE) {
  
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau <- tg$w
  
  v_ch <- se_ch^2
  
  # MAP prior:
  # mu_C | D_Ch, tau ~ N(xbar_ch, se_ch^2 + 2 * tau^2)
  v0_tau <- v_ch + 2 * tau^2
  sd0_tau <- sqrt(v0_tau)
  
  dens_mat <- sapply(seq_along(tau), function(j) {
    dnorm(mu_grid, mean = xbar_ch, sd = sd0_tau[j])
  })
  
  prior_density <- as.vector(dens_mat %*% w_tau)
  
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

rmap_prior_muC_weights <- function(mu_grid,
                                   xbar_ch,
                                   se_ch,
                                   sigma2,
                                   psi,
                                   w_rob = 0.8,
                                   vague_sd = 88,
                                   n_tau = 1001,
                                   normalise = TRUE) {
  
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau <- tg$w
  
  # MAP component:
  # mu_C | D_Ch, tau ~ N(xbar_ch, se_ch^2 + 2 * tau^2)
  v0_tau <- se_ch^2 + 2 * tau^2
  sd0_tau <- sqrt(v0_tau)
  
  dens_map_mat <- sapply(seq_along(tau), function(j) {
    dnorm(mu_grid, mean = xbar_ch, sd = sd0_tau[j])
  })
  
  map_density <- as.vector(dens_map_mat %*% w_tau)
  
  # Vague component:
  vague_density <- dnorm(
    mu_grid,
    mean = xbar_ch,
    sd = vague_sd
  )
  
  prior_density <- w_rob * map_density + (1 - w_rob) * vague_density
  
  if (normalise) {
    w <- prior_density / sum(prior_density)
  } else {
    w <- prior_density
  }
  
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
  
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau_prior <- tg$w
  
  ## MAP components
  v0_tau <- se_ch^2 + 2 * tau^2
  
  ## Vague component
  v0_vague <- vague_sd^2
  
  ## All prior components: MAP tau-components + vague component
  v0_all <- c(v0_tau, v0_vague)
  
  comp_prior_w <- c(w_rob * w_tau_prior, 1 - w_rob)
  
  ## predictive variance of xbar_C given component
  pred_var_all <- sigma2 / n_C + v0_all
  
  ## likelihood p(xbar_C | component)
  dens_all <- sapply(seq_along(v0_all), function(j) {
    dnorm(
      xbar_C_grid,
      mean = xbar_ch,
      sd = sqrt(pred_var_all[j])
    )
  })
  
  ## posterior weights p(component | xbar_C)
  w_comp_given_xbarC <- dens_all * matrix(
    comp_prior_w,
    nrow = length(xbar_C_grid),
    ncol = length(comp_prior_w),
    byrow = TRUE
  )
  
  w_comp_given_xbarC <- w_comp_given_xbarC / rowSums(w_comp_given_xbarC)
  
  ## posterior component distribution for mu_C | xbar_C, component
  prec_C <- n_C / sigma2
  prec_0_all <- 1 / v0_all
  
  v_C_all <- 1 / (prec_C + prec_0_all)
  
  m_C_all <- sapply(seq_along(v0_all), function(j) {
    v_C_all[j] * (
      xbar_C_grid * prec_C +
        xbar_ch * prec_0_all[j]
    )
  })
  
  ## posterior mixture mean
  m_mix <- rowSums(w_comp_given_xbarC * m_C_all)
  
  ## posterior mixture variance
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
  
  ## local ESS due to borrowing
  ESS_grid <- sigma2 / V_mix - n_C
  
  if (is.null(weights)) {
    weights <- rep(1 / length(xbar_C_grid), length(xbar_C_grid))
  } else {
    weights <- weights / sum(weights)
  }
  
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
    
    map_density = as.vector(dens_all[, seq_along(tau), drop = FALSE] %*% w_tau_prior),
    vague_density = dens_all[, length(v0_all)],
    prior_density = as.vector(dens_all %*% comp_prior_w)
  )
}

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
  
  tg <- make_tau_grid(psi, n_tau = n_tau)
  tau <- tg$tau
  w_tau_prior <- tg$w
  
  ## MAP components
  v0_tau <- se_ch^2 + 2 * tau^2
  
  ## Vague component
  v0_vague <- vague_sd^2
  
  ## All prior components: MAP tau components + vague component
  v0_all <- c(v0_tau, v0_vague)
  comp_prior_w <- c(w_rob * w_tau_prior, 1 - w_rob)
  
  ## function: local ESS for one observed xbar_C
  ESS_one_xbarC <- function(xbar_C) {
    
    pred_var_all <- sigma2 / n_C + v0_all
    
    log_w <- dnorm(
      xbar_C,
      mean = xbar_ch,
      sd = sqrt(pred_var_all),
      log = TRUE
    ) + log(comp_prior_w)
    
    log_w <- log_w - max(log_w)
    w_post <- exp(log_w)
    w_post <- w_post / sum(w_post)
    
    prec_C <- n_C / sigma2
    prec_0_all <- 1 / v0_all
    
    v_C_all <- 1 / (prec_C + prec_0_all)
    
    m_C_all <- v_C_all * (
      xbar_C * prec_C +
        xbar_ch * prec_0_all
    )
    
    m_mix <- sum(w_post * m_C_all)
    
    V_mix <- sum(w_post * (v_C_all + m_C_all^2)) - m_mix^2
    
    sigma2 / V_mix - n_C
  }
  
  ## weights for numerical integration over xbar_C | mu_C
  dz <- z_grid[2] - z_grid[1]
  z_w <- dnorm(z_grid) * dz
  z_w <- z_w / sum(z_w)
  
  ## E_{xbar_C | mu_C}[ESS(xbar_C)]
  ESS_grid <- sapply(mu_grid, function(mu_C) {
    
    xbar_C_grid <- mu_C + sqrt(sigma2 / n_C) * z_grid
    
    ESS_xbarC <- sapply(xbar_C_grid, ESS_one_xbarC)
    
    sum(z_w * ESS_xbarC)
  })
  
  if (is.null(weights)) {
    weights <- rep(1 / length(mu_grid), length(mu_grid))
  } else {
    weights <- weights / sum(weights)
  }
  
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

