######################################
#### Real-data example
#### Sample size determination
######################################

library(quadprog)
library(MCMCpack)
library(doParallel)
library(parallel)

## Source functions
source("case_study_funs.R")

## Load the calibrated analysis-prior results
load("Calibrated_analysis_priors.RData")

#### Settings

## No-borrowing control-arm sample size
nC_nob <- 20

## Save the current graphical settings
old_par <- par(no.readonly = TRUE)

## Specify the range of symmetric discrepancy regions
DEL_u_min <- 0
DEL_u_max <- 25

Del_grid <- seq(
  DEL_u_min,
  DEL_u_max,
  by = 0.01
)

#### Allocation ratio n_T:n_C = 2:1


#---- FPP ----

## Compute the maximum difference between the worst-case and target power
## over the borrowing parameter alpha_0
max_pow_dif_FPP <- function(DEL_l, DEL_u, n_C) {
  
  n_T <- 2 * n_C
  
  optimise(
    Calibration_f,
    interval = c(0, 1),
    DEL_l = DEL_l,
    DEL_u = DEL_u,
    alpha = alpha,
    beta = beta,
    n_T = n_T,
    n_C = n_C,
    se_ch = se_ch,
    xbar_ch = xbar_ch,
    sigma2 = sigma2,
    theta0 = theta0,
    theta1 = theta1,
    maximum = TRUE
  )
}


## Find the calibrated borrowing parameter alpha_0
find_alpha0_FPP <- function(DEL_l, DEL_u, n_C) {
  
  n_T <- 2 * n_C
  
  ## Identify the value of alpha_0 that maximises the power difference
  opt <- max_pow_dif_FPP(
    DEL_l,
    DEL_u,
    n_C
  )
  
  ## Evaluate the power difference under full borrowing
  pow_dif_1 <- Calibration_f(
    DEL_l = DEL_l,
    DEL_u = DEL_u,
    alpha_0 = 1,
    alpha = alpha,
    beta = beta,
    n_T = n_T,
    n_C = n_C,
    se_ch = se_ch,
    xbar_ch = xbar_ch,
    sigma2 = sigma2,
    theta0 = theta0,
    theta1 = theta1
  )
  
  ## Use full borrowing if the target power is already exceeded
  if (pow_dif_1 > 0) {
    return(1)
  }
  
  ## Otherwise, solve for the borrowing parameter that attains the target power
  uniroot(
    Calibration_f,
    interval = c(opt$maximum, 1),
    DEL_l = DEL_l,
    DEL_u = DEL_u,
    alpha = alpha,
    beta = beta,
    n_T = n_T,
    n_C = n_C,
    se_ch = se_ch,
    xbar_ch = xbar_ch,
    sigma2 = sigma2,
    theta0 = theta0,
    theta1 = theta1
  )$root
}


## Determine the minimum required sample size under the FPP
min_sample_size_FPP <- function(DEL_l, DEL_u, nC_nob) {
  
  n_C <- nC_nob
  
  ## Sequentially reduce n_C while the target power remains attainable
  repeat {
    
    pow_dif <- max_pow_dif_FPP(
      DEL_l,
      DEL_u,
      n_C
    )$objective
    
    if (pow_dif <= 0 || n_C <= 2) {
      break
    }
    
    n_C <- n_C - 1
  }
  
  ## Return to the smallest feasible control-arm sample size
  n_C <- n_C + 1
  n_T <- 2 * n_C
  
  ## Calibrate the borrowing parameter at the selected sample size
  alpha_0 <- find_alpha0_FPP(
    DEL_l,
    DEL_u,
    n_C
  )
  
  list(
    n_C = n_C,
    n_T = n_T,
    alpha_0 = alpha_0
  )
}


## Compute the difference between the minimum required n_C
## and a specified candidate sample size
Del_f <- function(Del, n_C, nC_nob) {
  
  min_sample_size_FPP(
    DEL_l = -Del,
    DEL_u = Del,
    nC_nob = nC_nob
  )$n_C - n_C
}


## Determine the minimum and maximum required control-arm sample sizes
## over the specified discrepancy range
n_C_min_FPP <- min_sample_size_FPP(
  DEL_l = -DEL_u_min,
  DEL_u = DEL_u_min,
  nC_nob = nC_nob
)$n_C

n_C_max_FPP <- min_sample_size_FPP(
  DEL_l = -DEL_u_max,
  DEL_u = DEL_u_max,
  nC_nob = nC_nob
)$n_C

## Sequence of candidate control-arm sample sizes
nC_seq <- n_C_min_FPP:n_C_max_FPP


## Find the largest discrepancy bound supported by a given n_C
find_Del_max <- function(n_C, Del_grid) {
  
  lo <- 1
  hi <- length(Del_grid)
  
  ## Return the largest grid value if it remains feasible
  if (Del_f(
    Del_grid[hi],
    n_C = n_C,
    nC_nob = nC_nob
  ) <= 0) {
    return(Del_grid[hi])
  }
  
  ## Return NA if even the smallest discrepancy value is infeasible
  if (Del_f(
    Del_grid[lo],
    n_C = n_C,
    nC_nob = nC_nob
  ) > 0) {
    return(NA_real_)
  }
  
  ## Use binary search to identify the largest feasible discrepancy value
  while (hi - lo > 1) {
    
    mid <- floor((lo + hi) / 2)
    
    if (Del_f(
      Del_grid[mid],
      n_C = n_C,
      nC_nob = nC_nob
    ) <= 0) {
      lo <- mid
    } else {
      hi <- mid
    }
  }
  
  Del_grid[lo]
}


## Determine the largest feasible discrepancy bound for each n_C
Del_max <- sapply(
  nC_seq,
  find_Del_max,
  Del_grid = Del_grid
)

## Summarise the sample size results
out_FPP_summary <- data.frame(
  n_C = nC_seq,
  n_T = 2 * nC_seq,
  Del_max = Del_max
)


#---- APP ----

## Find the first borrowing parameter delta0 that yields power above the target
find_plausible_delta0 <- function(DEL_l, DEL_u,
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
  # se_ch      -- standard error of the historical control mean
  # z_C        -- simulated standard normal values for the control arm
  # z_T        -- simulated standard normal values for the treatment arm
  # delta0_min -- lower bound of the search range for delta0
  # delta0_max -- upper bound of the search range for delta0
  # n_grid     -- number of grid points used in the search
  # ncore0     -- number of CPU cores used for parallel computation
  # tol        -- numerical tolerance for determining feasibility
  
  ## Construct the search grid for delta0
  delta_grid <- seq(
    delta0_min,
    delta0_max,
    length.out = n_grid
  )
  
  delta_grid <- sort(unique(delta_grid))
  n_delta <- length(delta_grid)
  
  if (n_delta < 1L) {
    stop("delta_grid must contain at least one value.")
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
  
  ## Store only the grid values that are actually evaluated
  delta_checked <- numeric(0)
  f_checked <- numeric(0)
  
  ## Search sequentially for the first delta0 with residual greater than tol
  for (i in seq_len(n_delta)) {
    
    delta_curr <- delta_grid[i]
    current_info <- evaluate_delta0(delta_curr)
    f_curr <- current_info$res
    
    if (!is.finite(f_curr)) {
      stop(
        "A non-finite objective value is obtained at delta0 = ",
        delta_curr,
        "."
      )
    }
    
    delta_checked <- c(delta_checked, delta_curr)
    f_checked <- c(f_checked, f_curr)
    
    ## Select the first grid value for which the target power is exceeded
    if (f_curr > tol) {
      
      root_case <- if (i == 1L) {
        "minimum_grid_value_is_positive"
      } else {
        "first_positive_grid_value"
      }
      
      return(list(
        found_positive = TRUE,
        root_case = root_case,
        delta0 = delta_curr,
        delta0_use = delta_curr,
        root_info = current_info,
        res_use = f_curr,
        f_at_delta0 = f_curr,
        p_0 = current_info$p_0,
        min_pow = current_info$min_pow,
        muC0_opt = current_info$muC0_opt,
        transition_index = i,
        delta_checked = delta_checked,
        f_checked = f_checked,
        delta_grid = delta_grid,
        message = paste0(
          "The first delta0-grid value satisfying res > tol ",
          "is selected: delta0 = ",
          round(delta_curr, 6),
          ", res = ",
          signif(f_curr, 6),
          "."
        )
      ))
    }
  }
  
  ## No feasible positive residual is found over the search range
  list(
    found_positive = FALSE,
    root_case = "no_positive_grid_value",
    delta0 = NA_real_,
    delta0_use = NA_real_,
    root_info = NULL,
    res_use = max(f_checked),
    f_at_delta0 = NA_real_,
    p_0 = NA_real_,
    min_pow = NA_real_,
    muC0_opt = NA_real_,
    transition_index = NA_integer_,
    delta_checked = delta_checked,
    f_checked = f_checked,
    delta_grid = delta_grid,
    message = paste(
      "No delta0-grid value satisfying res > tol was found",
      "over the specified delta0 range."
    )
  )
}


## Find the calibrated borrowing parameter delta0 for a given sample size
find_delta0_APP <- function(DEL_l, DEL_u, n_C) {
  
  n_T <- 2 * n_C
  
  out <- find_delta0_roots_cal_f_ad(
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
    z_C = z_C,
    z_T = z_T,
    delta0_min = 0,
    delta0_max = 3,
    n_grid = 201,
    ncore0 = 8,
    tol = 1e-6
  )
  
  out$delta0
}


## Check whether a feasible calibrated delta0 exists
is_feasible_APP <- function(DEL_l, DEL_u, n_C) {
  
  delta0 <- find_delta0_APP(
    DEL_l = DEL_l,
    DEL_u = DEL_u,
    n_C = n_C
  )
  
  length(delta0) > 0 && any(is.finite(delta0))
}


#---- Determine the maximum required n_C over the discrepancy range ----

n_C <- nC_nob
n_T <- 2 * n_C

## Evaluate the power difference at delta0 = 0
pow_diff_0 <- cal_f_ad(
  DEL_l = -DEL_u_max,
  DEL_u = DEL_u_max,
  n_T = n_T,
  n_C = n_C,
  alpha = alpha,
  beta = beta,
  theta0 = theta0,
  theta1 = theta1,
  sigma2 = sigma2,
  xbar_ch = xbar_ch,
  se_ch = se_ch,
  delta0 = 0,
  z_C = z_C,
  z_T = z_T,
  ncore0 = 8
)$res

## Reduce n_C until delta0 = 0 no longer attains the target power
while (pow_diff_0 >= 0) {
  
  n_C <- n_C - 1
  n_T <- 2 * n_C
  
  pow_diff_0 <- cal_f_ad(
    DEL_l = -DEL_u_max,
    DEL_u = DEL_u_max,
    n_T = n_T,
    n_C = n_C,
    alpha = alpha,
    beta = beta,
    theta0 = theta0,
    theta1 = theta1,
    sigma2 = sigma2,
    xbar_ch = xbar_ch,
    se_ch = se_ch,
    delta0 = 0,
    z_C = z_C,
    z_T = z_T,
    ncore0 = 8
  )$res
}

## Search for the minimum feasible n_C allowing some delta0 in the search range
n_T <- 2 * n_C

delta0 <- find_plausible_delta0(
  DEL_l = -DEL_u_max,
  DEL_u = DEL_u_max,
  n_T = n_T,
  n_C = n_C,
  alpha = alpha,
  beta = beta,
  theta0 = theta0,
  theta1 = theta1,
  sigma2 = sigma2,
  xbar_ch = xbar_ch,
  se_ch = se_ch,
  z_C = z_C,
  z_T = z_T,
  delta0_min = 0,
  delta0_max = 3,
  n_grid = 201,
  ncore0 = 8,
  tol = 1e-6
)

feasible <- length(delta0$delta0_use) > 0 &&
  any(is.finite(delta0$delta0_use))

if (feasible) {
  
  ## Reduce n_C until feasibility is lost
  while (feasible) {
    
    n_C <- n_C - 1
    n_T <- 2 * n_C
    
    delta0 <- find_plausible_delta0(
      DEL_l = -DEL_u_max,
      DEL_u = DEL_u_max,
      n_T = n_T,
      n_C = n_C,
      alpha = alpha,
      beta = beta,
      theta0 = theta0,
      theta1 = theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      z_C = z_C,
      z_T = z_T,
      delta0_min = 0,
      delta0_max = 3,
      n_grid = 201,
      ncore0 = 8,
      tol = 1e-6
    )
    
    feasible <- length(delta0$delta0_use) > 0 &&
      any(is.finite(delta0$delta0_use))
  }
  
  ## Return to the smallest feasible sample size
  n_C <- n_C + 1
  
} else {
  
  ## Increase n_C until feasibility is achieved
  while (!feasible) {
    
    n_C <- n_C + 1
    n_T <- 2 * n_C
    
    delta0 <- find_plausible_delta0(
      DEL_l = -DEL_u_max,
      DEL_u = DEL_u_max,
      n_T = n_T,
      n_C = n_C,
      alpha = alpha,
      beta = beta,
      theta0 = theta0,
      theta1 = theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      z_C = z_C,
      z_T = z_T,
      delta0_min = 0,
      delta0_max = 3,
      n_grid = 201,
      ncore0 = 8,
      tol = 1e-6
    )
    
    feasible <- length(delta0$delta0_use) > 0 &&
      any(is.finite(delta0$delta0_use))
  }
}

n_C_max_APP <- n_C


#---- Determine the minimum required n_C at DEL_u_min ----

n_C <- nC_nob
n_T <- 2 * n_C

pow_diff_0 <- cal_f_ad(
  DEL_l = -DEL_u_min,
  DEL_u = DEL_u_min,
  n_T = n_T,
  n_C = n_C,
  alpha = alpha,
  beta = beta,
  theta0 = theta0,
  theta1 = theta1,
  sigma2 = sigma2,
  xbar_ch = xbar_ch,
  se_ch = se_ch,
  delta0 = 0,
  z_C = z_C,
  z_T = z_T,
  ncore0 = 8
)$res

## Reduce n_C until the target power is no longer attained
while (pow_diff_0 >= 0) {
  
  n_C <- n_C - 1
  n_T <- 2 * n_C
  
  pow_diff_0 <- cal_f_ad(
    DEL_l = -DEL_u_min,
    DEL_u = DEL_u_min,
    n_T = n_T,
    n_C = n_C,
    alpha = alpha,
    beta = beta,
    theta0 = theta0,
    theta1 = theta1,
    sigma2 = sigma2,
    xbar_ch = xbar_ch,
    se_ch = se_ch,
    delta0 = 0,
    z_C = z_C,
    z_T = z_T,
    ncore0 = 8
  )$res
}

## Smallest feasible control-arm sample size
n_C_min_APP <- n_C + 1

## Sequence of feasible control-arm sample sizes
nC_seq_APP <- n_C_min_APP:n_C_max_APP

#---- Find the largest discrepancy bound for each feasible n_C ----

find_Del_max_APP <- function(n_C, Del_grid) {
  
  n_T <- 2 * n_C
  n_grid <- length(Del_grid)
  
  ## Evaluate the power difference at delta0 = 0
  pow_diff_0 <- function(i) {
    cal_f_ad(
      DEL_l = -Del_grid[i],
      DEL_u = Del_grid[i],
      n_T = n_T,
      n_C = n_C,
      alpha = alpha,
      beta = beta,
      theta0 = theta0,
      theta1 = theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      delta0 = 0,
      z_C = z_C,
      z_T = z_T,
      ncore0 = 8
    )$res
  }
  
  ## Check whether any delta0 in the search range yields a feasible design
  is_feasible <- function(i) {
    
    delta0 <- find_plausible_delta0(
      DEL_l = -Del_grid[i],
      DEL_u = Del_grid[i],
      n_T = n_T,
      n_C = n_C,
      alpha = alpha,
      beta = beta,
      theta0 = theta0,
      theta1 = theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      z_C = z_C,
      z_T = z_T,
      delta0_min = 0,
      delta0_max = 3,
      n_grid = 201,
      ncore0 = 8,
      tol = 1e-6
    )
    
    feasible <- length(delta0$delta0_use) > 0 &&
      any(is.finite(delta0$delta0_use))
    
    return(feasible)
  }
  
  ## Return NA if delta0 = 0 remains feasible at the largest discrepancy
  if (pow_diff_0(n_grid) > 0) {
    return(NA_real_)
  }
  
  ## Locate the first grid value at which delta0 = 0 becomes infeasible
  left <- 1L
  right <- n_grid
  
  while (right - left > 1L) {
    
    mid <- floor((left + right) / 2)
    
    if (pow_diff_0(mid) > 0) {
      left <- mid
    } else {
      right <- mid
    }
  }
  
  lo <- if (pow_diff_0(left) <= 0) left else right
  
  ## Check feasibility after allowing delta0 to vary
  feasible <- is_feasible(lo)
  
  if (feasible) {
    
    hi <- length(Del_grid)
    
    ## Return the largest grid value if it remains feasible
    if (is_feasible(hi)) {
      return(Del_grid[hi])
    }
    
    ## Binary search for the largest feasible discrepancy bound
    while (hi - lo > 1L) {
      
      mid <- floor((lo + hi) / 2L)
      
      if (is_feasible(mid)) {
        lo <- mid
      } else {
        hi <- mid
      }
    }
    
  } else {
    
    lo <- lo - 1L
  }
  
  ## Return the largest feasible discrepancy bound
  Del_grid[lo]
}



## Determine the largest feasible discrepancy bound for each n_C
Del_max_APP <- find_Del_max_APP(
  nC_seq_APP[1],
  Del_grid
)

for (i in 2:length(nC_seq_APP)) {
  
  Del_max_APP[i] <- find_Del_max_APP(
    nC_seq_APP[i],
    Del_grid[Del_grid > Del_max_APP[i - 1]]
  )
}

## Summarise the APP sample size results
out_APP_summary <- data.frame(
  n_C = nC_seq_APP,
  n_T = 2 * nC_seq_APP,
  Del_max = Del_max_APP
)

#---- MAP ----

## Find the first borrowing parameter psi that yields power above the target
find_plausible_psi <- function(DEL_l, DEL_u,
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
  # se_ch    -- standard error of the historical control mean
  # z_C      -- simulated standard normal values for the control arm
  # z_T      -- simulated standard normal values for the treatment arm
  # psi_grid -- grid of candidate values for psi
  # n_tau    -- number of grid points used for numerical integration over tau
  # ncore0   -- number of CPU cores used for parallel computation
  # tol      -- numerical tolerance for determining feasibility
  
  ## Prepare the search grid for psi
  psi_grid <- sort(unique(psi_grid))
  n_psi <- length(psi_grid)
  
  if (n_psi < 1L) {
    stop("psi_grid must contain at least one value.")
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
  
  ## Store only the grid values that are actually evaluated
  psi_checked <- numeric(0)
  res_checked <- numeric(0)
  
  ## Search sequentially for the first psi with residual greater than tol
  for (i in seq_len(n_psi)) {
    
    psi_curr <- psi_grid[i]
    out_curr <- evaluate_psi(psi_curr)
    res_curr <- out_curr$res
    
    if (!is.finite(res_curr)) {
      stop(
        "A non-finite residual is obtained at psi = ",
        psi_curr,
        "."
      )
    }
    
    psi_checked <- c(psi_checked, psi_curr)
    res_checked <- c(res_checked, res_curr)
    
    ## Select the first grid value for which the target power is exceeded
    if (res_curr > tol) {
      
      root_case <- if (i == 1L) {
        "minimum_grid_value_is_positive"
      } else {
        "first_positive_grid_value"
      }
      
      return(list(
        found_positive = TRUE,
        root_case = root_case,
        psi_use = psi_curr,
        res_use = res_curr,
        p_0 = out_curr$p_0,
        min_pow = out_curr$min_pow,
        muC0_opt = out_curr$muC0_opt,
        transition_index = i,
        psi_checked = psi_checked,
        res_checked = res_checked,
        psi_grid = psi_grid,
        message = paste0(
          "The first psi-grid value satisfying res > tol is selected: ",
          "psi = ", round(psi_curr, 6),
          ", res = ", signif(res_curr, 6), "."
        )
      ))
    }
  }
  
  ## No feasible positive residual is found over the search range
  list(
    found_positive = FALSE,
    root_case = "no_positive_grid_value",
    psi_use = NA_real_,
    res_use = max(res_checked),
    p_0 = NA_real_,
    min_pow = NA_real_,
    muC0_opt = NA_real_,
    transition_index = NA_integer_,
    psi_checked = psi_checked,
    res_checked = res_checked,
    psi_grid = psi_grid,
    message = paste(
      "No psi-grid value satisfying res > tol was found",
      "over the specified psi range."
    )
  )
}


## Find the calibrated borrowing parameter psi for a given sample size
find_psi_MAP <- function(DEL_l, DEL_u, n_C) {
  
  n_T <- 2 * n_C
  
  psi_res_map <- find_all_psi(
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
    z_C = z_C,
    z_T = z_T,
    psi_grid = seq(
      0.001,
      3 * sigma,
      length.out = 30
    ),
    n_tau = 401,
    ncore0 = 4
  )$psi_use
  
  psi_res_map
}


#---- Determine the maximum required n_C over the discrepancy range ----

n_C <- nC_nob
n_T <- 2 * n_C

## Evaluate the power difference at the smallest value of psi
pow_diff_MAP <- cal_f_map(
  DEL_l = -DEL_u_max,
  DEL_u = DEL_u_max,
  n_T = n_T,
  n_C = n_C,
  alpha = alpha,
  beta = beta,
  theta0 = theta0,
  theta1 = theta1,
  sigma2 = sigma2,
  xbar_ch = xbar_ch,
  se_ch = se_ch,
  psi = 0.001,
  z_C = z_C,
  z_T = z_T,
  n_tau = 401,
  ncore0 = 4
)$res

## Reduce n_C while the target power remains attainable
while (pow_diff_MAP >= 0) {
  
  n_C <- n_C - 1
  n_T <- 2 * n_C
  
  pow_diff_MAP <- cal_f_map(
    DEL_l = -DEL_u_max,
    DEL_u = DEL_u_max,
    n_T = n_T,
    n_C = n_C,
    alpha = alpha,
    beta = beta,
    theta0 = theta0,
    theta1 = theta1,
    sigma2 = sigma2,
    xbar_ch = xbar_ch,
    se_ch = se_ch,
    psi = 0.001,
    z_C = z_C,
    z_T = z_T,
    n_tau = 401,
    ncore0 = 4
  )$res
}

## Search for the minimum feasible n_C allowing psi to vary
psi_res_map <- find_plausible_psi(
  DEL_l = -DEL_u_max,
  DEL_u = DEL_u_max,
  n_T = n_T,
  n_C = n_C,
  alpha = alpha,
  beta = beta,
  theta0 = theta0,
  theta1 = theta1,
  sigma2 = sigma2,
  xbar_ch = xbar_ch,
  se_ch = se_ch,
  z_C = z_C,
  z_T = z_T,
  psi_grid = seq(
    0.001,
    3 * sigma,
    length.out = 30
  ),
  n_tau = 401,
  ncore0 = 4
)

feasible <- is.finite(psi_res_map$psi_use)

if (feasible) {
  
  ## Reduce n_C until feasibility is lost
  while (feasible) {
    
    n_C <- n_C - 1
    n_T <- 2 * n_C
    
    psi_res_map <- find_plausible_psi(
      DEL_l = -DEL_u_max,
      DEL_u = DEL_u_max,
      n_T = n_T,
      n_C = n_C,
      alpha = alpha,
      beta = beta,
      theta0 = theta0,
      theta1 = theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      z_C = z_C,
      z_T = z_T,
      psi_grid = seq(
        0.001,
        3 * sigma,
        length.out = 30
      ),
      n_tau = 401,
      ncore0 = 4
    )
    
    feasible <- is.finite(psi_res_map$psi_use)
  }
  
  ## Return to the smallest feasible sample size
  n_C <- n_C + 1
  
} else {
  
  ## Increase n_C until feasibility is achieved
  while (!feasible) {
    
    n_C <- n_C + 1
    n_T <- 2 * n_C
    
    psi_res_map <- find_plausible_psi(
      DEL_l = -DEL_u_max,
      DEL_u = DEL_u_max,
      n_T = n_T,
      n_C = n_C,
      alpha = alpha,
      beta = beta,
      theta0 = theta0,
      theta1 = theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      z_C = z_C,
      z_T = z_T,
      psi_grid = seq(
        0.001,
        3 * sigma,
        length.out = 30
      ),
      n_tau = 401,
      ncore0 = 4
    )
    
    feasible <- is.finite(psi_res_map$psi_use)
  }
}

n_C_max_MAP <- n_C


#---- Determine the minimum required n_C at DEL_u_min ----

## Start from the maximum required sample size
n_T <- 2 * n_C

pow_diff_MAP0 <- cal_f_map(
  DEL_l = -DEL_u_min,
  DEL_u = DEL_u_min,
  n_T = n_T,
  n_C = n_C,
  alpha = alpha,
  beta = beta,
  theta0 = theta0,
  theta1 = theta1,
  sigma2 = sigma2,
  xbar_ch = xbar_ch,
  se_ch = se_ch,
  psi = 0.001,
  z_C = z_C,
  z_T = z_T,
  n_tau = 401,
  ncore0 = 4
)$res

## Reduce n_C until the target power is no longer attained
while (pow_diff_MAP0 >= 0) {
  
  n_C <- n_C - 1
  n_T <- 2 * n_C
  
  pow_diff_MAP0 <- cal_f_map(
    DEL_l = -DEL_u_min,
    DEL_u = DEL_u_min,
    n_T = n_T,
    n_C = n_C,
    alpha = alpha,
    beta = beta,
    theta0 = theta0,
    theta1 = theta1,
    sigma2 = sigma2,
    xbar_ch = xbar_ch,
    se_ch = se_ch,
    psi = 0.001,
    z_C = z_C,
    z_T = z_T,
    n_tau = 401,
    ncore0 = 4
  )$res
}

## Smallest feasible control-arm sample size
n_C_min_MAP <- n_C + 1

## Sequence of feasible control-arm sample sizes
nC_seq_MAP <- n_C_min_MAP:n_C_max_MAP

#---- Find the largest discrepancy bound for each feasible n_C ----

find_Del_max_MAP <- function(n_C, Del_grid) {
  
  n_T <- 2 * n_C
  n_grid <- length(Del_grid)
  
  ## Evaluate the power difference at the smallest value of psi
  pow_diff_0 <- function(i) {
    cal_f_map(
      DEL_l = -Del_grid[i],
      DEL_u = Del_grid[i],
      n_T = n_T,
      n_C = n_C,
      alpha = alpha,
      beta = beta,
      theta0 = theta0,
      theta1 = theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      psi = 0.001,
      z_C = z_C,
      z_T = z_T,
      n_tau = 401,
      ncore0 = 4
    )$res
  }
  
  ## Check whether some psi in the search range yields a feasible design
  is_feasible <- function(i) {
    
    psi_res_map <- find_plausible_psi(
      DEL_l = -Del_grid[i],
      DEL_u = Del_grid[i],
      n_T = n_T,
      n_C = n_C,
      alpha = alpha,
      beta = beta,
      theta0 = theta0,
      theta1 = theta1,
      sigma2 = sigma2,
      xbar_ch = xbar_ch,
      se_ch = se_ch,
      z_C = z_C,
      z_T = z_T,
      psi_grid = seq(
        0.001,
        3 * sigma,
        length.out = 30
      ),
      n_tau = 401,
      ncore0 = 4
    )
    
    is.finite(psi_res_map$psi_use)
  }
  
  ## Return NA if the smallest psi remains feasible at the largest discrepancy
  if (pow_diff_0(n_grid) > 0) {
    return(NA_real_)
  }
  
  ## Locate the first grid value at which the smallest psi becomes infeasible
  left <- 1L
  right <- n_grid
  
  while (right - left > 1L) {
    
    mid <- floor((left + right) / 2)
    
    if (pow_diff_0(mid) > 0) {
      left <- mid
    } else {
      right <- mid
    }
  }
  
  lo <- if (pow_diff_0(left) <= 0) left else right
  
  ## Check feasibility after allowing psi to vary
  feasible <- is_feasible(lo)
  
  if (feasible) {
    
    hi <- length(Del_grid)
    
    ## Return the upper endpoint if it remains feasible
    if (is_feasible(hi)) {
      return(Del_grid[hi])
    }
    
    ## Binary search for the largest feasible discrepancy bound
    while (hi - lo > 1L) {
      
      mid <- floor((lo + hi) / 2L)
      
      if (is_feasible(mid)) {
        lo <- mid
      } else {
        hi <- mid
      }
    }
    
  } else {
    
    lo <- lo - 1L
  }
  
  ## Return the largest feasible discrepancy bound
  Del_grid[lo]
}

## Determine the largest feasible discrepancy bound for each n_C
Del_max_MAP <- find_Del_max_MAP(
  nC_seq_MAP[1],
  Del_grid
)

for (i in 2:length(nC_seq_MAP)) {
  
  Del_max_MAP[i] <- find_Del_max_MAP(
    nC_seq_MAP[i],
    Del_grid[Del_grid >= Del_max_MAP[i - 1]]
  )
}

#### Figure 4: Minimum required control-arm sample size

## Plot the minimum required n_C as a function of Delta_u
plot(
  c(0, Del_max),
  c(nC_seq, nC_seq[length(nC_seq)]),
  type = "s",
  lwd = 2,
  lty = 1,
  xlab = expression(Delta[u]),
  ylab = expression(n[C]),
  ylim = c(6, 22)
)

## Add the APP results
lines(
  c(0, Del_max_APP),
  c(nC_seq_APP, nC_seq_APP[length(nC_seq_APP)]),
  type = "s",
  col = 2,
  lty = 2,
  lwd = 2
)

## Add the MAP results
lines(
  c(0, Del_max_MAP),
  c(nC_seq_MAP, nC_seq_MAP[length(nC_seq_MAP)]),
  type = "s",
  col = 3,
  lty = 3,
  lwd = 2
)

## Add the legend
legend(
  "topleft",
  legend = c("FPP", "APP", "MAP"),
  col = c(1, 2, 3),
  lty = c(1, 2, 3),
  lwd = 2,
  bty = "n"
)
