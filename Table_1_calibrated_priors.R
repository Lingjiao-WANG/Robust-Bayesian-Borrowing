#### Table 1: Calibrated analysis priors

## Load the calibrated analysis-prior results
load("Rdatas/Calibrated_analysis_priors.RData")

## Construct Table 1 containing the calibrated borrowing parameters
## and posterior probability thresholds for each analysis prior
Table_1 <- rbind(
  
  ## Fixed power prior
  data.frame(
    Method = "FPP",
    Borrowing_parameter = paste0(
      "a0 = ",
      round(alpha_0, 3)
    ),
    Threshold = paste0(
      "p = ",
      round(p_0_fpp, 3)
    )
  ),
  
  ## Adaptive power prior
  data.frame(
    Method = "APP",
    Borrowing_parameter = paste0(
      "delta0 = ",
      round(delta0, 3)
    ),
    Threshold = paste0(
      "p = ",
      round(p_0_app, 3)
    )
  ),
  
  ## Meta-analytic predictive prior
  data.frame(
    Method = "MAP",
    Borrowing_parameter = paste0(
      "psi = ",
      round(psi_cal_map, 3)
    ),
    Threshold = paste0(
      "p = ",
      round(p_0_map, 3)
    )
  ),
  
  ## Robust MAP prior with fixed w and calibrated psi
  data.frame(
    Method = "rMAP_1",
    Borrowing_parameter = paste0(
      "w = ",
      round(w_rob, 3),
      ", psi = ",
      round(psi_cal_rmap, 3)
    ),
    Threshold = paste0(
      "p = ",
      round(p_0_rmap_1, 3)
    )
  ),
  
  ## Robust MAP prior with fixed psi and calibrated w
  data.frame(
    Method = "rMAP_2",
    Borrowing_parameter = paste0(
      "w = ",
      round(w_cal, 3),
      ", psi = ",
      round(psi_rmap_2, 3)
    ),
    Threshold = paste0(
      "p = ",
      round(p_0_rmap_2, 3)
    )
  )
)

## Display Table 1
print(
  Table_1,
  row.names = FALSE
)
