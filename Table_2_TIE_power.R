#### Table 2: Maximum TIE, maximum and minimum power,
####          and corresponding calibrated power gain and loss

## Load the Type I error and power results
load("Rdatas/Case_study_TIE_power.RData")

## Define the target power
target_power <- 1 - beta   # = 0.83

## Define the range of true current control-arm means
## used for evaluating the maximum and minimum power
theta_lower <- -71
theta_upper <- 151

#====================================================================
# Helper function
#====================================================================

make_table_row <- function(prior, TIE, Power) {
  
  ##----------------------------------------------------------
  ## Maximum TIE over the entire evaluated parameter range
  ##----------------------------------------------------------
  
  max_TIE <- max(
    TIE,
    na.rm = TRUE
  )
  
  
  ##----------------------------------------------------------
  ## Restrict power evaluation to:
  ##
  ##   theta_C in [-71, 151]
  ##
  ## and points where TIE <= alpha
  ##----------------------------------------------------------
  
  ind <- (
    mu_grid >= theta_lower &
      mu_grid <= theta_upper &
      TIE <= alpha
  )
  
  if (!any(ind)) {
    stop(
      paste(
        "No values satisfy TIE <= alpha for",
        prior,
        "over the specified power-evaluation range."
      )
    )
  }
  
  
  ##----------------------------------------------------------
  ## Maximum and minimum power
  ##----------------------------------------------------------
  
  max_power <- max(
    Power[ind],
    na.rm = TRUE
  )
  
  min_power <- min(
    Power[ind],
    na.rm = TRUE
  )
  
  
  ##----------------------------------------------------------
  ## Power gains relative to target power 1 - beta
  ##----------------------------------------------------------
  
  max_power_gain <- max_power - target_power
  min_power_gain <- min_power - target_power
  
  
  ##----------------------------------------------------------
  ## Construct table row
  ##----------------------------------------------------------
  
  data.frame(
    Prior = prior,
    
    max_TIE = sprintf(
      "%.3f",
      max_TIE
    ),
    
    max_power = sprintf(
      "%.3f (%.3f)",
      max_power,
      max_power_gain
    ),
    
    min_power = sprintf(
      "%.3f (%.3f)",
      min_power,
      min_power_gain
    ),
    
    stringsAsFactors = FALSE
  )
}


#====================================================================
# Construct rows for each borrowing prior
#====================================================================

tab_fpp <- make_table_row(
  prior = "FPP",
  TIE = TIE_fpp,
  Power = Power_fpp
)

tab_app <- make_table_row(
  prior = "APP",
  TIE = TIE_app_ps,
  Power = Power_app_ps
)

tab_map <- make_table_row(
  prior = "MAP",
  TIE = TIE_map,
  Power = Power_map
)

tab_rmap_1 <- make_table_row(
  prior = "rMAP_1",
  TIE = TIE_rmap_1,
  Power = Power_rmap_1
)

tab_rmap_2 <- make_table_row(
  prior = "rMAP_2",
  TIE = TIE_rmap_2,
  Power = Power_rmap_2
)


#====================================================================
# Combine all methods
#====================================================================

Table_2 <- rbind(
  tab_fpp,
  tab_app,
  tab_map,
  tab_rmap_1,
  tab_rmap_2
)


#====================================================================
# Rename columns
#====================================================================

names(Table_2) <- c(
  "Prior",
  "Max TIE",
  "Max Power (Gain)",
  "Min Power (Loss)"
)


#====================================================================
# Print Table 2
#====================================================================

print(
  Table_2,
  row.names = FALSE
)