#### Tables A2-A3: Average Type I error and power

## Load the average Type I error and power results
load("Rdatas/Case_study_ave_TIE_power.RData")

#### Table A2: Average operating characteristics over the full range

## Combine the results for all borrowing methods
Table_A2 <- rbind(
  tab_fpp,
  tab_app,
  tab_map,
  tab_rmap_1,
  tab_rmap_2
)

## Print Table A2
print(
  Table_A2,
  row.names = FALSE
)


#### Table A3: Average operating characteristics over the
#### prespecified discrepancy region

## Combine the results for all borrowing methods
Table_A3 <- rbind(
  tab_fpp_ps,
  tab_app_ps,
  tab_map_ps,
  tab_rmap_1_ps,
  tab_rmap_2_ps
)

## Print Table A3
print(
  Table_A3,
  row.names = FALSE
)
