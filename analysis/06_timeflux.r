
# ProducDiv
# Analysing changes in diversity and productivity over time

require(dplyr)

# SMALL SCALE LOG RATIO FLUXES =====================================================================

dt_time <- dt_variables %>% group_by(DataID, rareID, Year) %>%
  # calculate means for each time point across sampling events
  summarise(S_mean = mean(S), prod_mean = mean(prod_kg_yr)) %>%
  # sort it so that row order reflects time series order
  arrange(DataID, rareID, Year) %>% ungroup %>%
  group_by(DataID, rareID) %>%
  # time point column to guide the for loop later
  mutate(timepoint = 1:n()) 

dt_time$logratio_S <- NA
dt_time$logratio_prod <- NA

for (i in 1:nrow(dt_time)) {
  # ignore timepoint 1 because it doesn't have a t-1
  if (dt_time$timepoint[i] != 1) {
    dt_time$logratio_S[i] <- log(dt_time$S_mean[i]) - log(dt_time$S_mean[i-1])
    dt_time$logratio_prod[i] <- log(dt_time$prod_mean[i]) - log(dt_time$prod_mean[i-1])
  } else next
}

# LARGE SCALE FLUXES ===============================================================================

# calculate temporal stability of productivity and mean richness
dt_time_large <- dt_variables %>% group_by(DataID, rareID) %>%
  summarise(S_mean_series = mean(S), prod_stability = mean(prod_kg_yr) / sd(prod_kg_yr))
