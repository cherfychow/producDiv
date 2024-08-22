
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
    dt_time$logratio_S[i] <- log2(dt_time$S_mean[i]) - log2(dt_time$S_mean[i-1])
    dt_time$logratio_prod[i] <- log2(dt_time$prod_mean[i]) - log2(dt_time$prod_mean[i-1])
  } else next
}

# increase scale to every 2 years
dt_time$logratio2_S <- NA
dt_time$logratio2_prod <- NA
for (i in 1:nrow(dt_time)) {
  # ignore timepoint 1 because it doesn't have a t-1
  if (dt_time$timepoint[i] != 1 & dt_time$timepoint[i] %% 2 == 1) {
    dt_time$logratio2_S[i] <- log2(dt_time$S_mean[i]) - log2(dt_time$S_mean[i-2])
    dt_time$logratio2_prod[i] <- log2(dt_time$prod_mean[i]) - log2(dt_time$prod_mean[i-2])
  } else next
}

# LARGE SCALE FLUXES ===============================================================================

# calculate temporal stability of productivity and mean richness
dt_time_large <- dt_variables %>% group_by(DataID, rareID) %>%
  summarise(S_mean_series = mean(S), prod_stability = mean(prod_kg_yr) / sd(prod_kg_yr))

m_prodT <- as.list(rep(NA, nrow(dt_time_large)))
m_divT <- as.list(rep(NA, nrow(dt_time_large)))
for (i in 1:nrow(dt_time_large)) {
  # subset data
  m_prodT[[i]] <- lm(prod_kg_yr ~ Year, 
      data = dt_variables %>% filter(rareID == dt_time_large$rareID[i]))
  m_divT[[i]] <- lm(S ~ Year, 
      data = dt_variables %>% filter(rareID == dt_time_large$rareID[i]))
}

dt_time_large$S_slope <- NA
dt_time_large$S_se <- NA
dt_time_large$prod_slope <- NA
dt_time_large$prod_se <- NA

for (i in 1:nrow(dt_time_large)) {
  dt_time_large$S_slope[i] <- m_divT[[i]]$coefficients[2]
  dt_time_large$S_se[i] <- summary(m_divT[[i]])$coef[2,2]
  dt_time_large$prod_slope[i] <- m_prodT[[i]]$coefficients[2]
  dt_time_large$prod_se[i] <- summary(m_prodT[[i]])$coef[2,2]
}
