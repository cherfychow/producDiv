

# ProducDiv
# Calculating temporal change variables

require(tidytable)

# SMALL SCALE LOG RATIO FLUXES =====================================================================

vars_rare <- vars_rare %>% ungroup %>% 
  arrange(DataID, seriesID, Year) %>%
  group_by(DataID, seriesID) %>%
  # relative to time series start, inclusive
  mutate(Year_rel = Year - min(Year) + 1)
# also in order for each time series

# simpler calculation first for time series with one sampling event per year
# do it in a separate object just in case
vars_rare <- vars_rare %>% ungroup %>% 
  group_by(DataID, seriesID, Year) %>% 
  # temporary column to pick out these time series
  mutate(sample_n = n_distinct(SamplingEventID))

# need an object to help select time series samples
# productivity for all the time series, S only for some
vars_time <- vars_rare %>% distinct(DataID, seriesID, Year_rel, sample_n)


## ∆Prod ---------------------------------------------------------------------------------------

vars_time$logratio1_prod <- NA
# empties for model predicted 
vars_time$logratio1_prod_l <- NA
vars_time$logratio1_prod_u <- NA

# GLMs to estimate year to year difference
for (i in which(vars_time$Year_rel != 1)) {
  # ignore timepoint 1 because it doesn't have a t-1
  set.seed(24)
  years <- c(vars_time$Year_rel[i], (vars_time$Year_rel[i-1]))
  tempdata <- vars_rare %>% ungroup %>% 
    filter(DataID == vars_time$DataID[i],
           seriesID == vars_time$seriesID[i],
           Year_rel %in% years) %>% 
    select(Year_rel, starts_with('prod')) %>% 
    pivot_longer(.df = ., cols = starts_with('prod'), 
                 names_to = 'dump', values_to = 'prod') %>% 
    select(!dump) %>% mutate(prod = prod + 1)
  # translate by 1 so we still get proportional changes from negative/zero values
  tempP <- glm(data = tempdata, 
               prod ~ Year_rel, family = Gamma(link = 'log'))
  temppred <- predict.glm(tempP, newdata = dplyr::tibble(Year_rel = c(years[2], years[2]+1)), 
                          type = 'response', se.fit = T)
  vars_time$logratio1_prod[i] <- log2(temppred$fit[2]/temppred$fit[1])
  vars_time$logratio1_prod_l[i] <- log2(temppred$fit[2] - (1.96*temppred$se.fit[2])) - log2(temppred$fit[1] - (1.96*temppred$se.fit[1]))
  vars_time$logratio1_prod_u[i] <- log2(temppred$fit[2] + (1.96*temppred$se.fit[2])) - log2(temppred$fit[1] + (1.96*temppred$se.fit[1]))
  print(paste0('Time series ', i, ' complete of 14478'))
} 

rm(tempdata, years, tempP, temppred)

## ∆S ---------------------------------------------------------------------------------------

vars_time$logratio1_S <- NA
# empties for model predicted 
vars_time$logratio1_S_l <- NA
vars_time$logratio1_S_u <- NA

# only observed richness for now because productivity has upper and lowers
# the single sample ones don't need modelling
for (i in which(vars_time$sample_n == 1 & vars_time$Year_rel != 1)) {
  # ignore timepoint 1 because it doesn't have a t-1
    years <- c(vars_time$Year_rel[i], (vars_time$Year_rel[i-1]))
    temp <- vars_rare %>% ungroup %>% 
      filter(DataID == vars_time$DataID[i],
             seriesID == vars_time$seriesID[i],
             Year_rel %in% years) %>% pull(S)
      vars_time$logratio1_S[i] <- log2(temp[2] / temp[1])
}

rm(temp)

for (i in which(is.na(vars_time$logratio1_S) & vars_time$Year_rel != 1)) {
  set.seed(24)
  # ignore timepoint 1 because it doesn't have a t-1
  years <- c(vars_time$Year_rel[i], (vars_time$Year_rel[i-1]))
  tempdata <- vars_rare %>% ungroup %>% 
    filter(DataID == vars_time$DataID[i],
           seriesID == vars_time$seriesID[i],
           Year_rel %in% years) %>% 
    tidytable::select(Year_rel, S)
  # translate by 1 so we still get proportional changes from negative/zero values
  tempS <- glm(data = tempdata, 
               S ~ Year_rel, family = quasipoisson())
  temppred <- predict.glm(tempS, newdata = dplyr::tibble(Year_rel = c(years[2], years[2]+1)), 
                          type = 'response', se.fit = T)
  vars_time$logratio1_S[i] <- log2(temppred$fit[2]/temppred$fit[1])
  vars_time$logratio1_S_l[i] <- log2(temppred$fit[2] - (1.96*temppred$se.fit[2])) - log2(temppred$fit[1] - (1.96*temppred$se.fit[1]))
  vars_time$logratio1_S_u[i] <- log2(temppred$fit[2] + (1.96*temppred$se.fit[2])) - log2(temppred$fit[1] + (1.96*temppred$se.fit[1]))
  print(paste0('Time series ', i, ' complete of 14478'))
}

rm(tempS, tempdata, temppred)
summary(vars_time)


# LARGE SCALE TRAJECTORY ===========================================================================

# whole time series trajectory of productivity and diversity

# index object for time series assemblage level
vars_time_long <- vars_rare %>% distinct(DataID, seriesID)

# empty lists to store model objects
mod_prodt <- as.list(rep(NA, nrow(vars_time_long)))
mod_St <- as.list(rep(NA, nrow(vars_time_long)))
# empty fields to store model coefficients (time slope)
vars_time_long$S_slope <- NA
vars_time_long$S_se <- NA
vars_time_long$prod_slope <- NA
vars_time_long$prod_se <- NA

for (i in which(is.na(vars_time_long$S_se))) {
  set.seed(24)
  tempdata <- vars_rare %>% ungroup %>% 
    filter(DataID == vars_time_long$DataID[i],
           seriesID == vars_time_long$seriesID[i]) %>% 
    select(Year_rel, prod_kg_yr, S)
  
  mod_prodt[[i]] <- glm(data = tempdata, prod_kg_yr ~ Year_rel, family = Gamma(link = 'log'))
  mod_St[[i]] <- glm(data = tempdata, S ~ Year_rel, family = quasipoisson(link = 'log'))
  
  # store coefficients
  vars_time_long$S_slope[i] <- coefficients(mod_St[[i]])[2] # slope
  vars_time_long$S_se[i] <- summary(mod_St[[i]])$coefficients[2,2] # slope std error
  
  vars_time_long$prod_slope[i] <- coefficients(mod_prodt[[i]])[2] # slope
  vars_time_long$prod_se[i] <- summary(mod_prodt[[i]])$coefficients[2,2] # slope std error
  print(paste0('Time series ', i, ' complete of ', nrow(vars_time_long)))
}

summary(vars_time_long)

# GRAVITY COVARIATE ===============================================================================

# add in human gravity as well for temporal models

vars_time <- dt_rareass %>% distinct(seriesID, gridcell) %>% 
  left_join(vars_time, ., by = 'seriesID') %>% 
  left_join(., st_drop_geometry(grid_sf), by = 'gridcell')

vars_time_long <- dt_rareass %>% distinct(seriesID, gridcell) %>% 
  left_join(vars_time_long, ., by = 'seriesID') %>% 
  left_join(., st_drop_geometry(grid_sf), by = 'gridcell')
