

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
  print(paste0('Time series ', i, ' complete of 15115'))
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
  print(paste0('Time series ', i, ' complete of 15115'))
}

rm(tempS, tempdata, temppred)
summary(vars_time)


# LARGE SCALE TRAJECTORY ===========================================================================

# whole time series trajectory of productivity and diversity

# index object for time series assemblage level
vars_time_long <- vars_rare %>% ungroup %>% distinct(DataID, seriesID)

# empty lists to store model objects
mod_deltaprod <- as.list(rep(NA, nrow(vars_time_long)))
mod_deltaS <- as.list(rep(NA, nrow(vars_time_long)))
# empty fields to store model coefficients (time slope)
vars_time_long$S_slope <- NA
vars_time_long$S_se <- NA
vars_time_long$prod_slope <- NA
vars_time_long$prod_se <- NA

for (i in 1:nrow(vars_time_long)) {
  set.seed(24)
  tempdata <- vars_rare %>% ungroup %>%
    filter(DataID == vars_time_long$DataID[i],
           seriesID == vars_time_long$seriesID[i]) %>%
    select(Year_rel, prod_kg, S)

  mod_deltaprod[[i]] <- glm(data = tempdata, prod_kg ~ Year_rel, family = Gamma(link = 'log'))
  mod_deltaS[[i]] <- glm(data = tempdata, S ~ Year_rel, family = quasipoisson(link = 'log'))
  
  # store coefficients
  vars_time_long$S_slope[i] <- coefficients(mod_deltaS[[i]])[2] # slope
  vars_time_long$S_se[i] <- summary(mod_deltaS[[i]])$coefficients[2,2] # slope std error
  
  vars_time_long$prod_slope[i] <- coefficients(mod_deltaprod[[i]])[2] # slope
  vars_time_long$prod_se[i] <- summary(mod_deltaprod[[i]])$coefficients[2,2] # slope std error
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

## ∆SST ============================================================================================
dt_sst_t <- dt_rareass %>% select(DataID, seriesID, SamplingEventID, gridcell, Year, Month, date_res) %>% 
  left_join(., dt_sst, by = colnames(dt_sst)[-4])

vars_time_long$sst_slope <- NA
vars_time_long$sst_se <- NA
mod_sst <- as.list(rep(NA, nrow(vars_time_long)))

for (i in 1:nrow(vars_time_long)) {
  
  temp <- dt_sst_t %>%
    filter(DataID == vars_time_long$DataID[i], seriesID == vars_time_long$seriesID[i])

  mod_sst[[i]] <- glm(data = temp, sst ~ Year, family = Gamma(link = 'log'))
  
  # store coefficients
  vars_time_long$sst_slope[i] <- coefficients(mod_sst[[i]])[2] # slope
  vars_time_long$sst_se[i] <- summary(mod_sst[[i]])$coefficients[2,2] # slope std error
  print(paste0('Time series ', i, ' complete of ', nrow(vars_time_long)))
  
}

## ∆N ==============================================================================================
temp <- dt_rareass %>% select(DataID, seriesID, SamplingEventID, Abundance, Year) %>% 
  group_by(DataID, seriesID, SamplingEventID, Year) %>% 
  summarise(N = sum(Abundance))

vars_time_long$N_slope <- NA
vars_time_long$N_se <- NA
m_Nt <- as.list(rep(NA, nrow(vars_time_long)))

for (i in 1:nrow(vars_time_long)) {
  
  m_Nt[[i]] <- temp %>%
    filter(DataID == vars_time_long$DataID[i], seriesID == vars_time_long$seriesID[i]) %>% 
  glm(data = ., N ~ Year, family = Gamma(link = 'log'))
  
  # store coefficients
  vars_time_long$N_slope[i] <- coefficients(m_Nt[[i]])[2] # slope
  vars_time_long$N_se[i] <- summary(m_Nt[[i]])$coefficients[2,2] # slope std error
  print(paste0('Time series ', i, ' complete of ', nrow(vars_time_long)))
  
}


# CATEGORISE TEMPORAL TRENDS =======================================================================

vars_time <- vars_time %>% 
  mutate(delta_certain = case_when(.default = 0,
                                   is.na(logratio1_S) ~ NA,
                                   logratio1_S < 0 & logratio1_S_u < 0 ~ 1,
                                   logratio1_prod < 0 & logratio1_prod_u < 0  ~ 1,
                                   logratio1_S > 0 & logratio1_S_l > 0 ~ 1,
                                   logratio1_prod > 0 & logratio1_prod_l > 0  ~ 1
  ),
  Sdir = case_when(.default = NA, logratio1_S > 0 ~ 'inc',
                   logratio1_S < 0 ~ 'dec',
                   logratio1_S == 0 ~ 'zero'),
  Pdir = case_when(.default = NA,
                   logratio1_prod > 0 ~ 'inc',
                   logratio1_prod < 0 ~ 'dec',
                   logratio1_prod == 0 ~ 'zero'))

vars_time_long <- vars_time_long %>% 
  mutate(delta_certain = case_when(.default = 0,
                                   is.na(S_slope) ~ NA,
                                   S_slope > 0 & S_slope - 1.96*S_se > 0 ~ 1,
                                   S_slope < 0 & S_slope + 1.96*S_se < 0 ~ 1,
                                   prod_slope > 0 & prod_slope - 1.96*prod_se > 0 ~ 1,
                                   prod_slope < 0 & prod_slope + 1.96*prod_se < 0 ~ 1
  ),
  Sdir = case_when(.default = NA, S_slope > 0 ~ 'inc',
                   S_slope < 0 ~ 'dec',
                   S_slope == 0 ~ 'zero'),
  Pdir = case_when(.default = NA,
                   prod_slope > 0 ~ 'inc',
                   prod_slope < 0 ~ 'dec',
                   prod_slope == 0 ~ 'zero'))

dt_temp <- vars_time_long %>% 
  group_by(Sdir, delta_certain) %>% 
  summarise(n = n()) %>% mutate(group = 'full', metric = 'S') %>% 
  rename(direction = Sdir)

dt_temp <- vars_time_long %>% 
  group_by(Pdir, delta_certain) %>% 
  summarise(n = n()) %>% mutate(group = 'full', metric = 'Prod') %>% 
  rename(direction = Pdir) %>% 
  bind_rows(dt_temp, .)

dt_temp <- vars_time %>% 
  filter(!is.na(logratio1_prod)) %>% 
  group_by(Sdir, delta_certain) %>% 
  summarise(n = n()) %>% mutate(group = 'adjacent', metric = 'S') %>% 
  rename(direction = Sdir) %>% 
  bind_rows(dt_temp, .)

dt_temp <- vars_time %>% 
  filter(!is.na(logratio1_prod)) %>% 
  group_by(Pdir, delta_certain) %>% 
  summarise(n = n()) %>% mutate(group = 'adjacent', metric = 'Prod') %>% 
  rename(direction = Pdir) %>% 
  bind_rows(dt_temp, .)

dt_temp$delta_certain <- as.factor(dt_temp$delta_certain)

# CATEGORISE TREND PAIRS (parallel/opposite) =======================================================

vars_time <- vars_time %>%
  mutate(quad = case_when(.default = 'opp',
                          logratio1_S > 0 & logratio1_prod > 0 ~ 'par',
                          logratio1_S <= 0 & logratio1_prod <= 0 ~ 'par')) %>%
  mutate(quad = case_when(.default = quad,
                          logratio1_S > 0 & logratio1_prod > 0 & logratio1_S_l > 0 & logratio1_prod_l > 0 ~ 'par_certain',
                          logratio1_S <= 0 & logratio1_S_u <= 0 & logratio1_prod <= 0 & logratio1_prod_u <= 0  ~ 'par_certain',
                          logratio1_S > 0 & logratio1_S_l > 0 &  logratio1_prod <= 0 & logratio1_prod_u <= 0 ~ 'opp_certain',
                          logratio1_S <= 0 & logratio1_S_u <= 0 & logratio1_prod > 0 & logratio1_prod_l > 0 ~ 'opp_certain'))

vars_time_long <- vars_time_long %>% 
  mutate(quad = case_when(.default = 'opp',
                          S_slope > 0 & prod_slope > 0 ~ 'par',
                          S_slope <= 0 & prod_slope <= 0 ~ 'par')) %>% 
  mutate(quad = case_when(.default = quad,
                          S_slope > 0 & prod_slope > 0 & S_slope - 1.96*S_se > 0 & prod_slope - 1.96*prod_se > 0 ~ 'par_certain',
                          S_slope <= 0 & S_slope + 1.96*S_se <= 0 & prod_slope <= 0 & prod_slope + 1.96*prod_se <= 0  ~ 'par_certain',
                          S_slope > 0 & S_slope - 1.96*S_se > 0 &  prod_slope <= 0 & prod_slope + 1.96*prod_se <= 0 ~ 'opp_certain',
                          S_slope <= 0 & S_slope + 1.96*S_se <= 0 & prod_slope > 0 & prod_slope - 1.96*prod_se > 0 ~ 'opp_certain'))

# proportions in positive-negative for long and short term fluxes

dt_tchange <- vars_time_long %>% group_by(quad) %>% 
  summarise(n = n()) %>% mutate(group = 'full')
dt_tchange <- vars_time %>% group_by(quad) %>% 
  summarise(n = n()) %>% mutate(group = 'adjacent') %>% 
  bind_rows(dt_tchange, .)

dt_tchange <- dt_tchange %>% group_by(group) %>% mutate(total = sum(n), propn = n/total)


tempfirst <- vars_rare %>% filter(Year_rel == 1) %>% ungroup %>% group_by(DataID, seriesID) %>% 
  mutate(S1 = mean(S), P1 = mean(prod_kg)) %>% select(DataID, seriesID, S1, P1)
vars_time2 <- vars_rare %>% ungroup %>% group_by(DataID, seriesID) %>% 
  summarise(meanS = mean(S), meanP = mean(prod_kg)) %>% 
  left_join(., tempfirst, by = c('DataID', 'seriesID'))
vars_time2 <- vars_time_long %>% distinct(DataID, seriesID, quad, delta_certain) %>% 
  left_join(vars_time2, ., by = c('DataID', 'seriesID')) %>% mutate(quad = str_remove(quad, '\\_certain'))
rm(tempfirst)

