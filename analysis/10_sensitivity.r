

# ProducDiv
# Sensitivity tests

require(tidytable)
require(brms)
require(rfishprod)
require(ggplot2)
require(patchwork)
require(lme4)
require(dplyr)
require(MuMIn)

load('model_vars.RData')

# ==================================================================================================
# Length binning sensitivity
# ==================================================================================================

dt_length <- dt_rareass %>% ungroup %>%
  # match Kmax
  left_join(., dt_vbgf, 
            by = c('taxaID', 'DataID', 'seriesID', 'SamplingEventID')) %>% 
  filter(is.na(Length_lwr), Length <= MaxSizeTL) %>% sample_n(1000) %>%  
  select(Length, Species, a, b, MaxSizeTL, starts_with('Kmax'))

# get single ones
dt_length <- dt_single %>% ungroup %>%
  # match Kmax
  left_join(., prod_single, 
            by = c('taxaID', 'DataID', 'seriesID', 'SamplingEventID')) %>% 
  filter(is.na(Length_lwr), Length <= MaxSizeTL) %>% sample_n(1000) %>%  
  select(Length, Species, a, b, MaxSizeTL, starts_with('Kmax')) %>% 
  bind_rows(dt_length, .)

dt_length <- dt_length %>% 
  mutate(Length_bin = case_when(Length < 5 ~ 2.5,
                                between(Length, 5, 10) ~ 7.5,
                                Length >= 100 ~ floor(Length/50) * 50 + 25,
                                .default = floor(Length/10) * 10 + 5)) %>%  # artificially bin
  # constrain lengths like we did before
  mutate(Length_bin = case_when(Length_bin >= MaxSizeTL ~ MaxSizeTL - 0.1, .default = Length_bin),
         Length = case_when(Length >= MaxSizeTL ~ MaxSizeTL - 0.1, .default = Length))

# calculate gross productivity
dt_length <- dt_length %>% 
  mutate(bin_prod = somaGain(a, b, Lmeas = Length_bin, Lmax = MaxSizeTL, Kmax = Kmax),
         bin_prod_l = somaGain(a, b, Lmeas = Length_bin, Lmax = MaxSizeTL, Kmax = Kmax_lowq),
         bin_prod_u = somaGain(a, b, Lmeas = Length_bin, Lmax = MaxSizeTL, Kmax = Kmax_uppq),
         act_prod = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax),
         act_prod_l = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax_lowq),
         act_prod_u = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax_uppq))

dt_length$delta_prod <- with(dt_length, abs(bin_prod - act_prod))
dt_length$delta_length <- with(dt_length, abs(Length_bin - Length))
dt_length$delta_length_rel <- with(dt_length, Length/MaxSizeTL)

dt_length <- dt_length %>% mutate(across(contains('prod'), ~.x / 1000)) # to kg/day

# --------------------------------------------------------------------------------------------------
#  Sensitivity models
# --------------------------------------------------------------------------------------------------

# the relationship between actual and binned productivity
m_L1 <- glm(data = dt_length, bin_prod ~ log(act_prod), family = Gamma(link = 'log'))

# what drives binned length productivity differences?
# Species size? Individual size? Individual size relative to species size?
dt_length$rel_Length <- with(dt_length, Length/MaxSizeTL)
tempdt <- dt_length %>% select(delta_prod, MaxSizeTL, Length, rel_Length, delta_length) %>%
  filter_all(all_vars(. > 0))

set.seed(24)
m_L <- glm(data = tempdt, 
  delta_prod ~ log(MaxSizeTL) + log(Length) + rel_Length + log(delta_length),
  family = Gamma(link = 'log'), na.action = na.fail)

# dredge to find predictors
m_L_dredge <- dredge(m_L)
m_L_dredge$model <- rownames(m_L_dredge) %>% as.numeric
rownames(m_L_dredge) <- NULL
m_L_dredge <- m_L_dredge %>% arrange(model)

# best model didn't have length
set.seed(24)
m_L <- glm(data = tempdt, 
  delta_prod ~ log(MaxSizeTL) + rel_Length + log(delta_length),
  family = Gamma(link = 'log'))

# --------------------------------------------------------------------------------------------------
#  Generate predictions
# --------------------------------------------------------------------------------------------------

# actual vs binned prod
pred_L1 <- tidytable(act_prod = seq(min(dt_length$act_prod), max(dt_length$act_prod), length.out = 50))
temp <- predict(m_L1, newdata = pred_L1, se.fit = T, type = 'response')
pred_L1$fit <- temp$fit
pred_L1$se <- temp$se.fit

# drivers of productivity differences
pred_length <- tempdt %>% select(delta_length, rel_Length, MaxSizeTL)

temp <- tempdt %>% mutate(across(-delta_length, mean)) %>% 
  predict(m_L, newdata = ., se.fit = T, type = 'response')
pred_length$fit_deltaL <- temp$fit
pred_length$se_deltaL <- temp$se.fit

temp <- tempdt %>% mutate(across(-rel_Length, mean)) %>% 
  predict(m_L, newdata = ., se.fit = T, type = 'response')
pred_length$fit_relL <- temp$fit
pred_length$se_relL <- temp$se.fit

temp <- tempdt %>% mutate(across(-MaxSizeTL, mean)) %>% 
  predict(m_L, newdata = ., se.fit = T, type = 'response')
pred_length$fit_Lmax <- temp$fit
pred_length$se_Lmax <- temp$se.fit


# Sample completeness ==============================================================================

m_spatial_comp <- glmer(data = vars_single %>% filter(compness > 0.8),
                       prod_kg ~ log(S) + log(Grav_tot + .0000390505) + (1 + log(S)|biogeo) + (1|DataID) + (1|gridcell:DataID) + (1|grav_ID),
                       family = Gamma(link = 'log'))
save(m_spatial_comp, file = 'm_spatial_comp.RData')

pred_spatial_comp <- vars_single %>% filter(compness > 0.8) %>% 
  mutate(Grav_tot = mean(Grav_tot)) %>% select(S, Grav_tot) %>% 
  distinct
temp <- predict(m_spatial_comp, newdata = pred_spatial_comp, re.form = NA, se.fit = T)
pred_spatial_comp$prod_pred <- exp(temp$fit)
pred_spatial_comp$prod_lwr <- exp(temp$fit - 1.96*temp$se.fit)
pred_spatial_comp$prod_upr <- exp(temp$fit + 1.96*temp$se.fit)

# MPA effect on gravity ============================================================================

pr_spatial_mpa <- default_prior(prod_kg ~ log(S) + log(Grav_tot + .0000390505) + 
  level_fishing_protection + (log(Grav_tot + .0000390505) | level_fishing_protection) +
  (1 + log(S)|biogeo) + (1|DataID) + (1|gridcell:DataID) + (1|grav_ID),
                            data = vars_single_mpa,
                            family = Gamma(link = 'log'))
pr_spatial_mpa$prior[3] <- 'normal(-0.1, 3)' # human gravity, slight negative
pr_spatial_mpa$prior[4] <- 'normal(0.1, 3)' # richness effect, slight positive

set.seed(1)
m_spatial_mpa <- brm(data = vars_single_mpa,
                 formula = prod_kg ~ log(S) + log(Grav_tot + .0000390505) + 
                 level_fishing_protection + (log(Grav_tot + .0000390505) | level_fishing_protection) +
                 (1 + log(S)|biogeo) + (1|DataID) + (1|gridcell:DataID) + (1|grav_ID),
                 family = Gamma(link = 'log'),
                 prior = pr_spatial_mpa, iter = 5000)
saveRDS(m_spatial_mpa, file = 'outputs/m_spatial_mpa.rds')
print('Spatial MPA model saved')

# MPA effect on gravity N ============================================================================

# weakly informative priors
pr_spatialN_mpa <- default_prior(prod_kg ~ log(S) + log(N) + log(Grav_tot + .0000390505) + 
  level_fishing_protection + (log(Grav_tot + .0000390505) | level_fishing_protection) +
  (1 + log(S)|biogeo) + (1|DataID) + (1|gridcell:DataID) + (1|grav_ID),
                             data = vars_single_mpa,
                             family = Gamma(link = 'log'))
pr_spatialN_mpa$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative
pr_spatialN_mpa$prior[3] <- 'normal(0.1, 3)' # abundance effect, slight positive
pr_spatialN_mpa$prior[4] <- 'normal(0.1, 3)' # richness effect, slight positive

set.seed(1)
m_spatialN_mpa <- brm(data = vars_single_mpa,
                  formula = prod_kg ~ log(S) + log(N) + log(Grav_tot + .0000390505) + 
                  level_fishing_protection + (log(Grav_tot + .0000390505) | level_fishing_protection) +
                  (1 + log(S)|biogeo) + (1|DataID) + (1|gridcell:DataID) + (1|grav_ID),
                  family = Gamma(link = 'log'),
                  prior = pr_spatialN_mpa, iter = 5000)
saveRDS(m_spatialN_mpa, file = 'm_spatialN_mpa.rds')
print('Spatial with N + MPA model saved')