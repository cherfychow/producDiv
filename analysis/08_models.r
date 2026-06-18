
# ProducDiv
# Fitting Bayesian GLMMs

require(tidytable)
require(brms)
require(tidybayes)
require(stringr)

# MODEL 1: SPATIAL ==================================================================================
# what is the diversity-productivity relationship like across space

# weakly informative priors
pr_spatial <- default_prior(prod_kg ~ log(S) + log(Grav_tot + .0000390505) + (1 + log(S)|biogeo) + (1|DataID) + (1|gridcell:DataID) + (1|grav_ID),
                            data = vars_single,
                            family = Gamma(link = 'log'))
pr_spatial$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative
pr_spatial$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive


set.seed(1)
m_spatial <- brm(data = vars_single,
                 formula = prod_kg ~ log(S) + log(Grav_tot + .0000390505) + (1 + log(S)|biogeo) + (1|DataID) + (1|gridcell:DataID) + (1|grav_ID),
                 family = Gamma(link = 'log'),
                 prior = pr_spatial, iter = 5000)
# save(m_spatial, file = 'm_spatial2.RData')
print('Geographical model saved')

## POSTERIORS PREDICTIONS ==========================================================================

# get posterior draws
dr_spatial <- tidy_draws(m_spatial) %>% 
  select(starts_with(c('b', 'sd')), shape, Intercept) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')
unique(dr_spatial$vars)

dr_spatial <- dr_spatial %>% 
  mutate(var_type = case_when(str_detect(vars, 'b\\_log|sd\\_biogeo\\_\\_logS') ~ 'main', .default = 'randommodel'),
        prettyvar = case_when(vars == 'b_logS' ~ 'log(S)',
                              vars == 'b_logGrav_totP3.90505eM05' ~ 'log(total gravity)',
                              vars == 'sd_biogeo__logS' ~ '(S|biogeo)',
                              vars == 'sd_biogeo__Intercept' ~ '(1|biogeo)',
                              vars == 'sd_DataID__Intercept' ~ '(1|study)',
                              vars == 'sd_gridcell:DataID__Intercept' ~ '(1|location:study)',
                              vars == 'sd_gridcell__Intercept' ~ '(1|location)',
                              vars == 'sd_grav_ID__Intercept' ~ '(1|gravity cell)',
                              vars == 'b_Intercept' ~ 'β0',
                              .default = vars)) %>%
  filter(!vars == 'Intercept') # don't want to show everything that model generates
                                      

# get posteriors for each biogeographic region (level)

dr_biogeo <- ranef(m_spatial, summary = F, groups = 'biogeo')
dr_geo_b0 <- dr_biogeo[[1]][,,1]
dr_geo_S <- dr_biogeo[[1]][,,2]
rm(dr_biogeo)

# these only give the sigma for the random effects

dr_geo_b0 <- as_tidytable(dr_geo_b0) %>% 
  tidyr::pivot_longer(cols = everything(), names_to = 'biogeo', values_to = 'value') %>%
  mutate(value = -3.52 + value) # intercepts for each level

dr_geo_S <- as_tidytable(dr_geo_S) %>% 
  tidyr::pivot_longer(cols = everything(), names_to = 'biogeo', values_to = 'value') %>%
  mutate(value = 1.46 + value)

dr_geo_b0 %>% group_by(biogeo) %>% summarise(median = median(value))
dr_geo_S %>% group_by(biogeo) %>% summarise(median = median(value))

# model partial predictions
temp <- vars_single %>% ungroup %>% group_by(biogeo) %>% 
  summarise(minS = min(S), maxS = max(S),
            minG = min(Grav_tot), maxG = max(Grav_tot))
for (i in 1:nrow(temp)) {
  if (i > 1) {
    newdata <- tidytable(biogeo = temp$biogeo[i], 
                         S = runif(100, temp$minS[i], temp$maxS[i]) %>% floor(),
                         Grav_tot = mean(vars_single$Grav_tot)) %>% 
      bind_rows(newdata, .)
  } else {
    newdata <- tidytable(biogeo = temp$biogeo[i], 
                         S = runif(100, temp$minS[i], temp$maxS[i]) %>% floor(),
                         Grav_tot = mean(vars_single$Grav_tot))
  }
}

# partial effect of ∆S d
pred_biogeo <- newdata %>% mutate(Grav_tot = mean(vars_single$Grav_tot)) %>% 
  epred_draws(m_spatial, re_formula = prod_kg ~ log(S) + log(Grav_tot + .0000390505) + (1 + log(S) | biogeo), 
              newdata = ., value = 'prod_S') %>% ungroup %>% 
  select(!Grav_tot)

# pred_biogeo <- tidytable(Grav_tot = runif(n_distinct(vars_single$S), min(vars_single$Grav_tot, na.rm = T), 
#                                               max(vars_single$Grav_tot, na.rm = T)) %>% floor(),
#                              S = mean(vars_single$S, na.rm = T)) %>% 
#   epred_draws(m_spatial, re_formula = NA, newdata = ., value = 'prod_hgrav') %>% 
#   ungroup %>% select(Grav_tot, prod_hgrav) %>% 
#   bind_cols(pred_spatial, .) %>% ungroup

rm(newdata)

# partial effect of ∆S
pred_spatial <- vars_single %>% distinct(S) %>% 
  mutate(Grav_tot = mean(vars_single$Grav_tot)) %>% 
  epred_draws(m_spatial, re_formula = NA, newdata = ., value = 'prod_S') %>% 
  ungroup %>% select(!Grav_tot)

# residuals
spatial_resid <- vars_single %>% select(S, prod_kg)
spatial_resid <- pred_spatial %>% group_by(S) %>% summarise(prod_S = median(prod_S)) %>% 
  left_join(spatial_resid, ., by = 'S') %>% distinct
spatial_resid <- spatial_resid %>% mutate(resid = prod_kg - prod_S)

# SPATIAL ALTERNATE ================================================================================
# with abundance

# weakly informative priors
pr_spatialN <- default_prior(prod_kg ~ log(S) + log(Grav_tot + .0000390505) + log(N) + (1 + log(S)|biogeo) + (1|DataID) + (1|gridcell:DataID) + (1|grav_ID),
                             data = vars_single,
                             family = Gamma(link = 'log'))
pr_spatialN$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative
pr_spatialN$prior[3] <- 'normal(0.1, 3)' # abundance effect, slight positive
pr_spatialN$prior[4] <- 'normal(0.1, 3)' # richness effect, slight positive

set.seed(1)
m_spatialN <- brm(data = vars_single,
                  formula = prod_kg ~ log(S) + log(Grav_tot + .0000390505) + log(N) + (1 + log(S)|biogeo) + (1|DataID) + (1|gridcell:DataID) + (1|grav_ID),
                  family = Gamma(link = 'log'),
                  prior = pr_spatialN, iter = 5000)
print('Spatial with N model saved')



## POSTERIORS PREDICTIONS ==========================================================================

# get posterior draws
dr_spatialN <- tidy_draws(m_spatialN) %>% 
  select(starts_with(c('b', 'sd')), shape, Intercept) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')
unique(dr_spatialN$vars)

dr_spatialN <- dr_spatialN %>% 
  mutate(var_type = case_when(str_detect(vars, 'b\\_log|sd\\_biogeo\\_\\_logS') ~ 'main',
                              .default = 'randommodel'),
        prettyvar = case_when(vars == 'b_logS' ~ 'log(S)',
                              vars == 'b_logN' ~ 'log(N)',
                              vars == 'b_logGrav_totP3.90505eM05' ~ 'log(total gravity)',
                              vars == 'sd_biogeo__logS' ~ '(S|biogeo)',
                              vars == 'sd_biogeo__Intercept' ~ '(1|biogeo)',
                              vars == 'sd_DataID__Intercept' ~ '(1|study)',
                              vars == 'sd_gridcell:DataID__Intercept' ~ '(1|location:study)',
                              vars == 'sd_gridcell__Intercept' ~ '(1|location)',
                              vars == 'sd_grav_ID__Intercept' ~ '(1|gravity cell)',
                              vars == 'b_Intercept' ~ 'β0',
                              .default = vars)) %>%
  filter(!vars == 'Intercept')

pred_spatialN <- tidytable(S = seq(from = 1, to = 148, length.out = 100), 
  Grav_tot = mean(vars_single$Grav_tot), N = mean(vars_single$N)) %>% 
  epred_draws(m_spatialN, re_formula = NA, newdata = ., value = 'prod_S') %>% 
  ungroup %>% select(!c(Grav_tot, N))
temp <- tidytable(N = seq(1, to = max(vars_single$N), length.out = 100), 
  Grav_tot = mean(vars_single$Grav_tot), S = mean(vars_single$S)) %>% 
  epred_draws(m_spatialN, re_formula = NA, newdata = ., value = 'prod_N') %>% select(prod_N, N)
pred_spatialN <- bind_cols(pred_spatialN, temp)
rm(temp)

# MODEL 2: TEMPORAL ================================================================================
# does proportional diversity change result in proportional productivity change?

pr_time_full <- default_prior(prod_slope ~ S_slope + log(Grav_tot + .0000390505) + (1|DataID) + (1|gridcell:DataID) + (1|grav_ID),
                              data = vars_time_long, family = gaussian())
pr_time_full
pr_time_full$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative
pr_time_full$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive
set.seed(2)
m_time_full <- brm(data = vars_time_long,
                   formula = prod_slope ~ S_slope + log(Grav_tot + .0000390505) + (1|DataID) + (1|gridcell:DataID) + (1|grav_ID),
                   family = gaussian(),
                   prior = pr_time_full, iter = 8000)

## POSTERIORS PREDICTIONS ==========================================================================

# get posterior draws
dr_time_full <- tidy_draws(m_time_full) %>% 
  select(starts_with(c('b', 'sd', '.')), sigma) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

dr_time_full <- dr_time_full %>% 
  mutate(var_type = case_when(str_detect(vars, 'log|slope') ~ 'fixed',
                              .default = 'randommodel'),
        prettyvar = case_when(vars == 'b_S_slope' ~ '∆log S',
                              vars == 'b_logGrav_totP3.90505eM05' ~ 'log(total gravity)',
                              vars == 'sd_DataID__Intercept' ~ '(1|study)',
                              vars == 'sd_grav_ID__Intercept' ~ '(1|gravity cell)',
                              vars == 'sd_gridcell:DataID__Intercept' ~ '(1|location:study)',
                              vars == 'sigma' ~ 'σ',
                              vars == 'b_Intercept' ~ 'β0',
                              .default = vars))

# model partial predictions
# newdata <- tidytable(S_slope = runif(500, min(vars_time_long$S_slope, na.rm = T), 
#                                      max(vars_time_long$S_slope, na.rm = T)),
#                      Grav_tot = mean(vars_time_long$Grav_tot))

# partial effect of ∆S
# pred_time_full <- newdata %>% 
#   epred_draws(m_time_full, re_formula = NA, newdata = ., value = 'deltaProd_S') %>% ungroup %>% 
#   select(!Grav_tot)

# partial effect of human gravity
# pred_time_full <- newdata %>% mutate(S_slope = mean(vars_time_long$S_slope, na.rm = T)) %>% 
#   epred_draws(m_time_full, re_formula = NA, newdata = ., value = 'deltaProd_hgrav') %>% 
#   ungroup %>% select(Grav_tot, deltaProd_hgrav) %>% 
#   bind_cols(pred_time_full, .) %>% ungroup

# Temporal model: with N + SST ===========================================================

pr_timeN <- default_prior(prod_slope ~ S_slope + sst_slope + N_slope + log(Grav_tot + .0000390505) + (1|DataID) + (1|gridcell:DataID) + (1|grav_ID),
                             data = vars_time_long, family = gaussian())
pr_timeN
pr_timeN$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative
pr_timeN$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive
set.seed(2)
m_timeN <- brm(data = vars_time_long,
                  formula = prod_slope ~ S_slope + sst_slope + N_slope + log(Grav_tot + .0000390505) + (1|DataID) + (1|gridcell:DataID) + (1|grav_ID),
                  family = gaussian(),
                  prior = pr_timeN, iter = 8000)
save(m_timeN, file = 'outputs/m_timeN.rda')

# posteriors

dr_timeN <- tidy_draws(m_timeN) %>% 
  select(starts_with(c('b', 'sd', '.')), sigma) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')
unique(dr_timeN$vars)
dr_timeN <- dr_timeN %>% 
  mutate(var_type = case_when(str_detect(vars, 'log|slope') ~ 'fixed',
                              .default = 'randommodel'),
        prettyvar = case_when(vars == 'b_S_slope' ~ '∆log(S)',
                              vars == 'b_sst_slope' ~ '∆SST',
                              vars == 'b_N_slope' ~ '∆log(N)',
                              vars == 'b_logGrav_totP3.90505eM05' ~ 'log(total gravity)',
                              vars == 'sd_DataID__Intercept' ~ '(1|study)',
                              vars == 'sd_gridcell:DataID__Intercept' ~ '(1|location:study)',
                              vars == 'sd_grav_ID__Intercept' ~ '(1|gravity cell)',
                              vars == 'sigma' ~ 'σ',
                              vars == 'b_Intercept' ~ 'β0',
                              .default = vars))

# model partial predictions
newdata <- tidytable(S_slope = seq(min(vars_time_long$S_slope, na.rm = T), 
                                     max(vars_time_long$S_slope, na.rm = T), length.out = 50),
                    N_slope = seq(min(vars_time_long$N_slope, na.rm = T), 
                                     max(vars_time_long$N_slope, na.rm = T), length.out = 50),
                     Grav_tot = mean(vars_time_long$Grav_tot),
                     sst_slope = mean(vars_time_long$sst_slope))

# partial effect of ∆S
pred_time_full <- newdata %>% mutate(N_slope = mean(N_slope)) %>%
  epred_draws(m_timeN, re_formula = NA, ndraws = 2000, newdata = ., value = 'prod_S') %>%
  ungroup() %>% select(S_slope, prod_S) %>%
  group_by(S_slope) %>% 
  summarise(
    prod_S = quantile(prod_S, 0.5),
    prod_S_l = quantile(prod_S, 0.025),
    prod_S_u = quantile(prod_S, 0.975)
    )

temp <- newdata %>% mutate(S_slope = mean(S_slope)) %>%
  epred_draws(m_timeN, re_formula = NA, ndraws = 2000, newdata = ., value = 'prod_N') %>%
  ungroup() %>% select(N_slope, prod_N) %>%
  group_by(N_slope) %>% 
  summarise(
    prod_N = quantile(prod_N, 0.5),
    prod_N_l = quantile(prod_N, 0.025),
    prod_N_u = quantile(prod_N, 0.975)
    )

pred_time_full <- bind_cols(pred_time_full, temp)

# Alternative temporal model, fluxes ===============================================================
# does proportional diversity change result in proportional productivity change?

pr_time_flux <- default_prior(logratio1_prod ~ logratio1_S + log(Grav_tot + .0000390505) + (1|DataID) + (1|gridcell:DataID) + (1|grav_ID),
                              data = vars_time, family = gaussian())
pr_time_flux
pr_time_flux$prior[4] <- 'normal(-0.1, 3)' # human gravity, slight negative
pr_time_flux$prior[5] <- 'normal(0.1, 3)' # richness effect, slight positive

set.seed(2)
m_time_flux <- brm(data = vars_time,
                   formula = logratio1_prod ~ logratio1_S + log(Grav_tot + .0000390505) + (1|DataID) + (1|gridcell:DataID) + (1|grav_ID),
                   family = gaussian(),
                   prior = pr_time_flux, iter = 5000)

## POSTERIORS PREDICTIONS ==========================================================================
 
# get posterior draws
dr_time_flux <- tidy_draws(m_time_flux) %>% 
  select(starts_with(c('b', 'sd', '.')), sigma) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')
unique(dr_time_flux$vars)

dr_time_flux <- dr_time_flux %>% filter(!vars == 'Intercept') %>%
  mutate(var_type = case_when(str_detect(vars, 'log|delta') ~ 'fixed', .default = 'randommodel'),
        prettyvar = case_when(vars == 'b_logratio1_S' ~ '∆log (S)',
                              vars == 'b_logGrav_totP3.90505eM05' ~ 'log(total gravity)',
                              vars == 'sd_DataID__Intercept' ~ '(1|study)',
                              vars == 'sd_gridcell:DataID__Intercept' ~ '(1|location:study)',
                              vars == 'sd_grav_ID__Intercept' ~ '(1|gravity cell)',
                              vars == 'sigma' ~ 'σ',
                              vars == 'b_Intercept' ~ 'β0',
                              .default = vars))

# model partial predictions
newdata <- tidytable(logratio1_S = runif(500, min(vars_time$logratio1_S, na.rm = T), max(vars_time$logratio1_S, na.rm = T)),
                           Grav_tot = mean(vars_time$Grav_tot))

# partial effect of ∆S
pred_time_flux <- epred_draws(m_time_flux, re_formula = NA, newdata = newdata, value = 'deltaProd_S') %>% ungroup %>% 
  select(!Grav_tot)


# MODEL 5: QUANTILE SPATIAL ==================================================================================

require(quantreg)
m_quant <- as.list(rep(NA, 3))
quants = c(.01, .5, .99)
for (i in 1:3) {
  m_quant[[i]] <- rq(formula = log(prod_kg) ~ log(S) + log(Grav_tot + .0000390505), method = "fn",
                        tau = quants[i], data = vars_single)
}
rm(quants)

newdata <- tidytable(S = unique(vars_single$S), Grav_tot = mean(vars_single$S))

pred_quant <- as.list(rep(NA, 3))
names = c('lwr', 'med', 'upr')
for (i in 1:3) {
  pred_quant[[i]] <- predict(m_quant[[i]], newdata = newdata, se = 'boot', 
                      type = 'percentile', interval = 'confidence', level = 0.95) %>% as_tidytable
  pred_quant[[i]]$quantile = names[i]
  pred_quant[[i]]$S <- newdata$S
}
rm(newdata, names)
pred_quant <- bind_rows(pred_quant)
pred_quant <- pred_quant %>% mutate(across(c(fit, lower, higher), .fns = exp))
