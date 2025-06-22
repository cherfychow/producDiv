
# ProducDiv
# Fitting Bayesian GLMMs

require(tidytable)
require(brms)
require(tidybayes)
require(ggdist)
require(lqmm)

# MODEL 1: SPATIAL =================================================================================
# what is the diversity-productivity relationship like across space

# weakly informative priors
spatial_priors <- default_prior(prod_kg ~ log(S) + log(hgravity) + (1|DataID) + (1|gridcell),
              data = vars_single, family = Gamma(link = 'log'))
spatial_priors$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative
spatial_priors$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive

set.seed(1)
mod_spatial <- brm(data = vars_single,
                   formula = prod_kg_yr ~ log(S) + hgravity + (1|DataID) + (1|gridcell),
                   family = Gamma(link = 'log'),
                   prior = spatial_priors, iter = 5000)
summary(mod_spatial)
model_performance(mod_spatial)

## POSTERIORS PREDICTIONS ==========================================================================

# get posterior draws
spatial_draws <- tidy_draws(mod_spatial) %>% 
  select(starts_with(c('b', 'sd', '.')), shape, Intercept) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

spatial_draws <- spatial_draws %>% mutate(var_type = case_when(str_detect(vars, 'log') ~ 'fixed',
                                                   .default = 'randommodel'),
                              prettyvar = case_when(vars == 'b_logS' ~ 'S',
                                                    vars == 'b_loghgravityP1' ~ 'Human gravity',
                                                    vars == 'sd_DataID__Intercept' ~ '(1|DataID)',
                                                    vars == 'sd_gridcell__Intercept' ~ '(1|gridcell)',
                                                    vars == 'Intercept' ~ 'α',
                                                    vars == 'b_Intercept' ~ 'β0',
                                                    .default = vars))

# model partial predictions
newdata <- tidytable(S = runif(500, min(vars_single$S, na.rm = T), 
                                         max(vars_single$S, na.rm = T)) %>% floor(),
                     hgravity = runif(500, min(vars_single$hgravity), max(vars_single$hgravity)))

# partial effect of ∆S
spatial_partialpred <- newdata %>% mutate(hgravity = mean(vars_single$hgravity)) %>% 
  epred_draws(mod_spatial, re_formula = NA, newdata = ., value = 'prod_S') %>% ungroup %>% 
  select(!hgravity)
# partial effect of human gravity
spatial_partialpred <- newdata %>% mutate(S = mean(vars_single$S, na.rm = T)) %>% 
  epred_draws(mod_spatial, re_formula = NA, newdata = ., value = 'prod_hgrav') %>% 
  ungroup %>% select(hgravity, prod_hgrav) %>% 
  bind_cols(spatial_partialpred, .) %>% ungroup
rm(newdata)



# MODEL 2: TEMPORAL, ADJACENT YEARS ================================================================
# does proportional diversity change result in proportional productivity change?



## POSTERIORS PREDICTIONS ==========================================================================
 
# get posterior draws
t_draws <- tidy_draws(mod_t) %>% 
  select(starts_with(c('b', 'sd', '.')), sigma) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

t_draws <- t_draws %>% mutate(var_type = case_when(str_detect(vars, 'log|delta') ~ 'fixed',
                                                               .default = 'randommodel'),
                                          prettyvar = case_when(vars == 'b_logratio1_S' ~ '∆log (S)',
                                                                vars == 'b_loghgravityP1' ~ 'Human gravity',
                                                                vars == 'sd_DataID__Intercept' ~ '(1|study)',
                                                                vars == 'sd_gridcell:DataID__Intercept' ~ '(1|location:study)',
                                                                vars == 'sigma' ~ 'σ',
                                                                vars == 'b_Intercept' ~ 'β0',
                                                                .default = vars))

# model partial predictions
newdata <- tidytable(logratio1_S = runif(500, min(vars_time$logratio1_S, na.rm = T), max(vars_time$logratio1_S, na.rm = T)),
                           hgravity = runif(500, min(vars_time$hgravity), max(vars_time$hgravity)))

# partial effect of ∆S
t_partialpred <- newdata %>% mutate(hgravity = mean(vars_time$hgravity)) %>% 
  epred_draws(mod_t, re_formula = NA, newdata = ., value = 'deltaProd_S') %>% ungroup %>% 
  select(!hgravity)
# partial effect of human gravity
t_partialpred <- newdata %>% mutate(logratio1_S = mean(vars_time$logratio1_S)) %>% 
  epred_draws(mod_t, re_formula = NA, newdata = ., value = 'deltaProd_hgrav') %>% 
  ungroup %>% select(hgravity, deltaProd_hgrav) %>% bind_cols(t_partialpred, .)


# MODEL 2: TEMPORAL, FULL TIME SERIES ==============================================================
# does proportional diversity change result in proportional productivity change?



## POSTERIORS PREDICTIONS ==========================================================================

# get posterior draws
dr_time_full <- tidy_draws(m_time_full) %>% 
  select(starts_with(c('b', 'sd', '.')), sigma) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

dr_time_full <- dr_time_full %>% mutate(var_type = case_when(str_detect(vars, 'log|slope') ~ 'fixed',
                                                   .default = 'randommodel'),
                              prettyvar = case_when(vars == 'b_S_slope' ~ '∆log S',
                                                    vars == 'b_loghgravityP3.90505eM05' ~ 'Market gravity',
                                                    vars == 'sd_DataID__Intercept' ~ '(1|study)',
                                                    vars == 'sd_grav_ID__Intercept' ~ '(1|gravity cell)',
                                                    vars == 'sd_gridcell:DataID__Intercept' ~ '(1|location:study)',
                                                    vars == 'sigma' ~ 'σ',
                                                    vars == 'b_Intercept' ~ 'β0',
                                                    .default = vars))

# model partial predictions
newdata <- tidytable(S_slope = runif(500, min(vars_time_long$S_slope, na.rm = T), 
                                         max(vars_time_long$S_slope, na.rm = T)),
                     hgravity = runif(500, min(vars_time_long$hgravity), max(vars_time_long$hgravity)))

# partial effect of ∆S
pred_time_full <- newdata %>% 
  mutate(hgravity = mean(vars_time_long$hgravity)) %>% 
  epred_draws(m_time_full, re_formula = NA, newdata = ., value = 'deltaProd_S') %>% ungroup %>% 
  select(!hgravity)

# partial effect of human gravity
# pred_time_full <- newdata %>% mutate(S_slope = mean(vars_time_long$S_slope, na.rm = T)) %>% 
#   epred_draws(m_time_full, re_formula = NA, newdata = ., value = 'deltaProd_hgrav') %>% 
#   ungroup %>% select(hgravity, deltaProd_hgrav) %>% 
#   bind_cols(pred_time_full, .) %>% ungroup

# MODEL 4: BIOGEO ==================================================================================
# what is the diversity-productivity relationship like across space

# weakly informative priors
geo_priors <- default_prior(prod_kg ~ log(S) + log(hgravity + 1) + (1 + log(S)|biogeo) + (1|DataID) + (1|gridcell),
                            data = vars_single, 
                            family = Gamma(link = 'log'))
geo_priors$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative
geo_priors$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive


set.seed(1)
mod_geo <- brm(data = vars_single,
                   formula = prod_kg_yr ~ log(S) + log(hgravity + 1) + (1 + log(S)|biogeo) + (1|DataID) + (1|gridcell),
                   family = Gamma(link = 'log'),
                   prior = geo_priors, iter = 5000)


## POSTERIORS PREDICTIONS ==========================================================================

# get posterior draws
geo_draws <- tidy_draws(mod_geo) %>% 
  select(starts_with(c('b', 'sd')), shape, Intercept) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

geo_draws <- geo_draws %>% mutate(var_type = case_when(vars %in% c('b_logS', 'b_loghgravityP1', 'sd_biogeo__logS') ~ 'main',
                                                               .default = 'randommodel'),
                                          prettyvar = case_when(vars == 'b_logS' ~ 'S',
                                                                vars == 'b_loghgravityP1' ~ 'Human gravity',
                                                                vars == 'sd_biogeo__logS' ~ '(S|Biogeo)',
                                                                vars == 'sd_biogeo__Intercept' ~ '(1|Biogeo)',
                                                                vars == 'sd_DataID__Intercept' ~ '(1|DataID)',
                                                                vars == 'sd_gridcell__Intercept' ~ '(1|gridcell)',
                                                                vars == 'Intercept' ~ 'α',
                                                                vars == 'b_Intercept' ~ 'β0',
                                                                .default = vars))

# model partial predictions
temp <- vars_single %>% ungroup %>% group_by(biogeo) %>% 
  summarise(minS = min(S), maxS = max(S),
            minG = min(hgravity), maxG = max(hgravity))
for (i in 1:nrow(temp)) {
  if (i > 1) {
    newdata <- tidytable(biogeo = temp$biogeo[i], 
                         S = runif(100, temp$minS[i], temp$maxS[i]) %>% floor(),
                         hgravity = runif(100, temp$minG[i], temp$maxG[i]) %>% floor()) %>% 
      bind_rows(newdata, .)
  } else {
    newdata <- tidytable(biogeo = temp$biogeo[i], 
                         S = runif(100, temp$minS[i], temp$maxS[i]) %>% floor(),
                         hgravity = runif(100, temp$minG[i], temp$maxG[i]))
  }
}


# partial effect of ∆S d
geo_partialpred <- newdata %>% mutate(hgravity = mean(vars_single$hgravity)) %>% 
  epred_draws(mod_geo, re_formula = prod_kg ~ log(S) + log(hgravity + 1) + (1 + log(S) | biogeo), 
              newdata = ., value = 'prod_S') %>% ungroup %>% 
  select(!hgravity)

geo_partialpred <- tidytable(hgravity = runif(n_distinct(vars_single$S), min(vars_single$hgravity, na.rm = T), 
                                                max(vars_single$hgravity, na.rm = T)) %>% floor(),
                                 S = mean(vars_single$S, na.rm = T)) %>% 
  epred_draws(mod_geo, re_formula = NA, newdata = ., value = 'prod_hgrav') %>% 
  ungroup %>% select(hgravity, prod_hgrav) %>% 
  bind_cols(spatial_partialpred, .) %>% ungroup

rm(newdata)

# replace spatial preds with this model

# partial effect of ∆S
spatial_partialpred <- vars_single %>% distinct(S) %>% 
  mutate(hgravity = mean(vars_single$hgravity)) %>% 
  epred_draws(mod_geo, re_formula = NA, newdata = ., value = 'prod_S') %>% 
  ungroup %>% select(!hgravity)
# partial effect of human gravity
spatial_partialpred <- tidytable(hgravity = runif(121, min(vars_single$hgravity, na.rm = T), 
                                           max(vars_single$hgravity, na.rm = T)),
                                 S = mean(vars_single$S, na.rm = T)) %>% 
  epred_draws(mod_geo, re_formula = NA, newdata = ., value = 'prod_hgrav') %>% 
  ungroup %>% select(hgravity, prod_hgrav) %>% 
  bind_cols(spatial_partialpred, .) %>% ungroup

# residuals
spatial_resid <- vars_single %>% select(S, prod_kg)
spatial_resid <- spatial_partialpred %>% group_by(S) %>% summarise(prod_S = median(prod_S)) %>% 
  left_join(spatial_resid, ., by = 'S') %>% distinct
spatial_resid <- spatial_resid %>% mutate(resid = prod_S - prod_kg)

# MODEL 5: QUANTILE BIOGEO ==================================================================================

load('vars_single.RData')
require(lqmm)
vars_single$group <- with(vars_single, paste(gridcell, DataID, sep = '_'))
mod_quant <- lqmm(log(prod_kg) ~ log(S) + log(hgravity + 1), random = ~ 1, group = group, tau = c(0.025, 0.5, 0.975),
                  data = vars_single, control = list(LP_max_iter = 8000))
save(mod_quant, file = 'mod_quantile.RData')
summary(mod_quant)

# model predictions
quant_pred <- vars_single %>% mutate(hgravity = mean(hgravity), `log(S)` = log(S),
                                     `log(hgravity + 1)` = log(hgravity + 1))
quant_pred <- predict(mod_quant, level = 0, newdata = quant_pred) %>% 
  bind_cols(quant_pred, .)
quant_pred <- distinct(quant_pred, `log(S)`, `log(hgravity + 1)`, `0.025`, `0.500`, `0.975`)
# transform to response scale
quant_pred <- mutate(quant_pred, across(everything(), exp))
quant_pred <- quant_pred %>% rename(S = `log(S)`, hgravity = `log(hgravity + 1)`, ) %>% 
  mutate(hgravity = hgravity - 1)
quant_pred <- pivot_longer(quant_pred, starts_with('0'), names_to = 'quantile', values_to = 'fit')

# didn't converge!


require(quantreg)
quant_test <- as.list(rep(NA, 3))
quants = c(.01, .5, .99)
for (i in 1:3) {
  quant_test[[i]] <- rq(log(prod_kg) ~ log(S) + log(hgravity + 1), method = 'fn',
                        tau = quants[i], data = vars_single)
}
rm(quants)

newdata <- tidytable(S = unique(vars_single$S), hgravity = mean(vars_single$S))

testpred <- as.list(rep(NA, 3))
names = c('lwr', 'med', 'upr')
for (i in 1:3) {
  testpred[[i]] <- predict(quant_test[[i]], newdata = newdata, se = 'boot', 
                      type = 'percentile', interval = 'confidence', level = 0.95) %>% as_tidytable
  testpred[[i]]$quantile = names[i]
  testpred[[i]]$S <- newdata$S
}
rm(newdata, names)
testpred <- bind_rows(testpred)
testpred <- testpred %>% mutate(across(c(fit, lower, higher), .fns = exp))
