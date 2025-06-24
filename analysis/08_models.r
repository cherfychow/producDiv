
# ProducDiv
# Fitting Bayesian GLMMs

require(tidytable)
require(brms)
require(tidybayes)


# MODEL 1: SPATIAL BIOGEO ==================================================================================
# what is the diversity-productivity relationship like across space

# weakly informative priors
pr_spatial <- default_prior(prod_kg ~ log(S) + log(hgravity + .0000390505) + (1 + log(S)|biogeo) + (1|DataID) + (1|gridcell:DataID) + (1|grav_ID),
                            data = vars_single,
                            family = Gamma(link = 'log'))
pr_spatial$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative
pr_spatial$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive


set.seed(1)
m_spatial <- brm(data = vars_single,
                 formula = prod_kg ~ log(S) + log(hgravity + .0000390505) + (1 + log(S)|biogeo) + (1|DataID) + (1|gridcell:DataID) + (1|grav_ID),
                 family = Gamma(link = 'log'),
                 prior = pr_spatial, iter = 5000)
save(m_spatial, file = 'm_spatial2.RData')
print('Geographical model saved')

## POSTERIORS PREDICTIONS ==========================================================================

# get posterior draws
dr_spatial <- tidy_draws(m_spatial) %>% 
  select(starts_with(c('b', 'sd')), shape, Intercept) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

dr_spatial <- dr_spatial %>% mutate(var_type = case_when(str_detect(vars, 'b\\_log|sd\\_biogeo\\_\\_logS') ~ 'main',
                                                       .default = 'randommodel'),
                                  prettyvar = case_when(vars == 'b_logS' ~ 'S',
                                                        vars == 'b_loghgravityP3.90505eM05' ~ 'Market gravity',
                                                        vars == 'sd_biogeo__logS' ~ '(S|Biogeo)',
                                                        vars == 'sd_biogeo__Intercept' ~ '(1|Biogeo)',
                                                        vars == 'sd_DataID__Intercept' ~ '(1|study)',
                                                        vars == 'sd_gridcell:DataID__Intercept' ~ '(1|location:study)',
                                                        vars == 'sd_gridcell__Intercept' ~ '(1|location)',
                                                        vars == 'sd_grav_ID__Intercept' ~ '(1|gravity cell)',
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
                         hgravity = mean(vars_single$hgravity)) %>% 
      bind_rows(newdata, .)
  } else {
    newdata <- tidytable(biogeo = temp$biogeo[i], 
                         S = runif(100, temp$minS[i], temp$maxS[i]) %>% floor(),
                         hgravity = mean(vars_single$hgravity))
  }
}


# partial effect of ∆S d
pred_biogeo <- newdata %>% mutate(hgravity = mean(vars_single$hgravity)) %>% 
  epred_draws(m_spatial, re_formula = prod_kg ~ log(S) + log(hgravity + .0000390505) + (1 + log(S) | biogeo), 
              newdata = ., value = 'prod_S') %>% ungroup %>% 
  select(!hgravity)

# pred_biogeo <- tidytable(hgravity = runif(n_distinct(vars_single$S), min(vars_single$hgravity, na.rm = T), 
#                                               max(vars_single$hgravity, na.rm = T)) %>% floor(),
#                              S = mean(vars_single$S, na.rm = T)) %>% 
#   epred_draws(m_spatial, re_formula = NA, newdata = ., value = 'prod_hgrav') %>% 
#   ungroup %>% select(hgravity, prod_hgrav) %>% 
#   bind_cols(pred_spatial, .) %>% ungroup

rm(newdata)

# replace spatial preds with this model

# partial effect of ∆S
pred_spatial <- vars_single %>% distinct(S) %>% 
  mutate(hgravity = mean(vars_single$hgravity)) %>% 
  epred_draws(m_spatial, re_formula = NA, newdata = ., value = 'prod_S') %>% 
  ungroup %>% select(!hgravity)

# residuals
spatial_resid <- vars_single %>% select(S, prod_kg)
spatial_resid <- pred_spatial %>% group_by(S) %>% summarise(prod_S = median(prod_S)) %>% 
  left_join(spatial_resid, ., by = 'S') %>% distinct
spatial_resid <- spatial_resid %>% mutate(resid = prod_S - prod_kg)


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

# MODEL 5: QUANTILE BIOGEO ==================================================================================

require(quantreg)
m_quant <- as.list(rep(NA, 3))
quants = c(.01, .5, .99)
for (i in 1:3) {
  m_quant[[i]] <- rq(log(prod_kg) ~ log(S) + log(hgravity + min(grid_sf$hgravity)), method = 'fn',
                        tau = quants[i], data = vars_single)
}
rm(quants)

newdata <- tidytable(S = unique(vars_single$S), hgravity = mean(vars_single$S))

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
