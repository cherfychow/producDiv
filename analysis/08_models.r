
# ProducDiv
# Fitting Bayesian GLMMs

require(tidytable)
require(brms)
require(tidybayes)
require(emmeans)

# gotta scale
vars_single$hgrav2 <- vars_single$hgravity / 500
vars_time_long$hgrav2 <- vars_time_long$hgravity / 500

# MODEL 1: SPATIAL =================================================================================
# what is the diversity-productivity relationship like across space

# weakly informative priors
spatial_priors <- default_prior(prod_kg ~ log(S) + log(hgrav2) + (1|DataID) + (1|gridcell),
              data = vars_single, family = Gamma(link = 'log'))
spatial_priors$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative
spatial_priors$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive

set.seed(1)
mod_spatial <- brm(data = vars_single,
                   formula = prod_kg_yr ~ log(S) + hgravity + (1|DataID) + (1|gridcell),
                   family = Gamma(link = 'log'),
                   prior = spatial_priors, iter = 3000)



set.seed(24)
mod_spatial2 <- glmer(data = vars_single,
                      formula = prod_kg ~ log(S) + log(hgrav2 + 1) + (1|DataID),
                      family = Gamma(link = 'log'))
# partial S
set.seed(24)
spatial2_pred <- tibble(S = runif(500, 1, max(vars_single$S) %>% floor()), 
                        hgrav2 = runif(500, 0, max(vars_single$hgrav2)) %>% floor())
temp <- spatial2_pred %>% mutate(hgrav2 = mean(vars_single$hgrav2)) %>% 
  predict(mod_spatial2, newdata = ., se.fit = T, type = 'link', re.form = NA)
spatial2_pred <- spatial2_pred %>% mutate(prod_Spar = exp(temp$fit),
                                          prod_Spar_l = exp(temp$fit - 1.96*temp$se.fit),
                                          prod_Spar_u = exp(temp$fit + 1.96*temp$se.fit))
# partial grav
temp <- spatial2_pred %>% mutate(S = mean(vars_single$S)) %>% 
  predict(mod_spatial2, newdata = ., se.fit = T, type = 'link', re.form = NA)
spatial2_pred <- spatial2_pred %>% mutate(prod_gpar = exp(temp$fit),
                                          prod_gpar_l = exp(temp$fit - 1.96*temp$se.fit),
                                          prod_gpar_u = exp(temp$fit + 1.96*temp$se.fit))
rm(temp)
# save(mod_spatial, file = 'mod_spatial.RData')

## POSTERIORS PREDICTIONS ==========================================================================

spatial_draws <- tidy_draws(mod_spatial) %>% 
  select(starts_with(c('b', 'sd', '.')), shape, Intercept) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

spatial_draws <- spatial_draws %>% mutate(var_type = case_when(str_detect(vars, '^b\\_') ~ 'fixed',
                                             .default = 'randommodel'),
                        prettyvar = case_when(vars == 'b_logS' ~ 'logS',
                                              vars == 'b_hgravity' ~ 'Human gravity',
                                              vars == 'sd_DataID__Intercept' ~ '(1|DataID)',
                                              vars == 'sd_gridcell__Intercept' ~ '(1|gridcell)',
                                              .default = vars))

# marginal effect predictions
spatial_Spartial <- tidytable(S = runif(500, min(vars_single$S), max(vars_single$S)),
                              hgravity = mean(vars_single$hgravity),
                              gridcell = vars_single$gridcell[200],
                              DataID = vars_single$DataID[3000]) %>% 
  epred_draws(mod_spatial, re_formula = NULL, newdata = .)

spatial_Spartial <- spatial_Spartial %>% group_by(S) %>% 
  mutate(lwr = quantile(.epred, .025), upr = quantile(.epred, .975)) %>% 
  distinct(S, .keep_all = T)

# figure: posteriors, halfeye
ggplot(spatial_draws %>% filter(var_type == 'randommodel' | str_detect(vars, 'Intercept'))) +
  ggdist::stat_halfeye(aes(x = values, y = prettyvar), 
                       normalize = 'groups',
                       fill = '#a1b6c4', colour = '#22302c') +  style


# MODEL 2: TEMPORAL, FULL TIME SERIES ==============================================================
# does proportional diversity change result in proportional productivity change?

tlong_priors <- default_prior(prod_slope ~ S_slope + hgravity + (1|DataID) + (1|gridcell),
                              data = vars_time_long, family = gaussian())
tlong_priors
tlong_priors$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative
tlong_priors$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive
set.seed(2)
mod_tlong <- brm(data = vars_time_long,
                 formula = prod_slope ~ S_slope + hgravity + (1|DataID) + (1|gridcell),
                 family = gaussian(),
                 prior = tlong_priors, iter = 3000)
# save(mod_tlong, file = 'mod_tlong.RData')
summary(mod_tlong)
rm(tlong_priors)

mod_T2 <- lmer(data = vars_time_long,
                formula = prod_slope ~ S_slope + log(hgrav2 + 1) + (1|DataID) + (1|gridcell))

# partial S
set.seed(24)
T2_pred <- tibble(S_slope = runif(500, min(vars_time_long$S_slope), max(vars_time_long$S_slope)), 
                        hgrav2 = runif(500, 0, max(vars_time_long$hgrav2)))
temp <- T2_pred %>% mutate(hgrav2 = mean(vars_time_long$hgrav2)) %>% 
  predict(mod_T2, newdata = ., se.fit = T, type = 'response', re.form = NA)
T2_pred <- T2_pred %>% mutate(prod_Spar = temp$fit,
                              prod_Spar_l = temp$fit - 1.96*temp$se.fit,
                              prod_Spar_u = temp$fit + 1.96*temp$se.fit)
# partial grav
temp <- T2_pred %>% mutate(S_slope = mean(vars_time_long$S_slope)) %>% 
  predict(mod_T2, newdata = ., se.fit = T, type = 'response', re.form = NA)
T2_pred <- T2_pred %>% mutate(prod_gpar = temp$fit,
                              prod_gpar_l = temp$fit - 1.96*temp$se.fit,
                              prod_gpar_u = temp$fit + 1.96*temp$se.fit)
rm(temp)

## POSTERIORS PREDICTIONS ==========================================================================
# 
# # data frame of posterior draws
# tlong_post <- mod_tlong %>% spread_draws(r_DataID[DataID,], r_gridcell[gridcell,], drop = TRUE)