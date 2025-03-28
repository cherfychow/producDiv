
# ProducDiv
# Fitting Bayesian GLMMs

require(tidytable)
require(brms)
require(tidybayes)
require(performance)
require(ggdist)

# gotta scale
vars_single$hgrav2 <- vars_single$hgravity / 500
vars_time_long$hgrav2 <- vars_time_long$hgravity / 500
vars_time$hgrav2 <- vars_time$hgravity / 500

# MODEL 1: SPAtIAL =================================================================================
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
                   prior = spatial_priors, iter = 5000)
summary(mod_spatial)
model_performance(mod_spatial)

## POStERIORS PREDICtIONS ==========================================================================

# get posterior draws
spatial_draws <- tidy_draws(mod_spatial) %>% 
  select(starts_with(c('b', 'sd', '.')), shape, Intercept) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

spatial_draws <- spatial_draws %>% mutate(var_type = case_when(str_detect(vars, 'log') ~ 'fixed',
                                                   .default = 'randommodel'),
                              prettyvar = case_when(vars == 'b_logS' ~ 'S',
                                                    vars == 'b_loghgrav2P1' ~ 'Human gravity',
                                                    vars == 'sd_DataID__Intercept' ~ '(1|DataID)',
                                                    vars == 'sd_gridcell__Intercept' ~ '(1|gridcell)',
                                                    vars == 'Intercept' ~ 'α',
                                                    vars == 'b_Intercept' ~ 'β0',
                                                    .default = vars))

# model partial predictions
newdata <- tidytable(S = runif(500, min(vars_single$S, na.rm = T), 
                                         max(vars_single$S, na.rm = T)) %>% floor(),
                     hgrav2 = runif(500, min(vars_single$hgrav2), max(vars_single$hgrav2)))

# partial effect of ∆S
spatial_partialpred <- newdata %>% mutate(hgrav2 = mean(vars_single$hgrav2)) %>% 
  epred_draws(mod_spatial, re_formula = NA, newdata = ., value = 'prod_S') %>% ungroup %>% 
  select(!hgrav2)
# partial effect of human gravity
spatial_partialpred <- newdata %>% mutate(S = mean(vars_single$S, na.rm = T)) %>% 
  epred_draws(mod_spatial, re_formula = NA, newdata = ., value = 'prod_hgrav') %>% 
  ungroup %>% select(hgrav2, prod_hgrav) %>% 
  bind_cols(spatial_partialpred, .) %>% ungroup
rm(newdata)

# MODEL 2: tEMPORAL, FULL tIME SERIES ==============================================================
# does proportional diversity change result in proportional productivity change?

t_priors <- default_prior(prod_slope ~ S_slope + log(hgrav2 + 1) + (1|DataID) + (1|gridcell),
                              data = vars_time_long, family = gaussian())
t_priors
t_priors$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative
t_priors$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive
set.seed(2)
mod_t <- brm(data = vars_time_long,
                 formula = prod_slope ~ S_slope + log(hgrav2 + 1) + (1|DataID) + (1|gridcell),
                 family = gaussian(),
                 prior = t_priors, iter = 5000)

summary(mod_t)
model_performance(mod_t)

## POStERIORS PREDICtIONS ==========================================================================
 
# get posterior draws
t_draws <- tidy_draws(mod_t) %>% 
  select(starts_with(c('b', 'sd', '.')), sigma, Intercept) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

t_draws <- t_draws %>% mutate(var_type = case_when(str_detect(vars, 'slope|grav') ~ 'fixed',
                                                               .default = 'randommodel'),
                                          prettyvar = case_when(vars == 'b_S_slope' ~ '∆S',
                                                                vars == 'b_loghgrav2P1' ~ 'Human gravity',
                                                                vars == 'sd_DataID__Intercept' ~ '(1|DataID)',
                                                                vars == 'sd_gridcell__Intercept' ~ '(1|gridcell)',
                                                                vars == 'sigma' ~ 'σ',
                                                                vars == 'Intercept' ~ 'α',
                                                                vars == 'b_Intercept' ~ 'β0',
                                                                .default = vars))

# model partial predictions
newdata <- tidytable(S_slope = runif(500, min(vars_time_long$S_slope), max(vars_time_long$S_slope)),
                           hgrav2 = runif(500, min(vars_time_long$hgrav2), max(vars_time_long$hgrav2)))

# partial effect of ∆S
t_partialpred <- newdata %>% mutate(hgrav2 = mean(vars_time_long$hgrav2)) %>% 
  epred_draws(mod_t, re_formula = NA, newdata = ., value = 'deltaProd_S') %>% ungroup %>% 
  select(!hgrav2)
# partial effect of human gravity
t_partialpred <- newdata %>% mutate(S_slope = mean(vars_time_long$S_slope)) %>% 
  epred_draws(mod_t, re_formula = NA, newdata = ., value = 'deltaProd_hgrav') %>% 
  ungroup %>% select(hgrav2, deltaProd_hgrav) %>% bind_cols(t_partialpred, .)


# MODEL 2: tEMPORAL, FULL tIME SERIES ==============================================================
# does proportional diversity change result in proportional productivity change?

t_priors <- default_prior(logratio1_prod ~ logratio1_S + log(hgrav2 + 1) + (1|DataID) + (1|gridcell),
                          data = vars_time, family = gaussian())
t_priors
t_priors$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative
t_priors$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive
set.seed(2)
mod_t <- brm(data = vars_time,
             formula = logratio1_prod ~ logratio1_S + log(hgrav2 + 1) + (1|DataID) + (1|gridcell),
             family = gaussian(),
             prior = t_priors, iter = 5000)

## POStERIORS PREDICtIONS ==========================================================================

# get posterior draws
t_draws <- tidy_draws(mod_t) %>% 
  select(starts_with(c('b', 'sd', '.')), sigma, Intercept) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

t_draws <- t_draws %>% mutate(var_type = case_when(str_detect(vars, 'log|grav') ~ 'fixed',
                                                   .default = 'randommodel'),
                              prettyvar = case_when(vars == 'b_logratio1_S' ~ '∆log2(S)',
                                                    vars == 'b_loghgrav2P1' ~ 'Human gravity',
                                                    vars == 'sd_DataID__Intercept' ~ '(1|DataID)',
                                                    vars == 'sd_gridcell__Intercept' ~ '(1|gridcell)',
                                                    vars == 'sigma' ~ 'σ',
                                                    vars == 'Intercept' ~ 'α',
                                                    vars == 'b_Intercept' ~ 'β0',
                                                    .default = vars))

# model partial predictions
newdata <- tidytable(logratio1_S = runif(500, min(vars_time$logratio1_S, na.rm = T), 
                                         max(vars_time$logratio1_S, na.rm = T)),
                     hgrav2 = runif(500, min(vars_time$hgrav2), max(vars_time$hgrav2)))

# partial effect of ∆S
t_partialpred <- newdata %>% mutate(hgrav2 = mean(vars_time$hgrav2)) %>% 
  epred_draws(mod_t, re_formula = NA, newdata = ., value = 'deltaProd_S') %>% ungroup %>% 
  select(!hgrav2)
# partial effect of human gravity
t_partialpred <- newdata %>% mutate(logratio1_S = mean(vars_time$logratio1_S, na.rm = T)) %>% 
  epred_draws(mod_t, re_formula = NA, newdata = ., value = 'deltaProd_hgrav') %>% 
  ungroup %>% select(hgrav2, deltaProd_hgrav) %>% 
  bind_cols(t_partialpred, .) %>% ungroup

# MODEL 4: BIOGEO ==================================================================================
# what is the diversity-productivity relationship like across space

# weakly informative priors
geo_priors <- default_prior(prod_kg ~ log(S) + log(hgrav2) + (1 + log(S)|biogeo) + (1|DataID) + (1|gridcell),
                            data = vars_single, 
                            family = Gamma(link = 'log'))
geo_priors$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative
geo_priors$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive


set.seed(1)
mod_geo <- brm(data = vars_single,
                   formula = prod_kg_yr ~ log(S) + hgravity + (1 + log(S)|biogeo) + (1|DataID) + (1|gridcell),
                   family = Gamma(link = 'log'),
                   prior = geo_priors, iter = 5000)


## POStERIORS PREDICtIONS ==========================================================================

# get posterior draws
geo_draws <- tidy_draws(mod_geo) %>% 
  select(starts_with(c('b', 'sd', '.')), shape, Intercept) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

geo_draws <- geo_draws %>% mutate(var_type = case_when(str_detect(vars, 'log') ~ 'fixed',
                                                               .default = 'randommodel'),
                                          prettyvar = case_when(vars == 'b_logS' ~ 'S',
                                                                vars == 'b_loghgrav2P1' ~ 'Human gravity',
                                                                vars == 'sd_DataID__Intercept' ~ '(1|DataID)',
                                                                vars == 'sd_gridcell__Intercept' ~ '(1|gridcell)',
                                                                vars == 'Intercept' ~ 'α',
                                                                vars == 'b_Intercept' ~ 'β0',
                                                                .default = vars))

# model partial predictions
newdata <- tidytable(S = runif(500, min(vars_single$S, na.rm = T), 
                               max(vars_single$S, na.rm = T)) %>% floor(),
                     hgrav2 = runif(500, min(vars_single$hgrav2), max(vars_single$hgrav2)))

# partial effect of ∆S
geo_partialpred <- newdata %>% mutate(hgrav2 = mean(vars_single$hgrav2)) %>% 
  epred_draws(mod_geo, re_formula = NA, newdata = ., value = 'prod_S') %>% ungroup %>% 
  select(!hgrav2)
# partial effect of human gravity
geo_partialpred <- newdata %>% mutate(S = mean(vars_single$S, na.rm = T)) %>% 
  epred_draws(mod_geo, re_formula = NA, newdata = ., value = 'prod_hgrav') %>% 
  ungroup %>% select(hgrav2, prod_hgrav) %>% 
  bind_cols(geo_partialpred, .) %>% ungroup
rm(newdata)