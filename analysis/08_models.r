
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



# MODEL 2: TEMPORAL, ADJACENT YEARS ================================================================
# does proportional diversity change result in proportional productivity change?



## POSTERIORS PREDICTIONS ==========================================================================
 
# get posterior draws
t_draws <- tidy_draws(mod_t) %>% 
  select(starts_with(c('b', 'sd', '.')), sigma, Intercept) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

t_draws <- t_draws %>% mutate(var_type = case_when(str_detect(vars, 'log|delta') ~ 'fixed',
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
newdata <- tidytable(logratio1_S = runif(500, min(vars_time$logratio1_S, na.rm = T), max(vars_time$logratio1_S, na.rm = T)),
                           hgrav2 = runif(500, min(vars_time_long$hgrav2), max(vars_time_long$hgrav2)))

# partial effect of ∆S
t_partialpred <- newdata %>% mutate(hgrav2 = mean(vars_time$hgrav2)) %>% 
  epred_draws(mod_t, re_formula = NA, newdata = ., value = 'deltaProd_S') %>% ungroup %>% 
  select(!hgrav2)
# partial effect of human gravity
t_partialpred <- newdata %>% mutate(logratio1_S = mean(vars_time$logratio1_S)) %>% 
  epred_draws(mod_t, re_formula = NA, newdata = ., value = 'deltaProd_hgrav') %>% 
  ungroup %>% select(hgrav2, deltaProd_hgrav) %>% bind_cols(t_partialpred, .)


# MODEL 2: TEMPORAL, FULL TIME SERIES ==============================================================
# does proportional diversity change result in proportional productivity change?



## POSTERIORS PREDICTIONS ==========================================================================

# get posterior draws
T_draws <- tidy_draws(mod_T) %>% 
  select(starts_with(c('b', 'sd', '.')), sigma, Intercept) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

T_draws <- T_draws %>% mutate(var_type = case_when(str_detect(vars, 'log|grav|slope') ~ 'fixed',
                                                   .default = 'randommodel'),
                              prettyvar = case_when(vars == 'b_S_slope' ~ '∆S',
                                                    vars == 'b_loghgrav2P1' ~ 'Human gravity',
                                                    vars == 'b_N_slope' ~ '∆N',
                                                    vars == 'b_sst_slope' ~ '∆SST',
                                                    vars == 'sd_DataID__Intercept' ~ '(1|DataID)',
                                                    vars == 'sd_gridcell__Intercept' ~ '(1|gridcell)',
                                                    vars == 'sigma' ~ 'σ',
                                                    vars == 'Intercept' ~ 'α',
                                                    vars == 'b_Intercept' ~ 'β0',
                                                    .default = vars))

# model partial predictions
newdata <- tidytable(S_slope = runif(500, min(vars_time_long$S_slope, na.rm = T), 
                                         max(vars_time_long$S_slope, na.rm = T)),
                     hgrav2 = runif(500, min(vars_time_long$hgrav2), max(vars_time$hgrav2)))

# partial effect of ∆S
T_partialpred <- newdata %>% 
  mutate(hgrav2 = mean(vars_time_long$hgrav2)) %>% 
  epred_draws(mod_T, re_formula = NA, newdata = ., value = 'deltaProd_S') %>% ungroup %>% 
  select(!hgrav2)
# partial effect of human gravity
T_partialpred <- newdata %>% mutate(S_slope = mean(vars_time_long$S_slope, na.rm = T)) %>% 
  epred_draws(mod_T, re_formula = NA, newdata = ., value = 'deltaProd_hgrav') %>% 
  ungroup %>% select(hgrav2, deltaProd_hgrav) %>% 
  bind_cols(T_partialpred, .) %>% ungroup

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


## POSTERIORS PREDICTIONS ==========================================================================

# get posterior draws
geo_draws <- tidy_draws(mod_geo) %>% 
  select(starts_with(c('b', 'sd')), shape, Intercept) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

geo_draws <- geo_draws %>% mutate(var_type = case_when(vars %in% c('b_logS', 'b_loghgrav2P1', 'sd_biogeo__logS') ~ 'main',
                                                               .default = 'randommodel'),
                                          prettyvar = case_when(vars == 'b_logS' ~ 'S',
                                                                vars == 'b_loghgrav2P1' ~ 'Human gravity',
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
            minG = min(hgrav2), maxG = max(hgrav2))
for (i in 1:nrow(temp)) {
  if (i > 1) {
    newdata <- tidytable(biogeo = temp$biogeo[i], 
                         S = runif(100, temp$minS[i], temp$maxS[i]) %>% floor(),
                         hgrav2 = runif(100, temp$minG[i], temp$maxG[i]) %>% floor()) %>% 
      bind_rows(newdata, .)
  } else {
    newdata <- tidytable(biogeo = temp$biogeo[i], 
                         S = runif(100, temp$minS[i], temp$maxS[i]) %>% floor(),
                         hgrav2 = runif(100, temp$minG[i], temp$maxG[i]))
  }
}


# partial effect of ∆S
geo_partialpred <- newdata %>% mutate(hgrav2 = mean(vars_single$hgrav2)) %>% 
  epred_draws(mod_geo, re_formula = prod_kg ~ log(S) + log(hgrav2 + 1) + (1 + log(S) | biogeo), 
              newdata = ., value = 'prod_S') %>% ungroup %>% 
  select(!hgrav2)

geo_partialpred <- tidytable(hgrav2 = runif(n_distinct(vars_single$S), min(vars_single$hgrav2, na.rm = T), 
                                                max(vars_single$hgrav2, na.rm = T)) %>% floor(),
                                 S = mean(vars_single$S, na.rm = T)) %>% 
  epred_draws(mod_geo, re_formula = NA, newdata = ., value = 'prod_hgrav') %>% 
  ungroup %>% select(hgrav2, prod_hgrav) %>% 
  bind_cols(spatial_partialpred, .) %>% ungroup

rm(newdata)

# replace spatial preds with this model

# partial effect of ∆S
spatial_partialpred <- vars_single %>% distinct(S) %>% 
  mutate(hgrav2 = mean(vars_single$hgrav2)) %>% 
  epred_draws(mod_geo, re_formula = NA, newdata = ., value = 'prod_S') %>% 
  ungroup %>% select(!hgrav2)
# partial effect of human gravity
spatial_partialpred <- tidytable(hgrav2 = runif(121, min(vars_single$hgrav2, na.rm = T), 
                                           max(vars_single$hgrav2, na.rm = T)),
                                 S = mean(vars_single$S, na.rm = T)) %>% 
  epred_draws(mod_geo, re_formula = NA, newdata = ., value = 'prod_hgrav') %>% 
  ungroup %>% select(hgrav2, prod_hgrav) %>% 
  bind_cols(spatial_partialpred, .) %>% ungroup

# residuals
spatial_resid <- vars_single %>% select(S, prod_kg)
spatial_resid <- spatial_partialpred %>% group_by(S) %>% summarise(prod_S = median(prod_S)) %>% 
  left_join(spatial_resid, ., by = 'S') %>% distinct
spatial_resid <- spatial_resid %>% mutate(resid = prod_S - prod_kg)

# MODEL 5: QUANTILE BIOGEO ==================================================================================
# what is the diversity-productivity relationship like across space

# weakly informative priors
geoQ_priors <- default_prior(bf(prod_kg ~ log(S) + log(hgrav2 + 1) + (1|DataID) + (1|gridcell), quantile = 0.25),
                            data = vars_single, 
                            family = asym_laplace(link = 'log'))
geoQ_priors
geoQ_priors$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative
geoQ_priors$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive

mod_geoQ <- as.list(rep(NA, 3))
set.seed(1)
mod_geoQ[[1]] <- brm(bf(prod_kg ~ log(S) + log(hgrav2 + 1) + (1 + log(S)|biogeo) + (1|DataID) + (1|gridcell), quantile = 0.25),
               data = vars_single, 
               family = asym_laplace(link = 'log'),
               prior = geoQ_priors, iter = 5000)

set.seed(1)
mod_geoQ[[2]] <- brm(bf(prod_kg ~ log(S) + log(hgrav2 + 1) + (1 + log(S)|biogeo) + (1|DataID) + (1|gridcell), quantile = 0.5),
                     data = vars_single,
                     family = asym_laplace(link = 'log'),
                     prior = geoQ_priors, iter = 5000)

set.seed(1)
mod_geoQ[[3]] <- brm(bf(prod_kg ~ log(S) + log(hgrav2 + 1) + (1 + log(S)|biogeo) + (1|DataID) + (1|gridcell), quantile = 0.75),
                     data = vars_single, 
                     family = asym_laplace(link = 'log'),
                     prior = geoQ_priors, iter = 5000)

## POSTERIORS PREDICTIONS ==========================================================================

# get posterior draws
geoq_draws <- tidy_draws(mod_geoQ[[1]]) %>% 
  select(starts_with(c('b', 'sd')), sigma) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

geoq_draws <- geoq_draws %>% mutate(var_type = case_when(vars %in% c('b_logS', 'b_loghgrav2P1', 'sd_biogeoq__logS') ~ 'main',
                                                       .default = 'randommodel'),
                                  prettyvar = case_when(vars == 'b_logS' ~ 'S',
                                                        vars == 'b_loghgrav2P1' ~ 'Human gravity',
                                                        vars == 'sd_biogeoq__logS' ~ '(S|region)',
                                                        vars == 'sd_biogeoq__Intercept' ~ '(1|region)',
                                                        vars == 'sd_DataID__Intercept' ~ '(1|study)',
                                                        vars == 'sd_gridcell__Intercept' ~ '(1|location)',
                                                        vars == 'cor_biogeo__Intercept__logS' ~ 'cor(region, S, β0)',
                                                        vars == 'b_Intercept' ~ 'β0',
                                                        .default = vars))
# partial effect of ∆S
geoq_partialpred <- vars_single %>% distinct(S) %>% 
  mutate(hgrav2 = mean(vars_single$hgrav2)) %>% 
  epred_draws(mod_geoQ[[1]], re_formula = NA, newdata = ., value = 'prod_S', dpar = 'mu')
