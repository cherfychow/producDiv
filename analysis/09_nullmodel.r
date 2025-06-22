

# ProducDiv
# Fitting null model

require(tidytable)
require(brms)
require(rfishprod)
require(tidybayes)

dt_null <- vars_single %>% select(DataID, seriesID, SamplingEventID, S, hgrav2, gridcell)

# NULL PRODUCTIVITY/SIZE ===========================================================================

summary(prod_single)

# Kmax is estimated per species but now integrate with individual length records
dt_prod_null <- dt_single %>% 
  select(DataID, seriesID, fishbase_valid, SamplingEventID, Abundance) %>% 
  rename(Species = fishbase_valid)
dt_prod_null <- left_join(dt_prod_null, prod_single, by = c('Species', 'DataID', 'seriesID', 'SamplingEventID'), 
                            na_matches = 'na', unmatched = 'error', relationship = 'many-to-one') %>% 
  select(!c(Diet, Position))
summary(dt_prod_null)

# remove the records that don't qualify for productivity calculation
dt_prod_null <- dt_prod_null %>% filter(!is.na(Kmax)) %>% 
  mutate(Lmax = 0.9999 * MaxSizeTL,
         Lmin = 0.0001 * MaxSizeTL) # as close to zero as we can get productivity

summary(dt_prod_null)

# lowers and uppers refer to the Kmax value used!
# not necessarily the lower or upper of that metric
dt_prod_null <- dt_prod_null %>% 
  mutate(prod0_Lmin = somaGain(a, b, Lmeas = Lmin, Lmax = MaxSizeTL, Kmax = Kmax),
         prod0_Lmin_l = somaGain(a, b, Lmeas = Lmin, Lmax = MaxSizeTL, Kmax = Kmax_lowq),
         prod0_Lmin_u = somaGain(a, b, Lmeas = Lmin, Lmax = MaxSizeTL, Kmax = Kmax_uppq),
         prod0_Lmax = somaGain(a, b, Lmeas = Lmax, Lmax = MaxSizeTL, Kmax = Kmax),
         prod0_Lmax_l = somaGain(a, b, Lmeas = Lmax, Lmax = MaxSizeTL, Kmax = Kmax_lowq),
         prod0_Lmax_u = somaGain(a, b, Lmeas = Lmax, Lmax = MaxSizeTL, Kmax = Kmax_uppq)
  )

summary(dt_prod_null)

dt_prod_null <- dt_prod_null %>% select(!c(sstmean:b))

## Apply to assemblage data 

# scale by abundance
dt_prod_null <- dt_prod_null %>% 
  mutate(across(starts_with('prod0'), ~ Abundance * .x)) %>% 
  select(!Abundance)

# aggregate to sampling event level
dt_prod_null <- dt_prod_null %>% ungroup %>% 
  group_by(DataID, seriesID, SamplingEventID) %>% 
  summarise(across(starts_with('prod0'), ~ sum(.x))) %>% ungroup
summary(dt_prod_null)

# convert to kg/day
dt_prod_null <- dt_prod_null %>% mutate(across(starts_with('prod0'), ~ .x / 1000))

# kg/ha/day
dt_prod_null <- dt_meta %>% distinct(DataID, `Grain size`) %>% 
  mutate(Area = str_extract(`Grain size`, '[0-9\\.]+') %>% as.numeric) %>% 
  left_join(dt_prod_null, ., by = "DataID") %>% 
  select(!`Grain size`)

summary(dt_prod_null)
# unit conversion
dt_prod_null <- dt_prod_null %>% mutate(across(starts_with('prod0'), ~ .x * 10000 / Area))

# just add it onto vars_single
dt_null <- dt_prod_null %>% select(!Area) %>% 
  left_join(dt_null, ., by = c('DataID', 'seriesID', 'SamplingEventID'))

## NULL PRODUCTIVITY MODEL ==========================================================================
# weakly informative priors
null_priors <- default_prior(prod0_Lmin ~ log(S) + log(hgrav2 + 1) + (1|DataID) + (1|gridcell),
                                data = dt_null, family = Gamma(link = 'log'))
null_priors$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative
null_priors$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive

set.seed(1)
mod_null <- brm(data = dt_null,
                   formula = prod0_Lmin ~ log(S) + log(hgrav2 + 1) + (1|DataID) + (1|gridcell),
                   family = Gamma(link = 'log'),
                   prior = null_priors, iter = 5000)


## POSTERIORS PREDICtIONS 

# get posterior draws
null_draws <- tidy_draws(mod_null) %>% 
  select(starts_with(c('b', 'sd', '.')), shape, Intercept) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

null_draws <- null_draws %>% mutate(var_type = case_when(str_detect(vars, 'log') ~ 'fixed',
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
null_partialpred <- newdata %>% mutate(hgrav2 = mean(vars_single$hgrav2)) %>% 
  epred_draws(mod_null, re_formula = NA, newdata = ., value = 'prod_S') %>% ungroup %>% 
  select(!hgrav2)
# partial effect of human gravity
null_partialpred <- newdata %>% mutate(S = mean(vars_single$S, na.rm = T)) %>% 
  epred_draws(mod_null, re_formula = NA, newdata = ., value = 'prod_hgrav') %>% 
  ungroup %>% select(hgrav2, prod_hgrav) %>% 
  bind_cols(null_partialpred, .) %>% ungroup
rm(newdata)

## NULL MAX PRODUCTIVITY MODEL =====================================================================
# weakly informative priors

nullP2_priors <- default_prior(prod0_Lmin ~ log(S) + log(hgrav2 + 1) + (1|DataID) + (1|gridcell:DataID),
                             data = dt_null, family = Gamma(link = 'log'))
nullP2_priors$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative
nullP2_priors$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive

set.seed(1)
mod_nullP2 <- brm(data = dt_null,
                formula = prod0_Lmin ~ log(S) + log(hgrav2 + 1) + (1|DataID) + (1|gridcell:DataID),
                family = Gamma(link = 'log'),
                prior = nullP2_priors, iter = 5000)

# Model posteriors + predictions
nullP2_draws <- tidy_draws(mod_nullP2) %>% 
  select(starts_with(c('b', 'sd', '.')), shape) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

nullP2_draws <- nullP2_draws %>% mutate(var_type = case_when(str_detect(vars, 'log') ~ 'fixed',
                                                         .default = 'randommodel'),
                                    prettyvar = case_when(vars == 'b_logS' ~ 'S',
                                                          vars == 'b_loghgrav2P1' ~ 'Human gravity',
                                                          vars == 'sd_DataID__Intercept' ~ '(1|study)',
                                                          vars == 'sd_gridcell:DataID__Intercept' ~ '(1|location:study)',
                                                          vars == 'b_Intercept' ~ 'β0',
                                                          .default = vars))

# model partial predictions
newdata <- tidytable(S = runif(500, min(vars_single$S, na.rm = T), 
                               max(vars_single$S, na.rm = T)) %>% floor(),
                     hgrav2 = runif(500, min(vars_single$hgrav2), max(vars_single$hgrav2)))

# partial effect of ∆S
nullP2_partialpred <- newdata %>% mutate(hgrav2 = mean(vars_single$hgrav2)) %>% 
  epred_draws(mod_nullP2, re_formula = NA, newdata = ., value = 'prod0_nullP2') %>% ungroup %>% 
  select(!hgrav2)

# NULL SPECIES =====================================================================================

# Prep predictors 
dt_sp_null <- left_join(dt_single, dt_traits[c(1,3:5)], by = 'taxaID') %>% 
  filter(!is.na(MaxSizeTL)) %>% 
  select(DataID, seriesID, SamplingEventID, Abundance, Length, a, b, MaxSizeTL) %>% 
  mutate(Length = case_when(Length >= MaxSizeTL ~ .97*MaxSizeTL,
                            .default = (Length/MaxSizeTL)*MaxSizeTL)) # turn length relative to MaxSize

dt_sp_null <- prod_single %>% ungroup %>% 
  group_by(DataID, seriesID, SamplingEventID) %>% 
  summarise(Kmax_max = max(Kmax, na.rm = T), Kmax_min = min(Kmax, na.rm = T)) %>% # one kmax for the whole assemblage
  left_join(dt_sp_null, ., by = c('DataID', 'seriesID', 'SamplingEventID')) %>% 
  filter(!is.na(Kmax_min)) # some species don't qualify
summary(dt_sp_null)

# species null productivity
dt_sp_null <- dt_sp_null %>% 
  mutate(prod0_spmax = somaGain(a, b, Lmeas = Length, Kmax = Kmax_max, Lmax = MaxSizeTL),
         prod0_spmin = somaGain(a, b, Lmeas = Length, Kmax = Kmax_min, Lmax = MaxSizeTL)) %>% 
  mutate(prod0_spmax = Abundance * prod_maxsp / 1000,
         prod0_spmin = Abundance * prod_minsp / 1000) %>% 
  ungroup %>% group_by(DataID, seriesID, SamplingEventID) %>% 
  summarise(prod0_spmax = sum(prod0_spmax), prod0_spmin = sum(prod0_spmin))

# area standardised
dt_sp_null <- dt_meta %>% distinct(DataID, `Grain size`) %>% 
  mutate(Area = str_extract(`Grain size`, '[0-9\\.]+') %>% as.numeric) %>% 
  left_join(dt_sp_null, ., by = "DataID") %>% 
  select(!`Grain size`) %>% 
  mutate(across(starts_with('prod'), ~ .x * 10000 / Area)) %>% 
  select(!Area)

dt_null <- dt_sp_null %>% ungroup %>% 
  left_join(dt_null, ., by = c('DataID', 'seriesID', 'SamplingEventID'))

# species null model

# weakly informative priors
nullsp1_priors <- default_prior(prod_maxsp ~ log(S) + log(hgrav2 + 1) + (1|DataID) + (1|gridcell),
                                data = vars_single, family = Gamma(link = 'log'))
nullsp1_priors$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative
nullsp1_priors$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive

set.seed(1)
mod_nullsp1 <- brm(data = vars_single,
                   formula = prod_maxsp ~ log(S) + log(hgrav2 + 1) + (1|DataID) + (1|gridcell),
                   family = Gamma(link = 'log'),
                   prior = nullsp1_priors, iter = 5000)

# weakly informative priors
nullsp2_priors <- default_prior(prod_minsp ~ log(S) + log(hgrav2 + 1) + (1|DataID) + (1|gridcell),
                                data = vars_single, family = Gamma(link = 'log'))
nullsp2_priors$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative
nullsp2_priors$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive

set.seed(1)
mod_nullsp2 <- brm(data = vars_single,
                   formula = prod_minsp ~ log(S) + log(hgrav2 + 1) + (1|DataID) + (1|gridcell),
                   family = Gamma(link = 'log'),
                   prior = nullsp2_priors, iter = 5000)

# Model draws and predictions

# get posterior draws
nullsp1_draws <- tidy_draws(mod_nullsp1) %>% 
  select(starts_with(c('b', 'sd', '.')), shape, Intercept) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

nullsp1_draws <- nullsp1_draws %>% mutate(var_type = case_when(str_detect(vars, 'log') ~ 'fixed',
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
nullsp1_partial <- newdata %>% mutate(hgrav2 = mean(vars_single$hgrav2)) %>% 
  epred_draws(mod_nullsp1, re_formula = NA, newdata = ., value = 'prod_S') %>% ungroup %>% 
  select(!hgrav2)
# partial effect of human gravity
nullsp1_partial <- newdata %>% mutate(S = mean(vars_single$S, na.rm = T)) %>% 
  epred_draws(mod_nullsp1, re_formula = NA, newdata = ., value = 'prod_hgrav') %>% 
  ungroup %>% select(hgrav2, prod_hgrav) %>% 
  bind_cols(nullsp1_partial, .) %>% ungroup
rm(newdata)


# get posterior draws
nullsp2_draws <- tidy_draws(mod_nullsp2) %>% 
  select(starts_with(c('b', 'sd', '.')), shape, Intercept) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

nullsp2_draws <- nullsp2_draws %>% mutate(var_type = case_when(str_detect(vars, 'log') ~ 'fixed',
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
nullsp2_partial <- newdata %>% mutate(hgrav2 = mean(vars_single$hgrav2)) %>% 
  epred_draws(mod_nullsp2, re_formula = NA, newdata = ., value = 'prod_S') %>% ungroup %>% 
  select(!hgrav2)
# partial effect of human gravity
nullsp2_partial <- newdata %>% mutate(S = mean(vars_single$S, na.rm = T)) %>% 
  epred_draws(mod_nullsp2, re_formula = NA, newdata = ., value = 'prod_hgrav') %>% 
  ungroup %>% select(hgrav2, prod_hgrav) %>% 
  bind_cols(nullsp2_partial, .) %>% ungroup
rm(newdata)

# NULL ABUNDANCE ===================================================================================

# total N
N <- dt_single %>% select(DataID, seriesID, SamplingEventID, Abundance) %>% 
  group_by(DataID, seriesID, SamplingEventID) %>% summarise(totalN = sum(Abundance))

# all abundance = 1

dt_abund0 <- dt_single %>% 
  select(DataID, seriesID, fishbase_valid, SamplingEventID, Length) %>% 
  rename(Species = fishbase_valid) %>% distinct

dt_abund0 <- left_join(dt_abund0, prod_single, by = c('Species', 'DataID', 'seriesID', 'SamplingEventID'), 
                          na_matches = 'na', unmatched = 'error', relationship = 'many-to-one') %>% 
  select(!c(Diet, Position))
summary(dt_abund0)

dt_abund0 <- dt_abund0 %>% filter(!is.na(Kmax)) %>% 
  mutate(Length = case_when(.default = Length,
                            Length > MaxSizeTL ~ MaxSizeTL - 0.1))

# lowers and uppers refer to the Kmax value used!
# not necessarily the lower or upper of that metric
dt_abund0 <- dt_abund0 %>% 
  mutate(prod0_Amin = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax),
         prod0_Amin_l = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax_lowq),
         prod0_Amin_u = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax_uppq)
  )

summary(dt_abund0)

# aggregate to sampling event level
dt_abund0 <- dt_abund0 %>% ungroup %>% 
  group_by(DataID, seriesID, SamplingEventID) %>% 
  summarise(across(starts_with('prod0'), ~ sum(.x))) %>% ungroup
summary(dt_abund0)

# convert to kg/day
dt_abund0 <- dt_abund0 %>% mutate(across(starts_with('prod0'), ~ .x / 1000))

# kg/ha/day
dt_abund0 <- dt_meta %>% distinct(DataID, `Grain size`) %>% 
  mutate(Area = str_extract(`Grain size`, '[0-9\\.]+') %>% as.numeric) %>% 
  left_join(dt_abund0, ., by = "DataID") %>% 
  select(!`Grain size`)

summary(dt_abund0)
# unit conversion
dt_abund0 <- dt_abund0 %>% mutate(across(starts_with('prod0'), ~ .x * 10000 / Area))

dt_null <- dt_abund0 %>% select(!Area) %>% 
  left_join(dt_null, ., by = c('DataID', 'seriesID', 'SamplingEventID'))

## Fit model 

abund0_priors <- default_prior(prod0_Amin ~ log(S) + log(hgrav2 + 1) + (1|DataID) + (1|gridcell:DataID),
                               data = dt_null, family = Gamma(link = 'log'))
abund0_priors$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative
abund0_priors$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive

set.seed(1)
mod_abund0 <- brm(data = dt_null,
                  formula = prod0_Amin ~ log(S) + log(hgrav2 + 1) + (1|DataID) + (1|gridcell:DataID),
                  family = Gamma(link = 'log'),
                  prior = abund0_priors, iter = 5000)
# save(mod_abund0, file = 'mod_null_A.RData')

# posteriors and predictions

abund0_draws <- tidy_draws(mod_abund0) %>% 
  select(starts_with(c('b', 'sd', '.')), shape) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

abund0_draws <- abund0_draws %>% mutate(var_type = case_when(str_detect(vars, 'log') ~ 'fixed',
                                                               .default = 'randommodel'),
                                          prettyvar = case_when(vars == 'b_logS' ~ 'S',
                                                                vars == 'b_loghgrav2P1' ~ 'Human gravity',
                                                                vars == 'sd_DataID__Intercept' ~ '(1|study)',
                                                                vars == 'sd_gridcell:DataID__Intercept' ~ '(1|location:study)',
                                                                vars == 'b_Intercept' ~ 'β0',
                                                                .default = vars))

# model partial predictions
newdata <- tidytable(S = runif(500, min(vars_single$S, na.rm = T), 
                               max(vars_single$S, na.rm = T)) %>% floor(),
                     hgrav2 = mean(vars_single$hgrav2))

# partial effect of ∆S
abund0_partial <- epred_draws(mod_abund0, re_formula = NA, newdata = newdata, value = 'prod_S') %>% ungroup %>% 
  select(!hgrav2)
