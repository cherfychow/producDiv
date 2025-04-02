

# ProducDiv
# Fitting null model

require(tidytable)
require(brms)
require(rfishprod)

# Use the outputs of single-time data productivity to calculate a null productivity ================

summary(prod_single)

# Kmax is estimated per species but now integrate with individual length records
dt_prod_null <- dt_single %>% 
  select(DataID, seriesID, taxaID, SamplingEventID, Abundance)
dt_prod_null <- left_join(dt_prod_null, prod_single, by = c('taxaID', 'DataID', 'seriesID', 'SamplingEventID'), 
                            na_matches = 'na', unmatched = 'error', relationship = 'many-to-one') %>% 
  select(!c(Diet, Position))
summary(dt_prod_null)

# remove the records that don't qualify for productivity calculation
dt_prod_null <- dt_prod_null %>% filter(!is.na(Kmax)) %>% 
  mutate(Length = 0.9999 * MaxSizeTL) # as close to zero as we can get productivity

summary(dt_prod_null)

# lowers and uppers refer to the Kmax value used!
# not necessarily the lower or upper of that metric
dt_prod_null <- dt_prod_null %>% 
  mutate(biomass_gain = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax),
         biomass_gain_l = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax_lowq),
         biomass_gain_u = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax_uppq)
  )

summary(dt_prod_null)

dt_prod_null <- dt_prod_null %>% select(!c(Species:b))

## Apply to assemblage data =====================================================

# scale by abundance
dt_prod_null <- dt_prod_null %>% 
  mutate(across(starts_with('biomass'), ~ Abundance * .x)) %>% 
  select(!Abundance)

# aggregate to sampling event level
dt_prod_null <- dt_prod_null %>% ungroup %>% 
  group_by(DataID, seriesID, SamplingEventID) %>% 
  summarise(across(starts_with('biomass'), ~ sum(.x))) %>% ungroup
summary(dt_prod_null)

# convert to kg/day
dt_prod_null <- dt_prod_null %>% mutate(across(starts_with('biomass'), ~ .x / 1000))

colnames(dt_prod_null) <- str_replace(colnames(dt_prod_null), 'biomass\\_gain', 'prod_kg')

dt_meta$DataID <- str_to_lower(dt_meta$DataID)
dt_prod_null$DataID <- str_to_lower(dt_prod_null$DataID)
# kg/ha/day
dt_prod_null <- dt_meta %>% distinct(DataID, `Grain size`) %>% 
  mutate(Area = str_extract(`Grain size`, '[0-9\\.]+') %>% as.numeric) %>% 
  left_join(dt_prod_null, ., by = "DataID") %>% 
  select(!`Grain size`)

summary(dt_prod_null)
# unit conversion
dt_prod_null <- dt_prod_null %>% mutate(across(starts_with('prod'), ~ .x * 10000 / Area))

## Prep predictors =================================================================================

# just add it onto vars_single
colnames(dt_prod_null) <- str_replace(colnames(dt_prod_null), '\\_kg', '\\_null')
vars_single <- dt_prod_null %>% select(!Area) %>% 
  left_join(vars_single, ., by = c('DataID', 'seriesID', 'SamplingEventID'))
summary(vars_single)

# NULL SPATIAL MODEL ===============================================================================

# weakly informative priors
null_priors <- default_prior(prod_null ~ log(S) + log(hgrav2 + 1) + (1|DataID) + (1|gridcell),
                                data = vars_single, family = Gamma(link = 'log'))
null_priors$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative
null_priors$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive

set.seed(1)
mod_null <- brm(data = vars_single,
                   formula = prod_null ~ log(S) + log(hgrav2 + 1) + (1|DataID) + (1|gridcell),
                   family = Gamma(link = 'log'),
                   prior = null_priors, iter = 5000)


## POSTERIORS PREDICtIONS ==========================================================================

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
