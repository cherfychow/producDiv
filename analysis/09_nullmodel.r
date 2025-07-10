

# ProducDiv
# Fitting null model

require(tidytable)
require(brms)
require(rfishprod)
require(tidybayes)

dt_null <- vars_single %>% select(DataID, seriesID, SamplingEventID, S, hgravity, gridcell)

# NULL PRODUCTIVITY/SIZE ===========================================================================

summary(dt_Kmax)

# Kmax is estimated per species but now integrate with individual length records
dt_prod_null <- dt_single %>% 
  select(DataID, seriesID, fishbase_valid, SamplingEventID, Abundance) %>% 
  rename(Species = fishbase_valid)
dt_prod_null <- left_join(dt_prod_null, dt_Kmax, by = c('Species', 'DataID', 'seriesID', 'SamplingEventID'), 
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
  mutate(prod0_Lmax = somaGain(a, b, Lmeas = Lmax, Lmax = MaxSizeTL, Kmax = Kmax),
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

# NULL PRODUCTIVITY MODEL ==========================================================================

load('dt_null.RData')
require(brms)

# weakly informative priors
pr_prod0 <- default_prior(prod0_Lmax ~ log(S) + log(hgravity + .000039) + (1|DataID) + (1|gridcell) + (1|grav_ID),
                                data = dt_null, family = Gamma(link = 'log'))
pr_prod0$prior[2] <- 'normal(-0.1, 3)' # Market gravity, slight negative
pr_prod0$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive

set.seed(1)
m_prod0 <- brm(data = dt_null,
                   formula = prod0_Lmax ~ log(S) + log(hgravity + .000039) + (1|DataID) + (1|gridcell) + (1|grav_ID),
                   family = Gamma(link = 'log'),
                   prior = pr_prod0, iter = 5000)

save(m_prod0, file = 'm_prod0.RData')

## POSTERIORS PREDICtIONS 

# get posterior draws
dr_prod0 <- tidy_draws(m_prod0) %>% 
  select(starts_with(c('b', 'sd', '.')), shape, Intercept) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

dr_prod0 <- dr_prod0 %>% mutate(var_type = case_when(str_detect(vars, 'log') ~ 'fixed',
                                                               .default = 'randommodel'),
                                          prettyvar = case_when(vars == 'b_logS' ~ 'S',
                                                                vars == 'b_loghgravityP3.9eM05' ~ 'Market gravity',
                                                                vars == 'sd_DataID__Intercept' ~ '(1|DataID)',
                                                                vars == 'sd_gridcell__Intercept' ~ '(1|gridcell)',
                                                                vars == 'b_Intercept' ~ 'β0',
                                                                .default = vars))

# model partial predictions
newdata <- tidytable(S = runif(500, min(vars_single$S, na.rm = T), 
                               max(vars_single$S, na.rm = T)) %>% floor(),
                     hgravity = mean(vars_single$hgravity))

# partial effect of ∆S
pred_prod0 <- newdata %>% epred_draws(m_prod0, re_formula = NA, newdata = ., value = 'prod_S') %>% ungroup %>% 
  select(!hgravity)

rm(newdata)

# NULL SPECIES =====================================================================================

# Prep predictors 
temp <- left_join(dt_single, dt_Kmax, 
          by = c('gridcell', 'Year', 'Month', 'fishbase_valid' = 'Species'))

K_min <- temp %>% ungroup %>% arrange(DataID, seriesID, SamplingEventID, MaxSizeTL, Kmax) %>% 
  group_by(DataID, seriesID, SamplingEventID, Year) %>% 
  select(MaxSizeTL:Kmax) %>% 
  slice_head(n = 1) %>% ungroup # minimum

K_max <- temp %>% ungroup %>% arrange(DataID, seriesID, SamplingEventID, desc(MaxSizeTL), desc(Kmax)) %>% 
  group_by(DataID, seriesID, SamplingEventID, Year) %>% 
  select(MaxSizeTL:Kmax) %>% 
  slice_head(n = 1) %>% ungroup # maximum

dt_spmax <- left_join(dt_single, dt_Kmax, 
                           by = c('gridcell', 'Year', 'Month', 'fishbase_valid' = 'Species')) %>% 
  filter(!is.na(MaxSizeTL)) %>%  # remove records without max size information
  mutate(Length = case_when(Length >= MaxSizeTL ~ 0.99,
                            .default = Length/MaxSizeTL)) %>%  # turn length relative to MaxSize
  select(Abundance, Length, Year, DataID, seriesID, SamplingEventID) %>% 
  group_by(-Abundance) %>% # aggregate for simplicity since everyone will be one species anyway
  summarise(Abundance = sum(Abundance)) %>% ungroup %>% 
  left_join(., K_max, by = c('DataID', 'seriesID', 'SamplingEventID', 'Year')) %>% ungroup %>% 
  arrange(DataID, seriesID, SamplingEventID, Year, Length)

dt_spmax$Length <- with(dt_spmax, Length * MaxSizeTL) # scale back relative length

dt_spmin <- left_join(dt_single, dt_Kmax, 
                           by = c('gridcell', 'Year', 'Month', 'fishbase_valid' = 'Species')) %>% 
  filter(!is.na(MaxSizeTL)) %>%
  mutate(Length = case_when(Length >= MaxSizeTL ~ 0.99,
                            .default = Length/MaxSizeTL)) %>%  # turn length relative to MaxSize
  select(Abundance, Length, Year, DataID, seriesID, SamplingEventID) %>% 
  group_by(-Abundance) %>% # aggregate for simplicity since everyone will be one species anyway
  summarise(Abundance = sum(Abundance)) %>% ungroup %>% 
  left_join(., K_min, by = c('DataID', 'seriesID', 'SamplingEventID', 'Year')) %>% ungroup %>% 
  arrange(DataID, seriesID, SamplingEventID, Year, Length)
dt_spmin$Length <- with(dt_spmin, Length * MaxSizeTL) # scale back relative length

summary(dt_spmin)
summary(dt_spmax)

# species null productivity
dt_spmin <- dt_spmin %>% 
  mutate(prod0_spmin = somaGain(a, b, Lmeas = Length, Kmax = Kmax, Lmax = MaxSizeTL)) %>% 
  mutate(prod0_spmin = Abundance * prod0_spmin / 1000) %>%  # scale by abundance and convert to kg
  ungroup %>% group_by(DataID, seriesID, SamplingEventID) %>% 
  summarise(prod0_spmin = sum(prod0_spmin))

dt_spmax <- dt_spmax %>% 
  mutate(prod0_spmax = somaGain(a, b, Lmeas = Length, Kmax = Kmax, Lmax = MaxSizeTL)) %>% 
  mutate(prod0_spmax = Abundance * prod0_spmax / 1000) %>%  # scale by abundance and convert to kg
  ungroup %>% group_by(DataID, seriesID, SamplingEventID) %>% 
  summarise(prod0_spmax = sum(prod0_spmax))

# area standardised
dt_spmin <- dt_meta %>% distinct(DataID, `Grain size`) %>% 
  mutate(Area = str_extract(`Grain size`, '[0-9\\.]+') %>% as.numeric) %>% 
  left_join(dt_spmin, ., by = "DataID") %>% 
  select(!`Grain size`) %>% 
  mutate(across(starts_with('prod'), ~ .x * 10000 / Area)) %>% # standardise by hectare
  select(!Area)

dt_spmax <- dt_meta %>% distinct(DataID, `Grain size`) %>% 
  mutate(Area = str_extract(`Grain size`, '[0-9\\.]+') %>% as.numeric) %>% 
  left_join(dt_spmax, ., by = "DataID") %>% 
  select(!`Grain size`) %>% 
  mutate(across(starts_with('prod'), ~ .x * 10000 / Area)) %>% # standardise by hectare
  select(!Area)

dt_null <- dt_spmax %>% ungroup %>% 
  left_join(dt_null, ., by = c('DataID', 'seriesID', 'SamplingEventID')) %>% 
  left_join(., dt_spmin, by = c('DataID', 'seriesID', 'SamplingEventID'))

## Fit species null models =========================================================================

load('dt_null.RData')
require(brms)

# weakly informative priors
pr_sp0_max <- default_prior(prod0_spmax ~ log(S) + log(hgravity + .000039) + (1|DataID) + (1|gridcell) + (1|grav_ID),
                                data = dt_null, family = Gamma(link = 'log'))
pr_sp0_max$prior[2] <- 'normal(-0.1, 3)' # Market gravity, slight negative
pr_sp0_max$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive

set.seed(1)
m_sp0_max <- brm(data = dt_null,
                   formula = prod0_spmax ~ log(S) + log(hgravity + .000039) + (1|DataID) + (1|gridcell) + (1|grav_ID),
                   family = Gamma(link = 'log'),
                   prior = pr_sp0_max, iter = 5000)

save(m_sp0_max, file = 'm_sp0_max.RData')

load('dt_null.RData')
require(brms)

# weakly informative priors
pr_sp0_min <- default_prior(prod0_spmin ~ log(S) + log(hgravity + .000039) + (1|DataID) + (1|gridcell) + (1|grav_ID),
                                data = dt_null, family = Gamma(link = 'log'))
pr_sp0_min$prior[2] <- 'normal(-0.1, 3)' # Market gravity, slight negative
pr_sp0_min$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive

set.seed(1)
m_sp0_min <- brm(data = dt_null,
                   formula = prod0_spmin ~ log(S) + log(hgravity + .000039) + (1|DataID) + (1|gridcell) + (1|grav_ID),
                   family = Gamma(link = 'log'),
                   prior = pr_sp0_min, iter = 5000)
save(m_sp0_min, file = 'm_sp0_min.RData')

# Model draws and predictions

# get posterior draws
dr_sp0_max <- tidy_draws(m_sp0_max) %>% 
  select(starts_with(c('b', 'sd', '.')), shape, Intercept) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

dr_sp0_max <- dr_sp0_max %>% mutate(var_type = case_when(str_detect(vars, 'log') ~ 'fixed',
                                                         .default = 'randommodel'),
                                    prettyvar = case_when(vars == 'b_logS' ~ 'S',
                                                          vars == 'b_loghgravityP3.9eM05' ~ 'Market gravity',
                                                          vars == 'sd_DataID__Intercept' ~ '(1|DataID)',
                                                          vars == 'sd_gridcell__Intercept' ~ '(1|gridcell)',
                                                          vars == 'Intercept' ~ 'α',
                                                          vars == 'b_Intercept' ~ 'β0',
                                                          .default = vars))

# model partial predictions
newdata <- tidytable(S = runif(500, min(vars_single$S, na.rm = T), 
                               max(vars_single$S, na.rm = T)) %>% floor(),
                     hgravity = mean(vars_single$hgravity))

# partial effect of ∆S
pred_sp0_max <- newdata %>% mutate(hgravity = mean(vars_single$hgravity)) %>% 
  epred_draws(m_sp0_max, re_formula = NA, newdata = ., value = 'prod_S') %>% ungroup %>% 
  select(!hgravity)
# partial effect of Market gravity
pred_sp0_max <- newdata %>% mutate(S = mean(vars_single$S, na.rm = T)) %>% 
  epred_draws(m_sp0_max, re_formula = NA, newdata = ., value = 'prod_hgrav') %>% 
  ungroup %>% select(hgravity, prod_hgrav) %>% 
  bind_cols(pred_sp0_max, .) %>% ungroup
rm(newdata)


# get posterior draws
dr_sp0_min <- tidy_draws(m_sp0_min) %>% 
  select(starts_with(c('b', 'sd', '.')), shape, Intercept) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

dr_sp0_min <- dr_sp0_min %>% mutate(var_type = case_when(str_detect(vars, 'log') ~ 'fixed',
                                                               .default = 'randommodel'),
                                          prettyvar = case_when(vars == 'b_logS' ~ 'S',
                                                                vars == 'b_loghgravityP3.9eM05' ~ 'Market gravity',
                                                                vars == 'sd_DataID__Intercept' ~ '(1|DataID)',
                                                                vars == 'sd_gridcell__Intercept' ~ '(1|gridcell)',
                                                                vars == 'Intercept' ~ 'α',
                                                                vars == 'b_Intercept' ~ 'β0',
                                                                .default = vars))

# model partial predictions
newdata <- tidytable(S = runif(500, min(vars_single$S, na.rm = T), 
                               max(vars_single$S, na.rm = T)) %>% floor(),
                     hgravity = mean(vars_single$hgravity))

# partial effect of ∆S
pred_sp0_min <- newdata %>% mutate(hgravity = mean(vars_single$hgravity)) %>% 
  epred_draws(m_sp0_min, re_formula = NA, newdata = ., value = 'prod_S') %>% ungroup %>% 
  select(!hgravity)

ggplot(data = dt_null) + 
  # geom_point(data = vars_single, aes(x = S, y = prod_kg), 
  #            colour = '#666', shape = 21, alpha = 0.3) + 
  geom_point(aes(x = S, y = prod0_spmin), shape = 21, colour = 'blue', alpha = 0.2) + 
  # min Kmax
  geom_ribbon(data = pred_sp0_min %>% ungroup %>% group_by(S) %>% 
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)), 
              aes(x = S, ymin = lwr, ymax = upr), 
              fill = '#000', alpha = 0.1) +
  
  geom_line(data = pred_sp0_min %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = '#000', alpha = 0.8, linewidth = 0.75) +
  scale_x_log10() + scale_y_log10() + style

rm(newdata)

# NULL ABUNDANCE ===================================================================================

# total N
N <- dt_single %>% select(DataID, seriesID, SamplingEventID, Abundance) %>% 
  group_by(DataID, seriesID, SamplingEventID) %>% summarise(totalN = sum(Abundance))

# all abundance = 1

dt_abund0 <- dt_single %>% select(!Abundance) %>% distinct %>% 
  left_join(., dt_Kmax, by = c('gridcell', 'Year', 'Month', 'fishbase_valid' = 'Species')) %>% 
  filter(!is.na(MaxSizeTL))

dt_abund0 <- dt_abund0 %>% filter(!is.na(Kmax)) %>% 
  mutate(Length = case_when(.default = Length,
                            Length > MaxSizeTL ~ MaxSizeTL - 0.1))

# lowers and uppers refer to the Kmax value used!
# not necessarily the lower or upper of that metric
dt_abund0 <- dt_abund0 %>% 
  mutate(prod0_a = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax),
         prod0_a_l = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax_lowq),
         prod0_a_u = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax_uppq)
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

require(brms)
load('dt_null.RData')

pr_abund0 <- default_prior(prod0_a ~ log(S) + log(hgravity + .000039) + (1|DataID) + (1|gridcell:DataID) + (1|grav_ID),
                               data = dt_null, family = Gamma(link = 'log'))
pr_abund0$prior[2] <- 'normal(-0.1, 3)' # Market gravity, slight negative
pr_abund0$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive

set.seed(1)
m_abund0 <- brm(data = dt_null,
                  formula = prod0_a ~ log(S) + log(hgravity + .000039) + (1|DataID) + (1|gridcell:DataID) + (1|grav_ID),
                  family = Gamma(link = 'log'),
                  prior = pr_abund0, iter = 5000)
save(m_abund0, file = 'm_abund0.RData')

# posteriors and predictions

dr_abund0 <- tidy_draws(m_abund0) %>% 
  select(starts_with(c('b', 'sd', '.')), shape) %>% 
  tidyr::pivot_longer(cols = !starts_with('.'), names_to = 'vars', values_to = 'values')

dr_abund0 <- dr_abund0 %>% mutate(var_type = case_when(str_detect(vars, 'log') ~ 'fixed',
                                                               .default = 'randommodel'),
                                          prettyvar = case_when(vars == 'b_logS' ~ 'S',
                                                                vars == 'b_loghgravityP3.9eM05' ~ 'Market gravity',
                                                                vars == 'sd_DataID__Intercept' ~ '(1|study)',
                                                                vars == 'sd_gridcell:DataID__Intercept' ~ '(1|location:study)',
                                                                vars == 'b_Intercept' ~ 'β0',
                                                                .default = vars))

# model partial predictions
newdata <- tidytable(S = runif(500, min(vars_single$S, na.rm = T), 
                               max(vars_single$S, na.rm = T)) %>% floor(),
                     hgravity = mean(vars_single$hgravity))

# partial effect of ∆S
pred_abund0 <- epred_draws(m_abund0, re_formula = NA, newdata = newdata, value = 'prod_S') %>% ungroup %>% 
  select(!hgravity)
