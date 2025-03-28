
# ProducDiv
# Calculate productivity from fish traits + SST with observed lengths

require(tidytable)
require(stringr)
require(rfishprod)
require(readr)

# don't mix up time series and single time data

# TIME SERIES FIRST =========================================================================

# give predKmax data on all fishes observed
# then calculate productivity per rareID

# need to keep some sampling info since SST are specific to location + time
dt_rareass <- dt_sp %>% select(taxaID, ends_with('valid'), taxon_resolution) %>% 
  left_join(dt_assemblage, ., by = "taxaID") # join in
summary(dt_rareass)

# traits
# use taxaID as link
dt_traits <- dt_sp %>% select(taxaID, fishbase_valid) %>% 
  filter(!is.na(fishbase_valid)) %>% 
  left_join(., dt_traits, by = 'fishbase_valid', relationship = 'many-to-one')

# unique taxa at a given location and time
prodtrait <- dt_rareass %>% 
  distinct(DataID, seriesID, taxaID, gridcell, SamplingEventID, Year, Month, date_res) %>%
  inner_join(dt_traits, ., by = c('taxaID'))
# automatically filters out the uncertains that don't qualify for trait matching
# inner join somehow didn't quite work
prodtrait <- prodtrait %>% filter(!is.na(Diet))

# SST
prodtrait <- left_join(prodtrait, dt_sst, by = c('gridcell', 'Year', 'Month', 'date_res')) %>% # SST
  select(!c(gridcell, date_res, Year, Month)) %>% # no longer need these space time links
  rename(sstmean = sst, Species = fishbase_valid) # name matching
prodtrait <- prodtrait %>% relocate(ends_with('ID'), .after = sstmean)

# our own version of tidytrait
prodtrait <- prodtrait %>% mutate(
  Diet = factor(Diet, levels = sort(unique(dt_prodDB$Diet))),
  Position = factor(Position, levels = sort(unique(dt_prodDB$Position))),
  Method = factor('Otolth', levels = sort(unique(dt_prodDB$Method)))
)

prodtrait <- prodtrait %>% arrange(DataID, SamplingEventID)

dt_prodDB <- dt_prodDB %>% mutate(
  Diet = factor(Diet, levels = sort(unique(dt_prodDB$Diet))),
  Position = factor(Position, levels = sort(unique(dt_prodDB$Position))),
  Method = factor(Method, levels = sort(unique(dt_prodDB$Method)))
)
# make db match


## Kmax, growth coef ============================================================

# Formula from Morais and Bellwood (2018)
growthf <- formula(~ sstmean + MaxSizeTL + Diet + Position)

# small tweak to predKmax
predK <- function(traits, dataset, fmod, params = NULL, niter, nrounds = 150, 
                  verbose = 0, print_every = 1000, return = c("pred", "relimp", "models"), lowq = 0.25, uppq = 0.75) {
  if (missing(dataset)) {
    dataset <- db
  }
  modmatnd <- stats::model.matrix(fmod, data = traits)[, -1]
  modmat <- stats::model.matrix(fmod, dataset)[, -1]
  xgbmod <- list()
  rel_imp <- list()
  pred.grid <- as.data.frame(matrix(ncol = niter, nrow = nrow(traits), 
                                    dimnames = list(NULL, paste0("Boot", 1:niter))))
  for (i in 1:niter) {
    set.seed(i)
    xgbmod[[i]] <- xgboost::xgboost(modmat, label = dataset$Kmax, nrounds = nrounds, params = params, verbose = verbose, 
                                    print_every = print_every)
    rel_imp[[i]] <- xgboost::xgb.importance(colnames(modmat), 
                                            model = xgbmod[[i]])
    pred.grid[, i] <- stats::predict(xgbmod[[i]], newdata = modmatnd)
    cat(paste("Bootstrapping the model, round", i), "\n")
  }
  df <- cbind(traits, Kmax = apply(pred.grid, 1, stats::median), 
              Kmax_lowq = apply(pred.grid, 1, stats::quantile, probs = lowq), 
              Kmax_uppq = apply(pred.grid, 1, stats::quantile, probs = uppq))
  ret <- list(pred = df, relimp = rel_imp, models = xgbmod)
  return(ret)
}

pars  <- rfishprod:::xgboostparams # make sure the parameters are loaded
pars$max_depth <- 9

set.seed(24)

x = 20000
# make an empty list to populate
I = ceiling(nrow(prodtrait)/x)

# iterate through the dataset
for (i in 1:I) {
  if (i == I) { # how to handle the last iteration
    temp <- predK(traits = prodtrait[((i-1)*x+1):nrow(prodtrait),], 
                          dataset = dt_prodDB, fmod = growthf, params = pars, niter = 1000, 
                          lowq = 0.025, uppq = 0.975, verbose = T, return = 'pred')
  } else {
    temp <- predK(traits = prodtrait[((i-1)*x+1):(i*x),], 
                          dataset = dt_prodDB, fmod = growthf, params = pars, niter = 1000, 
                          lowq = 0.025, uppq = 0.975, verbose = T, return = 'pred')
  }
gc()
readr::write_csv(temp$pred, file = 'dt_vbgf.csv', append = T) # write to an output file instead of to the environment (BIG objects)
print(paste0('Run ', i, ' complete and saved'))
}

dt_vbgf <- read_csv('dt_vbgf.csv')
dt_vbgf <- dt_vbgf %>% distinct(DataID, seriesID, SamplingEventID, taxaID, .keep_all = T)
# should just be one Kmax estimate per species-location-time. if SST varies in a sample, it needs to be one value per sampling event id
summary(dt_vbgf)
dt_vbgf <- dt_vbgf %>% select(!ID)
  
  
## Gross somatic gain =============================================

# Kmax is estimated per species but now integrate with individual length records

dt_prod_rare <- dt_rareass %>% 
  select(ID, DataID, seriesID, taxaID, SamplingEventID, Abundance, Length, Length_lwr, Length_upr)
dt_prod_rare <- left_join(dt_prod_rare, dt_vbgf, by = c('taxaID', 'DataID', 'seriesID', 'SamplingEventID'), 
                     na_matches = 'na', unmatched = 'error', relationship = 'many-to-one') %>% 
  select(!c(Diet, Position, Method))
summary(dt_prod_rare)
# the NAs are from species that didn't qualify for trait/prod
# remove those records and constrain lengths
dt_prod_rare <- dt_prod_rare %>% filter(!is.na(Kmax)) %>% 
  mutate(Length = case_when(Length >= MaxSizeTL ~ MaxSizeTL - 0.2, .default = Length))

summary(dt_prod_rare)

# lowers and uppers refer to the metric calculated with the Kmax upper and lowers
# not necessarily the metric's own error propagated

dt_prod_rare <- dt_prod_rare %>% 
  mutate(biomass_gain = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax),
         biomass_gain_l = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax_lowq),
         biomass_gain_u = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax_uppq)
         )

summary(dt_prod_rare)

# dt_prod %>% filter(delta_g_day > 10) %>% View # make sure the large increases are proportional
# dt_prod %>% filter(delta_g_day < 0) %>% View
# don't need to keep traits
dt_prod_rare <- dt_prod_rare %>% select(!c(Species:sstmean))


## Apply to assemblage data =====================================================

# scale by abundance
dt_prod_rare <- dt_prod_rare %>% 
  mutate(across(starts_with('biomass_gain'), ~ Abundance * .x)) %>% 
  select(!Abundance)

# aggregate to sampling event level
dt_prod_rare <- dt_prod_rare %>% ungroup %>% 
  group_by(DataID, seriesID, SamplingEventID) %>% 
  summarise(across(starts_with('biomass_gain'), ~ sum(.x))) %>% ungroup
summary(dt_prod_rare)

# convert to kg/yr
dt_prod <- dt_prod %>% mutate(across(starts_with('biomass_gain'), ~ .x / 1000))

colnames(dt_prod_rare) <- str_replace(colnames(dt_prod_rare), 'biomass\\_gain', 'prod_kg')

dt_meta$DataID <- str_to_lower(dt_meta$DataID)
dt_prod_rare$DataID <- str_to_lower(dt_prod_rare$DataID)

# kg/ha/yr
dt_prod_rare <- dt_meta %>% distinct(DataID, `Grain size`) %>% 
  mutate(Area = str_extract(`Grain size`, '[0-9\\.]+') %>% as.numeric) %>% 
  left_join(dt_prod_rare, ., by = "DataID") %>% 
  select(!`Grain size`)

summary(dt_prod_rare)
# unit conversion
dt_prod_rare <- dt_prod_rare %>% mutate(across(starts_with('prod'), ~ .x * 10000 / Area))


# SINGLE TIME DATA =========================================================================

# give predKmax data on all fishes observed
# then calculate productivity per rareID

# need to keep some sampling info since SST are specific to location + time
# dt_single, because it's full data
dt_single <- dt_sp %>% select(taxaID, ends_with('valid'), taxon_resolution) %>% 
  left_join(dt_single, ., by = "taxaID") # join in
summary(dt_single)

# traits
# use taxaID as link
dt_traits <- dt_sp %>% select(taxaID, fishbase_valid) %>% 
  filter(!is.na(fishbase_valid)) %>% 
  left_join(., dt_traits, by = 'fishbase_valid', relationship = 'many-to-one')

# unique taxa at a given location and time
prodtrait <- dt_single %>% 
  distinct(DataID, seriesID, taxaID, gridcell, SamplingEventID, Year, Month, date_res) %>%
  inner_join(dt_traits, ., by = c('taxaID'))
# automatically filters out the uncertains that don't qualify for trait matching
# inner join somehow didn't quite work
prodtrait <- prodtrait %>% filter(!is.na(Diet))

# SST
prodtrait <- left_join(prodtrait, dt_sst, by = c('gridcell', 'Year', 'Month', 'date_res')) %>% # SST
  select(!c(gridcell, date_res, Year, Month)) %>% # no longer need these space time links
  rename(sstmean = sst, Species = fishbase_valid) # name matching
prodtrait <- prodtrait %>% relocate(ends_with('ID'), .after = sstmean)

# our own version of tidytrait
prodtrait <- prodtrait %>% mutate(
  Diet = factor(Diet, levels = sort(unique(dt_prodDB$Diet))),
  Position = factor(Position, levels = sort(unique(dt_prodDB$Position)))
)

prodtrait <- prodtrait %>% arrange(DataID, SamplingEventID)

dt_prodDB <- dt_prodDB %>% mutate(
  Diet = factor(Diet, levels = sort(unique(dt_prodDB$Diet))),
  Position = factor(Position, levels = sort(unique(dt_prodDB$Position))),
  Method = factor(Method, levels = sort(unique(dt_prodDB$Method)))
)
# make db match

samplelink <- prodtrait %>% select(ends_with('ID'), Species, sstmean)
prodtrait <- prodtrait %>% distinct(pick(Species:sstmean))

## Kmax, growth coef ============================================================

# Formula from Morais and Bellwood (2018)
growthf <- formula(~ sstmean + MaxSizeTL + Diet + Position)

# small tweak to predKmax
predK <- function(traits, dataset, fmod, params = NULL, niter, nrounds = 150, 
                  verbose = 0, print_every = 1000, return = c("pred", "relimp", "models"), lowq = 0.25, uppq = 0.75) {
  if (missing(dataset)) {
    dataset <- db
  }
  modmatnd <- stats::model.matrix(fmod, data = traits)[, -1]
  modmat <- stats::model.matrix(fmod, dataset)[, -1]
  xgbmod <- list()
  rel_imp <- list()
  pred.grid <- as.data.frame(matrix(ncol = niter, nrow = nrow(traits), 
                                    dimnames = list(NULL, paste0("Boot", 1:niter))))
  for (i in 1:niter) {
    set.seed(i)
    xgbmod[[i]] <- xgboost::xgboost(modmat, label = dataset$Kmax, nrounds = nrounds, params = params, verbose = verbose, 
                                    print_every = print_every)
    rel_imp[[i]] <- xgboost::xgb.importance(colnames(modmat), 
                                            model = xgbmod[[i]])
    pred.grid[, i] <- stats::predict(xgbmod[[i]], newdata = modmatnd)
    cat(paste("Bootstrapping the model, round", i), "\n")
  }
  df <- cbind(traits, Kmax = apply(pred.grid, 1, stats::median), 
              Kmax_lowq = apply(pred.grid, 1, stats::quantile, probs = lowq), 
              Kmax_uppq = apply(pred.grid, 1, stats::quantile, probs = uppq))
  ret <- list(pred = df, relimp = rel_imp, models = xgbmod)
  return(ret)
}

pars  <- rfishprod:::xgboostparams # make sure the parameters are loaded
pars$max_depth <- 9

set.seed(24)

temp <- predK(traits = prodtrait, 
              dataset = dt_prodDB, fmod = growthf, params = pars, niter = 1000, 
              lowq = 0.025, uppq = 0.975, verbose = T, return = 'pred')

prod_single <- samplelink %>% distinct(taxaID, DataID, seriesID, SamplingEventID, Species, .keep_all = T) %>% 
  left_join(., temp$pred, by = c('Species', 'sstmean'))


## Gross somatic gain =============================================================================

# Kmax is estimated per species but now integrate with individual length records

dt_prod_single <- dt_single %>% 
  select(DataID, seriesID, taxaID, SamplingEventID, Abundance, Length, Length_lwr, Length_upr)
dt_prod_single <- left_join(dt_prod_single, prod_single, by = c('taxaID', 'DataID', 'seriesID', 'SamplingEventID'), 
                            na_matches = 'na', unmatched = 'error', relationship = 'many-to-one') %>% 
  select(!c(Diet, Position))
summary(dt_prod_single)
# View(dt_prod %>% filter(is.na(Kmax))) # the NAs are from species that didn't qualify for trait/prod

# remove those records and constrain lengths
dt_prod_single <- dt_prod_single %>% filter(!is.na(Kmax)) %>% 
  mutate(Length = case_when(Length > MaxSizeTL ~ MaxSizeTL - 0.2, .default = Length))

summary(dt_prod_single)

# lowers and uppers refer to the Kmax value used!
# not necessarily the lower or upper of that metric
dt_prod_single <- dt_prod_single %>% 
  mutate(biomass_gain = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax),
         biomass_gain_l = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax_lowq),
         biomass_gain_u = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax_uppq)
  )

summary(dt_prod_single)

# dt_prod %>% filter(delta_g_day > 10) %>% View # make sure the large increases are proportional
# dt_prod %>% filter(delta_g_day < 0) %>% View
# don't need to keep traits
dt_prod_single <- dt_prod_single %>% select(!c(Species:sstmean))


## Apply to assemblage data =====================================================

# scale by abundance
dt_prod_single <- dt_prod_single %>% 
  mutate(across(starts_with('biomass'), ~ Abundance * .x)) %>% 
  select(!Abundance)

# aggregate to sampling event level
dt_prod_single <- dt_prod_single %>% ungroup %>% 
  group_by(DataID, seriesID, SamplingEventID) %>% 
  summarise(across(starts_with('biomass'), ~ sum(.x))) %>% ungroup
summary(dt_prod_single)

# convert to kg/yr
dt_prod_single <- dt_prod_single %>% mutate(across(starts_with('biomass'), ~ .x / 1000))

colnames(dt_prod_single) <- str_replace(colnames(dt_prod_single), 'biomass\\_gain', 'prod_kg')

dt_meta$DataID <- str_to_lower(dt_meta$DataID)
dt_prod_single$DataID <- str_to_lower(dt_prod_single$DataID)
# kg/ha/yr
dt_prod_single <- dt_meta %>% distinct(DataID, `Grain size`) %>% 
  mutate(Area = str_extract(`Grain size`, '[0-9\\.]+') %>% as.numeric) %>% 
  left_join(dt_prod_single, ., by = "DataID") %>% 
  select(!`Grain size`)

summary(dt_prod_single)
# unit conversion
dt_prod_single <- dt_prod_single %>% mutate(across(starts_with('prod'), ~ .x * 10000 / Area))

