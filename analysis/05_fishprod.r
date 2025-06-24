
# ProducDiv
# Calculate productivity from fish traits + SST with observed lengths

require(tidytable)
require(stringr)
require(rfishprod)
require(readr)

colnames(dt_rareass)
colnames(dt_single)
dt_assemblage <- bind_rows(dt_rareass, dt_single) %>% 
  select(ID, taxaID, Length, Abundance, SamplingEventID, seriesID, DataID, 
         Month, Year, gridcell, date_res) %>% 
  left_join(., dt_sp %>% select(taxaID, fishbase_valid, worms_valid), by = 'taxaID')

# match up traits
# use taxaID as link
dt_traits <- dt_sp %>% select(taxaID, fishbase_valid) %>% 
  filter(!is.na(fishbase_valid)) %>% 
  left_join(., dt_traits, by = 'fishbase_valid', relationship = 'many-to-one')

# gather trait and SST data ========================================================================

# unique taxa at a given location and time
prodtrait <- dt_assemblage %>% 
  distinct(gridcell, Year, Month, date_res, taxaID, fishbase_valid) %>%
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


# Estimate Kmax ====================================================================================


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

# Apply Kmax to assemblage data ====================================================================

# dt_vbgf <- dt_rareass %>% distinct(gridcell, Year, Month, SamplingEventID, seriesID) %>% 
#   left_join(dt_vbgf, ., by = c('SamplingEventID', 'seriesID'))
dt_vbgf <- dt_vbgf %>% select(Species:sstmean, starts_with('Kmax'), gridcell, Year, Month) %>% 
  ungroup %>% distinct
dt_Kmax <- bind_rows(dt_vbgf, prod_single) %>% distinct(Species, sstmean, gridcell, Year, Month, .keep_all = T)

dt_Kmax %>% ungroup %>% group_by(Species) %>% summarise(K_range = diff(range(Kmax))) %>% View

dt_prod_all <- left_join(dt_assemblage, dt_Kmax, by = c('fishbase_valid'='Species', 'gridcell', 'Year', 'Month'))

# calculate somatic gain
dt_prod_all %>% filter(Length >= MaxSizeTL) %>% nrow
dt_prod_all <- dt_prod_all %>% mutate(Length = case_when(Length >= MaxSizeTL ~ MaxSizeTL - 1, .default = Length))
summary(dt_prod_all)

dt_prod_all <- dt_prod_all %>% filter(!is.na(Kmax)) %>% 
mutate(biomass_gain = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax),
       biomass_gain_l = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax_lowq),
       biomass_gain_u = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax_uppq)
)

summary(dt_prod_all)

# scale by abundance
dt_prod_all <- dt_prod_all %>% ungroup %>% 
  mutate(across(starts_with('biomass_gain'), ~ Abundance * .x))


# aggregate to sampling event level
dt_prod_all <- dt_prod_all %>% ungroup %>% select(!Abundance) %>% 
  group_by(DataID, seriesID, SamplingEventID, Year, Month) %>% 
  summarise(biomass_gain = sum(biomass_gain),
            biomass_gain_l = sum(biomass_gain_l),
            biomass_gain_u = sum(biomass_gain_u)) %>% ungroup
summary(dt_prod_all)

# convert to kg/yr
dt_prod_all <- dt_prod_all %>% mutate(across(starts_with('biomass_gain'), ~ .x / 1000))

colnames(dt_prod_all) <- str_replace(colnames(dt_prod_all), 'biomass\\_gain', 'prod_kg')

dt_prod_all$DataID <- str_to_lower(dt_prod_all$DataID)

# kg/ha/yr
dt_prod_all <- dt_meta %>% distinct(DataID, `Grain size`) %>% 
  mutate(Area = str_extract(`Grain size`, '[0-9\\.]+') %>% as.numeric) %>% 
  left_join(dt_prod_all, ., by = "DataID") %>% 
  select(!`Grain size`)

summary(dt_prod_all)
# unit conversion
dt_prod_all <- dt_prod_all %>% mutate(across(starts_with('prod'), ~ .x * 10000 / Area))





