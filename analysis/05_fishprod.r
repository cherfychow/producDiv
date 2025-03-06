
# ProducDiv
# Calculate productivity from fish traits + SST with observed lengths

require(tidytable)
require(stringr)
require(rfishprod)

# PREP =========================================================================

# give predKmax data on all fishes observed
# then calculate productivity per rareID

# need to keep some sampling info since SST are specific to location + time
# dt_assemblage, because it's full data
dt_assemblage <- bind_rows(dt_rareass, dt_single)
dt_assemblage <- dt_sp %>% select(taxaID, ends_with('valid'), taxon_resolution) %>% 
  left_join(dt_assemblage, ., by = "taxaID") # join in
summary(dt_assemblage)

# traits
# use taxaID as link
dt_traits <- dt_sp %>% select(taxaID, fishbase_valid) %>% left_join(., dt_traits, by = 'fishbase_valid')

# just IDs, time, trait taxa name
prodtrait <- dt_assemblage %>% 
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
sort(unique(prodtrait$Diet)) == sort(unique(dt_prodDB$Diet))
sort(unique(prodtrait$Position)) == sort(unique(dt_prodDB$Position))
# have to fit the db's factor ordering
prodtrait <- prodtrait %>% 
  mutate(Method = factor('Otolth', levels = sort(unique(dt_prodDB$Method))),
         Diet = factor(Diet, levels = sort(unique(dt_prodDB$Diet))),
         Position = factor(Position, levels = levels(dt_prodDB$Position)))
dt_prodDB <- dt_prodDB %>% 
  mutate(Method = factor(Method, levels = sort(unique(dt_prodDB$Method))),
         Diet = factor(Diet, levels = sort(unique(dt_prodDB$Diet))),
         Position = factor(Position, levels = levels(dt_prodDB$Position)))
summary(dt_prodDB)
summary(prodtrait)
# make db match

# Formula from Morais and Bellwood (2018)
growthf <- formula(~ sstmean + MaxSizeTL + Diet + Position)

# Kmax, growth coef ============================================================

# small tweak to predKmax force pars override
predKmax <- function(traits, dataset, fmod, params = NULL, niter, nrounds = 150, 
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


# Kmax estimate per species-gridcell-year-dataID / species-samplingeventid-dataID
# recommended boostrapts on 100s to 1000s scale
set.seed(24)
prod_par <- rfishprod:::xgboostparams # make sure the parameters are loaded
dt_prodDB <- dt_prodDB %>% select(!starts_with('trophic'))


x = 20000
# nrow(prodtrait)/30000
dt_kmax <- as.list(rep(NA, (ceiling(nrow(prodtrait)/x))))
for (i in 1:length(dt_kmax)) {
  if (i == length(dt_kmax)) { # how to handle the last iteration
    dt_kmax[[i]] <- predK(traits = prodtrait[((i-1)*x+1):nrow(prodtrait),],
      dataset = dt_prodDB, fmod = growthf, params = prod_par, niter = 1000, lowq = 0.025, uppq = 0.975)
  } else {
      dt_kmax[[i]] <- predKmax(traits = prodtrait[((i-1)*x+1):(i*x),],
        dataset = dt_prodDB, fmod = growthf, params = prod_par, niter = 1000, lowq = 0.025, uppq = 0.975)
  }
  print(paste0('Run ', i, ' of ', length(dt_kmax), ' complete'))
}

summary(vbgf_par$pred)
rm(prodtrait, growthf)

# Mortality and gross somatic gain =============================================

# mortality
vbgf_par <- vbgf_par %>% mutate(
  M = predM(Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax, method = 'Lorenzen'),
  M_l = predM(Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax_lowq, method = 'Lorenzen'),
  M_u = predM(Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax_uppq, method = 'Lorenzen'))

vbgf_par <- vbgf_par %>% filter(Kmax > 0) %>% # leave out the nonsensicals for now
  mutate(
    Size_t1 = applyVBGF(Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax), # rate of change day^-1
    size_delta = Size_t1 - Size
  ) # length change

# lowers and uppers refer to the Kmax value used!
# not necessarily the lower or upper of that metric
vbgf_par <- vbgf_par %>% 
  mutate(biomass_gain = somaGain(a, b, Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax),
         biomass_loss = somaLoss(M, Lmeas = Size, a, b),
         delta_g_day = biomass_gain - biomass_loss,
         biomass_gain_l = somaGain(a, b, Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax_lowq),
         biomass_loss_l = somaLoss(M_l, Lmeas = Size, a, b),
         delta_g_day_l = biomass_gain_l - biomass_loss_l,
         biomass_gain_u = somaGain(a, b, Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax_uppq),
         biomass_loss_u = somaLoss(M_u, Lmeas = Size, a, b),
         delta_g_day_u = biomass_gain_u - biomass_loss)

summary(vbgf_par)
# don't need to keep traits
vbgf_par <- vbgf_par %>% mutate(rareID = paste(DataID, gridcell, sep = '_')) %>% 
  select(rareID, SamplingEventID, Species, Size, Kmax:delta_g_day_u)

# Apply to assemblage data =====================================================

colnames(dt_rareass)
dt_rareass <- dt_rareass %>% 
  # known uncertain blind species names for productivity calculations
  mutate(prod_genus_sp = case_when(str_detect(genus_sp, '(?<=idae)\\ssp\\d*$') ~ str_remove(genus_sp, '\\ssp\\d*$'),
                                   .default = str_remove(genus_sp, '\\d+')))

# match up
dt_prod <- vbgf_par %>% select(rareID:Size, starts_with('delta')) %>% 
  left_join(dt_rareass, ., by = c('rareID', 'SamplingEventID', 
                                  'prod_genus_sp' = 'Species',
                                  'Length' = 'Size'))
# aggregate to species level
dt_prod <- dt_prod %>% 
  group_by(pick(!c(Abundance, starts_with('delta')))) %>%
  # multiply by count
  summarise(delta_g_day = Abundance * delta_g_day,
            delta_g_day_l = Abundance * delta_g_day_l,
            delta_g_day_u = Abundance * delta_g_day_u) %>% 
  ungroup %>% group_by(pick(!c(Length, delta_g_day, delta_g_day_l, delta_g_day_u))) %>% 
  # aggregate total length and biomass change for the species
  summarise(delta_g_day = sum(delta_g_day, na.rm = T),
            delta_g_day_l = sum(delta_g_day_l, na.rm = T),
            delta_g_day_u = sum(delta_g_day_u, na.rm = T)) %>% ungroup

# aggregate to site level

dt_prod <- dt_prod %>% ungroup %>% 
  group_by(rareID, Year, SamplingEventID) %>%
  summarise(delta_g_day = sum(delta_g_day, na.rm = T),
            delta_g_day_l = sum(delta_g_day_l, na.rm = T),
            delta_g_day_u = sum(delta_g_day_u, na.rm = T)) %>% 
  # convert to kg/yr
  mutate(delta_g_day = delta_g_day * 365 / 1000,
         delta_g_day_l = delta_g_day_l * 365 / 1000,
         delta_g_day_u = delta_g_day_u * 365 / 1000) %>% ungroup

colnames(dt_prod) <- str_replace(colnames(dt_prod), 'delta\\_g\\_day', 'prod_kg_yr')

# kg/ha/yr
dt_prod <- dt_rareass %>% distinct(rareID, DataID) %>% 
  left_join(dt_meta[c(1,7)], ., by = "DataID") %>% 
  mutate(Area = str_extract(`Grain size`, '\\d{2,}')) %>% 
  select(DataID, rareID, Area) %>% 
  left_join(dt_prod, ., by = "rareID")
dt_prod$Area <- as.numeric(dt_prod$Area)
# unit conversion
dt_prod <- dt_prod %>% mutate(across(starts_with('prod'), ~ .x * 10000 / Area))

## Figures ===================================================================


  

