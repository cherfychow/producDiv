
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

# just IDs, time, trait taxa name, length
prodtrait <- dt_assemblage %>% 
  distinct(ID, DataID, seriesID, taxaID, gridcell, SamplingEventID, Year, Month, date_res, Length) %>%
  inner_join(dt_traits, ., by = c('taxaID'))
# automatically filters out the uncertains that don't qualify for trait matching
# inner join somehow didn't quite work
prodtrait <- prodtrait %>% filter(!is.na(Length))
prodtrait <- prodtrait %>% filter(!is.na(Diet))

# SST
prodtrait <- left_join(prodtrait, dt_sst, by = c('gridcell', 'Year', 'Month', 'date_res')) %>% # SST
  select(!c(gridcell, date_res, Year, Month)) %>% # no longer need these space time links
  rename(sstmean = sst, Size = Length, Species = fishbase_valid) # name matching
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

# check if any fish are larger than their max size
prodtrait %>% filter(Size > MaxSizeTL) %>% nrow
prodtrait %>% filter(Size > MaxSizeTL) %>% View
# 109k rows

# limit observed size to be below max size
prodtrait <- mutate(prodtrait, 
                    Size = case_when(Size > MaxSizeTL ~ MaxSizeTL - 0.1, 
                                     .default = Size))

# Formula from Morais and Bellwood (2018)
growthf <- formula(~ sstmean + MaxSizeTL + Diet + Position + Method)

# Kmax, growth coef ============================================================

# Kmax estimate per species-gridcell-year-dataID / species-samplingeventid-dataID
# recommended boostrapts on 100s to 1000s scale
set.seed(24)
prod_par <- rfishprod:::xgboostparams # make sure the parameters are loaded
dt_prodDB <- dt_prodDB %>% select(!starts_with('trophic'))
vbgf_par <- predKmax(traits = prodtrait, dataset = dt_prodDB, params = prod_par,
           fmod = growthf, niter = 1000, return = 'pred')

summary(vbgf_par$pred)
rm(prodtrait, growthf)

summary(vbgf_par$pred)

# negative Kmax values
# vbgf_par$pred %>% filter(Kmax <= 0) %>% View

# try rerun
negK <- vbgf_par$pred %>% filter(Kmax <= 0) %>% select(!starts_with('Kmax'))
vbgf_par$pred <- vbgf_par$pred %>% filter(Kmax > 0)
sort(unique(negK$Species))
negK <- predKmax(negK, dataset = dt_prodDB, params = mod_params,
                 fmod = growthf, niter = 1000, return = 'pred')

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


  

