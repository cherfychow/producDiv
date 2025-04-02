

# ProducDiv
# Calculate asymptotic diversity and final prep of diversity predictors

require(iNEXT.3D)
require(tidytable)
require(dplyr)
require(ggplot2)

# requires 01_gridding for the data source

# CALCULATE OBSERVED RICHNESS ======================================================================

dt_div <- dt_sp %>% select(taxaID, worms_valid) %>% 
  left_join(dt_rareass, ., by = 'taxaID') %>% 
  group_by(DataID, seriesID, SamplingEventID, Year) %>% 
  summarise(S = n_distinct(worms_valid)) %>% ungroup

# ASYMPTOTIC RICHNESS ==============================================================================
# use dt_div to index since it's already structured by rareID + samplingevents
# time intensive loop

# for all the sampling events in each rareID
# we use timeseries to index by seriesID, so 

samplelink <- dt_rareass %>% ungroup %>% distinct(DataID, seriesID, SamplingEventID, Year) %>% 
  arrange(DataID, seriesID, Year, SamplingEventID)
samplelink$sampleID <- paste0('rare_',1:nrow(samplelink))
dt_rareass <- left_join(dt_rareass, samplelink, by = colnames(samplelink)[-5], relationship = 'many-to-one')

set.seed(24)

for (i in 1:n_distinct(dt_rareass$DataID)) { # run per study
  
  # prep data input for iNEXT
  # wide format. rows = species. cols = time/space
  study <- sort(unique(dt_rareass$DataID))[i]
  abund <- dt_rareass %>% 
    filter(DataID == study) %>% 
    mutate(grainsize = dt_meta$`Grain size`[which(dt_meta$DataID == study)] %>% 
             str_extract(., '\\d{2,}') %>% as.numeric) %>% 
    select(sampleID, Abundance, taxaID, grainsize) %>% ungroup %>% 
    group_by_at(vars(-Abundance)) %>%
    # collapse length
    summarise(Abundance = sum(Abundance*grainsize) %>% round(., 0)) %>% ungroup # density back to counts
  
  abund <- dt_sp %>% select(taxaID, worms_valid) %>% 
    left_join(abund, ., by = 'taxaID') %>% 
    select(!c(taxaID, grainsize)) %>% 
    ungroup %>% 
    group_by(sampleID, worms_valid) %>% 
    summarise(Abundance = sum(Abundance))
  
  
  # wide format
  # drop year. we can match that back from eventID
  abund <- abund %>%
    tidyr::pivot_wider(names_from = sampleID, values_from = Abundance)
  na <- as.list(rep(0, ncol(abund)-1))
  names(na) <- colnames(abund)[-1]
  abund <- tidyr::replace_na(abund, na)
  
  # calculate asymptotic q0
  errorhandle <- tryCatch(temp <- iNEXT.3D::estimate3D(diversity = 'TD',
                                             data = as.data.frame(abund[-1]), 
                                             q = 0, datatype = 'abundance',
                                             base = 'coverage'),
                          error = function(e)e)
  
  # fill in NAs for error fail
  if (inherits(errorhandle, 'error')) {
    if (exists('S_asy') == F) {
      S_asy <- tibble(SamplingEventID = colnames(abund)[-1], 
                      seriesID = dt_div$seriesID[R],
                      metric = 'S_asy',
                      value = NA)
    } else {
      S_asy <- tibble(SamplingEventID = colnames(abund)[-1], 
                      seriesID = dt_div$seriesID[R],
                      metric = 'S_asy',
                      value = NA) %>% bind_rows(S_asy, .)
    }
    print(paste0('Iteration ', R, ' failed'))
    next
  }
  
  # store in a tibble if it worked
  if (exists('S_asy') == F) { # assign new
    S_asy <- temp$AsyEst %>% filter(Diversity == 'Species richness') %>%
      select(Estimator, LCL, UCL) %>% 
      mutate(seriesID = dt_div$seriesID[R], 
             SamplingEventID = colnames(abund)[-1],
             metric = 'S_asy') %>%
      rename(value = Estimator, lwr = LCL, upr = UCL)
  } else { # append
    S_asy <- temp$AsyEst %>% filter(Diversity == 'Species richness') %>%
      select(Estimator, LCL, UCL) %>% 
      mutate(seriesID = dt_div$seriesID[R], 
             SamplingEventID = colnames(abund)[-1],
             metric = 'S_asy') %>%
      rename(value = Estimator, lwr = LCL, upr = UCL) %>%
      bind_rows(S_asy, .)
  }
  print(paste0("Iteration ", R, ' of ', nrow(dt_div), ' done'))
  
  save(S_asy, file = 'expS.RData')
}

rm(abund, temp, na, rare_id, R, errorhandle)
str(S_asy)
summary(S_asy)

# isolate, rerun fix
fix <- S_asy %>% filter(is.na(value)) %>% select(seriesID)
S_asy_clean <- S_asy %>% filter(!seriesID %in% fix$seriesID) # remove from main
# it was because of decimal points and 0.5s. round densities
# fix data won't converge. omit
rm(fix)

# join up with observed species richness
dt_div <- S_asy_clean %>% rename(S_asy = value, asy_lwr = lwr, asy_upr = upr) %>% select(!metric) %>%
  left_join(dt_div,., by = c('seriesID', 'SamplingEventID'))

