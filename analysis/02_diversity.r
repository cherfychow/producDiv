

# ProducDiv
# Calculate asymptotic diversity

require(iNEXT)
require(dplyr)
require(ggplot2)
require(beepr)

# requires 01_gridding for the data source

# asymptotic species richness
# use dt_div to index since it's already structured by rareID + samplingevents
# time intensive loop

# for studies where there is no location, make seriesID the same as rareID
dt_rareass$seriesID[is.na(dt_rareass$seriesID)] <- dt_rareass$rareID[is.na(dt_rareass$seriesID)]
timeseries <- dt_rareass %>% distinct(DataID, rareID, seriesID)

# for all the sampling events in each rareID
# we use timeseries to index by seriesID, so 
for (R in 1:nrow(timeseries)) {
  
  # prep data input for iNEXT
  # wide format. rows = species. cols = time/space

  abund <- dt_rareass %>%
    filter(seriesID == timeseries$seriesID[R]) %>%
    select(Year, SamplingEventID, Abundance, worms_valid) %>% ungroup %>% 
    group_by_at(vars(-Abundance)) %>%
    # collapse length
    summarise(Abundance = sum(Abundance)) %>% ungroup
  
  # wide format
  # drop year. we can match that back from eventID
  abund <- abund %>% select(!Year) %>%  
    tidyr::pivot_wider(names_from = SamplingEventID, values_from = Abundance)
  na <- as.list(rep(0, ncol(abund)-1))
  names(na) <- colnames(abund)[-1]
  abund <- tidyr::replace_na(abund, na)
  
  # calculate asymptotic q0
  errorhandle <- tryCatch(temp <- iNEXT(x = as.data.frame(abund[-1]), q = 0, datatype = 'abundance', se = T),
           error = function(e)e)
  
  # fill in NAs for error fail
  if (inherits(errorhandle, 'error')) {
    if (exists('S_asy') == F) {
      S_asy <- tibble(SamplingEventID = colnames(abund)[-1], 
                      seriesID = timeseries$seriesID[R],
                      metric = 'S_asy',
                      value = NA)
    } else {
      S_asy <- tibble(SamplingEventID = colnames(abund)[-1], 
                      seriesID = timeseries$seriesID[R],
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
      mutate(seriesID = timeseries$seriesID[R], 
             SamplingEventID = colnames(abund)[-1],
             metric = 'S_asy') %>%
      rename(value = Estimator, lwr = LCL, upr = UCL)
  } else { # append
    S_asy <- temp$AsyEst %>% filter(Diversity == 'Species richness') %>%
      select(Estimator, LCL, UCL) %>% 
      mutate(seriesID = timeseries$seriesID[R], 
             SamplingEventID = colnames(abund)[-1],
             metric = 'S_asy') %>%
      rename(value = Estimator, lwr = LCL, upr = UCL) %>%
      bind_rows(S_asy, .)
  }
  print(paste0("Iteration ", R, ' of ', nrow(timeseries), ' done'))
  
  save(S_asy, file = 'asydiv.RData')
}

beep()


rm(abund, temp, na, rare_id, R)
str(S_asy)
summary(S_asy)
# join up with observed species richness
S_asy <- S_asy %>% rename(S_asy = value, asy_lwr = lwr, asy_upr = upr) %>% select(!metric) %>% 
  left_join(., timeseries[-1], by = "seriesID")
S_asy <- dt_rareass %>% distinct(rareID, seriesID, SamplingEventID, Year) %>% 
  left_join(., S_asy, by = c('rareID', 'seriesID', 'SamplingEventID'))
summary(S_asy)

# some sampling events with NAs in asymptotic richness? try again
missing <- S_asy %>% filter(is.na(S_asy)) %>% distinct(SamplingEventID) %>% pull

missing <- dt_rareass %>%
  filter(SamplingEventID %in% missing) %>%
  select(SamplingEventID, Abundance, worms_valid) %>% ungroup %>% 
  group_by_at(vars(-Abundance)) %>%
  # collapse length
  summarise(Abundance = sum(Abundance)) %>% ungroup

# wide format
# drop year. we can match that back from eventID
missing <- missing %>% 
  tidyr::pivot_wider(names_from = SamplingEventID, values_from = Abundance)
na <- as.list(rep(0, ncol(missing)-1))
names(na) <- colnames(missing)[-1]
missing <- tidyr::replace_na(missing, na)
missing_asy <- iNEXT(x = as.data.frame(missing[-1]), q = 0, datatype = 'incidence_freq', se = T)


View(S_asy, )
rm(S_asy)

# visualise
source('https://gist.githubusercontent.com/cherfychow/e9ae890fd16f4c86730748c067feee2b/raw/b2db138ab8164c237129308ea643de78cd54b252/cherulean.R')

f_div <- ggplot(data = dt_div) +
    geom_errorbar(aes(x = Year, ymax = upr, ymin = lwr, color = rareID)) +
    geom_point(aes(x = Year, y = value, color = rareID)) +
    guides(fill = 'none', color = 'none') +
    facet_wrap(vars(metric), scales = 'free') +
    scale_color_cherulean(palette = 'cheridis', discrete = T) +
    scale_fill_cherulean(palette = 'cheridis', discrete = T)

plotly::ggplotly(f_div)
