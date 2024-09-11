

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

# for all the sampling events in each rareID
for (R in 1:nrow(dt_div)) {
  
  # prep data input for iNEXT
  # wide format. rows = species. cols = time/space
  

  abund <- dt_rareass %>% filter(seriesID == dt_div$seriesID[R] | rareID == dt_div$seriesID[R]) %>%
    select(Year, SamplingEventID, Abundance, worms_valid) %>%
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
  
  if (inherits(errorhandle, 'error')) {
    # fill in NAs for error fail
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
    next
  }
  
  # store in a tibble
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
  print(paste0("Iteration", R, ' done'))
}
save.image('02_div.RData')
# beep()
# 
# 
# rm(abund, temp, na, rare_id, R)
# str(S_asy)
# summary(S_asy)
# S_asy <- S_asy[-6] %>% 
#   # match up back with years
#   left_join(., dt_div[1:3], by = c('SamplingEventID', 'rareID'))
# S_asy$metric <- 'S_asy'
# 
# dt_div <- bind_rows(dt_div, S_asy)  # join up with species richness
# 
# rm(S_asy)
# 
# # visualise
# source('https://gist.githubusercontent.com/cherfychow/e9ae890fd16f4c86730748c067feee2b/raw/b2db138ab8164c237129308ea643de78cd54b252/cherulean.R')
# 
# f_div <- ggplot(data = dt_div) +
#     geom_errorbar(aes(x = Year, ymax = upr, ymin = lwr, color = rareID)) +
#     geom_point(aes(x = Year, y = value, color = rareID)) +
#     guides(fill = 'none', color = 'none') +
#     facet_wrap(vars(metric), scales = 'free') +
#     scale_color_cherulean(palette = 'cheridis', discrete = T) +
#     scale_fill_cherulean(palette = 'cheridis', discrete = T)
# 
# plotly::ggplotly(f_div)
