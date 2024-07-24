

# ProducDiv
# Calculate diversity

require(iNEXT)
require(dplyr)
require(ggplot2)
require(beepr)

# requires 01_gridding for the data source
# iter should also be an inherited parameter
dt_div <- dt_rareass %>% filter(iter >= 1, iter <= 5) %>%
    group_by(iter, rareID, Year) %>%
    # calculate richness for each iteration
    summarise(S = n_distinct(paste0(Family, Genus, Species))) %>% 
    ungroup() %>% group_by(rareID, Year) %>%
    # summarise richness
    summarise(metric = 'S', value = quantile(S, probs = 0.5),
    lwr = quantile(S, probs = 0.025), upr = quantile(S, probs = 0.975))
    

head(dt_div)

# delete this for a full run
iter = 1
# asymptotic species richness
# not quite sure how to integrate iterative dataset in this yet
for (i in 1:iter) { # up to iteration 5 for now
    for (R in 1:n_distinct(dt_rareass$rareID)) {
    rare_id <- unique(dt_rareass$rareID)[R]
    # prep data input for iNEXT
    abund <- dt_rareass %>% filter(rareID == rare_id, iter == i) %>%
        select(Year, Abundance, Family, Genus, Species) %>%
        group_by_at(vars(-Abundance)) %>%
        # collapse length
        summarise(Abundance = sum(Abundance)) %>% ungroup %>% 
        mutate(Taxa = paste(Family, Genus, Species, sep = " "))
    # wide format    
    abund <- abund %>% select(!c(Family:Species)) %>% 
        tidyr::pivot_wider(names_from = Year, values_from = Abundance)
    na <- as.list(rep(0, ncol(abund)-1))
    names(na) <- colnames(abund)[-1]
    abund <- tidyr::replace_na(abund, na)

    # calculate asymptotic q0
    temp <- iNEXT(x = as.data.frame(abund[-1]), q = 0, datatype = 'abundance', se = T)
    # store in a tibble
    if (exists('S_asy') == F) { # assign new
        S_asy <- temp$AsyEst %>% filter(Diversity == 'Species richness') %>%
        select(Estimator, LCL, UCL) %>% 
        mutate(rareID = rare_id, Year = as.numeric(colnames(abund)[-1]), iter = i, metric = 'S_asy') %>%
        rename(value = Estimator, lwr = LCL, upr = UCL)
    } else { # append
        S_asy <- temp$AsyEst %>% filter(Diversity == 'Species richness') %>%
        select(Estimator, LCL, UCL) %>% 
        mutate(rareID = rare_id, Year = as.numeric(colnames(abund)[-1]), iter = i, metric = 'S_asy') %>%
        rename(value = Estimator, lwr = LCL, upr = UCL) %>%
        bind_rows(S_asy, .)
    }
    }
}
beep()
rm(abund, temp, na, rare_id, R)
str(S_asy)
summary(S_asy)
dt_div <- bind_rows(dt_div, S_asy[-6]) # join up with species richness
tail(dt_div)
rm(S_asy)

# visualise
source('https://gist.githubusercontent.com/cherfychow/e9ae890fd16f4c86730748c067feee2b/raw/b2db138ab8164c237129308ea643de78cd54b252/cherulean.R')

f_div <- ggplot(data = dt_div) +
    geom_ribbon(aes(x = Year, ymin = lwr, ymax = upr, fill = rareID),
        alpha = 0.2) +
    geom_line(aes(x = Year, y = value, color = rareID)) +
    guides(fill = 'none', color = 'none') +
    facet_wrap(vars(metric), scales = 'free') +
    scale_color_cherulean(palette = 'cheridis', discrete = T) +
    scale_fill_cherulean(palette = 'cheridis', discrete = T)

plotly::ggplotly(f_div)
