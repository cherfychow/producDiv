

# ProducDiv
# Calculate diversity

require(iNEXT)

# requires 01_gridding for the data source
# iter should also be an inherited parameter
dt_div <- dt_rareass %>% group_by(iter, rareID, Year) %>%
    # calculate richness for each iteration
    summarise(S = n_distinct(paste0(Family, Genus, Species))) %>% 
    ungroup() %>% group_by(rareID, Year) %>%
    # summarise richness
    summarise(metric = 'S', value = quantile(S, probs = 0.5),
    lwr = quantile(S, probs = 0.025), upr = quantile(S, probs = 0.975))
    

head(dt_div)

# asymptotic species richness
# not quite sure how to integrate iterative dataset in this yet
S_asy <- as.list(rep('', n_distinct(dt_rareass$rareID)))
for (i in 1:5) { # up to iteration 5 for now
    for (R in 1:n_distinct(dt_rareass$rareID)) {
    rare_id <- unique(dt_rareass$rareID)[R]
    # prep data input for iNEXT
    temp <- dt_rareass %>% filter(rareID == rare_id, iter == i) %>%
        select(Year, Abundance, Family, Genus, Species) %>%
        group_by_at(vars(-Abundance)) %>%
        # collapse length
        summarise(Abundance = sum(Abundance)) %>% ungroup %>% 
        mutate(Taxa = paste(Family, Genus, Species, sep = " "))
    # wide format    
    temp <- temp %>% select(!c(Family:Species)) %>% 
        tidyr::pivot_wider(names_from = Year, values_from = Abundance)
    na <- as.list(rep(0, ncol(temp)-1))
    names(na) <- colnames(temp)[-1]
    temp <- tidyr::replace_na(temp, na)

    # calculate asymptotic q0
    S_asy[[R]] <- iNEXT(x = as.data.frame(temp[-1]), q = 0, datatype = 'abundance', se = T)
    S_asy[[R]] <- S_asy[[R]]$AsyEst %>% filter(Diversity == 'Species richness') %>%
        select(Estimator, LCL, UCL) %>% 
        mutate(rareID = rare_id, Year = colnames(temp)[-1], iter = i, 
        metric = 'S_asy', .before = Estimator) %>%
        rename(value = Estimator, lwr = LCL, upr = UCL)
    }
}
beep()
rm(temp, na, rare_id, R)
sasy <- bind_rows(S_asy)
sasy$Year <- as.numeric(sasy$Year)
dt_div <- bind_rows(dt_div, sasy[-3]) # join up with species richness
rm(S_asy, sasy)

# visualise
ggplot(data = dt_div) +
    geom_ribbon(aes(x = Year, ymin = lwr, ymax = upr, fill = rareID),
        alpha = 0.3) +
    geom_line(aes(x = Year, y = value, color = rareID)) +
    guides(fill = 'none', color = 'none') +
    facet_wrap(vars(metric))

