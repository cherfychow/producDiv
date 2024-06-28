

# ProducDiv
# Calculate diversity

require(iNEXT)

# requires 01_gridding for the data source
dt_rareass$rareID <- with(dt_rareass, paste(DataID, gridcell))
dt_div <- dt_rareass %>% group_by(rareID, Year) %>%
    summarise(S = n_distinct(paste0(Family, Genus, Species)))
head(dt_div)

S_asy <- as.list(rep('', n_distinct(dt_rareass$rareID)))
for (R in 1:n_distinct(dt_rareass$rareID)) {
    rare_id <- unique(dt_rareass$rareID)[R]
    temp <- dt_rareass %>% filter(rareID == rare_id) %>%
        select(Year, Abundance, Family, Genus, Species) %>%
        group_by_at(vars(-Abundance)) %>%
        # collapse length
        summarise(Abundance = sum(Abundance)) %>% ungroup %>% 
        mutate(Taxa = paste(Family, Genus, Species, sep = " "))
    temp <- temp %>% select(!c(Family:Species)) %>% 
        tidyr::pivot_wider(names_from = Year, values_from = Abundance)
}