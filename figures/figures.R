
# ProducDiv
# Figure code

require(ggplot2)
require(patchwork)
require(sf)
require(rnaturalearth)
require(ggdist)

# colour palettes
source('https://gist.githubusercontent.com/cherfychow/e9ae890fd16f4c86730748c067feee2b/raw/b2db138ab8164c237129308ea643de78cd54b252/cherulean.R')
style <- theme_classic(base_size = 15) + theme(panel.grid = element_blank(),
                                          axis.ticks = element_line(linewidth = 0.3),
                                          axis.line = element_line(linewidth = 0.3))


# Dataset map ---------------------------------------------------------------------------------

proj <- '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'

world <- ne_coastline('medium', returnclass = 'sf') %>%
  st_make_valid(.) %>% st_break_antimeridian(., lon_0 = 160)

# make a geo point object with the data
point_rare <- dt_rareass %>% group_by(DataID, rareID) %>% 
  mutate(Duration = max(Year) - min(Year) + 1) %>%
  distinct(Longitude, Latitude, DataID, rareID, Year, Duration) %>%
  filter(is.na(Latitude) == F) %>%
  st_as_sf(., crs = 4326, coords = c("Longitude", "Latitude")) %>%
  st_make_valid(.) %>%
  st_cast(., to = "MULTIPOINT")

point_single <- dt_single %>% group_by(DataID) %>% 
  distinct(Longitude, Latitude, DataID) %>%
  filter(is.na(Latitude) == F) %>%
  st_as_sf(., crs = 4326, coords = c("Longitude", "Latitude")) %>%
  st_make_valid %>% st_cast(., to = "MULTIPOINT")

# by dataset ID
f_map_time <- ggplot() +
  geom_sf(data = world, color='grey75', fill='grey95') +
  geom_sf(data = point_rare, 
          aes(fill = DataID, size = Duration), shape = 21) +
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) +
  theme_minimal() + guides(fill = 'none') +
  scale_fill_cherulean(palette = 'cheridis', discrete = T, alpha = 0.4) +
  scale_size_continuous(breaks = c(2,5,10,15,20))
f_map_single <- ggplot() +
  geom_sf(data = world, color='grey75', fill='grey95') +
  geom_sf(data = point_single, 
          aes(fill = DataID), shape = 21, size = 2) +
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) +
  theme_minimal() + guides(fill = 'none') +
  scale_fill_cherulean(palette = 'cheridis', discrete = T, alpha = 0.4)
(f_map_time + theme(axis.text.x = element_blank())) / f_map_single

# all together
points <- dt_assemblage %>% group_by(DataID) %>% distinct(Longitude, Latitude, DataID) %>% 
  filter(is.na(Latitude) == F) %>% 
  st_as_sf(., crs = 4326, coords = c("Longitude", "Latitude")) %>%
  st_make_valid(.) %>%
  st_cast(., to = "MULTIPOINT")

f_map <- ggplot() +
  geom_sf(data = world, color='grey75', fill='grey95') +
  geom_sf(data = points, 
          aes(fill = DataID), shape = 21, size = 2.5) +
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) +
  theme_minimal(base_size = 13, base_family = 'Helvetica') +
  guides(fill = 'none') +
  scale_fill_cherulean(palette = 'cheridis', discrete = T, alpha = 0.4)

f_map_eco <- ggplot() +
  geom_sf(data = world, color='grey75', fill='grey95') +
  geom_sf(data = ecoR, aes(fill = region)) +
  geom_sf(data = points, shape = 21, size = 1.5, fill = 'white') +
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) +
  theme_minimal(base_size = 13, base_family = 'Helvetica') +
  theme(legend.position = 'bottom') +
  scale_fill_cherulean(palette = 'cheridis', discrete = T)

f_map_eco / f_map

# NULL MODEL  ======================================================================================

# posterior distributions

f_null_post1 <- ggplot(null_draws %>% filter(var_type == 'fixed')) +
  stat_slab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  scale_fill_manual(values = '#DBC539', na.value = 'grey90') +
  stat_spike(aes(x = values, y = prettyvar, ), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_null_post2 <- ggplot(null_draws %>% filter(var_type == 'randommodel')) +
  stat_slab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  scale_fill_manual(values = '#DBC539', na.value = 'grey90') +
  stat_spike(aes(x = values, y = prettyvar), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_null_post <- f_null_post1 / f_null_post2 + plot_layout(heights = c(2, 5))
f_null_post

f_null <- ggplot(data = vars_single) +
  labs(x = 'S', y = bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  geom_point(aes(x = S, y = prod_null), size = 1, shape = 21, 
             colour = '#666', alpha = 0.2) +
  
  # geom_point(aes(x = S, y = prod_kg, colour = biogeo, fill = biogeo), 
  #            shape = 21, size = 2, alpha = 0.5) +
  # scale_colour_cherulean(palette = 'cheridis', discrete = T) +
  # scale_fill_cherulean(palette = 'cheridis', discrete = T) +
  # 
  # geom_ribbon(data = spatial_partialpred %>% ungroup %>% group_by(S) %>% 
  #               summarise(lwr = quantile(prod_S, 0.025),
  #                         upr = quantile(prod_S, 0.975),
  #                         prod_S = quantile(prod_S, 0.5)), 
  #             aes(x = S, ymin = lwr, ymax = upr), 
  #             fill = 'transparent', colour = '#000', linetype = 'dashed') +
  # 
  # geom_line(data = spatial_partialpred %>% ungroup %>% group_by(S) %>% 
  #             summarise(lwr = quantile(prod_S, 0.025),
  #                       upr = quantile(prod_S, 0.975),
  #                       prod_S = quantile(prod_S, 0.5)), 
  #           aes(x = S, y = prod_S), 
  #           colour = '#000', linewidth = 0.75) +
  
  geom_ribbon(data = null_partialpred %>% ungroup %>% group_by(S) %>%
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)),
              aes(x = S, ymin = lwr, ymax = upr),
              fill = 'transparent', colour = '#000', linetype = 'dashed') +

  geom_line(data = null_partialpred %>% ungroup %>% group_by(S) %>%
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)),
            aes(x = S, y = prod_S),
            colour = '#000', linewidth = 0.75) +
  scale_x_log10() + scale_y_log10() + style +
  guides(colour = 'none', fill = 'none')
f_null

# Spatial single time point model ==================================================================

## Posterior distributions =========================================================================

f_space_post1 <- ggplot(spatial_draws %>% filter(var_type == 'fixed')) +
  stat_slab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  scale_fill_manual(values = '#DBC539', na.value = 'grey90') +
  stat_spike(aes(x = values, y = prettyvar, ), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_space_post2 <- ggplot(spatial_draws %>% filter(var_type == 'randommodel')) +
  stat_slab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  scale_fill_manual(values = '#DBC539', na.value = 'grey90') +
  stat_spike(aes(x = values, y = prettyvar), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_space_post <- f_space_post1 / f_space_post2 + plot_layout(heights = c(2, 5))
f_space_post

## Partial predictions =============================================================================

f_space_S <- ggplot(data = vars_single) +
  labs(x = 'S', y = bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  geom_linerange(aes(x = S, ymin = prod_kg_l, ymax = prod_kg_u),
                linewidth = 0.2, colour = '#888') +
  
  geom_point(aes(x = S, y = prod_kg, colour = biogeo, fill = biogeo), 
             shape = 21, size = 2, alpha = 0.5) +
  scale_colour_cherulean(palette = 'cheridis', discrete = T) +
  scale_fill_cherulean(palette = 'cheridis', discrete = T) +
  
  geom_ribbon(data = spatial_partialpred %>% ungroup %>% group_by(S) %>% 
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)), 
              aes(x = S, ymin = lwr, ymax = upr), 
              fill = 'transparent', colour = '#000', linetype = 'dashed') +
  
  geom_line(data = spatial_partialpred %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = '#000', linewidth = 1) +
  
  # null model lines
  geom_ribbon(data = nullsp1_partial %>% ungroup %>% group_by(S) %>% 
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)), 
              aes(x = S, ymin = lwr, ymax = upr), 
              fill = '#000', alpha = 0.1) +
  
  geom_line(data = nullsp1_partial %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = '#000', alpha = 0.5, linewidth = 0.75) +
  
  geom_ribbon(data = nullsp2_partial %>% ungroup %>% group_by(S) %>% 
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)), 
              aes(x = S, ymin = lwr, ymax = upr), 
              fill = '#000', alpha = 0.1) +
  
  geom_line(data = nullsp2_partial %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = '#000', alpha = 0.5, linewidth = 0.75) +
  
  # null productivity
  geom_ribbon(data = null_partialpred %>% ungroup %>% group_by(S) %>%
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)),
              aes(x = S, ymin = lwr, ymax = upr),
              fill = '#444444', alpha = 0.2) +
  
  geom_line(data = null_partialpred %>% ungroup %>% group_by(S) %>%
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)),
            aes(x = S, y = prod_S),
            colour = 'grey30', alpha = 0.5, linewidth = 0.75) +
  
  scale_x_continuous(breaks = c(1,5,10,50,100,300,500)) + 
  scale_y_continuous(breaks = c(0,0.001,.01, .1, 1, 10, 50, 100)) +
  coord_trans(y = "log10", x = 'log10') + style +
  guides(colour = 'none', fill = 'none')

# with human gravity
f_space_g <- ggplot(data = vars_single) +
  labs(x = 'Human gravity', y =  bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  geom_errorbar(aes(x = hgrav2 + 1, ymin = prod_kg_l, ymax = prod_kg_u),
                linewidth = 0.2, colour = '#666') +
  geom_point(aes(x = hgrav2 + 1, y = prod_kg, colour = biogeo, fill = biogeo), 
             shape = 21, size = 2, alpha = 0.5) +
  scale_colour_cherulean(palette = 'cheridis', discrete = T) +
  scale_fill_cherulean(palette = 'cheridis', discrete = T) +
  
  geom_ribbon(data = spatial_partialpred %>% ungroup %>% group_by(hgrav2) %>% 
                summarise(lwr = quantile(prod_hgrav, 0.025),
                          upr = quantile(prod_hgrav, 0.975),
                          prod_hgrav = quantile(prod_hgrav, 0.5)), 
              aes(x = hgrav2 + 1, ymin = lwr, ymax = upr), 
              fill = 'transparent', colour = '#000', linetype = 'dashed') +
  geom_line(data = spatial_partialpred %>% ungroup %>% group_by(hgrav2) %>% 
              summarise(lwr = quantile(prod_hgrav, 0.025),
                        upr = quantile(prod_hgrav, 0.975),
                        prod_hgrav = quantile(prod_hgrav, 0.5)), 
            aes(x = hgrav2 + 1, y = prod_hgrav), 
            colour = '#333660', linewidth = 0.75) +
  
  scale_y_continuous(breaks = c(0,0.001,.01, .1, 1, 10, 50, 100)) +
  coord_trans(y = "log10", x = 'log10') + style +
  guides(colour = 'none', fill = 'none')

(f_space_S/f_space_g)

# Biogeographic single time point model ============================================================

## Posterior distributions =========================================================================

f_geo_post1 <- ggplot(geo_draws %>% filter(var_type == 'main')) +
  stat_slab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95, colour = '#666', linewidth = 0.2) +  
  scale_fill_manual(values = '#BBB', na.value = '#FFF') +
  stat_spike(aes(x = values, y = prettyvar, ), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_geo_post2 <- ggplot(geo_draws %>% filter(var_type == 'randommodel')) +
  stat_slab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95, colour = '#666', linewidth = 0.2) +  
  scale_fill_manual(values = '#BBB', na.value = '#FFF') +
  stat_spike(aes(x = values, y = prettyvar), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_geo_post <- f_geo_post1 / f_geo_post2 + plot_layout(heights = c(2, 5))
f_geo_post

## Partial predictions =============================================================================

# split by bioregion
f_geo <- ggplot(data = vars_single) +
  labs(x = 'S', y = bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  geom_point(aes(x = S, y = prod_kg, colour = biogeo), 
             shape = 21, size = 1.5, alpha = 0.5, fill = 'transparent') +
  
  geom_ribbon(data = geo_partialpred %>% ungroup %>% group_by(biogeo, S) %>% 
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)), 
              aes(x = S, ymin = lwr, ymax = upr, group = biogeo), 
              fill = 'transparent', colour = '#000', linetype = 'dashed') +
  
  # overall S effect line
  geom_line(data = spatial_partialpred %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = '#000', linewidth = 0.75, alpha = 0.4) +
  
  geom_line(data = geo_partialpred %>% ungroup %>% group_by(biogeo, S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S, group = biogeo), 
            colour = '#000', linewidth = 0.75) +

  scale_colour_cherulean(palette = 'cheridis', discrete = T) +
  scale_fill_cherulean(palette = 'cheridis', discrete = T) +
  scale_x_log10() + scale_y_log10() +
  guides(colour = 'none', fill = 'none') +
  style + facet_wrap(~biogeo)

ggplot(data = vars_single) +
  labs(x = 'log(S)', y = bquote("log Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  
  geom_line(data = geo_partialpred %>% ungroup %>% group_by(biogeo, S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S, colour = biogeo), linewidth = 1.2, alpha = 0.3) +
  
  scale_x_log10() + scale_y_log10() + style



## Residuals =======================================================================================

geo_resid1 <- ggplot(data = spatial_resid) + style +
  labs(x = 'Data', y = 'Model predicted') +
  geom_point(aes(x = prod_kg, y = prod_S), 
             shape = 21, colour = '#000', alpha = 0.2) +
  scale_x_log10() + scale_y_log10()

geo_resid2 <- ggplot(data = spatial_resid) + style +
  labs(x = 'S', y = 'Model residuals') +
  geom_point(aes(x = S, y = resid), 
             shape = 21, colour = '#000', alpha = 0.2) +
  scale_x_log10()

geo_resid1 / geo_resid2

# ∆Prod ~ ∆S, long term trajectory =================================================================

## Posterior distributions =========================================================================

f_T_post1 <- ggplot(T_draws %>% filter(var_type == 'fixed')) +
  stat_slab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  scale_fill_manual(values = '#DBC539', na.value = 'grey90') +
  stat_spike(aes(x = values, y = prettyvar, ), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_T_post2 <- ggplot(T_draws %>% filter(var_type == 'randommodel')) +
  stat_slab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  scale_fill_manual(values = '#DBC539', na.value = 'grey90') +
  stat_spike(aes(x = values, y = prettyvar), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_T_post <- f_T_post1 / f_T_post2 + plot_layout(heights = c(2, 5))
f_T_post


## Partial prediction plot =========================================================================

vars_time_long <- vars_time_long %>% 
  mutate(quad = case_when(.default = 'neg',
                          S_slope > 0 & prod_slope > 0 ~ 'pos',
                          S_slope <= 0 & prod_slope <= 0 ~ 'pos')) %>% 
  mutate(quad = case_when(.default = quad,
                          S_slope > 0 & prod_slope > 0 & S_slope - 1.96*S_se > 0 & prod_slope - 1.96*prod_se > 0 ~ 'pos_certain',
                          S_slope <= 0 & S_slope + 1.96*S_se <= 0 & prod_slope <= 0 & prod_slope + 1.96*prod_se <= 0  ~ 'pos_certain',
                          S_slope > 0 & S_slope - 1.96*S_se > 0 &  prod_slope <= 0 & prod_slope + 1.96*prod_se <= 0 ~ 'neg_certain',
                          S_slope <= 0 & S_slope + 1.96*S_se <= 0 & prod_slope > 0 & prod_slope - 1.96*prod_se > 0 ~ 'neg_certain'))

f_T_S1 <- ggplot(data = vars_time_long) +
  labs(x = "∆S", y = bquote("∆Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = '#AAA') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#AAA') +
  geom_linerange(aes(xmin = S_slope - 1.96*S_se, 
                     xmax = S_slope + 1.96*S_se, 
                     y = prod_slope, colour = quad), size = 0.2, alpha = 0.3) +
  geom_linerange(aes(ymin = prod_slope - 1.96*prod_se, 
                     ymax = prod_slope + 1.96*prod_se, 
                     x = S_slope, colour = quad), size = 0.2, alpha = 0.3) +
  scale_colour_manual(values = c('#455EBA','#455EBA','#DBC539','#DBC539')) +
  # uncertain points
  geom_point(data = vars_time_long %>% filter(!str_detect(quad, 'certain')),
             aes(x = S_slope, y = prod_slope), 
             size = 2, shape = 21, alpha = 0.3, colour = '#444', fill = '#FFF') +
  # certain slopes
  geom_point(data = vars_time_long %>% filter(str_detect(quad, 'certain')),
             aes(x = S_slope, y = prod_slope, fill = quad), 
             size = 2, shape = 21, colour = '#444444AA') +
  scale_fill_manual(values = c('#455EBA','#DBC539')) +
  geom_ribbon(data = T_partialpred %>% ungroup %>% group_by(S_slope) %>% 
                    summarise(lwr = quantile(deltaProd_S, 0.025),
                              upr = quantile(deltaProd_S, 0.975),
                              deltaProd_S = quantile(deltaProd_S, 0.5)), 
                  aes(x = S_slope, ymin = lwr, ymax = upr), 
                  colour = '#000', fill = 'transparent', linetype = 'dashed') +
  geom_line(data = T_partialpred %>% ungroup %>% group_by(S_slope) %>% 
              summarise(lwr = quantile(deltaProd_S, 0.025),
                        upr = quantile(deltaProd_S, 0.975),
                        deltaProd_S = quantile(deltaProd_S, 0.5)), 
            aes(x = S_slope, y = deltaProd_S), 
            colour = '#000', linewidth = 0.5) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 'dashed', colour = '#444', alpha = 0.6) +
  guides(colour = 'none', fill = 'none') +
  scale_x_continuous(limits = c(-2.5, 1.5)) +
  scale_y_continuous(limits = c(-3.5, 3)) + style

f_T_Sx <- ggplot(data = vars_time_long) +
  geom_density(aes(x = S_slope), fill = '#CCC', colour = '#666') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#666') +
  theme_void() + scale_x_continuous(limits = c(-2.5, 1.5))

f_T_Sy <- ggplot(data = vars_time_long) +
  geom_density(aes(y = prod_slope), fill = '#CCC', colour = '#666') +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = '#666') +
  theme_void() + scale_y_continuous(limits = c(-3.5, 3))

f_T_S <- f_T_Sx + plot_spacer() + f_T_S1 + f_T_Sy + 
  plot_layout(ncol = 2, nrow = 2, widths = c(5,1), heights = c(1,6))
f_T_S

f_T_g <- ggplot(data = vars_time_long) +
  labs(x = "Human gravity", y = bquote("∆Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  geom_linerange(aes(ymin = prod_slope - 1.96*prod_se, 
                     ymax = prod_slope + 1.96*prod_se, 
                     x = hgrav2 + 1, colour = quad), size = 0.2) +
  scale_colour_manual(values = c('#455EBA','#455EBA','#DBC539','#DBC539')) +
  geom_point(aes(x = hgrav2 + 1, y = prod_slope, fill = quad), 
             size = 2, shape = 21, colour = '#444') +
  scale_fill_manual(values = c('#455EBA','#455EBA','#DBC539','#DBC539')) +
  
  geom_ribbon(data = T_partialpred %>% ungroup %>% group_by(hgrav2) %>% 
                summarise(lwr = quantile(deltaProd_hgrav, 0.025),
                          upr = quantile(deltaProd_hgrav, 0.975),
                          deltaProd_hgrav = quantile(deltaProd_hgrav, 0.5)), 
              aes(x = hgrav2 + 1, ymin = lwr, ymax = upr), 
              colour = '#000', fill = 'transparent', linetype = 'dashed') +
  geom_line(data = T_partialpred %>% ungroup %>% group_by(hgrav2) %>% 
              summarise(lwr = quantile(deltaProd_hgrav, 0.025),
                        upr = quantile(deltaProd_hgrav, 0.975),
                        deltaProd_hgrav = quantile(deltaProd_hgrav, 0.5)), 
            aes(x = hgrav2 + 1, y = deltaProd_hgrav), 
            colour = '#000', linewidth = 0.5) +
  
  guides(fill = 'none', colour = 'none') +
  scale_x_log10() + style
f_T_g


# ∆Prod ~ ∆S, log ratio adjacent years =============================================================

f_t_post1 <- ggplot() +
  # log ratio fluxes
  stat_slab(data = t_draws %>% filter(var_type == 'fixed'),
            aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  scale_fill_manual(values = '#DBC539', na.value = 'grey90') +
  stat_spike(data = t_draws %>% filter(var_type == 'fixed'),
             aes(x = values, y = prettyvar, ), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  # full time series
  stat_slab(data = T_draws %>% filter(var_type == 'fixed'),
            aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  stat_spike(data = T_draws %>% filter(var_type == 'fixed'),
             aes(x = values, y = prettyvar, ), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_t_post2 <- ggplot() +
  stat_slab(data = t_draws %>% filter(var_type == 'randommodel'),
            aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  scale_fill_manual(values = '#DBC539', na.value = 'grey90') +
  stat_spike(data = t_draws %>% filter(var_type == 'randommodel'),
             aes(x = values, y = prettyvar), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  
  # full time series
  stat_slab(data = T_draws %>% filter(var_type == 'randommodel'),
            aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  stat_spike(data = T_draws %>% filter(var_type == 'randommodel'),
             aes(x = values, y = prettyvar, ), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_t_post <- f_t_post1 / f_t_post2 + plot_layout(heights = c(1,2))


## Partial prediction plot =========================================================================
# 
# vars_time <- vars_time %>% 
#   mutate(quad = case_when(.default = 'neg',
#                           logratio1_S > 0 & logratio1_prod > 0 ~ 'pos',
#                           logratio1_S <= 0 & logratio1_prod <= 0 ~ 'pos')) %>% 
#   mutate(quad = case_when(.default = quad,
#                           logratio1_S > 0 & logratio1_prod > 0 & logratio1_S_l > 0 & logratio1_prod_l > 0 ~ 'pos_certain',
#                           logratio1_S <= 0 & logratio1_S_u <= 0 & logratio1_prod <= 0 & logratio1_prod_u <= 0  ~ 'pos_certain',
#                           logratio1_S > 0 & logratio1_S_l > 0 &  logratio1_prod <= 0 & logratio1_prod_u <= 0 ~ 'neg_certain',
#                           logratio1_S <= 0 & logratio1_S_u <= 0 & logratio1_prod > 0 & logratio1_prod_l > 0 ~ 'neg_certain'))

f_t_S1 <- ggplot(data = vars_time) +
  labs(x = bquote("∆log"[2]~"S"), y = bquote("∆log"[2]~"Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = '#AAA') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#AAA') +
  geom_linerange(aes(xmin = logratio1_S_l, 
                     xmax = logratio1_S_u, 
                     y = logratio1_prod, colour = quad), size = 0.2, alpha = 0.3) +
  geom_linerange(aes(ymin = logratio1_prod_l, 
                     ymax = logratio1_prod_u, 
                     x = logratio1_S, colour = quad), size = 0.2, alpha = 0.3) +
  scale_colour_manual(values = c('#455EBA','#455EBA','#DBC539','#DBC539')) +
  # uncertain points
  geom_point(data = vars_time %>% filter(!str_detect(quad, 'certain')),
             aes(x = logratio1_S, y = logratio1_prod), 
             size = 2, shape = 21, alpha = 0.3, colour = '#444', fill = '#FFF') +
  # certain slopes
  geom_point(data = vars_time %>% filter(str_detect(quad, 'certain')),
             aes(x = logratio1_S, y = logratio1_prod, fill = quad), 
             size = 2, shape = 21, colour = '#444444AA') +
  scale_fill_manual(values = c('#455EBA','#DBC539')) +
  geom_ribbon(data = t_partialpred %>% ungroup %>% group_by(logratio1_S) %>% 
                summarise(lwr = quantile(deltaProd_S, 0.025),
                          upr = quantile(deltaProd_S, 0.975),
                          deltaProd_S = quantile(deltaProd_S, 0.5)), 
              aes(x = logratio1_S, ymin = lwr, ymax = upr), 
              colour = '#000', fill = 'transparent', linetype = 'dashed') +
  geom_line(data = t_partialpred %>% ungroup %>% group_by(logratio1_S) %>% 
              summarise(lwr = quantile(deltaProd_S, 0.025),
                        upr = quantile(deltaProd_S, 0.975),
                        deltaProd_S = quantile(deltaProd_S, 0.5)), 
            aes(x = logratio1_S, y = deltaProd_S), 
            colour = '#000', linewidth = .75) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 'dashed', colour = '#444', alpha = 0.6) +
  guides(colour = 'none', fill = 'none') +
  coord_cartesian(ylim = c(-4,7.5), xlim = c(-3.2, 4)) + style

f_t_Sx <- ggplot(data = vars_time) +
  geom_density(aes(x = logratio1_S), fill = '#CCC', colour = '#666') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#666') +
  theme_void() + scale_x_continuous(limits = c(-3.2, 4))

f_t_Sy <- ggplot(data = vars_time) +
  geom_density(aes(y = logratio1_prod), fill = '#CCC', colour = '#666') +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = '#666') +
  theme_void() + scale_y_continuous(limits = c(-4, 7.5))

f_t_S <- f_t_Sx + plot_spacer() + f_t_S1 + f_t_Sy + 
  plot_layout(ncol = 2, nrow = 2, widths = c(5,1), heights = c(1,6))
f_t_S

f_t_g <- ggplot(data = vars_time) +
  labs(x = "Human gravity",
       y = bquote("∆log"[2]~"Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  geom_linerange(aes(ymin = logratio1_prod - 1.96*prod_se, 
                     ymax = logratio1_prod + 1.96*prod_se, 
                     x = hgrav2 + 1, colour = quad), size = 0.2) +
  scale_colour_manual(values = c('#455EBA','#455EBA','#DBC539','#DBC539')) +
  geom_point(aes(x = hgrav2 + 1, y = logratio1_prod, fill = quad), 
             size = 2, shape = 21, colour = '#444') +
  scale_fill_manual(values = c('#455EBA','#455EBA','#DBC539','#DBC539')) +
  geom_ribbon(data = t2_pred, aes(ymin = prod_gpar_l, ymax = prod_gpar_u, x = hgrav2),
              fill = '#333660', alpha = 0.3) +
  geom_line(data = t2_pred, aes(x = hgrav2, y = prod_gpar), colour = '#333660', linewidth = 0.5) +
  scale_x_log10() + style
f_t_g

# proportions in positive-negative for long and short term fluxes

dt_tchange <- vars_time_long %>% group_by(quad) %>% 
  summarise(n = n()) %>% mutate(group = 'full')
dt_tchange <- vars_time %>% group_by(quad) %>% 
  summarise(n = n()) %>% mutate(group = 'adjacent') %>% 
  bind_rows(dt_tchange, .)

ggplot(data = dt_tchange) +
  geom_col(aes(y = n, fill = quad, x = group), position = 'fill') + 
  scale_fill_manual(values = c('#9CB2FF','#455EBA','#fff29b','#DBC539'), name = NULL,
                    labels = c('negative', 'negative 95% CI', 'positive', 'positive 95% CI')) +
  theme_minimal() + labs(x = NULL, y = NULL) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(panel.grid = element_blank(), legend.position = 'bottom')

# Spatiotemporal figure ============================================================================

# get gridcell centroid as a point

grid_sf <- grid_sf %>% mutate(centroid = st_centroid(geometry))
grid_sf$cell_lat <- st_coordinates(grid_sf$centroid)[,2]
grid_sf$cell_lon <- st_coordinates(grid_sf$centroid)[,1]
vars_time_long <- grid_sf %>% select(gridcell, centroid) %>% 
  left_join(vars_time_long, ., by = 'gridcell')

ggplot() +
  geom_sf(data = world, color='grey75', fill='grey95') +
  geom_sf(data = vars_time_long, aes(geometry = centroid, colour = quad),
              alpha = 0.3) +
  scale_colour_manual(values = c('#455EBA','#455EBA','#DBC539','#DBC539')) +
  theme_minimal(base_size = 13) +
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76))


# NULL SPECIES MODEL ===============================================================================

## Posterior distributions =========================================================================

f_nullsp1_post1 <- ggplot(nullsp1_draws %>% filter(var_type == 'fixed')) +
  stat_slab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  scale_fill_manual(values = '#DBC539', na.value = 'grey90') +
  stat_spike(aes(x = values, y = prettyvar, ), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_nullsp1_post2 <- ggplot(nullsp1_draws %>% filter(var_type == 'randommodel')) +
  stat_slab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  scale_fill_manual(values = '#DBC539', na.value = 'grey90') +
  stat_spike(aes(x = values, y = prettyvar), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_nullsp1_post <- f_nullsp1_post1 / f_nullsp1_post2 + plot_layout(heights = c(2, 5))
f_nullsp1_post

## Partial predictions =============================================================================

f_nullsp <- ggplot(data = vars_single) +
  labs(x = 'S', y = bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  
  geom_point(aes(x = S, y = prod_maxsp), 
             shape = 21, alpha = 0.2, colour = 'grey60') +
  
  geom_point(aes(x = S, y = prod_minsp), 
             shape = 21, alpha = 0.2, colour = 'grey80') +
  
  geom_ribbon(data = nullsp1_partial %>% ungroup %>% group_by(S) %>% 
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)), 
              aes(x = S, ymin = lwr, ymax = upr), 
              fill = 'transparent', colour = 'darkred', linetype = 'dashed') +
  
  geom_line(data = nullsp1_partial %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = 'darkred', linewidth = 0.75) +
  
  geom_ribbon(data = nullsp2_partial %>% ungroup %>% group_by(S) %>% 
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)), 
              aes(x = S, ymin = lwr, ymax = upr), 
              fill = 'transparent', colour = 'darkred', linetype = 'dashed') +
  
  geom_line(data = nullsp2_partial %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = 'darkred', linewidth = 0.75) +
  
  # geom_ribbon(data = geo_partialpred %>% ungroup %>% group_by(biogeo, S) %>% 
  #               summarise(lwr = quantile(prod_S, 0.025),
  #                         upr = quantile(prod_S, 0.975),
  #                         prod_S = quantile(prod_S, 0.5)), 
  #             aes(x = S, ymin = lwr, ymax = upr, group = biogeo), 
  #             fill = 'transparent', colour = '#000', linetype = 'dashed') +
  
  geom_line(data = spatial_partialpred %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = '#000', linewidth = 0.75) +
  
  scale_x_log10() + scale_y_log10() + style +
  guides(colour = 'none', fill = 'none')


# Geospatial quantile with nulls ===================================================================

ggplot(data = vars_single) + style +
  labs(x = 'S', y = bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  geom_linerange(aes(x = S, ymin = prod_kg_l, ymax = prod_kg_u),
                 linewidth = 0.2, colour = '#888') +
  
  # data
  geom_point(aes(x = S, y = prod_kg, colour = biogeo), 
             shape = 21, size = 2, alpha = 0.3, fill = 'transparent') +
  scale_colour_cherulean(palette = 'cheridis', discrete = T) +
  
  # mod_geo S
  geom_ribbon(data = spatial_partialpred %>% ungroup %>% group_by(S) %>% 
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)), 
              aes(x = S, ymin = lwr, ymax = upr), 
              fill = 'transparent', colour = '#000', linetype = 'dashed') +
  
  geom_line(data = spatial_partialpred %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = '#000', linewidth = 1) +
  
  # 0.25 quantile
  geom_ribbon(data = geoq_partialpred %>% ungroup %>% group_by(S) %>% 
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)), 
              aes(x = S, ymin = lwr, ymax = upr), 
              fill = 'transparent', colour = '#000', linetype = 'dashed') +
  
  geom_line(data = geoq_partialpred %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = '#000', linewidth = 1) +
  
  # # null model lines
  # geom_ribbon(data = nullsp1_partial %>% ungroup %>% group_by(S) %>% 
  #               summarise(lwr = quantile(prod_S, 0.025),
  #                         upr = quantile(prod_S, 0.975),
  #                         prod_S = quantile(prod_S, 0.5)), 
  #             aes(x = S, ymin = lwr, ymax = upr), 
  #             fill = '#000', alpha = 0.1) +
  # 
  # geom_line(data = nullsp1_partial %>% ungroup %>% group_by(S) %>% 
  #             summarise(lwr = quantile(prod_S, 0.025),
  #                       upr = quantile(prod_S, 0.975),
  #                       prod_S = quantile(prod_S, 0.5)), 
  #           aes(x = S, y = prod_S), 
  #           colour = '#000', alpha = 0.5, linewidth = 0.75) +
  # 
  # geom_ribbon(data = nullsp2_partial %>% ungroup %>% group_by(S) %>% 
  #               summarise(lwr = quantile(prod_S, 0.025),
  #                         upr = quantile(prod_S, 0.975),
  #                         prod_S = quantile(prod_S, 0.5)), 
  #             aes(x = S, ymin = lwr, ymax = upr), 
  #             fill = '#000', alpha = 0.1) +
  # 
  # geom_line(data = nullsp2_partial %>% ungroup %>% group_by(S) %>% 
  #             summarise(lwr = quantile(prod_S, 0.025),
  #                       upr = quantile(prod_S, 0.975),
  #                       prod_S = quantile(prod_S, 0.5)), 
  #           aes(x = S, y = prod_S), 
  #           colour = '#000', alpha = 0.5, linewidth = 0.75) +
  # 
  # # null productivity
  # geom_ribbon(data = null_partialpred %>% ungroup %>% group_by(S) %>%
  #               summarise(lwr = quantile(prod_S, 0.025),
  #                         upr = quantile(prod_S, 0.975),
  #                         prod_S = quantile(prod_S, 0.5)),
  #             aes(x = S, ymin = lwr, ymax = upr),
  #             fill = '#444444', alpha = 0.2) +
  # 
  # geom_line(data = null_partialpred %>% ungroup %>% group_by(S) %>%
  #             summarise(lwr = quantile(prod_S, 0.025),
  #                       upr = quantile(prod_S, 0.975),
  #                       prod_S = quantile(prod_S, 0.5)),
  #           aes(x = S, y = prod_S),
  #           colour = 'grey30', alpha = 0.5, linewidth = 0.75) +
  
  scale_x_log10() + scale_y_log10() +
  guides(colour = 'none', fill = 'none')
