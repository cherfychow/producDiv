
# ProducDiv
# Figure code

require(ggplot2)
require(patchwork)
require(sf)
require(rnaturalearth)
require(ggMarginal)

# colour palettes
source('https://gist.githubusercontent.com/cherfychow/e9ae890fd16f4c86730748c067feee2b/raw/b2db138ab8164c237129308ea643de78cd54b252/cherulean.R')
style <- theme_classic(base_size = 13) + theme(panel.grid = element_blank(),
                                          axis.ticks = element_line(linewidth = 0.3),
                                          axis.line = element_line(linewidth = 0.3))


# Dataset map ---------------------------------------------------------------------------------

proj <- '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'

world <- ne_coastline('medium', returnclass = 'sf') %>%
  st_make_valid(.) %>% st_break_antimeridian(., lon_0 = 160)
world <- st_cast(world, to = 'MULTIPOLYGON')

# make a spatial point object with the data
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
            colour = '#000', linewidth = 0.75) +
  
  scale_x_continuous(breaks = c(1,5,10,50,100,300,500)) + 
  scale_y_continuous(breaks = c(0,0.001,.01, .1, 1, 10, 50, 100)) +
  coord_trans(y = "log10", x = 'log10') + style +
  guides(colour = 'none', fill = 'none')

# with human gravity
f_space_g <- ggplot(data = vars_single) +
  labs(x = 'Human hgravity', y =  bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
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

f_space_S/f_space_g

# split by bioregion
ggplot(data = vars_single) +
  labs(x = 'log(S)', y = bquote("log Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  geom_linerange(aes(x = S, ymin = prod_kg_l, ymax = prod_kg_u),
                linewidth = 0.2, colour = '#666') +
  geom_point(aes(x = S, y = prod_kg, colour = biogeo, fill = biogeo), 
             shape = 21, size = 2, alpha = 0.5) +
  
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
            colour = '#000', linewidth = 0.75) +

  scale_colour_cherulean(palette = 'cheridis', discrete = T) +
  scale_fill_cherulean(palette = 'cheridis', discrete = T) +
  coord_trans(x = 'log10', y = 'log10') +
  guides(colour = 'none', fill = 'none') +
  style + facet_wrap(~biogeo)


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
  geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = '#666') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#666') +
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
                  fill = '#333660', alpha = 0.3) +
  geom_line(data = T_partialpred %>% ungroup %>% group_by(S_slope) %>% 
              summarise(lwr = quantile(deltaProd_S, 0.025),
                        upr = quantile(deltaProd_S, 0.975),
                        deltaProd_S = quantile(deltaProd_S, 0.5)), 
            aes(x = S_slope, y = deltaProd_S), 
            colour = '#333660', linewidth = 0.5) +
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

f_T_g <- ggplot(data = vars_time_long) +
  labs(x = "Human gravity", y = bquote("∆Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  geom_linerange(aes(ymin = prod_slope - 1.96*prod_se, 
                     ymax = prod_slope + 1.96*prod_se, 
                     x = hgrav2 + 1, colour = quad), size = 0.2) +
  scale_colour_manual(values = c('#455EBA','#455EBA','#DBC539','#DBC539')) +
  geom_point(aes(x = hgrav2 + 1, y = prod_slope, fill = quad), 
             size = 2, shape = 21, colour = '#444') +
  scale_fill_manual(values = c('#455EBA','#455EBA','#DBC539','#DBC539')) +
  geom_ribbon(data = T2_pred, aes(ymin = prod_gpar_l, ymax = prod_gpar_u, x = hgrav2),
              fill = '#333660', alpha = 0.3) +
  geom_line(data = T2_pred, aes(x = hgrav2, y = prod_gpar), colour = '#333660', linewidth = 0.5) +
  scale_x_log10() + style

f_T_S
f_T_g


# ∆Prod ~ ∆S, log ratio adjacent years =============================================================

f_t_post1 <- ggplot(t_draws %>% filter(var_type == 'fixed')) +
  stat_slab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  scale_fill_manual(values = '#DBC539', na.value = 'grey90') +
  stat_spike(aes(x = values, y = prettyvar, ), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_t_post2 <- ggplot(t_draws %>% filter(var_type == 'randommodel')) +
  stat_slab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  scale_fill_manual(values = '#DBC539', na.value = 'grey90') +
  stat_spike(aes(x = values, y = prettyvar), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_t_post <- f_t_post1 / f_t_post2 + plot_layout(heights = c(2, 5))
f_t_post


## Partial prediction plot =========================================================================

vars_time <- vars_time %>% 
  mutate(quad = case_when(.default = 'neg',
                          logratio1_S > 0 & logratio1_prod > 0 ~ 'pos',
                          logratio1_S <= 0 & logratio1_prod <= 0 ~ 'pos')) %>% 
  mutate(quad = case_when(.default = quad,
                          logratio1_S > 0 & logratio1_prod > 0 & logratio1_S_l > 0 & logratio1_prod_l > 0 ~ 'pos_certain',
                          logratio1_S <= 0 & logratio1_S_u <= 0 & logratio1_prod <= 0 & logratio1_prod_u <= 0  ~ 'pos_certain',
                          logratio1_S > 0 & logratio1_S_l > 0 &  logratio1_prod <= 0 & logratio1_prod_u <= 0 ~ 'neg_certain',
                          logratio1_S <= 0 & logratio1_S_u <= 0 & logratio1_prod > 0 & logratio1_prod_l > 0 ~ 'neg_certain'))

f_T_S1 <- ggplot(data = vars_time_long) +
  labs(x = "∆S", y = bquote("∆Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = '#666') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#666') +
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
              fill = '#333660', alpha = 0.3) +
  geom_line(data = T_partialpred %>% ungroup %>% group_by(S_slope) %>% 
              summarise(lwr = quantile(deltaProd_S, 0.025),
                        upr = quantile(deltaProd_S, 0.975),
                        deltaProd_S = quantile(deltaProd_S, 0.5)), 
            aes(x = S_slope, y = deltaProd_S), 
            colour = '#333660', linewidth = 0.5) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 'dashed', colour = '#444', alpha = 0.6) +
  guides(colour = 'none', fill = 'none') +
  scale_x_continuous(limits = c(-2.5, 1.5)) +
  scale_y_continuous(limits = c(-3.5, 3)) + style

f_t_Sx <- ggplot(data = vars_time) +
  geom_density(aes(x = logratio1_S), fill = '#CCC', colour = '#666') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#666') +
  theme_void() + scale_x_continuous(limits = c(-2.5, 1.5))

f_t_Sy <- ggplot(data = vars_time) +
  geom_density(aes(y = logratio1_prod), fill = '#CCC', colour = '#666') +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = '#666') +
  theme_void() + scale_y_continuous(limits = c(-3.5, 3))

f_t_S <- f_t_Sx + plot_spacer() + f_t_S1 + f_t_Sy + 
  plot_layout(ncol = 2, nrow = 2, widths = c(5,1), heights = c(1,6))

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

f_t_S
f_t_g

# proportions in positive-negative for long and short term fluxes

dt_tchange <- vars_time_long %>% group_by(quad) %>% 
  summarise(n = n()) %>% mutate(group = 'full')
dt_tchange <- vars_time %>% group_by(quad) %>% 
  summarise(n = n()) %>% mutate(group = 'adjacent') %>% 
  bind_rows(dt_tchange, .)

ggplot(data = dt_tchange) +
  geom_col(aes(y = n, fill = quad, x = group), position = 'fill') + 
  guides(fill = 'none') +
  scale_fill_manual(values = c('#114c54', '#d989e5')) +
  theme_minimal() + labs(x = NULL, y = NULL) +
  theme(panel.grid = element_blank())
