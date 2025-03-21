
# ProducDiv
# Figure code

require(ggplot2)
require(patchwork)
require(sf)
require(rnaturalearth)
require(ggMarginal)

# colour palettes
source('https://gist.githubusercontent.com/cherfychow/e9ae890fd16f4c86730748c067feee2b/raw/b2db138ab8164c237129308ea643de78cd54b252/cherulean.R')

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
  scale_fill_cherulean(palette = 'quoyi', discrete = T, alpha = 0.4) +
  scale_size_continuous(breaks = c(2,5,10,15,20))
f_map_single <- ggplot() +
  geom_sf(data = world, color='grey75', fill='grey95') +
  geom_sf(data = point_single, 
          aes(fill = DataID), shape = 21, size = 2) +
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) +
  theme_minimal() + guides(fill = 'none') +
  scale_fill_cherulean(palette = 'quoyi', discrete = T, alpha = 0.4)
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
  scale_fill_cherulean(palette = 'quoyi', discrete = T, alpha = 0.4)

f_map_eco <- ggplot() +
  geom_sf(data = world, color='grey75', fill='grey95') +
  geom_sf(data = ecoR, aes(fill = region)) +
  geom_sf(data = points, shape = 21, size = 1.5, fill = 'white') +
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) +
  theme_minimal(base_size = 13, base_family = 'Helvetica') +
  theme(legend.position = 'bottom') +
  scale_fill_cherulean(palette = 'quoyi', discrete = T)

f_map_eco / f_map


# Spatial prod ~ S

style <- theme_bw(base_size = 13) + theme(panel.grid = element_blank())

ggplot(data = vars_single) +
  labs(x = 'S', y = bquote("Productivity (kg⋅ha"^-1~"⋅yr"^-1~")")) +
  geom_errorbar(aes(x = S, ymin = prod_kg_yr_l, ymax = prod_kg_yr_u),
                size = 0.2, colour = '#666') +
  geom_point(aes(x = S, y = prod_kg_yr, fill = biogeo), shape = 21, size = 2, 
             colour = '#000', alpha = 0.5) +
  scale_fill_cherulean(palette = 'quoyi', discrete = T) +
  style + 
  scale_y_log10()


# ∆Prod ~ ∆S, long term trajectory

# slope scatter
vars_time_long <- vars_time_long %>% mutate(quad = case_when(.default = 'neg',
                                                             S_slope > 0 & prod_slope > 0 ~ 'pos',
                                                             S_slope <= 0 & prod_slope <= 0 ~ 'pos'))
f_Timelong1 <- ggplot(data = vars_time_long) +
  labs(x = '∆log(S)', y = bquote("∆log(Productivity) (kg⋅ha"^-1~"⋅yr"^-2~")")) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = '#666') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#666') +
  geom_linerange(aes(xmin = S_slope - 1.96*S_se, 
                     xmax = S_slope + 1.96*S_se, 
                     y = prod_slope, colour = quad), size = 0.2, alpha = 0.25) +
  geom_linerange(aes(ymin = prod_slope - 1.96*prod_se, 
                     ymax = prod_slope + 1.96*prod_se, 
                     x = S_slope, colour = quad), size = 0.2, alpha = 0.25) +
  geom_point(aes(x = S_slope, y = prod_slope, colour = quad), size = 2, alpha = 0.3) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 'dashed', colour = '#444', alpha = 0.6) +
  guides(colour = 'none') +
  scale_colour_manual(values = c('#114c54', '#d989e5')) +
  scale_x_continuous(limits = c(-3, 3)) +
  scale_y_continuous(limits = c(-4, 3)) + style

f_Timelongx <- ggplot(data = vars_time_long) +
  geom_density(aes(x = S_slope), fill = '#CCC', colour = '#666') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#666') +
  theme_void() + scale_x_continuous(limits = c(-3, 3))

f_Timelongy <- ggplot(data = vars_time_long) +
  geom_density(aes(y = prod_slope), fill = '#CCC', colour = '#666') +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = '#666') +
  theme_void() + scale_y_continuous(limits = c(-4, 3))

f_Timelongx + plot_spacer() + f_Timelong1 + f_Timelongy + 
  plot_layout(ncol = 2, nrow = 2, widths = c(5,1), heights = c(1,6))

# density grid
ggplot(data = vars_time_long) +
  labs(x = '∆S', y = bquote("∆Productivity (kg⋅ha"^-1~"⋅yr"^-2~")"), 
       title = 'Diversity-productivity trajectories (time series slope)') +
  stat_bin2d(bins = 100, aes(x = S_slope, y = prod_slope)) +
  scale_fill_cherulean(palette = 'quoyi', discrete = F) +
  theme_bw() + theme(panel.grid = element_blank())


# ∆Prod ~ ∆S, log ratio adjacent years
vars_time <- vars_time %>% mutate(quad = case_when(.default = 'neg',
                                                   logratio1_S > 0 & logratio1_prod > 0 ~ 'pos',
                                                   logratio1_S <= 0 & logratio1_prod <= 0 ~ 'pos'))

f_Timeadj1 <- ggplot(data = vars_time) +
  labs(x = '∆log(S)', y = bquote("∆log(Productivity) (kg⋅ha"^-1~"⋅yr"^-2~")")) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = '#666') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#666') +
  geom_linerange(aes(xmin = logratio1_S_l, 
                     xmax = logratio1_S_u, 
                     y = logratio1_prod, colour = quad), size = 0.2, alpha = 0.25) +
  geom_linerange(aes(ymin = logratio1_prod_l, 
                     ymax = logratio1_prod_u, 
                     x = logratio1_S, colour = quad), size = 0.2, alpha = 0.25) +
  geom_point(aes(x = logratio1_S, y = logratio1_prod, colour = quad), size = 2, alpha = 0.3) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 'dashed', colour = '#444', alpha = 0.6) +
  guides(colour = 'none') +
  scale_colour_manual(values = c('#114c54', '#d989e5')) +
  scale_x_continuous(limits = c(-4, 4)) +
  scale_y_continuous(limits = c(-6, 8)) + style

f_Timeadjx <- ggplot(data = vars_time) +
  geom_density(aes(x = logratio1_S), fill = '#CCC', colour = '#666') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#666') +
  theme_void() + scale_x_continuous(limits = c(-4, 4))

f_Timeadjy <- ggplot(data = vars_time) +
  geom_density(aes(y = logratio1_prod), fill = '#CCC', colour = '#666') +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = '#666') +
  theme_void() + scale_y_continuous(limits = c(-6, 8))

f_Timeadjx + plot_spacer() + f_Timeadj1 + f_Timeadjy + 
  plot_layout(ncol = 2, nrow = 2, widths = c(5,1), heights = c(1,6))

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
