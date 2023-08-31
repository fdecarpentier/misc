### General purpose packages ###
library(tidyverse)
library(cowplot)
library(ggblend)
library(units)
library(viridis)
library(ggnewscale)

### Spatial analysis packages ###
library(sf)
library(stars) 
library(terra)
library(tidyterra)
library(gstat)
library(osmdata)
library(elevatr)

### Obtain the basic features ###
# Features are obtained from OpenStreetMap
keys <- c("name", "natural", "highway", "building", "man_made")
features <- lapply(seq_along(keys), function (i) {
  matrix(c(-123.015, 37.6934, -122.9966, 37.7052), nrow = 2, dimnames = list(c("x","y"), c("min","max"))) |>
    opq() |> 
    add_osm_feature(key = keys[i]) |>
    osmdata_sf()}) |> 
  set_names(keys)

# For an easy use of the boundaries and bounding box of the SE island (+ Maintop or not)
boundaries <- list(two_isl = features$natural$osm_polygons |> filter(name %in% c("Maintop Island", "Southeast Farallon")),
                   main_isl = features$natural$osm_polygons |> filter(name %in% c("Southeast Farallon")))
bbox <- lapply(boundaries, st_bbox)

### Create plots with soil metrics as colored points ###
points_large <- read.csv("farallon_soil.csv")
points_large_sf <- st_as_sf(points_large, crs = "EPSG:4326", coords = c("Long", "Lat"))

points_long_sf <- points_large_sf |> 
  pivot_longer(cols = !c(Plot, geometry), names_to = "Analysis", values_to = "Value") %>% 
  mutate(Analysis_full = case_when(
    Analysis == "Anthr"~ "Anthropogenic disturb.",
    Analysis == "Burr" ~ "Burrow density_/m2", 
    Analysis == "Elev" ~ "Elevation_m", 
    Analysis == "Ec" ~ "Elec. conductivity 1:1_dS/m", 
    Analysis == "Bare" ~ "Bare ground_%", 
    Analysis == "Litter" ~ "Litter_%", 
    Analysis == "Rock" ~ "Rocky substrate_%", 
    Analysis == "Solar" ~ "Solar_MJ/cm2/yr", 
    Analysis == "Depth" ~ "Depth_cm", 
    Analysis == "Slope" ~ "Slope_°", 
    T ~ Analysis),
    Value = ifelse(Analysis %in% c("Bare", "Litter", "Rock"), Value*100, Value)) |>
  separate(Analysis_full, into = c("Analysis_disp", "Unit"), sep = "_", remove = F) |>
  filter(Analysis %in% c("Rock", "Litter", "Depth", "pH", "Ec", "Bare"))

points_plots <- lapply(split(points_long_sf, points_long_sf$Analysis), function(i) {
  ggplot() + 
    geom_sf(data = boundaries$main_isl, fill = "grey90", color = NA) + 
    geom_sf(data = i, aes(color = Value)) +
    scale_color_viridis(option = "B", direction = -1) +
    labs(title = unique(i$Analysis_disp),
         color = unique(i$Unit) %>% ifelse(is.na(.), "", .))})

# points_grid <- plot_grid(plotlist = points_plots)

### A few plots (not used in the final figure) to select what sites to keep ###
sites_plot <- ggplot(points_large, aes(x = Long, y = Lat, color = as.character(Plot))) +
  geom_point() +
  geom_text(aes(label = Plot), nudge_x = .0003) +
  geom_sf(data = features$natural$osm_polygons$geometry, inherit.aes = F, fill = NA) +
  coord_sf() +
  theme_void() +
  theme(legend.position = "none")
sites_plot

sites_selection <- c(2, 4, 6, 8, 9, 12, 15, 16, 19, 20, 22, 24, 28, 32, 34, 38, 41, 42)
length(sites_selection)

sites_density <- ggplot(points_long_sf, aes(x=Value)) +
  geom_density() +
  geom_text(aes(y = 0, label = Plot, color = as.character(Plot))) +
  geom_point(data = points_long_sf |> filter(Plot %in% sites_selection), aes(y = .02, color = as.character(Plot))) +
  facet_wrap(~Analysis, scales = "free") +
  theme_classic() +
  theme(legend.position = "none")

sites_final <- points_large |> 
  select(Plot, Lat, Long) |>
  mutate(New = "Established site") |>
  filter(Plot %in% sites_selection) |>
  rbind(read.csv("farallon_extra.csv") |> mutate(New = "New site")) |> # add new sampling plots
  mutate(Meta = ifelse(Plot %in% c(8, 15, 41, 42), "Culture + metagenomic", "Culture only")) |>
  st_as_sf(crs = "EPSG:4326", coords = c("Long", "Lat"))

### Interpolation (when relevant) with the kriging method ###
# Determination of a buffered convex hull = contour of the points cloud
hull <- points_large_sf |>
  summarise() |>
  st_convex_hull() |> 
  st_buffer(50)

# Creation of an empty raster within the coastline and hull boundaries
grid <- st_as_stars(bbox$main_isl, dx = 0.000025) |>
  st_crop(boundaries$main_isl) |>
  st_crop(hull)

# Addition of manual cutoff for the variograms when necessary
points_long_sf_krige <- filter(points_long_sf, Analysis %in% c("Depth", "Rock", "Ec", "Bare")) |>
  mutate(Cutoff = case_when(
    Analysis == "Depth" ~ .402,
    Analysis == "Bare" ~ .6,
    # Analysis == "Rock" ~ .25,
    # Analysis == "Litter" ~ .807,
    # Analysis == "Burr" ~ .055,
    T ~ NA)) # If NA, will compute auto variogram 

# Calculation of variograms with auto cutoff
varios_auto <- points_long_sf_krige |>
  filter(is.na(Cutoff)) %>% 
  split(., .$Analysis_full) |>
  lapply(function (i){
    variogram(Value~1, i) |>
      mutate(id = unique(i$Analysis_full))})

# Calculation of variograms with manual cutoff
varios_cut <- points_long_sf_krige |>
  filter(!is.na(Cutoff)) %>% 
  split(., .$Analysis_full) |>
  lapply(function (i){
    variogram(Value~1, i, cutoff = unique(i$Cutoff), width = unique(i$Cutoff/15)) |>
      mutate(id = unique(i$Analysis_full))})

# Fitting of exponential models to variograms
varios <- c(varios_auto, varios_cut)
varios_model <- lapply(varios, function (i){
  fit.variogram(i, vgm(c("Exp"), fit.kappa = T))}) #, "Gau", "Sph", "Mat"
varios_model

# Visualization of variograms and fitted models
varios_model_line <- lapply(seq_along(varios_model), function (i){
  variogramLine(varios_model[[i]], maxdist = max(varios[[i]]$dist))}) |>
  set_names(names(varios_model)) |>
  bind_rows(.id = "id")

varios_plot <- bind_rows(varios) |> 
  ggplot(aes(dist, gamma)) +
  geom_line(data = varios_model_line, color = "deeppink3") +
  geom_point() +
  geom_text(aes(label = np), position = position_nudge(x = 0.015)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Distance h (m)", y = expression(gamma(h))) +
  facet_wrap(~id, scales = "free") +
  theme_classic()

# Calculation of interpolated raster with the kriging method
kriges <- lapply(varios, function (i){
  points_long_sf %>%
    filter(Analysis_full == unique(i$id)) %>% 
    krige(Value~1, ., grid, fit.variogram(i, vgm(c("Exp"), fit.kappa = T)))}) #, "Gau", "Sph", "Mat"

# Creation of interpolated raster
kriges_plots <- lapply(seq_along(kriges), function (i){
  ggplot() +
    geom_sf(data = boundaries$main_isl, fill = "grey90", color = NA) +
    geom_stars(data = kriges[[i]], aes(fill = var1.pred, x = x, y = y)) + 
    geom_sf(data = points_long_sf, color = "grey70", shape = 21, stroke = .5) +
    scale_fill_viridis(option = "B", direction = -1, na.value = NA) +
    labs(title = names(kriges)[[i]] |> word(1, sep = "_"),
         fill = names(kriges)[[i]] |> word(2, sep = "_"))}) |>
  set_names(names(kriges))

# kriges_grid <- plot_grid(plotlist = kriges_plots, align = "hv")

### Interpolation with inverse distance weight (when krige is not possible) ###
idws <- filter(points_long_sf, Analysis %in% c("Litter", "pH")) %>% 
  split(., .$Analysis_full) |> 
  lapply(function (i){
    idw(i$Value~1, i, grid, idp = 1.8)})

idws_plots <- lapply(seq_along(idws), function (i){
  ggplot() +
    geom_sf(data = boundaries$main_isl, fill = "grey90", color = NA) +
    geom_stars(data = idws[[i]], aes(fill = var1.pred, x = x, y = y)) + 
    geom_sf(data = points_long_sf, color = "grey70", shape = 21, stroke = .5) +
    scale_fill_viridis(option = "B", direction = -1, na.value = NA) +
    labs(title = names(idws)[[i]] |> word(1, sep = "_"),
         fill = names(idws)[[i]] |> word(2, sep = "_") %>% ifelse(is.na(.), "", .))}) |>
  set_names(names(idws))

# idws_grid <- plot_grid(plotlist = idws_plots, align = "hv")

### Create an optimal zone for Chlamydomonas reinhardtii ###
# Create a list of raster with a cutoff
# The code could be prettier, I know ... but I don't have time right now ...
zones <- list(`Elec. conductivity 1:1_dS/m` = kriges$`Elec. conductivity 1:1_dS/m` |>
                mutate(var1.pred = ifelse(var1.pred < 1, var1.pred, NA)), # EC (salinity) below 1
              pH = idws$pH |>
                 mutate(var1.pred = ifelse(var1.pred > 4.8, var1.pred, NA)), # pH above 4.8
              `Rocky substrate_%` = kriges$`Rocky substrate_%` |>
                mutate(var1.pred = ifelse(var1.pred < 25, var1.pred, NA)), # Rock content below 25%
              `Litter_%` = idws$`Litter_%` |>
                mutate(var1.pred = ifelse(var1.pred > 1, var1.pred, NA))) # Litter content above 1%

# Transforms stars rasters in sf polygons
zones_sf <- lapply(zones, function (i){
  i |>
    select(var1.pred) |>
    st_as_sf(merge = T) |>
    st_union()}) # Very slow, would be nice to find a better method

# Creates an intersection of all zones
zones_final <- st_intersection(zones_sf$`Elec. conductivity 1:1_dS/m`, 
                            zones_sf$`Rocky substrate_%`) |>
  st_intersection(st_intersection(zones_sf$`Litter_%`, zones_sf$pH)) |>
  smoothr::smooth(method = "ksmooth", smoothness = 20)

# ggplot() +
#   geom_sf(data = boundaries$two_isl, fill = "grey90", color = "grey10") +
#   geom_sf(data = zones_final, fill = "grey", alpha = .5) +
#   geom_sf(data = zones_sf$`Rocky substrate_%`, color = "orange", fill = NA) +
#   geom_sf(data = zones_sf$pH, color = "pink3", fill = NA) +
#   geom_sf(data = zones_sf$`Elec. conductivity 1:1_dS/m`, color = "blue", fill = NA) +
#   geom_sf(data = zones_sf$Litter, color = "red3", fill = NA) +
#   geom_text(data = points_large, aes(x = Long, y = Lat, label = Plot), nudge_x = -.00025) +
#   geom_sf(data = points_long_sf, color = "grey10") +
#   theme_void()
# ggsave("sites_zone.pdf", width = 10, height = 5)

### Create an elevation maps ###
# Get elevation raster with elevatr package
elev <- get_elev_raster(locations = boundaries$two_isl, z = 14, clip = "locations") %>% 
  st_as_stars() %>% 
  st_crop(boundaries$two_isl) %>% 
  setNames("elevation") 

# Calculate contour lines
elev_contour <- st_contour(elev, contour_lines = F, breaks = 10*-2:8+2)

# Select the labels to display
labels <- st_as_sf(features$name$osm_points) |>
  filter(name %in% c("Southeast Farallon", "Maintop Island", "Maintop Bay", "Fisherman Bay"))

# Moves the labels
labels$geometry[labels$name == "Maintop Island"] <- st_point(c(-123.0124, 37.6998)) |> st_sfc(crs = "EPSG:4326")                          
labels$geometry[labels$name == "Maintop Bay"] <- st_point(c(-123.0085, 37.7005)) |> st_sfc(crs = "EPSG:4326")                          
labels$geometry[labels$name == "Fisherman Bay"] <- st_point(c(-123.0017, 37.7017)) |> st_sfc(crs = "EPSG:4326")                          
labels$geometry[labels$name == "Southeast Farallon"] <- st_point(c(-123.0019, 37.69855)) |> st_sfc(crs = "EPSG:4326")                          

# Create the elevation plot
elev_plot <- ggplot() +
  geom_sf(data = boundaries$two_isl, fill = "grey90", color = "grey10") + 
  geom_sf(data = zones_final, aes(fill = ""), color = "#cfff00", alpha = .5) +
  geom_sf(data = features$highway$osm_lines$geometry, color = "grey30", linewidth = .1) +
  geom_sf(data = features$building$osm_polygons$geometry, color = "seashell4", fill = "seashell1", alpha = .4) +
  geom_sf(data = features$man_made$osm_polygons$geometry, color = "seashell4", fill = "seashell1", alpha = .4) +
  geom_sf(data = elev_contour, color = "grey30", fill = NA, linetype = "dotted") + 
  geom_sf_text(data = labels, aes(label = name)) +
  geom_sf_text(data = sites_final, aes(label = Plot, color = Meta), nudge_x = -.0001, hjust = 1) +
  geom_sf(data = sites_final, aes(shape = New, color = Meta)) +
  scale_fill_manual(values = c("#cfff04"), labels = expression(paste(italic("C. reinhardtii"), " optimal zone"))) +
  scale_color_manual(values = c("#fa2f84", "#000bdc")) +
  guides(shape = guide_legend(order = 1), color = guide_legend(order = 2), fill = guide_legend(order = 3)) +
  coord_sf(xlim = c(bbox$two_isl$xmin, bbox$two_isl$xmax), ylim = c(bbox$two_isl$ymin, bbox$two_isl$ymax)) +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.justification = c(0.5, 0),
        legend.position = c(.4, .08))

### Create a plot with hillshade effect ###
# Get elevation raster with elevatr on a buffered zone
elev_star <- get_elev_raster(locations = st_buffer(boundaries$two_isl, 20), z = 14, clip = "locations") |>
  st_as_stars(elev_buff) |> 
  setNames("elevation") 

# Resample the elevation raster to increase resolution
elev_grid <- st_bbox(elev_star) |> st_as_stars(dx = 0.000015)
elev_warp <- st_warp(elev_star, elev_grid, use_gdal = T, method = "cubicspline") |> #Don't pay attention to the warning (cf. documentation)
  setNames("elevation") |>
  rast()

# Calculate slope and aspect rasters based on elevation
sl <- terrain(elev_warp, "slope", unit = "radians")
asp <- terrain(elev_warp, "aspect", unit = "radians")

# Creates a hillshade effect raster with several directions
hill <- map(seq(180, 300, length.out = 4), function(dir){ #Be carefull, does not work properly if angles around 360
  shade(sl, asp, angle = 30, direction = dir, normalize= T)}) |>
  rast() |>
  sum()

# Crops the hillshade and elevation raster to actual boundaries
hill_crop <- mask(hill, vect(boundaries$two_isl), touches = F) |> st_as_stars()
elev_crop <- mask(elev_warp, vect(boundaries$two_isl), touches = F) |> st_as_stars()

# Created a final hillshade plot
hill_plot <- ggplot() +
  list(
    geom_stars(data = hill_crop, show.legend = F, alpha = 1), 
    scale_fill_distiller(palette = "Greys", na.value = NA), #"lightsteelblue3" is okay-ish for sea
    new_scale_fill(),
    geom_stars(data = elev_crop, show.legend = F, alpha = .8),
    scale_fill_hypso_tint_c(limits = c(-10, 75), palette = "arctic_hypso")) %>% #"dem_screen" and "arctic_hypso" are nice
  blend("multiply") +
  geom_sf(data = boundaries$two_isl, fill = NA) +
  new_scale_fill() +
  geom_sf(data = zones_final, aes(fill = ""), color = "#cfff00", alpha = .5) +
  geom_sf(data = features$highway$osm_lines$geometry, linewidth = .1, alpha = 1) +
  geom_sf(data = features$building$osm_polygons$geometry, color = "seashell4", fill = "seashell1", alpha = .4) +
  geom_sf(data = features$man_made$osm_polygons$geometry, color = "seashell4", fill = "seashell1", alpha = .4) +
  geom_sf_text(data = labels, aes(label = name)) +
  geom_sf_text(data = sites_final, aes(label = Plot, color = Meta), nudge_x = -.0001, hjust = 1) +
  geom_sf(data = sites_final, aes(shape = New, color = Meta)) +
  scale_color_manual(values = c("#fa2f84", "#000bdc")) +
  scale_fill_manual(values = c("#cfff04"), labels = expression(paste(italic("C. reinhardtii"), " optimal zone"))) +
  guides(shape = guide_legend(order = 1), color = guide_legend(order = 2), fill = guide_legend(order = 3)) +
  coord_sf(xlim = c(bbox$two_isl$xmin, bbox$two_isl$xmax), ylim = c(bbox$two_isl$ymin, bbox$two_isl$ymax)) +
  theme_void() +
  theme(legend.title = element_blank(),
      legend.direction = "horizontal",
      legend.justification = c(0.5, 0),
      legend.position = c(.4, .08))

### Produces a final arrangement plot ###
plots_selection <- c(points_plots[c("pH")], idws_plots[c("Litter_%")], kriges_plots) |>
  lapply(function (i){
    i +
    geom_sf(data = boundaries$main_isl, fill = NA, color = "grey10") +
    coord_sf(xlim = c(bbox$main_isl$xmin, -122.997), expand = F) +
    theme_void()+
    theme(legend.justification = c(0, .5), 
          legend.position = c(.85, .5))})

selection_grid <- plot_grid(plotlist = plots_selection, ncol = 3)

full_grid_hillshade <- plot_grid(selection_grid, hill_plot, labels = c("A", "B"), nrow = 2, rel_heights = c(1, 1.4))
ggsave("farallon_hillshade.pdf", full_grid_hillshade, width = 8, height = 8, scale = 1.1)
ggsave("farallon_hillshade.png", full_grid_hillshade, width = 8, height = 8, scale = 1.1, device = png, type = "cairo", bg = "white")

full_grid_contour <- plot_grid(selection_grid, elev_plot, labels = c("A", "B"), nrow = 2, rel_heights = c(1, 1.4))
ggsave("farallon_contour.pdf", full_grid_contour, width = 8, height = 8, scale = 1.1)
ggsave("farallon_contour.jpg", full_grid_contour, width = 8, height = 8, scale = 1.1, dpi = 1200)


