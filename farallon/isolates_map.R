library(tidyverse)
library(data.table)
library(maps)
library(ggrepel)

World_map <- map_data("world") %>% tibble()

Cre_coord <- fread("isolates_coord.csv") %>% 
  mutate(type = factor(type, levels = c("Lab", "Historical", "New isolate", "Candidate")))

Cre_region <- Cre_coord %>% 
  filter(is.na(lat)) %>% 
  select(name, region, type)

Cre_outline <- World_map %>% 
  filter(region %in% Cre_region$region) %>% 
  right_join(Cre_region, by = "region", relationship = "many-to-many")

Cre_label <- Cre_outline %>% 
  filter(is.na(subregion)) %>% 
  group_by(name) %>% 
  filter(ifelse(region == "France", lat == max(lat), lat == min(lat))) %>%
  select(name, long, lat, type) %>% 
  rbind(Cre_coord %>% filter(!is.na(lat))) %>% 
  group_by(type, lat, long) %>% 
  summarise(name = paste(name, collapse = "\n"))

Cre_map <- ggplot(data = Cre_coord, aes(x = long, y = lat)) +
  geom_polygon(data = World_map, aes(x = long, y = lat, group = group), 
               fill = "seashell2", colour = "white", linewidth = .2) +  
  geom_point(aes(color = type), shape = 16, alpha = .8) +
  geom_polygon(data = Cre_outline, aes(x = long, y = lat, group = group, color = type), 
               fill = NA, linewidth = .5) +
  geom_text_repel(data = Cre_label, aes(label = name, color = type),
                  min.segment.length = 0, max.overlaps = 40, hjust = 0, seed = 132, force = 2, force_pull = 1) +
  scale_color_manual(values = c("tomato2", "deeppink2", "steelblue3", "green4")) +
  coord_equal() +
  theme_void() +
  theme(legend.position = c(.68, .4),
        legend.title = element_blank(),
        legend.background = element_blank())
Cre_map

ggsave("isolates_map.pdf", Cre_map, width = 20, height = 10)
ggsave("isolates_map.jpg", Cre_map, width = 20, height = 10, dpi = 800)

