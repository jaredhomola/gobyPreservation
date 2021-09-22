## Load data and libraries
# Jared's setwd
library(ggmap)
library(tidyverse)
library(ggsn)
library(ggrepel)

# Need to create a free API key

# Create reference tibbles
samplingSites.tib <- tibble(
  site = c("Lake Michigan",
           "Cedar Creek"),
  lat = c(43.223175,
          43.305793),
  lon = c(-86.340194,
          -86.115011)
)

mapExtent.tib <- tibble(
  lat = c(43.125,
          43.42),
  long = c(-86.45,
           -86.05)
)

#### Create inset map
inset <-  ggplotGrob(
  ggmap(
    get_googlemap(
      center = c(lon = -85.82, lat = 43.773),
      zoom = 5,
      scale = 2,
      maptype = "satellite",
      color = "color")) +
    theme(
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = unit(c(0, 0,-1,-1), 'lines')) +
    geom_point(
      aes(x = -86.2, y = 43.25),
      shape = 18,
      color = "red",
      size = 4)
  #geom_rect(
  #  aes(
  #    xmin = -86.5,
  #    xmax = -86.0,
  #    ymin = 43.1,
  #    ymax = 43.4
  #  ),
  #  color = "red",
  #  fill = NA,
  #  size = 2
  #) 
)

#### Create site map
siteMap <-
    ggmap(get_googlemap(
      center = c(lon = -86.245769, lat = 43.262737),
      zoom = 11,
      scale = 2,
      maptype = 'terrain',
      color = 'color'
    )) +
    geom_point(data = samplingSites.tib,
               aes(x = lon, y = lat),
               color = 'red',
               size = 4) +
    geom_label_repel(
      aes(lon, lat, label = site),
      data = samplingSites.tib,
      size = 5,
      point.padding = 0.3,
      segment.color = 'grey50'
    ) +
    scalebar(
      x.min = -86.45,
      x.max = -86.05,
      y.min = 43.125,
      y.max = 43.42,
      transform = TRUE,
      dist_unit = "km",
      location = "bottomright",
      height = 0.03,
      st.dist = 0.03,
      dist = 5,
      model = 'WGS84'
    ) +
  xlab("Longitude") +
  ylab("Latitude") +
  inset(grob = inset, xmin = -86.47, xmax = -86.25, ymin = 43.3, ymax = 43.41) +
  theme_bw(base_size = 14)

ggsave("siteMap.png", units = "in", dpi = 300, width = 7, height = 7)













#### Mapping DAPC ####
library(ggmap)
library(rnaturalearth)
library(maps)
library(mapdata)
library(tidyverse)

mapInfo <- read_delim("G:/My Drive/Gobies/Project data/Ethanol Experiment/JMorpho/Lorencen et al figures and tables/latLong.csv",
                      delim = ",")

state_prov <- ne_states(c("united states of america", "canada"))
greatLakes <- subset(state_prov, name %in% c("Michigan", "Wisconsin", "Ontario", "Illinois", "Indiana", "Ohio", "Minnesota", "Pennsylvania", "New York"))

# compute the bounding box
bbox <- make_bbox(lat = latitude, lon = longitude, data = mapInfo)
bbox

# grab the maps from google
base <- get_map(location = bbox, 
                source = "stamen", 
                maptype = "toner-hybrid")

ggmap(base) + 
  geom_point(data = mapInfo, 
             mapping = aes(x = longitude, y = latitude))



ggplot() + 
  geom_polygon(data = greatLakes, aes(x = long, 
                               y = lat, 
                               group = group)) + 
  coord_fixed(xlim = c(-87, -85),  
              ylim = c(42, 44),
              ratio = 1.3)
