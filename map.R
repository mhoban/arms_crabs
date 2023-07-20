library(tidyverse)
library(sf)
library(mapdata)
library(marmap)
library(patchwork)
library(here)
library(rnaturalearth)
library(rnaturalearthhires)


if (!exists("crab_data")) {
  crab_data <- read_csv(here("data","metadata.csv"),col_types = cols())
}

crab_sites <- crab_data %>%
  distinct(island,island_group,lat,lon,.keep_all = FALSE) %>%
  mutate( island = fct_reorder(island,-lat) ) %>%
  st_as_sf(coords=c("lon","lat"),crs=4326, remove=F) #%>%
# st_shift_longitude()

hawaii <- ne_states("united states of america",returnclass = "sf") %>%
  filter(name == "Hawaii")

# hawaii <- rnaturalearthhires::states10 %>%
#   st_as_sf() %>%
#   filter(name == "Hawaii")

# run this section --------------------------------------------------------
trans_box <- function(b,crs) {
  box_before <- st_sfc(
    st_point(b[c(1,2)]),
    st_point(b[c(3,4)]),
    crs=4326
  )
  box_before %>%
    st_transform(crs) %>%
    st_bbox()
}


box_before <- st_bbox(st_buffer(crab_sites,units::as_units(50,"km")))
box_before <- st_bbox(crab_sites)
# box_before[2] <- 18
# box_before[4] <- 29
bath <- getNOAA.bathy(lon1 = box_before[1], lon2 = box_before[3], lat1 = box_before[2], lat2 = box_before[4], resolution = 5) 
transformer <- str_glue("+proj=leac +lat_1={box_before[2]} +lat_2={box_before[4]} +lon_0={mean(box_before[c(1,3)])}")

box_after <- box_before %>%
  trans_box(transformer)
mappr <- hawaii %>%
  st_transform(transformer)
trans <- crab_sites %>% 
  st_transform(transformer)

bath_latlon <- as.raster(bath) %>% 
  as.bathy() %>% 
  as.xyz() %>%
  as_tibble() %>%
  drop_na(V3)

bath_proj <- as.raster(bath) %>% 
  raster::projectRaster(crs = transformer) %>% 
  as.bathy() %>% 
  as.xyz() %>%
  as_tibble() %>%
  drop_na(V3)

pmnm <- st_read(here("data","gis","pmnm.shp")) %>%
  st_transform(transformer)

nwp <- ggplot() + 
  geom_contour(data=bath_proj,mapping=aes(x=V1,y=V2,z = V3),  color="#eeeeee") +
  geom_sf(data=mappr) + 
  geom_sf(data=trans,aes(color=island),size=3) +
  geom_sf(data=pmnm,fill="#6ed8f6",alpha=0.2) + 
  # theme_minimal() + 
  coord_sf(xlim=box_after[c(1,3)],ylim=box_after[c(2,4)],expand=FALSE) +
  labs(x="",y="")
nwp
