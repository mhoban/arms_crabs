library(tidyverse)
library(sf)
library(marmap)
library(fs)
library(here)
library(rnaturalearth)

if (!exists("crab_data")) {
  crab_data <- read_csv(here("data","metadata.csv"),col_types = cols())
}

crab_sites <- crab_data %>%
  distinct(island,island_group,lat,lon,.keep_all = FALSE) %>%
  mutate( island = fct_reorder(island,-lat) ) %>%
  st_as_sf(coords=c("lon","lat"),crs=4326, remove=F) #%>%

hawaii <- ne_states("united states of america",returnclass = "sf") %>%
  filter(name == "Hawaii")

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
bath <- getNOAA.bathy(lon1 = box_before[1], lon2 = box_before[3], lat1 = box_before[2], lat2 = box_before[4], resolution = 3) 
transformer <- str_glue("+proj=leac +lat_1={box_before[2]} +lat_2={box_before[4]} +lon_0={mean(box_before[c(1,3)])}")

box_after <- box_before %>%
  trans_box(transformer)

bath_proj <- as.raster(bath) %>% 
  raster::projectRaster(crs = transformer) %>% 
  as.bathy() %>% 
  as.xyz() %>%
  as_tibble() %>%
  drop_na(V3) %>%
  filter(V3 < 0)

bpr <- as.raster(bath) %>% 
  raster::projectRaster(crs = transformer) %>%
  as.bathy() %>%
  fortify.bathy()

pmnm <- st_read(here("data","gis","pmnm.shp")) %>%
  st_transform(transformer)


# section break -----------------------------------------------------------

island_transform <- tribble(
  ~island,~new_name,
  "Kure Atoll","Hōlanikū",
  "Pearl and Hermes Reef","Manawai",
  "Lisianski","Kapou",
  "French Frigate Shoals","Lalo",
  "Kauai","Kaua‘i",
  "Oahu","O‘ahu",
  "Maui","Maui",
  "Hawaii Island","Hawai‘i Island"
)

cs <- crab_sites %>%
  left_join(island_transform,by="island") %>%
  mutate(
    new_name = fct_reorder(new_name,-lat)
  )

plotz <- ggplot() + 
  # geom_contour(data=bath_proj,mapping=aes(x=V1,y=V2,z = V3),  color="#cecece") +
  # geom_raster(data=bpr,aes(x=x,y=y,fill=z)) +
  # geom_tile(data=bath_proj,mapping=aes(x=V1,y=V2,fill = V3)) +
  # scale_fill_gradient2(low="dodgerblue4", mid="gainsboro", high="darkgreen") +
  geom_sf(data=hawaii %>% st_transform(transformer), color="black", fill="green", linewidth=0.4) + 
  scale_fill_brewer(palette="Accent",name="Island") + 
  geom_sf(data=pmnm,fill="#6ed8f6",alpha=0.1) + 
  geom_sf(data=cs %>% st_transform(transformer),aes(fill=new_name),size=3,shape=21,color="black") + 
  coord_sf(xlim=box_after[c(1,3)],ylim=box_after[c(2,4)]) +
  scale_x_continuous(breaks = -seq(180, 155, -5)) +
  theme_minimal() +
  xlab("Longitude") + 
  ylab("Latitude")

labela <- tribble(
  ~lon,~lat,~text,
  -170,25.4,"Paphānaumokuākea Marine\nNational Monument\n(NWHI)"
) %>%
  st_as_sf(coords=c("lon","lat"),crs=4326, remove=F)  %>%
  st_transform(transformer)

plotz <- plotz + 
  geom_sf_text(aes(label=text),data=labela,size=3,angle=-17.0)
plotz

dir_create(here("output","figures"),recurse=TRUE)
ggsave(here("output","figures","sample_map.svg"),plotz,width=10,height=6,units="in",device=svg)

