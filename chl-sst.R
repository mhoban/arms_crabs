library(tidyverse)
library(sf)
library(raster)
library(here)
library(lubridate)
library(httr)
library(fs)

# set random seed for consistency
set.seed(31337)


## download SST and CHL data from NASA and aggregate it by ARMS sampling site


# pre-parse metadata and get deployment & recovery years as well as deployment range
arms <- read_csv(here("data","metadata.csv")) %>%
  rename(sample=unit) %>%
  st_as_sf(coords=c("lon","lat"),crs=4326,remove=FALSE) %>%
  separate_wider_regex(sample,patterns=c(recovery_year="^[0-9]+","_",sample=".*"),cols_remove=FALSE)  %>%
  mutate(recovery_year=as.numeric(recovery_year)) %>%
  dplyr::select(sample,everything()) %>%
  unite("range",deployment_year,recovery_year,sep="--",remove=FALSE) 


# sample 30 months from each date range
dates <- arms %>%
  group_by(range) %>%
  group_map(\(data,grp) {
    start_date <- str_glue("{min(data$deployment_year)}-01-01") %>%
      as.Date()
    end_date <- str_glue("{max(data$recovery_year)}-12-31") %>%
      as.Date()
    
    tibble(
      range = grp$range,
      date =  seq(start_date,end_date,by="month") %>%
        sample(30) %>%
        as.character()
    )
  }) %>%
  bind_rows()

# chlorophyll and sst datasets
datasets <- c(chl="MY1DMM_CHLORA", sst="MYD28M")


# url templates for getting image IDs and GeoTIFF data
id_url="https://neo.gsfc.nasa.gov/lib/imageSizes.php?datasetId={dataset}&curDate={timestamp}"
tiff_url="https://neo.gsfc.nasa.gov/servlet/RenderData?si={tiff_id}&cs=rgb&format=FLOAT.TIFF&width=3600&height=1800"

# get unix timestamp
ts <- function(d) {
  as.numeric(as.POSIXct(d,tz="HST"))
}

# retrieve a URL
get_url <- function(u) {
  r <- GET(u)
  if (r$status_code == 200) {
    return(content(r,as="text",encoding="UTF-8"))
  } else {
    return(NULL)
  }
}

# make sure directories exist
datasets %>%
  iwalk(\(ds,n) {
    dir_create(here("data","chlorophyll",n))
  })

# download tiffs to appropriate directories for unique dates
dates$date %>%
  unique() %>%
  walk(\(d) {
    datasets %>% 
      iwalk(\(dataset,ds_name) {
        timestamp <- ts(d)           # get timestamp
        idu <- str_glue(id_url)      # format url
        tiff_id <- get_url(idu)      # get tiff ID
        if (!is.null(tiff_id)) {     # sanity check
          tu <- str_glue(tiff_url)   # format url
          # download tiff
          download.file(tu,destfile = here("data","chlorophyll",ds_name,str_glue("{d}.tiff")))
        }
      })
  })

# cc$crabs_untransformed %>%
#   sample_tibble() %>%
#   st_as_sf(coords=c("lon","lat"),crs=4326) %>%
#   separate_wider_regex(sample,patterns=c(recovery_year="^[0-9]+","_",sample=".*"))  %>%
#   mutate(recovery_year=as.numeric(recovery_year)) %>%
#   dplyr::select(sample,everything()) %>%
#   unite("range",deployment_year,recovery_year,sep="--",remove=FALSE) %>%


bbox <- c(
  xmin = min(arms$lon),
  ymin = min(arms$lat),
  xmax = max(arms$lon),
  ymax = max(arms$lat)
)

in_box <- function(x,y,b) {
  (x >= b['xmin'] & x <= b['xmax']) & (y >= b['ymin'] & y <= b['ymax'])
}

get_raster <- function(f,lat,lon,buffer=300) {
  tryCatch(
    {
      rr <- raster(f)
      geo <- tibble(lat=lat,lon=lon) %>%
        st_as_sf(coords=c("lon","lat"),crs=4326)
      raster::extract(rr,geo,buffer=buffer) %>%
        as.numeric()
    },
    error = function(e) return(NA_real_)
  )
}

files <- dir_ls(
  here("data","chlorophyll"),
  regexp = "[0-9]{4}-[0-9]{2}-[0-9]{2}\\.tiff",
  recurse = TRUE) %>%
  map_dfr(\(fn) {
    d <- path_ext_remove(path_file(fn))
    m <- path_file(path_dir(fn))
    if (!str_detect(fn,"old_tiff")) {
      if (m %in% c("sst","chl")) {
        c(date = d, file = fn, measurement = m)
      }
    }
  })

dates <- files %>%
  pivot_wider(names_from = measurement,values_from = file, values_fill=NA)

# load sst data from rasters into a big data frame
raster_data <- dates %>%
  pmap_dfr(\(...) {
    row <- list(...)
    d <- row$date
    
    arms %>%
      dplyr::select(sample,deployment_year,recovery_year,lon,lat) %>%
      filter(year(d) >= deployment_year & year(d) <= recovery_year) %>%
      mutate(
        day=d,
        sst = get_raster(row$sst,lat,lon),
        chl = get_raster(row$chl,lat,lon),
        sst = na_if(sst,99999),
        chl = na_if(chl,99999)
      ) %>%
      dplyr::select(sample,day,sst,chl)
  })


rsum <- raster_data %>%
  group_by(sample) %>%
  summarise(sst=mean(sst,na.rm=T),chl=mean(chl,na.rm=T))

mm <- read_csv("data/metadata.csv")

aa <- arms %>%
  left_join(rsum,by="sample",suffix=c("","_sat")) %>%
  rename(unit=sample) %>%
  # unite("unit",recovery_year,sample) %>%
  dplyr::select(all_of(names(mm)),ends_with("_sat"))
  
write_csv(aa,here("data","metadata.csv"))


