library(tidyverse)
library(sf)
library(raster)
library(here)
library(lubridate)
library(httr)
library(fs)

# set random seed for consistency
set.seed(31337)

get_raster <- function(f,lat,lon,buffer=300,multi=c("first","mean"),naval=99999) {
  tryCatch(
    {
      rr <- raster(f)
      NAvalue(rr) <- naval
      geo <- tibble(lat=lat,lon=lon) %>%
        st_as_sf(coords=c("lon","lat"),crs=4326)
      rr <- raster::extract(rr,geo,buffer=buffer)
      switch(
        match.arg(multi),
        first = rr %>%
          map_dbl(~na.omit(.x)[1]),
        mean = rr %>%
          map_dbl(~mean(.x,na.rm=TRUE))
      )      
    },
    error = function(e) return(NA_real_)
  )
}

# pre-parse metadata and get deployment & recovery years as well as deployment range
arms <- read_csv(here("data","metadata.csv")) %>%
  rename(sample=unit) %>%
  separate_wider_regex(sample,patterns=c(recovery_year="^[0-9]+","_",sample=".*"),cols_remove=FALSE)  %>%
  mutate(recovery_year=as.numeric(recovery_year)) %>%
  dplyr::select(sample,everything()) %>%
  unite("range",deployment_year,recovery_year,sep="--",remove=FALSE) 

ranges <- arms %>%
  distinct(range) %>%
  pull(range)

range_means_sst <- arms %>%
  distinct(deployment_year,recovery_year) %>%
  array_branch(margin=1) %>%
  set_names(ranges) %>%
  map(\(row) {
    start <- row[1]
    end <- row[2]
    yrs <- seq(start,end)
    files <- here("data","chlorophyll","means",str_glue("{yrs}_mean_sst.tif"))
    r <- files %>%
      map(raster)
    s <- stack(r)
    calc(s,mean,na.rm=T)
  },.progress=TRUE)

range_means_chl <- arms %>%
  distinct(deployment_year,recovery_year) %>%
  array_branch(margin=1) %>%
  set_names(ranges) %>%
  map(\(row) {
    start <- row[1]
    end <- row[2]
    yrs <- seq(start,end)
    files <- here("data","chlorophyll","means",str_glue("{yrs}_mean_chl.tif"))
    r <- files %>%
      map(raster)
    s <- stack(r)
    calc(s,mean,na.rm=T)
  },.progress=TRUE)

arms_sf <- read_csv("data/metadata.csv") %>%
  st_as_sf(coords=c("lon","lat"),crs=4326,remove=FALSE)
  

# chl_data <- raster("data/chlorophyll/oracle/Present.Surface.Chlorophyll.Mean.tif") 

aa <- arms %>%
  group_by(range) %>%
  group_modify(\(data,grp) {
    srast <- range_means_sst[[grp$range]] 
    crast <- range_means_chl[[grp$range]]
    f <- arms_sf %>%
      filter(unit %in% data$sample)
    tdata <- srast %>%
      raster::extract(f,buffer=1000) %>%
      map_dbl(~mean(.x,na.rm=T))
    cdata <- crast %>%
      raster::extract(f,buffer=1000) %>%
      map_dbl(~mean(.x,na.rm=T))
    data %>%
      mutate(sst_new = tdata,chl_new = cdata)
  }) %>%
  ungroup()


mm <- read_csv(here("data","metadata.csv"))

f <- aa %>%
  rename(unit=sample) %>%
  dplyr::select(all_of(names(mm)),ends_with("_new")) 

f %>%
  write_csv(here("data","metadata.csv"))


