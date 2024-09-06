library(tidyverse)
library(sf)
library(raster)
library(here)
library(lubridate)
library(httr)
library(fs)
library(furrr)



# set random seed for consistency
set.seed(31337)

# setup parallel workers
if (rstudioapi::isAvailable()) {
  plan(multisession,workers=8)
} else {
  plan(multicore,workers=8)
}

# output directory
outdir <- here("output","modis")
dir_create(outdir)

# create monthly directory
monthly <- here("data","chlorophyll","monthly")
dir_create(monthly)


if (!exists("skip_dl")) {
  skip_dl <<- TRUE
}
# download code -----------------------------------------------------------
library(httr)
library(raster)

if (!skip_dl) {
    # netrc file must exist with login credentials to nasa data warehouse
  netrc <- here("data","chlorophyll","netrc")
  # load data urls
  # urls <- dir_ls(here("data","chlorophyll"),regex=r'(monthly_9km.*\.txt)') %>%
  urls <- dir_ls(here("data","chlorophyll"),regex=r'(nasa_.*\.txt)') %>%
    map(~read_lines(.x)) %>%
    unlist()
    # download files in parallel
  # show progress bar
  urls %>%
    future_walk(\(url) {
      GET(url, write_disk(here("data","chlorophyll","monthly",path_file(url)), overwrite = TRUE), 
          config(netrc = TRUE, netrc_file = netrc, httpauth=1), set_cookies("LC" = "cookies"))
    },.options = furrr_options(seed = 31337),.progress = TRUE)
}



# do calculations -------------------------------------------------------------------------------------------------

yrs <- seq(2009,2016)
resolutions <- c("4km", "9km")

# average years for chl-A
yearmeans_chl <- resolutions %>%
  set_names() %>%
  map(\(res) {
    cat(str_glue("\naveraging years for {res} chl-A\n"))
    yrs %>%
      set_names() %>%
      future_map(\(yr) {
        chl <- dir_ls(monthly,regex=str_glue("AQUA_MODIS\\.{yr}.*CHL.*\\.{res}\\.nc")) %>%
          map(\(nc) raster(nc,varname="chlor_a"))
        s <- stack(chl)
        ml <- calc(s,mean,na.rm=T) 
        return(ml)
      },.progress=TRUE,.options = furrr_options(seed=31337))
    })
cat("\nsaving averaged years\n")

# save files
yearmeans_chl %>%
  iwalk(\(res,rr) {
    res %>%
      future_iwalk(\(r,yr) {
        # filename <- str_glue("output/{yr}_mean_chl.tif")
        filename <- here(outdir,str_glue("{yr}_mean_chl_{rr}.tif"))
        writeRaster(r,filename,format="GTiff",overwrite=TRUE)
      },.progress=TRUE,.options = furrr_options(seed=31337))
  })

cat("\ncalculating average for all years and saving\n")

# average all years
yearmeans_chl %>%
  future_iwalk(\(res,r) {
    y <- stack(res)
    allyears <- calc(y,mean,na.rm=T)
    writeRaster(allyears,here(outdir,str_glue("allyears_chl_{r}.tif")),format="GTiff",overwrite=TRUE)
  },.progress=TRUE,.options = furrr_options(seed=31337))

# average years for SSt
yearmeans_sst <- resolutions %>%
  set_names() %>%
  map(\(res) {
    cat(str_glue("\naveraging years for {res} SST\n"))
    yrs %>%
      set_names() %>%
      future_map(\(yr) {
        sst <- dir_ls(monthly,regex=str_glue("AQUA_MODIS\\.{yr}.*SST.*\\.{res}\\.nc")) %>%
          map(\(nc) raster(nc,varname="sst"))

        s <- stack(sst)
        ml <- calc(s,mean,na.rm=T) 
        return(ml)
      },.progress=TRUE,.options = furrr_options(seed=31337))
  })

cat("\nsaving averaged years\n")

# save files
yearmeans_sst %>%
  iwalk(\(res,rr) {
    res %>%
      future_iwalk(\(r,yr) {
        filename <- here(outdir,str_glue("{yr}_mean_sst_{rr}.tif"))
        writeRaster(r,filename,format="GTiff",overwrite=TRUE)
    },.progress=TRUE,.options = furrr_options(seed=31337))
  })

cat("\ncalculating average for all years and saving\n")
# average all years
yearmeans_sst %>%
  future_iwalk(\(res,r) {
    y <- stack(res)
    allyears <- calc(y,mean,na.rm=T)
    writeRaster(allyears,here(outdir,str_glue("allyears_sst_{r}.tif")),format="GTiff",overwrite=TRUE)
  },.progress=TRUE,.options = furrr_options(seed=31337))
cat("\n")


# get satellite-derived data --------------------------------------------------------------------------------------

# load metadata and get deployment range
arms <- read_csv(here("data","metadata.csv")) %>%
  rename(sample=unit) %>%
  unite("range",deployment_year,recovery_year,sep="--",remove=FALSE) 

ranges <- arms %>%
  distinct(range) %>%
  pull(range)

# calculate means across deployment year ranges for sst
range_means_sst <- arms %>%
  distinct(deployment_year,recovery_year) %>%
  mutate(across(everything(),as.numeric)) %>%
  array_branch(margin=1) %>%
  set_names(ranges) %>%
  map(\(row) {
    start <- row[1]
    end <- row[2]
    yrs <- seq(start,end)
    resolutions %>%
      set_names() %>%
      map(\(res) {
        files <- here(outdir,str_glue("{yrs}_mean_sst_{res}.tif"))
        r <- files %>%
          future_map(\(r) raster(r,varname="sst"),.options = furrr_options(seed = 31337))
        s <- stack(r)
        calc(s,mean,na.rm=T)
      })
  },.progress=TRUE)

  # calculate means across deployment year ranges for chl
range_means_chl <- arms %>%
  distinct(deployment_year,recovery_year) %>%
  mutate(across(everything(),as.numeric)) %>%
  array_branch(margin=1) %>%
  set_names(ranges) %>%
  map(\(row) {
    start <- row[1]
    end <- row[2]
    yrs <- seq(start,end)
    resolutions %>%
      set_names() %>%
      map(\(res) {
        files <- here(outdir,str_glue("{yrs}_mean_chl_{res}.tif"))
        r <- files %>%
          future_map(\(r) raster(r,varname="chlor_a"),.options = furrr_options(seed = 31337))
        s <- stack(r)
        calc(s,mean,na.rm=T)
      })
  },.progress=TRUE)

# create sf object for arms deployment sites
arms_sf <- read_csv("data/metadata.csv") %>%
  st_as_sf(coords=c("lon","lat"),crs=4326,remove=FALSE)
  

# modify arms metadata with new chl/sst values
aa <- arms %>%
  group_by(range) %>%
  group_modify(\(data,grp) {
    srast <- range_means_sst[[grp$range]] 
    crast <- range_means_chl[[grp$range]]
    f <- arms_sf %>%
      filter(unit %in% data$sample)
    tdata <- srast %>%
      map(\(r) {
        r %>%
          raster::extract(f,buffer=1000) %>%
          map_dbl(~mean(.x,na.rm=T))
      })
    cdata <- crast %>%
      map(\(r) {
        r %>% 
          raster::extract(f,buffer=1000) %>%
          map_dbl(~mean(.x,na.rm=T))
      })
    data %>%
      mutate(
        sst_4km = tdata$`4km`, sst_9km = tdata$`9km`,
        chl_4km = cdata$`4km`, chl_9km = cdata$`9km`
      )
      # mutate(sst_new = tdata,chl_new = cdata)
  }) %>%
  ungroup()

# save new metadata
aa %>%
  rename(unit=sample) %>%
  dplyr::select(-range) %>%
  write_csv(here("data","metadata.csv"))


