# prep/load libraries ---------------------------------------------------------------------------------------------
library(tidyverse)
library(sf)
library(here)
library(lubridate)
library(httr)
library(fs)
library(furrr)
library(httr)
library(terra)
source("extract_raster.R")

# set random seed for consistency
set.seed(31337)

# output directory
outdir <- here("output","modis")
dir_create(outdir)

# create monthly directory
monthly <- here("data","chlorophyll","monthly")
dir_create(monthly)

# human impacts directory/file
impacts_data <- here("data","impact")
impacts_dir <- here("output","impact")
dir_create(impacts_data)
dir_create(impacts_dir)


if (!exists("skip_dl")) {
  skip_dl <<- TRUE
}

# load metadata and get deployment range
arms <- read_csv(here("data","metadata.csv")) %>%
  rename(sample=unit) %>%
  unite("range",deployment_year,recovery_year,sep="--",remove=FALSE) 

ranges <- arms %>%
  distinct(range) %>%
  pull(range)

# create sf object for arms deployment sites
arms_sf <- read_csv("data/metadata.csv") %>%
  st_as_sf(coords=c("lon","lat"),crs=4326,remove=FALSE)
  

# initial point search radius
initial_radius <- 10

# download chl/sst datasets ---------------------------------------------------------------------------------------
# setup parallel workers
if (rstudioapi::isAvailable()) {
  plan(multisession,workers=8)
} else {
  plan(multicore,workers=8)
}

if (!skip_dl) {
  # netrc file must exist with login credentials to nasa data warehouse
  netrc <- here("data","chlorophyll","netrc")
  # load data urls
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
          rast(subds="chlor_a") %>%
          app(mean,na.rm=T)
        names(chl) <- "chl"
        return(wrap(chl))
      },.options = furrr_options(seed=31337))
    },.progress=TRUE)

cat("\nsaving averaged years\n")

# save files
yearmeans_chl %>%
  iwalk(\(res,rr) {
    res %>%
      iwalk(\(r,yr) {
        filename <- here(outdir,str_glue("{yr}_mean_chl_{rr}.tif"))
        writeRaster(unwrap(r),filename,filetype="GTiff",overwrite=TRUE)
        gc()
      },.progress=TRUE)
  },.progress=TRUE)

cat("\ncalculating average for all years and saving\n")

# average all years
yearmeans_chl %>%
  iwalk(\(res,r) {
    y <- res %>%
      map(unwrap) %>%
      rast()
    allyears <- app(y,mean,na.rm=T)
    writeRaster(allyears,here(outdir,str_glue("allyears_chl_{r}.tif")),filetype="GTiff",overwrite=TRUE)
    gc()
  },.progress=TRUE)

# average years for SSt
yearmeans_sst <- resolutions %>%
  set_names() %>%
  map(\(res) {
    cat(str_glue("\naveraging years for {res} SST\n"))
    yrs %>%
      set_names() %>%
      future_map(\(yr) {
        sst <- dir_ls(monthly,regex=str_glue("AQUA_MODIS\\.{yr}.*SST.*\\.{res}\\.nc")) %>%
          rast(subds="sst") %>%
          app(mean,na.rm=T)
        names(sst) <- "sst"
        return(wrap(sst))
      },.progress=TRUE,.options = furrr_options(seed=31337))
  })

cat("\nsaving averaged years\n")

# save files
yearmeans_sst %>%
  iwalk(\(res,rr) {
    res %>%
      iwalk(\(r,yr) {
        filename <- here(outdir,str_glue("{yr}_mean_sst_{rr}.tif"))
        writeRaster(unwrap(r),filename,filetype="GTiff",overwrite=TRUE)
        gc()
    },.progress=TRUE)
  },.progress=TRUE)

cat("\ncalculating average for all years and saving\n")
# average all years
yearmeans_sst %>%
  iwalk(\(res,r) {
    y <- res %>%
      map(unwrap) %>%
      rast()
    allyears <- app(y,mean,na.rm=T)
    writeRaster(allyears,here(outdir,str_glue("allyears_sst_{r}.tif")),filetype="GTiff",overwrite=TRUE)
  },.progress=TRUE)
cat("\n")

rm(yearmeans_sst)
rm(yearmeans_chl)
gc()



# integrate chl/sst data ------------------------------------------------------------------------------------------

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
      future_map(\(res) {
        files <- here(outdir,str_glue("{yrs}_mean_sst_{res}.tif"))
         rast(files) %>%
          app(mean,na.rm=T) %>%
          setNames("sst") %>%
          wrap()
      },.options = furrr_options(seed = 31337))
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
      future_map(\(res) {
        files <- here(outdir,str_glue("{yrs}_mean_sst_{res}.tif"))
        rast(files) %>%
         app(mean,na.rm=T) %>%
         setNames("chl") %>%
         wrap()
      },.options = furrr_options(seed = 31337))
  },.progress=TRUE)


# modify arms metadata with new chl/sst values
aa <- arms %>%
  group_by(range) %>%
  group_modify(\(data,grp) {
    srast <- range_means_sst[[grp$range]] 
    crast <- range_means_chl[[grp$range]]
    arms_subset <- arms_sf %>%
      filter(unit %in% data$sample)
    tdata <- srast %>%
      map(\(r) {
        f <- r %>%
          unwrap() %>%
          nearExtract(arms_subset,b=initial_radius)
      }) %>%
      list_cbind()
    cdata <- crast %>%
      map(\(r) {
        r %>% 
          unwrap() %>%
          nearExtract(arms_subset,b=initial_radius)
      }) %>%
      list_cbind()
    data %>%
      mutate(
        sst_4km = tdata$`4km`, sst_9km = tdata$`9km`,
        chl_4km = cdata$`4km`, chl_9km = cdata$`9km`
      )
  }) %>%
  ungroup()

# save new metadata
aa %>%
  rename(unit=sample) %>%
  dplyr::select(-range) %>%
  write_csv(here("data","metadata.csv"))



# download human impacts dataset ----------------------------------------------------------------------------------

# data DOI: doi:10.5063/F12B8WBS
# source of raster data: https://oceanhealthindex.org/resources/data/cumulative-human-impacts/

# urls of impacts for each year
impacts <- list(
  `2009` = "https://knb.ecoinformatics.org/knb/d1/mn/v2/object/urn%3Auuid%3A5a1c71ec-a441-4592-9021-e05308f0e525",
  `2010` = "https://knb.ecoinformatics.org/knb/d1/mn/v2/object/urn%3Auuid%3Ac2fe2f40-0e48-4575-9f0b-572b29245e83",
  `2011` = "https://knb.ecoinformatics.org/knb/d1/mn/v2/object/urn%3Auuid%3Ad3534e9b-b906-4e22-bbd6-7d8fdadee142",
  `2012` = "https://knb.ecoinformatics.org/knb/d1/mn/v2/object/urn%3Auuid%3Ac6a1b528-f957-4225-91d4-84cfba7696a8",
  `2013` = "https://knb.ecoinformatics.org/knb/d1/mn/v2/object/urn%3Auuid%3Ac2fe2f40-0e48-4575-9f0b-572b29245e83"
)
# years we have data for
all_years <- names(impacts)

# load metadata and get deployment range
arms <- read_csv(here("data","metadata.csv")) %>%
	rename(sample=unit) %>%
	unite("range",deployment_year,recovery_year,sep="--",remove=FALSE) 

# get distinct year ranges
ranges <- arms %>%
	distinct(range) %>%
	pull(range)

# download impact rasters in parallel
if (!skip_dl) {
	impacts %>%
		future_iwalk(\(url,yr) {
			GET(url,write_disk(here(impacts_data,str_glue("impacts_{yr}.tif")),overwrite = TRUE))
		},.options = furrr_options(seed = 31337),.progress=TRUE)
}


# calculate human impacts means -----------------------------------------------------------------------------------

# save means for ranges (serially because they're big and will use a lot of memory)
arms %>% 
	distinct(deployment_year,recovery_year) %>%
	array_branch(1) %>%
	walk(\(yrs) {
		if (all(yrs %in% all_years)) {
		  # get sequence of years
			yy <- seq(min(yrs),max(yrs))
			# load rasters and calculate means
			files <- path(impacts_data,str_glue("impacts_{yy}.tif")) 
			rr <- rast(files) %>%
				app(mean,na.rm=T) %>%
				setNames("impacts")
			# save new raster and do garbage collection to free up memory
			writeRaster(rr,path(impacts_dir,str_glue("impacts_{min(yrs)}_{max(yrs)}.tif")),overwrite=TRUE,filetype="GTiff")
			gc()
		}
	},.progress=TRUE)

# calculate means for all years, because we don't have data beyond 2013
files <- path(impacts_data,str_glue("impacts_{all_years}.tif")) 
rr <- rast(files) %>%
	app(mean,na.rm=T) %>%
	setNames("impacts")
writeRaster(rr,path(impacts_dir,str_glue("impacts_{min(all_years)}_{max(all_years)}.tif")),overwrite=TRUE,filetype="GTiff")
gc()


# integrate human impacts data ------------------------------------------------------------------------------------
all_years <- seq(2009,2013)
aa <- arms %>%
  # modify each date range
  group_by(range) %>%
  group_modify(\(data,grp) {
    # get date range
    y1 <- min(data$deployment_year)
    y2 <- max(data$recovery_year)
    # if our range is outside of the available range,
    # take the mean across the whole time period we have data for
    if (!all(c(y1,y2) %in% all_years)) {
      y1 <- min(all_years)
      y2 <- max(all_years)
    }
    ff <- path(impacts_dir,str_glue("impacts_{y1}_{y2}.tif"))
    if (file_exists(ff)) {
      # load the raster
      impacts <- rast(ff)
      # convert the current grouped data to an sf object and extract
      # the impacts, sequentially increasing radius if we don't find it right away
      asf <- data %>%
        st_as_sf(coords=c("lon","lat"),crs=4326,remove=FALSE) %>%
        st_transform(st_crs(impacts))
      iid <- nearExtract(impacts,asf,b=initial_radius) 
      data %>%
        mutate(cumulative_impacts = iid$impacts, impact_radius = iid$radius)
    } else {
      return(data %>% mutate(cumulative_impacts = NA_real_))
    }
  }) %>%
  ungroup()

aa <- aa %>%
  group_by(island) %>%
  mutate(human_impact = mean(cumulative_impacts)) %>%
  ungroup() %>%
  select(-range,unit=sample) %>%
  select(unit,everything())

write_csv(aa,"data/metadata.csv")

