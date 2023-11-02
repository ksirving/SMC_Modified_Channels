## get gage data for inventory

library(tidyverse)
library(tidylog)
library(sf)
library(mapview)
library(glue)
library(here)
library(lubridate)
library(beepr) # to tell us when stuff is done

devtools::install_github("USGS-R/nhdplusTools", force=T)
library(nhdplusTools)

## channel engineering data

BioEng <- read.csv("ignore/02_chan_eng.csv") %>%
  select(-c(X,channel_engineering_classification_date, channel_engineering_personnel, channel_engineering_comments)) %>%
  mutate(Class2 = ifelse(channel_engineering_class =="NAT", "Natural", "Modified"))

BioEng


# upload spatial gage data ------------------------------------------------


gagesCA <- st_read("/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/SpatialData/USGS_California_Stream_Gages/USGS_California_Stream_Gages.shp")

  

 # HUC12s
load("input_data/huc12_sf.rda") # CA h12s
# check size:
pryr::object_size(h12)

## map it

basemapsList <- c("Esri.WorldTopoMap", "Esri.WorldImagery",
                  "Esri.NatGeoWorldMap",
                  "OpenTopoMap", "OpenStreetMap", 
                  "CartoDB.Positron", "Stamen.TopOSMFeatures")

mapviewOptions(basemaps=basemapsList, fgb = FALSE)



# this map of all sites in same HUC 12
m1 <- mapview(gagesCA, cex=2, col.regions=c("orange"),
              layer.name="USGS Gages")


m1
# m1@map %>% leaflet::addMeasure(primaryLengthUnit = "meters")
# 
# mapshot(m1, url = paste0(getwd(), "/output_data/01_bio_sites_socal_counties_mapview.html"),
#         file = paste0(getwd(), "/ignore/01_bio_sites_socal_counties_mapview.png"))

getwd()
# Upload bio sites --------------------------------------------------------

bioSites <- st_read("input_data/01_bio_sites_all.shp")
head(bioSites)
st_crs(bioSites)

## bugs data
load(file = "/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/Cannabis_Eflows/ignore/SMC_bmi_cali_July2023.RData")
# csciScores <- read.csv("ignore/01_csci_comp_mets_comid_socal.csv")
head(bug_tax_ca)

bugSites <- bug_tax_ca %>%
  select(masterid, latitude, longitude, county, comid) %>%
  mutate(BioType= "BMI") %>%
  distinct()

### algae scores
load(file="/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/Cannabis_Eflows/ignore/SMC_algae_cali_July2023.RData")
head(alg_tax_ca)

algSites <- alg_tax_ca %>%
  select(masterid, latitude, longitude, county, comid) %>%
  mutate(BioType= "Algae") %>%
  distinct()

## join and make spatial
allSites <- bind_rows(bugSites, algSites) %>% distinct() %>%
  st_as_sf(coords=c("longitude", "latitude"), crs=4326, remove=F)
head(allSites)

# this map of all sites in same HUC 12
m1 <- mapview(gagesCA, cex=2, col.regions=c("orange"),
              layer.name="USGS Gages") +
  mapview(allSites, cex=2, col.regions=c("green"),
          layer.name="Bioassessment Sites")


m1


# HUC  --------------------------------------------------------------------

sf::sf_use_s2(FALSE) ## switch off spherical geom

## add HUC
## match CRS
h12 <- st_transform(h12, st_crs(gagesCA)) 
h12
# Add H12 to  Gages (adds ATTRIBUTES, retains ALL pts if left=TRUE)
all_gages_h12 <- st_join(gagesCA,left = TRUE, h12[c("HUC_12")])
all_gages_h12

#### bio
## match CRS
allSites <- st_transform(allSites, st_crs(gagesCA)) 

bio_h12 <- st_join(allSites,left = TRUE, h12[c("HUC_12")])
bio_h12

## make bio a df
bio_h12_df <- as.data.frame(bio_h12) %>% select(-geometry)

## lat, long are bio sites
## gemoetry is gage site

# now join based on H12: what  stations share same H12 as gage?

########### bio gages with same H12 -
## both coords for  and gages included (Lat, Lon & geometry = gages, latitude/longitude = bio)
bugs_gages_h12 <- inner_join(all_gages_h12, bio_h12_df, by = "HUC_12") %>% 
  distinct() %>% filter(!is.na(masterid))

bugs_gages_h12

# number of unique?
length(unique(factor(bugs_gages_h12$HUC_12))) # h12=929
length(unique(bugs_gages_h12$masterid)) # bug sites = 4663
length(unique(bugs_gages_h12$SiteNum)) # gages = 3428

# Subset to socal ---------------------------------------------------------

## upload shapefile
socal <- st_read("/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/SpatialData/SoCal_Counties/SoCal_Counties.shp")
socal <- st_transform(socal, crs=(st_crs(bugs_gages_h12)))

counties <- socal$NAME_PCASE

bugs_gages_h12 <- bugs_gages_h12 %>%
  filter(county %in% counties)
bugs_gages_h12

### make bio sites separate and spatial
## (i.e. change geometry to bug site)
bio_gages_h12 <- bugs_gages_h12 %>%
  select(SiteNum:BioType) %>% as.data.frame() %>% select(-geometry) %>% 
  st_as_sf(coords=c("longitude", "latitude"), crs=4326, remove=F) %>%
  rename(COMID = comid)

bio_gages_h12 ## has matched bio and gages in HUCs, geom id=s bio, sitename of gages for id

## same with gages
gages_bio_h12 <- bugs_gages_h12 %>%
  select(SiteName:masterid, comid) %>% as.data.frame() %>% select(-geometry) %>% 
  st_as_sf(coords=c("Lon", "Lat"), crs=4326, remove=F) %>%
  rename(COMID = comid)

gages_bio_h12

## same with HUC
mhucs <- unique(bugs_gages_h12$HUC_12)

h12Red <- h12 %>%
  filter(HUC_12 %in% mhucs)
h12Red
## map

# set background basemaps:
basemapsList <- c("Esri.WorldTopoMap", "Esri.WorldImagery",
                  "Esri.NatGeoWorldMap",
                  "OpenTopoMap", "OpenStreetMap", 
                  "CartoDB.Positron", "Stamen.TopOSMFeatures")

mapviewOptions(basemaps=basemapsList, fgb = FALSE)

# this map of all sites in same HUC 12
m1 <- mapview(bio_gages_h12, cex=6, col.regions="orange",
              layer.name="Bio Stations") +
  mapview(gages_bio_h12, col.regions="skyblue", cex=2, color="blue2",
          layer.name="USGS Gages") +
  mapview(h12Red, col.regions="dodgerblue", alpha.region=0.1,
          color="darkblue", legend=FALSE, layer.name="HUC12") 



m1
m1@map %>% leaflet::addMeasure(primaryLengthUnit = "meters")

# Check Ryan data ---------------------------------------------------------

# get ALL bug data (distinct stations)      
load("/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/bmi_ffm_links/data_output/01_bmi_stations_distinct.rda")
# algae_stations_distinct
load("/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/bmi_ffm_links/data_output/01a_algae_stations_distinct.rdata") 


bio_ffm<- read_rds("https://github.com/ryanpeek/flow_seasonality/blob/main/output/10_ffc_filtered_final_combined_rev.rds?raw=true")
bio_ffm


soHuc <- unique(bugs_gages_h12$HUC_12)# %in% unique(allRdata$HUC_12)
## filter by huc 12

bio_ffm <- bio_ffm %>%
  filter(HUC_12 %in% soHuc) #%>%
  # select(-Scores) %>%
  # distinct()

length(unique(bio_ffm$StationCode))
# asci (N=150)
asci_sites <- bio_ffm %>% filter(bioindicator=="ASCI") %>% 
  select(gageid:csci) %>% distinct(.keep_all=TRUE) %>% 
  left_join(., algae_stations_distinct) %>% 
  distinct(StationCode, .keep_all=TRUE) %>% 
  st_as_sf(coords=c("Longitude","Latitude"), crs=4269, remove=FALSE)

# csci (N=124)
csci_sites <- bio_ffm %>% filter(bioindicator=="CSCI") %>% 
  select(gageid:csci) %>% distinct(.keep_all=TRUE) %>% 
  left_join(., bmi_stations_distinct) %>% 
  distinct(StationCode, .keep_all=TRUE) %>% 
  st_as_sf(coords=c("longitude","latitude"), crs=4269, remove=FALSE)

# gages (N=103)
gages_sites <- bio_ffm %>% 
  select(gageid:csci) %>% distinct(gageid, .keep_all=TRUE) %>% 
  st_as_sf(coords=c("usgs_lon","usgs_lat"), crs=4269, remove=FALSE)

# HUCs
h12Red <- h12 %>%
  filter(HUC_12 %in% soHuc)
h12Red
## map

# set background basemaps:
basemapsList <- c("Esri.WorldTopoMap", "Esri.WorldImagery",
                  "Esri.NatGeoWorldMap",
                  "OpenTopoMap", "OpenStreetMap", 
                  "CartoDB.Positron", "Stamen.TopOSMFeatures")

mapviewOptions(basemaps=basemapsList, fgb = FALSE)

# this map of all sites in same HUC 12
m1 <- mapview(csci_sites, cex=6, col.regions="orange",
              layer.name="CSCI Sites") +
  mapview(asci_sites, cex=6, col.regions="green",
          layer.name="ASCI Sites") +
  mapview(gages_sites, col.regions="skyblue", cex=2, color="blue2",
          layer.name="USGS Gages") 
  # mapview(h12Red, col.regions="dodgerblue", alpha.region=0.1,
  #         color="darkblue", legend=FALSE, layer.name="HUC12") 



m1
m1@map %>% leaflet::addMeasure(primaryLengthUnit = "meters")

names(bio_ffm)

bio_ffm_red <- bio_ffm %>%
  select(comid, gageid, StationCode, bioindicator, REACHCODE:GNIS_NAME) %>%
  distinct()

write.csv(bio_ffm_red, "output_data/crosswalk_bio_gages.csv")

## join with mod channels

## channel engineering data

BioEng <- read.csv("ignore/02_chan_eng.csv") %>%
  select(-c(X,channel_engineering_classification_date, channel_engineering_personnel, channel_engineering_comments)) %>%
  mutate(Class2 = ifelse(channel_engineering_class =="NAT", "Natural", "Modified"))

BioEng
socaleng
socaleng <- left_join(bio_ffm_red, BioEng, by = c("StationCode" = "masterid"))

socaleng %>% group_by(channel_engineering_class, bioindicator) %>% tally() %>% as.data.frame() 

# channel_engineering_class bioindicator  n
# 1                         HB         ASCI 23
# 2                         HB         CSCI 29
# 3                        NAT         ASCI 85
# 4                        NAT         CSCI 66
# 5                        SB0         ASCI 12
# 6                        SB0         CSCI  9
# 7                        SB1         ASCI 13
# 8                        SB1         CSCI 10
# 9                        SB2         ASCI 13
# 10                       SB2         CSCI  6

# Measure distance --------------------------------------------------------



# Upstream flowlines from Gage --------------------------------------------
## throwing errors and very slow!!!!!!

## TRANSFORM TO UTM datum for flowlines
h12Red <- st_transform(h12Red, crs = 3310) ## huc 12
bio_gages_h12 <- st_transform(bio_gages_h12, crs = 3310) ## bio sites
gages_bio_h12 <- st_transform(gages_bio_h12, crs = 3310) ## gages

# use a list of comids to make a list to pass to the nhdplusTools function
# first do gages, then see how many bio sites are on lines, then same with transects 
str(gages_bio_h12)
gages_bio_h12 <- gages_bio_h12 %>%
  mutate(COMID = as.numeric(COMID))

# Use the gage com_list
coms_list <- map(gages_bio_h12$COMID, ~list(featureSource = "COMID", featureID=.x))
coms_list_test <- coms_list[1:5]
coms_list_test
# check
coms_list[3] # should list feature source and featureID


# Get upstream mainstem streamlines (10 km limit) from gages
# coms_list can be a dataframe (then may need to change `coms_list` to `coms_list$comid` or just a list of comids. 
mainstemsUS <- map(coms_list_test, ~navigate_nldi(nldi_feature = .x,
                                             mode="UM", # upstream main
                                             distance_km = 10))
mainstemsUS_c
# transform the sf layer to match mainstems crs (4326)
gages_bio_h12 <- gages_bio_h12 %>% st_transform(4326)

# check length (for NAs?)
mainstemsUS %>% 
  purrr::map_lgl(~ length(.x)>1) %>% table()

# drop NA/empty elements
mainstemsUS_c <- mainstemsUS %>% purrr::compact()
# mainstemsUS_c 
# testList <- mainstemsUS_c[172]

# testList
class(mainstemsUS_c)
# make a single flat layer
# this accesses each item in our list of items...nhdplus returns a list of 2 dataframes that include $UM_flowlines, $origin. 
# the name UM_flowlines can change depending on the "mode" above (may be DS or DD_flowlines).

mainstems_flat_us <- map_df(mainstemsUS_c, ~mutate(.x$UM_flowlines, comid_origin=.x$origin$comid, .after=nhdplus_comid))

head(mainstems_flat_us)

length(unique(mainstems_flat_us$nhdplus_comid)) ## 974
length(unique(mainstems_flat_us$comid_origin)) ## 215

# bind together
mainstems_us <- sf::st_as_sf(mainstems_flat_us, use.names = TRUE, fill = TRUE)

# add direction to gage col
mainstems_us <- mainstems_us %>% 
  mutate(from_gage = "UM")

length(unique(mainstems_us$nhdplus_comid)) ## 974
length(unique(mainstems_us$comid_origin)) ## 215

# rm temp files
rm(mainstems_flat_us, mainstemsUS)

st_crs(mainstems_us)

# set background basemaps:
basemapsList <- c("Esri.WorldTopoMap", "Esri.WorldImagery",
                  "Esri.NatGeoWorldMap",
                  "OpenTopoMap", "OpenStreetMap", 
                  "CartoDB.Positron", "Stamen.TopOSMFeatures")

mapviewOptions(basemaps=basemapsList, fgb = FALSE)

# this map of all sites in same HUC 12 and upstream 10km
m1 <- mapview(mainstems_us, color = "navyblue") +
  mapview(bio_gages_h12, cex=6, col.regions="orange",
              layer.name="Bio Stations") +
  mapview(gages_bio_h12, col.regions="skyblue", cex=2, color="blue2",
          layer.name="USGS Gages") +
  mapview(h12Red, col.regions="dodgerblue", alpha.region=0.1,
          color="darkblue", legend=FALSE, layer.name="HUC12") 



m1
m1@map %>% leaflet::addMeasure(primaryLengthUnit = "meters")


