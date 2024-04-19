## level 2: site specific analysis,  how many of mod chans have issues with flow metric, keep 4 outcomes 

library(gam)
library(tidylog)
library(tidyverse)
library(sf)

## directory for figures
out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/SMC_Modified_Channels/Figures/"

## workflow
## use 4 outcome categories to count how many mod channels have issues with flow metrics
## which flow metrics - visualise
## which flow metrics also issue with bio - visualise


# Upload Data -------------------------------------------------------------

## within limits doc
imps <- read.csv("output_data/03_impact_ffm_bio.csv") %>%
  select( -X)

## tally of categories
tal <- read.csv("output_data/03_count_impact.csv")


# Format and explore ------------------------------------------------------

head(imps)

## NA's in result when ffm limits columns are NA
## Rules to improve:
## if ffm column NA, bio is within limits = indeterminant bio/hydro
## if ffm column NA, bio is not within limits = bio impact, indetermiant bio/hydro

impx <- imps %>%
  #select(-c(CorObs_Pred, sampleyear, Lower, Upper, longitude:comid )) %>%
  drop_na(IndexValue) %>% ## remove NAs as all sites should have a score
  filter(!Threshold == "SB1") %>%
  mutate(Result2 = case_when(is.na(WithinHydroLimits) & WithinBioLimits == "NotWithin"  ~ "BioImpact_HydroIndeterminant", ## bio is impacted but hydro is NA
         is.na(WithinHydroLimits) & WithinBioLimits == "Within" ~ "HydroIndeterminant")) %>% ## Bio not impacted and hydro is NA
  mutate(Result = ifelse(is.na(Result), Result2, Result))

dim(impx)
## where are the NAs in bio
# ind <- which(is.na(impx$WithinBioLimits))
# ind
# 
# test <-impx[ind,] ## all SB1 = remove as not needed


## Tally Result categories

## tally of impact per ffm - with the indeterminant categories

tallyImpact <- impx %>%
  group_by(Index, Hydro_endpoint, Flow.Metric.Name,Threshold, Result, SmoothingFunc) %>%
  distinct() %>%
  tally() %>%
  drop_na(Result) %>%
  mutate(PercChans = (n/sum(n)*100))

tallyImpact

## same as above but assign only the 4 categories - make sure to note the caveat

## NA's in result when ffm limits columns are NA
## Rules to improve:
## if ffm column NA, other is within limits = NA
## if ffm column NA, other is not within limits = bio impact

impx1 <- imps %>%
  #select(-c(CorObs_Pred, sampleyear, Lower, Upper, longitude:comid )) %>%
  drop_na(IndexValue) %>%
  filter(!Threshold == "SB1") %>%
  mutate(Result2 = case_when(is.na(WithinHydroLimits) & WithinBioLimits == "NotWithin"  ~ "BioImpact", ## bio is impacted but hydro is NA
                             is.na(WithinHydroLimits) & WithinBioLimits == "Within" ~ NA)) %>%
  mutate(Result = ifelse(is.na(Result), Result2, Result))

names(impx1)

## Tally Result categories

## tally of impact per ffm - with the indeterminant categories

tallyImpact1 <- impx1 %>%
  group_by(Index, Hydro_endpoint, Flow.Metric.Name,Threshold, Result, SmoothingFunc) %>%
  distinct() %>%
  tally() %>%
  drop_na(Result) %>%
  mutate(PercChans = (n/sum(n)*100))

tallyImpact1

## NAs in hydro limit come from when the curve doesn't reach the bio limit?
## so no flow targets

## same as above but assign only the 4 categories - make sure to note the caveat

## NA's in result when ffm limits columns are NA
## Rules to improve:
## if ffm column NA, indeterminant

impx2 <- imps %>%
  #select(-c(CorObs_Pred, sampleyear, Lower, Upper, longitude:comid)) %>%
  drop_na(IndexValue) %>%
  filter(!Threshold == "SB1") %>%
  mutate(Result2 = case_when(is.na(WithinHydroLimits)   ~ "Indeterminant")) %>% ## bio is impacted but hydro is NA
  mutate(Result = ifelse(is.na(Result), Result2, Result))

impx2
names(impx2)
## Tally Result categories

## tally of impact per ffm - with the indeterminant categories

tallyImpact2 <- impx2 %>%
  group_by(Index, Hydro_endpoint, Flow.Metric.Name,Threshold, Result, SmoothingFunc) %>%
  distinct() %>%
  tally() %>%
  drop_na(Result) %>%
  mutate(PercChans = (n/sum(n)*100))

class(tallyImpact2$Threshold)


# Clean Tables ------------------------------------------------------------
unique(tallyImpact2$Threshold)

## change order of cols, rename cols, change names of modified streams calss
FinalTablex <- tallyImpact1 %>% ## using option with least NAs
  ungroup() %>%
  select(Index, Flow.Metric.Name, SmoothingFunc, Threshold, Result, PercChans, n) %>%
  rename(FlowMetric = Flow.Metric.Name, ModifiedClass = Threshold, Impact = Result, NumberOfSites = n, PercentageOfSites = PercChans) %>%
  mutate(ModifiedClass = case_when((ModifiedClass == "HB") ~ "Hard Bottom",
                                      (ModifiedClass == "NAT") ~ "Natural",
                                      (ModifiedClass == "NATHigh") ~ "Overall 1st",
                                      (ModifiedClass == "NATLow") ~ "Overall 30th",
                                      (ModifiedClass == "NATMed") ~ "Overall 10th",
                                      (ModifiedClass == "SB0") ~ "Soft Bottom (no hard sides)",
                                      (ModifiedClass == "SB2") ~ "Soft Bottom (two hard sides)"))
         
FinalTablex    
## get total sites per class
sums <- FinalTablex %>%
  group_by(Index, FlowMetric, ModifiedClass, SmoothingFunc) %>%
  summarise(NumberSitesPerClass = sum(NumberOfSites))

sums
  
## pivot impact results wider
FinalTable <- FinalTablex %>%  
  select(-NumberOfSites) %>%
  filter(SmoothingFunc == 6) %>% ## change later!
  pivot_wider(names_from = Impact, values_from = PercentageOfSites) %>%
  mutate(across(everything(), .fns = ~replace_na(.,0))) %>%
  select(Index:ModifiedClass, NoImpact, BioImpact, HydroImpact, BothImpact)

## join number of sites

FinalTable <- full_join(FinalTable, sums, by = c("Index", "FlowMetric", "ModifiedClass", "SmoothingFunc"))
FinalTable

## save out

write.csv(FinalTable, "output_data/03_percent_impacts_each_Class.csv")

## table with columns as channel class

FinalTableWide <- FinalTablex %>%
  select(-NumberOfSites) %>%
  pivot_wider(names_from = ModifiedClass, values_from = PercentageOfSites) %>%
  mutate(across(everything(), .fns = ~replace_na(.,0))) #%>%
  # select(-)

  FinalTableWide
  
  write.csv(FinalTableWide, "output_data/03_percent_impacts_each_Class_wide.csv")

# Map results -------------------------------------------------------------

## upload boundary

socal <- st_read("/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/SpatialData/SoCal_Counties/SoCal_Counties.shp")
socal

plot(socal)

## data per site

head(impx2)

## make spatial

imps_sf <- impx1 %>% 
  # st_as_sf(coords=c( "longitude", "latitude"), crs=4326, remove=F) %>%
  # st_transform(crs = st_crs(socal)) %>%
  select(Index, Hydro_endpoint, Threshold, BioThresh, SmoothingFunc, masterid, COMID, Flow.Metric.Name, Flow.Component, Result, longitude, latitude)  %>%
  mutate(Threshold = as.factor(Threshold), Result = as.factor(Result))

## for now remove indeterminant
imps_sf <- imps_sf %>% 
  filter(!Result == "Indeterminant", Index == "csci") 
imps_sf

## use ggmap to get google 
library(ggmap)
library("ggsci")

## register google key 
## in notes

# Philadelphia Lat 39.95258 and Lon is -75.16522
# basemap <- get_map(location=c(lon = -119.2742, lat = 34.41183), zoom=11, maptype = 'terrain-background', source = 'stamen')
# basemap <- ggmap(get_googlemap(center = c(lon = -117.803056, lat = 33.595774), zoom=8, maptype = "terrain", color = "color" ))
# ## c("terrain", "satellite", "roadmap", "hybrid")
# ?get_googlemap
# print(basemap)

# save(basemap, file = "Figures/03_basemap_socal_ggplot.RData")

load(file = "Figures/03_basemap_socal_ggplot.RData")

# Plot per flow metric and result

metric <- unique(imps_sf$Hydro_endpoint)
metric

## defione impact acategory
results <- unique(imps_sf$Result)
r=2
m=2
results

imps_sf
  
  for(m in 1:length(metric)) {
    
    imps_sf1 <- imps_sf %>%
      filter(Hydro_endpoint == metric[m], SmoothingFunc == 6 )
    
    # for(r in 1:length(results)) {
      
      ## extract - impact result
      # imps_sfx <- imps_sf1 %>%
      #   filter(Result == results[r])
    
    m1 <- basemap + 
      geom_point(data = imps_sf1, aes(x = longitude, y = latitude, colour = Threshold), size = 5) +
      scale_color_discrete(name = "Modified Class") +
      guides(size = "none") +
      facet_wrap(~Result)
    m1
    
    file.name1 <- paste0("Figures/03_", metric[m],"_", results[r], "_map_4_cats_per_impact.jpg")
    file.name1
    ggsave(m1, filename=file.name1, dpi=500, height=10, width=15)
    
  }
  
}

# Plot per flow metric and channel type

metric <- unique(imps_sf$Hydro_endpoint)
metric

## defione impact acategory
type <- unique(imps_sf$Threshold)
r=2
m=2
results

for(m in 1:length(metric)) {
  
  imps_sf1 <- imps_sf %>%
    filter(Hydro_endpoint == metric[m], SmoothingFunc == 6)
  
  # for(r in 1:length(results)) {
  #   
  #   ## extract - impact result
  #   imps_sfx <- imps_sf1 %>%
  #     filter(Result == results[r])
  #   
    m1 <- basemap + 
      geom_point(data = imps_sf1, aes(x = longitude, y = latitude, colour = Result), size = 5) +
      scale_color_jco(name = "Impact Type") + ## colours re-order
      guides(size = "none") +
      facet_wrap(~Threshold)
    
    m1
    file.name1 <- paste0("Figures/03_", metric[m], "_map_4_cats_per_channel_type.jpg")
    file.name1
    ggsave(m1, filename=file.name1, dpi=500, height=10, width=15)
    
  
  
}


## plot
# m1 <- ggplot() +
#   geom_sf(data = socal) +
#   geom_sf(data = imps_sfx, aes(colour = Threshold)) +
#   facet_wrap(~Result) + basemap
#   
#   # scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1", 3, type = "discrete")) +
#   # scale_color_gradientn(colours = rev(inferno(4))) +
#   # scale_color_viridis(name = "", discrete=TRUE)+
#   # facet_wrap(~Scenario) +
#   # theme(legend.title=element_blank()) +
#   # geom_sf(data)
#   theme_bw()+
#   # scale_fill_viridis_d("Probability Difference")+
#   ggtitle("")
# 
# m1

file.name1 <- "Figures/03_map_4_cats_per_impact.jpg"
ggsave(m1, filename=file.name1, dpi=500, height=10, width=15)


## map

basemapsList <- c("Esri.WorldTopoMap", "Esri.WorldImagery",
                  "Esri.NatGeoWorldMap",
                  "OpenTopoMap", "OpenStreetMap", 
                  "CartoDB.Positron", "Stamen.TopOSMFeatures")

mapviewOptions(basemaps=basemapsList, fgb = FALSE)

## filter dry season and make spatial
impsx <- impx1 %>%
  filter(Hydro_endpoint == "DS_Mag_50", Index == "csci",
         Result == "BothImpact", !Threshold %in% c("NATLow", "NATMed", "NATHigh")) %>%
  st_as_sf(coords=c("longitude", "latitude"), crs=4326, remove=F)
names(impsx)
unique(impsx$masterid)

class(impsx)
st_write(impsx, "output_data/03_bothImpact_dry_season_csci.shp")
?st_write
m1 <- mapview(impsx, zcol = "Threshold",  col.regions=c("red", "green", "orange", "blue"),
              layer.name="Channel Type")


m1@map %>% leaflet::addMeasure(primaryLengthUnit = "meters")
# 
mapshot(m1, url = paste0(getwd(), "/output_data/01_dry_season_bfl_channel_types_flow_problem.html"),
        file = paste0(getwd(), "/ignore/01_dry_season_bfl_channel_types_flow_problem.png"))





