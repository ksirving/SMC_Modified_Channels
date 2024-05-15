## level 2: site specific analysis,  how many of mod chans have issues with flow metric, keep 4 outcomes 

library(gam)
library(tidyverse)
library(sf)
library(tidylog)

## directory for figures
out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/SMC_Modified_Channels/Figures/"

## workflow
## use 4 outcome categories to count how many mod channels have issues with flow metrics
## which flow metrics - visualise
## which flow metrics also issue with bio - visualise


# Upload Data -------------------------------------------------------------

## within limits doc
imps <- read.csv("output_data/05_impact_ffm_bio.csv") %>%
  select( -X)

## tally of categories
tal <- read.csv("output_data/05_count_impact.csv")
head(tal)
## get percentages

tal <- tal %>%
  group_by(Index, Hydro_endpoint, Flow.Metric.Name, Threshold, Quantile, metric) %>%
  mutate(TotalPerCat = sum(n)) %>%
  mutate(PercChans = (n/TotalPerCat)*100)

# Clean Tables ------------------------------------------------------------


## change order of cols, rename cols, change names of modified streams calss
FinalTablex <- tal %>% ## using option with least NAs
  ungroup() %>%
  select(Index, Flow.Metric.Name,  Threshold, Result,  n, PercChans, Quantile) %>% # 
  rename(FlowMetric = Flow.Metric.Name, ModifiedClass = Threshold, Impact = Result, PercentageOfSites = PercChans, NumberOfSites = n,) %>% #  
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
  group_by(Index, FlowMetric, ModifiedClass, Quantile) %>%
  summarise(NumberSitesPerClass = sum(NumberOfSites))

sums

write.csv(sums, "output_data/06_number_of_sites_per_FFM.csv")

  
## pivot impact results wider
FinalTable <- FinalTablex %>%  
  select(-NumberOfSites) %>%
  # filter(SmoothingFunc == 6) %>% ## change later!
  pivot_wider(names_from = Impact, values_from = PercentageOfSites) %>%
  mutate(across(everything(), .fns = ~replace_na(.,0))) %>%
  select(Index:Quantile, NoImpact, BioImpact, HydroImpact, BothImpact)

## join number of sites

FinalTable <- full_join(FinalTable, sums, by = c("Index", "FlowMetric", "ModifiedClass", "Quantile"))
FinalTable

## save out

write.csv(FinalTable, "output_data/06_percent_impacts_each_Class.csv")

## table with columns as channel class

FinalTableWide <- FinalTablex %>%
  select(-NumberOfSites) %>%
  pivot_wider(names_from = ModifiedClass, values_from = PercentageOfSites) %>%
  mutate(across(everything(), .fns = ~replace_na(.,0))) #%>%
  # select(-)

  FinalTableWide
  
  write.csv(FinalTableWide, "output_data/06_percent_impacts_each_Class_wide.csv")

# Map results -------------------------------------------------------------

## upload boundary

socal <- st_read("/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/SpatialData/SoCal_Counties/SoCal_Counties.shp")
socal

plot(socal)

## for now remove natlow,med,high
removes <- unique(imps$Threshold)[c(1,5,7)]

## make spatial

imps_sf <- imps %>% 
  # st_as_sf(coords=c( "longitude", "latitude"), crs=4326, remove=F) %>%
  # st_transform(crs = st_crs(socal)) %>%
  select(Index, Hydro_endpoint, Threshold, BioThresh, SmoothingFunc, masterid, COMID, Flow.Metric.Name, Flow.Component, Result, longitude, latitude)  %>%
  mutate(Threshold = as.factor(Threshold), Result = as.factor(Result)) %>%
  filter(!Threshold %in% removes)


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

# save(basemap, file = "Figures/06_basemap_socal_ggplot.RData")

load(file = "Figures/03_basemap_socal_ggplot.RData")

# Plot per flow metric and result

metrics <- unique(imps_sf$Hydro_endpoint)
metrics

## define index

ind <- unique(imps_sf$Index)
ind

## loop over metrics and index
i=2
m=2
imps_sf
  
  for(m in 1:length(metrics)) {
    
    imps_sf1 <- imps_sf %>%
      filter(Hydro_endpoint == metrics[m])
    
    for(i in 1:length(ind)) {
      
      ## extract - index
      imps_sfx <- imps_sf1 %>%
        filter(Index == ind[i])
    
    m1 <- basemap + 
      geom_point(data = imps_sfx, aes(x = longitude, y = latitude, colour = Threshold), size = 5) +
      scale_color_jco(name = "Modified Class") +
      guides(size = "none") +
      facet_grid(rows = vars(Result), cols = vars(Quantile))
    m1
    
    file.name1 <- paste0("Figures/06_", ind[i], "_", metrics[m], "_map_4_cats_per_impact.jpg")

    ggsave(m1, filename=file.name1, dpi=500, height=10, width=15)
    
  }
  
}

# Plot per flow metric and channel type

## define impact acategory
type <- unique(imps_sf$Threshold)

for(m in 1:length(metrics)) {
  
  imps_sf1 <- imps_sf %>%
    filter(Hydro_endpoint == metrics[m])
  
  for(i in 1:length(ind)) {
    
    ## extract - index
    imps_sfx <- imps_sf1 %>%
      filter(Index == ind[i])
    
    ## map 
    m1 <- basemap + 
      geom_point(data = imps_sfx, aes(x = longitude, y = latitude, colour = Result), size = 5) +
      scale_color_jco(name = "Impact Type") + ## colours re-order
      guides(size = "none") +
      # facet_wrap(~Threshold)
      facet_grid(rows = vars(Threshold), cols = vars(Quantile))
    
    m1
    file.name1 <- paste0("Figures/06_", ind[i], "_", metrics[m], "_map_4_cats_per_channel_type.jpg")
    file.name1
    ggsave(m1, filename=file.name1, dpi=500, height=10, width=15)
    
  }
  
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

file.name1 <- "Figures/06_map_4_cats_per_impact.jpg"
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
st_write(impsx, "output_data/06_bothImpact_dry_season_csci.shp")
?st_write
m1 <- mapview(impsx, zcol = "Threshold",  col.regions=c("red", "green", "orange", "blue"),
              layer.name="Channel Type")


m1@map %>% leaflet::addMeasure(primaryLengthUnit = "meters")
# 
mapshot(m1, url = paste0(getwd(), "/output_data/01_dry_season_bfl_channel_types_flow_problem.html"),
        file = paste0(getwd(), "/ignore/01_dry_season_bfl_channel_types_flow_problem.png"))





