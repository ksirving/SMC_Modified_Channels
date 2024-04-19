## exploring possible appraoches to flow eco in modified channels

library(tidylog)
library(tidyverse)
library(sf)
library(mapview)

out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/SMC_Modified_Channels/Figures/"

getwd()


# Flow data ---------------------------------------------------------------

dh_data <- read.csv("/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/SMC_Modified_Channels/ignore/2023-07-20_RFpred_output_alldata_biositesCOMIDs_med_dlt_FFM_test12_test2.csv")
head(dh_data)

dh_median <- dh_data %>%
  pivot_longer(d_ds_mag_50:delta_q99, names_to = "flow_metric", values_to = "deltah_final") 
  # select(site, flow_metric, deltah_final)
dim(dh_median) # 7236
head(dh_median)

## full names for labels
labels <- read.csv("input_data/ffm_names.csv")
labels <- labels[1:24, ]
labels
labels <- labels %>% rename(hydro.endpoints = Flow.Metric.Code)
labels[25, 1] <- "Magnitude of largest annual storm"
labels[25, 2] <- "Q99"
labels[25, 3] <- "Peak flow"
labels


labels <- labels %>%
  mutate(flow_metric = case_when(hydro.endpoints == "DS_Mag_50" ~ "d_ds_mag_50",
                         hydro.endpoints == "FA_Mag" ~ "d_fa_mag",
                         hydro.endpoints == "Peak_10" ~ "d_peak_10",
                         hydro.endpoints == "Peak_2" ~ "d_peak_2",
                         hydro.endpoints == "Peak_5" ~ "d_peak_5",
                         hydro.endpoints == "SP_Mag" ~ "d_sp_mag",
                         hydro.endpoints == "Wet_BFL_Mag_10" ~ "d_wet_bfl_mag_10",
                         hydro.endpoints == "Wet_BFL_Mag_50" ~ "d_wet_bfl_mag_50",
                         hydro.endpoints == "Q99" ~ "delta_q99")) %>%
  drop_na(flow_metric)

## join with data

## sites with flow data

flowSites <- unique(dh_median$masterid)

flowSites

## channel engineering data

BioEng <- read.csv("ignore/02_chan_eng.csv") %>%
  select(-c(X,channel_engineering_classification_date, channel_engineering_personnel, channel_engineering_comments)) %>%
  mutate(Class2 = ifelse(channel_engineering_class =="NAT", "Natural", "Modified"))

BioEng

tallyFFMClass <- BioEng %>%
  group_by(channel_engineering_class) %>%
  # select(masterid, COMID) %>%
  distinct() %>%
  drop_na(channel_engineering_class) %>%
  tally()

tallyFFMClass

## count number of asci/csci sites per channel class
## boxplots of delta h for all channel types
## test mod streams with derived delta h limits

# Bio data ----------------------------------------------------------------

##  sites only

bioSites <- st_read("input_data/01_bio_sites_all.shp")
head(bioSites)
dim(bioSites)

# set background basemaps:
basemapsList <- c("Esri.WorldTopoMap", "Esri.WorldImagery",
                  "Esri.NatGeoWorldMap",
                  "OpenTopoMap", "OpenStreetMap", 
                  "CartoDB.Positron", "Stamen.TopOSMFeatures")

mapviewOptions(basemaps=basemapsList, fgb = FALSE)

# this map of all sites 
m1 <- mapview(bioSites, cex=2, col.regions="black",
              layer.name="Bio Sites") 



m1
m1@map %>% leaflet::addMeasure(primaryLengthUnit = "meters")


## bugs data
load(file = "/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/Cannabis_Eflows/ignore/SMC_csci_cali_Mar2023_v2.RData")
# csciScores <- read.csv("ignore/01_csci_comp_mets_comid_socal.csv")
head(bug_tax_ca)

### algae scores
load(file="/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/Cannabis_Eflows/ignore/SMC_asci_cali_Mar2023_v3.RData")
head(alg_tax_ca)
# asciScores <- read.csv("ignore/01_asci_comp_mets_comid_socal.csv")

# Join bio and flow -------------------------------------------------------

## how many flow in bug sites
sum(flowSites %in% bug_tax_ca$masterid) ## 472

sum(flowSites %in% alg_tax_ca$masterid) ## 354

# ## filter bug sites to flow sites
# bugflowSites <- bioSites %>%
#   filter(masterid %in% flowSites)

## filter bug data using masterid ### remove reps - remove 2nd rep for now, change later!!!!
csciScores <- bug_tax_ca %>%
  # select(-X, -stationcode) %>%
  filter(fieldreplicate == 1 ) %>%
  # separate(sampledate, into = c("sampledate", "Time"), sep= "T", remove = F) %>%
  # separate(sampledate, into = c("year", "Month", "Day"), sep= "-", remove = F) %>%
  # mutate(year = as.numeric(year)) %>%
  mutate(Metric = "csci", csci = as.numeric(csci)) %>%
  rename(MetricValue = csci, COMID = comid) %>%
  select(masterid, sampleyear, Metric, MetricValue, longitude, latitude, COMID)

length(unique(csciScores$masterid)) ## 4413 sites (all ca)


## filter bug data using masterid ### remove reps - remove 2nd rep for now, change later!!!!
asciScores <- alg_tax_ca %>%
  # select(-X, -stationcode) %>%
  filter(replicate == 1 ) %>%
  separate(sampledate, into = c("sampledate", "Time"), sep= " ", remove = F) %>%
  separate(sampledate, into = c("sampleyear", "Month", "Day"), sep= "-", remove = F) %>%
  mutate(sampleyear = as.numeric(sampleyear))  %>%
  mutate(Metric = "asci") %>%
  rename(MetricValue = result, COMID = comid) %>%
  mutate(MetricValue = as.numeric(MetricValue)) %>%
  select(masterid, sampleyear, Metric, MetricValue, longitude, latitude,  COMID)


str(asciScores)
length(unique(asciScores$masterid)) ## 2306 sites 


names(asciScores)
names(csciScores)

# Join bio sites to flow data ---------------------------------------------

## join asci and csci
scoresSites <- bind_rows(asciScores, csciScores)

## join channel class to bio sites
engsites <- inner_join(scoresSites, BioEng, by = "masterid")
engsites

length(unique(engsites$masterid))

## tally asci csci by class 

tallyclass <- engsites %>%
  group_by(Metric, channel_engineering_class) %>%
  select(-sampleyear, - MetricValue) %>% distinct() %>%
  tally()
tallyclass
head(dh_median)
flowsites <- dh_median %>%
  select(masterid) %>%
  mutate(HasFFM = "Yes")

flowsites

## join to flow 
AllData1 <- full_join(engsites, flowsites, by = c("masterid"), relationship = "many-to-many") %>%
  mutate(HasFFM = replace_na(HasFFM, "No"))
AllData1 


## count FFM per channel class

tallyFFM <- AllData1 %>%
  group_by(channel_engineering_class, HasFFM) %>%
  select(-sampleyear, - MetricValue, -Metric) %>% distinct() %>%
  drop_na(channel_engineering_class) %>%
  tally()

tallyFFM

## join all data
 
AllData <- right_join(engsites, dh_median, by = c("masterid", "longitude", "latitude"), relationship = "many-to-many")%>%
     left_join(labels, by = "flow_metric")
AllData

## save out
save(AllData, file = "output_data/00_bugs_algae_flow_joined_by_masterid.RData")

# AllDataA <- left_join(asciScoresLA,  dh_median, by = c("masterid")) 
# head(AllDataA)
# dim(AllDataA) ## 
# 
# ## save out
# save(AllDataA, file = "output_data/00_algae_flow_joined_by_masterid.RData")


## add points to figures
## scam package
## add all soft bottom together, hard, soft

# Ranges of Delta H -------------------------------------------------------

load(file = "output_data/00_bugs_algae_flow_joined_by_masterid.RData")
head(AllData)


AllDatax <- AllData %>%
  drop_na(channel_engineering_class) %>%
  mutate(channel_engineering_class = recode_factor(channel_engineering_class, NAT = "Natural", SB0 = "Soft: Zero Hard Sides",
                                                   SB1 = "Soft: One Hard Side", SB2 = "Soft: Two Hard Sides", 
                                                   HB = "Hard Bottom"))

unique(AllDatax$channel_engineering_class)
## box plot of Ranges
m=1
mets <- unique(na.omit(AllDatax$Flow.Metric.Name))
mets

for(m in 1:length(mets)) {
  
  T1 <- (ggplot(subset(AllDatax, Flow.Metric.Name == mets[m]),  aes(x=channel_engineering_class, y=deltah_final)) +
           geom_boxplot() +
           scale_x_discrete(name=paste("")) +
           scale_y_continuous(name = paste0("Delta: ", mets[m])))

  
  T1
  
  file.name1 <- paste0(out.dir, "00_", mets[m], "_boxplot_delta_range.jpg")
  ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
}

## jitter plot: height to be 0, unchanges, changed below a therehold, changed above a threshold. delta h limits below critical 


# GLM Models:  combined channel types------------------------------------------------------------------

AllDataLong2 <- AllData %>%
  filter(Metric %in% c("csci", "ASCI_Hybrid"))

## tally of pos/negs per metric and channel type

names(AllDataLong2)

tally_3 <-
  AllDataLong2 %>% drop_na(channel_engineering_class) %>%
  # filter(Metric == "csci") %>%
  mutate(Type = ifelse(deltah_final < 0, "Negative", "Positive")) %>%
  group_by(flow_metric, channel_engineering_class, Metric, Type) %>% 
  tally() %>% 
  arrange(flow_metric, channel_engineering_class, Metric, Type) %>%
  pivot_wider(names_from = Type, values_from = n) %>%
  mutate(Total = Negative + Positive)

print.data.frame(tally_3)

## save
write.csv(tally_3, "output_data/00_sites_per_category_flow_metric.csv")

## tally per category

tally_2 <-
  AllDataLong2 %>% drop_na(channel_engineering_class) %>%
  # filter(Metric == "csci") %>%
  mutate(Type = ifelse(deltah_final < 0, "Negative", "Positive")) %>%
  group_by(channel_engineering_class, Metric, Type) %>% 
  tally() %>% 
  arrange(channel_engineering_class, Metric, Type) %>%
  pivot_wider(names_from = Type, values_from = n) %>%
  mutate(Total = Negative + Positive)

print.data.frame(tally_2)

## save
write.csv(tally_2, "output_data/00_sites_per_category.csv")


unique(AllDataLong2$flow_metric)
## create df of model configurations

## bio 
biol.endpoints<-unique(AllDataLong2$Metric)
biol.endpoints

## flow
flow.endpoints<- unique(na.omit(AllDataLong2$flow_metric))
flow.endpoints

# direction of alteration

Direction.alt <- c("Pos", "Neg")
Direction.alt

# channel type

channel.type.combined <- unique(na.omit(AllDataLong2$Class2))
channel.type.sep <- unique(na.omit(AllDataLong2$channel_engineering_class))

### make grid
bio_h_summary<-  expand.grid(biol.endpoints=biol.endpoints,flow.endpoints=flow.endpoints, 
                             Direction.alt= Direction.alt, channel.type = channel.type.combined, stringsAsFactors = F)

## add bio thresholds: standard
bio_h_summary <- bio_h_summary %>%
  mutate(BioThresholds = ifelse(biol.endpoints == "csci", 0.79, 0.86))

bio_h_summary
# bio_h_summary <- bio_h_summary[-c(1:9),]
bio_h_summary
head(AllDataLong2)
i=1
## model of each configuration
log.lm <-lapply(1:nrow(bio_h_summary), function(i)
{
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  dmet<-as.character(bio_h_summary[i,"Direction.alt"])
  cmet<-as.character(bio_h_summary[i,"channel.type"])
  thresh<-(bio_h_summary[i,"BioThresholds"])
  
  if(dmet == "Pos") {
    
    mydat<-AllDataLong2 %>%
      filter(Metric == bmet,
             flow_metric == tmet, Class2 == cmet) %>%
      select(Metric, MetricValue, deltah_final, Class2) %>% ## only metrics needed
      drop_na( deltah_final) %>%
      filter(deltah_final > 0) %>%
      filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
      distinct()
  } else {
    
    mydat<-AllDataLong2 %>%
      filter(Metric == bmet,
             flow_metric == tmet,  Class2 == cmet) %>%
      select(Metric, MetricValue, deltah_final, Class2) %>% ## only metrics needed
      drop_na( deltah_final) %>%
      filter(deltah_final < 0) %>%
      filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
      distinct()
  }
  


  mydat$Condition<-ifelse(mydat$MetricValue< thresh ,0, 1) ## convert to binary
  # mydat <- mydat %>% drop_na(Condition)
  
  mydat<-mydat[order(mydat$MetricValue),] ## order by csci value
  try(glm(Condition~deltah_final, family=binomial(link="logit"), data=mydat), silent=T) ### glm
  
  
})
log.lm
## save models
save(log.lm, file = "output_data/00_csci_asci_all_glm_flow_combined.RData")

### get rsqds and pvals
for(i in 1:length(log.lm)) {
  if (class(log.lm[[i]]) == "try-error") {
    
    # mod <- summary(log.lm[[i]])
    bio_h_summary$AIC[i] <- NA ##1-mod$deviance/mod$null.deviance ## mcfaddens r2
    bio_h_summary$PValue[i] <- NA
    bio_h_summary$McFaddensR2[i] <- NA
    bio_h_summary$n[i] <- NA
    
  } else {
    
    mod <- summary(log.lm[[i]])
    bio_h_summary$AIC[i] <- mod$aic ##1-mod$deviance/mod$null.deviance ## mcfaddens r2
    bio_h_summary$PValue[i] <- mod$coefficients[8]
    bio_h_summary$McFaddensR2[i] <- 1-mod$deviance/mod$null.deviance
    bio_h_summary$n[i] <- mod$df[2]+1
  }
  
}
## save configs and r sqds
save(bio_h_summary, file="output_data/00_csci_asci_glm_rsqds_combined.RData")
bio_h_summary
## make df of predicted values to predict on - need to be different for each flow metric

## blank df
DF <- NULL
DF <- as.data.frame(DF)
AllDataLong2
i=1
### get predictions and fitted values
for(i in 1:length(log.lm)) {
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  dmet<-as.character(bio_h_summary[i,"Direction.alt"])
  cmet<-as.character(bio_h_summary[i,"channel.type"])
  tmet
  
  if (class(log.lm[[i]]) == "try-error") {
    ## get NAs where no data/model
    DFX <- cbind(bmet, as.numeric(NA), NA, cmet, NA, NA, dmet, tmet) %>% as.data.frame() %>%
      rename(Metric = bmet, MetricValue = V2, Value = V3, ChannelType = cmet, Condition = V5, 
             predictedVals = V6, DirectionAlt = dmet, Variable = tmet) 
    
  } else {
    
    ## get model, predict, extract all data and categories
    mod <- log.lm[[i]]
    predictedVals <- predict.glm(mod,  type = "response")

    DFX <- cbind(na.omit(mod$data), as.data.frame(predictedVals)) %>%
      rename(ChannelType = Class2, Value = deltah_final) %>%
      mutate(DirectionAlt = dmet, Variable = tmet) %>%
      mutate(MetricValue = as.character(MetricValue), Value = as.character(Value),
             Condition = as.character(Condition), predictedVals = as.character(predictedVals))
    
  }
  
  DF <- bind_rows(DF, DFX)
  
}

## change back to numeric
DF <- DF %>%
  mutate(Value = as.numeric(Value), MetricValue = as.numeric(MetricValue),
         predictedVals = as.numeric(predictedVals), Condition = as.numeric(Condition))
str(DF)

### predicted figures
bio <- unique(DF$Metric)
bio
mets <- unique(DF$Variable)
DF

for(b in 1:length(bio)) {
  
  DF1 <- DF %>%
    filter(Metric == bio[b])
  
  for(m in 1:length(mets)) {
    
    T1 <- (ggplot(subset(DF1, Variable == mets[m]), aes(y=predictedVals, x=Value, group = ChannelType, color = ChannelType)) +
             # geom_point(size=0.2) +
             geom_line( linewidth = 1)+
             # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
             # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
             facet_wrap(~DirectionAlt, scales = "free") +
             scale_x_continuous(name=paste(mets[m])) +
             scale_y_continuous(name = paste0("Probability of Good ", bio[b]))) 
    # theme(legend.position = "none"))
    
    T1
    
    file.name1 <- paste0(out.dir, "00_", bio[b], "_", mets[m], "_flow_response_predicted_glm_combined.jpg")
    ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
  
}


}



# GLM Models: separate channel types------------------------------------------------------------------

AllDataLong2 <- AllData %>%
  filter(Metric %in% c("csci", "ASCI_Hybrid"))

unique(AllDataLong2$flow_metric)



## create df of model configurations

## bio 
biol.endpoints<- unique(AllDataLong2$Metric)
biol.endpoints


## flow
flow.endpoints<- unique(na.omit(AllDataLong2$flow_metric))
# flow.endpoints <- flow.endpoints[-c(14:16)]

# direction of alteration

Direction.alt <- c("Pos", "Neg")
Direction.alt

# channel type

channel.type.combined <- unique(na.omit(AllDataLong2$Class2))
channel.type.sep <- unique(na.omit(AllDataLong2$channel_engineering_class))


bio_h_summary<-  expand.grid(biol.endpoints=biol.endpoints,flow.endpoints=flow.endpoints, 
                             Direction.alt= Direction.alt, channel.type = channel.type.sep, stringsAsFactors = F)

## add bio thresholds: standard
bio_h_summary <- bio_h_summary %>%
  mutate(BioThresholds = ifelse(biol.endpoints == "csci", 0.79, 0.86))

# bio_h_summary <- bio_h_summary[-c(1:9),]
bio_h_summary
head(AllDataLong2)

## model of each configuration
log.lm <-lapply(1:nrow(bio_h_summary), function(i)
{
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  dmet<-as.character(bio_h_summary[i,"Direction.alt"])
  cmet<-as.character(bio_h_summary[i,"channel.type"])
  thresh<-(bio_h_summary[i,"BioThresholds"])

  
  if(dmet == "Pos") {
    
    mydat<-AllDataLong2 %>%
      filter(Metric == bmet,
             flow_metric == tmet, channel_engineering_class == cmet) %>%
      select(Metric, MetricValue, deltah_final, channel_engineering_class) %>% ## only metrics needed
      drop_na( deltah_final) %>%
      filter(deltah_final > 0) %>%
      filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
      distinct()
    
  } else {
    
    mydat<-AllDataLong2 %>%
      filter(Metric == bmet,
             flow_metric == tmet,  channel_engineering_class == cmet) %>%
      select(Metric, MetricValue, deltah_final, channel_engineering_class) %>% ## only metrics needed
      drop_na( deltah_final) %>%
      filter(deltah_final < 0) %>%
      filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
      distinct()
  }
  
  

  mydat$Condition<-ifelse(mydat$MetricValue< thresh ,0, 1) ## convert to binary
  # mydat <- mydat %>% drop_na(Condition)
  
  mydat<-mydat[order(mydat$MetricValue),] ## order by csci value
  try(glm(Condition~deltah_final, family=binomial(link="logit"), data=mydat), silent=T) ### glm

  
})

log.lm
## save models
save(log.lm, file = "output_data/00_csci_asci_all_glm_flow_separate.RData")

### get rsqds and pvals
for(i in 1:length(log.lm)) {
  
  if (class(log.lm[[i]]) == "try-error") {
    
    # mod <- summary(log.lm[[i]])
    bio_h_summary$AIC[i] <- NA ##1-mod$deviance/mod$null.deviance ## mcfaddens r2
    bio_h_summary$PValue[i] <- NA
    bio_h_summary$McFaddensR2[i] <- NA
    bio_h_summary$n[i] <- NA
    
  } else {
    
    mod <- summary(log.lm[[i]])
    bio_h_summary$AIC[i] <- mod$aic ##1-mod$deviance/mod$null.deviance ## mcfaddens r2
    bio_h_summary$PValue[i] <- mod$coefficients[8]
    bio_h_summary$McFaddensR2[i] <- 1-mod$deviance/mod$null.deviance
    bio_h_summary$n[i] <- mod$df[2]+1
  }
  
}
## save configs and r sqds
save(bio_h_summary, file="output_data/00_csci_asci_glm_rsqds_separate.RData")
bio_h_summary
## make df of predicted values to predict on - need to be different for each flow metric

## blank df
DF <- NULL
DF <- as.data.frame(DF)

### get predictions and fitted values
for(i in 1:length(log.lm)) {

  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  dmet<-as.character(bio_h_summary[i,"Direction.alt"])
  cmet<-as.character(bio_h_summary[i,"channel.type"])

  if (class(log.lm[[i]]) == "try-error") {
    ## get NAs where no data/model
    DFX <- cbind(bmet, as.numeric(NA), NA, cmet, NA, NA, dmet, tmet) %>% as.data.frame() %>%
      rename(Metric = bmet, MetricValue = V2, Value = V3, ChannelType = cmet, Condition = V5, 
             predictedVals = V6, DirectionAlt = dmet, Variable = tmet) 
    
  } else {
    
    ## get model, predict, extract all data and categories
    mod <- log.lm[[i]]
    predictedVals <- predict.glm(mod,  type = "response")
    
    DFX <- cbind(na.omit(mod$data), as.data.frame(predictedVals)) %>%
      rename(ChannelType = channel_engineering_class, Value = deltah_final) %>%
      mutate(DirectionAlt = dmet, Variable = tmet) %>%
      mutate(MetricValue = as.character(MetricValue), Value = as.character(Value),
             Condition = as.character(Condition), predictedVals = as.character(predictedVals))
    
  }
  
  DF <- bind_rows(DF, DFX)
  
}
dmet
## change back to numeric
DF <- DF %>%
  mutate(Value = as.numeric(Value), MetricValue = as.numeric(MetricValue),
         predictedVals = as.numeric(predictedVals), Condition = as.numeric(Condition))
str(DF)

### predicted figures
bio <- unique(DF$Metric)
bio
mets <- unique(DF$Variable)
DF

b=1

for(b in 1:length(bio)) {
  
  DF1 <- DF %>%
    filter(Metric == bio[b])
  


for(m in 1:length(mets)) {
  
  T1 <- (ggplot(subset(DF1, Variable == mets[m]), aes(y=predictedVals, x=Value, group = ChannelType, color = ChannelType)) +
           # geom_point(size=0.2) +
           geom_line( linewidth = 1)+
           # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
           # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
           facet_wrap(~DirectionAlt, scales = "free") +
           scale_x_continuous(name=paste(mets[m])) +
           scale_y_continuous(name = paste0("Probability of Good ", bio[b]))) 
  # theme(legend.position = "none"))
  

  
  file.name1 <- paste0(out.dir, "00_", bio[b], "_", mets[m], "_flow_response_predicted_glm_separate.jpg")
  ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
}
}

# GLM Models: separate channel types & thresholds------------------------------------------------------------------

AllDataLong2 <- AllData %>%
  filter(Metric %in% c("csci", "ASCI_Hybrid"))

unique(AllDataLong2$flow_metric)



## create df of model configurations

## bio 
biol.endpoints<- unique(AllDataLong2$Metric)
biol.endpoints


## flow
flow.endpoints<- unique(na.omit(AllDataLong2$flow_metric))
# flow.endpoints <- flow.endpoints[-c(14:16)]

# direction of alteration

Direction.alt <- c("Pos", "Neg")
Direction.alt

# channel type

channel.type.combined <- unique(na.omit(AllDataLong2$Class2))
channel.type.sep <- unique(na.omit(AllDataLong2$channel_engineering_class))


bio_h_summary<-  expand.grid(biol.endpoints=biol.endpoints,flow.endpoints=flow.endpoints, 
                             Direction.alt= Direction.alt, channel.type = channel.type.sep, stringsAsFactors = F)

## add bio thresholds: standard
bio_h_summary <- bio_h_summary %>%
  mutate(BioThresholds = case_when(biol.endpoints == "csci" & channel.type == "SB0" ~ 0.78,
                                   biol.endpoints == "csci" & channel.type == "SB1" ~ 0.78,
                                   biol.endpoints == "csci" & channel.type == "SB2" ~ 0.75,
                                   biol.endpoints == "csci" & channel.type == "HB" ~ 0.67,
                                   biol.endpoints == "csci" & channel.type == "NAT" ~ 0.79,
                                   biol.endpoints == "ASCI_Hybrid" & channel.type == "SB0" ~ 0.79,
                                   biol.endpoints == "ASCI_Hybrid" & channel.type == "SB1" ~ 0.79,
                                   biol.endpoints == "ASCI_Hybrid" & channel.type == "SB2" ~ 0.76,
                                   biol.endpoints == "ASCI_Hybrid" & channel.type == "HB" ~ 0.87,
                                   biol.endpoints == "ASCI_Hybrid" & channel.type == "NAT" ~ 0.86))

# bio_h_summary <- bio_h_summary[-c(1:9),]
bio_h_summary
head(AllDataLong2)

## model of each configuration
log.lm <-lapply(1:nrow(bio_h_summary), function(i)
{
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  dmet<-as.character(bio_h_summary[i,"Direction.alt"])
  cmet<-as.character(bio_h_summary[i,"channel.type"])
  thresh<-(bio_h_summary[i,"BioThresholds"])
  
  
  if(dmet == "Pos") {
    
    mydat<-AllDataLong2 %>%
      filter(Metric == bmet,
             flow_metric == tmet, channel_engineering_class == cmet) %>%
      select(Metric, MetricValue, deltah_final, channel_engineering_class) %>% ## only metrics needed
      drop_na( deltah_final) %>%
      filter(deltah_final > 0) %>%
      filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
      distinct()
    
  } else {
    
    mydat<-AllDataLong2 %>%
      filter(Metric == bmet,
             flow_metric == tmet,  channel_engineering_class == cmet) %>%
      select(Metric, MetricValue, deltah_final, channel_engineering_class) %>% ## only metrics needed
      drop_na( deltah_final) %>%
      filter(deltah_final < 0) %>%
      filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
      distinct()
  }
  
  
  
  mydat$Condition<-ifelse(mydat$MetricValue< thresh ,0, 1) ## convert to binary
  # mydat <- mydat %>% drop_na(Condition)
  
  mydat<-mydat[order(mydat$MetricValue),] ## order by csci value
  try(glm(Condition~deltah_final, family=binomial(link="logit"), data=mydat), silent=T) ### glm
  
  
})

log.lm
## save models
save(log.lm, file = "output_data/00_csci_asci_all_glm_flow_separate_mod_thresh.RData")

### get rsqds and pvals
for(i in 1:length(log.lm)) {
  
  if (class(log.lm[[i]]) == "try-error") {
    
    # mod <- summary(log.lm[[i]])
    bio_h_summary$AIC[i] <- NA ##1-mod$deviance/mod$null.deviance ## mcfaddens r2
    bio_h_summary$PValue[i] <- NA
    bio_h_summary$McFaddensR2[i] <- NA
    bio_h_summary$n[i] <- NA
    
  } else {
    
    mod <- summary(log.lm[[i]])
    bio_h_summary$AIC[i] <- mod$aic ##1-mod$deviance/mod$null.deviance ## mcfaddens r2
    bio_h_summary$PValue[i] <- mod$coefficients[8]
    bio_h_summary$McFaddensR2[i] <- 1-mod$deviance/mod$null.deviance
    bio_h_summary$n[i] <- mod$df[2]+1
  }
  
}
## save configs and r sqds
save(bio_h_summary, file="output_data/00_csci_asci_glm_rsqds_separate_mod_thresh.RData")
bio_h_summary
## make df of predicted values to predict on - need to be different for each flow metric

## blank df
DF <- NULL
DF <- as.data.frame(DF)

### get predictions and fitted values
for(i in 1:length(log.lm)) {
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  dmet<-as.character(bio_h_summary[i,"Direction.alt"])
  cmet<-as.character(bio_h_summary[i,"channel.type"])
  
  if (class(log.lm[[i]]) == "try-error") {
    ## get NAs where no data/model
    DFX <- cbind(bmet, as.numeric(NA), NA, cmet, NA, NA, dmet, tmet) %>% as.data.frame() %>%
      rename(Metric = bmet, MetricValue = V2, Value = V3, ChannelType = cmet, Condition = V5, 
             predictedVals = V6, DirectionAlt = dmet, Variable = tmet) 
    
  } else {
    
    ## get model, predict, extract all data and categories
    mod <- log.lm[[i]]
    predictedVals <- predict.glm(mod,  type = "response")
    
    DFX <- cbind(na.omit(mod$data), as.data.frame(predictedVals)) %>%
      rename(ChannelType = channel_engineering_class, Value = deltah_final) %>%
      mutate(DirectionAlt = dmet, Variable = tmet) %>%
      mutate(MetricValue = as.character(MetricValue), Value = as.character(Value),
             Condition = as.character(Condition), predictedVals = as.character(predictedVals))
    
  }
  
  DF <- bind_rows(DF, DFX)
  
}
dmet
## change back to numeric
DF <- DF %>%
  mutate(Value = as.numeric(Value), MetricValue = as.numeric(MetricValue),
         predictedVals = as.numeric(predictedVals), Condition = as.numeric(Condition))
str(DF)

### predicted figures
bio <- unique(DF$Metric)
bio
mets <- unique(DF$Variable)
DF

b=1

for(b in 1:length(bio)) {
  
  DF1 <- DF %>%
    filter(Metric == bio[b])
  
  
  
  for(m in 1:length(mets)) {
    
    T1 <- (ggplot(subset(DF1, Variable == mets[m]), aes(y=predictedVals, x=Value, group = ChannelType, color = ChannelType)) +
             # geom_point(size=0.2) +
             geom_line( linewidth = 1)+
             # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
             # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
             facet_wrap(~DirectionAlt, scales = "free") +
             scale_x_continuous(name=paste(mets[m])) +
             scale_y_continuous(name = paste0("Probability of Good ", bio[b]))) 
    # theme(legend.position = "none"))
    
    
    
    file.name1 <- paste0(out.dir, "00_", bio[b], "_", mets[m], "_flow_response_predicted_glm_separate_mod_thresh.jpg")
    ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
  }
}


# Exceedance to delta limits ----------------------------------------------

## delta limits upload
limits <- read.csv("/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/SOC_Irving_et_al_data/output_data/Manuscript/07_ALL_delta_thresholds_scaled.csv")
limits

## filter out only standard thresholds
limits <- limits %>%
  filter(Bio_threshold %in% c(0.79, 0.86))

limits

## clean up df, make longer
thresholds <- limits %>%
  select(-X.1, -X, -n) %>%
  rename(Hydro_Metric = Hydro_endpoint) %>% 
  mutate(metric = paste0(Bio_endpoint, "_", Hydro_Metric, "_", Bio_threshold)) %>%
  pivot_longer(Threshold25:Threshold75, names_to = "Threshold") %>%
  rename(DeltaH = value) 

## make wider with type - pos/neg
thresholdsall <- thresholds %>%
  pivot_wider(names_from = Type, values_from = DeltaH)

thresholdsall

## delta H 
head(AllDataLong2)
delta <- AllDataLong2 %>%
  rename(Hydro_Metric = flow_metric) #%>%

## join limits with delta data
delta_df <- full_join(delta, thresholdsall, by="Hydro_Metric") %>%
  drop_na(Bio_threshold)
head(delta_df)

## define alteration per subbasin, per year - within limits
delta_dfx <- delta_df %>%
  group_by(masterid, Hydro_Metric, Biol, Threshold, Bio_threshold) %>%
  mutate(Alteration_Current = ifelse(deltah_final <= Positive & deltah_final >= Negative, "Unaltered", "Altered")) 

head(delta_dfx)

### percentage alteration per class

PercAlt <- delta_dfx %>%
  ungroup() %>%
  drop_na(channel_engineering_class) %>%
  filter(Threshold == "Threshold50") %>%
  group_by(Hydro_Metric, Biol, Bio_threshold, channel_engineering_class) %>%
  count(Alteration_Current) %>%
  mutate(Percentage = n/sum(n)*100) %>%
  select(-n) %>%
  pivot_wider(names_from = Alteration_Current, values_from = Percentage) #%>%
  
  PercAlt
  
  

##  direction per site
metric_tally_current_dir <- delta_dfx %>%
  drop_na(channel_engineering_class) %>%
  filter(Threshold == "Threshold50") %>%
  group_by(masterid, Hydro_Metric, Biol, Threshold, Bio_threshold, channel_engineering_class) %>%
  filter(Alteration_Current == "Altered") %>%
  mutate(Alteration_Current_direction = ifelse(deltah_final > Positive, "High", "Low")) %>%
  ungroup(masterid, Threshold) %>%
  select(masterid,  channel_engineering_class, Hydro_Metric, Biol, Alteration_Current, Alteration_Current_direction) %>%
  count(Alteration_Current_direction) %>%
  mutate(Percentage = n/sum(n)*100) %>%
  select(-n) %>%
  pivot_wider(names_from = Alteration_Current_direction, values_from = Percentage) #%>%


alldata <- full_join(PercAlt, metric_tally_current_dir, by = c("Hydro_Metric", "Biol", "Bio_threshold", "channel_engineering_class")) %>%
  mutate(channel_engineering_class = recode_factor(channel_engineering_class, NAT = "NAT", SB0 = "SB0",
                                                   SB1 = "SB1", SB2 = "SB2", 
                                                   HB = "HB")) %>%
  mutate(Direction = ifelse(High > 50, "High", "Low"))

head(alldata)

alldata
## save

write.csv(PercAlt, "output_data/00_percent_alteration_per_channel_type.csv")
write.csv(metric_tally_current_dir, "output_data/00_alteration_direction_per_channel_type.csv")
write.csv(alldata, "output_data/00_alteration_per_channel_type")

read.csv("output_data/00_percent_alteration_per_channel_type.csv")


## plot

p <- ggplot(alldata,  aes(fill = Direction, x=channel_engineering_class, y=Altered)) +
    geom_bar(position = "dodge", stat="identity") +
    facet_grid(vars(Biol), vars(Hydro_Metric)) +
    scale_x_discrete(name=paste("")) 
    # scale_y_continuous(name = paste0("Delta: ", mets[m]))

p
file.name1 <- paste0(out.dir, "00_alteration_per_channelType.jpg")
ggsave(p, filename=file.name1, dpi=300, height=5, width=7.5)

### stacked bars to illustrate high/low
## jitter plot with delta h limits


# GAMs: All combined channel types ----------------------------------------

# GAMs: combined channel types --------------------------------------------
# install.packages("gam")
library(gam)

AllDataLong2 <- AllData %>%
  filter(Metric %in% c("csci", "asci"))


unique(AllDataLong2$flow_metric)
## create df of model configurations

## bio 
biol.endpoints<-unique(AllDataLong2$Metric)
biol.endpoints

## flow
flow.endpoints<- unique(na.omit(AllDataLong2$flow_metric))
flow.endpoints

### make grid
bio_h_summary<-  expand.grid(biol.endpoints=biol.endpoints,flow.endpoints=flow.endpoints, stringsAsFactors = F)

bio_h_summary
# bio_h_summary <- bio_h_summary[-c(1:9),]
bio_h_summary
head(AllDataLong2)
i=1
## model of each configuration
gam.lm <-lapply(1:nrow(bio_h_summary), function(i)
{
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  # dmet<-as.character(bio_h_summary[i,"Direction.alt"])
  # cmet<-as.character(bio_h_summary[i,"channel.type"])
  # thresh<-(bio_h_summary[i,"BioThresholds"])
  
  
  mydat<-AllDataLong2 %>%
    filter(Metric == bmet,
           flow_metric == tmet) %>%
    select(Metric, MetricValue, deltah_final, Class2) %>% ## only metrics needed
    drop_na( deltah_final) %>%
    # filter(deltah_final > 0) %>%
    filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
    distinct()
  
  
  mydat<-mydat[order(mydat$MetricValue),] ## order by csci value
  # try(glm(Condition~deltah_final, family=binomial(link="logit"), data=mydat), silent=T) ### glm
  
  gam(MetricValue~s(deltah_final), family = Gamma(link = "log"), data = mydat)
  
  
})
gam.lm
## save models
save(gam.lm, file = "output_data/00_csci_asci_all_gam_flow_combined_all.RData")

### get rsqds and pvals

for(i in 1:length(gam.lm)) {
  
  # if (class(gam.lm[[i]]) == "try-error") {
  #   
  #   # mod <- summary(log.lm[[i]])
  #   bio_h_summary$AIC[i] <- NA ##1-mod$deviance/mod$null.deviance ## mcfaddens r2
  #   bio_h_summary$PValue[i] <- NA
  #   bio_h_summary$McFaddensR2[i] <- NA
  #   bio_h_summary$n[i] <- NA
  #   
  
  mod <- summary(gam.lm[[i]])
  
  bio_h_summary$AIC[i] <- mod$aic ##1-mod$deviance/mod$null.deviance ## mcfaddens r2
  bio_h_summary$PValue[i] <- mod$anova$`Pr(F)`[2]
  bio_h_summary$McFaddensR2[i] <- 1-mod$deviance/mod$null.deviance
  # bio_h_summary$n[i] <- mod$df[2]+1
  
  
}


## save configs and r sqds
save(bio_h_summary, file="output_data/00_csci_asci_gam_rsqds_combined_all.RData")
bio_h_summary
## make df of predicted values to predict on - need to be different for each flow metric

## blank df
DF <- NULL
DF <- as.data.frame(DF)
AllDataLong2
i=1
### get predictions and fitted values
for(i in 1:length(gam.lm)) {
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  # dmet<-as.character(bio_h_summary[i,"Direction.alt"])
  # cmet<-as.character(bio_h_summary[i,"channel.type"])
  
  
  
  ## get model, predict, extract all data and categories
  mod <- gam.lm[[i]]
  predictedVals <- predict(mod,  type = "response")
  predictedVals
  ?predict
  DFX <- cbind(na.omit(mod$data), as.data.frame(predictedVals)) %>%
    rename(Value = deltah_final) %>%
    mutate(Variable = tmet) 
  
  DF <- bind_rows(DF, DFX)
  
}
DF

## change back to numeric
DF <- DF %>%
  mutate(Value = as.numeric(Value), MetricValue = as.numeric(MetricValue),
         predictedVals = as.numeric(predictedVals))
str(DF)

### predicted figures
bio <- unique(DF$Metric)
bio
mets <- unique(DF$Variable)
DF
m=1
for(b in 1:length(bio)) {
  
  DF1 <- DF %>%
    filter(Metric == bio[b])
  
  for(m in 1:length(mets)) {
    
    T1 <- (ggplot(subset(DF1, Variable == mets[m]), aes(y=predictedVals, x=Value)) +
             # geom_point(size=0.2) +
             geom_line( linewidth = 1)+
             # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
             # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
             # facet_wrap(~DirectionAlt, scales = "free") +
             scale_x_continuous(name=paste(mets[m])) +
             scale_y_continuous(name = paste0("Probability of Good ", bio[b]))) 
    # theme(legend.position = "none"))
    
    T1
    
    file.name1 <- paste0(out.dir, "00_", bio[b], "_", mets[m], "_flow_response_predicted_gam_combined_all.jpg")
    ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
    
  }
  
  
}




# GAMs: separate channel types --------------------------------------------
# install.packages("gam")
library(gam)

AllDataLong2 <- AllData %>%
  filter(Metric %in% c("csci", "asci"))


unique(AllDataLong2$flow_metric)
## create df of model configurations

## bio 
biol.endpoints<-unique(AllDataLong2$Metric)
biol.endpoints

## flow
flow.endpoints<- unique(na.omit(AllDataLong2$flow_metric))
flow.endpoints


# channel type

channel.type.combined <- unique(na.omit(AllDataLong2$Class2))
channel.type.sep <- unique(na.omit(AllDataLong2$channel_engineering_class))

### make grid
bio_h_summary<-  expand.grid(biol.endpoints=biol.endpoints,flow.endpoints=flow.endpoints, 
                              channel.type = channel.type.combined, stringsAsFactors = F)

bio_h_summary
# bio_h_summary <- bio_h_summary[-c(1:9),]
bio_h_summary
head(AllDataLong2)
i=1
## model of each configuration
gam.lm <-lapply(1:nrow(bio_h_summary), function(i)
{
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  # dmet<-as.character(bio_h_summary[i,"Direction.alt"])
  cmet<-as.character(bio_h_summary[i,"channel.type"])
  # thresh<-(bio_h_summary[i,"BioThresholds"])
  
  
    mydat<-AllDataLong2 %>%
      filter(Metric == bmet,
             flow_metric == tmet, Class2 == cmet) %>%
      select(Metric, MetricValue, deltah_final, Class2) %>% ## only metrics needed
      drop_na( deltah_final) %>%
      # filter(deltah_final > 0) %>%
      filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
      distinct()
    

  mydat<-mydat[order(mydat$MetricValue),] ## order by csci value
  # try(glm(Condition~deltah_final, family=binomial(link="logit"), data=mydat), silent=T) ### glm

  gam(MetricValue~s(deltah_final), family = Gamma(link = "log"), data = mydat)

  
})
gam.lm
## save models
save(gam.lm, file = "output_data/00_csci_asci_all_gam_flow_combined.RData")

### get rsqds and pvals

for(i in 1:length(gam.lm)) {
  
  # if (class(gam.lm[[i]]) == "try-error") {
  #   
  #   # mod <- summary(log.lm[[i]])
  #   bio_h_summary$AIC[i] <- NA ##1-mod$deviance/mod$null.deviance ## mcfaddens r2
  #   bio_h_summary$PValue[i] <- NA
  #   bio_h_summary$McFaddensR2[i] <- NA
  #   bio_h_summary$n[i] <- NA
  #   

    mod <- summary(gam.lm[[i]])

    bio_h_summary$AIC[i] <- mod$aic ##1-mod$deviance/mod$null.deviance ## mcfaddens r2
    bio_h_summary$PValue[i] <- mod$anova$`Pr(F)`[2]
    bio_h_summary$McFaddensR2[i] <- 1-mod$deviance/mod$null.deviance
    # bio_h_summary$n[i] <- mod$df[2]+1
  
  
}


## save configs and r sqds
save(bio_h_summary, file="output_data/00_csci_asci_gam_rsqds_combined.RData")
bio_h_summary
## make df of predicted values to predict on - need to be different for each flow metric

## blank df
DF <- NULL
DF <- as.data.frame(DF)
AllDataLong2
i=1
### get predictions and fitted values
for(i in 1:length(gam.lm)) {
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  # dmet<-as.character(bio_h_summary[i,"Direction.alt"])
  cmet<-as.character(bio_h_summary[i,"channel.type"])


    
    ## get model, predict, extract all data and categories
    mod <- gam.lm[[i]]
    predictedVals <- predict(mod,  type = "response")
    
    DFX <- cbind(na.omit(mod$data), as.data.frame(predictedVals)) %>%
      rename(ChannelType = Class2, Value = deltah_final) %>%
      mutate(Variable = tmet) 
  
  DF <- bind_rows(DF, DFX)
  
}
DF
## change back to numeric
DF <- DF %>%
  mutate(Value = as.numeric(Value), MetricValue = as.numeric(MetricValue),
         predictedVals = as.numeric(predictedVals))
str(DF)

### predicted figures
bio <- unique(DF$Metric)
bio
mets <- unique(DF$Variable)
DF
m=1
for(b in 1:length(bio)) {
  
  DF1 <- DF %>%
    filter(Metric == bio[b])
  
  for(m in 1:length(mets)) {
    
    T1 <- (ggplot(subset(DF1, Variable == mets[m]), aes(y=predictedVals, x=Value, group = ChannelType, color = ChannelType)) +
             # geom_point(size=0.2) +
             geom_line( linewidth = 1)+
             # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
             # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
             # facet_wrap(~DirectionAlt, scales = "free") +
             scale_x_continuous(name=paste(mets[m])) +
             scale_y_continuous(name = paste0("Probability of Good ", bio[b]))) 
    # theme(legend.position = "none"))
    
    T1
    
    file.name1 <- paste0(out.dir, "00_", bio[b], "_", mets[m], "_flow_response_predicted_gam_combined.jpg")
    ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
    
  }
  
  
}



# GAM Models: separate channel types------------------------------------------------------------------

AllDataLong2 <- AllData %>%
  filter(Metric %in% c("csci", "ASCI_Hybrid"))

unique(AllDataLong2$flow_metric)

## create df of model configurations

## bio 
biol.endpoints<- unique(AllDataLong2$Metric)
biol.endpoints


## flow
flow.endpoints<- unique(na.omit(AllDataLong2$flow_metric))
# flow.endpoints <- flow.endpoints[-c(14:16)]


# channel type

channel.type.combined <- unique(na.omit(AllDataLong2$Class2))
channel.type.sep <- unique(na.omit(AllDataLong2$channel_engineering_class))


bio_h_summary<-  expand.grid(biol.endpoints=biol.endpoints,flow.endpoints=flow.endpoints, 
                             channel.type = channel.type.sep, stringsAsFactors = F)

# bio_h_summary <- bio_h_summary[-c(1:9),]
bio_h_summary
head(AllDataLong2)

## model of each configuration
gam.lm <-lapply(1:nrow(bio_h_summary), function(i)
{
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  # dmet<-as.character(bio_h_summary[i,"Direction.alt"])
  cmet<-as.character(bio_h_summary[i,"channel.type"])
  # thresh<-(bio_h_summary[i,"BioThresholds"])
  
  
  # if(dmet == "Pos") {
  #   
  #   mydat<-AllDataLong2 %>%
  #     filter(Metric == bmet,
  #            flow_metric == tmet, channel_engineering_class == cmet) %>%
  #     select(Metric, MetricValue, deltah_final, channel_engineering_class) %>% ## only metrics needed
  #     drop_na( deltah_final) %>%
  #     # filter(deltah_final > 0) %>%
  #     filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
  #     distinct()
  #   
  # } else {
    
    mydat<-AllDataLong2 %>%
      filter(Metric == bmet,
             flow_metric == tmet,  channel_engineering_class == cmet) %>%
      select(Metric, MetricValue, deltah_final, channel_engineering_class) %>% ## only metrics needed
      drop_na( deltah_final) %>%
      # filter(deltah_final < 0) %>%
      filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
      distinct()

  # mydat$Condition<-ifelse(mydat$MetricValue< thresh ,0, 1) ## convert to binary
  # mydat <- mydat %>% drop_na(Condition)
  
  mydat<-mydat[order(mydat$MetricValue),] ## order by csci value
  # try(glm(Condition~deltah_final, family=binomial(link="logit"), data=mydat), silent=T) ### glm
  
  try(gam(MetricValue~s(deltah_final), family = Gamma(link = "log"), data = mydat), silent = T)
  
})

gam.lm
## save models
save(gam.lm, file = "output_data/00_csci_asci_all_gam_flow_separate.RData")

### get rsqds and pvals
for(i in 1:length(gam.lm)) {
  
  if (class(gam.lm[[i]]) == "try-error") {
    
    # mod <- summary(log.lm[[i]])
    bio_h_summary$AIC[i] <- NA ##1-mod$deviance/mod$null.deviance ## mcfaddens r2
    bio_h_summary$PValue[i] <- NA
    bio_h_summary$McFaddensR2[i] <- NA
    # bio_h_summary$n[i] <- NA
    
  } else {
    
    mod <- summary(gam.lm[[i]])
    bio_h_summary$AIC[i] <- mod$aic ##1-mod$deviance/mod$null.deviance ## mcfaddens r2
    bio_h_summary$PValue[i] <- mod$anova$`Pr(F)`[2]
    bio_h_summary$McFaddensR2[i] <- 1-mod$deviance/mod$null.deviance
    # bio_h_summary$n[i] <- mod$df[2]+1
  }
  
}
## save configs and r sqds
save(bio_h_summary, file="output_data/00_csci_asci_gam_rsqds_separate.RData")
bio_h_summary
## make df of predicted values to predict on - need to be different for each flow metric

## blank df
DF <- NULL
DF <- as.data.frame(DF)

### get predictions and fitted values
for(i in 1:length(gam.lm)) {
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  dmet<-as.character(bio_h_summary[i,"Direction.alt"])
  cmet<-as.character(bio_h_summary[i,"channel.type"])
  
  if (class(gam.lm[[i]]) == "try-error") {
    ## get NAs where no data/model
    DFX <- cbind(bmet, as.numeric(NA), NA, cmet, NA,  tmet) %>% as.data.frame() %>%
      rename(Metric = bmet, MetricValue = V2, Value = V3, ChannelType = cmet,  
             predictedVals = V5,  Variable = tmet) 
  
  } else {
    
    ## get model, predict, extract all data and categories
    mod <- gam.lm[[i]]
    predictedVals <- predict(mod,  type = "response")
    
    DFX <- cbind(na.omit(mod$data), as.data.frame(predictedVals)) %>%
      rename(ChannelType = channel_engineering_class, Value = deltah_final) %>%
      mutate(Variable = tmet) %>%
      mutate(MetricValue = as.character(MetricValue), Value = as.character(Value),
             predictedVals = as.character(predictedVals))
    
  }
  
  DF <- bind_rows(DF, DFX)
  
}

DF
## change back to numeric
DF <- DF %>%
  mutate(Value = as.numeric(Value), MetricValue = as.numeric(MetricValue),
         predictedVals = as.numeric(predictedVals))
str(DF)

### predicted figures
bio <- unique(DF$Metric)
bio
mets <- unique(DF$Variable)
DF

b=1

for(b in 1:length(bio)) {
  
  DF1 <- DF %>%
    filter(Metric == bio[b])
  
  
  for(m in 1:length(mets)) {
    
    T1 <- (ggplot(subset(DF1, Variable == mets[m]), aes(y=predictedVals, x=Value, group = ChannelType, color = ChannelType)) +
             # geom_point(size=0.2) +
             geom_line( linewidth = 1)+
             # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
             # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
             # facet_wrap(~DirectionAlt, scales = "free") +
             scale_x_continuous(name=paste(mets[m])) +
             scale_y_continuous(name = paste0("Probability of Good ", bio[b]))) 
    # theme(legend.position = "none"))
    
    T1
    
    file.name1 <- paste0(out.dir, "00_", bio[b], "_", mets[m], "_flow_response_predicted_gam_separate.jpg")
    ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
  }
}


