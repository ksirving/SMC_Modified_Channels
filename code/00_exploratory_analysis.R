## exploring possible appraoches to flow eco in modified channels

library(tidylog)
library(tidyverse)
library(sf)

out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/SMC_Modified_Channels/Figures/"

getwd()


# Flow data ---------------------------------------------------------------

dh_data <- read.csv("/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/SMC_Modified_Channels/ignore/SoCal_bio_deltaH_summary_supp_final.csv")
head(dh_data)

dh_median <- dh_data %>%
  filter(summary.statistic =="median") %>%
  select(site, flow_metric, deltah_final) %>%
  rename(masterid= site)
dim(dh_median) # 7236
head(dh_median)

## full names for labels
labels <- read.csv("input_data/ffm_names.csv")
labels <- labels[1:24, ]
labels <- labels %>% rename(Hydro_endpoint = Flow.Metric.Code)
labels[25, 1] <- "Magnitude of largest annual storm"
labels[25, 2] <- "Q99"
labels[25, 3] <- "Peak Flow"
labels

limits <- left_join(limits, labels, by = "Hydro_endpoint")

## sites with flow data

flowSites <- unique(dh_median$masterid)

flowSites

## channel engineering data

BioEng <- read.csv("ignore/02_chan_eng.csv") %>%
  select(-c(X,channel_engineering_classification_date, channel_engineering_personnel, channel_engineering_comments)) %>%
  mutate(Class2 = ifelse(channel_engineering_class =="NAT", "Natural", "Modified"))

BioEng

## count number of asci/csci sites per channel class
## boxplots of delta h for all channel types
## test mod streams with derived delta h limits

# Bio data ----------------------------------------------------------------

##  sites only

bioSites <- st_read("input_data/01_bio_sites_all.shp")
head(bioSites)

## bugs data

csciScores <- read.csv("ignore/01_csci_comp_mets_comid_socal.csv")
head(csciScores)

### algae scores

asciScores <- read.csv("ignore/01_asci_comp_mets_comid_socal.csv")
head(asciScores)


# Join bio and flow -------------------------------------------------------

## how many flow in bug sites
sum(flowSites %in% bioSites$masterid) ## 422

sum(flowSites %in% csciScores$masterid) ## 402

## filter bug sites to flow sites
bugflowSites <- bioSites %>%
  filter(masterid %in% flowSites)

## filter bug data using masterid ### remove reps - remove 2nd rep for now, change later!!!!
csciScoresLA <- csciScores %>%
  select(-X, -stationcode) %>%
  filter(masterid %in% bugflowSites$masterid, fieldreplicate == 1 ) %>%
  separate(sampledate, into = c("sampledate", "Time"), sep= "T", remove = F) %>%
  separate(sampledate, into = c("year", "Month", "Day"), sep= "-", remove = F) %>%
  mutate(year = as.numeric(year))

length(unique(csciScoresLA$masterid)) ## 705 sites in LA region with flow


## filter bug data using masterid ### remove reps - remove 2nd rep for now, change later!!!!
asciScoresLA <- asciScores %>%
  select(-X, -stationcode) %>%
  filter(masterid %in% bugflowSites$masterid,replicate == 1 ) %>%
  separate(sampledate, into = c("sampledate", "Time"), sep= "T", remove = F) %>%
  separate(sampledate, into = c("year", "Month", "Day"), sep= "-", remove = F) %>%
  mutate(year = as.numeric(year))

length(unique(asciScoresLA$masterid)) ## 336 sites in LA region with flow


# Join bio sites to flow data ---------------------------------------------

scoresflowSites <- csciScoresLA 


AllData <- left_join(scoresflowSites, dh_median, by = c("masterid"))  %>%
  left_join(BioEng, by = "masterid")
head(AllData)
dim(AllData) #

meancsci <- AllData %>%
  filter(Metric == "csci") %>%
  summarise(meancsci = mean(na.omit(MetricValue)))
meancsci

## save out
save(AllData, file = "output_data/00_bugs_flow_joined_by_masterid.RData")

# AllDataA <- left_join(asciScoresLA,  dh_median, by = c("masterid")) 
# head(AllDataA)
# dim(AllDataA) ## 
# 
# ## save out
# save(AllDataA, file = "output_data/00_algae_flow_joined_by_masterid.RData")


# GLM Models: CSCI ------------------------------------------------------------------

AllDataLong2 <- AllData %>%
  filter(!Metric == "count")

unique(AllDataLong2$flow_metric)
## create df of model configurations

## bio 
biol.endpoints<-"csci"

## flow
flow.endpoints<- unique(na.omit(AllDataLong2$flow_metric))
flow.endpoints

# Thresholds for index
index.thresholds <- c(0.79)

# direction of alteration

Direction.alt <- c("Pos", "Neg")
Direction.alt

# channel type

channel.type <- unique(na.omit(AllDataLong2$Class2))
channel.type


bio_h_summary<-  expand.grid(biol.endpoints=biol.endpoints,flow.endpoints=flow.endpoints, Direction.alt= Direction.alt, channel.type = channel.type, stringsAsFactors = F)

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
  


  mydat$Condition<-ifelse(mydat$MetricValue< 0.79 ,0, 1) ## convert to binary
  # mydat <- mydat %>% drop_na(Condition)
  
  mydat<-mydat[order(mydat$MetricValue),] ## order by csci value
  glm(Condition~deltah_final, family=binomial(link="logit"), data=mydat) ### glm
  
  
})
log.lm
## save models
save(log.lm, file = "output_data/00_csci_all_glm_flow.RData")

### get rsqds and pvals
for(i in 1:length(log.lm)) {
  
  mod <- summary(log.lm[[i]])
  bio_h_summary$AIC[i] <- mod$aic ##1-mod$deviance/mod$null.deviance ## mcfaddens r2
  bio_h_summary$PValue[i] <- mod$coefficients[8]
  bio_h_summary$McFaddensR2[i] <- 1-mod$deviance/mod$null.deviance
  bio_h_summary$n[i] <- mod$df[2]+1
  
}
## save configs and r sqds
save(bio_h_summary, file="output_data/00_csci_glm_rsqds.RData")
bio_h_summary
## make df of predicted values to predict on - need to be different for each flow metric

## blank df
DF <- NULL
DF <- as.data.frame(DF)
AllDataLong2
### get predictions and fitted values
for(i in 1:length(log.lm)) {
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  dmet<-as.character(bio_h_summary[i,"Direction.alt"])
  cmet<-as.character(bio_h_summary[i,"channel.type"])
  tmet
  
  ## new data - plus and minus 10% as a start
  # flowvalues <- seq(range(data$Curflow)[1]-10,range(data$Curflow)[2]+10,0.05)
  
  ## get model, predict, extract all data and categories
  mod <- log.lm[[i]]
  predictedVals <- predict.glm(mod,  type = "response")
  
  DFX <- cbind(mod$data, as.data.frame(predictedVals)) %>%
    rename(ChannelType = Class2, Value = deltah_final) %>%
    mutate(DirectionAlt = dmet, Variable = tmet)

 
  
  DF <- bind_rows(DF, DFX)
  
}
dmet

DF
### predicted figures
mets <- unique(DF$Variable)
mets
m=1

## facet labels



for(m in 1:length(mets)) {
  
  T1 <- (ggplot(subset(DF, Variable == mets[m]), aes(y=predictedVals, x=Value, group = ChannelType, color = ChannelType)) +
           # geom_point(size=0.2) +
           geom_line( linewidth = 1)+
           # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
           # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
           facet_wrap(~DirectionAlt, scales = "free") +
           scale_x_continuous(name=paste(mets[m])) +
           scale_y_continuous(name = paste0("Probability of 0.79 CSCI"))) 
           # theme(legend.position = "none"))
  
  T1
  
  file.name1 <- paste0(out.dir, "00_", mets[m], "_csci_flow_response_predicted_glm.jpg")
  ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
}

### predicted figures with index thresholds on one figure
mets <- unique(DF$Variable)
mets[m]

for(m in 1:length(mets)) {
  
  T1 <- (ggplot(subset(DF, Variable == mets[m]), aes(y=predictedVals, x=Value, group = BioThreshold, color = BioThreshold)) +
           # geom_point(size=0.2) +
           geom_line(linewidth = 1)+
           # scale_color_manual(values=c('blue','red')) +
           # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
           # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
           # facet_wrap(~Variable, labeller = as_labeller(supp.labs),
           #            scales = "free") +
           # geom_vline(data=filter(DF, !Variable %in% c("Max_Wkl_Max_StreamT_grt_30_","Max_Wkl_Rng_StreamT", "Mean_Wkl_Rng_StreamT")), 
           #            aes(xintercept = 86), linetype="dashed", color = "red", linewidth=0.5, show.legend = T) +
           geom_vline(data=filter(DF, !Variable %in% c("Max_Wkl_Max_StreamT_grt_30_","Max_Wkl_Rng_StreamT", "Mean_Wkl_Rng_StreamT")),
                      aes(xintercept = 80), linetype="dashed", color = "blue", linewidth=0.5) +
           scale_x_continuous(name="Water flow (°F)") +
           scale_y_continuous(name = paste0("Probability of Good CSCI")))
  
  T1
  
  file.name1 <- paste0(out.dir, "03_", mets[m], "_csci_flow_response_predicted_glm.jpg")
  ggsave(T1, filename=file.name1, dpi=300, height=6, width=6)
}

## get flows for different probabilities

head(DF)

MaxDFCSCI <- DF %>%
  filter(Variable == "Max_Wkl_Max_StreamT")

MaxDF

# Models: ASCI -----------------------------------------------------------------

AllDataLong2 <- AllDataA %>%
  filter(!Metric == "count")

unique(AllDataLong2$CurMetric)
## create df of model configurations
## bio 
biol.endpoints<-"ASCI_Hybrid"

## flow
flow.endpoints<- unique(na.omit(AllDataLong2$CurMetric))

# Thresholds for index
index.thresholds <- c(0.75, 0.86)

bio_h_summary<-  expand.grid(biol.endpoints=biol.endpoints,flow.endpoints=flow.endpoints, 
                             index.thresholds= index.thresholds, stringsAsFactors = F)
bio_h_summary
## model of each configuration
log.lm <-lapply(1:nrow(bio_h_summary), function(i)
{
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  imet<-as.character(bio_h_summary[i,"index.thresholds"])
  
  mydat<-AllDataLong2 %>%
    filter(Metric == bmet,
           CurMetric == tmet) %>%
    select(Curflow, MetricValue, COMID) %>% ## only metrics needed
    drop_na(Curflow, MetricValue) %>%
    filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
    distinct()
  
  mydat$Condition<-ifelse(mydat$MetricValue< imet ,0,1) ## convert to binary
  # mydat <- mydat %>% drop_na(Condition)
  
  mydat<-mydat[order(mydat$MetricValue),] ## order by csci value
  glm(Condition~Curflow, family=binomial(link="logit"), data=mydat) ### glm
  
  
})

## save models
save(log.lm, file = "output_data/03_asci_glm_currentflow.RData")

### get rsqds and pvals
for(i in 1:length(log.lm)) {
  
  mod <- summary(log.lm[[i]])
  bio_h_summary$AIC[i] <- mod$aic ##1-mod$deviance/mod$null.deviance ## mcfaddens r2
  bio_h_summary$PValue[i] <- mod$coefficients[8]
  bio_h_summary$McFaddensR2[i] <- 1-mod$deviance/mod$null.deviance
  bio_h_summary$n[i] <- mod$df[2]+1
  
}
## save configs and r sqds
save(bio_h_summary, file="output_data/03_asci_glm_rsqds.RData")
bio_h_summary
## make df of predicted values to predict on - need to be different for each flow metric

## blank df
DF <- NULL
DF <- as.data.frame(DF)

### get predictions and fitted values
for(i in 1:length(log.lm)) {
  
  bio <- bio_h_summary[i,"biol.endpoints"]
  flow <- bio_h_summary[i,"flow.endpoints"]
  ind <- bio_h_summary[i,"index.thresholds"]
  
  data <- AllDataLong2 %>%
    filter(Metric == bio,
           CurMetric == flow) %>%
    select(Curflow, MetricValue, COMID) %>% ## only metrics needed
    drop_na(Curflow, MetricValue) %>%
    filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
    distinct()
  
  
  ## new data - plus and minus 10% as a start
  flowvalues <- seq(range(data$Curflow)[1]-10,range(data$Curflow)[2]+10,0.05)
  
  ## get model, predict, extract all data and categories
  mod <- log.lm[[i]]
  predictedVals <- predict.glm(mod,list(Curflow = flowvalues),  type = "response")
  DFX <- as.data.frame(predictedVals)
  DFX$Value <- flowvalues
  DFX$Bio <- bio
  DFX$Variable <- flow
  DFX$BioThreshold <- ind
  DFX$MinVal <-  range(data$Curflow)[1]
  DFX$MaxVal <-  range(data$Curflow)[2]
  
  DF <- bind_rows(DF, DFX)
  
}

DF$BioThreshold <- as.factor(DF$BioThreshold)
# plot(DF$predictedVals, DF$BioVals)

### predicted figures
mets <- unique(DF$BioThreshold)

## facet labels
supp.labs <- c(
  "Max_Wkl_Max_StreamT_grt_30_"="Weeks greater than 86F",
  "Max_Wkly_Mean_StreamT" = "Max Weekly Mean",
  "Max_Wkl_Max_StreamT" = "Weekly Maximum",
  "Min_Wkl_Min_StreamT" = "Weekly Minimum",
  "Max_Wkl_Rng_StreamT" = " Max Weekly Range",
  "Mean_Wkl_Rng_StreamT" =  "Av Weekly Range" 
)

for(m in 1:length(mets)) {
  
  T1 <- (ggplot(subset(DF, BioThreshold == mets[m]), aes(y=predictedVals, x=Value, group = Variable, color = Variable)) +
           # geom_point(size=0.2) +
           geom_line( linewidth = 1)+
           # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
           # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
           facet_wrap(~Variable, labeller = as_labeller(supp.labs),
                      scales = "free") +
           # geom_vline(data=filter(DF, !Variable %in% c("Max_Wkl_Max_StreamT_grt_30_","Max_Wkl_Rng_StreamT", "Mean_Wkl_Rng_StreamT")), 
           #            aes(xintercept = 86), linetype="dashed", color = "red", linewidth=0.5, show.legend = T) +
           geom_vline(data=filter(DF, !Variable %in% c("Max_Wkl_Max_StreamT_grt_30_","Max_Wkl_Rng_StreamT", "Mean_Wkl_Rng_StreamT")), 
                      aes(xintercept = 80), linetype="dashed", color = "blue", linewidth=0.5, show.legend = T) +
           scale_x_continuous(name="Water flow (°F)") +
           scale_y_continuous(name = paste0("Probability of ", mets[m], " CSCI")) +
           theme(legend.position = "none"))
  
  T1
  
  file.name1 <- paste0(out.dir, "03_", mets[m], "_asci_flow_response_predicted_glm.jpg")
  ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
}

### predicted figures with index thresholds on one figure
mets <- unique(DF$Variable)
mets[m]

for(m in 1:length(mets)) {
  
  T1 <- (ggplot(subset(DF, Variable == mets[m]), aes(y=predictedVals, x=Value, group = BioThreshold, color = BioThreshold)) +
           # geom_point(size=0.2) +
           geom_line( linewidth = 1, aes( color = BioThreshold) )+
           # scale_color_manual(values=c('blue','red')) +
           # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
           # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
           # facet_wrap(~Variable, labeller = as_labeller(supp.labs),
           #            scales = "free") +
           # geom_vline(data=filter(DF, !Variable %in% c("Max_Wkl_Max_StreamT_grt_30_","Max_Wkl_Rng_StreamT", "Mean_Wkl_Rng_StreamT")), 
           #            aes(xintercept = 86), linetype="dashed", color = "red", linewidth=0.5, show.legend = T) +
           geom_vline(data=filter(DF, !Variable %in% c("Max_Wkl_Max_StreamT_grt_30_","Max_Wkl_Rng_StreamT", "Mean_Wkl_Rng_StreamT")),
                      aes(xintercept = 80), linetype="dashed", color = "blue", linewidth=0.5) +
           scale_x_continuous(name="Water flow (°F)") +
           scale_y_continuous(name = paste0("Probability of Good ASCI")) +
           theme(legend.position = "none"))
  
  T1
  
  file.name1 <- paste0(out.dir, "03_", mets[m], "_asci_flow_response_predicted_glm.jpg")
  ggsave(T1, filename=file.name1, dpi=300, height=6, width=6)
}

## get flows for different probabilities

head(DF)

MaxDFASCI <- DF %>%
  filter(Variable == "Max_Wkl_Max_StreamT")

MaxDFASCI

## join with csci

MaxDF <- rbind(MaxDFASCI, MaxDFCSCI)

write.csv(MaxDF, "ignore/03_probs_per_flow.csv")


# Alteration models: CSCI -------------------------------------------------

AllDataLong2 <- AllData %>%
  filter(!Metric == "count") #%>%
# pivot_longer(AltflowMinus80:AltflowDivide80, names_to="AltMetric", values_to = "AltValues")

names(AllDataLong2)

## create df of model configurations

## bio 
biol.endpoints<-"csci"

## flow
flow.endpoints<- unique(na.omit(AllDataLong2$CurMetric))[2:4]
flow.endpoints

## Alt flow
alt.endpoints<- unique(na.omit(AllDataLong2$AltMetric))
alt.endpoints


bio_h_summary<-  expand.grid(biol.endpoints=biol.endpoints,flow.endpoints=flow.endpoints,  alt.endpoints= alt.endpoints, stringsAsFactors = F)

# bio_h_summary <- bio_h_summary[-c(1:9),]

## model of each configuration
log.glm <-lapply(1:nrow(bio_h_summary), function(i)
{
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  amet<-as.character(bio_h_summary[i,"alt.endpoints"])
  
  mydat<-AllDataLong2 %>%
    filter(Metric == bmet,
           CurMetric == tmet, 
           AltMetric == amet) %>%
    select(Altflow, MetricValue, COMID) %>% ## only metrics needed
    drop_na(Altflow, MetricValue) %>%
    filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
    distinct() #%>%
  # filter(Altflow > 0) ## take only positive altered flows
  
  mydat$Condition<-ifelse(mydat$MetricValue< 0.60 ,0,1) ## convert to binary
  # mydat <- mydat %>% drop_na(Condition)
  
  mydat<-mydat[order(mydat$MetricValue),] ## order by csci value
  glm(Condition~Altflow, family=binomial(link="logit"), data=mydat) ### glm
  
  # View(mydat)
})


## save models
save(log.glm, file = "output_data/03_csci_glm_altered.RData")

### get rsqds and pvals
for(i in 1:length(log.glm)) {
  
  mod <- summary(log.glm[[i]])
  bio_h_summary$AIC[i] <- mod$aic 
  bio_h_summary$PValue[i] <- mod$coefficients[8]
  bio_h_summary$McFaddensR2[i] <- 1-mod$deviance/mod$null.deviance
  # mod
}
## save configs and r sqds
save(bio_h_summary, file="output_data/03_csci_altered_glm_rsqds.RData")


## make df of predicted values to predict on - need to be different for each flow metric
## blank df
DF <- NULL
DF <- as.data.frame(DF)
bio_h_summary
### get predictions and fitted values
for(i in 1:length(log.glm)) {
  
  bio <- bio_h_summary[i,"biol.endpoints"]
  flow <- bio_h_summary[i,"flow.endpoints"]
  alt <- bio_h_summary[i,"alt.endpoints"]
  
  data <- AllDataLong2 %>%
    filter(Metric == bio,
           CurMetric == flow,
           AltMetric == alt) %>%
    select(Altflow, MetricValue, COMID) %>% ## only metrics needed
    drop_na(Altflow, MetricValue) %>%
    filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
    distinct() 
  
  ## new data - plus and minus 10 as a start
  
  if (alt == "AltflowDivide80" ) {
    flowvalues <- seq(range(data$Altflow)[1]*0.9,range(data$Altflow)[2]*1.1,0.05)
  } else {
    flowvalues <- seq(range(data$Altflow)[1]-10,range(data$Altflow)[2]+10,0.05)
  }
  
  
  ## get model, predict, extract all data and categories
  mod <- log.glm[[i]]
  predictedVals <- predict.glm(mod,list(Altflow = flowvalues),  type = "response")
  DFX <- as.data.frame(predictedVals)
  DFX$Value <- flowvalues
  DFX$Bio <- bio
  DFX$Variable <- flow
  DFX$AltVar <- alt
  DFX$MinVal <-  range(data$Altflow)[1]
  DFX$MaxVal <-  range(data$Altflow)[2]
  
  DF <- bind_rows(DF, DFX)
  
}


## remove mean weekly as not an expected relationship
## and scale probability?

DF <- DF %>%
  filter(!Variable == "Max_Wkly_Mean_StreamT")


### predicted figures

mets <- unique(DF$AltVar)

## facet labels
supp.labs <- c(
  "Max_Wkl_Max_StreamT_grt_30_"="Weeks greater than 86F",
  "Max_Wkly_Mean_StreamT" = "Max Weekly Mean",
  "Max_Wkl_Max_StreamT" = "Weekly Maximum",
  "Min_Wkl_Min_StreamT" = "Weekly Minimum",
  "Max_Wkl_Rng_StreamT" = " Max Weekly Range",
  "Mean_Wkl_Rng_StreamT" =  "Av Weekly Range" 
)

for(m in 1:length(mets)) {
  
  T1 <- (ggplot(subset(DF, AltVar == mets[m]), aes(y=predictedVals, x=Value, group = Variable, color = Variable)) +
           # geom_point(size=0.2) +
           geom_line( linewidth = 1)+
           # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
           # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
           facet_wrap(~Variable, labeller = as_labeller(supp.labs),
                      scales = "free") +
           # geom_vline(data=filter(DF, !Variable %in% c("Max_Wkl_Max_StreamT_grt_30_","Max_Wkl_Rng_StreamT", "Mean_Wkl_Rng_StreamT")), 
           #            aes(xintercept = 86), linetype="dashed", color = "red", linewidth=0.5, show.legend = T) +
           # geom_vline(data=filter(DF, !Variable %in% c("Max_Wkl_Max_StreamT_grt_30_","Max_Wkl_Rng_StreamT", "Mean_Wkl_Rng_StreamT")), 
           #            aes(xintercept = 80), linetype="dashed", color = "blue", linewidth=0.5, show.legend = T) +
           scale_x_continuous(name=paste0(mets[m],": flowerature Alteration (°F)")) +
           scale_y_continuous(name = "Probability of good CSCI") +
           theme(legend.position = "none"))
  
  T1
  
  file.name1 <- paste0(out.dir, "03_", mets[m], "_csci_flow_response_predicted_glm.jpg")
  ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
}



# Alteration Models: ASCI -------------------------------------------------

AllDataLong2 <- AllDataA %>%
  filter(!Metric == "count") #%>%
# pivot_longer(AltflowMinus80:AltflowDivide80, names_to="AltMetric", values_to = "AltValues")

## create df of model configurations

## bio 
biol.endpoints<-"ASCI_Hybrid"

## flow
flow.endpoints<- unique(na.omit(AllDataLong2$CurMetric))[2:4]
flow.endpoints

## Alt flow
alt.endpoints<- unique(na.omit(AllDataLong2$AltMetric))
alt.endpoints


bio_h_summary<-  expand.grid(biol.endpoints=biol.endpoints,flow.endpoints=flow.endpoints,  alt.endpoints= alt.endpoints, stringsAsFactors = F)

# bio_h_summary <- bio_h_summary[-c(1:9),]
bio_h_summary
i=1
i
## model of each configuration
log.glm <-lapply(1:nrow(bio_h_summary), function(i)
{
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  amet<-as.character(bio_h_summary[i,"alt.endpoints"])
  
  mydat<-AllDataLong2 %>%
    filter(Metric == bmet,
           CurMetric == tmet, 
           AltMetric == amet) %>%
    select(Altflow, MetricValue, COMID) %>% ## only metrics needed
    drop_na(Altflow, MetricValue) %>%
    filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
    distinct() #%>%
  # filter(Altflow > 0) ## take only positive altered flows
  
  mydat$Condition<-ifelse(mydat$MetricValue< 0.75 ,0,1) ## convert to binary
  # mydat <- mydat %>% drop_na(Condition)
  
  mydat<-mydat[order(mydat$MetricValue),] ## order by csci value
  glm(Condition~Altflow, family=binomial(link="logit"), data=mydat) ### glm
  
  # View(mydat)
})


## save models
save(log.glm, file = "output_data/03_asci_glm_altered.RData")

### get rsqds and pvals
for(i in 1:length(log.glm)) {
  
  mod <- summary(log.glm[[i]])
  bio_h_summary$AIC[i] <- mod$aic 
  bio_h_summary$PValue[i] <- mod$coefficients[8]
  bio_h_summary$McFaddensR2[i] <- 1-mod$deviance/mod$null.deviance
  # mod
}
## save configs and r sqds
save(bio_h_summary, file="output_data/03_asci_altered_glm_rsqds.RData")


## make df of predicted values to predict on - need to be different for each flow metric
## blank df
DF <- NULL
DF <- as.data.frame(DF)

### get predictions and fitted values
for(i in 1:length(log.glm)) {
  
  bio <- bio_h_summary[i,"biol.endpoints"]
  flow <- bio_h_summary[i,"flow.endpoints"]
  alt <- bio_h_summary[i,"alt.endpoints"]
  
  data <- AllDataLong2 %>%
    filter(Metric == bio,
           CurMetric == flow,
           AltMetric == alt) %>%
    select(Altflow, MetricValue, COMID) %>% ## only metrics needed
    drop_na(Altflow, MetricValue) %>%
    filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
    distinct() 
  
  ## new data - plus and minus 10% as a start
  if (alt == "AltflowDivide80" ) {
    flowvalues <- seq(range(data$Altflow)[1]*0.9,range(data$Altflow)[2]*1.1,0.05)
  } else {
    flowvalues <- seq(range(data$Altflow)[1]-10,range(data$Altflow)[2]+10,0.05)
  }
  ## get model, predict, extract all data and categories
  mod <- log.glm[[i]]
  predictedVals <- predict.glm(mod,list(Altflow = flowvalues),  type = "response")
  DFX <- as.data.frame(predictedVals)
  DFX$Value <- flowvalues
  DFX$Bio <- bio
  DFX$Variable <- flow
  DFX$AltVar <- alt
  DFX$MinVal <-  range(data$Altflow)[1]
  DFX$MaxVal <-  range(data$Altflow)[2]
  
  DF <- bind_rows(DF, DFX)
  
}


## remove mean weekly as not an expected relationship
## and scale probability?

DF <- DF %>%
  filter(!Variable == "Max_Wkly_Mean_StreamT")


### predicted figures

mets <- unique(DF$AltVar)

## facet labels
supp.labs <- c(
  "Max_Wkl_Max_StreamT_grt_30_"="Weeks greater than 86F",
  "Max_Wkly_Mean_StreamT" = "Max Weekly Mean",
  "Max_Wkl_Max_StreamT" = "Weekly Maximum",
  "Min_Wkl_Min_StreamT" = "Weekly Minimum",
  "Max_Wkl_Rng_StreamT" = " Max Weekly Range",
  "Mean_Wkl_Rng_StreamT" =  "Av Weekly Range" 
)

for(m in 1:length(mets)) {
  
  T1 <- (ggplot(subset(DF, AltVar == mets[m]), aes(y=predictedVals, x=Value, group = Variable, color = Variable)) +
           # geom_point(size=0.2) +
           geom_line( linewidth = 1)+
           # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
           # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
           facet_wrap(~Variable, labeller = as_labeller(supp.labs),
                      scales = "free") +
           # geom_vline(data=filter(DF, !Variable %in% c("Max_Wkl_Max_StreamT_grt_30_","Max_Wkl_Rng_StreamT", "Mean_Wkl_Rng_StreamT")), 
           #            aes(xintercept = 86), linetype="dashed", color = "red", linewidth=0.5, show.legend = T) +
           # geom_vline(data=filter(DF, !Variable %in% c("Max_Wkl_Max_StreamT_grt_30_","Max_Wkl_Rng_StreamT", "Mean_Wkl_Rng_StreamT")), 
           #            aes(xintercept = 80), linetype="dashed", color = "blue", linewidth=0.5, show.legend = T) +
           scale_x_continuous(name=paste0(mets[m],": flowerature Alteration (°F)")) +
           scale_y_continuous(name = "Probability of median ASCI") +
           theme(legend.position = "none"))
  
  T1
  
  file.name1 <- paste0(out.dir, "03_", mets[m], "_asci_flow_response_predicted_glm.jpg")
  ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
}



# Channel engineering -----------------------------------------------

## join engineering data to bio

AllDataLong2 <- AllData %>%
  filter(!Metric == "count") %>%
  left_join(BioEng, by = "masterid")

head(AllDataLong2)

unique(AllDataLong2$CurMetric)
## create df of model configurations

## bio 
biol.endpoints<-"csci"

## flow
flow.endpoints<- unique(na.omit(AllDataLong2$CurMetric))
flow.endpoints

## engineering
eng.endpoints <- unique(na.omit(AllDataLong2$Class2))

bio_h_summary<-  expand.grid(biol.endpoints=biol.endpoints,flow.endpoints=flow.endpoints, eng.endpoints=eng.endpoints,stringsAsFactors = F)

# bio_h_summary <- bio_h_summary[-c(1:9),]
bio_h_summary
head(AllDataLong2)
## model of each configuration
log.lm <-lapply(1:nrow(bio_h_summary), function(i)
{
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  emet<-as.character(bio_h_summary[i,"eng.endpoints"])
  
  mydat<-AllDataLong2 %>%
    filter(Metric == bmet,
           CurMetric == tmet,
           Class2 == emet) %>%
    select(Curflow, MetricValue, COMID) %>% ## only metrics needed
    drop_na(Curflow, MetricValue) %>%
    filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
    distinct()
  mydat
  mydat$Condition<-ifelse(mydat$MetricValue< 0.79 ,0,1) ## convert to binary
  # mydat <- mydat %>% drop_na(Condition)
  
  mydat<-mydat[order(mydat$MetricValue),] ## order by csci value
  glm(Condition~Curflow, family=binomial(link="logit"), data=mydat) ### glm
  
  
})

## save models
save(log.lm, file = "output_data/03_csci_glm_currentflow_chan_eng.RData")

### get rsqds and pvals
for(i in 1:length(log.lm)) {
  
  mod <- summary(log.lm[[i]])
  bio_h_summary$AIC[i] <- mod$aic ##1-mod$deviance/mod$null.deviance ## mcfaddens r2
  bio_h_summary$PValue[i] <- mod$coefficients[8]
  bio_h_summary$McFaddensR2[i] <- 1-mod$deviance/mod$null.deviance
  bio_h_summary$n[i] <- mod$df[2]+1
  
}
## save configs and r sqds
save(bio_h_summary, file="output_data/03_csci_glm_rsqds_chan_eng.RData")
bio_h_summary
## make df of predicted values to predict on - need to be different for each flow metric

## blank df
DF <- NULL
DF <- as.data.frame(DF)

### get predictions and fitted values
for(i in 1:length(log.lm)) {
  
  bio <- bio_h_summary[i,"biol.endpoints"]
  flow <- bio_h_summary[i,"flow.endpoints"]
  eng <- bio_h_summary[i,"eng.endpoints"]
  
  data <- AllDataLong2 %>%
    filter(Metric == bio,
           CurMetric == flow,
           Class2 == eng) %>%
    select(Curflow, MetricValue, COMID) %>% ## only metrics needed
    drop_na(Curflow, MetricValue) %>%
    filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
    distinct()
  
  
  ## new data - plus and minus 10% as a start
  flowvalues <- seq(range(data$Curflow)[1]-10,range(data$Curflow)[2]+10,0.05)
  
  ## get model, predict, extract all data and categories
  mod <- log.lm[[i]]
  predictedVals <- predict.glm(mod,list(Curflow = flowvalues),  type = "response")
  DFX <- as.data.frame(predictedVals)
  DFX$Value <- flowvalues
  DFX$Bio <- bio
  DFX$Variable <- flow
  DFX$Eng <- eng
  DFX$MinVal <-  range(data$Curflow)[1]
  DFX$MaxVal <-  range(data$Curflow)[2]
  
  DF <- bind_rows(DF, DFX)
  
}
DF
# plot(DF$predictedVals, DF$BioVals)

### predicted figures
mets <- unique(DF$Eng)
mets
m=1
## facet labels
supp.labs <- c(
  "Max_Wkl_Max_StreamT_grt_30_"="Weeks greater than 86F",
  "Max_Wkly_Mean_StreamT" = "Max Weekly Mean",
  "Max_Wkl_Max_StreamT" = "Weekly Maximum",
  "Min_Wkl_Min_StreamT" = "Weekly Minimum",
  "Max_Wkl_Rng_StreamT" = " Max Weekly Range",
  "Mean_Wkl_Rng_StreamT" =  "Av Weekly Range" 
)


for(m in 1:length(mets)) {
  
  T1 <- (ggplot(subset(DF, Eng == mets[m]), aes(y=predictedVals, x=Value, group = Eng, color = Eng)) +
           # geom_point(size=0.2) +
           geom_line( linewidth = 1)+
           # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
           # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
           facet_wrap(~Variable, labeller = as_labeller(supp.labs),
                      scales = "free") +
           # geom_vline(data=filter(DF, !Variable %in% c("Max_Wkl_Max_StreamT_grt_30_","Max_Wkl_Rng_StreamT", "Mean_Wkl_Rng_StreamT")), 
           #            aes(xintercept = 86), linetype="dashed", color = "red", linewidth=0.5) +
           geom_vline(data=filter(DF, !Variable %in% c("Max_Wkl_Max_StreamT_grt_30_","Max_Wkl_Rng_StreamT", "Mean_Wkl_Rng_StreamT")), 
                      aes(xintercept = 80), linetype="dashed", color = "blue", linewidth=0.5) +
           scale_x_continuous(name="Water flow (°F)") +
           scale_y_continuous(name =  "Probability of median CSCI")) +
    theme(legend.title=element_blank())
  
  T1
  
  file.name1 <- paste0(out.dir, "03_csci_", mets[m], "_flow_response_predicted_glm_ModChannels.jpg")
  ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
  
}




# Models: ASCI -----------------------------------------------------------------

AllDataLong2 <- AllDataA %>%
  filter(!Metric == "count")%>%
  left_join(BioEng, by = "masterid")

head(AllDataLong2)

unique(AllDataLong2$CurMetric)
## create df of model configurations

## bio 
biol.endpoints<-"ASCI_Hybrid"

## flow
flow.endpoints<- unique(na.omit(AllDataLong2$CurMetric))
flow.endpoints

## engineering
eng.endpoints <- unique(na.omit(AllDataLong2$Class2))

bio_h_summary<-  expand.grid(biol.endpoints=biol.endpoints,flow.endpoints=flow.endpoints, eng.endpoints=eng.endpoints,stringsAsFactors = F)

# bio_h_summary <- bio_h_summary[-c(1:9),]
bio_h_summary
head(AllDataLong2)
## model of each configuration
log.lm <-lapply(1:nrow(bio_h_summary), function(i)
{
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  emet<-as.character(bio_h_summary[i,"eng.endpoints"])
  
  mydat<-AllDataLong2 %>%
    filter(Metric == bmet,
           CurMetric == tmet,
           Class2 == emet) %>%
    select(Curflow, MetricValue, COMID) %>% ## only metrics needed
    drop_na(Curflow, MetricValue) %>%
    filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
    distinct()
  
  mydat$Condition<-ifelse(mydat$MetricValue< 0.86 ,0,1) ## convert to binary
  # mydat <- mydat %>% drop_na(Condition)
  
  mydat<-mydat[order(mydat$MetricValue),] ## order by csci value
  glm(Condition~Curflow, family=binomial(link="logit"), data=mydat) ### glm
  
  
})

## save models
save(log.lm, file = "output_data/03_asci_glm_currentflow_chan_eng.RData")

### get rsqds and pvals
for(i in 1:length(log.lm)) {
  
  mod <- summary(log.lm[[i]])
  bio_h_summary$AIC[i] <- mod$aic ##1-mod$deviance/mod$null.deviance ## mcfaddens r2
  bio_h_summary$PValue[i] <- mod$coefficients[8]
  bio_h_summary$McFaddensR2[i] <- 1-mod$deviance/mod$null.deviance
  
}
## save configs and r sqds
save(bio_h_summary, file="output_data/03_asci_glm_rsqds_chan_eng.RData")
bio_h_summary
## make df of predicted values to predict on - need to be different for each flow metric

## blank df
DF <- NULL
DF <- as.data.frame(DF)

### get predictions and fitted values
for(i in 1:length(log.lm)) {
  
  bio <- bio_h_summary[i,"biol.endpoints"]
  flow <- bio_h_summary[i,"flow.endpoints"]
  eng <- bio_h_summary[i,"eng.endpoints"]
  
  data <- AllDataLong2 %>%
    filter(Metric == bio,
           CurMetric == flow,
           Class2 == eng) %>%
    select(Curflow, MetricValue, COMID) %>% ## only metrics needed
    drop_na(Curflow, MetricValue) %>%
    filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
    distinct()
  
  
  ## new data - plus and minus 10% as a start
  flowvalues <- seq(range(data$Curflow)[1]-10,range(data$Curflow)[2]+10,0.05)
  
  ## get model, predict, extract all data and categories
  mod <- log.lm[[i]]
  predictedVals <- predict.glm(mod,list(Curflow = flowvalues),  type = "response")
  DFX <- as.data.frame(predictedVals)
  DFX$Value <- flowvalues
  DFX$Bio <- bio
  DFX$Variable <- flow
  DFX$Eng <- eng
  DFX$MinVal <-  range(data$Curflow)[1]
  DFX$MaxVal <-  range(data$Curflow)[2]
  
  DF <- bind_rows(DF, DFX)
  
}

# plot(DF$predictedVals, DF$BioVals)

### predicted figures
mets <- unique(DF$Eng)
## facet labels
supp.labs <- c(
  "Max_Wkl_Max_StreamT_grt_30_"="Weeks greater than 86F",
  "Max_Wkly_Mean_StreamT" = "Max Weekly Mean",
  "Max_Wkl_Max_StreamT" = "Weekly Maximum",
  "Min_Wkl_Min_StreamT" = "Weekly Minimum",
  "Max_Wkl_Rng_StreamT" = " Max Weekly Range",
  "Mean_Wkl_Rng_StreamT" =  "Av Weekly Range" 
)

for(m in 1:length(mets)) {
  
  T1 <- (ggplot(subset(DF, Eng == mets[m]), aes(y=predictedVals, x=Value, group = Eng, color = Eng)) +
           # geom_point(size=0.2) +
           geom_line( linewidth = 1)+
           # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
           # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
           facet_wrap(~Variable, labeller = as_labeller(supp.labs),
                      scales = "free") +
           # geom_vline(data=filter(DF, !Variable %in% c("Max_Wkl_Max_StreamT_grt_30_","Max_Wkl_Rng_StreamT", "Mean_Wkl_Rng_StreamT")), 
           #            aes(xintercept = 86), linetype="dashed", color = "red", linewidth=0.5) +
           geom_vline(data=filter(DF, !Variable %in% c("Max_Wkl_Max_StreamT_grt_30_","Max_Wkl_Rng_StreamT", "Mean_Wkl_Rng_StreamT")), 
                      aes(xintercept = 80), linetype="dashed", color = "blue", linewidth=0.5) +
           scale_x_continuous(name="Water flow (°F)") +
           scale_y_continuous(name =  "Probability of median CSCI")) +
    theme(legend.title=element_blank())
  
  T1
  
  file.name1 <- paste0(out.dir, "03_asci_", mets[m], "flow_response_predicted_glm_modChannels.jpg")
  ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
  
}

# Alteration models: CSCI -------------------------------------------------

AllDataLong2 <- AllData %>%
  filter(!Metric == "count") #%>%
# pivot_longer(AltflowMinus80:AltflowDivide80, names_to="AltMetric", values_to = "AltValues")

names(AllDataLong2)

## create df of model configurations

## bio 
biol.endpoints<-"csci"

## flow
flow.endpoints<- unique(na.omit(AllDataLong2$CurMetric))[2:4]
flow.endpoints

## Alt flow
alt.endpoints<- unique(na.omit(AllDataLong2$AltMetric))
alt.endpoints


bio_h_summary<-  expand.grid(biol.endpoints=biol.endpoints,flow.endpoints=flow.endpoints,  alt.endpoints= alt.endpoints, stringsAsFactors = F)

# bio_h_summary <- bio_h_summary[-c(1:9),]

## model of each configuration
log.glm <-lapply(1:nrow(bio_h_summary), function(i)
{
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  amet<-as.character(bio_h_summary[i,"alt.endpoints"])
  
  mydat<-AllDataLong2 %>%
    filter(Metric == bmet,
           CurMetric == tmet, 
           AltMetric == amet) %>%
    select(Altflow, MetricValue, COMID) %>% ## only metrics needed
    drop_na(Altflow, MetricValue) %>%
    filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
    distinct() #%>%
  # filter(Altflow > 0) ## take only positive altered flows
  
  mydat$Condition<-ifelse(mydat$MetricValue< 0.60 ,0,1) ## convert to binary
  # mydat <- mydat %>% drop_na(Condition)
  
  mydat<-mydat[order(mydat$MetricValue),] ## order by csci value
  glm(Condition~Altflow, family=binomial(link="logit"), data=mydat) ### glm
  
  # View(mydat)
})


## save models
save(log.glm, file = "output_data/03_csci_glm_altered.RData")

### get rsqds and pvals
for(i in 1:length(log.glm)) {
  
  mod <- summary(log.glm[[i]])
  bio_h_summary$AIC[i] <- mod$aic 
  bio_h_summary$PValue[i] <- mod$coefficients[8]
  bio_h_summary$McFaddensR2[i] <- 1-mod$deviance/mod$null.deviance
  # mod
}
## save configs and r sqds
save(bio_h_summary, file="output_data/03_csci_altered_glm_rsqds.RData")


## make df of predicted values to predict on - need to be different for each flow metric
## blank df
DF <- NULL
DF <- as.data.frame(DF)
bio_h_summary
### get predictions and fitted values
for(i in 1:length(log.glm)) {
  
  bio <- bio_h_summary[i,"biol.endpoints"]
  flow <- bio_h_summary[i,"flow.endpoints"]
  alt <- bio_h_summary[i,"alt.endpoints"]
  
  data <- AllDataLong2 %>%
    filter(Metric == bio,
           CurMetric == flow,
           AltMetric == alt) %>%
    select(Altflow, MetricValue, COMID) %>% ## only metrics needed
    drop_na(Altflow, MetricValue) %>%
    filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
    distinct() 
  
  ## new data - plus and minus 10 as a start
  
  if (alt == "AltflowDivide80" ) {
    flowvalues <- seq(range(data$Altflow)[1]*0.9,range(data$Altflow)[2]*1.1,0.05)
  } else {
    flowvalues <- seq(range(data$Altflow)[1]-10,range(data$Altflow)[2]+10,0.05)
  }
  
  
  ## get model, predict, extract all data and categories
  mod <- log.glm[[i]]
  predictedVals <- predict.glm(mod,list(Altflow = flowvalues),  type = "response")
  DFX <- as.data.frame(predictedVals)
  DFX$Value <- flowvalues
  DFX$Bio <- bio
  DFX$Variable <- flow
  DFX$AltVar <- alt
  DFX$MinVal <-  range(data$Altflow)[1]
  DFX$MaxVal <-  range(data$Altflow)[2]
  
  DF <- bind_rows(DF, DFX)
  
}


## remove mean weekly as not an expected relationship
## and scale probability?

DF <- DF %>%
  filter(!Variable == "Max_Wkly_Mean_StreamT")


### predicted figures

mets <- unique(DF$AltVar)

## facet labels
supp.labs <- c(
  "Max_Wkl_Max_StreamT_grt_30_"="Weeks greater than 86F",
  "Max_Wkly_Mean_StreamT" = "Max Weekly Mean",
  "Max_Wkl_Max_StreamT" = "Weekly Maximum",
  "Min_Wkl_Min_StreamT" = "Weekly Minimum",
  "Max_Wkl_Rng_StreamT" = " Max Weekly Range",
  "Mean_Wkl_Rng_StreamT" =  "Av Weekly Range" 
)

for(m in 1:length(mets)) {
  
  T1 <- (ggplot(subset(DF, AltVar == mets[m]), aes(y=predictedVals, x=Value, group = Variable, color = Variable)) +
           # geom_point(size=0.2) +
           geom_line( linewidth = 1)+
           # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
           # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
           facet_wrap(~Variable, labeller = as_labeller(supp.labs),
                      scales = "free") +
           # geom_vline(data=filter(DF, !Variable %in% c("Max_Wkl_Max_StreamT_grt_30_","Max_Wkl_Rng_StreamT", "Mean_Wkl_Rng_StreamT")), 
           #            aes(xintercept = 86), linetype="dashed", color = "red", linewidth=0.5, show.legend = T) +
           # geom_vline(data=filter(DF, !Variable %in% c("Max_Wkl_Max_StreamT_grt_30_","Max_Wkl_Rng_StreamT", "Mean_Wkl_Rng_StreamT")), 
           #            aes(xintercept = 80), linetype="dashed", color = "blue", linewidth=0.5, show.legend = T) +
           scale_x_continuous(name=paste0(mets[m],": flowerature Alteration (°F)")) +
           scale_y_continuous(name = "Probability of good CSCI") +
           theme(legend.position = "none"))
  
  T1
  
  file.name1 <- paste0(out.dir, "03_", mets[m], "_csci_flow_response_predicted_glm.jpg")
  ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
}



# Alteration Models: ASCI -------------------------------------------------

AllDataLong2 <- AllDataA %>%
  filter(!Metric == "count") #%>%
# pivot_longer(AltflowMinus80:AltflowDivide80, names_to="AltMetric", values_to = "AltValues")

## create df of model configurations

## bio 
biol.endpoints<-"ASCI_Hybrid"

## flow
flow.endpoints<- unique(na.omit(AllDataLong2$CurMetric))[2:4]
flow.endpoints

## Alt flow
alt.endpoints<- unique(na.omit(AllDataLong2$AltMetric))
alt.endpoints


bio_h_summary<-  expand.grid(biol.endpoints=biol.endpoints,flow.endpoints=flow.endpoints,  alt.endpoints= alt.endpoints, stringsAsFactors = F)

# bio_h_summary <- bio_h_summary[-c(1:9),]
bio_h_summary
i=1
i
## model of each configuration
log.glm <-lapply(1:nrow(bio_h_summary), function(i)
{
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  amet<-as.character(bio_h_summary[i,"alt.endpoints"])
  
  mydat<-AllDataLong2 %>%
    filter(Metric == bmet,
           CurMetric == tmet, 
           AltMetric == amet) %>%
    select(Altflow, MetricValue, COMID) %>% ## only metrics needed
    drop_na(Altflow, MetricValue) %>%
    filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
    distinct() #%>%
  # filter(Altflow > 0) ## take only positive altered flows
  
  mydat$Condition<-ifelse(mydat$MetricValue< 0.86 ,0,1) ## convert to binary
  # mydat <- mydat %>% drop_na(Condition)
  
  mydat<-mydat[order(mydat$MetricValue),] ## order by csci value
  glm(Condition~Altflow, family=binomial(link="logit"), data=mydat) ### glm
  
  # View(mydat)
})


## save models
save(log.glm, file = "output_data/03_asci_glm_altered.RData")

### get rsqds and pvals
for(i in 1:length(log.glm)) {
  
  mod <- summary(log.glm[[i]])
  bio_h_summary$AIC[i] <- mod$aic 
  bio_h_summary$PValue[i] <- mod$coefficients[8]
  bio_h_summary$McFaddensR2[i] <- 1-mod$deviance/mod$null.deviance
  # mod
}
## save configs and r sqds
save(bio_h_summary, file="output_data/03_asci_altered_glm_rsqds.RData")


## make df of predicted values to predict on - need to be different for each flow metric
## blank df
DF <- NULL
DF <- as.data.frame(DF)

### get predictions and fitted values
for(i in 1:length(log.glm)) {
  
  bio <- bio_h_summary[i,"biol.endpoints"]
  flow <- bio_h_summary[i,"flow.endpoints"]
  alt <- bio_h_summary[i,"alt.endpoints"]
  
  data <- AllDataLong2 %>%
    filter(Metric == bio,
           CurMetric == flow,
           AltMetric == alt) %>%
    select(Altflow, MetricValue, COMID) %>% ## only metrics needed
    drop_na(Altflow, MetricValue) %>%
    filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
    distinct() 
  
  ## new data - plus and minus 10% as a start
  if (alt == "AltflowDivide80" ) {
    flowvalues <- seq(range(data$Altflow)[1]*0.9,range(data$Altflow)[2]*1.1,0.05)
  } else {
    flowvalues <- seq(range(data$Altflow)[1]-10,range(data$Altflow)[2]+10,0.05)
  }
  ## get model, predict, extract all data and categories
  mod <- log.glm[[i]]
  predictedVals <- predict.glm(mod,list(Altflow = flowvalues),  type = "response")
  DFX <- as.data.frame(predictedVals)
  DFX$Value <- flowvalues
  DFX$Bio <- bio
  DFX$Variable <- flow
  DFX$AltVar <- alt
  DFX$MinVal <-  range(data$Altflow)[1]
  DFX$MaxVal <-  range(data$Altflow)[2]
  
  DF <- bind_rows(DF, DFX)
  
}


## remove mean weekly as not an expected relationship
## and scale probability?

DF <- DF %>%
  filter(!Variable == "Max_Wkly_Mean_StreamT")


### predicted figures

mets <- unique(DF$AltVar)

## facet labels
supp.labs <- c(
  "Max_Wkl_Max_StreamT_grt_30_"="Weeks greater than 86F",
  "Max_Wkly_Mean_StreamT" = "Max Weekly Mean",
  "Max_Wkl_Max_StreamT" = "Weekly Maximum",
  "Min_Wkl_Min_StreamT" = "Weekly Minimum",
  "Max_Wkl_Rng_StreamT" = " Max Weekly Range",
  "Mean_Wkl_Rng_StreamT" =  "Av Weekly Range" 
)

for(m in 1:length(mets)) {
  
  T1 <- (ggplot(subset(DF, AltVar == mets[m]), aes(y=predictedVals, x=Value, group = Variable, color = Variable)) +
           # geom_point(size=0.2) +
           geom_line( linewidth = 1)+
           # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
           # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
           facet_wrap(~Variable, labeller = as_labeller(supp.labs),
                      scales = "free") +
           # geom_vline(data=filter(DF, !Variable %in% c("Max_Wkl_Max_StreamT_grt_30_","Max_Wkl_Rng_StreamT", "Mean_Wkl_Rng_StreamT")), 
           #            aes(xintercept = 86), linetype="dashed", color = "red", linewidth=0.5, show.legend = T) +
           # geom_vline(data=filter(DF, !Variable %in% c("Max_Wkl_Max_StreamT_grt_30_","Max_Wkl_Rng_StreamT", "Mean_Wkl_Rng_StreamT")), 
           #            aes(xintercept = 80), linetype="dashed", color = "blue", linewidth=0.5, show.legend = T) +
           scale_x_continuous(name=paste0(mets[m],": flowerature Alteration (°F)")) +
           scale_y_continuous(name = "Probability of good ASCI") +
           theme(legend.position = "none"))
  
  T1
  
  file.name1 <- paste0(out.dir, "03_", mets[m], "_asci_flow_response_predicted_glm.jpg")
  ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
}



