### GAM model on all combined channel types

library(gam)
library(tidylog)
library(tidyverse)
library(sf)

## directory for figures
out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/SMC_Modified_Channels/Figures/"


# Upload & Format Data -------------------------------------------------------------

load(file = "output_data/00_bugs_algae_flow_joined_by_masterid.RData")
head(AllData)

## take only acsi h and csci
AllDataLong2 <- AllData %>%
  filter(Metric %in% c("csci", "asci"))

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
save(gam.lm, file = "output_data/01_csci_asci_all_gam_flow_combined_all.RData")

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
save(bio_h_summary, file="output_data/01_csci_asci_gam_rsqds_combined_all.RData")
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

  DFX <- cbind(na.omit(mod$data), as.data.frame(predictedVals)) %>%
    rename(Value = deltah_final) %>%
    mutate(Variable = tmet) 
  
  DF <- bind_rows(DF, DFX)
  
}
DF

## change back to numeric
DF <- DF %>%
  mutate(HydroValue = as.numeric(Value), BioValue = as.numeric(MetricValue),
         predictedVals = as.numeric(predictedVals)) %>%
  rename(hydro.endpoints = Variable) %>%
  select(-Value, -MetricValue)

str(DF)

### predicted figures
bio <- unique(DF$Metric)
bio
mets <- unique(DF$hydro.endpoints)
DF
m=1

for(b in 1:length(bio)) {
  
  DF1 <- DF %>%
    filter(Metric == bio[b])

  for(m in 1:length(mets)) {
    
    T1 <- (ggplot(subset(DF1, hydro.endpoints == mets[m]), aes(y=predictedVals, x=HydroValue)) +
             # geom_point(size=0.2) +
             geom_line( linewidth = 1)+
             # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
             # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
             # facet_wrap(~DirectionAlt, scales = "free") +
             scale_x_continuous(name=paste(mets[m])) +
             scale_y_continuous(name = paste0(bio[b], " Score"))) +
            theme_bw()
    # theme(legend.position = "none"))

    file.name1 <- paste0(out.dir, "01_", bio[b], "_", mets[m], "_flow_response_predicted_gam_combined_all.jpg")
    ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
    
  }
  
  
}


# format for thresholds ----------------------------------------------------------

## function to find value in curve
load("/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/Flow_curves_new/Code/functions/root_interpolation_function.Rdata")

## full names for labels
labels <- read.csv("/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/Flow_curves_new/Data/ffm_names.csv")
labels <- labels[1:24, ]
labels <- labels %>% rename(hydro.endpoints = Flow.Metric.Code)
labels[25, 1] <- "Peak Flow Magnitude (Q99, cfs)"
labels[25, 2] <- "Q99"
labels[25, 3] <- "Peak Flow Magnitude"
labels


## FIX NAMES TO MATCH LABELS AND LIMITS - Rachel 9/6
dfx <- DF %>% 
  mutate(hydro.endpoints = case_when(hydro.endpoints == "d_ds_mag_50" ~ "DS_Mag_50",            
                                     hydro.endpoints == "d_fa_mag" ~ "FA_Mag",
                                     hydro.endpoints == "d_peak_10" ~ "Peak_10",
                                     hydro.endpoints == "d_peak_2" ~ "Peak_2",
                                     hydro.endpoints == "d_peak_5" ~ "Peak_5",
                                     hydro.endpoints == "d_sp_mag" ~ "SP_Mag",
                                     hydro.endpoints == "d_wet_bfl_mag_10" ~ "Wet_BFL_Mag_10",
                                     hydro.endpoints == "d_wet_bfl_mag_50" ~ "Wet_BFL_Mag_50", 
                                     hydro.endpoints == "delta_q99" ~ "Q99")) %>%
  mutate(comb_code_type = paste(Metric, "_", hydro.endpoints, sep=""))
dfx

all_data <- left_join(dfx, labels, by ="hydro.endpoints")

head(all_data)


# calculate thresholds ----------------------------------------------------

### ASCI
## create df
df <- as.data.frame(matrix(ncol=13))
colnames(df) <- c("metric", "ThresholdNAT_Lower", "ThresholdNAT_Upper", "ThresholdHB_Lower", "ThresholdHB_Upper",
                  "ThresholdSB0_Lower", "ThresholdSB0_Upper", "ThresholdSB2_Lower", "ThresholdSB2_Upper",
                  "n", "Index", "Hydro_endpoint", "CorObs_Pred")
df

## define metrics
metrics <- unique(all_data$comb_code_type)

metrics

i = 1
## loop through metrics
for(i in 1: length(metrics)) {
  
  ## define metric
  met <- metrics[i]
  ## filter by metric
  hydroxx <- all_data %>%
    filter(comb_code_type == met)
  
  ## define thresholds - different for asci and csci
  if (hydroxx$Metric[1] == "asci") {
    
    ## get curves values at different probabilities
    threshNAT <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$predictedVals, 0.86) 
    # threshNAT <- ifelse(length(threshNAT) == 0, NA, threshNAT)

    threshHB <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$predictedVals, 0.87) 
    # threshHB <- ifelse(length(threshHB) == 0, NA, threshHB)
    
    # threshSB1 <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$predictedVals, 0.9) 
    # threshSB1 <- ifelse(length(threshSB1) == 0, NA, thresh90)
    
    threshSB0 <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$predictedVals, 0.79) 
    # threshSB0 <- ifelse(length(threshSB0) == 0, NA, threshSB0)
    
    threshSB2 <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$predictedVals, 0.76) 
    # threshSB2 <- ifelse(length(threshSB2) == 0, NA, threshSB2)
    
  } else {
    
    ## get curves values at different probabilities
    threshNAT <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$predictedVals, 0.79) 
    # threshNAT <- ifelse(length(threshNAT) == 0, NA, threshNAT)
    
    threshHB <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$predictedVals, 0.67) 
    # threshHB <- ifelse(length(threshHB) == 0, NA, threshHB)
    
    # threshSB1 <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$predictedVals, 0.9) 
    # threshSB1 <- ifelse(length(threshSB1) == 0, NA, thresh90)
    
    threshSB0 <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$predictedVals, 0.78) 
    # threshSB0 <- ifelse(length(threshSB0) == 0, NA, threshSB0)
    
    threshSB2 <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$predictedVals, 0.75) 
    # threshSB2 <- ifelse(length(threshSB2) == 0, NA, threshSB2)
    
  }


  ## add info to df
  df[i, 1] <- met
  df[i, 2] <- threshNAT[1]
  df[i, 3] <- threshNAT[2]
  df[i, 4] <- threshHB[1]
  df[i, 5] <- threshHB[2]
  df[i, 6] <- threshSB0[1]
  df[i, 7] <- threshSB0[2]
  df[i, 8] <- threshSB2[1]
  df[i, 9] <- threshSB2[2]
  df[i, 10] <- length(hydroxx$predictedVals)
  df[i ,11] <- hydroxx$Metric[1]
  df[i, 12] <- hydroxx$hydro.endpoints[1]
  df[i, 13] <- cor(hydroxx$predictedVals, hydroxx$BioValue)
  

  
}

hydroxx
df
write.csv(df, "output_data/01_delta_thresholds_GAMs.csv")


# format to apply data ----------------------------------------------------

# add best observed thresholds
## apply one ffm - raw data
## how many modified channels are within thresholds

head(df) ## thresholds

## format 
limits <- df %>%
  dplyr::select( -n) %>%
  mutate(metric = paste0(Index, "_", Hydro_endpoint)) %>%  # make code 
  pivot_longer(ThresholdNAT_Lower:ThresholdSB2_Upper, names_to = "Threshold", values_to = "DeltaH") %>% # make threshold names longer
  separate(Threshold, into=c("Threshold", "Type")) %>%
  mutate(Threshold = gsub("Threshold", "", Threshold)) %>%
  mutate(BioThresh = case_when((Index == "asci" & Threshold == "NAT") ~ 0.86, ## ad bio thresholds
                               (Index == "asci" & Threshold == "HB") ~ 0.87,
                               (Index == "asci" & Threshold == "SB0") ~ 0.79,
                               (Index == "asci" & Threshold == "SB2") ~ 0.76,
                               (Index == "csci" & Threshold == "NAT") ~ 0.79,
                               (Index == "csci" & Threshold == "HB") ~ 0.67,
                               (Index == "csci" & Threshold == "SB0") ~ 0.78,
                               (Index == "csci" & Threshold == "SB2") ~ 0.75))

limits

## make wider with type - lower/upper
limits <- limits %>%
  pivot_wider(names_from = Type, values_from = DeltaH)    

# joining labels table to limits - by hydroendpoints 
limits <- left_join(limits, labels, by = c("Hydro_endpoint" = "hydro.endpoints"))
limits

write.csv(limits, "output_data/01_deltaH_limits.csv")

## data
load(file = "output_data/00_bugs_algae_flow_joined_by_masterid.RData")
head(AllData)

names(AllData)

## take only acsi h and csci
AllDataLong2 <- AllData %>%
  select(-X) %>%
  filter(Metric %in% c("csci", "asci")) %>%
  drop_na(deltah_final) %>%
  rename(Hydro_endpoint = hydro.endpoints,
         Index = Metric,
         Threshold = channel_engineering_class,
         IndexValue = MetricValue)

## count ffm per class
tallyFFM <- AllData %>%
  group_by(channel_engineering_class) %>%
  select(masterid, COMID) %>%
  distinct() %>%
  drop_na(channel_engineering_class) %>%
  tally()

tallyFFM

## join limits to data
names(AllDataLong2)
names(limits)

allLims <- full_join(limits, AllDataLong2, by = c("Index", "Hydro_endpoint", "Threshold"))
head(allLims)
names(imps)

## define if hydro is within mod limits & bio above threshold
imps <- allLims %>%
  group_by(metric, Index, Hydro_endpoint, Threshold, BioThresh, masterid, COMID) %>%
  mutate(WithinHydroLimits = ifelse(deltah_final <= Upper & deltah_final >= Lower, "Within", "NotWithin")) %>%  ## within hydro limits
  mutate(WithinBioLimits = ifelse(IndexValue >= BioThresh, "Within", "NotWithin")) %>%## above bio thresholds
  mutate(Result = case_when((WithinHydroLimits == "Within" & WithinBioLimits == "Within") ~ "NoImpact", ## add impact
                            (WithinHydroLimits == "Within" & WithinBioLimits == "NotWithin") ~ "BioImpact", ##
                            (WithinHydroLimits == "NotWithin" & WithinBioLimits == "Within") ~ "HydroImpact",
                            (WithinHydroLimits == "NotWithin" & WithinBioLimits == "NotWithin") ~ "BothImpact")) %>%
  mutate(Result = factor(Result, levels = c("NoImpact", "BioImpact", "HydroImpact", "BothImpact")))


  ## tally of imopact per ffm

tallyImpact <- imps %>%
  group_by(Index, Hydro_endpoint, Flow.Metric.Name,Threshold, Result) %>%
  distinct() %>%
  tally() %>%
  drop_na(Result) %>%
  mutate(PercChans = (n/sum(n)*100))
  
tallyImpact      

### make table for a couple of ffms

## csci and ds mag 50
cscids
cscids <- tallyImpact %>%
  filter(Index == "csci", Hydro_endpoint == "DS_Mag_50") %>%
  select(-n) %>%
  pivot_wider(names_from = Result, values_from = PercChans) %>%
  select(Index, Flow.Metric.Name, Threshold, NoImpact, BioImpact:BothImpact)

## save out
write.csv(cscids, "01_csci_ds_mag_50.csv")

## csci & wet season
csciws <- tallyImpact %>%
  filter(Index == "csci", Hydro_endpoint == "Wet_BFL_Mag_50") %>%
  select(-n) %>%
  pivot_wider(names_from = Result, values_from = PercChans) %>%
  select(Index, Flow.Metric.Name, Threshold, NoImpact, BioImpact, HydroImpact, BothImpact)

## save out
write.csv(csciws, "01_csci_ws_bfl_mag_50.csv")
         