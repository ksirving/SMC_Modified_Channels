## Gams with Gam package

library(tidylog)
library(tidyverse)
library(sf)
library(gam)
library(mapview)

out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/SMC_Modified_Channels/Figures/"

getwd()

setwd("/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/SMC_Modified_Channels")

# Upload & Format Data -------------------------------------------------------------

load(file = "output_data/00_bugs_algae_flow_joined_by_masterid.RData")
head(AllData)

# dfx <- AllData %>%
#   select(masterid, COMID, longitude, latitude) %>%
#   distinct()
# 
# write.csv(dfx, "/Users/katieirving/Library/CloudStorage/OneDrive-SharedLibraries-SCCWRP/Temp Modelling - General/Data/Raw/LAR/SoCal_Bio_sites.csv")

## how many bio sites 
length(unique(AllData$masterid))

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

## smoothing functions
smooth_funcs <- c(3,6,9,12,15,18,20)

### make grid
bio_h_summary<-  expand.grid(biol.endpoints=biol.endpoints,flow.endpoints=flow.endpoints, smooth_funcs = smooth_funcs, stringsAsFactors = F)

bio_h_summary

## create df of model configurations

head(AllDataLong2)
i=1
## model of each configuration
gam.lm <-lapply(1:nrow(bio_h_summary), function(i)
{
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  smet<-as.numeric(bio_h_summary[i,"smooth_funcs"])

  
  mydat<-AllDataLong2 %>%
    filter(Metric == bmet,
           flow_metric == tmet) %>%
    select(Metric, MetricValue, deltah_final, channel_engineering_class) %>% ## only metrics needed
    drop_na( deltah_final) %>%
    # filter(deltah_final > 0) %>%
    filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
    distinct()
  
  
  mydat<-mydat[order(mydat$MetricValue),] ## order by csci value
  # try(glm(Condition~deltah_final, family=binomial(link="logit"), data=mydat), silent=T) ### glm
  
  gam(MetricValue~s(deltah_final, smet), family = Gamma(link = "log"), data = mydat)
  
  
})

## save models
save(gam.lm, file = "output_data/03_csci_asci_all_gam_flow_combined_all_Gam_pack.RData")

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
save(bio_h_summary, file="output_data/03_csci_asci_gam_rsqds_combined_all.RData")
bio_h_summary

## plot k vals with r2 

bio_h_summary

P1 <- ggplot(data = bio_h_summary, aes(y=McFaddensR2, x=smooth_funcs, group = flow.endpoints, colour = flow.endpoints)) +
         # geom_point(size=0.2) +
         geom_line( linewidth = 1)+
         # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
         # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
         facet_wrap(~biol.endpoints, scales = "free") +
         scale_x_continuous(name= "Smoothing Function (k)") +
         scale_y_continuous(name = "McFaddens Rq") 
# theme(legend.position = "none"))

P1

file.name1 <- paste0(out.dir, "03_smoothing_ks_rsq.jpg")
ggsave(P1, filename=file.name1, dpi=300, height=5, width=7.5)

## plot k vals with AIC

P2 <- ggplot(data = bio_h_summary, aes(y=AIC, x=smooth_funcs, group = flow.endpoints, colour = flow.endpoints)) +
  # geom_point(size=0.2) +
  geom_line( linewidth = 1)+
  # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
  # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
  facet_wrap(~biol.endpoints, scales = "free") +
  scale_x_continuous(name= "Smoothing Function (k)") +
  scale_y_continuous(name = "AIC") 
# theme(legend.position = "none"))

P2

file.name1 <- paste0(out.dir, "03_smoothing_ks_AIC.jpg")
ggsave(P2, filename=file.name1, dpi=300, height=5, width=7.5)


## make df of predicted values to predict on - need to be different for each flow metric

## blank df
DF <- NULL
DF <- as.data.frame(DF)

### get predictions and fitted values
for(i in 1:length(gam.lm)) {
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  smet<-as.numeric(bio_h_summary[i,"smooth_funcs"])
  
  ## get model, predict, extract all data and categories
  mod <- gam.lm[[i]]
  predictedVals <- predict(mod,  type = "response")
  
  ## data from model, remove NAs
  modData <- mod$data %>%
    drop_na(deltah_final, MetricValue)
  modData
  ## combine with predictions
  DFX <- cbind(modData, as.data.frame(predictedVals)) %>%
    # rename(Value = deltah_final) %>%
    mutate(Variable = tmet, SmoothingFunc = smet) 
  
  DF <- bind_rows(DF, DFX)
  
}
DF

## change back to numeric
DF <- DF %>%
  mutate(deltah_final = as.numeric(deltah_final), MetricValue = as.numeric(MetricValue),
         predictedVals = as.numeric(predictedVals))
str(DF)

### predicted figures
bio <- unique(DF$Metric)
mets <- unique(DF$Variable)



for(b in 1:length(bio)) {
  
  DF1 <- DF %>%
    filter(Metric == bio[b])
  
  head(DF1)
  
  for(m in 1:length(mets)) {
    
    T1 <- (ggplot(subset(DF1, Variable == mets[m]), aes(y=predictedVals, x=deltah_final)) +
             geom_smooth( linewidth = 1, aes())+
             geom_point( aes(x=deltah_final, y = MetricValue, colour = channel_engineering_class)) +
             # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
             # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
             geom_vline(xintercept =0) +
             # facet_wrap(~SmoothingFunc, scales = "free") +
             scale_x_continuous(name=paste(mets[m])) +
             scale_y_continuous(name = paste0("Probability of Good ", bio[b]))) 
    # theme(legend.position = "none"))
    
    T1
    
    file.name1 <- paste0(out.dir, "03_", bio[b], "_", mets[m], "_flow_response_predicted_gam_combined_all_smoothingKs.jpg")
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

DF
## FIX NAMES TO MATCH LABELS AND LIMITS 
dfx <- DF %>% 
  mutate(hydro.endpoints = case_when(Variable == "d_ds_mag_50" ~ "DS_Mag_50",            
                                     Variable == "d_fa_mag" ~ "FA_Mag",
                                     Variable == "d_peak_10" ~ "Peak_10",
                                     Variable == "d_peak_2" ~ "Peak_2",
                                     Variable == "d_peak_5" ~ "Peak_5",
                                     Variable == "d_sp_mag" ~ "SP_Mag",
                                     Variable == "d_wet_bfl_mag_10" ~ "Wet_BFL_Mag_10",
                                     Variable == "d_wet_bfl_mag_50" ~ "Wet_BFL_Mag_50", 
                                     Variable == "delta_q99" ~ "Q99")) %>%
  mutate(comb_code_type = paste(Metric, "_", hydro.endpoints, "_", SmoothingFunc, sep="")) %>%
  rename(HydroValue = deltah_final, BioValue = MetricValue)
dfx

all_data <- left_join(dfx, labels, by ="hydro.endpoints")

head(all_data)


# calculate thresholds ----------------------------------------------------

## create df
df <- as.data.frame(matrix(ncol=20))
colnames(df) <- c("metric", "ThresholdNAT_Lower", "ThresholdNAT_Upper", "ThresholdNATMed_Lower", "ThresholdNATMed_Upper",
                  "ThresholdNATLow_Lower", "ThresholdNATLow_Upper",
                  "ThresholdNATHigh_Lower", "ThresholdNATHigh_Upper", "ThresholdHB_Lower", "ThresholdHB_Upper",
                  "ThresholdSB0_Lower", "ThresholdSB0_Upper", "ThresholdSB2_Lower", "ThresholdSB2_Upper",
                  "n", "Index", "Hydro_endpoint","SmoothingFunc", "CorObs_Pred")
df

## define metrics
metrics <- unique(all_data$comb_code_type)
metrics
## loop through metrics
for(i in 1: length(metrics)) {
  
  ## define metric
  met <- metrics[i]
  
  ## filter by metric
  hydroxx <- all_data %>%
    filter(comb_code_type == met)
  
  hydroxx
  
  ## define thresholds - different for asci and csci
  if (hydroxx$Metric[1] == "asci") {
    
    ## get curves values at different probabilities
    threshNAT <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$predictedVals, 0.86) 
    # threshNAT <- ifelse(length(threshNAT) == 0, NA, threshNAT)
    threshNATMed <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$predictedVals, 0.86)## sep ref threshold for overall tally
    
    threshNATLow <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$predictedVals, 0.75) 
    
    threshNATHigh <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$predictedVals, 0.94) 
    
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
    
    threshNATMed <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$predictedVals, 0.79) ## sep ref threshold for overall tally
    
    threshNATLow <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$predictedVals, 0.63) 
    
    threshNATHigh <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$predictedVals, 0.92) 
    
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
  df[i, 4] <- threshNATMed[1]
  df[i, 5] <- threshNATMed[2]
  df[i, 6] <- threshNATLow[1]
  df[i, 7] <- threshNATLow[2]
  df[i, 8] <- threshNATHigh[1]
  df[i, 9] <- threshNATHigh[2]
  df[i, 10] <- threshHB[1]
  df[i, 11] <- threshHB[2]
  df[i, 12] <- threshSB0[1]
  df[i, 13] <- threshSB0[2]
  df[i, 14] <- threshSB2[1]
  df[i, 15] <- threshSB2[2]
  df[i, 16] <- length(hydroxx$predictedVals)
  df[i ,17] <- hydroxx$Metric[1]
  df[i, 18] <- hydroxx$hydro.endpoints[1]
  df[i, 19] <- hydroxx$SmoothingFunc[1]
  df[i, 20] <- cor(hydroxx$predictedVals, hydroxx$BioValue)
  
  
  
}

hydroxx
df ## NAs - where curve doesn't reach threshold
write.csv(df, "output_data/03_delta_thresholds_GAMs.csv")


# format to apply data to categories----------------------------------------------------
df <- read.csv("output_data/03_delta_thresholds_GAMs.csv")
# add best observed thresholds
## apply one ffm - raw data
## how many modified channels are within thresholds

## number of sites in each ffm
## csci = 304
## asci = 171

head(df) ## thresholds

## format 
limits <- df %>%
  dplyr::select( -n) %>%
  mutate(metric = paste0(Index, "_", Hydro_endpoint)) %>%  # make code 
  pivot_longer(ThresholdNAT_Lower:ThresholdSB2_Upper, names_to = "Threshold", values_to = "DeltaH") %>% # make threshold names longer
  separate(Threshold, into=c("Threshold", "Type")) %>%
  mutate(Threshold = gsub("Threshold", "", Threshold)) %>%
  mutate(BioThresh = case_when((Index == "asci" & Threshold == "NAT") ~ 0.86, ## ad bio thresholds
                               (Index == "asci" & Threshold == "NATMed") ~ 0.86,
                               (Index == "asci" & Threshold == "NATLow") ~ 0.75,
                               (Index == "asci" & Threshold == "NATHigh") ~ 0.94,
                               (Index == "asci" & Threshold == "HB") ~ 0.87,
                               (Index == "asci" & Threshold == "SB0") ~ 0.79,
                               (Index == "asci" & Threshold == "SB2") ~ 0.76,
                               (Index == "csci" & Threshold == "NAT") ~ 0.79,
                               (Index == "csci" & Threshold == "NATMed") ~ 0.79,
                               (Index == "csci" & Threshold == "NATLow") ~ 0.63,
                               (Index == "csci" & Threshold == "NATHigh") ~ 0.92,
                               (Index == "csci" & Threshold == "HB") ~ 0.67,
                               (Index == "csci" & Threshold == "SB0") ~ 0.78,
                               (Index == "csci" & Threshold == "SB2") ~ 0.75))

limits

## make wider with type - lower/upper
limits <- limits %>%
  pivot_wider(names_from = Type, values_from = DeltaH)    
limits
# joining labels table to limits - by hydroendpoints 
# limits <- left_join(limits, labels, by = c("Hydro_endpoint" = "hydro.endpoints"))
# limits
# 
# write.csv(limits, "output_data/01_deltaH_limits_mod_classes.csv")

## data
load(file = "output_data/00_bugs_algae_flow_joined_by_masterid.RData")
head(AllData)

names(AllData)

## take only acsi h and csci and rename to match columns
AllDataLong2 <- AllData %>%
  select(-X, -Class2) %>%
  filter(Metric %in% c("csci", "asci")) %>%
  drop_na(deltah_final) %>%
  rename(Hydro_endpoint = hydro.endpoints,
         Index = Metric,
         Threshold = channel_engineering_class,
         IndexValue = MetricValue)

AllDataLong2

### sep df for overall thresholds - nat low/high

AllDataLongOverall <- AllDataLong2 %>%
  mutate(ThresholdNatLow = "NATLow", ThresholdNatMed = "NATMed", ThresholdNatHigh = "NATHigh") %>%
  select(-c(Threshold)) %>%
  pivot_longer(ThresholdNatLow:ThresholdNatHigh, names_to = "Check", values_to = "Threshold") %>%
  select(-Check) 

## join categorised and overall dfs together

AllDF <- bind_rows(AllDataLongOverall, AllDataLong2) %>% distinct()
class(AllDF)
names(AllDF)

## count site per class with FFM
tallyFFM <- AllDF %>%
  group_by(Threshold) %>%
  select(masterid, COMID) %>%
  distinct() %>%
  drop_na(Threshold) %>%
  tally()

tallyFFM

## join limits to data
names(AllDF)
names(limits)

allLims <- full_join(limits, AllDF, by = c("Index", "Hydro_endpoint", "Threshold"), relationship = "many-to-many")
names(allLims)
names(imps)

## define if hydro is within mod limits & bio above threshold
imps <- allLims %>%
  group_by(metric, Index, Hydro_endpoint, Threshold, BioThresh, masterid, COMID, SmoothingFunc) %>%
  mutate(WithinHydroLimits = ifelse(deltah_final <= Upper & deltah_final >= Lower, "Within", "NotWithin")) %>%  ## within hydro limits
  mutate(WithinBioLimits = ifelse(IndexValue >= BioThresh, "Within", "NotWithin")) %>%## above bio thresholds
  mutate(Result = case_when((WithinHydroLimits == "Within" & WithinBioLimits == "Within") ~ "NoImpact", ## add impact
                            (WithinHydroLimits == "Within" & WithinBioLimits == "NotWithin") ~ "BioImpact", ##
                            (WithinHydroLimits == "NotWithin" & WithinBioLimits == "Within") ~ "HydroImpact",
                            (WithinHydroLimits == "NotWithin" & WithinBioLimits == "NotWithin") ~ "BothImpact")) %>%
  mutate(Result = factor(Result, levels = c("NoImpact", "BioImpact", "HydroImpact", "BothImpact"))) # %>%
# select(-c(Flow.Metric.Name.y, Flow.Component.y)) %>%
# rename(Flow.Metric.Name = Flow.Metric.Name.x, Flow.Component = Flow.Component.x)

write.csv(imps, "output_data/03_impact_ffm_bio.csv")

unique(imps$Threshold)
## tally of impact per ffm

tallyImpact <- imps %>%
  group_by(Index, Hydro_endpoint, Flow.Metric.Name,Threshold, Result, SmoothingFunc) %>%
  distinct() %>%
  tally() %>%
  drop_na(Result) %>%
  mutate(PercChans = (n/sum(n)*100))

tallyImpact

write.csv(tallyImpact, "output_data/03_count_impact.csv")

tallyImpact <- read.csv("output_data/03_count_impact.csv")

# Make tables for slides --------------------------------------------------

## use ggmap to get google 
library(ggmap)
library("ggsci")

# tallyImpact <- read.csv("output_data/03_percent_impacts_each_Class.csv")
names(tallyImpact1)

### make table for a couple of ffms

tallyImpact1 <- tallyImpact %>% ## using option with least NAs
  ungroup() %>%
  select(Index, Flow.Metric.Name, Threshold, Result, SmoothingFunc, PercChans, n) %>%
  rename(FlowMetric = Flow.Metric.Name, ModifiedClass = Threshold, Impact = Result, NumberOfSites = n, PercentageOfSites = PercChans) %>%
  mutate(ModifiedClass = case_when((ModifiedClass == "HB") ~ "Hard Bottom",
                                   (ModifiedClass == "NAT") ~ "Natural",
                                   (ModifiedClass == "NATHigh") ~ "Overall 1st",
                                   (ModifiedClass == "NATLow") ~ "Overall 30th",
                                   (ModifiedClass == "NATMed") ~ "Overall 10th",
                                   (ModifiedClass == "SB0") ~ "Soft Bottom (no hard sides)",
                                   (ModifiedClass == "SB2") ~ "Soft Bottom (two hard sides)"))

## csci and ds mag 50
cscids
cscids <- tallyImpact1 %>%
  filter(Index == "csci", FlowMetric == "Dry-season median baseflow", SmoothingFunc == 6) %>%
  select( -NumberOfSites) %>%
  pivot_wider(names_from = Impact, values_from = PercentageOfSites) %>%
  select(Index, FlowMetric, ModifiedClass, SmoothingFunc, NoImpact, BioImpact:BothImpact) %>%
  filter(ModifiedClass %in% c("Natural", "Hard Bottom", "Soft Bottom (no hard sides)", "Soft Bottom (two hard sides)")) %>%
  mutate(ModifiedClass = factor(ModifiedClass, 
                                levels = c("Natural", "Hard Bottom", "Soft Bottom (no hard sides)", "Soft Bottom (two hard sides)"),
                                labels = c("Natural", "Hard Bottom", "Soft Bottom (0)", "Soft Bottom (2)")))

## save out
# write.csv(cscids, "01_csci_ds_mag_50.csv")

## wet season

csciws <- tallyImpact1 %>%
  filter(Index == "csci", FlowMetric == "Wet-season median baseflow", SmoothingFunc == 6) %>%
  select( -NumberOfSites) %>%
  # pivot_wider(names_from = Result, values_from = PercChans) %>%
  pivot_wider(names_from = Impact, values_from = PercentageOfSites) %>%
  select(Index, FlowMetric, ModifiedClass, SmoothingFunc, NoImpact, BioImpact:BothImpact) %>%
  filter(ModifiedClass %in% c("Natural", "Hard Bottom", "Soft Bottom (no hard sides)", "Soft Bottom (two hard sides)")) %>%
  mutate(ModifiedClass = factor(ModifiedClass, 
                                levels = c("Natural", "Hard Bottom", "Soft Bottom (no hard sides)", "Soft Bottom (two hard sides)"),
                                labels = c("Natural", "Hard Bottom", "Soft Bottom (0)", "Soft Bottom (2)")))
csciws
cscids

## dry season
ds <- ggplot(cscids) +
  geom_col(mapping = aes(y=BothImpact, x=ModifiedClass, col = ModifiedClass, fill = ModifiedClass)) +
  scale_color_jco(name = "Modified Class") +
  labs(x="Modified Channel Type", y= "Percent Channels") +
  scale_y_continuous(limits = c(0,100)) +
  theme(axis.title = element_text(size = 15,
                                  face = "bold")) +
  theme(axis.text = element_text(size = 12)) +
  theme(legend.position = "none")
ds
file.name1 <- paste0(out.dir, "03_drySeason_column_perc_bothImpact_allTypesIncNat.jpg")
ggsave(ds, filename=file.name1, dpi=600, height=7, width=10)

## wet season
ws <- ggplot(csciws) +
  geom_col(mapping = aes(y=BothImpact, x=ModifiedClass, col = ModifiedClass, fill = ModifiedClass)) +
  scale_color_jco(name = "Modified Class") +
  labs(x="Modified Channel Type", y= "Percent Channels") +
  scale_y_continuous(limits = c(0,100)) +
  theme(axis.title = element_text(size = 15,
                                  face = "bold")) +
  theme(axis.text = element_text(size = 12)) +
  theme(legend.position = "none")
ws
file.name1 <- paste0(out.dir, "03_wetSeason_column_perc_bothImpact_allTypesIncNat.jpg")
ggsave(ws, filename=file.name1, dpi=600, height=7, width=10)

### stacked columns

## join ds and ws together

head(tallyImpact)

tallyImpactx <- tallyImpact %>%
  select(-NumberSitesPerClass, -X) %>%
  filter(Index == "csci") %>%
  filter(ModifiedClass %in% c("Natural", "Hard Bottom", "Soft Bottom (no hard sides)", "Soft Bottom (two hard sides)")) %>%
  mutate(ModifiedClass = factor(ModifiedClass, 
                                levels = c("Natural", "Hard Bottom", "Soft Bottom (no hard sides)", "Soft Bottom (two hard sides)"),
                                labels = c("Natural", "Hard Bottom", "Soft Bottom (0)", "Soft Bottom (2)"))) %>%
  pivot_longer(NoImpact:BothImpact, names_to = "Result", values_to = "PercChans") %>%
  mutate(Result = factor(Result, levels = c("NoImpact", "BioImpact", "HydroImpact", "BothImpact"), 
                         labels = c("Bio & Flow Within Range", "Low Bio", "Altered Flow", "Low Bio & Altered Flow")))
tallyImpactx


a1 <- ggplot(tallyImpactx, aes(fill=Result, y=PercChans, x=ModifiedClass)) + 
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~FlowMetric) +
  scale_fill_jco(name = "Category") +
  scale_y_continuous(name = "Sites (%)")

a1

file.name1 <- paste0(out.dir, "03_all_hydro_stacked_perc.jpg")
ggsave(a1, filename=file.name1, dpi=600, height=7, width=10)


tallyImpactxLong
## map

basemapsList <- c("Esri.WorldTopoMap", "Esri.WorldImagery",
                  "Esri.NatGeoWorldMap",
                  "OpenTopoMap", "OpenStreetMap", 
                  "CartoDB.Positron", "Stamen.TopOSMFeatures")

mapviewOptions(basemaps=basemapsList, fgb = FALSE)

## filter dry season and make spatial
impsx <- imps %>%
  filter(Hydro_endpoint == "DS_Mag_50", Index == "csci",
         Result == "BothImpact") %>%
  st_as_sf(coords=c("longitude", "latitude"), crs=4326, remove=F)
names(impsx)

m1 <- mapview(impsx, zcol = "Threshold",  col.regions=c("red", "green", "orange", "blue"),
              layer.name="Channel Type")


m1@map %>% leaflet::addMeasure(primaryLengthUnit = "meters")
# 
mapshot(m1, url = paste0(getwd(), "/output_data/01_dry_season_bfl_channel_types_flow_problem.html"),
        file = paste0(getwd(), "/ignore/01_dry_season_bfl_channel_types_flow_problem.png"))

