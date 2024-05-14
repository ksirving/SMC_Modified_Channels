## get GAM thresholds

library(gam)
library(tidyverse)
library(sf)
library(tidylog)
## load DF 

load(file = "ignore/04_quantGams_smooths_predictions_noPeakAugs.RData")
head(DF)

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
  mutate(comb_code_type = paste0(Metric, "_", hydro.endpoints, "_", Smooths, "_", Quant)) %>%
  rename(HydroValue = deltah_final, BioValue = predictedVals)

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
                  "n", "Index", "Hydro_endpoint","SmoothingFunc", "Quantile")
df

## define metrics
metrics <- unique(all_data$comb_code_type)
metrics

i=74
## loop through metrics
for(i in 1: length(metrics)) {
  
  ## define metric
  met <- metrics[i]
  
  ## filter by metric
  hydroxx <- all_data %>%
    filter(comb_code_type == met)
  
  head(hydroxx)
  
  ## define thresholds - different for asci and csci
  if (hydroxx$Metric[1] == "asci") {
    
    ## get curves values at different probabilities
    threshNAT <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.86) 
    # threshNAT <- ifelse(length(threshNAT) == 0, NA, threshNAT)
  
    
    threshNATMed <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.86)## sep ref threshold for overall tally
    
    threshNATLow <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.75) 
    
    threshNATHigh <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.94) 
    
    threshHB <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.87) 
    # threshHB <- ifelse(length(threshHB) == 0, NA, threshHB)
    
    # threshSB1 <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.9) 
    # threshSB1 <- ifelse(length(threshSB1) == 0, NA, thresh90)
    
    threshSB0 <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.79) 
    # threshSB0 <- ifelse(length(threshSB0) == 0, NA, threshSB0)
    
    threshSB2 <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.76) 
    # threshSB2 <- ifelse(length(threshSB2) == 0, NA, threshSB2)
    
  } else {
    
    ## get curves values at different probabilities
    threshNAT <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.79) 
    threshNAT
    # threshNAT <- ifelse(length(threshNAT) == 0, NA, threshNAT)
    
    threshNATMed <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.79) ## sep ref threshold for overall tally
    
    threshNATLow <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.63) 
    
    threshNATHigh <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.92) 
    
    threshHB <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.67) 
    # threshHB <- ifelse(length(threshHB) == 0, NA, threshHB)
    
    # threshSB1 <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.9) 
    # threshSB1 <- ifelse(length(threshSB1) == 0, NA, thresh90)
    
    threshSB0 <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.78) 
    # threshSB0 <- ifelse(length(threshSB0) == 0, NA, threshSB0)
    threshSB0
    
    threshSB2 <- RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.75) 
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
  df[i, 16] <- length(hydroxx$BioValue)
  df[i ,17] <- hydroxx$Metric[1]
  df[i, 18] <- hydroxx$hydro.endpoints[1]
  df[i, 19] <- hydroxx$Smooths[1]
  df[i, 20] <- paste(hydroxx$Quant[1])

  
  
}

hydroxx
df ## NAs - where curve doesn't reach threshold
write.csv(df, "output_data/05_delta_thresholds_GAMs.csv")

df <- read.csv("output_data/05_delta_thresholds_GAMs.csv")

head(df)

### metrics and smooths

# csci
## DS_Mag_50 - 6
## FA_Mag - 6

sm6csci <- c("DS_Mag_50", "FA_Mag")
sm6csci

## Peak_2 - 3
## Peak_5 - 3
## Peak_10 - 3
## Q99 - 3
## Sp_mag - 3
## wet_bfl_mag_10 - 3
## wet_bfl_mag_10 - 3

sm3csci <- c("Peak_10", "Peak_2", "Peak_5", "SP_Mag", "Wet_BFL_Mag_10", "Wet_BFL_Mag_50", "Q99")
sm3csci
## asci
## DS_Mag_50 - 6
## FA_Mag - 6
## Sp_mag - 6

sm6asci <- c("DS_Mag_50", "FA_Mag", "SP_Mag")
sm6asci

## Peak_2 - 3
## Peak_5 - 3
## Peak_10 - 3
## Q99 - 3
## wet_bfl_mag_10 - 3
## wet_bfl_mag_10 - 3

sm3asci <- c("Peak_10", "Peak_2", "Peak_5", "Wet_BFL_Mag_10", "Wet_BFL_Mag_50", "Q99")
sm3asci

unique(df$Hydro_endpoint)

## subset models/thresholds to only the ones we're using
##  quants 0.9 & 0.5 for all metrics
dfsub <- df %>%
  # filter(Quantile == c(0.9, 0.5)) %>%
  mutate(ChosenSmth = case_when((Index == "csci" & Hydro_endpoint %in% sm6csci) ~ 6,
                               (Index == "csci" & Hydro_endpoint %in% sm3csci) ~ 3,
                               (Index == "asci" & Hydro_endpoint %in% sm6asci) ~ 6,
                               (Index == "asci" & Hydro_endpoint %in% sm3asci) ~ 3)) %>%
  mutate(ModelstoUse = ifelse(SmoothingFunc == ChosenSmth, "Yes", "No")) %>%
  filter(ModelstoUse == "Yes") %>%
  filter(Quantile %in% c(0.90, 0.50))
  
write.csv(dfsub, "output_data/05_chosen_models_thresholds.csv")

## only these metrics

# "DS_Mag_50", "Peak_10", "FA_Mag", "Wet_BFL_Mag_50"

# format to apply data to categories----------------------------------------------------

df <- read.csv("output_data/05_chosen_models_thresholds.csv")
# add best observed thresholds
## apply one ffm - raw data
## how many modified channels are within thresholds


## number of sites in each ffm
## csci = 304
## asci = 171

head(df) ## thresholds

## format 
limits <- df %>%
  dplyr::select( -n, -X.1, -X) %>%
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
  pivot_wider(names_from = Type, values_from = DeltaH)  %>%
  mutate(MetricCurve = paste0(metric,"_", SmoothingFunc, "_", Quantile))

unique(limits$MetricCurve)

## organise limits by type of curve: AboveCurve, BelowCurve

AboveCurve <- c("csci_SP_Mag_3_0.9", "csci_Wet_BFL_Mag_50_3_0.9")
BelowCurve <- c("asci_DS_Mag_50_6_0.5", "asci_FA_Mag_6_0.5", "asci_Peak_10_3_0.5", "asci_Peak_2_3_0.5", 
                "asci_Peak_5_3_0.5", "asci_Wet_BFL_Mag_10_3_0.5", "asci_Wet_BFL_Mag_50_3_0.5", "asci_Q99_3_0.5")

## assign labels to know wheich curves are above or below thresholds
## some thresholds hit curve, so make sure it's NA too
limits <- limits %>%
  mutate(CurveStatus = ifelse(is.na(Upper) & MetricCurve %in% AboveCurve, "Above", "Passes")) %>%
  mutate(CurveStatus = ifelse(is.na(Lower) & MetricCurve %in% BelowCurve, "Below", CurveStatus)) 

## format thresholds to follow interpretation
## is there are 2 negative values, take one closest to zero
## make sure it's the lower threshold
## a negative and a positive is fine
## a positive in the lower column should be an upper value

limitsx <- limits %>%
  mutate(Lower2 = ifelse(Upper < 0, Upper, Lower)) %>% ## 2 negs, take closest to 0
  mutate(Upper2 = ifelse(Lower > 0, Lower, Upper)) %>% ## change +ves in lower value to upper
  mutate(Lower = ifelse(Lower > 0, NA, Lower)) %>% ## change +ve ones to NA in lower column
  mutate(Lower2 = ifelse(is.na(Lower2), Lower, Lower2)) ## switch over remaining lower values


## add in minimum delta H values, where no limits can be found
## here: alteration isn't an issue as curve always within limits

MinDel <- min(na.omit(AllData$deltah_final)) ## -44933.65
MaxDel <- max(na.omit(AllData$deltah_final)) ## 7536.166

## remove negative values from upper - residual from swapping closest to zeros above
## also add min max for curves above thresholds

limitsx1 <- limitsx %>%
  mutate(Upper2 = ifelse(Upper2 < 0, NA, Upper2)) %>% ## remove -ve vals from upper
  mutate(Lower2 = ifelse(CurveStatus == "Above", MinDel, Lower2),
         Upper2 = ifelse(CurveStatus == "Above", MaxDel, Upper2)) %>% ## all curves above add min and max values
  mutate(MetricCurve2 = paste0(MetricCurve, "_", Threshold)) 

## next issue: some curves don't intersect with some thresholds
## extract these instances to take a look

cheCur <- limitsx1 %>%
  filter(is.na(Lower2) & is.na(Upper2) & CurveStatus == "Passes")
head(cheCur)
## all 0.5s should be below, and 0.9s above

cheCur2 <- cheCur %>%
  # select(MetricCurve2, Quantile, Threshold, BioThresh, CurveStatus) %>% ## take only some columns
  mutate(CurveStatus = ifelse(Quantile == 0.5, "Below", "Above")) %>% 
  mutate(Lower2 = ifelse(CurveStatus == "Above", MinDel, Lower2),
         Upper2 = ifelse(CurveStatus == "Above", MaxDel, Upper2)) ## all curves above add min and max values

## define these curves
ToRemove <- unique(cheCur2$MetricCurve2)

## remove from main DF
limitsx2 <- limitsx1 %>%
  filter(!MetricCurve2 %in% ToRemove)

## add others back in
## all curve below threshold change to 0, should result in "outside limits"
limitsx3 <- bind_rows(cheCur2, limitsx2) %>%
  mutate(Lower2 = ifelse(CurveStatus == "Below", 0, Lower2),
         Upper2 = ifelse(CurveStatus == "Below", 0, Upper2)) 

## finally: if a curve has only 1 value, add min or max to missing values

limitsx3 <- limitsx3 %>%
  mutate(Lower2 = ifelse(is.na(Lower2), MinDel, Lower2)) %>%
  mutate(Upper2 = ifelse(is.na(Upper2), MaxDel, Upper2)) %>%
  select(-c(Lower,Upper, MetricCurve, MetricCurve2, ModelstoUse))

head(limitsx3)

## save out
write.csv(limitsx3, "output_data/05_QGAM_FFM_Ranges_With_Adjustments.csv")


# Add to predicted value --------------------------------------------------

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
names(limitsx3)

allLims <- full_join(limitsx3, AllDF, by = c("Index", "Hydro_endpoint", "Threshold"), relationship = "many-to-many")
names(allLims)

# Get categories ----------------------------------------------------------

## define if hydro is within mod limits & bio above threshold
imps <- allLims %>%
  group_by(metric, Index, Hydro_endpoint, Threshold, BioThresh, masterid, COMID, SmoothingFunc, Quantile) %>%
  mutate(WithinHydroLimits = ifelse(deltah_final <= Upper2 & deltah_final >= Lower2, "Within", "NotWithin")) %>%  ## within hydro limits
  mutate(WithinBioLimits = ifelse(IndexValue >= BioThresh, "Within", "NotWithin"))  %>%## above bio thresholds
  mutate(Result = case_when((WithinHydroLimits == "Within" & WithinBioLimits == "Within") ~ "NoImpact", ## add impact
                            (WithinHydroLimits == "Within" & WithinBioLimits == "NotWithin") ~ "BioImpact", ##
                            (WithinHydroLimits == "NotWithin" & WithinBioLimits == "Within") ~ "HydroImpact",
                            (WithinHydroLimits == "NotWithin" & WithinBioLimits == "NotWithin") ~ "BothImpact")) %>%
  mutate(Result = factor(Result, levels = c("NoImpact", "BioImpact", "HydroImpact", "BothImpact"))) # %>%


write.csv(imps, "output_data/05_impact_ffm_bio.csv")

unique(imps$Threshold)

## NAs are from where we have data in raw DF but not channel engineering - remove
## others are when we have no asci value - remove

imps <- imps %>%
  drop_na(metric, WithinBioLimits)

## count NAs in result
sum(is.na(imps$Result)) ## 0

## tally of impact per ffm

tallyImpact <- imps %>%
  group_by(Index, Hydro_endpoint, Flow.Metric.Name,Threshold, Result) %>%
  distinct() %>%
  tally() %>%
  # drop_na(Result) %>%
  mutate(PercChans = (n/sum(n)*100))

tallyImpact

write.csv(tallyImpact, "output_data/05_count_impact.csv")

tallyImpact <- read.csv("output_data/05_count_impact.csv")

# Ma