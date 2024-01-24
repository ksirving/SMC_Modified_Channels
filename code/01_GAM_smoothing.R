## testing GAm smoothing/wigglyness
## a significant smooth term is one where you can not draw a horizontal line through the 95% confidence interval.

library(gam)
library(tidylog)
library(tidyverse)
library(sf)

library(mgcv)

mcycle <- MASS::mcycle
gam_mod <- mgcv::gam(accel ~ s(times), data=mcycle)

# Extract the model coefficients
coef(gam_mod)

## directory for figures
out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/SMC_Modified_Channels/Figures/"


# Upload & Format Data -------------------------------------------------------------

load(file = "output_data/00_bugs_algae_flow_joined_by_masterid.RData")
head(AllData)

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
# bio_h_summary <- bio_h_summary[-c(1:9),]
bio_h_summary
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
    select(Metric, MetricValue, deltah_final, Class2, channel_engineering_class) %>% ## only metrics needed
    drop_na(deltah_final, MetricValue) %>%
    # filter(deltah_final > 0) %>%
    filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
    distinct()
  
  # mydat
  mydat<-mydat[order(mydat$MetricValue),] ## order by csci value
  # try(glm(Condition~deltah_final, family=binomial(link="logit"), data=mydat), silent=T) ### glm
  
  # gam(MetricValue~s(deltah_final, k=9), family = Gamma(link = "log"), data = mydat)
  
   mgcv::gam(MetricValue~s(deltah_final, k=smet), method = "REML", data = mydat)
  
  
})

mod <- summary(gam.lm[[i]])


## save models
save(gam.lm, file = "output_data/01_csci_asci_all_gam_flow_combined_all_smoothing_functions.RData")

### get rsqds and pvals

for(i in 1:length(gam.lm)) {
  
  mod <- gam.lm[[i]]
  modx <- summary(mod)

  bio_h_summary$AIC[i] <- mod$aic 
  bio_h_summary$EDF[i] <-  modx$s.table[1] ## effective degrees of freedom
  
  bio_h_summary$PValue[i] <- modx$s.table[4] ## pvalue
  bio_h_summary$McFaddensR2[i] <- 1-mod$deviance/mod$null.deviance ## r2
  # bio_h_summary$n[i] <- mod$df[2]+1
  
  
}


## save configs and r sqds
save(bio_h_summary, file="output_data/01_csci_asci_gam_rsqds_combined_all_smoothing_functions.RData")
bio_h_summary
## make df of predicted values to predict on - need to be different for each flow metric

## blank df
DF <- NULL
DF <- as.data.frame(DF)


### get predictions and fitted values
for(i in 1:length(gam.lm)) {
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  smet<-as.numeric(bio_h_summary[i,"smooth_funcs"])
  # dmet<-as.character(bio_h_summary[i,"Direction.alt"])
  # cmet<-as.character(bio_h_summary[i,"channel.type"])
  
  
  
  ## get model, predict, extract all data and categories
  mod <- gam.lm[[i]]
  predictedVals <- predict(mod,  type = "response")

  DFX <- cbind(mod$model, as.data.frame(predictedVals)) %>%
    rename(Value = deltah_final) %>%
    mutate(Variable = tmet, SmoothingK = smet, Metric = bmet) 
  
  DF <- bind_rows(DF, DFX)
  
}
head(DFX)

## change back to numeric
DF <- DF %>%
  mutate(HydroValue = as.numeric(Value), BioValue = as.numeric(MetricValue),
         predictedVals = as.numeric(predictedVals)) %>%
  rename(hydro.endpoints = Variable) %>%
  select(-Value, -MetricValue)

str(DF)
names(DF)
### predicted figures
bio <- unique(DF$Metric)
bio
mets <- unique(DF$hydro.endpoints)
b=1
m=1

for(b in 1:length(bio)) {
  
  DF1 <- DF %>%
    filter(Metric == bio[b])
  
  # dataRaw <- AllDataLong2 %>%
  #   filter(Metric == bio[b])
    

  for(m in 1:length(mets)) {
    
    # dataRawx <- dataRaw %>%
    #   filter(hydro.endpoints == mets[m])
    
    T1 <- (ggplot(subset(DF1, hydro.endpoints == mets[m]), aes(y=predictedVals, x=HydroValue)) +
             # geom_point( data = dataRawx, aes(x=deltah_final, y = MetricValue, colour = channel_engineering_class)) +
             geom_line( linewidth = 1)+
             # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
             # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
             facet_wrap(~SmoothingK, scales = "free") +
             scale_x_continuous(name=paste(mets[m])) +
             scale_y_continuous(name = paste0(bio[b], " Score"))) +
      theme_bw()
    # theme(legend.position = "none"))
    T1
    
    file.name1 <- paste0(out.dir, "01_", bio[b], "_", mets[m], "_flow_response_predicted_gam_combined_all_smoothing.jpg")
    ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
    
  }
  
  
}
