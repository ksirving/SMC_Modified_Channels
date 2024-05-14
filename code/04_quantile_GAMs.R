## Quantile Gams

## tutorial
# install.packages("qgam")

library(qgam); library(MASS)

# if( suppressWarnings(require(RhpcBLASctl)) ){ blas_set_num_threads(1) } # Optional
# 
# fit <- qgam(accel~s(times, k=20, bs="ad"), 
#             data = mcycle, 
#             qu = 0.8)
# 
# # Plot the fit
# xSeq <- data.frame(cbind("accel" = rep(0, 1e3), "times" = seq(2, 58, length.out = 1e3)))
# pred <- predict(fit, newdata = xSeq, se=TRUE)
# plot(mcycle$times, mcycle$accel, xlab = "Times", ylab = "Acceleration", ylim = c(-150, 80))
# lines(xSeq$times, pred$fit, lwd = 1)
# lines(xSeq$times, pred$fit + 2*pred$se.fit, lwd = 1, col = 2)
# lines(xSeq$times, pred$fit - 2*pred$se.fit, lwd = 1, col = 2)
# 
# check(fit$calibr, 2)


# Modified channel data ---------------------------------------------------

library(tidylog)
library(tidyverse)
library(sf)
library(gam)
library(mapview)

out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/SMC_Modified_Channels/Figures/"

getwd()

# Upload & Format Data -------------------------------------------------------------

load(file = "output_data/00_bugs_algae_flow_joined_by_masterid.RData")
head(AllData)

# test <- AllData %>%
#   filter(Metric == "csci", flow_metric == "delta_q99")
# 
# head(test)
# 
# ggplot(test, aes(x = deltah_final, y = MetricValue)) +
#   geom_point() +
#   # xlim(0,1000) +
#   geom_smooth() 


## how many bio sites 
length(unique(AllData$masterid)) ## 479

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

## quantiles
quants <- c(0.99, 0.90, 0.70, 0.50)
# quantsL <- c(0.05, 0.1, 0.15, 0.20, 0.25,0.25)

### make grid
bio_h_summary<-  expand.grid(biol.endpoints=biol.endpoints,flow.endpoints=flow.endpoints, 
                              quants = quants, smooth_funcs = smooth_funcs, stringsAsFactors = F)

bio_h_summary

## create df of model configurations

head(AllDataLong2)
i=43
## model of each configuration
gam.lm <-lapply(1:nrow(bio_h_summary), function(i)
{
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  qmet<-as.numeric(bio_h_summary[i,"quants"])
  smet<-as.numeric(bio_h_summary[i,"smooth_funcs"])
  
  mydat<-AllDataLong2 %>%
    filter(Metric == bmet,
           flow_metric == tmet) %>%
    dplyr::select(Metric, MetricValue, deltah_final, channel_engineering_class) %>% ## only metrics needed
    drop_na( deltah_final) %>%
    # filter(deltah_final > 0) %>%
    filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
    distinct()
  
  
  mydat<-mydat[order(mydat$MetricValue),] ## order by csci value
  mydat
  # try(glm(Condition~deltah_final, family=binomial(link="logit"), data=mydat), silent=T) ### glm
  # gam(MetricValue~s(deltah_final, smet), family = Gamma(link = "log"), data = mydat)
  
  qgam(MetricValue~s(deltah_final, k=smet), 
       data = mydat, 
        qu = qmet)
  

})

# mod1 <- qgam(MetricValue~s(deltah_final,k=smet, bs="ad"), 
#      data = mydat, 
#      qu = qmet)
# 
# plot(mod1)
# 
# # Error in seq.default(xl - dx * (m[1] + 1), xu + dx * (m[1] + 1), length = nk +  : 
# #                        'from' must be a finite number
# 
# qgam
# 
# mod2 <- gam(MetricValue~s(deltah_final), family = Gamma(link = "log"), data = mydat)
# 
# plot(mod2)

## save models
save(gam.lm, file = "output_data/04_csci_asci_quantile_gam.RData")
gam.lm
## load
load(file = "output_data/04_csci_asci_quantile_gam.RData")

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

  
  bio_h_summary$DevExplained[i] <- mod$dev.expl ##1-mod$deviance/mod$null.deviance ## mcfaddens r2
  bio_h_summary$PValue[i] <- mod$s.table[4]
  bio_h_summary$R2[i] <- mod$r.sq
  bio_h_summary$n[i] <- mod$n
  bio_h_summary$EDF[i] <- mod$edf
  
}

## save configs and r sqds
save(bio_h_summary, file="output_data/04_csci_asci_quantile_gam_rsqds.RData")

## load
load(file="output_data/04_csci_asci_quantile_gam_rsqds.RData")
bio_h_summary

## plot k vals with r2 

head(bio_h_summary)

P1 <- ggplot(data = bio_h_summary, aes(y=R2, x=quants, group = flow.endpoints, colour = flow.endpoints)) +
  # geom_point(size=0.2) +
  geom_line( linewidth = 1)+
  # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
  # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
  facet_wrap(~biol.endpoints, scales = "free") +
  facet_grid(rows = vars(smooth_funcs), cols = vars(biol.endpoints), scales = "free") +
  scale_x_continuous(name= "Quantiles") +
  scale_y_continuous(name = "McFaddens (?) Rsq") 
# theme(legend.position = "none"))

P1

file.name1 <- paste0(out.dir, "04_quants_rsq_smooths.jpg")
ggsave(P1, filename=file.name1, dpi=300, height=5, width=7.5)

## plot k vals with AIC

P2 <- ggplot(data = bio_h_summary, aes(y=DevExplained, x=quants, group = flow.endpoints, colour = flow.endpoints)) +
  # geom_point(size=0.2) +
  geom_line( linewidth = 1)+
  # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
  # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
  facet_wrap(~biol.endpoints, scales = "free") +
  facet_grid(rows = vars(smooth_funcs), cols = vars(biol.endpoints), scales = "free") +
  scale_x_continuous(name= "Quantiles") +
  scale_y_continuous(name = "Deviance Explained") 
# theme(legend.position = "none"))

P2

file.name1 <- paste0(out.dir, "04_Quantiles_DevExpl_smooths.jpg")
ggsave(P2, filename=file.name1, dpi=300, height=5, width=7.5)

## make df of predicted values to predict on - need to be different for each flow metric

## blank df
DF <- NULL
DF <- as.data.frame(DF)

## shorten to test predictions
# bio_h_summary <- bio_h_summary[1:10,]
# bio_h_summary
# dim(bio_h_summary)

### get predictions and fitted values
for(i in 1:length(gam.lm)) {
  
  print(paste0("Model ", i))
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  qmet<-as.numeric(bio_h_summary[i,"quants"])
  smet<-as.numeric(bio_h_summary[i,"smooth_funcs"])
  
  ## get input data
  mydat<-AllDataLong2 %>%
    filter(Metric == bmet,
           flow_metric == tmet) %>%
    select(Metric, MetricValue, deltah_final, Class2, channel_engineering_class) %>% ## only metrics needed
    drop_na(deltah_final, MetricValue) %>%
    # filter(deltah_final > 0) %>%
    filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
    distinct()
  
  incr <- mean(mydat$deltah_final)/500

  ## get new data from range of iunput data
  newdata <- as.data.frame(seq(min(mydat$deltah_final), max(mydat$deltah_final), abs(incr)))
  colnames(newdata) <- "deltah_final" ## change name to match original
  
  ## get model, predict, extract all data and categories
  mod <- gam.lm[[i]]
  
  predictedVals <- predict(mod, newdata, type = "response") ## add newdata if needed

  ## join
  DFX <- cbind(newdata, as.data.frame(predictedVals)) %>%
    # rename(Value = "s(deltah_final)") %>%
    mutate(Variable = tmet,  Metric = bmet, Quant = qmet, Smooths = smet) %>%
    distinct()
  
  # head(DFX)
  
  DF <- bind_rows(DF, DFX)
  
}

head(DF)

## change back to numeric
DF <- DF %>%
  mutate(deltah_final = as.numeric(deltah_final),
         predictedVals = as.numeric(predictedVals),
         Quant = as.factor(Quant))

save(DF, file = "ignore/04_quantGams_smooths_predictions.RData")
str(DF)

### predicted figures
bio <- unique(DF$Metric)
mets <- unique(DF$Variable)
mets

b=2
m=3
## plot curve with points for each quant
## get se e.g., plot(mcycle$times, mcycle$accel, xlab = "Times", ylab = "Acceleration", ylim = c(-150, 80))
# lines(xSeq$times, pred$fit, lwd = 1)
# lines(xSeq$times, pred$fit + 2*pred$se.fit, lwd = 1, col = 2)
# lines(xSeq$times, pred$fit - 2*pred$se.fit, lwd = 1, col = 2)  

## plot per quant
## use 6 k smooth to check
## 99, 90, 70, 50

for(b in 1:length(bio)) {
  
  ## filter to index
  DF1 <- DF %>%
    filter(Metric == bio[b], Smooths == 3)

  ## get thresholds for ref and HB channels
  if(bio[b] == "asci") {
    
    refT <- 0.86
    hbT <- 0.87
    
  } else {
    
    refT <- 0.79
    hbT <- 0.67
  }
  
  # names(AllDataLong2)
  # head(DF1)

  ### get observed data points for overlay
  ptsbio <- AllDataLong2 %>%
    filter(Metric == bio[b])
  
# str(ptsbio)
  
  for(m in 1:length(mets)) {
    
    ptsbiox <- ptsbio %>%
      filter(flow_metric == mets[m]) %>%
      mutate(channel_engineering_class = as.factor(channel_engineering_class))
    
    # ptsbiox

    T1 <- ggplot() +
             geom_smooth(data = subset(DF1, Variable == mets[m]), aes(y=predictedVals, x=deltah_final, col = Quant), linewidth = 1)+
             geom_point(data = ptsbiox, aes(x=deltah_final, y = MetricValue,
                                            col = channel_engineering_class)) +
             # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
             geom_hline(yintercept = refT,  linetype="dashed", linewidth=0.5, color = "grey50") +
             geom_hline(yintercept = hbT,  linetype="dashed", linewidth=0.5, color = "red") +
             geom_vline(xintercept = 0) +
             # facet_wrap(~Smooths, scales = "free") +
             scale_x_continuous(name=paste(mets[m])) +
             scale_y_continuous(name = paste0("Index Score  (", bio[b], ")"))
             
    
    # theme(legend.position = "none"))
    
    T1
    
    file.name1 <- paste0(out.dir, "04_", bio[b], "_", mets[m], "_flow_response_predicted_gam_combined_all_quants.jpg")
    ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
    
  }
  
  
  
  
}


# test <- subset(DF1, Variable == mets[m])
# test
## smoothing per lowest quants

# for(b in 1:length(bio)) {
# 
# DF1 <- DF %>%
#   filter(Metric == bio[b], Quant == 0.7)
# 
# # head(DF1)
# 
# for(m in 1:length(mets)) {
#   
#   T1 <- (ggplot(subset(DF1, Variable == mets[m]) %>% arrange(predictedVals), aes(y=predictedVals, x=deltah_final)) +
#            geom_smooth( linewidth = 1)+ 
#            # geom_point() +
#            # geom_point( aes(x=deltah_final, y = MetricValue, colour = channel_engineering_class)) +
#            # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
#            # geom_hline(yintercept = 0.6,  linetype="dashed", linewidth=0.5, color = "grey50") +
#            facet_wrap(~SmoothFunc, scales = "free") +
#            scale_x_continuous(name=paste(mets[m])) +
#            scale_y_continuous(name = paste0("Probability of Good ", bio[b]))) 
#   # theme(legend.position = "none"))
#   
#   T1
#   
#   file.name1 <- paste0(out.dir, "04_", bio[b], "_", mets[m], "_flow_response_predicted_gam_combined_all_smoothingKs_at_low_quantile.jpg")
#   ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
#   
# }
# 
# 
# 
# 
# }

### david code

met.plot<-met.df %>% 
  ggplot(., aes(x=comp_val, y=csci))+
  theme_bw()+
  theme(panel.grid = element_blank(),
  )+
  geom_point(shape=21, size=1.5, fill="#CCCCCC")+
  geom_quantile(aes(color=as_factor(after_stat(quantile)), linetype=as_factor(after_stat(quantile))),method="rq", 
                quantiles=c(0.9, 0.8,0.5,0.2,0.1), linewidth=1.25, key_glyph="rect")+
  geom_smooth(method = "lm", color="#990000", linewidth=1.25,se=FALSE)+
  labs(x=xxx, y="CSCI Score", color="Estimate\nType", linetype="Estimate\nType")+
  scale_color_manual(limits=c("0.9", "0.8","0.5","0.49","0.2","0.1"),
                     values=c("#3399FF", "#663399", "#339900","#990000","#663399", "#3399FF"), 
                     breaks=c(0.9, 0.8,0.5,0.49,0.2,0.1), 
                     labels=c("90th", "80th", "Median","Mean", "20th", "10th"))+
  scale_linetype_manual(limits=c("0.9", "0.8","0.5","0.49", "0.2","0.1"),
                        values=c(2,4,6,1,4,2),
                        breaks=c(0.9, 0.8,0.5, 0.49,0.2,0.1),
                        labels=c("90th", "80th", "Median","Mean", "20th", "10th"))+
  ggtitle(paste("All Comparators", xxx, "Quantile Regression", sep=" "))+
  NULL


# Testing Peak augmentation -----------------------------------------------


