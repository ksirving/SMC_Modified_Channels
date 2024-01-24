## level 1 - regional analysis:  Flow targets derived for standard csci thresholds - do they up or down by category

library(gam)
library(tidylog)
library(tidyverse)
library(sf)

## directory for figures
out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/SMC_Modified_Channels/Figures/"

## workflow
## targets for standard threshold all streams
## targets for modified threshold all streams
## all 4 outcomes for each
## How much does it change per category?

# Upload Data -------------------------------------------------------------

## thresholds for each FFM and index
df <- read.csv("output_data/01_delta_thresholds_GAMs.csv")
df

## all data: index score, ffm by site
load(file = "output_data/00_bugs_algae_flow_joined_by_masterid.RData")
head(AllData)
