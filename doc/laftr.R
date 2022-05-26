## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----import, warning=FALSE,message=FALSE--------------------------------------
library(readxl)
library(dplyr)

library(laftr)

# load the U concentrations
sample <- read_xlsx('../data/example.xlsx', sheet = 'sample')  %>%
    rename(U = UCa, sU = sUCa)

# load the calibration measurements for zeta factor
zeta.measurements <- read_xlsx('../data/example.xlsx', sheet = 'standard') %>%
    rename(U = UCa, sU = sUCa)

## ----zeta1--------------------------------------------------------------------
session.zeta <- zeta_ICP(zeta.measurements, spotsize = 40, mineral = 'apatite', standard = 'Dur') 

## ----zeta2--------------------------------------------------------------------
head(session.zeta$data) %>% as_tibble()

## ----zeta3--------------------------------------------------------------------
head(session.zeta$results) %>% as_tibble()

## ----zeta4--------------------------------------------------------------------
session.zeta$zeta %>% as_tibble()

## ----res----------------------------------------------------------------------
result <- age_ICP(sample, spotsize = 40, zeta = session.zeta$zeta)

## ----res.data-----------------------------------------------------------------
head(result$data) %>% as_tibble()

## ----ages---------------------------------------------------------------------
head(result$ages) %>% as_tibble()

## ----ages.pooled--------------------------------------------------------------
result$age.pooled

## ----chisq--------------------------------------------------------------------
result$stats

## ----isoplotr-----------------------------------------------------------------
#library(IsoplotR)
data.iso <- as.fissiontracks(sample, spotsize = 40, zeta = session.zeta$zeta)
ages.iso <- as.fissiontracks(result$ages, ages = TRUE)

## ----kde----------------------------------------------------------------------
IsoplotR::kde(ages.iso)

## ----radial-------------------------------------------------------------------
IsoplotR::radialplot(ages.iso)

## ----central------------------------------------------------------------------
IsoplotR::central(ages.iso)

## ----w.mean-------------------------------------------------------------------
IsoplotR::weightedmean(ages.iso)

## ----peaks--------------------------------------------------------------------
IsoplotR::peakfit(ages.iso)

