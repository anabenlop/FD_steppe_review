# Functional Diversity analysis for steppe birds
# Date: 05/08/2024
# Author: Ana Benítez


library(dplyr)
library(tidyr)
library(tibble)
library(mFD)
library(ggplot2)

#Import data

sp_data <- read.csv("data/TraitsSteppebirdsFull.csv")

biomass_matrix <- fauna %>% group_by(SITE, REFERENCE_TAXA) %>% 
  summarize(BIOMASS=sum(BIOMASS_KG)) %>%
  spread(key = REFERENCE_TAXA,value = BIOMASS)%>% #biomass reported in kg
  replace(is.na(.), 0) %>% #replace NAs with 0
  tibble::column_to_rownames(var = "SITE")

rm(fauna) #clean environment

#Create matrix depicting vertebrate's traits
