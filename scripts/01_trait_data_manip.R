# Calculate functional space of Iberian steppe birds 

# load libraries
# library(sf)
library(dplyr)
library(mFD)
# library(picante)
library(tidyverse)
library(corrplot)
library(caret)
library(missForest)
library(tibble)
# library(biscale)
# library(ggspatial)
# library(cowplot)
library(FD)
library(ggplot2)

# clean environment
rm(list = ls())

# Load steppe birds trait data:
traits <- read.csv("data/TraitsSteppebirdsFull.csv", stringsAsFactors = F)
traits <- traits[-c(42:1012),]
str(traits)
names(traits)

# Removing columns related to common names, family, order, references and some traits that wont be worth to analyse 
# such as: Trophic level (Trophic niche is more specific), Maximum longevity (most values are the same as longevity), 
# relative_brain_size (removing now to calculate after imputation of brain size NAs values)
traits <- traits[,c("binomial", "ESP_LAT", "Beak.Length_Culmen", "Beak.Length_Nares",
                                    "Beak.Width", "Beak.Depth", "Tarsus.Length", "Wing.Length", "Kipps.Distance",
                                    "Hand.Wing.Index", "Tail.Length", "Mass", "Habitat", "Habitat.Density",
                                    "Migration", "Trophic.Level", "Trophic.Niche", "Primary.Lifestyle", "Range.Size",                
                                    "hab_breadth", "adult_length_cm", "female_maturity_d", "clutch_size_n", "clutches_per_y",
                                    "maximum_longevity_y", "egg_mass_g", "incubation_d", "fledging_age_d", "longevity_y", 
                                    "birth_or_hatching_weight_g", "annual_survival", "PopulationDensity_ind_km2",
                                    "brain_size_g", "eye_axial_length", "eye_tranverse_diameter", 
                                    "activity", "Nest_type", "Foraging", "Degree_of_development", "Development_continuum",
                                    "Gregariousness", "Mating_system")]

# check relationship between development categorical and continuous index
ggplot(traits) + geom_boxplot((aes(Degree_of_development, Development_continuum)))
ggplot(traits) + geom_histogram((aes(Development_continuum)))

traits$Habitat <- as.factor(traits$Habitat)
traits$Habitat.Density <- as.factor(traits$Habitat.Density)
traits$Migration <- as.factor(traits$Migration)
traits$Trophic.Niche <- as.factor(traits$Trophic.Niche)
traits$Primary.Lifestyle <- as.factor(traits$Primary.Lifestyle)
traits$activity <- as.factor(traits$activity)
traits$Nest_type <- as.factor(traits$Nest_type)
traits$Foraging <- as.factor(traits$Foraging)
traits$Degree_of_development <- as.factor(traits$Degree_of_development)
traits$Gregariousness <- as.factor(traits$Gregariousness)
traits$Mating_system <- as.factor(traits$Mating_system)

# Counting the number of NAs per each variable
nasTraits = colSums(is.na(traits))

nasTraits
# Just a few NAs: PopDensity_ind_km2(3), birth_or_hatching_weight_g(10), brain size(8), hab_breadth(2), 
# degree of development(1)               

# Imputation for NAs
set.seed(42)
imp <- missForest(traits[,c(2:33)], verbose=FALSE)
traits_imp <- cbind(traits$ESP_LAT, imp$ximp)

names (traits_imp)[1] = "ESP_LAT"

# Now calculating relative brain size to body mass using residuals. Positive values relate to brains that are
# bigger than expected
model <- lm(brain_size_g ~ Mass, data = traits_imp)
residuals  <-  model$residuals
traits_imp$relative_brain <- residuals
traits_imp$brain_size_g  <- NULL

