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


# Load steppe birds trait data:
Species_traits <- read.csv("data/TraitsSteppebirdsFull.csv", stringsAsFactors = F)
str(Species_traits)
names(Species_traits)

# Removing columns related to common names, family, order, references and some traits that wont be worth to analyse 
# such as: Trophic level (Trophic niche is more specific), Maximum longevity (most values are the same as longevity), 
# relative_brain_size (removing now to calculate after imputation of brain size NAs values)
Species_traits <- Species_traits[,c("binomial", "ESP_LAT", "Beak.Length_Culmen", "Beak.Length_Nares",
                                    "Beak.Width", "Beak.Depth", "Tarsus.Length", "Wing.Length", "Kipps.Distance",
                                    "Secondary1", "Hand.Wing.Index", "Tail.Length", "Mass", "Habitat", "Habitat.Density",
                                    "Migration", "Trophic.Level", "Trophic.Niche", "Primary.Lifestyle", "Range.Size",                
                                    "hab_breadth", "adult_length_cm", "female_maturity_d", "clutch_size_n", "clutches_per_y",
                                    "maximum_longevity_y", "egg_mass_g", "incubation_d", "fledging_age_d", "longevity_y", 
                                    "birth_or_hatching_weight_g", "brood_value", "annual_survival", "PopulationDensity_ind_km2",
                                    "brain_size_g", "relative_brain_size", "eye_axial_length", "eye_tranverse_diameter", 
                                    "activity", "Nest_type", "Foraging", "Degree_of_development", "Development_continuum",
                                    "Gregariousness", "Mating_system")]

Species_traits$Habitat <- as.factor(Species_traits$Habitat)
Species_traits$Habitat.Density <- as.factor(Species_traits$Habitat.Density)
Species_traits$Migration <- as.factor(Species_traits$Migration)
Species_traits$Trophic.Niche <- as.factor(Species_traits$Trophic.Niche)
Species_traits$Primary.Lifestyle <- as.factor(Species_traits$Primary.Lifestyle)
Species_traits$activity <- as.factor(Species_traits$activity)
Species_traits$Nest_type <- as.factor(Species_traits$Nest_type)
Species_traits$Foraging <- as.factor(Species_traits$Foraging)
Species_traits$Degree_of_development <- as.factor(Species_traits$Degree_of_development)
Species_traits$Gregariousness <- as.factor(Species_traits$Gregariousness)
Species_traits$Mating_system <- as.factor(Species_traits$Mating_system)

# Counting the number of NAs per each variable
nasTraits = colSums(is.na(Species_traits))

nasTraits
# Just a few NAs: PopDensity_ind_km2(3), birth_or_hatching_weight_g(10), brain size(8), hab_breadth(2), 
# degree of development(1)               

# Imputation for NAs
set.seed(42)
imp <- missForest(Species_traits[,c(2:33)], verbose=FALSE)
Species_traits_imp <- cbind(Species_traits$ESP_LAT, imp$ximp)

names (Species_traits_imp)[1] = "ESP_LAT"

# Now calculating relative brain size to body mass using residuals. Positive values relate to brains that are
# bigger than expected
model <- lm(brain_size_g ~ Mass, data = Species_traits_imp)
residuals  <-  model$residuals
Species_traits_imp$relative_brain <- residuals
Species_traits_imp$brain_size_g  <- NULL

