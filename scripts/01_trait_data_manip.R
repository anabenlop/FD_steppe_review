# Manipulate, explore and impute trait data of Iberian steppe birds 

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
data <- read.csv("data/TraitsSteppebirdsFull.csv", stringsAsFactors = F)
data <- data[-c(42:1012),]
str(data)
names(data)

data$Type <- ifelse(data$Type == "Strict Esteparias Estrictas", "Strict", "NonStrict")

# Create matrix depicting strict and non-strict steppe bird species
sp_matrix <- data %>% group_by(Type, binomial) %>% 
  summarize(n = n()) %>%
  spread(key = binomial,value = n) %>% 
  replace(is.na(.), 0) %>% #replace NAs with 0
  tibble::column_to_rownames(var = "Type")

# Removing columns related to common names, family, order, references 
traits <- data[,c("binomial", "ESP_LAT", "Beak.Length_Culmen", "Beak.Length_Nares","Beak.Width", "Beak.Depth", 
                                    "Tarsus.Length", "Wing.Length", "Kipps.Distance", "Hand.Wing.Index", "Tail.Length", "Mass", 
                                    "Habitat", "Habitat.Density","Migration", "Trophic.Level", "Trophic.Niche", "Primary.Lifestyle", "Range.Size",                
                                    "hab_breadth", "adult_length_cm", "female_maturity_d", "clutch_size_n", "clutches_per_y",
                                    "maximum_longevity_y", "egg_mass_g", "incubation_d", "fledging_age_d", "longevity_y", 
                                    "birth_or_hatching_weight_g", "annual_survival", "PopulationDensity_ind_km2",
                                    "brain_size_g", "eye_axial_length", "eye_tranverse_diameter", 
                                    "activity", "Nest_type", "Foraging", "Degree_of_development", "Development_continuum",
                                    "Gregariousness", "Mating_system")]


# change character to factor for categorical variables
traits$Habitat <- as.factor(traits$Habitat)
traits$Habitat.Density <- as.factor(traits$Habitat.Density)
traits$Migration <- as.factor(traits$Migration)
traits$Trophic.Level <- as.factor(traits$Trophic.Level)
traits$Trophic.Niche <- as.factor(traits$Trophic.Niche)
traits$Primary.Lifestyle <- as.factor(traits$Primary.Lifestyle)
traits$activity <- as.factor(traits$activity)
traits$Nest_type <- as.factor(traits$Nest_type)
traits$Foraging <- as.factor(traits$Foraging)
traits$Degree_of_development <- as.factor(traits$Degree_of_development)
traits$Gregariousness <- as.factor(traits$Gregariousness)
traits$Mating_system <- as.factor(traits$Mating_system)

# explore trait distribution
ggplot(traits) + geom_histogram(aes(Beak.Length_Culmen))
ggplot(traits) + geom_histogram(aes(Beak.Length_Nares))
ggplot(traits) + geom_histogram(aes(Beak.Width))
ggplot(traits) + geom_histogram(aes(Beak.Depth))
ggplot(traits) + geom_histogram(aes(Tarsus.Length))
ggplot(traits) + geom_histogram(aes(Wing.Length))
ggplot(traits) + geom_histogram(aes(Kipps.Distance))
ggplot(traits) + geom_histogram(aes(Hand.Wing.Index))
ggplot(traits) + geom_histogram(aes(Tail.Length))
ggplot(traits) + geom_histogram(aes(Mass))
ggplot(traits) + geom_histogram(aes(hab_breadth))
ggplot(traits) + geom_histogram(aes(Range.Size))
ggplot(traits) + geom_histogram(aes(adult_length_cm))
ggplot(traits) + geom_histogram(aes(female_maturity_d))
ggplot(traits) + geom_histogram(aes(clutch_size_n))
ggplot(traits) + geom_histogram(aes(clutches_per_y))
ggplot(traits) + geom_histogram(aes(maximum_longevity_y))
ggplot(traits) + geom_histogram(aes(egg_mass_g))
ggplot(traits) + geom_histogram(aes(incubation_d))
ggplot(traits) + geom_histogram(aes(fledging_age_d))
ggplot(traits) + geom_histogram(aes(longevity_y))
ggplot(traits) + geom_histogram(aes(birth_or_hatching_weight_g))
ggplot(traits) + geom_histogram(aes(annual_survival))
ggplot(traits) + geom_histogram(aes(PopulationDensity_ind_km2))
ggplot(traits) + geom_histogram(aes(brain_size_g))
ggplot(traits) + geom_histogram(aes(eye_axial_length))
ggplot(traits) + geom_histogram(aes(eye_tranverse_diameter))
ggplot(traits) + geom_histogram(aes(Development_continuum))


# check relationship between development categorical and continuous index
ggplot(traits) + geom_boxplot((aes(Degree_of_development, Development_continuum))) # makes more sense to use the continuum?

# Counting the number of NAs per each variable
nasTraits = colSums(is.na(traits))

nasTraits
# Just a few NAs: PopDensity_ind_km2(1), birth_or_hatching_weight_g(17), brain size(15), hab_breadth(2), eye_axial_length (20), eye_tranverse_diameter (25),

# # Missing data plot
Amelia::missmap(traits[,-c(1:2)])

# Imputation for NAs
set.seed(42)
imp <- missForest(traits[,c(3:42)], verbose=TRUE)
traits_imp <- cbind(traits$binomial, imp$ximp)

names (traits_imp)[1] = "binomial"

# Now calculating relative brain size to body mass using residuals. Positive values relate to brains that are
# bigger than expected
model <- lm(log10(brain_size_g) ~ log10(Mass), data = traits_imp)
summary(model)
residuals  <-  model$residuals
traits_imp$relative_brain <- residuals

ggplot(traits_imp) + geom_histogram(aes(relative_brain))

# explore trait distribution with imputed values
ggplot(traits_imp) + geom_histogram(aes(Beak.Length_Culmen))
ggplot(traits_imp) + geom_histogram(aes(Beak.Length_Nares))
ggplot(traits_imp) + geom_histogram(aes(Beak.Width))
ggplot(traits_imp) + geom_histogram(aes(Beak.Depth))
ggplot(traits_imp) + geom_histogram(aes(Tarsus.Length))
ggplot(traits_imp) + geom_histogram(aes(Wing.Length))
ggplot(traits_imp) + geom_histogram(aes(Kipps.Distance))
ggplot(traits_imp) + geom_histogram(aes(Hand.Wing.Index))
ggplot(traits_imp) + geom_histogram(aes(Tail.Length))
ggplot(traits_imp) + geom_histogram(aes(Mass))
ggplot(traits_imp) + geom_histogram(aes(hab_breadth))
ggplot(traits_imp) + geom_histogram(aes(Range.Size))
ggplot(traits_imp) + geom_histogram(aes(adult_length_cm))
ggplot(traits_imp) + geom_histogram(aes(female_maturity_d))
ggplot(traits_imp) + geom_histogram(aes(clutch_size_n))
ggplot(traits_imp) + geom_histogram(aes(clutches_per_y))
ggplot(traits_imp) + geom_histogram(aes(maximum_longevity_y))
ggplot(traits_imp) + geom_histogram(aes(egg_mass_g))
ggplot(traits_imp) + geom_histogram(aes(incubation_d))
ggplot(traits_imp) + geom_histogram(aes(fledging_age_d))
ggplot(traits_imp) + geom_histogram(aes(longevity_y))
ggplot(traits_imp) + geom_histogram(aes(birth_or_hatching_weight_g))
ggplot(traits_imp) + geom_histogram(aes(annual_survival))
ggplot(traits_imp) + geom_histogram(aes(PopulationDensity_ind_km2))
ggplot(traits_imp) + geom_histogram(aes(brain_size_g))
ggplot(traits_imp) + geom_histogram(aes(relative_brain))
ggplot(traits_imp) + geom_histogram(aes(eye_axial_length))
ggplot(traits_imp) + geom_histogram(aes(eye_tranverse_diameter))
ggplot(traits_imp) + geom_histogram(aes(Development_continuum))


# Standardize continuous traits
continuous <- c(2:11, 18:34, 39, 42) #vector of columns with continuous traits
traits_imp <- BAT::standard(data.frame(traits_imp), method = "standard", convert = continuous) 

# traits_imp$Beak.Length_Culmen<- as.numeric(scale(log10(traits_imp$Beak.Length_Culmen), center=TRUE, scale=TRUE))
# traits_imp$Beak.Length_Nares<- as.numeric(scale(log10(traits_imp$Beak.Length_Nares), center=TRUE, scale=TRUE))
# traits_imp$Beak.Width<- as.numeric(scale(log10(traits_imp$Beak.Width), center=TRUE, scale=TRUE))
# traits_imp$Beak.Depth<- as.numeric(scale(log10(traits_imp$Beak.Depth), center=TRUE, scale=TRUE))
# traits_imp$Tarsus.Length<- as.numeric(scale(log10(traits_imp$Tarsus.Length), center=TRUE, scale=TRUE))
# traits_imp$Wing.Length<- as.numeric(scale(log10(traits_imp$Wing.Length), center=TRUE, scale=TRUE))
# traits_imp$Kipps.Distance <- as.numeric(scale(log10(traits_imp$Kipps.Distance), center=TRUE, scale=TRUE))
# traits_imp$Hand.Wing.Index <- as.numeric(scale(log10(traits_imp$Hand.Wing.Index), center=TRUE, scale=TRUE))
# traits_imp$Tail.Length <- as.numeric(scale(log10(traits_imp$Tail.Length), center=TRUE, scale=TRUE))
# traits_imp$Mass <- as.numeric(scale(log10(traits_imp$Mass), center=TRUE, scale=TRUE))
# traits_imp$hab_breadth <- as.numeric(scale(log10(traits_imp$hab_breadth), center=TRUE, scale=TRUE))
# traits_imp$Range.Size <- as.numeric(scale(log10(traits_imp$Range.Size), center=TRUE, scale=TRUE))
# traits_imp$adult_length_cm <- as.numeric(scale(log10(traits_imp$adult_length_cm), center=TRUE, scale=TRUE))
# traits_imp$female_maturity_d <- as.numeric(scale(log10(traits_imp$female_maturity_d), center=TRUE, scale=TRUE))
# traits_imp$clutch_size_n <- as.numeric(scale(log10(traits_imp$clutch_size_n), center=TRUE, scale=TRUE))
# traits_imp$clutches_per_y <- as.numeric(scale(log10(traits_imp$clutches_per_y), center=TRUE, scale=TRUE))
# traits_imp$maximum_longevity_y <- as.numeric(scale(log10(traits_imp$maximum_longevity_y), center=TRUE, scale=TRUE))
# traits_imp$egg_mass_g <- as.numeric(scale(log10(traits_imp$egg_mass_g), center=TRUE, scale=TRUE))
# traits_imp$incubation_d <- as.numeric(scale(log10(traits_imp$incubation_d), center=TRUE, scale=TRUE))
# traits_imp$fledging_age_d <- as.numeric(scale(log10(traits_imp$fledging_age_d), center=TRUE, scale=TRUE))
# traits_imp$longevity_y <- as.numeric(scale(log10(traits_imp$longevity_y), center=TRUE, scale=TRUE))
# traits_imp$birth_or_hatching_weight_g <- as.numeric(scale(log10(traits_imp$birth_or_hatching_weight_g), center=TRUE, scale=TRUE))
# traits_imp$annual_survival <- as.numeric(scale(log10(traits_imp$annual_survival), center=TRUE, scale=TRUE))
# traits_imp$PopulationDensity_ind_km2 <- as.numeric(scale(log10(traits_imp$PopulationDensity_ind_km2), center=TRUE, scale=TRUE))
# traits_imp$brain_size_g <- as.numeric(scale(log10(traits_imp$brain_size_g), center=TRUE, scale=TRUE))
# traits_imp$eye_tranverse_diameter <- as.numeric(scale(log10(traits_imp$eye_tranverse_diameter), center=TRUE, scale=TRUE))

# Check Collinearity between continuous traits
M = cor(traits_imp[,continuous])
corrplot(M, method = 'number') # color
dev.off()

# Check Collinearity
psych::pairs.panels(traits_imp[,colnames(traits_imp[continuous])])

# remove redundant traits
traits_imp$brain_size_g  <- NULL #we already have relative brain
traits_imp$longevity_y  <- NULL #max longevity included
traits_imp$Trophic.Level  <- NULL #Trophic.niche more informative
traits_imp$Degree_of_development  <- NULL #Development continuum more informative

#Loading metadata: name of each variable and type of variable: Q (quantitative), N (nominal)
traits_cat <- read.csv("data/Traits_categories.csv", stringsAsFactors = F)

# remove redundant traits here as well
traits_cat <- traits_cat[-c(31,27, 14, 37),]

# save trait matrix and trait categories for further analyses in script 02
write.csv(traits_imp, "clean_data//traits_clean.csv", row.names = F)
write.csv(traits_cat, "clean_data/traits_cat_clean.csv", row.names = F)

# End of script ------
