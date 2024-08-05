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
traits <- traits[,c("binomial", "ESP_LAT", "Beak.Length_Culmen", "Beak.Length_Nares","Beak.Width", "Beak.Depth", 
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
# traits_imp$brain_size_g  <- NULL

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


#log transform and scale quantitative traits
traits_imp$Beak.Length_Culmen<- as.numeric(scale(log10(traits_imp$Beak.Length_Culmen), center=TRUE, scale=TRUE))
traits_imp$Beak.Length_Nares<- as.numeric(scale(log10(traits_imp$Beak.Length_Nares), center=TRUE, scale=TRUE))
traits_imp$Beak.Width<- as.numeric(scale(log10(traits_imp$Beak.Width), center=TRUE, scale=TRUE))
traits_imp$Beak.Depth<- as.numeric(scale(log10(traits_imp$Beak.Depth), center=TRUE, scale=TRUE))
traits_imp$Tarsus.Length<- as.numeric(scale(log10(traits_imp$Tarsus.Length), center=TRUE, scale=TRUE))
traits_imp$Wing.Length<- as.numeric(scale(log10(traits_imp$Wing.Length), center=TRUE, scale=TRUE))
traits_imp$Kipps.Distance <- as.numeric(scale(log10(traits_imp$Kipps.Distance), center=TRUE, scale=TRUE))
traits_imp$Hand.Wing.Index <- as.numeric(scale(log10(traits_imp$Hand.Wing.Index), center=TRUE, scale=TRUE))
traits_imp$Tail.Length <- as.numeric(scale(log10(traits_imp$Tail.Length), center=TRUE, scale=TRUE))
traits_imp$Mass <- as.numeric(scale(log10(traits_imp$Mass), center=TRUE, scale=TRUE))
traits_imp$hab_breadth <- as.numeric(scale(log10(traits_imp$hab_breadth), center=TRUE, scale=TRUE))
traits_imp$Range.Size <- as.numeric(scale(log10(traits_imp$Range.Size), center=TRUE, scale=TRUE))
traits_imp$adult_length_cm <- as.numeric(scale(log10(traits_imp$adult_length_cm), center=TRUE, scale=TRUE))
traits_imp$female_maturity_d <- as.numeric(scale(log10(traits_imp$female_maturity_d), center=TRUE, scale=TRUE))
traits_imp$clutch_size_n <- as.numeric(scale(log10(traits_imp$clutch_size_n), center=TRUE, scale=TRUE))
traits_imp$clutches_per_y <- as.numeric(scale(log10(traits_imp$clutches_per_y), center=TRUE, scale=TRUE))
traits_imp$maximum_longevity_y <- as.numeric(scale(log10(traits_imp$maximum_longevity_y), center=TRUE, scale=TRUE))
traits_imp$egg_mass_g <- as.numeric(scale(log10(traits_imp$egg_mass_g), center=TRUE, scale=TRUE))
traits_imp$incubation_d <- as.numeric(scale(log10(traits_imp$incubation_d), center=TRUE, scale=TRUE))
traits_imp$fledging_age_d <- as.numeric(scale(log10(traits_imp$fledging_age_d), center=TRUE, scale=TRUE))
traits_imp$longevity_y <- as.numeric(scale(log10(traits_imp$longevity_y), center=TRUE, scale=TRUE))
traits_imp$birth_or_hatching_weight_g <- as.numeric(scale(log10(traits_imp$birth_or_hatching_weight_g), center=TRUE, scale=TRUE))
traits_imp$annual_survival <- as.numeric(scale(log10(traits_imp$annual_survival), center=TRUE, scale=TRUE))
traits_imp$PopulationDensity_ind_km2 <- as.numeric(scale(log10(traits_imp$PopulationDensity_ind_km2), center=TRUE, scale=TRUE))
traits_imp$brain_size_g <- as.numeric(scale(log10(traits_imp$brain_size_g), center=TRUE, scale=TRUE))
traits_imp$eye_tranverse_diameter <- as.numeric(scale(log10(traits_imp$eye_tranverse_diameter), center=TRUE, scale=TRUE))

# assess correlation between quantitative traits
M = cor(traits_imp[,c(2:11, 18:34, 39, 42)])
corrplot(M, method = 'number') # color

#Loading metadata: name of each variable and type of variable: Q (quantitative), N (nominal)
traits_cat <- read.csv("data/Traits_categories.csv", stringsAsFactors = F)

#Setting species names as row names
rownames(traits_imp) <- traits_imp$binomial
traits_imp$binomial <- NULL 

# Species traits summary:
species_traits_summ <- mFD::sp.tr.summary(
  tr_cat     = traits_cat,   
  sp_tr      = traits_imp, 
  stop_if_NA = TRUE)

#calculate functional distance
func_dist <- mFD::funct.dist(sp_tr = traits_imp, tr_cat = traits_cat, metric="gower")

# build up functional spaces: calculate with different weighting of deviation and
# with or without scaling functional distance.
func_space <- mFD::quality.fspaces(sp_dist= func_dist, maxdim_pcoa = 9,
                                   deviation_weighting = c("absolute","squared"), fdist_scaling = c(TRUE,FALSE),
                                   fdendro= "average")

#check quality metrics 
func_space$quality_fspaces #check quality metrics
apply(func_space$quality_fspaces, 2, which.min) # PCoA 9D seems to be the best

#plot quality metrics across PCoAs
library(magrittr)
func_space$"quality_fspaces" %>%
  tibble::as_tibble(rownames = "Functional_space") %>%
  tidyr::pivot_longer(cols =! Functional_space, names_to = "quality_metric", values_to = "Quality") %>%
  ggplot2::ggplot(ggplot2::aes(x = Functional_space, y = Quality, fill=quality_metric,
                               colour= quality_metric, shape = quality_metric, size=3)) +
  ggplot2::geom_point() + ggplot2::theme_minimal() + ggplot2::scale_shape_manual(values=c(1,2,6,5)) + 
  viridis::scale_color_viridis(discrete = TRUE, option="D") 


#plot quality of functional spaces with Mean Absolute Deviation 
mFD::quality.fspaces.plot(fspaces_quality= func_space,
                          quality_metric  = "mad_scaled",
                          fspaces_plot = c("tree_average", "pcoa_1d", "pcoa_2d", "pcoa_3d",
                                           "pcoa_4d", "pcoa_5d", "pcoa_6d",
                                           "pcoa_7d", "pcoa_8d", "pcoa_9d"),
                          gradient_deviation = c(neg = "#030091", nul = "#66B79C", pos = "#C6FFB7"),
                          gradient_deviation_quality = c(low = "#B7E3FF", high = "#2800B2"),
                          x_lab= "Trait-based distance")

# From XXX
# convex hull-based indices require a space with less axes than the number of species number, and their computation time increases with the number of axes 
# (to the point that functional beta-diversity indices are hardly computable in more than five dimensions). So if, for example, the best space is the one with 
# six axes while the quality index of the 4D and 5D spaces are close, keeping the 4D space will be a pragmatic choice.

# As FD indices will eventually be computed on coordinates on space (raw distance), 
# we hereafter will consider only the mean absolute-deviation metric

# Plot relation between functional axes and traits
sp_coords <- func_space$details_fspaces$sp_pc_coord

# group traits by type
morpho <- c("Beak.Length_Culmen", "Beak.Length_Nares","Beak.Width", "Beak.Depth", 
                        "Tarsus.Length", "Wing.Length", "Kipps.Distance", "Hand.Wing.Index", "Mass", # not Tail.length
                        "adult_length_cm")

demog <- c("female_maturity_d", "clutch_size_n", "clutches_per_y", "maximum_longevity_y", "egg_mass_g", "incubation_d", "fledging_age_d",  
                       "birth_or_hatching_weight_g", "annual_survival", "PopulationDensity_ind_km2") # no logenvity_y or Range.Size

develop <-  c("brain_size_g", "relative_brain", "eye_axial_length", "eye_tranverse_diameter", 
                          "Degree_of_development", "Development_continuum", "activity")

behav <-  c("Habitat", "Habitat.Density","Migration", "Trophic.Niche", "Primary.Lifestyle", "hab_breadth", "Nest_type", # no trophic level
                        "Foraging",  "Gregariousness", "Mating_system")

mFD::traits.faxes.cor(sp_tr=traits_imp[morpho], 
                      sp_faxes_coord = sp_coords[ , c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9")], 
                      plot= TRUE)

mFD::traits.faxes.cor(sp_tr=traits_imp[,demog], 
                      sp_faxes_coord = sp_coords[ , c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9")], 
                      plot= TRUE)

mFD::traits.faxes.cor(sp_tr=traits_imp[,develop], 
                      sp_faxes_coord = sp_coords[ , c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9")], 
                      plot= TRUE)

mFD::traits.faxes.cor(sp_tr=traits_imp[,behav], 
                      sp_faxes_coord = sp_coords[ , c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9")], 
                      plot= TRUE)


# convex hulls
