# Functional Diversity analysis for steppe birds
# Date: 05/08/2024
# Author: Ana Benítez


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

# Load steppe birds trait data (clean):

#Import data
trait_data <- read.csv("clean_data/traits_clean.csv", stringsAsFactors = T)
trait_cat <- read.csv("clean_data/traits_cat_clean.csv")
sp_data <- read.csv("data/Species_trends.csv")
sp_tax <- read.csv("data/TraitsSteppebirdsFull.csv")
sp_tax <- sp_tax[,c("binomial", "Order", "Family")]

sp_data$Type <- ifelse(sp_data$Type == "Strict steppe birds", "Strict", "NonStrict")

sp_data$pop_size <- (sp_data$Population.size.EBBA2..min. + sp_data$Population.size.EBBA2..max.)/2

# change character to factor for categorical variables
trait_data$Habitat <- as.factor(trait_data$Habitat)
trait_data$Habitat.Density <- as.factor(trait_data$Habitat.Density)
trait_data$Migration <- as.factor(trait_data$Migration)
trait_data$Trophic.Niche <- as.factor(trait_data$Trophic.Niche)
trait_data$Primary.Lifestyle <- as.factor(trait_data$Primary.Lifestyle)
trait_data$activity <- as.factor(trait_data$activity)
trait_data$Nest_type <- as.factor(trait_data$Nest_type)
trait_data$Foraging <- as.factor(trait_data$Foraging)
trait_data$Gregariousness <- as.factor(trait_data$Gregariousness)
trait_data$Mating_system <- as.factor(trait_data$Mating_system)

# Create matrix depicting strict and non-strict steppe bird species, with presence/absence and with abundance data
sp_matrix <- sp_data %>% group_by(Type, binomial) %>% 
  summarize(n = n()) %>%
  spread(key = binomial,value = n) %>% 
  replace(is.na(.), 0) %>% #replace NAs with 0
  tibble::column_to_rownames(var = "Type")

sp_ab_matrix <- sp_data %>% group_by(Type, binomial) %>% 
  summarize(pop_size = pop_size) %>%
  spread(key = binomial,value = pop_size)%>% #population size (total number individuals)
  replace(is.na(.), 0) %>% #replace NAs with 0
  tibble::column_to_rownames(var = "Type")

#Setting species names as row names
rownames(trait_data) <- trait_data$binomial
trait_data$binomial <- NULL 

# Species traits summary:
species_traits_summ <- mFD::sp.tr.summary(
  tr_cat     = trait_cat,   
  sp_tr      = trait_data, 
  stop_if_NA = TRUE)

# Functional distance analyses ####
#calculate functional distance matrix 
func_dist <- mFD::funct.dist(sp_tr = trait_data, tr_cat = trait_cat, metric="gower")

# build up functional spaces: calculate with different weighting of deviation and
# with or without scaling functional distance.
func_space <- mFD::quality.fspaces(sp_dist= func_dist, maxdim_pcoa = 8,
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
                          fspaces_plot = c("tree_average", "pcoa_2d", "pcoa_3d",
                                           "pcoa_4d", "pcoa_5d", "pcoa_6d",
                                           "pcoa_7d", "pcoa_8d"),
                          gradient_deviation = c(neg = "#030091", nul = "#66B79C", pos = "#C6FFB7"),
                          gradient_deviation_quality = c(low = "#B7E3FF", high = "#2800B2"),
                          x_lab= "Trait-based distance")

# From Magneville et al. 2022
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
           "birth_or_hatching_weight_g", "annual_survival", "PopulationDensity_ind_km2") # no longevity_y or Range.Size

develop <-  c("relative_brain", "eye_axial_length", "eye_tranverse_diameter", "Development_continuum", "activity")

behav <-  c("Habitat", "Habitat.Density","Migration", "Trophic.Niche", "Primary.Lifestyle", "hab_breadth", "Nest_type", # no trophic level
            "Foraging",  "Gregariousness", "Mating_system")

trait_faxes1 <- mFD::traits.faxes.cor(sp_tr=trait_data[morpho], 
                                      sp_faxes_coord = sp_coords[ , c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")], 
                                      plot= TRUE)

trait_faxes2 <- mFD::traits.faxes.cor(sp_tr=trait_data[,demog], 
                                      sp_faxes_coord = sp_coords[ , c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")], 
                                      plot= TRUE)

trait_faxes3 <- mFD::traits.faxes.cor(sp_tr=trait_data[,develop], 
                                      sp_faxes_coord = sp_coords[ , c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")], 
                                      plot= TRUE)

trait_faxes4 <- mFD::traits.faxes.cor(sp_tr=trait_data[,behav], 
                                      sp_faxes_coord = sp_coords[ , c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")], 
                                      plot= TRUE)


trait_faxes1$"tr_faxes_stat"[which(trait_faxes1$"tr_faxes_stat"$"p.value" < 0.05), ]

trait_faxes2$"tr_faxes_stat"[which(trait_faxes2$"tr_faxes_stat"$"p.value" < 0.05), ]

trait_faxes3$"tr_faxes_stat"[which(trait_faxes3$"tr_faxes_stat"$"p.value" < 0.05), ]

trait_faxes4$"tr_faxes_stat"[which(trait_faxes4$"tr_faxes_stat"$"p.value" < 0.05), ]

# Plotting 
trait_faxes1$"tr_faxes_plot"
trait_faxes2$"tr_faxes_plot"
trait_faxes3$"tr_faxes_plot"
trait_faxes4$"tr_faxes_plot"

# convex hulls

big_plot <- mFD::funct.space.plot(
  sp_faxes_coord  = sp_coords[ ,c("PC1", "PC2", "PC3", "PC4")],
  faxes           = c("PC1", "PC2", "PC3", "PC4"),
  name_file       = NULL,
  faxes_nm        = NULL,
  range_faxes     = c(NA, NA),
  color_bg        = "grey95",
  color_pool      = "darkgreen",
  fill_pool       = "white",
  shape_pool      = 21,
  size_pool       = 1,
  plot_ch         = TRUE,
  color_ch        = "black",
  fill_ch         = "white",
  alpha_ch        = 0.5,
  plot_vertices   = TRUE,
  color_vert      = "blueviolet",
  fill_vert       = "blueviolet",
  shape_vert      = 23,
  size_vert       = 1,
  plot_sp_nm      = NULL,
  nm_size         = 3,
  nm_color        = "black",
  nm_fontface     = "plain",
  check_input     = TRUE)

big_plot

# Alpha FD for strict and non-strict steppe birds ##############################

FD_indices <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_coords[ , c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")],
  asb_sp_w         = as.matrix(sp_matrix), # here change matrix
  ind_vect         = c("fdis", "fmpd", "fnnd", "feve", "fric", "fdiv", "fori", 
                       "fspe", "fide"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

FD_indices_ab <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_coords[ , c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")],
  asb_sp_w         = as.matrix(sp_ab_matrix), # here change matrix
  ind_vect         = c("fdis", "fmpd", "fnnd", "feve", "fric", "fdiv", "fori", 
                       "fspe", "fide"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

FD_df <- FD_indices$"functional_diversity_indices" #transform to df
FD_df$Type <- row.names(FD_df)

FD_ab_df <- FD_indices_ab$"functional_diversity_indices" #transform to df
FD_ab_df$Type <- row.names(FD_ab_df)


#Plot relationship between FD indices and type of steppe bird
ggplot(FD_df) + geom_point(aes(x = Type, y = fdiv))
ggplot(FD_df) + geom_point(aes(x = Type, y = fdis))
ggplot(FD_df) + geom_point(aes(x = Type, y = fspe))

ggplot(FD_ab_df) + geom_point(aes(x = Type, y = fdiv))
ggplot(FD_ab_df) + geom_point(aes(x = Type, y = fdis))
ggplot(FD_ab_df) + geom_point(aes(x = Type, y = fspe))

# plots convex hull all, strict steppe and non-strict steppe species

plots_alpha <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = FD_indices,
  plot_asb_nm              = c("Strict", "NonStrict"),
  ind_nm                   = c("fdis", "fide", "fnnd", "feve", "fric", 
                               "fdiv", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "white",
  shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
  size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
  color_sp                 = c(pool = "grey50", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  color_vert               = c(pool = "grey50", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  fill_sp                  = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  fill_vert                = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  color_ch                 = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  fill_ch                  = c(pool = "grey90", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fdiv      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE) 

plots_alpha$"fric"$"patchwork" # convex hull
plots_alpha$"fric"$PC1_PC2 #
# plots_alpha$"fdiv"$"patchwork"
# plots_alpha$"fspe"$"patchwork"
# plots_alpha$"fdis"$"patchwork"
# plots_alpha$"fori"$"patchwork"


### Plot trait space Mammola et al. 2022 ----
source("scripts/get_position.R")
myCol<-viridis::viridis(n=6, option="B") #Make Color palette 
pie(seq_along(myCol),myCol,col= myCol)

ord <- cmdscale(func_dist)
coordinates <- data.frame(ord)
colnames(coordinates) <- c("PC1", "PC2")

centroid <- coordinates %>% as_tibble() %>%
  add_column(order = sp_tax$Order) %>%
  group_by(order) %>%
  summarise(cen.1 = mean(PC1), cen.2 = mean(PC2))

fit <- vegan::envfit(ord = ord, env = trait_data, w = NULL, na.rm = TRUE)

plot(ord)
trait_position <- get_position(fit, add = TRUE)

trait_pos <- data.frame(trait_position) %>%
  rownames_to_column("Trait") %>%
  mutate_at(
    "Trait",
    .funs = forcats::fct_recode,
    "Beak culmen" = "Beak.Length_Culmen",
    "Beak nares" = "Beak.Length_Nares",
    "Beak width" = "Beak.Width",
    "Beak depth" = "Beak.Depth",
    "Tarsus length" = "Tarsus.Length",
    "Wing length" = "Wing.Length",
    "Kipps dist." = "Kipps.Distance",
    "Hand wing index" = "Hand.Wing.Index",
    "Tail length" = "Tail.length",
    "Mass" = "Mass",
    "Habitat density" = "Habitat.Density",
    "Trophic niche" = "Trophic.Niche",
    "Lifestyle" = "Primary.Lifestyle",
    "Range size" = "Range.Size",
    "Habitat breadth" = "hab_breadth",
    "Length" = "adult_length_cm",
    "Fem. maturity" = "female_maturity_d",
    "Clutch size" = "clutch_size_n",
    "Clutches per year" = "clutches_per_y",
    "Max. longevity" =  "maximum_longevity_y",
    "Egg mass" = "egg_mass_g",
    "Incubation period" = "incubation_d",
    "Fledging age" = "fledging_age_d",
    "Hatching weight" = "birth_or_hatching_weight_g",
    "Survival" =  "annual_survival",
    "Pop. density" = "PopulationDensity_ind_km2",
    "Eye axial length" = "eye_axial_length",
    "Eye tranverse diameter" = "eye_tranverse_diameter",
    "Activity" =  "activity",
    "Nest type" = "Nest_type",
    "Development continuum" = "Development_continuum",
    "Mating system" = "Mating_system",
    "Relative brain size" =  "relative_brain"

  ) 
# %>%
#   filter(
#     Trait != c(
#       "AME_typeAbsent_Adaptation",
#       "AME_typeAbsent_Ontology",
#       "AME_typePresent"
#     )
#   )


theme_set(theme_minimal()) #Setting the theme

plot_space_trait <-  ggplot(data = coordinates, aes(PC1, PC2)) +
  stat_density_2d(
    aes(fill = after_stat(level)),
    geom = "polygon",
    colour = NA,
    alpha = .5,
    h = .25
  ) +
  geom_hline(aes(yintercept = 0), linetype = 3, colour = "gray70") +
  geom_vline(aes(xintercept = 0), linetype = 3, colour = "gray70") +
  geom_point(shape = 19,
             size = .5,
             colour = "black") +
  scale_fill_gradientn(colours = rev(myCol)) +
  theme(
    panel.background = element_rect(
      fill = NA,
      colour = "black",
      linewidth  = 1,
      linetype = "solid"
    ),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "PCoA 1 (48%)", y = "PCoA 2 (29%)") +
  ylim(-.46, .36) + xlim(-.46, .45) +
  coord_fixed()

plot_space_trait
