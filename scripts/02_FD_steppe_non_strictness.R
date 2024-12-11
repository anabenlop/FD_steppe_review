# Functional Diversity analysis for steppe birds
# Date: 10/09/2024
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

# sp_data$Type <- ifelse(sp_data$Type == "Strict steppe birds", "Strict", "NonStrict")

sp_data$pop_size <- (sp_data$Population.size.EBBA2..min. + sp_data$Population.size.EBBA2..max.)/2

# Threat status in Europe --- Alaudala rufescens no category in Europe. From Cornell: "Decreases reported in Iberia, where species 
# is locally regarded as “near-threatened” (race apetzii), and in Canary Is, where considered “endangered” on Gran Canaria 
# (polatzeki) and “critical” on Tenerife (nominate). Considered Non-threatened for now. Although pop trends are Decreasing
# López-Jiménez, N., Editor. (2021). Libro Rojo de las Aves de España. Sociedad Española de Ornitología, Madrid.
# https://seo.org/wp-content/uploads/2022/09/LIbro-Rojo-web-3_01.pdf

sp_data$threat <- ifelse(sp_data$Status.protection.Europe.IUCN == "Vulnerable" | 
                           sp_data$Status.protection.Europe.IUCN == "Endangered", "Threatened", "Non-threatened")

sp_data[sp_data$binomial == "Alaudala rufescens",]$Population.trend.Europe.IUCN <- "Decreasing"


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
sp_matrix <- sp_data %>% group_by(threat, binomial) %>% 
  summarize(n = n()) %>%
  spread(key = binomial,value = n) %>% 
  replace(is.na(.), 0) %>% #replace NAs with 0
  tibble::column_to_rownames(var = "threat")

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
func_space <- mFD::quality.fspaces(sp_dist = func_dist, maxdim_pcoa = 10,
                                   deviation_weighting = c("absolute","squared"), fdist_scaling = c(TRUE,FALSE),
                                   fdendro= "average")

#check quality metrics 
func_space$quality_fspaces #check quality metrics
apply(func_space$quality_fspaces, 2, which.min) # PCoA 9D seems to be the best, but the gain is little after PCoA 4D

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
                          quality_metric  = "rmsd_scaled",
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

func_space$details_fspaces$pc_eigenvalues #get explained variance by each PC

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
                                      sp_faxes_coord = sp_coords[ , c("PC1", "PC2", "PC3", "PC4")], 
                                      plot = TRUE)

trait_faxes2 <- mFD::traits.faxes.cor(sp_tr=trait_data[,demog], 
                                      sp_faxes_coord = sp_coords[ , c("PC1", "PC2", "PC3", "PC4")], 
                                      plot = TRUE)

trait_faxes3 <- mFD::traits.faxes.cor(sp_tr=trait_data[,develop], 
                                      sp_faxes_coord = sp_coords[ , c("PC1", "PC2", "PC3", "PC4")], 
                                      plot = TRUE)

trait_faxes4 <- mFD::traits.faxes.cor(sp_tr=trait_data[,behav], 
                                      sp_faxes_coord = sp_coords[ , c("PC1", "PC2", "PC3", "PC4")], 
                                      plot = TRUE)

trait_faxes4a <- mFD::traits.faxes.cor(sp_tr=trait_data[,behav[1:5]], 
                                      sp_faxes_coord = sp_coords[ , c("PC1", "PC2", "PC3", "PC4")], 
                                      plot = TRUE)

trait_faxes4b <- mFD::traits.faxes.cor(sp_tr=trait_data[,behav[6:10]], 
                                       sp_faxes_coord = sp_coords[ , c("PC1", "PC2", "PC3", "PC4")], 
                                       plot = TRUE)


trait_faxes1$"tr_faxes_stat"[which(trait_faxes1$"tr_faxes_stat"$"p.value" < 0.05), ]

trait_faxes2$"tr_faxes_stat"[which(trait_faxes2$"tr_faxes_stat"$"p.value" < 0.05), ]

trait_faxes3$"tr_faxes_stat"[which(trait_faxes3$"tr_faxes_stat"$"p.value" < 0.05), ]

trait_faxes4$"tr_faxes_stat"[which(trait_faxes4$"tr_faxes_stat"$"p.value" < 0.05), ]

# Plotting 
trait_faxes1$"tr_faxes_plot"
trait_faxes2$"tr_faxes_plot"
trait_faxes3$"tr_faxes_plot"
trait_faxes4$"tr_faxes_plot"
trait_faxes4a$"tr_faxes_plot"
trait_faxes4b$"tr_faxes_plot"

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

# calculate indices for 4 axes
FD_indices <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_coords[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = as.matrix(sp_matrix), # here change matrix
  ind_vect         = c("fdis", "fmpd", "fnnd", "feve", "fric", "fdiv", "fori", 
                       "fspe", "fide"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

FD_df <- FD_indices$"functional_diversity_indices" #transform to df
FD_df$threat <- row.names(FD_df)

#Plot relationship between FD indices and type of steppe bird
ggplot(FD_df) + geom_point(aes(x = threat, y = fric))
ggplot(FD_df) + geom_point(aes(x = threat, y = fdiv))
ggplot(FD_df) + geom_point(aes(x = threat, y = fdis))
ggplot(FD_df) + geom_point(aes(x = threat, y = fspe))

# calculate indices for 2 axes
FD_indices_2axes <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_coords[ , c("PC1", "PC2")],
  asb_sp_w         = as.matrix(sp_matrix), # here change matrix
  ind_vect         = c("fdis", "fmpd", "fnnd", "feve", "fric", "fdiv", "fori", 
                       "fspe", "fide"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

FD_2ax_df <- FD_indices_2axes$"functional_diversity_indices" #transform to df
FD_2ax_df$Type <- row.names(FD_2ax_df)
FD_2ax_df
# 
# plots convex hull all, strict steppe and non-strict steppe species

plots_alpha <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = FD_indices,
  plot_asb_nm              = c("Non-threatened", "Threatened"),
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
plots_alpha$"fric"$PC3_PC4 #
# plots_alpha$"fdiv"$"patchwork"
# plots_alpha$"fspe"$"patchwork"
# plots_alpha$"fdis"$"patchwork"
# plots_alpha$"fori"$"patchwork"


# Plot trait space Mammola et al. 2022 ----
source("scripts/get_position.R")
myCol <- viridis::viridis(n=6, option="B") #Make Color palette 
myCol2 <- viridis::viridis(n=6, option="D") #Make Color palette 
pie(seq_along(myCol),myCol,col= myCol)
pie(seq_along(myCol2),myCol2,col= myCol2)

# ord <- cmdscale(func_dist, k = 4, eig = TRUE) MAMMOLA CALCULATION

# we get species coordinates from the functional dimensions calculated above
func_space$details_fspaces$pc_eigenvalues #get explained variance by each PC

coordinates <- data.frame(sp_coords[,1:4])
colnames(coordinates) <- c("PC1", "PC2", "PC3", "PC4")
coordinates$binomial <- rownames(coordinates)
coordinates <- inner_join(coordinates, sp_tax[,c("binomial", "Order")], by = "binomial")
coordinates <- inner_join(coordinates, sp_data[,c("binomial", "Type", "threat")], by = "binomial")

centroid <- coordinates %>% as_tibble() %>%
  # add_column(order = sp_tax$Order) %>%
  group_by(Order) %>%
  summarise(cen.1 = mean(PC1), cen.2 = mean(PC2))

centroid_threat12 <- coordinates %>% as_tibble() %>%
  # add_column(order = sp_tax$Order) %>%
  group_by(threat) %>%
  summarise(cen.1 = mean(PC1), cen.2 = mean(PC2))

centroid_threat34 <- coordinates %>% as_tibble() %>%
  # add_column(order = sp_tax$Order) %>%
  group_by(threat) %>%
  summarise(cen.3 = mean(PC3), cen.4 = mean(PC4))


fit <- vegan::envfit(ord = sp_coords, env = trait_data, w = NULL, na.rm = TRUE)

plot(sp_coords)
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
    "Tail length" = "Tail.Length",
    "Mass" = "Mass",
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
    "Development continuum" = "Development_continuum",
    "Monogomaus" = "Mating_systemmonogamous",
    "Polyandrous" = "Mating_systempolyandrous",
    "Polyginic" = "Mating_systempolyginic",
    "Monog-polyg." = "Mating_systemmonogamous-polyginic",
    "Seasonal gregarious" = "Gregariousnessseasonal",
    "All-year gregarious" = "Gregariousnessall-year round",
    "Solitary" = "Gregariousnessno",
    "Relative brain size" =  "relative_brain",
    "Rock" = "HabitatRock",
    "Grassland" = "HabitatGrassland",
    "Shrubland" = "HabitatShrubland",
    "Desert" = "HabitatDesert",
    "Human-made" = "HabitatHuman Modified",
    "Coast" = "HabitatCoastal",
    "Woodland" = "HabitatWoodland",
    "Wetland" = "HabitatWetland", 
    "Semi-open habitat" = "Habitat.Density2",  
    "Open habitat" = "Habitat.Density3",
    "Sedendary" = "Migration1",
    "Partial mig." = "Migration2",
    "Migratory" = "Migration3",
    "Omnivore" = "Trophic.NicheOmnivore",
    "Invertivore" = "Trophic.NicheInvertivore",
    "Granivore" = "Trophic.NicheGranivore",
    "Herbivore" = "Trophic.NicheHerbivore terrestrial",
    "Carnivore" = "Trophic.NicheVertivore",
    "Terrestrial" = "Primary.LifestyleTerrestrial",
    "Aerial" = "Primary.LifestyleAerial",
    "Insessorial" = "Primary.LifestyleInsessorial",
    "Generalist" = "Primary.LifestyleGeneralist",
    "Diurnal" =  "activitydiurnal",
    "Nocturnal" =  "activitynocturnal",
    "Ground-nesting" = "Nest_typeExposedGround",
    "Cavity-nesting" = "Nest_typeCavity",
    "Elevated-nesting" = "Nest_typeExposedElevated",
    "Generalist foraging" = "ForagingForagingGeneralist",
    "Ground foraging" = "ForagingGroundForaging",
    "Arboreal gleaning" = "ForagingArborealGleaning",
    "Aerial sallying" = "ForagingAerialSallying",
    "Aerial screening" = "ForagingAerialScreening",
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
  scale_fill_gradientn(colours = rev(myCol2)) +
  theme(
    panel.background = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 1,
      linetype = "solid"
    ),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)) +
  labs(x = "PCoA 1 (48%)", y = "PCoA 2 (21%)") +
  ylim(-.46, .46) + xlim(-.46, .46) +
  coord_fixed()

plot_space_trait

plot_space_trait_34 <-  ggplot(data = coordinates, aes(PC3, PC4)) +
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
  scale_fill_gradientn(colours = rev(myCol2)) +
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
  labs(x = "PCoA 3 (7.1%)", y = "PCoA 4 (6.4%)") +
  ylim(-.40, .40) + xlim(-.45, .45) +
  coord_fixed()

plot_space_trait_34

# plot functional space for threatened and non-threatened bird species ------
plot_threat <-
  ggplot(coordinates, aes(PC1, PC2)) +
  stat_density_2d(
    aes(fill = after_stat(level)),
    geom = "polygon",
    colour = NA,
    alpha = .5,
    h = .25
  ) +
  geom_hline(aes(yintercept = 0), linetype = 3, colour = "gray70") +
  geom_vline(aes(xintercept = 0), linetype = 3, colour = "gray70") +
  geom_point(
    data = centroid_threat12,
    aes(x = cen.1, y = cen.2),
    shape = 21,
    fill = "white",
    size = 2.5
  ) +
  geom_point(
    data = centroid_threat,
    aes(x = cen.1, y = cen.2),
    shape = 19,
    colour = "black",
    size = 1
  ) +
  scale_fill_gradientn(colours = rev(myCol2)) +
  ggrepel::geom_text_repel(data = centroid_threat12,
                           aes(x = cen.1,
                               y = cen.2,
                               label = threat)) +
  theme(
    panel.background = element_rect(
      fill = NA,
      colour = "black",
      linewidth  = 1,
      linetype = "solid"
    ),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)) +
  labs(x = "PCoA 1 (48%)", y = "PCoA 2 (21%)") +
  ylim(-.46, .46) + xlim(-.46, .46) +
  coord_fixed()

plot_threat

plot_threat2 <-
  ggplot(coordinates, aes(PC1, PC2)) +
  stat_density_2d(
    aes(fill = after_stat(level)),
    geom = "polygon",
    colour = NA,
    alpha = .5,
    h = .25
  ) +
  geom_hline(aes(yintercept = 0), linetype = 3, colour = "gray70") +
  geom_vline(aes(xintercept = 0), linetype = 3, colour = "gray70") +
  geom_point(
    data = centroid_threat12,
    aes(x = cen.1, y = cen.2),
    shape = 21,
    fill = "white",
    size = 2.5
  ) +
  geom_point(
    data = centroid_threat12,
    aes(x = cen.1, y = cen.2),
    shape = 19,
    colour = "black",
    size = 1
  ) +
  facet_wrap(vars(threat), ncol = 2) +
  scale_fill_gradientn(colours = rev(myCol2)) +
  theme(
    panel.background = element_rect(
      fill = NA,
      colour = "black",
      linewidth  = 1,
      linetype = "solid"
    ),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)) +
  labs(x = "PCoA 1 (48%)", y = "PCoA 2 (21%)")+
  ylim(-.46, .46) + xlim(-.46, .46) +
  coord_fixed()

plot_threat2

plot_threat3 <-
  ggplot(coordinates, aes(PC3, PC4)) +
  stat_density_2d(
    aes(fill = after_stat(level)),
    geom = "polygon",
    colour = NA,
    alpha = .5,
    h = .25
  ) +
  geom_hline(aes(yintercept = 0), linetype = 3, colour = "gray70") +
  geom_vline(aes(xintercept = 0), linetype = 3, colour = "gray70") +
  geom_point(
    data = centroid_threat34,
    aes(x = cen.3, y = cen.4),
    shape = 21,
    fill = "white",
    size = 2.5
  ) +
  geom_point(
    data = centroid_threat34,
    aes(x = cen.3, y = cen.4),
    shape = 19,
    colour = "black",
    size = 1
  ) +
  scale_fill_gradientn(colours = rev(myCol2)) +
  ggrepel::geom_text_repel(data = centroid_threat34,
                           aes(x = cen.3,
                               y = cen.4,
                               label = threat)) +
  theme(
    panel.background = element_rect(
      fill = NA,
      colour = "black",
      linewidth  = 1,
      linetype = "solid"
    ),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)) +
  labs(x = "PCoA 3 (7.1%)", y = "PCoA 4 (6.4%)") +
  ylim(-.46, .46) + xlim(-.46, .46) +
  coord_fixed()

plot_threat3

plot_threat4 <-
  ggplot(coordinates, aes(PC3, PC4)) +
  stat_density_2d(
    aes(fill = after_stat(level)),
    geom = "polygon",
    colour = NA,
    alpha = .5,
    h = .25
  ) +
  geom_hline(aes(yintercept = 0), linetype = 3, colour = "gray70") +
  geom_vline(aes(xintercept = 0), linetype = 3, colour = "gray70") +
  geom_point(
    data = centroid_threat34,
    aes(x = cen.3, y = cen.4),
    shape = 21,
    fill = "white",
    size = 2.5
  ) +
  geom_point(
    data = centroid_threat34,
    aes(x = cen.3, y = cen.4),
    shape = 19,
    colour = "black",
    size = 1
  ) +
  facet_wrap(vars(threat), ncol = 2) +
  scale_fill_gradientn(colours = rev(myCol2)) +
  theme(
    panel.background = element_rect(
      fill = NA,
      colour = "black",
      linewidth  = 1,
      linetype = "solid"
    ),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)) +
  labs(x = "PCoA 3 (7.1%)", y = "PCoA 4 (6.4%)") +
  ylim(-.46, .46) + xlim(-.46, .46) +
  coord_fixed()

plot_threat4

# plot functional space type of steppe bird with silhouttes ----
library(png)
library(grid)

# LOAD SILLOUETES

#Otididae
# img <- readPNG("Silhouette/Otididae.png")
# pholcidae <- rasterGrob(img, interpolate = TRUE)
# pholcidae$width <- unit(.6, "npc")
# pholcidae$height <- unit(.6, "npc")
# 
# #Pteroclidae
# img <- readPNG("Silhouette/Pteroclidae.png")
# dysderidae <- rasterGrob(img, interpolate = TRUE)
# dysderidae$width <- unit(.6, "npc")
# dysderidae$height <- unit(.6, "npc")
# 
# #Nesticidae
# img <- readPNG("Silhouette/Nesticidae.png")
# nesticidae <- rasterGrob(img, interpolate = TRUE)
# nesticidae$width <- unit(.6, "npc")
# nesticidae$height <- unit(.6, "npc")

# plot traits onto ordination space
plot_traits <-   ggplot(coordinates, aes(PC1, PC2)) +
  stat_density_2d(
    aes(fill = after_stat(level)),
    geom = "polygon",
    colour = NA,
    alpha = .5,
    h = .25
  ) +
  geom_hline(aes(yintercept = 0), linetype = 3, colour = "gray70") +
  geom_vline(aes(xintercept = 0), linetype = 3, colour = "gray70") +
  geom_point(
    data = trait_pos,
    aes(x = PC1, y = PC2),
    shape = 21,
    fill = "white",
    size = 2.5
  ) +
  geom_point(
    data = trait_pos,
    aes(x = PC1, y = PC2),
    shape = 19,
    colour = "black",
    size = 1
  ) +
  scale_fill_gradientn(colours = rev(myCol2)) +
  ggrepel::geom_text_repel(data = trait_pos, aes(x = PC1, y = PC2, label = Trait), max.overlaps = 20) +
  theme(
    panel.background = element_rect(
      fill = NA,
      colour = "black",
      linewidth  = 1,
      linetype = "solid"
    ),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)) +
  labs(x = "PCoA 1 (48%)", y = "PCoA 2 (21%)") +
  ylim(-.46, .46) + xlim(-.46, .46) +
  coord_fixed()

plot_traits

## plot family silhouettes

plot_traits2 <-   ggplot(coordinates, aes(PC1, PC2)) +
  stat_density_2d(
    aes(fill = ..level..),
    geom = "polygon",
    colour = NA,
    alpha = .5,
    h = .25
  ) +
  geom_hline(aes(yintercept = 0), linetype = 3, colour = "gray70") +
  geom_vline(aes(xintercept = 0), linetype = 3, colour = "gray70") +
  geom_point(
    data = trait_pos,
    aes(x = PC1, y = PC2),
    shape = 21,
    fill = "white",
    size = 2.5
  ) +
  geom_point(
    data = trait_pos,
    aes(x = PC1, y = PC2),
    shape = 19,
    colour = "black",
    size = 1
  ) +
  scale_fill_gradientn(colours = rev(myCol2)) +
  ggrepel::geom_text_repel(data = trait_pos, aes(x = PC1, y = PC2, label = Trait), max.overlaps = 15) +
  theme(
    panel.background = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 1,
      linetype = "solid"
    ),
    panel.grid = element_blank(),
    legend.position = "none"
  )  +
  labs(x = "PCoA 1 (43%)", y = "PCoA 2 (35%)") +
  ylim(-.30, .30) + xlim(-.50, .50) +
  coord_fixed() +
  facet_wrap(vars(Type), ncol = 2) 

plot_traits2

plot_traits3 <-   ggplot(coordinates, aes(PC1, PC2)) +
  stat_density_2d(
    aes(fill = ..level..),
    geom = "polygon",
    colour = NA,
    alpha = .5,
    h = .25
  ) +
  geom_hline(aes(yintercept = 0), linetype = 3, colour = "gray70") +
  geom_vline(aes(xintercept = 0), linetype = 3, colour = "gray70") +
  # geom_point(
  #   data = trait_pos,
  #   aes(x = PC1, y = PC2),
  #   shape = 21,
  #   fill = "white",
  #   size = 2.5
  # ) +
  # geom_point(
  #   data = trait_pos,
  #   aes(x = PC1, y = PC2),
  #   shape = 19,
  #   colour = "black",
  #   size = 1
  # ) +
  scale_fill_gradientn(colours = rev(myCol2)) +
  # ggrepel::geom_text_repel(data = trait_pos, aes(x = PC1, y = PC2, label = Trait), max.overlaps = 15) +
  theme(
    panel.background = element_rect(
      fill = NA,
      colour = "black",
      size = 1,
      linetype = "solid"
    ),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)) +
  labs(x = "PCoA 1 (48%)", y = "PCoA 2 (21%)") +
  ylim(-.50, .50) + xlim(-.50, .50) +
  coord_fixed() +
  facet_wrap(vars(Type, threat), ncol = 2, nrow = 2) 

plot_traits3

## check relationship between PCs and poptrends----
head(sp_data)
class(sp_coords)

sp_PC <- data.frame(binomial = row.names(sp_coords), sp_coords[,1:2])

sp_trends <- inner_join(sp_data, sp_PC, by = "binomial")

hist(sp_trends$Change.index)

#check with change index
ggplot(sp_trends) + geom_point(aes(PC1, Change.index)) +
                    geom_smooth(aes(PC1, Change.index), method = "lm")

ggplot(sp_trends) + geom_point(aes(PC2, Change.index)) +
  geom_smooth(aes(PC2, Change.index), method = "lm")

summary(lm(Change.index ~ PC1, data = sp_trends))
summary(lm(Change.index ~ PC2, data = sp_trends))

#check with EBBA2 trends
sp_trends$Pop.trend.EBBA2_recalc <- ifelse(is.na(sp_trends$Change.index), "Unknown",
                                           ifelse(sp_trends$Change.index > 0, "Increase", "Loss"))
                                           
table(sp_trends$Pop.trend.EBBA2_recalc)
ggplot(sp_trends) + geom_boxplot(aes(Pop.trend.EBBA2_recalc, PC1))

ggplot(sp_trends) + geom_boxplot(aes(Pop.trend.EBBA2_recalc, PC2))


#check with Population.trend.Europe.IUCN
table(sp_trends$Population.trend.Europe.IUCN)
ggplot(sp_trends) + geom_boxplot(aes(Population.trend.Europe.IUCN, PC1))
 
ggplot(sp_trends) + geom_boxplot(aes(Population.trend.Europe.IUCN, PC2))
