# Functional Diversity analysis for steppe birds
# Date: 10/09/2024
# Author: Ana Benítez

# Here we remove species that are not really steppe-bird species according to our criteria
# Charadrius morinellus, Gelochelidon nilotica, Oenanthe leucura, Falco tinnunculus, Bucanetes githagineus 


# load libraries
library(dplyr)
library(mFD)
library(tidyverse)
library(corrplot)
library(caret)
library(missForest)
library(tibble)
library(FD)
library(ggplot2)
library(picante)
library(tictoc)

# clean environment
rm(list = ls())

# Load steppe birds trait data (clean):

#Import data
trait_data <- read.csv("clean_data/traits_clean.csv", stringsAsFactors = T)
trait_cat <- read.csv("clean_data/traits_cat_clean.csv")
sp_data <- read.csv("data/Species_trends.csv")
sp_tax <- read.csv("data/TraitsSteppebirdsFull.csv")
sp_tax <- sp_tax[,c("binomial", "Order", "Family")]
# sp_tax <- sp_tax[1:41,]

# remove non steppe bird species

# Charadrius morinellus, Gelochelidon nilotica, Oenanthe leucura, Saxicola dacotiae

trait_data <- trait_data %>% 
              filter(!binomial %in% c("Eudromias morinellus", "Gelochelidon nilotica","Oenanthe leucura", 
                                      "Falco tinnunculus", "Bucanetes githagineus"))

sp_tax <- sp_tax %>% 
          filter(!binomial %in% c("Eudromias morinellus", "Gelochelidon nilotica","Oenanthe leucura", 
                          "Falco tinnunculus", "Bucanetes githagineus"))

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
trait_data$Nest_type <- as.factor(trait_data$Nest_type)
trait_data$Foraging <- as.factor(trait_data$Foraging)
trait_data$Gregariousness <- as.factor(trait_data$Gregariousness)

# Create matrix depicting threatened and non-threatened steppe bird species, with presence/absence and with abundance data
sp_matrix <- sp_data %>% group_by(threat, binomial) %>% 
  summarize(n = n()) %>%
  spread(key = binomial,value = n) %>% 
  replace(is.na(.), 0) %>% #replace NAs with 0
  tibble::column_to_rownames(var = "threat")

# All_sp <- rep(1, length(sp_matrix))
# sp_matrix <- rbind(sp_matrix,All_sp)
# rownames(sp_matrix)[3] <- "All"

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
func_space <- mFD::quality.fspaces(sp_dist = func_dist, maxdim_pcoa = 4,
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
                                           "pcoa_4d"),
                          gradient_deviation = c(neg = "#030091", nul = "#66B79C", pos = "#C6FFB7"),
                          gradient_deviation_quality = c(low = "#B7E3FF", high = "#2800B2"),
                          x_lab= "Trait-based distance")


# Plot relation between functional axes and traits
sp_coords <- func_space$details_fspaces$sp_pc_coord

func_space$details_fspaces$pc_eigenvalues #get explained variance by each PC
 # PC1 50.0%; PC2 22.8% --- # PC1 50%; PC2 23% Total: 73%


# group traits by type
morpho <- c("Beak.Length_Culmen", "Beak.Length_Nares","Beak.Width", "Beak.Depth", 
            "Tarsus.Length", "Wing.Length", "Kipps.Distance", "Hand.Wing.Index", "Mass", # not Tail.length
            "adult_length_cm")

demog <- c("female_maturity_d", "clutch_size_n", "clutches_per_y", "maximum_longevity_y", "egg_mass_g", "incubation_d", "fledging_age_d",  
           "birth_or_hatching_weight_g", "annual_survival", "PopulationDensity_ind_km2") # no longevity_y or Range.Size

develop <-  c("relative_brain", "eye_axial_length", "eye_tranverse_diameter", "Development_continuum")

behav <-  c("Habitat", "Habitat.Density","Migration", "Trophic.Niche", "Primary.Lifestyle", "hab_breadth", "Nest_type", # no trophic level
            "Foraging",  "Gregariousness")

trait_faxes1 <- mFD::traits.faxes.cor(sp_tr=trait_data[morpho], 
                                      sp_faxes_coord = sp_coords[ , c("PC1", "PC2")], 
                                      plot = TRUE)

trait_faxes2 <- mFD::traits.faxes.cor(sp_tr=trait_data[,demog], 
                                      sp_faxes_coord = sp_coords[ , c("PC1", "PC2")],  
                                      plot = TRUE)

trait_faxes3 <- mFD::traits.faxes.cor(sp_tr=trait_data[,develop], 
                                      sp_faxes_coord = sp_coords[ , c("PC1", "PC2")], 
                                      plot = TRUE)

trait_faxes4 <- mFD::traits.faxes.cor(sp_tr=trait_data[,behav], 
                                      sp_faxes_coord = sp_coords[ , c("PC1", "PC2")], 
                                      plot = TRUE)

trait_faxes4a <- mFD::traits.faxes.cor(sp_tr=trait_data[,behav[1:5]], 
                                       sp_faxes_coord = sp_coords[ , c("PC1", "PC2")], 
                                       plot = TRUE)

trait_faxes4b <- mFD::traits.faxes.cor(sp_tr=trait_data[,behav[6:9]], 
                                       sp_faxes_coord = sp_coords[ , c("PC1", "PC2")], 
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

# Alpha FD for strict and non-strict steppe birds ##############################
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

fric_mFD <- data.frame(Type = FD_2ax_df$Type, sp_ric = FD_2ax_df$sp_richn, fric = FD_2ax_df$fric)

# we might need to control for uneven number of species with null models - calculate SESFD ----
# 

#########################################################################
#####################CALCULATING FD (FRic) WITH FD PACKAGE ############# 
#########################################################################
trait_data_ord <- trait_data[order(row.names(trait_data)), ] # order trait data by species name

FDindexes_FD <- dbFD(
  trait_data_ord, as.matrix(sp_matrix), w.abun = FALSE, 
  corr ="sqrt", calc.FRic=TRUE, m = 2, 
  stand.FRic = TRUE, scale.RaoQ = FALSE, calc.FGR = FALSE,
  calc.CWM = FALSE, print.pco = TRUE
)

# FDindexes_FD$x.values
# FD_package  <- as.data.frame(FDindexes_FD)

unique(length(FDindexes_FD$FRic))

# # Converting row names into a column:
# FD_package_0 <- rownames_to_column(FD_package, "CUADRICULA")

######################################################
###Are FRic values obtained with mFD and FD different?###
############################################################


fric_mFD
#             Type sp_ric      fric
# 1 Non-threatened     30 0.9573449
# 2     Threatened      7 0.4546568

FDindexes_FD$FRic 
# Non-threatened     Threatened 
# 0.9639259      0.5615138 


######################################################################
###  Null model permutations to calculate SESFRic#####################
######################################################################
# SES values can be positive or negative, and their statistical significance at P < 0.05 is assumed when the 
# SES value falls outside the range of −1.96 to 1.96, assuming a normal distribution of deviations9,20. 
# On one hand, a SES value greater than 1.96 indicates that there is functional over-dispersion due to a limiting 
# similarity mechanism. On the other hand, SES values lower than −1.96 indicate that there is functional 
# under-dispersion, which means that an environmental filtering mechanism occurs. When SES values fall within 
# the range of −1.96 to 1.96, mechanistic processes cannot be proved, and then stochastic processes are assumed 
# as the observed values cannot be distinguished from random values.

# Positive SES values indicate that the considered variable (i.e., FD) is higher than expected considering 
# the number of species present in a community, and vice versa
# Observed values cannot be distinguished from random values.

# mFD Package
set.seed(123)
null.model.type <- "independentswap" # available methods: "independentswap", "frequency", "richness", "trialswap"
nperm <- 999
fric <- matrix(NaN, nrow = nperm, ncol = length(fric_mFD$fric))
abund.null <- vector(mode="list", nperm)

tic("Run 999 permutations to calculate SESFRic")
for(i in seq(nperm)){ # loop for each permutation
  # Create permutated abundance matrix
  abund.i <- randomizeMatrix(as.matrix(sp_matrix), null.model = null.model.type)
  traits.i <- trait_data_ord
  
  abund.null[[i]] <- abund.i
  
  # Record results
  res.i <- mFD::alpha.fd.multidim(
    sp_faxes_coord   = sp_coords[ , c("PC1", "PC2")],
    asb_sp_w         = abund.i, # here change matrix
    ind_vect         = c("fric"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  
  # Save result
  fric[i,] <- res.i$functional_diversity_indices$fric
  
  rm(abund.i); rm(res.i)
  print(paste0(round(i/nperm*100), "% of permutations completed"))
}
toc() # 6.3 sec elapsed

### Summary stats
SESFRic_mfd = (fric_mFD$fric - apply(fric, 2, mean)) / apply(fric, 2, sd) # standardized effect size
SESFRic_mfd
# Non-threatened     Threatened 
#  0.4599247       0.2677439

# check if they are outside expected distribution
q_fric <- NaN*fric_mFD$fric
for(i in seq(q_fric)){
  q_fric[i] <- sum(fric_mFD$fric[i] > fric[,i]) / length(fric[,i])
}
sigfric <- q_fric < 0.05 | q_fric > 0.95 # test if outside distribution
sigfric
# Non-threatened     Threatened 
# FALSE          FALSE 

# Observed values cannot be distinguished from random values.


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
coordinates <- inner_join(coordinates, sp_tax[,c("binomial", "Order", "Family")], by = "binomial")
coordinates <- inner_join(coordinates, sp_data[,c("binomial", "threat")], by = "binomial")

centroid <- coordinates %>% as_tibble() %>%
  # add_column(order = sp_tax$Order) %>%
  group_by(Order) %>%
  summarise(cen.1 = mean(PC1), cen.2 = mean(PC2))

centroid_fam <- coordinates %>% as_tibble() %>%
  # add_column(order = sp_tax$Order) %>%
  group_by(Family) %>%
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
    # "Monogomaus" = "Mating_systemmonogamous",
    # "Polyandrous" = "Mating_systempolyandrous",
    # "Polyginic" = "Mating_systempolyginic",
    # "Monog-polyg." = "Mating_systemmonogamous-polyginic",
    "Seasonal gregarious" = "Gregariousnessseasonal",
    "All-year gregarious" = "Gregariousnessall-year round",
    "Solitary" = "Gregariousnessno",
    "Relative brain size" =  "relative_brain",
    "Rock" = "HabitatRock",
    "Grassland" = "HabitatGrassland",
    "Shrubland" = "HabitatShrubland",
    "Desert" = "HabitatDesert",
    "Human-made" = "HabitatHuman Modified",
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
    # "Diurnal" =  "activitydiurnal",
    # "Nocturnal" =  "activitynocturnal",
    "Ground-nesting" = "Nest_typeExposedGround",
    "Cavity-nesting" = "Nest_typeCavity",
    "Elevated-nesting" = "Nest_typeExposedElevated",
    "Generalist foraging" = "ForagingForagingGeneralist",
    "Ground foraging" = "ForagingGroundForaging",
    "Arboreal gleaning" = "ForagingArborealGleaning",
    "Aerial sallying" = "ForagingAerialSallying",
    "Aerial screening" = "ForagingAerialScreening"
    ) 

trait_pos_short <- trait_pos %>% filter(
    Trait %in%
      c("Wing length", 
      "Length",
      "Tarsus length",
      "Hatching weight",
      "Max. longevity",
      "Fledging age",
      "Semi-open habitat",
      "Open habitat",
      "Terrestrial",
      "Ground-nesting",
      "Ground foraging",
      "Relative brain size",
      "Elevated-nesting",
      "Carnivore"
    )
)
   
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
  labs(x = "PCoA 1 (50%)", y = "PCoA 2 (23%)") +
  ylim(-.48, .48) + xlim(-.48, .48) +
  coord_fixed()

plot_space_trait

ggsave(
  plot = plot_space_trait,
  filename = "output/funcspace.pdf",
  width = 14,
  height = 10,
  units = "in")

ggsave(
  plot = plot_space_trait,
  filename = "output/funcspace.png",
  width = 7,
  height = 5,
  units = "in")

# plot traits onto ordination space
plot_trait_pos <-   ggplot(coordinates, aes(PC1, PC2)) +
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
    data = trait_pos_short,
    aes(x = PC1, y = PC2),
    shape = 21,
    fill = "white",
    size = 2.5
  ) +
  geom_point(
    data = trait_pos_short,
    aes(x = PC1, y = PC2),
    shape = 19,
    colour = "black",
    size = 1
  ) +
  scale_fill_gradientn(colours = rev(myCol2)) +
  ggrepel::geom_text_repel(data = trait_pos_short, aes(x = PC1, y = PC2, label = Trait), max.overlaps = 20) +
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
  labs(x = "PCoA 1 (50%)", y = "PCoA 2 (23%)") +
  ylim(-.48, .48) + xlim(-.48, .48) +
  coord_fixed()

plot_trait_pos

ggsave(
  plot = plot_trait_pos,
  filename = "output/funcspace_trait.png",
  width = 7,
  height = 5,
  units = "in")

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
    data = centroid_threat12,
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
  theme(strip.text = element_text(size=14),
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
  labs(x = "PCoA 1 (50%)", y = "PCoA 2 (23%)") +
  ylim(-.48, .48) + xlim(-.48, .48) +
  coord_fixed()

plot_threat

ggsave(
  plot = plot_threat,
  filename = "output/funcspace_threat.pdf",
  width = 7,
  height = 5,
  units = "in")


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
  theme(strip.text = element_text(size=14),
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
  labs(x = "PCoA 1 (50%)", y = "PCoA 2 (23%)") +
  ylim(-.48, .48) + xlim(-.48, .48) +
  coord_fixed()

plot_threat2

# ggsave(
#   plot = plot_threat2,
#   filename = "output/funcspace_threat_split.pdf",
#   width = 14,
#   height = 5,
#   units = "in")

ggsave(
  plot = plot_threat2,
  filename = "output/funcspace_threat_split.png",
  width = 10,
  height = 5,
  units = "in")


## plot order-family centroids

plot_traits_order <-   ggplot(coordinates, aes(PC1, PC2)) +
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
    data = centroid,
    aes(x = cen.1, y = cen.2),
    shape = 21,
    fill = "white",
    size = 2.5
  ) +
  geom_point(
    data = centroid,
    aes(x = cen.1, y = cen.2),
    shape = 19,
    colour = "black",
    size = 1
  ) +
  scale_fill_gradientn(colours = rev(myCol2)) +
  ggrepel::geom_text_repel(data = centroid, aes(x = cen.1, y = cen.2, label = Order), max.overlaps = 15) +
  theme(strip.text = element_text(size=14),
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
  labs(x = "PCoA 1 (50%)", y = "PCoA 2 (23%)") +
  ylim(-.48, .48) + xlim(-.48, .48) +
  coord_fixed() 

plot_traits_order

ggsave(
  plot = plot_traits_order,
  filename = "output/funcspace_order.png",
  width = 7,
  height = 5,
  units = "in")

plot_traits_family <-   ggplot(coordinates, aes(PC1, PC2)) +
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
    data = centroid_fam,
    aes(x = cen.1, y = cen.2),
    shape = 21,
    fill = "white",
    size = 2.5
  ) +
  geom_point(
    data = centroid_fam,
    aes(x = cen.1, y = cen.2),
    shape = 19,
    colour = "black",
    size = 1
  ) +
  scale_fill_gradientn(colours = rev(myCol2)) +
  ggrepel::geom_text_repel(data = centroid_fam, aes(x = cen.1, y = cen.2, label = Family), max.overlaps = 15) +
  theme(strip.text = element_text(size=14),
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
  labs(x = "PCoA 1 (50%)", y = "PCoA 2 (23%)") +
  ylim(-.48, .48) + xlim(-.48, .48) +
  coord_fixed() 

plot_traits_family

ggsave(
  plot = plot_traits_family,
  filename = "output/funcspace_fam.png",
  width = 7,
  height = 5,
  units = "in")

library(ggpubr)
fig1 <- ggarrange(plot_trait_pos, plot_traits_family, labels = "AUTO", hjust = -1, vjust = 2)
fig1

ggsave(
  plot = fig1,
  filename = "output/figure1.png",
  width = 10,
  height = 5,
  units = "in")

## check relationship between PCs and poptrends----
head(sp_data)
class(sp_coords)

sp_PC <- data.frame(binomial = row.names(sp_coords), sp_coords[,1:2])

sp_trends <- inner_join(sp_data, sp_PC, by = "binomial")

hist(sp_trends$Change.index)

#check with change index
# it does not include Chlamydotis undulata subsp. fuertaventurae, Cursorius cursor, Anthus berthelotii  
ggplot(sp_trends) + geom_point(aes(PC1, Change.index)) +
  geom_smooth(aes(PC1, Change.index), method = "lm") + theme_bw()

# ggplot(sp_trends) + geom_point(aes(PC1, Change.index)) +
#   geom_smooth(aes(PC1, Change.index), method = "loess")

ggplot(sp_trends) + geom_point(aes(PC2, Change.index)) +
  geom_smooth(aes(PC2, Change.index), method = "lm") + theme_bw()

# ggplot(sp_trends) + geom_point(aes(PC2, Change.index)) +
#   geom_smooth(aes(PC2, Change.index), method = "loess")

# check relationship between trends and PC scores
m1 <- lm(Change.index ~ PC1, data = sp_trends)
m2 <- lm(Change.index ~ PC2, data = sp_trends)
m3 <- lm(Change.index ~ PC1 + PC2, data = sp_trends)
m4 <- lm(Change.index ~ PC1*PC2, data = sp_trends)
m1b <- lm(Change.index ~ poly(PC1,2), data = sp_trends)
m2b <- lm(Change.index ~ poly(PC2,2), data = sp_trends)
m3b <- lm(Change.index ~ poly(PC1,2) + poly(PC2,2),, data = sp_trends)
m4b <- lm(Change.index ~ poly(PC1,2)*poly(PC2,2), data = sp_trends)

AIC(m1,m2,m3,m4,m1b,m2b,m3b,m4b)

library(effects)
plot(allEffects(m1))
plot(allEffects(m2b))

summary(m1)
summary(m2b)

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

## check relationship between PCs and scoring
sp_rank <- read.csv("data/sp_rank.csv", stringsAsFactors = F)

sp_trends <- inner_join(sp_trends, sp_rank, by = "binomial")

write.csv(sp_trends, "clean_data/sp_trends.csv", row.names = F)

# End of script -----

