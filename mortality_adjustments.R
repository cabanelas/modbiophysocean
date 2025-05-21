################################################################################
#############      Mod Biol Phys Ocean Project     #############################
#############             APR-2025                 #############################
#############      Adjusting Natural Mortality     #############################
## by: Alexandra Cabanelas 
################################################################################
# This script shows how I adjusted mortality: GLIM and Charnov M

# weight = grams 
# lengths = centimetres
# time = years 
# scaleModel is useful to change the units im using

# NES  ----- NEUSCS_params
#purely size-based predation, i.e., no species-specific interactions
#Set up with three fishing gears targeting small, medium and large species. 
#Vulnerabilities are represented by changing the clearance rate constant (gamma) 
#between species. Calibrated to efforts effort = c(small = 0.4, medium = 0.3, large = 0.25)

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
#install.packages("mizer")
#devtools::install_github("sizespectrum/mizerExamples")
library(tidyverse)
library(mizer)
library(mizerExamples)
library(mizerExperimental)
library(RColorBrewer)
library(plotly)
library(scales)
library(patchwork)
library(cowplot)

#help(package="mizer")
#help(package="mizerExamples")
#help(package="mizerExperimental")

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
NEUSCS_params <- validParams(NEUSCS_params)
spCharnov <- read.csv("sp_params_Charnov.csv")

## ------------------------------------------ ##
#           Mortality adjustments 
## ------------------------------------------ ##
species_params(NEUSCS_params)$z0 
# in the original model, external mortality is constant = 0.15

## ------------------------- ##
#   1) GLIM (Lorenzen 2022)
## ------------------------- ##

NEUSCS_params2 <- NEUSCS_params
# GLIM external mortality: M(w) = z0 / w^c + s * (w / w_inf)
glim_ext_mort <- outer(
  species_params(NEUSCS_params2)$w_inf,
  w(NEUSCS_params2),
  FUN = function(w_inf, w) {
    z0 <- 0.02  # base constant
    c <- 1
    s <- 0.1 #senescence factor
    (z0 / (w^c)) + s * (w / w_inf)
  }
)

# assign dimnames and class 
dimnames(glim_ext_mort) <- list(
  species = species_params(NEUSCS_params2)$species,
  w = signif(w(NEUSCS_params2), 4)
)
class(glim_ext_mort) <- "mizer_array"

# create updated params object 
NEUSCS_glim <- NEUSCS_params2
# set ext_mort directly
ext_mort(NEUSCS_glim) <- glim_ext_mort

# simulate w GLIM mortality - no fishing 
#NEUSCS_glim <- projectToSteady(NEUSCS_glim, t_max = 100, effort = 0)
sim_glim <- project(NEUSCS_glim, t_max = 100, effort = 0)
sim_glim_fished <- project(NEUSCS_glim, t_max = 100)

getEffort(sim_glim_fished)

#saveRDS(sim_glim, file = "NEUSCS_glim_unfished.rds")
#saveRDS(sim_glim_fished, file = "NEUSCS_glim_fished.rds")

## ------------------------- ##
#   2) Gislason-Charnov mortality = Charnov M
## ------------------------- ##
# allometric or size-structured natural mortality function

# Gislason et al. (2010) empirically derived a size-dependent natural mortality 
#function based on observed data across species
# Charnov et al. (2012) provided a mechanistic underpinning, linking it to 
#life-history theory and optimal reproduction strategies

mizer_species <- species_params(NEUSCS_params)$species

# double check order matches
all.equal(spCharnov$COMNAME, species_params(NEUSCS_params)$species)

spCharnov <- spCharnov[match(species_params(NEUSCS_params)$species, 
                             spCharnov$COMNAME), ]

# Charnov mortality: M(w) = (w / w_inf)^(-1.5 / b) * K

K_vec <- spCharnov$median_K
b_vec <- spCharnov$b
w_inf_vec <- species_params(NEUSCS_params)$w_inf
#w_inf_vec <- spCharnov$w_max # use if wanting to use data derived vals
w_grid <- w(NEUSCS_params)

mort_list <- mapply(function(K, b, w_inf) {
  (w_grid / w_inf)^(-1 / b) * K
}, K_vec, b_vec, w_inf_vec, SIMPLIFY = FALSE)

charnov_ext_mort <- do.call(rbind, mort_list)

charnov_ext_mort[charnov_ext_mort > 4] <- 4 

dimnames(charnov_ext_mort) <- list(
  species = species_params(NEUSCS_params)$species,
  w = signif(w_grid, 4)
)

class(charnov_ext_mort) <- "mizer_array"

# assign to a new model
NEUSCS_charnov <- NEUSCS_params
ext_mort(NEUSCS_charnov) <- charnov_ext_mort

sim_charnov <- project(NEUSCS_charnov, effort = 0, t_max = 100)
sim_charnov_fished <- project(NEUSCS_charnov, t_max = 100)

#saveRDS(sim_charnov, file = "sim_charnov_unfished.rds")
#saveRDS(sim_charnov_fished, file = "sim_charnov_fished.rds")
summary(charnov_ext_mort)
range(charnov_ext_mort)

## ------------ ##
#   Summary Results
## ------------ ##

dim(N(sim_glim))
finalN(sim_glim) # results at the final time step
NResource(sim_glim)
finalNResource(sim_glim)
getBiomass(sim_glim)
getMort(NEUSCS_glim)
getPredMort(NEUSCS_glim)
getCommunitySlope(sim_glim)

ext_mort(NEUSCS_glim)
dim(ext_mort(NEUSCS_glim))  # Should be [species x weight]
species_params(NEUSCS_glim)$z0 

## ------------ ##
#   Simulate Plots
## ------------ ##
plot(sim_glim) 
plot(sim_glim_fished) 
plot(sim_charnov)
plot(sim_charnov_fished)

plotBiomass(sim_glim)
plotBiomass(sim_glim_fished)
plotBiomass(sim_charnov)
plotBiomass(sim_charnov_fished)

plotSpectra(sim_glim)      # size spectrum
plotSpectra(sim_glim_fished) 
plotSpectra(sim_charnov) 
plotSpectra(sim_charnov_fished) 




plotFeedingLevel(sim_glim) # feeding level by size

plotGrowthCurves(sim_glim) 

## ------------------------------------------ ##
#           Plots
## ------------------------------------------ ##

## ------------------------- ##
#   1) GLIM (Lorenzen 2022)
## ------------------------- ##

## ------------ ##
#   Biomass,Pred,Spectra
## ------------ ##
p1 <- plotSpectra(sim_glim) +
  theme_bw() +
  labs(title = "Size-Spectra") +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 11, color = "black")
  )

p2 <- plotBiomass(sim_glim) +
  theme_bw() +
  labs(title = "Biomass") +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 11, color = "black")
  )

plots_combined <- p1 | p2
legend <- get_legend(
  p1 + theme(legend.position = "right", 
             legend.margin = margin(0, 0, 0, 0),
             legend.title = element_blank()) +
    guides(colour = guide_legend(ncol = 1))
)

(p1fin <- plot_grid(
  plots_combined,
  legend,
  ncol = 2,
  rel_widths = c(4, 1)
))
#ggsave("glim_nofishing2.png", plot = p1fin, width = 12, height = 6, dpi = 300, bg = "white")


## ------------ ##
#   Number density
## ------------ ##
# extract final time step
n_array <- sim_glim@n  # 3D array: sp × size × time
# final time step (t = 100)
n_final <- n_array[dim(n_array)[1], , ]  # Final time step, all species, all sizes
w <- as.numeric(dimnames(n_array)$w) # or sim_orig@params@w ; size class in g
sp <- dimnames(n_array)$sp

df <- as.data.frame(t(n_final))
df$weight <- w

df_long <- df %>%
  pivot_longer(cols = -weight, names_to = "species", values_to = "abundance") %>%
  filter(abundance > 0)

min_w <- min(df_long$weight)
max_w <- max(df_long$weight)
log_breaks <- seq(floor(log2(min_w)), ceiling(log2(max_w)), by = 1)
bin_breaks <- 2^log_breaks

# to understand concept of bin size; do it by hand
data_with_bins <- df_long |>
  mutate(bin = cut(weight, breaks = bin_breaks, right = FALSE,
                   labels = FALSE))
head(data_with_bins)

# group the data by bin and calculate the number of fish in each bin
binned_numbers <- data_with_bins |> 
  group_by(bin) |> 
  summarise(Number = n())
binned_numbers %>% print(n = 29)

# average number density
binned_numbers <- binned_numbers |> 
  mutate(bin_start = bin_breaks[-length(bin_breaks)],
         bin_end = bin_breaks[-1],
         bin_width = bin_end - bin_start,
         Number_density = Number / bin_width)
binned_numbers %>% print(n = 29)

binned_numbers <- binned_numbers |>
  mutate(bin_midpoint = exp((log(bin_start) + log(bin_end)) / 2))

model <- lm(log(Number_density) ~ log(bin_midpoint), data = binned_numbers)
model
slope <- round(coef(model)[2], 2)
# the straight-line approximation has a slope of about -1.127  

# log(N(w)) = 4.7 - 1.1 log (w)
# N(w) = exp(14.6)w^-1.1 = N(1)w^-lambda

# shallower slope = large fish are relatively more abundant 
# maybe makes sense? more data for larger fish here 

(numdensp <- ggplot(binned_numbers, aes(x = bin_midpoint, y = Number_density)) +
    geom_line(color = "black", linewidth = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "gray", linetype = "dashed") +
    scale_x_log10(name = "Weight [g]", 
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(name = "Number density [1/g]", 
                  labels = trans_format("log10", math_format(10^.x))) +
    annotate("text", x = min(binned_numbers$bin_midpoint) * 40, 
             y = max(binned_numbers$Number_density) / 6, 
             label = paste("Slope =", slope), 
             hjust = 0, size = 5) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = 0.3),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(size = 15, color="black"),
      axis.text = element_text(size = 13, color="black")
    )
)
#ggsave("origNEUSCS_nofishing_numDensity.png", plot = numdensp, width = 9, height = 8, dpi = 300, bg = "white")

## ------------ ##
#   FISHED VS UNFISHED PLOTS
## ------------ ##
 
# extract final year biomass for each simulation
biomass_unfished <- getBiomass(sim_glim) %>%
  tail(1) %>%
  t() %>%
  as.data.frame()

biomass_fished <- getBiomass(sim_glim_fished) %>%
  tail(1) %>%
  t() %>%
  as.data.frame()

colnames(biomass_unfished) <- "Unfished"
colnames(biomass_fished) <- "Fished"
biomass_unfished$Species <- rownames(biomass_unfished)
biomass_fished$Species   <- rownames(biomass_fished)

biomass_df <- biomass_unfished %>%
  left_join(biomass_fished, by = "Species") %>%
  mutate(pct_change = 100 * (Fished - Unfished) / Unfished)


ggplot(biomass_df, aes(x = reorder(Species, -pct_change), y = pct_change)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Fished vs Unfished NEUSCS Simulation",
    y = "% Change in Biomass"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.title.x = element_blank())

size_class_df <- data.frame(
  Species = c(
    "Butterfish", "Herring", "Atlantic menhaden", "Yellowtail Flounder", "Scup",
    "Acadian Redfish", "Witch Flounder", "Croaker",
    "Winter Flounder", "Black sea bass", "American Plaice", "Weakfish", "Spiny dogfish",
    "Bluefish", "Summer flounder", "Haddock",
    "Monkfish", "Cusk", "Tilefish", "Pollock", "White hake", "Striped bass", "Cod", "Halibut"
  ),
  SizeGroup = c(
    rep("Small", 8),
    rep("Medium", 8),
    rep("Large", 8)
  )
)

biomass_df <- biomass_df %>%
  left_join(size_class_df, by = "Species")

# plot with size group color
(fvu <- ggplot(biomass_df, aes(x = reorder(Species, -pct_change), y = pct_change, fill = SizeGroup)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c(Small = "#56B4E9", Medium = "#F0E442", Large = "#E69F00")) +
  labs(
    title = "Percent Change in Biomass (Fished vs Unfished)",
    x = "Species", y = "% Change in Biomass"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.title.x = element_blank(),
        legend.title = element_blank())
)
#ggsave("glim_fishedvunfished.png", plot = fvu, width = 9, height = 8, dpi = 300, bg = "white")

## ------------------------------------------ ##
#           Plots
## ------------------------------------------ ##

## ------------------------- ##
#   2) Gislason-Charnov mortality = Charnov M
## ------------------------- ##
## ------------ ##
#   Biomass,Pred,Spectra
## ------------ ##
p1 <- plotSpectra(sim_charnov) +
  theme_bw() +
  labs(title = "Size-Spectra") +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 11, color = "black")
  )

p2 <- plotBiomass(sim_charnov) +
  theme_bw() +
  labs(title = "Biomass") +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 11, color = "black")
  )

plots_combinedsim_charnov <- p1 | p2
legend <- get_legend(
  p1 + theme(legend.position = "right", 
             legend.margin = margin(0, 0, 0, 0),
             legend.title = element_blank()) +
    guides(colour = guide_legend(ncol = 1))
)

(p1fin <- plot_grid(
  plots_combinedsim_charnov,
  legend,
  ncol = 2,
  rel_widths = c(4, 1)
))
#ggsave("charnov_nofishingv2.png", plot = p1fin, width = 12, height = 6, dpi = 300, bg = "white")

## ------------ ##
#   Number density
## ------------ ##
# extract final time step
n_array <- sim_glim@n  # 3D array: sp × size × time
# final time step (t = 100)
n_final <- n_array[dim(n_array)[1], , ]  # Final time step, all species, all sizes

w <- as.numeric(dimnames(n_array)$w) # or sim_orig@params@w ; size class in g
sp <- dimnames(n_array)$sp

df <- as.data.frame(t(n_final))
df$weight <- w

df_long <- df %>%
  pivot_longer(cols = -weight, names_to = "species", values_to = "abundance") %>%
  filter(abundance > 0)

min_w <- min(df_long$weight)
max_w <- max(df_long$weight)
log_breaks <- seq(floor(log2(min_w)), ceiling(log2(max_w)), by = 1)
bin_breaks <- 2^log_breaks

# to understand concept of bin size; do it by hand
data_with_bins <- df_long |>
  mutate(bin = cut(weight, breaks = bin_breaks, right = FALSE,
                   labels = FALSE))
head(data_with_bins)

# group the data by bin and calculate the number of fish in each bin
binned_numbers <- data_with_bins |> 
  group_by(bin) |> 
  summarise(Number = n())
binned_numbers %>% print(n = 29)

# average number density
binned_numbers <- binned_numbers |> 
  mutate(bin_start = bin_breaks[-length(bin_breaks)],
         bin_end = bin_breaks[-1],
         bin_width = bin_end - bin_start,
         Number_density = Number / bin_width)
binned_numbers %>% print(n = 29)

binned_numbers <- binned_numbers |>
  mutate(bin_midpoint = exp((log(bin_start) + log(bin_end)) / 2))

model <- lm(log(Number_density) ~ log(bin_midpoint), data = binned_numbers)
model
slope <- round(coef(model)[2], 2)
# the straight-line approximation has a slope of about -1.127  

# log(N(w)) = 4.7 - 1.1 log (w)
# N(w) = exp(14.6)w^-1.1 = N(1)w^-lambda

# shallower slope = large fish are relatively more abundant 
# maybe makes sense? more data for larger fish here 

(numdensp <- ggplot(binned_numbers, aes(x = bin_midpoint, y = Number_density)) +
    geom_line(color = "black", linewidth = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "gray", linetype = "dashed") +
    scale_x_log10(name = "Weight [g]", 
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(name = "Number density [1/g]", 
                  labels = trans_format("log10", math_format(10^.x))) +
    annotate("text", x = min(binned_numbers$bin_midpoint) * 40, 
             y = max(binned_numbers$Number_density) / 6, 
             label = paste("Slope =", slope), 
             hjust = 0, size = 5) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = 0.3),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(size = 15, color="black"),
      axis.text = element_text(size = 13, color="black")
    )
)
#ggsave("origNEUSCS_nofishing_numDensity.png", plot = numdensp, width = 9, height = 8, dpi = 300, bg = "white")

## ------------ ##
#   FISHED VS UNFISHED PLOTS
## ------------ ##

# extract final year biomass for each simulation
biomass_unfished <- getBiomass(sim_charnov) %>%
  tail(1) %>%
  t() %>%
  as.data.frame()

biomass_fished <- getBiomass(sim_charnov_fished) %>%
  tail(1) %>%
  t() %>%
  as.data.frame()

colnames(biomass_unfished) <- "Unfished"
colnames(biomass_fished) <- "Fished"
biomass_unfished$Species <- rownames(biomass_unfished)
biomass_fished$Species   <- rownames(biomass_fished)

biomass_df <- biomass_unfished %>%
  left_join(biomass_fished, by = "Species") %>%
  mutate(pct_change = 100 * (Fished - Unfished) / Unfished)


ggplot(biomass_df, aes(x = reorder(Species, -pct_change), y = pct_change)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Fished vs Unfished NEUSCS Simulation",
    y = "% Change in Biomass"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.title.x = element_blank())

size_class_df <- data.frame(
  Species = c(
    "Butterfish", "Herring", "Atlantic menhaden", "Yellowtail Flounder", "Scup",
    "Acadian Redfish", "Witch Flounder", "Croaker",
    "Winter Flounder", "Black sea bass", "American Plaice", "Weakfish", "Spiny dogfish",
    "Bluefish", "Summer flounder", "Haddock",
    "Monkfish", "Cusk", "Tilefish", "Pollock", "White hake", "Striped bass", "Cod", "Halibut"
  ),
  SizeGroup = c(
    rep("Small", 8),
    rep("Medium", 8),
    rep("Large", 8)
  )
)

biomass_df <- biomass_df %>%
  left_join(size_class_df, by = "Species")

# plot with size group color
(fvu <- ggplot(biomass_df, aes(x = reorder(Species, -pct_change), y = pct_change, fill = SizeGroup)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c(Small = "#56B4E9", Medium = "#F0E442", Large = "#E69F00")) +
    labs(
      title = "Percent Change in Biomass (Fished vs Unfished)",
      x = "Species", y = "% Change in Biomass"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          axis.title.x = element_blank(),
          legend.title = element_blank())
)
#ggsave("charnov_fishedvunfished.png", plot = fvu, width = 9, height = 8, dpi = 300, bg = "white")









#####
# trash
biomass_ts <- getBiomass(NEUSCS_glim)
cv_fished <- apply(biomass_ts, 2, function(x) sd(x)/mean(x))
# Higher CV = more unstable

ext_mort <- outer(
  species_params(NEUSCS_params3)$w_inf,
  w(NEUSCS_params3),
  FUN = function(w_inf, w) {
    zeta1 <- 0.2         # mortality constant 0.27
    A0 <- 3.66            # growth constant
    b <- 0.6787           # allometric growth exponent
    h <- 0.4          # mortality scaling exponent 0.4641
    a_lambda_T <- 1       # metabolic temp effect (1 = not modeling temp)
    
    lambda <- exp(zeta1 * (A0 / 3)) * a_lambda_T
    lambda * w^(-h) * w_inf^(h + b - 1)
  }
)
