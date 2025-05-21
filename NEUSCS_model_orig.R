################################################################################
#############      Mod Biol Phys Ocean Project     #############################
#############             APR-2025                 #############################
#############      Exploring NEUSCS model          #############################
## by: Alexandra Cabanelas 
################################################################################
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
library(patchwork)
library(cowplot)
help(package="mizer")
help(package="mizerExamples")
help(package="mizerExperimental")

NEUSCS_params <- validParams(NEUSCS_params)

## ------------------------------------------ ##
#           data exploration
## ------------------------------------------ ##
summary(NEUSCS_params)
slotNames(NEUSCS_params)

species_params(NEUSCS_params) # per-species biological and ecological parameters
gear_params(NEUSCS_params)

NEUSCS_params@rates_funcs
names(getRateFunction(NEUSCS_params))
getRateFunction(NEUSCS_params)

NEUSCS_params@w                # size grid
NEUSCS_params@w_full           # extended size grid
NEUSCS_params@interaction      # predator-prey interaction matrix
NEUSCS_params@sc               # abundance scaling
NEUSCS_params@rates_funcs$RDD  # recruitment function
NEUSCS_params@species_params$z0 # a constant background mortality of 0.15 is being used
species_params(NEUSCS_params)$erepro # reproductive efficiency

getEncounter(NEUSCS_params)    # predator-prey encounter rate
getSearchVolume(NEUSCS_params) # predation power of the predator. interpreted as a search volume

interaction_matrix(NEUSCS_params)
getFeedingLevel(NEUSCS_params)
getPredMort(NEUSCS_params) # predation mortality rate 
getExtMort(NEUSCS_params)  # can set it with setExtMort()
getFMort(NEUSCS_params)    
getMort(NEUSCS_params)     # total mortality rate; external+predation+fishing
getInitialEffort(NEUSCS_params) #fishing effort for each

growth_rateNES <- getEGrowth(NEUSCS_params) #instantaneous per-capita growth rate

select(species_params(NEUSCS_params), beta, sigma)
# beta = PPMR = 100
# sigma = width of the feeding kernel = 1.3

species_params(NEUSCS_params)$gamma

ext_mort(NEUSCS_params)
dim(ext_mort(NEUSCS_params))  # Should be [species x weight]
species_params(NEUSCS_params)$z0 

## ------------ ##
#   Data explore plots - steady state
## ------------ ##

plot(NEUSCS_params)

plotGrowthCurves(NEUSCS_params) # von Bertalanffy growth
plotFeedingLevel(NEUSCS_params, # feeding level by size
                 species = NEUSCS_params@species_params$species) 
plotlyPredMort(NEUSCS_params) # predation mortality rate of each species against size
plotlyFMort(NEUSCS_params)    # total fishing mortality of each species by size

# biomass density as a function of weight = num density as a function of log weight 
plotSpectra(NEUSCS_params)            # biomass density 
plotSpectra(NEUSCS_params, power = 0) # number density 
plotSpectra(NEUSCS_params, power = 2) # biomass density log weight = Sheldon spectrum

n_NES <- initialN(NEUSCS_params) #number density 
dimnames(n_NES) # weight in grams at the start of each size class

# Percent of total biomass by weight by species
numbers <- n_NES * dw(NEUSCS_params) #number of fish in each size class
biomass_density <- n_NES * w(NEUSCS_params) # biomass density
biomass <- biomass_density * dw(NEUSCS_params) # biomass 
 
# anitialise an array 
cumulative_biomass <- biomass
# calculate the cumulative sum of all biomasses in previous bins
cumulative_biomass[] <- cumsum(biomass)
# normalise this so that it is given as a percentage of the total biomass
cdf <- cumulative_biomass / cumulative_biomass[1, 100] * 100

ggplot(melt(cdf), aes(x = w, y = value)) +
  geom_line() + 
  labs(x = "Weight [g]",
       y = "% of total biomass")  +
  facet_wrap(~sp) +
  scale_x_log10()


cumulative_biomass[] <- t(apply(biomass, 1, cumsum))
total_biomass_per_species <- rowSums(biomass)
cdf <- sweep(cumulative_biomass, 1, total_biomass_per_species, "/") * 100
matplot(t(cdf), type = "l", log = "x", lty = 1, col = 1:nrow(cdf))

growth_rate_frame <- melt(growth_rateNES)
ggplot(growth_rate_frame) +
  geom_line(aes(x = w, y = value)) +
  scale_x_log10(name = "Weight [g]") +
  scale_y_log10(name = "Growth rate [g/year]") +
  facet_wrap(~sp)

mort_rate_frame <- melt(getMort(NEUSCS_params))
ggplot(mort_rate_frame) +
  geom_line(aes(x = w_prey, y = value)) +
  scale_x_log10(name = "Weight [g]") +
  scale_y_log10(name = "Mortality rate [1/year]") +
  facet_wrap(~prey)
mm <- lm(log(mort_rate_frame$value) ~ log(mort_rate_frame$w_prey))
mm
m0 <- exp(coef(mm)[[1]])
m0 # mortality rate (at 1g?)


nf <- melt(n_NES) %>% 
  filter(w <= 10)

ggplot(nf) +
  geom_line(aes(x = w, y = value)) +
  scale_x_log10(name = "Weight [g]") +
  scale_y_log10(name = "Number density [1/g]") +
  facet_wrap(~sp)

growth_rate_df <- melt(growth_rateNES)
colnames(growth_rate_df)

ggplot(growth_rate_df) +
  geom_line(aes(x = w, y = value, color = sp)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Weight [g]",
       y = "Growth rate [g/year]")

## ------------ ##
#   Simulate
## ------------ ##

#simulating - no fishing 
sim_orig <- project(NEUSCS_params, t_max = 100, effort = 0) #effort = c(small = 0.5, medium = 0.2, large = 0.5

plot(sim_orig) # can add power=2
# Feeding level plot = how satiated an individual is, with 0 being unfed; 1 fully satiated
# Biomass plot = can see when species reach stable equilibrium 
# other plots show last time step of the simulation (in this case equilibrium)

## ------------ ##
#   Plot
## ------------ ##

## ------- ##
#  1)
## ------- ##
p1 <- plotBiomass(sim_orig) +
  theme_bw() +
  labs(title = "Biomass") +
  theme(legend.position = "none")

p2 <- plotPredMort(sim_orig) +
  theme_bw() +
  labs(title = "Predation Mortality") +
  theme(legend.position = "none")

p3 <- plotSpectra(sim_orig) +
  theme_bw() +
  labs(title = "Size-Spectra") +
  theme(legend.title = element_blank())

legend <- get_legend(
  p3 + theme(legend.position = "right", legend.margin = margin(0, 0, 0, 0)) +
    guides(colour = guide_legend(ncol = 1))
)

p3_c <- plotSpectra(sim_orig) +
  theme_bw() +
  labs(title = "Size-Spectra") +
  theme(legend.position = "none")

plots_combined <- p1 / p2 / p3_c

(p1fin_biom_pred_ss <- plot_grid(
  plots_combined,
  legend,
  ncol = 2,
  rel_widths = c(4, 1)
))
#ggsave("origNEUSCS_nofishing.png", plot = p1fin_biom_pred_ss, width = 9, height = 8, dpi = 300, bg = "white")

## ------- ##
#  2)
## ------- ##

p1 <- plotBiomass(sim_orig) +
  theme_bw() +
  labs(title = "Biomass") +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 11, color = "black")
  )
plots_combined <- p3 | p1

legend <- get_legend(
  p3 + theme(legend.position = "right", 
             legend.margin = margin(0, 0, 0, 0),
             legend.title = element_blank()) +
    guides(colour = guide_legend(ncol = 1))
)

(p1fin_ss_biom <- plot_grid(
  plots_combined,
  legend,
  ncol = 2,
  rel_widths = c(4, 1)
))
#ggsave("origNEUSCS_nofishing_biomass_sizespectra2.png", plot = p1fin_ss_biom, width = 12, height = 6, dpi = 300, bg = "white")

## ------------ ##
#   Simulate Plots
## ------------ ##

plotSpectra(sim_orig)      # size spectrum
plotPredMort(sim_orig)     # predation mortality
plotFeedingLevel(sim_orig) # feeding level by size
plotGrowthCurves(sim_orig) 

getBiomass(sim_orig)         # total biomass per species
getYield(sim_orig)           # catch per species
getPredMort(sim_orig)        # predation mortality at size
getFeedingLevel(sim_orig)    # how full each species is


## ------------ ##
#   Number density
## ------------ ##
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

(p_number_density <- ggplot(binned_numbers) +
    geom_line(aes(x = bin_midpoint, y = Number_density)) +
    labs(x = "Weight [g]", y = "Number density [1/g]")+
    scale_x_log10() +
    scale_y_log10()
)

model <- lm(log(Number_density) ~ log(bin_midpoint), data = binned_numbers)
model
slope <- round(coef(model)[2], 2)
# the straight-line approximation has a slope of about -1.127  
# log(N(w)) = 4.7 - 1.1 log (w)
# N(w) = exp(14.6)w^-1.1 = N(1)w^-lambda

# shallower slope = large fish are relatively more abundant 
# maybe makes sense? more data for larger fish here 

library(scales)
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

plotSpectra(NEUSCS_params, power = 0) # number density 

n <- nrow(df_long)
w_min <- min(df_long$weight)
lambda <- 1 + n / sum(log(df_long$weight / w_min)) #Max Likelihood
lambda
# this gives a lambda of 1.1 


## ------------------------------------------ ##           
## ------------------------------------------ ##
## ------------------------------------------ ##
#           Fishing 
## ------------------------------------------ ##
## ------------------------------------------ ##
sim_orig_fished <- project(NEUSCS_params, t_max = 100)
getEffort(sim_orig_fished)

sim_fished_transition <- project(
  object = sim_orig@params,
  effort = c(small = 0.4, medium = 0.3, large = 0.25),
  t_max = 100,
  initial_n = sim_orig@n[100, , ],
  initial_n_pp = sim_orig@n_pp[100, ]
)

plotBiomass(sim_orig_fished)
plotBiomass(sim_fished_transition)

#saveRDS(sim_orig, file = "sim_orig_unfished.rds")
#saveRDS(sim_orig_fished, file = "sim_orig_fished.rds")

## ------------ ##
#   Summary Results
## ------------ ##

dim(N(sim_orig_fished))
finalN(sim_orig_fished) # results at the final time step
NResource(sim_orig_fished)
finalNResource(sim_orig_fished)
getBiomass(sim_orig_fished)
getCommunitySlope(sim_orig_fished)

## ------------ ##
#   Simulate Plots
## ------------ ##

plot(sim_orig_fished)  
plotBiomass(sim_orig_fished)
plotSpectra(sim_orig_fished)      # size spectrum
plotSpectra(sim_orig_fished, power = 2) 
plotPredMort(sim_orig_fished)     # predation mortality
plotFeedingLevel(sim_orig_fished) # feeding level by size
plotGrowthCurves(sim_orig_fished) 


biomass <- getBiomass(sim_orig_fished)
biomass_frame <- melt(biomass)
str(biomass_frame)
pp <- plot_ly(biomass_frame) %>% 
  add_lines(x = ~time, y = ~value, color = ~sp)


#visualize differences in size spectra
plotSpectra2(sim_orig_fished, name1 = "Fished",
             sim_orig, name2 = "Unfished", 
             power = 2)


#We can get the proportion of Herrings in terms of biomass that have a weight above 50g in each of the 10 years
getProportionOfLargeFish(sim_orig, 
                         species = "Herring", 
                         threshold_w = 50, 
                         biomass_proportion = TRUE)

#abundances in final time step = can compare 2 simulations 
relative_abundance <- finalN(sim_orig) / finalN(sim_orig_fished)

total_abund0 <- colSums(finalN(sim_orig))
total_abund1 <- colSums(finalN(sim_orig_fished))
relative_abundance <- total_abund1 / total_abund0
plot(x = w(NEUSCS_params), y = relative_abundance, log = "xy", type = "n", 
     xlab = "Size (g)", ylab = "Relative abundance", ylim = c(0.1, 10))
lines(x = w(NEUSCS_params), y = relative_abundance)
lines(x = c(min(w(NEUSCS_params)), max(w(NEUSCS_params))), y = c(1, 1), lty = 2)


m2_mort <- getPredMort(NEUSCS_params, finalN(sim_orig))[1, ]
m2_fished <- getPredMort(NEUSCS_params, finalN(sim_orig_fished))[1, ]

plot(x = w(NEUSCS_params), y = m2_fished, log = "x", type = "n", 
     xlab = "Size [g]", ylab = "Predation Mortality [1/year]")
lines(x = w(NEUSCS_params), y = m2_mort, lty = 2)
lines(x = w(NEUSCS_params), y = m2_fished)

## ------------ ##
#   Fishing Plots
## ------------ ##
p1 <- plotBiomass(sim_orig_fished) +
  theme_bw() +
  labs(title = "Biomass") +
  theme(legend.position = "none")

p2 <- plotPredMort(sim_orig_fished) +
  theme_bw() +
  labs(title = "Predation Mortality") +
  theme(legend.position = "none")

p3 <- plotSpectra(sim_orig_fished) +
  theme_bw() +
  labs(title = "Size-Spectra") +
  theme(legend.title = element_blank())

legend <- get_legend(
  p3 + theme(legend.position = "right", legend.margin = margin(0, 0, 0, 0)) +
    guides(colour = guide_legend(ncol = 1))
)

p3_c <- plotSpectra(sim_orig_fished) +
  theme_bw() +
  labs(title = "Size-Spectra") +
  theme(legend.position = "none")

plots_combined <- p1 / p2 / p3_c

(p1fin <- plot_grid(
  plots_combined,
  legend,
  ncol = 2,
  rel_widths = c(4, 1)
))
#ggsave("origNEUSCS_withfishing.png", plot = p1fin, width = 9, height = 8, dpi = 300, bg = "white")


# extract final time step
n_array <- sim_orig_fished@n  # 3D array: sp × size × time
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
#   Fished vs unfished plots
## ------------ ##

# final year biomass for each simulation
biomass_unfished <- getBiomass(sim_orig) %>%
  tail(1) %>%
  t() %>%
  as.data.frame()

biomass_fished <- getBiomass(sim_orig_fished) %>%
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

#plot
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

#plot
fvu <- ggplot(biomass_df, aes(x = reorder(Species, -pct_change), y = pct_change, fill = SizeGroup)) +
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
#ggsave("origNEUSCS_fishedvunfished.png", plot = fvu, width = 9, height = 8, dpi = 300, bg = "white")


## ------------ ##
#   to check when reaches stability
## ------------ ##

biomass <- getBiomass(sim_orig)

# percent change year to year for each species
rel_change <- abs(diff(as.matrix(biomass)) / biomass[-nrow(biomass), ]) * 100

# first year where all species change < 1% for N consecutive years (e.g., 10)
threshold <- 1  # percent
window <- 10    # number of consecutive years

# logical matrix where TRUE = stable
stable_matrix <- rel_change < threshold

# earliest year where all species are stable for N consecutive years
stabilized_year <- NA
for (i in 1:(nrow(stable_matrix) - window + 1)) {
  window_slice <- stable_matrix[i:(i + window - 1), ]
  if (all(colSums(window_slice) == window)) {
    stabilized_year <- as.numeric(rownames(biomass)[i + 1])  # +1 since diff shortens by 1
    break
  }
}

stabilized_year


## ------------------------------------------ ##           
## ------------------------------------------ ##
## ------------------------------------------ ##
## ------------------------------------------ ##           
## ------------------------------------------ ##
## ------------------------------------------ ##
## ------------ ##
#   Extra Plots
## ------------ ##

## ------- ##
#  Extra 1)
## ------- ##
# array that contains the predation mortality at time by species by size
pred_mort_orig <- getPredMort(sim_orig) 
idxFinalT(sim_orig)
t_final <- idxFinalT(sim_orig)

# predation mortality at size in the final time step
pred_mort_final <- pred_mort_orig[t_final, , ]

pred_kernel <- getPredKernel(NEUSCS_params)
dimnames(n_NES)
pred_kernel_32 <- pred_kernel[, 55, , drop = FALSE] #32 g bin

ggplot(melt(pred_kernel_32)) +
  geom_line(aes(x = w_prey, y = value)) +
  scale_x_log10() +
  facet_wrap(~sp)
# a predator size 32g likes to fed on prey 1e+01

plotlyDeath(NEUSCS_params, species = "Cod") #who eats cod
plotlyDeath(NEUSCS_params, species = "Butterfish")
plotlyDeath(NEUSCS_params, species = "Cod", proportion = F)

plotDiet(NEUSCS_params, species = "Cod") #who is eaten by cod
plotDiet(NEUSCS_params, species = "Butterfish")

plotDiet(NEUSCS_params, species = 1:9)

psiNES <- melt(getMaturityProportion(NEUSCS_params) * getReproductionProportion(NEUSCS_params))
maturation_lines <- species_params(NEUSCS_params) %>%
  select(species, w_mat, w_mat25) %>%
  rename(sp = species)

ggplot() +
  geom_line(data = psiNES, aes(x = w, y = value)) +
  labs(x = "Weight [g]",
       y = "Proportion invested into reproduction") + 
  facet_wrap(~sp, scales = "free_x") + 
  geom_vline(data = maturation_lines, aes(xintercept = w_mat, group = sp), 
             lty = 2) +
  geom_vline(data = maturation_lines, aes(xintercept = w_mat25, group = sp), 
             lty = 2, col = "grey") +
  theme(legend.position = "none") +
  scale_x_log10()

## ------- ##
#  Extra 2)
## ------- ##

# extract final time step
n_array <- sim_orig@n  # 3D array: sp × size × time
# final time step (t = 100)
n_final <- n_array[dim(n_array)[1], , ]  # Final time step, all species, all sizes
w <- as.numeric(dimnames(n_array)$w) # or sim_orig@params@w ; size class in g
sp <- dimnames(n_array)$sp

df <- as.data.frame(t(n_final))
df$weight <- w

df_long <- df %>%
  pivot_longer(cols = -weight, names_to = "species", values_to = "abundance") %>%
  filter(abundance > 0)

(hist_plot <- ggplot(df_long, aes(x = weight, weight = abundance)) +
    geom_histogram(bins = 40, fill = "blue", color = "black") +
    labs(x = "Weight [g]", y = "Number of fish") +
    scale_x_log10() +
    theme_minimal()
)

hist_plot + facet_wrap(~species, scales = "free")

range(df_long$weight)

log_breaks <- seq(from = -10, to = 11, by = 1)
log_breaks
bin_breaks <- 2 ^ log_breaks # octave scaling 
bin_breaks

ggplot(df_long) +
  geom_histogram(aes(weight, weight = abundance), fill = "blue", colour = "black",
                 breaks = bin_breaks) +
  labs(x = "Weight [g]",
       y = "Number of fish") +
  scale_y_log10()

log_breaks_double <- seq(from = -10, to = 11, by = 0.5)  # now 0.5 instead of 1
bin_breaks_double <- 2 ^ log_breaks_double
bin_breaks_double

ggplot(df_long) +
  geom_histogram(aes(weight, weight = abundance), fill = "blue", colour = "black",
                 breaks = bin_breaks_double) +
  labs(x = "Weight [g]",
       y = "Number of fish") +
  scale_y_log10() + 
  scale_x_log10()

## ------- ##
#  Extra 3)Biomass density
## ------- ##

#average biomass density in each bin = biomass of all the fish in each weight bin / divide that by the width of each bin
# sum weight not numbers 
binned_biomass <- data_with_bins |> 
  group_by(bin) |> 
  summarise(Biomass = sum(weight)) |>
  mutate(bin_start = bin_breaks[-length(bin_breaks)],
         bin_end = bin_breaks[-1],
         bin_width = bin_end - bin_start,
         bin_midpoint = exp((log(bin_start) + log(bin_end)) / 2),
         Biomass_density = Biomass / bin_width)

ggplot(binned_biomass, aes(x = bin_midpoint, y = Biomass_density)) +
  geom_line() + 
  scale_x_log10(name = "Weight [g]") + 
  scale_y_log10(name = "Biomass density") +
  geom_smooth(method = 'lm')

modelbiom <- lm(log(Biomass_density) ~ log(bin_midpoint), data = binned_biomass)
slopebiom <- round(coef(modelbiom)[2], 2)

(biomdensp <- ggplot(binned_biomass, aes(x = bin_midpoint, y = Biomass_density)) +
    geom_line(color = "black", linewidth = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "gray", linetype = "dashed") +
    scale_x_log10(name = "Weight [g]", 
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(name = "Biomass density", 
                  labels = trans_format("log10", math_format(10^.x))) +
    annotate("text", x = min(binned_biomass$bin_midpoint) * 5000, 
             y = max(binned_biomass$Biomass_density) / 6, 
             label = paste("Slope =", slopebiom), 
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
# biomass increases with size (due to flatter slope than -1)

plotSpectra(sim_orig)            # biomass density 
plotSpectra(sim_orig, power = 2) # biomass density log weight = Sheldon spectrum

## ------- ##
#  Extra 4)
## ------- ##

binned_data <- 
  left_join(binned_numbers, binned_biomass,
            by = c("bin", "bin_start", "bin_end", "bin_width", 
                   "bin_midpoint")) |>
  mutate(bin_width_log_w = log10(bin_end) - log10(bin_start),
         Number_density_log_w = Number / bin_width_log_w,
         Biomass_density_log_w = Biomass / bin_width_log_w)
long_data <- binned_data |>
  pivot_longer(cols = contains("density"),
               names_to = "Type", values_to = "Density")

ggplot(long_data, aes(x = bin_midpoint, y = Density, colour = Type)) +
  geom_line() + 
  scale_x_log10(name = "Weight [g]") + 
  scale_y_log10(name = "Density") +
  geom_smooth(method = 'lm', se = FALSE) +
  labs(title = "Alternative ways of representing the size spectrum")

# N(w) = 1/grams = (slope = -2)
# B(w) and N_logw(w) = dimensionless = (slope = -1)
# Biomass density in log weight = B_logw(w) = grams = Sheldon's density
#= (slope = 0)

lm(log(Biomass_density_log_w) ~ log(bin_midpoint), data = binned_data)

p <- ggplot(df_long) +
  geom_density(aes(weight, stat(count), colour = species), adjust = 4) +
  geom_density(aes(weight, stat(count)), colour = "black", lwd = 1.2, adjust = 4) +
  scale_x_continuous(trans = "log10", name = "Weight [g]") +
  scale_y_continuous(trans = "log10", limits = c(1, NA), name = "Number density in log w")

plotly::ggplotly(p)













#########################################################################
#########################################################################
#########################################################################
### END
## ------------ ##
#   More Plots
## ------------ ##

gear_map <- gear_params(NEUSCS_params)[, c("species", "gear")]
biomass_grouped <- biomass_frame %>%
  left_join(gear_map, by = c("sp" = "species"))

palette8 <- brewer.pal(8, "Dark2")
biomass_grouped <- biomass_grouped %>%
  mutate(sp_panel = paste(gear, sp, sep = ": "))
biomass_grouped <- biomass_grouped %>%
  group_by(gear) %>%
  mutate(color_id = as.integer(as.factor(sp))) %>%
  ungroup()


ggplot(biomass_grouped) +
  geom_line(aes(x = time, y = value, colour = sp),
            linewidth = 0.8) +
  facet_wrap(~gear, scales = "free_y") +
  scale_colour_manual(values = rep(palette8, 3)) +
  theme_bw() +
  theme(legend.position = "right") +
  guides(colour = guide_legend(title = "Species (gear group)"))

plot_by_gear <- function(gear_group) {
  data_subset <- biomass_grouped %>% filter(gear == gear_group)
  
  ggplot(data_subset) +
    geom_line(aes(x = time, y = value, colour = sp),
              linewidth = 0.8) +
    scale_colour_manual(values = palette8) +
    labs(title = paste(gear_group, "species"), colour = "Species") +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_blank())
}

plot_by_gear <- function(gear_group, y_axis = TRUE, x_axis = TRUE) {
  data_subset <- biomass_grouped %>% filter(gear == gear_group)
  
  ggplot(data_subset) +
    geom_line(aes(x = time, y = value, colour = sp), linewidth = 0.8) +
    scale_colour_manual(values = palette8) +
    labs(
      title = stringr::str_to_title(gear_group),  # Capitalize: "Small", etc.
      x = if (x_axis) "Time [years]" else NULL,
      y = if (y_axis) "Biomass [g]" else NULL,
      colour = NULL
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.title.y = if (!y_axis) element_blank() else element_text(),
      axis.title.x = if (!x_axis) element_blank() else element_text()
    ) +
    guides(colour = guide_legend(nrow = 4))
}

# Generate individual plots
p_small  <- plot_by_gear("small")
p_medium <- plot_by_gear("medium")
p_large  <- plot_by_gear("large")

library(patchwork)
p_small  <- plot_by_gear("small", y_axis = TRUE,  x_axis = TRUE)
p_medium <- plot_by_gear("medium", y_axis = FALSE, x_axis = TRUE)
p_large  <- plot_by_gear("large", y_axis = FALSE, x_axis = TRUE)

p_small | p_medium | p_large

(p_small + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 4))) | 
  (p_medium + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 4))) | 
  (p_large + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 4)))

(p_small | p_medium | p_large) + plot_layout(ncol = 3)

two_species_biomass <- filter(biomass_frame, sp %in% c("White hake", 
                                                       "Monkfish"))

ggplot(two_species_biomass) +
  geom_line(aes(x = time, y = value, color = sp))


##plotting spectra by hand
final_n <- N(sim_orig)[idxFinalT(sim_orig), , , drop = FALSE]
str(final_n)
nf <- melt(final_n)
str(nf)
nf <- filter(nf, value > 0)

ggplot(nf) +
  geom_line(aes(x = w, y = value, color = sp)) +
  scale_x_log10() +
  scale_y_log10()


maturity <- getMaturityProportion(NEUSCS_params)["Cod", ]
repro_prop <- getReproductionProportion(NEUSCS_params)["Cod", ]
df <- data.frame(Size = w(NEUSCS_params), 
                 Reproduction = repro_prop, 
                 Maturity = maturity, 
                 Total = maturity * repro_prop)
dff <- melt(df, id.vars = "Size", 
            variable.name = "Type", 
            value.name = "Proportion")

ggplot(dff) + geom_line(aes(x = Size, y = Proportion, colour = Type))

## ------------ ##
#   unstable species
## ------------ ##
plotlyBiomass(sim_glim)
tail(getBiomass(sim_glim)) #for stability the last numbers should barely be changing
plotBiomass(sim_glim, species = c("Acadian Redfish", "Scup", "Spiny dogfish"))

ext_vals <- ext_mort(NEUSCS_glim)["Acadian Redfish", ]
plot(w(NEUSCS_glim), ext_vals, type = "l", log = "xy", 
     main = "GLIM Mortality for Acadian Redfish", 
     xlab = "Weight (g)", ylab = "M(w)")

animateSpectra(sim_glim_fish, total = TRUE, power = 2, 
               ylim = c(1e-6, NA), wlim = c(1e-3, NA))