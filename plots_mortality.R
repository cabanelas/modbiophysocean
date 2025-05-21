################################################################################
#############      Mod Biol Phys Ocean Project     #############################
#############             APR-2025                 #############################
#############               PLOTS                  #############################
## by: Alexandra Cabanelas 
################################################################################

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
#install.packages("mizer")
#devtools::install_github("sizespectrum/mizerExamples")
library(tidyverse)
library(mizer)
library(mizerExamples)
library(mizerExperimental)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

## ------------------------------------------ ##
#           Data
## ------------------------------------------ ##
sim_NEUSCS <- readRDS("output/sim_orig_unfished.rds")
sim_NEUSCS_fished <- readRDS("output/sim_orig_fished.rds")

sim_glim <- readRDS("output/NEUSCS_glim_unfished.rds")
sim_glim_fished <- readRDS("output/NEUSCS_glim_fished.rds")

sim_charnov <- readRDS("output/sim_charnov_unfished.rds")
sim_charnov_fished <- readRDS("output/sim_charnov_fished.rds")

## ------------------------------------------ ##
#           Plots
## ------------------------------------------ ##

## ------------ ##
#  1) predation mortality
## ------------ ##

label_unfished <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Unfished", 
           size = 5, fontface = "bold") +
  theme_void()

label_fished <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Fished", 
           size = 5, fontface = "bold") +
  theme_void()

no_axes <- theme(axis.title.x = element_blank(),
                 axis.title.y = element_blank())

p1 <- plotPredMort(sim_NEUSCS) +
  annotate("text", x = Inf, y = Inf, 
           label = "NEUSCS", 
           hjust = 1.1, vjust = 1.5, size = 5) +
  no_axes

p2 <- plotPredMort(sim_NEUSCS_fished) +
  annotate("text", x = Inf, y = Inf, 
           label = "NEUSCS", 
           hjust = 1.1, vjust = 1.5, size = 5) +
  no_axes

p3 <- plotPredMort(sim_glim) +
  labs(y = "Predation mortality") +
  annotate("text", x = Inf, y = Inf, 
           label = "GLIM", 
           hjust = 1.1, vjust = 1.5, size = 5) +
  theme(axis.title.x = element_blank(),
        axis.title.y= element_text(size = 14.2))

p4 <- plotPredMort(sim_glim_fished) +
  annotate("text", x = Inf, y = Inf, 
           label = "GLIM", 
           hjust = 1.1, vjust = 1.5, size = 5) +
  no_axes

p5 <- plotPredMort(sim_charnov) +
  labs(x = "Size (g)") +
  annotate("text", x = Inf, y = Inf, 
           label = "Charnov M", 
           hjust = 1.1, vjust = 1.5, size = 5) +
  theme(axis.title.y = element_blank(),
        axis.title.x= element_text(size = 16))

p6 <- plotPredMort(sim_charnov_fished) +
  labs(x = "Size (g)") +
  annotate("text", x = Inf, y = Inf, 
           label = "Charnov M", 
           hjust = 1.1, vjust = 1.5, size = 5) +
  theme(axis.title.y = element_blank(),
        axis.title.x= element_text(size = 16)) 

(predmortp <- (label_unfished | label_fished) /
  ((p1 | p2) / (p3 | p4) / (p5 | p6)) +
  plot_layout(heights = c(0.05, 1), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key = element_rect(fill = "white", colour = NA),
    legend.key.width = unit(2, "lines"),
    legend.key.height = unit(1, "lines"),
    legend.text = element_text(size = 11)
  )
)
#ggsave("figures/plotsPreMortallmods.png", plot = predmortp, width = 9, height = 9, dpi = 300, bg = "white")

## ------------ ##
#  2) final biomass per species
## ------------ ##

total_biomass <- tibble(
  Model = rep(c("Constant", "GLIM", "Charnov"), each = 2),
  Fishing = rep(c("Unfished", "Fished"), 3),
  Biomass = c(
    sim_NEUSCS %>% getBiomass() %>% tail(10) %>% rowSums() %>% mean(),
    sim_NEUSCS_fished %>% getBiomass() %>% tail(10) %>% rowSums() %>% mean(),
    sim_glim %>% getBiomass() %>% tail(10) %>% rowSums() %>% mean(),
    sim_glim_fished %>% getBiomass() %>% tail(10) %>% rowSums() %>% mean(),
    sim_charnov %>% getBiomass() %>% tail(10) %>% rowSums() %>% mean(),
    sim_charnov_fished %>% getBiomass() %>% tail(10) %>% rowSums() %>% mean()
  )
)

(plotbiomassall <- ggplot(total_biomass, aes(x = Model, y = Biomass, fill = Fishing)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::comma) +
  labs(y = "Total Biomass (g)") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        axis.text = element_text(size = 13, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 14))
)
#ggsave("figures/barplottotbiomass.png", plot = plotbiomassall, width = 9, height = 8, dpi = 300, bg = "white")


## ------------ ##
#  3) 
## ------------ ##

# total biomass density by size class 
get_size_biomass <- function(sim) {
  N <- getN(sim, drop = FALSE)        # w_full × species
  w <- sim@params@w_full              
  dw <- sim@params@dw_full
  
  biomass_density <- rowSums(N) * w   # total biomass per size class
  tibble(
    size = w,
    biomass = biomass_density
  )
}

df_constant <- get_size_biomass(sim_NEUSCS_fished) %>%
  rename(biomass_fished = biomass) %>%
  mutate(biomass_unfished = get_size_biomass(sim_NEUSCS)$biomass,
         model = "Constant")

df_glim <- get_size_biomass(sim_glim_fished) %>%
  rename(biomass_fished = biomass) %>%
  mutate(biomass_unfished = get_size_biomass(sim_glim)$biomass,
         model = "GLIM")

df_charnov <- get_size_biomass(sim_charnov_fished) %>%
  rename(biomass_fished = biomass) %>%
  mutate(biomass_unfished = get_size_biomass(sim_charnov)$biomass,
         model = "Charnov")

plot_df <- bind_rows(df_constant, df_glim, df_charnov) %>%
  mutate(log_ratio = log10(biomass_fished / biomass_unfished))

(fishp <- ggplot(plot_df, aes(x = size, y = log_ratio)) +
  geom_line(size = 1) +
  facet_wrap(~model, scales = "free") +
  scale_x_log10() +
  labs(
    x = "Size (g)",
    y = expression(log[10]*"(Fished / Unfished biomass)"),
  ) +
  theme_minimal(base_size = 14)
)

#ggsave("figures/fishing.png", plot = fishp, width = 12, height = 6, dpi = 300, bg = "white")

## ------------ ##
#  4) mortality
## ------------ ##
species_name <- "Cod"

species_index <- which(species_params(NEUSCS_params)$species == species_name)

# weight grid
w_vals <- w(NEUSCS_params)

mort_orig <- ext_mort(sim_NEUSCS@params)[species_index, ]
mort_glim <- ext_mort(sim_glim@params)[species_index, ]

mort_df <- data.frame(
  weight = w_vals,
  Constant = mort_orig,
  GLIM = mort_glim
) %>%
  pivot_longer(cols = c("Constant", "GLIM"), 
               names_to = "Type", values_to = "Mortality")

ggplot(mort_df, aes(x = weight, y = Mortality, color = Type)) +
  geom_line(linewidth = 1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Body weight (g)",
    y = "External mortality rate (1/year)",
    title = paste("External mortality for", species_name),
    color = "Mortality type"
  ) +
  theme_minimal(base_size = 14)

## ------------ ##
#  5) mortality
## ------------ ##
# sp names and weight grid
species <- species_params(NEUSCS_params)$species
w_vals <- w(NEUSCS_params)

# external mortality arrays
mort_orig_all <- ext_mort(sim_NEUSCS@params)
mort_glim_all <- ext_mort(sim_glim@params)

# long df
mort_df2 <- bind_rows(
  as.data.frame(mort_orig_all) %>%
    mutate(Species = rownames(mort_orig_all), Type = "Constant"),
  as.data.frame(mort_glim_all) %>%
    mutate(Species = rownames(mort_glim_all), Type = "GLIM")
) %>%
  pivot_longer(
    cols = -c(Species, Type),
    names_to = "Weight",
    values_to = "Mortality"
  ) %>%
  mutate(
    Weight = as.numeric(Weight)
  )

ggplot(mort_df2, aes(x = Weight, y = Mortality, color = Species, linetype = Type)) +
  geom_line(size = 1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Body weight (g)",
    y = "External mortality rate (1/year)",
    color = "Species",
    linetype = "Mortality Type"
  ) +
  theme_minimal(base_size = 13)

## ------------ ##
#  6) mortality
## ------------ ##

color_palette <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(unique(mort_df2$Species)))
glim_df <- subset(mort_df2, Type == "GLIM")
const_df <- subset(mort_df2, Type == "Constant")

ggplot() +
  geom_line(data = glim_df, aes(x = Weight, y = Mortality, 
                                color = Species), size = 1) +
  geom_line(data = const_df, aes(x = Weight, y = Mortality), 
            color = "black", linetype = "dashed", size = 1) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(
    values = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(glim_df$Species)))
  ) +
  labs(
    x = "Body weight (g)",
    y = "External mortality rate (1/year)",
    linetype = "Mortality Type"
  ) +
  theme_minimal(base_size = 13) +
  guides(color = "none")  

(mortplots1 <- ggplot() +
  geom_line(data = glim_df,
            aes(x = Weight, y = Mortality, color = Species, linetype = "GLIM"),
            size = 1) +
  geom_line(data = const_df,
            aes(x = Weight, y = Mortality, linetype = "Constant"),
            color = "black",
            size = 1) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(
    values = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(glim_df$Species)))
  ) +
  scale_linetype_manual(
    name = "Mortality Type",
    values = c("GLIM" = "solid", "Constant" = "dashed")
  ) +
  labs(
    x = "Body weight (g)",
    y = "External mortality rate (1/year)"
  ) +
  theme_minimal(base_size = 13) +
  guides(color = "none")  
)
#ggsave("figures/glim_mortality_plot.png", plot = mortplots1, width = 7, height = 5, dpi = 300, bg = "white")


## ------------ ##
#  7) biomass
## ------------ ##

# biomass data
biomass_orig <- getBiomass(sim_NEUSCS)
biomass_glim <- getBiomass(sim_glim)
biomass_Charnov <- getBiomass(sim_charnov)

biomass_df <- data.frame(
  Year = rep(1:101, times = ncol(biomass_orig)),
  Species = rep(colnames(biomass_orig), each = nrow(biomass_orig)),
  Biomass_Orig = as.vector(biomass_orig),
  Biomass_GLIM = as.vector(biomass_glim),
  Biomass_Charn = as.vector(biomass_Charnov)
)

ggplot(biomass_df, aes(x = Year)) +
  geom_line(aes(y = Biomass_Orig, color = "Constant M")) +
  geom_line(aes(y = Biomass_GLIM, color = "GLIM M")) +
  geom_line(aes(y = Biomass_Charn, color = "Charnov M")) +
  facet_wrap(~Species, scales = "free_y") +
  labs(title = "Biomass Trends Under Different Mortality Models",
       y = "Biomass", x = "Year") +
  theme_minimal()

## ------------ ##
#  8) biomass
## ------------ ##

(p1sim <- plotBiomass(sim_NEUSCS) +
  labs(title = "Constant Mortality"))

(p1glim <- plotBiomass(sim_glim) +
  labs(title = "GLIM Mortality"))
#ggsave("figures/biomass_sim.png", plot = p1sim, width = 8, height = 6, dpi = 300, bg = "white")
#ggsave("figures/biomass_glim.png", plot = p1glim, width = 8, height = 6, dpi = 300, bg = "white")

## ------------ ##
#  9) biomass
## ------------ ##

(pspectrasim <- plotSpectra(sim_NEUSCS, time_range = 100) +
  labs(title = "Constant M"))

(pspectraglim <- plotSpectra(sim_glim, time_range = 100) +
  labs(title = "GLIM M"))
#ggsave("figures/spectra_sim.png", plot = pspectrasim, width = 8, height = 6, dpi = 300, bg = "white")
#ggsave("figures/spectra_glim.png", plot = pspectraglim, width = 8, height = 6, dpi = 300, bg = "white")

## ------------ ##
#  10) biomass
## ------------ ##

biomass_orig_AR <- getBiomass(sim_NEUSCS)[, "Acadian Redfish"]
biomass_glim_AR <- getBiomass(sim_glim)[, "Acadian Redfish"]
biomass_df <- data.frame(
  Year = 1:101,
  Biomass_Orig_AR = biomass_orig_AR,
  Biomass_GLIM_AR = biomass_glim_AR
)
ggplot(biomass_df, aes(x = Year)) +
  geom_line(aes(y = Biomass_Orig_AR, color = "Constant Mortality"), linewidth = 1.2) +
  geom_line(aes(y = Biomass_GLIM_AR, color = "GLIM Mortality"), linewidth = 1.2) +
  labs(y = "Biomass (American Plaice)", x = "Year", color = NULL) +
  theme_minimal()

## ------------ ##
#   to check when reaches stability
## ------------ ##

biomass <- getBiomass(sim_glim)

# % change year to year for each species
rel_change <- abs(diff(as.matrix(biomass)) / biomass[-nrow(biomass), ]) * 100

# first year where all species change < 1% for N consecutive years
threshold <- 1  # %
window <- 10    # # of consecutive years

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


