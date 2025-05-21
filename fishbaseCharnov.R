################################################################################
#############      Mod Biol Phys Ocean Project     #############################
#############             May-2025                 #############################
#############    Fish Base parameters Charnov M    #############################
## by: Alexandra Cabanelas 
################################################################################

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
#remotes::install_github("ropensci/rfishbase")
library(rfishbase)
library(dplyr)
library(stringr)
library(TropFishR)
library(ggplot2)

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
max_size <- read.csv("raw/max_lengths.csv") %>%
  mutate(latin_name = str_to_sentence(SCINAME)) %>%
  dplyr::select(-c(SCINAME,X)) %>%
  mutate(
    l_max = ifelse(latin_name == "Sebastes fasciatus", 49, l_max)
  )
# this contains length data to get asymptotic size w_max (size in g)
# these are the largest fish of each species recorded in spring+fall bottom trawl surveys NOAA

# the l max for acadian redfish seems unreasonable, so changed above
length_frequency <- read.csv("raw/length_frequency_trawls.csv")


#herring_latin <- common_to_sci("Butterfish")
#arrange(herring_latin, SpecCode)

# selecting model species
sp <- data.frame(latin_name = c(
  "Peprilus triacanthus",       # Butterfish
  "Clupea harengus",            # Herring
  "Brevoortia tyrannus",        # Atlantic menhaden
  "Limanda ferruginea",         # Yellowtail Flounder - NOT ON FISH BASE!!
  "Stenotomus chrysops",        # Scup
  "Sebastes fasciatus",         # Acadian Redfish
  "Glyptocephalus cynoglossus", # Witch Flounder
  "Micropogonias undulatus",    # Croaker
  "Pseudopleuronectes americanus", # Winter Flounder
  "Centropristis striata",      # Black sea bass
  "Hippoglossoides platessoides", # American Plaice
  "Cynoscion regalis",          # Weakfish
  "Squalus acanthias",          # Spiny dogfish
  "Pomatomus saltatrix",        # Bluefish
  "Paralichthys dentatus",      # Summer flounder
  "Melanogrammus aeglefinus",   # Haddock
  "Lophius americanus",         # Monkfish
  "Brosme brosme",              # Cusk
  "Lopholatilus chamaeleonticeps",  # Tilefish 
  "Pollachius virens",          # Pollock
  "Urophycis tenuis",           # White hake
  "Morone saxatilis",           # Striped bass
  "Gadus morhua",               # Cod
  "Hippoglossus hippoglossus"   # Halibut
)
)

available_tables()

# get K
growth_data <- popgrowth(sp$latin_name)

# lots of entries for same species
# take the median values 
median_K <- growth_data |>
  group_by(Species) |>
  summarise(median_K = median(K))

# convert from length to weight we use the allometric length-weight 
lw_data <- length_weight(sp$latin_name)
#length_weight <- estimate(sp$latin_name, fields = c("Species", "a", "b"))

lw_median <- lw_data |>
  group_by(Species) |>
  summarize(
    a = median(a, na.rm = TRUE),
    b = median(b, na.rm = TRUE),
    .groups = "drop"
  )

###################################################################
# check l max from trawls vs l inifinity from fishbase
Linf_summary <- growth_data %>%
  group_by(Species) %>%
  summarise(L_inf = median(Loo, na.rm = TRUE), .groups = "drop")

Linf_comparison <- max_size %>%
  rename(Species = latin_name) %>%  # to match FishBase naming
  left_join(Linf_summary, by = "Species")

Linf_comparison <- Linf_comparison %>%
  mutate(
    ratio_lmax_Linf = l_max / L_inf,
    diff_lmax_Linf = l_max - L_inf
  )

ggplot(Linf_comparison, aes(x = L_inf, y = l_max)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "FishBase L∞ (cm)",
    y = "Observed l_max from surveys (cm)",
    title = "Comparison of L∞ from FishBase vs l_max from surveys"
  )
###################################################################


sp <- sp |>
  left_join(lw_median, by = c("latin_name" = "Species")) |> 
  left_join(max_size) |>
  mutate(w_max = a * l_max ^ b)

sp <- sp |>
  left_join(median_K, by = c("latin_name" = "Species"))

sp <- sp %>%
  mutate(COMNAME = str_to_sentence(COMNAME))

# since yellowtail flounder is not on fishbase
# https://apps-nefsc.fisheries.noaa.gov/saw/sasi.php
# got its data from here
# one val per region

# L infinity
mean(c(38.07, 42.27, 42.42))  # 40.9 cm 
# K
mean(c(0.477, 0.506, 0.489))  # 0.491 
# b
mean(c(3.103541, 3.025541, 3.107285))  # 3.0788
# a (in log)
exp(mean(c(-12.07506, -11.8009, -12.04101))) # 5.77e-06

sp <- sp %>%
  mutate(
    a = ifelse(latin_name == "Limanda ferruginea", 5.77e-06, a),
    b = ifelse(latin_name == "Limanda ferruginea", 3.0788, b),
    median_K = ifelse(latin_name == "Limanda ferruginea", 0.491, median_K),
    w_max = ifelse(latin_name == "Limanda ferruginea", 5.77e-06 * l_max^3.0788, w_max)
  )

# to match mizer names
name_map <- c(
  "Atlantic herring"    = "Herring",
  "Yellowtail flounder" = "Yellowtail Flounder",
  "Acadian redfish"     = "Acadian Redfish",
  "Witch flounder"      = "Witch Flounder",
  "Atlantic croaker"    = "Croaker",
  "Winter flounder"     = "Winter Flounder",
  "American plaice"     = "American Plaice",
  "Goosefish"           = "Monkfish",
  "Atlantic cod"        = "Cod",
  "Atlantic halibut"    = "Halibut"
)

sp <- sp %>%
  mutate(COMNAME = ifelse(COMNAME %in% names(name_map),
                          name_map[COMNAME],
                          COMNAME))

#write.csv(sp, "output/sp_params_Charnov.csv")






## ------------------------------------------ ##
#            extra -----
## ------------------------------------------ ##

## Growth parameters
# how fast a species grows is determined by the maximum intake rate and the 
#feeding level. 
#Mizer chooses a sensible default for the feeding level and you only need to 
#give the coefficient h of the maximum intake rate. 

maturity_tbl <- rfishbase::maturity(sp$latin_name)

# lots of entries for same species
# take the median values 
median_maturity <- maturity_tbl |>
  group_by(Species) |>
  filter(!is.na(tm), !is.na(Lm)) |>
  summarise(age_mat = median(tm),
            l_mat = median(Lm))

sp <- sp |>
  left_join(median_maturity, by = c("latin_name" = "Species")) |>
  mutate(w_mat = a * l_mat ^ b)

sp$beta <- 100
sp$sigma <- 1.3
sp$biomass_observed <- NA #ideally get from stock assessment reports
sp$biomass_cutoff <- sp$w_mat