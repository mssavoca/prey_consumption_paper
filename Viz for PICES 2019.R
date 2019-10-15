################################################################
# Lunge rates and prey consumption visualizations for PICES 2019----
################################################################

# load data and packages and fuctions 
library(data.table)
library(readxl)
library(forcats)
library(tidyverse)
library(ggstance)
library(mgcv)
library(lme4)
library(lmerTest)
library(MuMIn)
library(ggpubr)


SE = function(x){sd(x)/sqrt(sum(!is.na(x)))}

# Abbreviate a binomial e.g. Balaenoptera musculus -> B. musculus
abbr_binom <- function(binom) {
  paste(str_sub(binom, 1, 1), 
        str_extract(binom, " .*"), 
        sep = ".")
}

# For barplot of what other's have found
prior_prey_est <- tribble(
  ~Species,                 ~prey_wt_kg,   ~Metric,   ~Source,
  "Megaptera novaeangliae", 528, "R", "Trites et al. 1997",
  "Megaptera novaeangliae", 1927, "MDC", "Trites et al. 1997",
  "Megaptera novaeangliae", 532, "R", "Croll et al. 2006",
  "Megaptera novaeangliae", 1612, "MDC", "Croll et al. 2006",
  "Megaptera novaeangliae", 628, "R", "Barlow et al. 2008",
  "Megaptera novaeangliae", 1902, "MDC", "Barlow et al. 2008",
  "Balaenoptera physalus", 709, "R", "Trites et al. 1997",
  "Balaenoptera physalus", 2744, "MDC", "Trites et al. 1997",
  "Balaenoptera physalus", 901, "R", "Croll et al. 2006",
  "Balaenoptera physalus", 2730, "MDC", "Croll et al. 2006",
  "Balaenoptera physalus", 914, "R", "Barlow et al. 2008",
  "Balaenoptera physalus", 2769, "MDC", "Barlow et al. 2008",
  "Balaenoptera musculus", 1230, "R", "Trites et al. 1997",
  "Balaenoptera musculus", 4490, "MDC", "Trites et al. 1997",
  "Balaenoptera musculus", 1120, "R", "Croll et al. 2006",
  "Balaenoptera musculus", 3393, "MDC", "Croll et al. 2006",
  "Balaenoptera musculus", 1501, "R", "Barlow et al. 2008",
  "Balaenoptera musculus", 4547, "MDC", "Barlow et al. 2008"
)

prior_prey_fig <- prior_prey_est %>% 
  filter(Source == "Croll et al. 2006") %>% 
  ggplot(aes(fct_relevel(abbr_binom(Species), "M. novaeangliae", "B. physalus"), prey_wt_kg, 
                                             fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(.~Source) + 
  labs(x = "Species",
       y = bquote('Prey ingested'~(kg~d^-1))) + 
  ylim(0,4000) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(face = "italic"))
prior_prey_fig

# Raincloud plot code----
### This script creates an R function to generate raincloud plots, then simulates
### data for plots. If using for your own data, you only need lines 1-80.
### It relies largely dddon code previously written by David Robinson
### (https://gist.github.com/dgrtwo/eb7750e74997891d7c20)
### and the package ggplot2 by Hadley Wickham

# Check if required packages are installed  
packages <- c("cowplot", "readr", "ggplot2", "dplyr", "lavaan", "smooth", "Hmisc")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}

# Load packages
library(ggplot2)
# Defining the geom_flat_violin (Raincloud) function ----
# Note: the below code modifies the
# existing github page by removing a parenthesis in line 50

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(
                ymin = min(y),
                ymax = max(y),
                xmin = x,
                xmax = x + width / 2
              )
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data,
                              xminv = x,
                              xmaxv = x + violinwidth * (xmax - x)
            )
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(
              plyr::arrange(transform(data, x = xminv), y),
              plyr::arrange(transform(data, x = xmaxv), -y)
            )
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1, ])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(
            weight = 1, colour = "grey20", fill = "white", size = 0.5,
            alpha = NA, linetype = "solid"
          ),
          
          required_aes = c("x", "y")
  )








# # feeding rates generated by Max----

# for CATS data
load("~/Documents/Research Data/Whale research/ShirelCh2/rates_CATS.RData")

lunge_rates_CATS <- lunge_rates %>% 
  mutate_at(vars(Rate), ~replace(., is.nan(.), 0)) %>% 
  mutate(tag_type = "CATS")


# for DTAG data

load("~/Documents/Research Data/Whale research/ShirelCh2/rates_DTAG only_9.17.19_long.RData")

lunge_rates_DTAG <- dtag_lunges_long %>% 
  mutate_at(vars(Rate), ~replace(., is.nan(.), 0)) %>% 
  rename(`Time (hours)` = Hours) %>% 
  mutate(prey_general = "Krill",
         Study_Area = "SoCal",
         Region = "Eastern North Pacific",
         tag_type = "DTAG")


load("~/Documents/Research Data/Whale research/ShirelCh2/rates_DTAG_problems addressed.RData")

lunge_rates_DTAG_problems_addressed <- dtag_lunges_long %>% 
  mutate_at(vars(Rate), ~replace(., is.nan(.), 0)) %>% 
  rename(`Time (hours)` = Hours) %>% 
  mutate(prey_general = "Krill",
         Study_Area = "SoCal",
         Region = "Eastern North Pacific",
         tag_type = "DTAG")

lunge_rates_all <- lunge_rates_DTAG_problems_addressed %>% 
  bind_rows(lunge_rates_CATS) %>% 
  mutate(SpeciesCode = substr(ID, 1, 2))

## REMOVE BAD DEPLOYMENTS ###
bad_ID <- read_csv("Bad ID.csv")

#put all DTAG deployments together 
lunge_rates_all <- lunge_rates_DTAG_problems_addressed %>% 
  bind_rows(lunge_rates_CATS) %>% 
  mutate(SpeciesCode = substr(ID, 1, 2)) %>% 
  anti_join(bad_ID, by = "ID")



# READ IN and COMBINE DATA
# Species-specific average length data (from Shirel's paper) to recalculate average engulfment capacity----
v_data <- read_csv("MWmeasurements.csv") %>% 
  group_by(Species) %>% 
  dplyr::summarize(med_TLm = median(TLm)) %>% 
  rename(CommonName = Species) %>% 
  mutate(SpeciesCode = case_when(CommonName ==  "Blue Whale" ~ "bw",
                                 CommonName == "Fin Whale"~ "bp",
                                 CommonName ==  "Humpback Whale" ~ "mn",
                                 CommonName == "Minke Whale" ~ "bb", 
                                 CommonName == "Bryde's Whale" ~ "be",
                                 CommonName == "Sei Whale" ~ "bs"))
# # v_data$L <- v_data$MW*0.9766  #creates column that converts MW (kg) to liters

# Allometric equations from Shirel's paper----
# creating fucntions from Shirel's paper for MW (in kg) for engulfment capacity in liters for each species where we have a known length
engulf_allo <- tribble(
  ~SpeciesCode, ~slope,   ~intercept,
  "bb",     3.10910,  0.69969,
  "be",     3.1453,   0.5787,
  "bs",     3.05591,   0.77219,
  "bp",     3.54883,  0.15604,
  "bw",     3.667316, -0.014078,
  "mn",     3.24603,  0.85934
)


# Species specific krill prey weight per m^3----
#from DEC's file BaleenWhaleForagingDistBigKrill100Bins, I think this is the NULL distribution
# Gives mean biomass of krill in kg per m^3----
prey_dist <- tribble(
  ~SpeciesCode, ~ln_mean,          ~ln_sd,
  "bb",     log(10^-0.302053984),  log(10^(0.402095747^2))^0.5,    # E. superba
  "bp",     log(10^-0.207001968),  log(10^(0.397598067^2))^0.5,    # T.spin
  "bw",     log(10^-0.200903497),  log(10^(0.391739563^2))^0.5,    # T.spin
  "mn",     log(10^-0.21552581),  log(10^(0.400348263^2))^0.5      # T.spin
)

prey_dist_ENP <- tribble(
  ~SpeciesCode,     ~ln_mean,             ~ln_sd,
  "Best_lower",     log(10^0.021189299),  log(10^(0.255272505^2))^0.5,    # T.spin, ALL FROM MRY
  "Best_upper",     log(10^0.204119983),  log(10^(0.113943352^2))^0.5,    # T.spin, ALL FROM MRY
  "Hypothetical_low",     log(10^-0.74836178),  0.7374944     #log(10^(-0.514278574)^2)^0.5      # had to calculate manually or NaN is produced
)


prey_dist_ENP <- tribble(
  ~SpeciesCode,     ~ln_mean, ~ln_sd,
  "Best_lower",     log(1.05),  log(1.8),    # T.spin, ALL FROM MRY
  "Best_upper",     log(1.6),  log(1.3),    # T.spin, ALL FROM MRY
  "Hypothetical_low",  log(0.1785),  log(0.306)     #log(10^(-0.514278574)^2)^0.5      # had to calculate manually or NaN is produced
)

# example of what a blue whale m3 krill density distribution looks like ----
index = 1:10000
MRY_hyp_low_prey_dens_m3 <- rlnorm(1e4, log(1.05) + log(0.17),   log(1.8))  # this hypothetical low distribution used extremely large ENP krill (28mm); AND target strength for E. superba, which is bigger than 28mm 
MRY_best_lower_prey_dens_m3 <- rlnorm(1e4, log(1.05),  log(1.8)) 
MRY_best_upper_prey_dens_m3 <- rlnorm(1e4, log(1.6),  log(1.3)) 
ex_prey_dens <- as.data.frame(cbind(index, MRY_hyp_low_prey_dens_m3, MRY_best_lower_prey_dens_m3, MRY_best_upper_prey_dens_m3))

ex_prey_dens <- ggplot(ex_prey_dens) + 
  #geom_histogram(binwidth = 0.1, fill = "gray80") +
  #geom_density(color = "blue") +
  #geom_freqpoly(aes(MRY_hyp_low_prey_dens_m3), binwidth = 0.1, color = "gray60") +
  geom_freqpoly(aes(MRY_hyp_low_prey_dens_m3), binwidth = 0.1, color = "gray20") +
  geom_freqpoly(aes(MRY_best_lower_prey_dens_m3), binwidth = 0.1, color = "blue") +
  geom_freqpoly(aes(MRY_best_upper_prey_dens_m3), binwidth = 0.1, color = "red") +
  labs(x = bquote('krill biomass'~(kg~per~m^3))) +
  theme_classic(base_size = 18)
ex_prey_dens




# Whale population data---- 
#Current data from IUCN 2019 Redlist, whaling data compiled in Rocha et al. 2014
pop_data <- read_excel("Filtration project_Whale population and feeding information.xlsx", sheet = 1)


#Read in tag guide----
tag_guide <- read_excel("TAG GUIDE_9.5.19.xlsx", skip = 2) %>%  # Skips first two rows
  rename(Study_Area = `Study_Area     _`,
         SpeciesCode = `Spec      _`,
         whaleLength = `Drone  _`) 

tag_guide_ENP <- tag_guide %>% 
  filter(Study_Area %in% c("Monterey", "Cordell Bank", "SoCal", "San Diego", "WA Coast", "Alaska"))
#View(tag_guide_ENP)  # this is for CATS tags only, need to add in DTags for full dataset

# Master filtration rate file----
filtration_master <- lunge_rates_all %>%
  left_join(select(tag_guide, ID, SpeciesCode, Study_Area, whaleLength), by = c("ID", "SpeciesCode")) %>%
  mutate(Study_Area = coalesce(Study_Area.x, Study_Area.y)) %>% 
  select(-c(Study_Area.x, Study_Area.y)) %>% 
  filter(SpeciesCode %in% c("bb", "be", "bp", "bw", "mn", "bs")) %>% 
  mutate(
    Year = substr(ID, 3, 4),
    whaleLength = parse_number(whaleLength),
    Study_Area = replace_na(Study_Area, "SoCal"),
    CommonName = case_when(
      SpeciesCode == "bw" ~ "Blue Whale",
      SpeciesCode == "bp" ~ "Fin Whale",
      SpeciesCode == "bs" ~ "Sei Whale",
      SpeciesCode == "mn" ~ "Humpback Whale",
      SpeciesCode %in% c("ba", "bb") ~ "Minke Whale", 
      TRUE ~ "Bryde's Whale"),
    Species = case_when(
      SpeciesCode == "bw" ~ "Balaenoptera musculus",
      SpeciesCode == "bp" ~ "Balaenoptera physalus",
      SpeciesCode == "mn" ~ "Megaptera novaeangliae",
      SpeciesCode %in% c("bb","ba") ~ "Balaenoptera bonaerensis", 
      SpeciesCode == "be" ~ "Balaenoptera edeni",
      SpeciesCode == "bs" ~ "Balaenoptera borealis"),
    Region = case_when(
      Study_Area %in% c("Monterey", "SoCal", "Cordell Bank", "San Diego", "WA Coast", "Alaska", "Everett")  ~ "Eastern North Pacific",
      Study_Area %in% c("Stellwagen", "Norway", "Azores", "Greenland") ~ "North Atlantic",
      Study_Area == "South Africa" ~ "South Africa",
      Study_Area == "Antarctic" ~ "Antarctic",
      Study_Area == "Chile" ~ "Chile",
      Study_Area == "Falklands" ~ "Falklands")) %>%  
  # Join population numbers
  left_join(pop_data, by = "Species") %>% 
  # Join species-median size data
  left_join(select(v_data, -CommonName), by = "SpeciesCode") %>% 
  # Join allometric equations
  left_join(engulf_allo, by = "SpeciesCode") %>% 
  # Calculate engulfment volumes
  mutate(
    BestLengthEst = coalesce(whaleLength, med_TLm),
    Engulfment_L = BestLengthEst ^ slope * 10 ^ intercept * 0.9766,
    EngulfVolPerHr = Engulfment_L*Rate) %>% 
  left_join(prey_dist, by = "SpeciesCode") %>% 
  ungroup %>% 
  mutate(prey_hyp_low_lnmean = -1.72316668,
         prey_hyp_low_lnsd = 0.7374944,
         prey_best_low_lnmean = 0.04879016,
         prey_best_low_lnsd = 0.3873574,
         prey_best_upper_lnmean = 0.47000363,
         prey_best_upper_lnsd = 0.1729007)


filtration_master %>%
  group_by(Species) %>%
  summarise(mean_length = median(whaleLength))




# plots of raw data----

pal <- c("B. bonaerensis" = "firebrick3", "B. borealis" = "goldenrod2", "B. edeni" = "darkorchid3",  "M. novaeangliae" = "gray30", "B. physalus" = "chocolate3", "B. musculus" = "dodgerblue2")

filter(Region == "Eastern North Pacific", ID != "bw180904-44") %>% 
  group_by(ID, Species) %>% 
  summarise(whaleLength = first(whaleLength)) %>% 
  ungroup %>% 
  mutate(sp_lbl = factor(abbr_binom(Species)) %>% fct_rev) %>% 
  
  
All_droned_lengths_ENP <- filtration_master %>% 
  filter(Region == "Eastern North Pacific", 
         ID != "bw180904-44") %>% 
  group_by(ID, Species) %>% 
  summarise(whaleLength = first(whaleLength)) %>% 
  ungroup %>% 
  mutate(sp_lbl = factor(abbr_binom(Species)) %>% fct_rev) %>% 
  ggplot(aes(x = sp_lbl, 
             y = whaleLength,
             fill = abbr_binom(Species))) +
  geom_flat_violin(position = position_nudge(x = 0.1, y = 0), alpha = .8) +
  geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
  geom_point(position = position_jitter(width = .05), alpha = 0.6) +
  coord_flip() +
  scale_fill_manual(values = pal) +
  labs(x = "Species",
       y = "Measured length (m)") +
  theme_classic(base_size = 24) +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"))
All_droned_lengths_ENP

dev.copy2pdf(file="All_droned_lengths_ENP.pdf", width=12, height=8)


Engulfment_capacity_ENP <- filtration_master %>% 
  filter(Region == "Eastern North Pacific", 
         ID != "bw180904-44") %>% 
  drop_na(whaleLength) %>% 
  group_by(ID, Species, Engulfment_L) %>% 
  summarise(whaleLength = first(whaleLength)) %>%
  mutate(Engulfment_m3 = Engulfment_L/1000) %>%
  ungroup %>% 
  ggplot(aes(x = fct_reorder(abbr_binom(Species), Engulfment_m3), y = Engulfment_m3,
             fill = abbr_binom(Species))) +
  geom_flat_violin(position = position_nudge(x = 0.075, y = 0), alpha = .8) +
  geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
  geom_point(position = position_jitter(width = .05), alpha = 0.6) +
  coord_flip() +
  scale_fill_manual(values = pal) +
  labs(x = "Species",
       y = bquote('Engulfment capacity'~(m^3))) +
  scale_y_log10(labels = scales::comma) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"))
Engulfment_capacity_ENP

dev.copy2pdf(file="Engulfment_capacity_ENP.pdf", width=12, height=8)


#plot of feeding rates per hour for ENP species 

Feeding_rate_h <- filtration_master %>% 
  filter(prey_general %in% c("Fish", "Krill"), 
         Phase == "Total",
         `Time (hours)` > 1) %>%       # %in% c("Day", "Twilight", "Night")) %>%      
  group_by(ID, Species, Engulfment_L, Phase, Rate, prey_general, Year) %>% 
  summarise(whaleLength = first(whaleLength)) %>%
  ungroup %>%
  mutate(prey_general = factor(prey_general),
         prey_general = recode_factor(prey_general, 
                                      Fish = "Fish-feeding", 
                                      Krill = "Krill-feeding"),
         Engulfment_m3 = Engulfment_L/1000) %>% 
  ungroup %>%
  ggplot(aes(x = fct_reorder(abbr_binom(Species), Rate), y = Rate,
             fill = abbr_binom(Species))) +
  geom_flat_violin(position = position_nudge(x = 0.1, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  geom_point(position = position_jitter(width = .05), alpha = 0.6) +
  coord_flip() +
  scale_fill_manual(values = pal) +
  facet_grid(prey_general~., scales = "free", space = "free") +  
  labs(x = "Species",
       y = bquote('Feeding rate'~(lunges~h^-1))) + 
  theme_classic(base_size = 24) +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"))
Feeding_rate_h

dev.copy2pdf(file="Feeding_rate_h.pdf", width=14, height=6.5)

prey_consumed_h <- filtration_master %>% 
  filter(Phase == "Total") %>%       # %in% c("Day", "Twilight", "Night")) %>%      
  group_by(ID, Species, Engulfment_L, Phase, Rate, prey_general, Year) %>% 
  summarise(whaleLength = first(whaleLength)) %>%
  mutate(Engulfment_m3 = Engulfment_L/1000) %>%
  ungroup %>% 
  ggplot(aes(x = fct_reorder(abbr_binom(Species), Rate), y = Rate,
             fill = abbr_binom(Species))) +
  geom_flat_violin(position = position_nudge(x = 0.1, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  geom_point(position = position_jitter(width = .05), alpha = 0.6) +
  coord_flip() +
  scale_fill_manual(values = pal) +
  #facet_grid(.~prey_general, scales = "free", space = "free") +  # freeing scales doesn't work with flat_violin plots
  labs(x = "Species",
       y = bquote('Feeding rate'~(lunges~h^-1))) + 
  theme_classic(base_size = 24) +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"))
prey_consumed_h


dev.copy2pdf(file="Feeding_rate_h.pdf", width=14, height=6.5)



################################
## CREATE SAMPLING FUNCTION----
################################

#determining daily feeding rates ----

# Get a sample of whale lengths
sample_length <- function(sp_code, size) {
  # Check species code is valid
  if (!sp_code %in% c("bw", "bp", "mn", "bb")) {
    stop(sprintf("%s is an invalid species code", sp_code))
  }
  # Filter to species
  lengths <- filtration_master %>% 
    filter(SpeciesCode == sp_code) %>% 
    drop_na(whaleLength)
  # Verify available whale lengths
  if (nrow(lengths) == 0) {
    stop(sprintf("No lengths available for %s", sp_code))
  }
  # Draw sample
  sample(lengths$whaleLength, size, replace = TRUE)
}

# Seems to be working
sample_rates <- function(rows, keys) {
  rows_wide <- rows %>% 
    select(ID, Phase, Rate, `Time (hours)`) %>% 
    pivot_wider(names_from = Phase, values_from = c(Rate, `Time (hours)`))
  
  # Skip groups with one animal
  if(nrow(rows_wide) == 1) return(tibble())
  
  sample2 <- function(x, prob, size, replace) {
    tryCatch({
      if (all(prob == 0)) {
        0
      } else {
        sample(x, prob = prob, size = size, replace = replace)
      }
    }, error = function(e) browser())
  }
  
  mass_per_day <- function(rate, capacity, lnmean, lnsd) {
    densities <- rlnorm(rate, lnmean, lnsd)
    sum(capacity * densities)
  }
  
  # Sample size
  sz <- 10000
  result <- tibble(day_rate = sample2(rows_wide$Rate_Day, prob = rows_wide$`Time (hours)_Day`, size = sz, replace = TRUE),
                   tw_rate = sample2(rows_wide$Rate_Twilight, prob = rows_wide$`Time (hours)_Twilight`, size = sz, replace = TRUE),
                   night_rate = sample2(rows_wide$Rate_Night, prob = rows_wide$`Time (hours)_Night`, size = sz, replace = TRUE),
                   daily_rate = floor(day_rate*14 + tw_rate*3 + night_rate *7),
                   length_distrib = sample_length(keys$SpeciesCode, sz),
                   slope = rows$slope[1],
                   intercept = rows$intercept[1],
                   measured_engulfment_cap_m3 = (length_distrib ^ slope * 10 ^ intercept * 0.9766)/1000,
                   prey_mass_per_day_hyp_low_kg = map2_dbl(daily_rate, measured_engulfment_cap_m3, mass_per_day, 
                                                           lnmean = rows$prey_hyp_low_lnmean[1], 
                                                           lnsd = rows$prey_hyp_low_lnsd[1]),
                   prey_mass_per_day_best_low_kg = map2_dbl(daily_rate, measured_engulfment_cap_m3, mass_per_day, 
                                                            lnmean = rows$prey_best_low_lnmean[1], 
                                                            lnsd = rows$prey_best_low_lnsd[1]),
                   prey_mass_per_day_best_upper_kg = map2_dbl(daily_rate, measured_engulfment_cap_m3, mass_per_day, 
                                                              lnmean = rows$prey_best_upper_lnmean[1], 
                                                              lnsd = rows$prey_best_upper_lnsd[1]),
                   prey_mass_per_lunge = prey_mass_per_day_best_low_kg / daily_rate)
}

d_strapped <- filtration_master %>% 
  #filter(prey_general == "Krill") %>% 
  group_by(Species, SpeciesCode, prey_general, Year, Region) %>% 
  group_modify(sample_rates) %>% 
  ungroup 


d_strapped_daily_rates_ENP <- d_strapped %>% 
  filter(daily_rate > 5, Region == "Eastern North Pacific",
         prey_general == "Krill") %>% 
  group_by(Species) %>% 
  summarize(avg_daily_lunges = mean(daily_rate))


Daily_rate_ENP <- d_strapped %>% 
  filter(daily_rate > 5,
    Region == "Eastern North Pacific",
         prey_general == "Krill") %>%   #%in% c("Fish", "Krill")) %>%  
  ggplot(aes(x = fct_reorder(abbr_binom(Species), daily_rate), y = daily_rate, 
             fill = abbr_binom(Species))) +
  geom_flat_violin(position = position_nudge(x = 0.1, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  coord_flip() + 
  scale_fill_manual(values = pal) +
  #scale_y_continuous(breaks = seq(from = 0, to = 2500, by = 500)) +
  # facet_grid(prey_general~., scales = "free", space = "free",         # freeing scales doesn't work with flat_violin plots
  #            labeller = labeller(prey_general = names(c("Fish-feeding", "Krill-feeding")))) +  #Labeller not working
  labs(x = "Species",
       y = bquote('Estimated feeding rate'~(lunges~d^-1))) + 
  #ylim(0, 2000) +
  #scale_y_log10(labels = scales::comma, breaks = c(10,100,500,1000)) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none",
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic"))
Daily_rate_ENP

dev.copy2pdf(file="Daily_rate_ENP.pdf", width=13, height=8)


Daily_biomass_ingested_ENP <- d_strapped %>% 
  filter(daily_rate > 5,
         Region == "Eastern North Pacific",
         prey_general == "Krill") %>%  
  ggplot() +
  #hyp-low
  geom_flat_violin(aes(x = fct_reorder(abbr_binom(Species), prey_mass_per_day_hyp_low_kg), y = prey_mass_per_day_hyp_low_kg/1000,
                       fill = abbr_binom(Species)), color = NA, position = position_nudge(x = 0.2, y = 0), alpha = .4, adjust = 2) +
  geom_boxplot(aes(x = fct_reorder(abbr_binom(Species), prey_mass_per_day_hyp_low_kg), y = prey_mass_per_day_hyp_low_kg/1000,
                   fill = abbr_binom(Species)), color = "gray20", width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.1,
               position = position_nudge(x = -0.12, y = 0)) +
  #best-low *our best estimate*
  geom_flat_violin(aes(x = fct_reorder(abbr_binom(Species), prey_mass_per_day_best_low_kg), y = prey_mass_per_day_best_low_kg/1000, 
                       fill = abbr_binom(Species)), position = position_nudge(x = 0.2, y = 0), alpha = 1, adjust = 2) +
  geom_boxplot(aes(x = fct_reorder(abbr_binom(Species), prey_mass_per_day_best_low_kg), y = prey_mass_per_day_best_low_kg/1000, 
                   fill = abbr_binom(Species)), color = "blue", width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.8,
               position = position_nudge(x = 0.12, y = 0)) +
  #best-upper
  geom_flat_violin(aes(x = fct_reorder(abbr_binom(Species), prey_mass_per_day_best_upper_kg), y = prey_mass_per_day_best_upper_kg/1000,
                       fill = abbr_binom(Species)), color = NA, position = position_nudge(x = 0.2, y = 0), alpha = .4, adjust = 2) +
  geom_boxplot(aes(x = fct_reorder(abbr_binom(Species), prey_mass_per_day_best_upper_kg), y = prey_mass_per_day_best_upper_kg/1000,
                   fill = abbr_binom(Species)), color = "red", width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.1) +
  #ylim(0,100) +
  coord_flip() +
  scale_fill_manual(values = pal) +
  scale_y_log10(labels = scales::comma, limits = c(0.1,100), breaks = c(0,1,10,25,100)) +
  geom_hline(yintercept = c(1.814,2.747,4.143), color = c("gray30", "chocolate3", "dodgerblue2"), linetype = "dashed") +
  labs(x = "Species",
       y = bquote('Estimated prey consumed'~(tonnes~d^-1))) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none",
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic"))
Daily_biomass_ingested_ENP

dev.copy2pdf(file="Daily_biomass_ingested_ENP.pdf", width=13, height=8)


Irvine_data <- read_csv(file = "Irvine et al. 2019_med_dur_tag_data.csv")

Irvine_feeding_rates <- Irvine_data %>% 
  group_by(ID, Species, Tag_duration_d) %>% 
  summarise(Total_lunges_fulltag = sum(Total_lunges)) %>% 
  ungroup %>% 
  mutate(daily_rate = Total_lunges_fulltag / Tag_duration_d) %>% 
  ggplot(aes(x = fct_reorder(abbr_binom(Species), daily_rate), y = daily_rate,
             fill = abbr_binom(Species))) +
  geom_flat_violin(position = position_nudge(x = 0.1, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  geom_point(position = position_jitter(width = .05), alpha = 0.6) +
  coord_flip() +
  scale_fill_manual(values = pal) +
  labs(x = "Species",
       y = bquote('Feeding rate'~(lunges~d^-1))) + 
  theme_classic(base_size = 28) +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"))
Irvine_feeding_rates

dev.copy2pdf(file="Irvine_feeding_rates.pdf", width=13, height=8)
