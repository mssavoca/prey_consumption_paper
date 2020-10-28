# Analysis and Figure code for submission Savoca et al. ----
## High-resolution foraging measurements reveal baleen whales as global ecosystem engineers


library(tidyverse)
library(maptools)
library(readxl)
library(ggpubr)
library(purrr)

SE = function(x){sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x)))}

# color palette for figures
pal <- c("B. bonaerensis" = "firebrick3", "B. borealis" = "goldenrod2", "B. brydei" = "darkorchid3",  "M. novaeangliae" = "gray30", "B. physalus" = "chocolate3", "B. musculus" = "dodgerblue2")

# Abbreviate a binomial e.g. Balaenoptera musculus -> B. musculus
abbr_binom <- function(binom) {
  paste(str_sub(binom, 1, 1), 
        str_extract(binom, " .*"), 
        sep = ".")
}

# summarize columns in a dataframe with an additional last row


col_summ <- function(tbl, f) {
  rbind(
    tbl,
    purrr::map(tbl, ~ if(is.numeric(.x)) f(.x) else NA)
  )
}


# Raincloud GEOM_FLAT_VIOLIN plot code----
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




# Feeding rates from MFC ----

load("rates_MFC_627.RData") # Updated CATS deployment summaries

load("rates_MedTermTag.RData") # Medium-term deployment summary

load("rates[1].RData") # DTAG deployment summary



#Read in tag guide----
tag_guide <- read_excel("TAG GUIDE_5.27.20.xlsx", skip = 2) %>%  # Skips first two rows
  rename(Study_Area = `Study_Area     _`,
         SpeciesCode = `Spec      _`,
         whaleLength = `Drone  _`)  #loads the tag guide 


dtag_lunges_long <- dtag_lunges_long %>% 
  mutate(prey_general = "Krill",
         tag_type = "DTAG")

rates_longer <- function(depid, day_lunges, night_lunges, day_hours, night_hours, day_rate, night_rate, species) {
  tibble(
    ID = depid,
    SpeciesCode = species,
    Phase = c("Day", "Night", "Total"),
    Lunges = c(day_lunges, night_lunges, day_lunges + night_lunges),
    `Time (hours)` = c(day_hours, night_hours, day_hours + night_hours),
    Rate = c(day_rate, night_rate, (day_lunges + night_lunges) / (day_hours + night_hours))
  )
}
medterm_lunges_long <- pmap_dfr(medterm_lunge_rates, rates_longer) %>%   
  mutate(prey_general = "Krill", 
         tag_type = "MedDur")

  
All_rorqual_deployments <- dtag_lunges_long %>% 
  bind_rows(lunge_rates) %>% 
  mutate(`Time (hours)` = coalesce(Hours, `Time (hours)`),
         tag_type = replace_na(tag_type, "CATS"),
         SpeciesCode = substr(ID, 1, 2)) %>% 
  dplyr::select(-Hours) %>% 
  bind_rows(medterm_lunges_long) %>% 
  left_join(
    dplyr::select(tag_guide, ID, Study_Area, whaleLength), by = "ID") %>% 
  mutate(whaleLength = parse_number(whaleLength)) %>%  
  # manually remove deployments from Azores, Chile, and Norway
  filter(!Study_Area %in% c("Chile", "Azores", "Norway")) %>%  
  filter(!ID %in% c("mn161107-36b", "mn161106-36b", "mn161105-36", "mn151103-7", "mn151103-4", "mn151103-3")) %>%   #includes DEC's CB 2016 whales
  filter(!ID %in% c("mn151031-3", "mn151031-4")) # Removes DEC's CB 2016 whales; ALL SA humpbacks taken out if left on



summ_table <- All_rorqual_deployments %>% filter(Phase == "Total") %>% 
  col_summ(sum)

#write.csv(summ_table, "All_rorqual_deployments.csv")



# droned whales----
# non-tagged unique individuals from KC
KC_nontagged_lengths <- read_xlsx("Whales_for_Savoca_untagged_fromKC.xlsx", skip = 1) %>% 
  rename(length = "Max",
         name = "Animal_ID") %>% 
  mutate(SpeciesCode = case_when(Species == "Antarctic minke" ~ "bb",
                                 Species == "blue" ~ "bw",
                                 Species == "humpback" ~ "mn"))



# bringing in whale length data from Paolo; using allometric equations from Shirel to convert to engulfment capacity
# Allometric equations from Shirel UPDATED from her paper
# creating fucntions using new data from Shirel and DC for for engulfment capacity in m3 for each species where we have a known length

engulf_allo <- tribble(
  ~SpeciesCode, ~slope,   ~intercept,
  "bb",     3.11151,  -2.31328,
  "be",     3.1206,   -2.4039,
  "bs",     3.05563,   -2.23809,
  "bp",     3.54883,  -2.85446,
  "bw",     3.667115, -3.024580,
  "mn",     3.24603,  -2.15150
)

All_rorqual_deployments_forjoin <- All_rorqual_deployments %>% 
  rename(name = ID)

tag_guide_forjoin <- tag_guide %>% 
  rename(name = ID)

whale_lengths <- read.csv("whale_masses.csv") %>% 
  left_join(select(All_rorqual_deployments_forjoin, name, Study_Area), by = "name") %>% 
  left_join(select(tag_guide_forjoin, name, Study_Area), by = "name") %>% 
  distinct() %>% 
  rename(SpeciesCode = "species") %>% 
  full_join(select(KC_nontagged_lengths, name, length, SpeciesCode, Region)) %>% 
  left_join(engulf_allo, by = "SpeciesCode") %>% 
  mutate(
    Study_Area = coalesce(Study_Area.x, Study_Area.y),
    Region = coalesce(Study_Area, Region),
    Region = case_when(Region == "Antarctic" ~ "WAP",
                       Study_Area %in% c("SoCal", "Monterey") ~ "CCE",
                       Study_Area == "South Africa" ~ "SA",
                       SpeciesCode == "bw" ~ "CCE",
                       SpeciesCode == "bb" ~ "WAP"),
    Engulfment_m3 = length ^ slope * 10 ^ intercept,
         #numbers from Lockyer 1976
         Lockyer_a = case_when(SpeciesCode == "bb" ~ 0.0496,   # From Table 1, the rest from Table 2
                               SpeciesCode == "be" ~ 0.0122,
                               SpeciesCode == "mn" ~ 0.0158,
                               SpeciesCode == "bs" ~ 0.0242,
                               SpeciesCode == "bp" ~ 0.0015,  # using Northern Hemisphere fin whale numbers 
                               SpeciesCode == "bw" ~ 0.0004   # using B.m. brevicauda numbers as they are most similar to ENP blues
         ),
         Lockyer_b = case_when(SpeciesCode == "bb" ~ 2.31,    # From Table 1, the rest from Table 2
                               SpeciesCode == "be" ~ 2.74,
                               SpeciesCode == "mn" ~ 2.95,
                               SpeciesCode == "bs" ~ 2.43,
                               SpeciesCode == "bp" ~ 3.46,   # using Northern Hemisphere fin whale numbers
                               SpeciesCode == "bw" ~ 3.97    # using B.m. brevicauda numbers as they are most similar to ENP blues
         ),
         Lockyer_mass_t = (Lockyer_a*length ^Lockyer_b)*1.06,  # Accounting for fluid losses; Lockyer table 2
         Species = case_when(
           SpeciesCode == "bw" ~ "Balaenoptera musculus",
           SpeciesCode == "bp" ~ "Balaenoptera physalus",
           SpeciesCode == "mn" ~ "Megaptera novaeangliae",
           SpeciesCode %in% c("bb","ba") ~ "Balaenoptera bonaerensis", 
           SpeciesCode == "be" ~ "Balaenoptera brydei",
           SpeciesCode == "bs" ~ "Balaenoptera borealis"))

# print out, add missing Regions, put back in
# write_csv(whale_lengths,"whale_lengths2.csv")
# whale_lengths2 <- read_csv("whale_lengths2.csv")

measured_length_wt <- whale_lengths %>% 
  group_by(SpeciesCode) %>% 
  summarize(sample_size = n(),
            avg_wt = mean(mass, na.rm = TRUE),
            SE_wt = SE(mass),
            avg_length = mean(length),
            SE_length = SE(length))
#write_csv(measured_length,"whale_lengths_summary.csv")

# Dave's KRILL data (updated 10/20/20)----
#natural log transformed, which is what we were doing in the prior iteration

pooled_log_mean <- function(m1, m2) {
  a = (log(m1) + log(m2))/2
  exp(a)  #turns values back into biomass units (kg/m3)
}


# pooled sd
pooled_sd_mean <- function(sd1, sd2, n1, n2) {
  
  a=log(sd1)
  b=log(sd2)
  
  pooled_sd = sqrt(((n1-1)*a^2 + (n2-1)*b^2)/(n1+n2-2))  # https://www.statisticshowto.datasciencecentral.com/pooled-standard-deviation/ 
  
  exp(pooled_sd)   #turns values back into biomass units (kg/m3)
}


# # read in, clean, combine data
# MRY_krill_data <-   read_csv("MontereyKrillData_10.21.20.csv") %>% 
#   #read_csv("MontereyKrillDataDiveMeans.csv") %>% 
#   filter(!Species %in% c("bigBW", "bb")) %>% 
#   select(Species:BiomassTop50sd) %>%    # remember: biomass is in kg/m3
#   mutate(
#     Biomass_hyp_low = exp(log(Biomass) + log(0.17)),
#     Study_Area = "Monterey",
#     Region = "Eastern North Pacific") %>% 
#   rename(SpeciesCode = "Species") 
# 
# SoCal_krill_data <-   read_csv("SoCalKrillData_10.21.20.csv") %>% 
#   #read_csv("SoCalKrillDataDiveMeans.csv") %>% 
#   filter(!Species %in% c("bigBW", "bb")) %>% 
#   select(Species:BiomassTop50sd) %>%    # remember: biomass is in kg/m3
#   mutate(
#     Biomass_hyp_low = exp(log(Biomass) + log(0.17)),
#     Study_Area = "SoCal",
#     Region = "Eastern North Pacific") %>% 
#   rename(SpeciesCode = "Species") 
# 
# # combining Monterey and SoCal prey data
# ENP_krill_data <- rbind(MRY_krill_data, SoCal_krill_data) %>% 
#   #mutate(Study_Area2 = Study_Area) %>% 
#   pivot_wider(names_from = Study_Area, values_from = c(`Num Days`:Biomass_hyp_low)) %>% 
#   #rename(Study_Area = Study_Area2) %>% 
#   mutate(
#     `Num Days` = `Num Days_Monterey` + `Num Days_SoCal`, 
#     Biomass_hyp_low = pooled_log_mean(Biomass_hyp_low_Monterey, Biomass_hyp_low_SoCal),
#     Biomass = pooled_log_mean(Biomass_Monterey, Biomass_SoCal),
#     `Biomass sd` = pooled_sd_mean(`Biomass sd_Monterey`, `Biomass sd_SoCal`, `Num Days_Monterey`, `Num Days_SoCal`),
#     BiomassTop50 = pooled_log_mean(BiomassTop50_Monterey, BiomassTop50_SoCal),
#     BiomassTop50sd = pooled_sd_mean(`BiomassTop50sd_Monterey`, `BiomassTop50sd_SoCal`, `Num Days_Monterey`, `Num Days_SoCal`),
#     Study_Area = NA
#     # 
#     # `Num Days` = coalesce(`Num Days_Monterey`, `Num Days_SoCal`),
#     # Biomass_hyp_low = coalesce(Biomass_hyp_low_Monterey, Biomass_hyp_low_SoCal),
#     # Biomass = coalesce(Biomass_Monterey, Biomass_SoCal),
#     # `Biomass sd` = coalesce(`Biomass sd_Monterey`, `Biomass sd_SoCal`),
#     # BiomassTop50 = coalesce(BiomassTop50_Monterey, BiomassTop50_SoCal),
#     # BiomassTop50sd = coalesce(`BiomassTop50sd_Monterey`, `BiomassTop50sd_SoCal`)
#   ) %>% 
#   select(-c(`Num Days_Monterey`:Biomass_hyp_low_SoCal)) %>% 
#   mutate_at(vars(Biomass_hyp_low:BiomassTop50sd) , ~log(10^.x))


krill_data_Scaling_paper <- read_csv("KrillData_Scaling_paper.csv") %>% 
  mutate(Biomass_hyp_low = ifelse(region == "Polar", NA, log(Biomass) + log(0.17))) %>% 
  rename(SpeciesCode = "Species") %>% 
  mutate_at(vars(Biomass, `Biomass sd`), log)


# WAP_krill_data <-   read_csv("AntarcticKrillData_10.21.20.csv") %>% 
#   #read_csv("AntarcticKrillDataDiveMeans.csv") %>% 
#   filter(!Species %in% c("bigBW")) %>% 
#   select(Species:BiomassTop50sd) %>%    # remember: biomass is in kg/m3
#   mutate(
#     Biomass_hyp_low = NA,
#     Study_Area = "Antarctic",
#     Region = "Antarctic") %>% 
#   rename(SpeciesCode = "Species")
# 
# # combining ENP prey data
# ENP_krill_data <- rbind(MRY_krill_data, SoCal_krill_data) %>% 
#   #mutate(Study_Area2 = Study_Area) %>% 
#   pivot_wider(names_from = Study_Area, values_from = c(`Num Days`:Biomass_hyp_low)) %>% 
#   #rename(Study_Area = Study_Area2) %>% 
#   mutate(
#     `Num Days` = `Num Days_Monterey` + `Num Days_SoCal`, 
#     Biomass_hyp_low = pooled_log_mean(Biomass_hyp_low_Monterey, Biomass_hyp_low_SoCal),
#     Biomass = pooled_log_mean(Biomass_Monterey, Biomass_SoCal),
#     `Biomass sd` = pooled_sd_mean(`Biomass sd_Monterey`, `Biomass sd_SoCal`, `Num Days_Monterey`, `Num Days_SoCal`),
#     BiomassTop50 = pooled_log_mean(BiomassTop50_Monterey, BiomassTop50_SoCal),
#     BiomassTop50sd = pooled_sd_mean(`BiomassTop50sd_Monterey`, `BiomassTop50sd_SoCal`, `Num Days_Monterey`, `Num Days_SoCal`),
#     Study_Area = NA
# 
#   ) %>% 
#   select(-c(`Num Days_Monterey`:Biomass_hyp_low_SoCal)) 
# 
# 
# All_krill_data <- rbind(MRY_krill_data, SoCal_krill_data, WAP_krill_data) 
# 
# All_krill_data_ENPcombined <- rbind(ENP_krill_data, WAP_krill_data)

# Whale population data---- 
#Current data from IUCN 2019 Redlist, whaling data compiled in Rocha et al. 2014
pop_data <- read_excel("Filtration project_Whale population and feeding information.xlsx", sheet = 1) %>%
  mutate(SpeciesCode = case_when(Species == "Balaenoptera musculus" ~ "bw",
                                 Species == "Balaenoptera physalus" ~ "bp",
                                 Species == "Megaptera novaeangliae" ~ "mn",
                                 Species == "Balaenoptera borealis" ~ "bs",
                                 Species == "Balaenoptera brydei" ~ "be",
                                 Species == "Balaenoptera acutorostrata" ~ "ba",
                                 Species == "Balaenoptera bonaerensis" ~ "bb",
                                 Species == "Eubalaena glacialis" ~ "eg",
                                 Species == "Eubalaena japonica" ~ "ej",
                                 Species == "Eubalaena australis" ~ "ea",
                                 Species == "Balaena mysticetus" ~ "bm"
  ))



# data projected for Antarctic blue and fin whales
Krill_data_Ant_projection <- read_excel("AntarcticKrillData_10.21.20.xlsx", sheet = 3) %>% 
  filter(Species %in% c("bp","bw")) %>% 
  rename(SpeciesCode = "Species") %>% 
  mutate(region = "Polar")

# All_krill_data_ENPcombined_noSA <- rbind(ENP_krill_data, WAP_krill_data) %>% 
#   mutate(region = ifelse(Region == "Antarctic",  "Polar","Temperate")) %>% 
#   filter(!SpeciesCode %in% c("bp","bw") | region != "Polar") %>% 
#   bind_rows(Krill_data_Ant_projection) %>%        # bringing in projected data for Antarctic blue and fin whales
#   select(-c(Region, Sv:SvTop50sd)) %>%
#   mutate(Study_Area = case_when(region == "Temperate" ~ "Eastern North Pacific",
#                                 region == "Polar" ~ "Antarctic")) %>% 
#   mutate_at(vars(Biomass_hyp_low:BiomassTop50sd) , ~log(.x))  # natural log transformed, which is what we were doing in the prior iteration


# Generate daily rates----
# combine day and twilight
combine_rates <- function(df, keys) {
  day_lunges <- df$Lunges[df$Phase == "Day"]
  twilight_lunges <- df$Lunges[df$Phase == "Twilight"]
  night_lunges <- df$Lunges[df$Phase == "Night"]
  tibble(
    day_lunges = sum(df$Lunges[df$Phase %in% c("Day", "Twilight")]),
    day_dur = sum(df$`Time (hours)`[df$Phase %in% c("Day", "Twilight")]),
    day_rate = day_lunges / day_dur,
    night_lunges = df$Lunges[df$Phase == "Night"],
    night_dur = df$`Time (hours)`[df$Phase == "Night"],
    night_rate = df$Rate[df$Phase == "Night"]
  )
}
daynight_rates <- All_rorqual_deployments %>% 
  left_join(dplyr::select(tag_guide, ID), by = "ID") %>% 
  mutate(
    region = case_when(
      Study_Area == "Antarctic" ~ "Polar",
      TRUE ~ "Temperate"
    )
  ) %>% 
  group_by(ID, region, prey_general, SpeciesCode) %>% 
  group_modify(combine_rates) %>% 
  ungroup() %>% 
  mutate(prey_general = ifelse(is.na(prey_general), "Krill", prey_general)) 

#summary table of deployments by species and regions
ss_table_dn <- daynight_rates %>%
  group_by(SpeciesCode, region) %>%
  summarise(sample_size = n())






# Weighting function for tag duration in daily feeding rate model
dur_weight <- function(dur) {
  # Coefficients for asymptotic model. Tells us the error for daily feeding
  # rate estimates from short-duration tags.
  Asym <- 0.7348501
  R0 <- 6.3862214
  lrc <- -1.7281635
  
  # Any duration longer than 10 hours has the same weight
  dur[dur > 10] <- 10
  
  # Weight duration by the inverse of the predicted error
  pred_err <- Asym + (R0 - Asym) * exp(-exp(lrc) * dur)
  1 / pred_err
}
plot(dur_weight, 
     xlim = c(0, 24), 
     ylim = c(0, 0.6),
     xlab = "Tag duration (hours)",
     ylab = "weight")

dev.copy2pdf(file="weighting fuction.pdf", width=4.5, height=4)

# Generate feeding rate distributions for each species
iter <- 1000
estimate_rate <- function(prey, min_lunges, min_dur) {
  daynight_rates %>%
    filter(prey_general == prey,
           day_lunges + night_lunges >= min_lunges,
           day_dur + night_dur >= min_dur) %>% 
    group_by(SpeciesCode, region) %>% 
    group_modify(
      function(data, keys) {
        tibble(i = 1:iter) %>% 
          mutate(
            day_rate = sample(
              data$day_rate, 
              size = iter,
              replace = TRUE,
              prob = dur_weight(data$day_dur)
            ),
            night_rate = sample(
              data$night_rate[data$night_dur > 0], 
              size = iter,
              replace = TRUE,
              prob = data$night_dur[data$night_dur > 0]
            )
          )
      }
    )
}
krill_rate_estimates <- estimate_rate("Krill", 1, 1)

plot(dur, pred_err, type = "l")


dur_sd <- function(dur) {
  Asym <- -1.846
  R0 <- 12.67
  lrc <- -2.057
  
  pred_err <- Asym + (R0 - Asym) * exp(-exp(lrc) * dur)
  pred_err_10 <- Asym + (R0 - Asym) * exp(-exp(lrc) * 10)
  ifelse(dur > 10, pred_err_10, pred_err)/2 
}


# Krill-feeding consumption estimates ----
krill_biomass_estimates <- krill_rate_estimates %>% 
  group_by(SpeciesCode) %>% 
  mutate(Engulf_cap_m3 = sample(whale_lengths$Engulfment_m3[whale_lengths$SpeciesCode == SpeciesCode[1]],
                                size = n(),
                                replace = TRUE),
         Mass_est_t = sample(whale_lengths$Lockyer_mass_t[whale_lengths$SpeciesCode == SpeciesCode[1]],
                             size = n(),
                             replace = TRUE),
         season_center = ifelse(
           region == "Temperate",
           213, # August 1 in northern hemisphere
           30 # January 30 in southern hemisphere
         )) %>% 
  ungroup %>% 
  left_join(krill_data_Scaling_paper, by = c("SpeciesCode", "region"))



# quick look at species average rates; Assuming Sep 15, 12.5 hours of daylight
krill_rate_estimates %>% 
  filter(SpeciesCode == "bp") %>% 
  mutate(daily_rate = day_rate * 12.5 + night_rate * 11.5) %>% 
  pull(daily_rate) %>% 
  quantile(c(0.25, 0.5, 0.75, 0.95))



latitudes <- tribble(
  ~ region,    ~ latitude,
  "Polar",     65,
  "Temperate", 36
)

# Overall daily rate, Modified by MATT
estimate_daily <- function(rate_estimates, latitude, season_len) {
  yday_start <- floor(yday_center - season_len / 2)
  yday_end <- floor(yday_center + season_len / 2)
  year_start <- as.POSIXct("2019-12-31", tz = "UTC")
  feeding_days <- year_start + lubridate::days(yday_start:yday_end)
  
  daylengths <- crossing(day = feeding_days, latitude) %>% 
    mutate(
      sunrise = maptools::sunriset(cbind(0, latitude), day, direction = "sunrise"),
      sunset = maptools::sunriset(cbind(0, latitude), day, direction = "sunset"),
      daylen = 24 * (sunset - sunrise),
      nightlen = 24 - daylen,
      yday = lubridate::yday(day)
    )
  
  # This function allows for NAs in logmean (i.e. Antarctic hypothetical low) and truncation
  # to the upper 50% (i.e. to account for selectivity)
  get_daily_mean <- function(logmean, logsd, lunges, trunc) {
    rand <- function(n, meanlog, sdlog) {
      if (trunc) {
        # Truncate values below median
        EnvStats::rlnormTrunc(n, meanlog, sdlog, min = exp(logmean))
      } else {
        rlnorm(n, meanlog, sdlog)
      }
    }
    result <- pmap_dbl(list(lunges, 
                            ifelse(is.na(logmean), 0, logmean), 
                            logsd), 
                       ~ mean(rand(floor(..1), ..2, ..3)))
    result <- result[is.na(logmean)] <- NA
    result
  }
  
  # MATT: change rate_estimates to biomass_estimates
  krill_biomass_estimates %>% 
    left_join(daylengths, by = "region") %>% 
    # MATT: change daily_rate to daily_biomass
    group_by(SpeciesCode, region) %>% 
    # mutate(biomass_best_low_m3 = rlnorm(n(), Biomass[1], `Biomass sd`[1]),
    #        biomass_best_high_m3 = rlnorm(n(), BiomassTop50[1], BiomassTop50sd[1]),
    mutate(daily_rate = day_rate * daylen + night_rate * nightlen) %>% 
    ungroup() %>% 
    mutate(
      daily_mean_hyp_low_kg = get_daily_mean(Biomass_hyp_low, `Biomass sd`, daily_rate, trunc = FALSE),
      daily_mean_kg = get_daily_mean(Biomass, `Biomass sd`, daily_rate, trunc = FALSE),
      daily_mean_hyp_high_kg = get_daily_mean(Biomass, `Biomass sd`, daily_rate, trunc = TRUE),
      
      daily_consumption_hyp_low_kg = daily_mean_hyp_low_kg * Engulf_cap_m3 * daily_rate,
      daily_consumption_kg = daily_mean_kg * Engulf_cap_m3 * daily_rate, 
      daily_consumption_hyp_high_kg = daily_mean_hyp_high_kg * Engulf_cap_m3 * daily_rate)
}

# MATT: Change krill_rate_estimates to krill_biomass_estimates. Should have
# columns day_biomass and night_biomass.


krill_daily <- estimate_daily(krill_biomass_estimates, latitudes, peak_day, 120)  %>% 
  mutate(Total_energy_intake_best_low_kJ = case_when(region == "Polar" ~ daily_consumption_low_kg*4575,
                                                     region == "Temperate" ~ daily_consumption_low_kg*3628),
         Total_energy_intake_best_high_kJ = case_when(region == "Polar" ~ daily_consumption_high_kg*4575,
                                                      region == "Temperate" ~ daily_consumption_high_kg*3628),
         Mass_specifc_energy_intake_best_high_kJ = Total_energy_intake_best_high_kJ/(Mass_est_t*1000),
         Mass_specifc_energy_intake_best_low_kJ = Total_energy_intake_best_low_kJ/(Mass_est_t*1000),
         Species = case_when(
           SpeciesCode == "bw" ~ "Balaenoptera musculus",
           SpeciesCode == "bp" ~ "Balaenoptera physalus",
           SpeciesCode == "mn" ~ "Megaptera novaeangliae",
           SpeciesCode == "bb" ~ "Balaenoptera bonaerensis", 
           SpeciesCode == "be" ~ "Balaenoptera brydei",
           SpeciesCode == "bs" ~ "Balaenoptera borealis")
  )


#save(krill_daily, file = "daynights_krill.RData") # THIS SAVES DATA WITH SA, Chile, ETC
#load("daynights_krill.RData") # THIS LOADS DATA WITH SA, Chile, ETC

#save(krill_daily, file = "daynights_krill_v2.RData") # THIS SAVES DATA WITHOUT SA, Chile, Azores, Norway 
#load("daynights_krill_v2.RData") # THIS LOADS DATA WITHOUT SA, Chile, Azores, Norway 

#save(krill_daily, file = "daynights_krill_v3.RData") # THIS SAVES *DIVEMEANS* DATA From DEC on 10.20.20 DATA; WITHOUT SA, Chile, Azores, Norway 
#load("daynights_krill_v3.RData") # This loads the DIVEMEANS data from DEC 10.20.20


# check to see prey consumption rates 
krill_daily %>%  
  group_by(SpeciesCode) %>% 
  filter(SpeciesCode == "bw", region == "Temperate") %>%
  pull(daily_consumption_high_kg) %>%
  summary()





# Fish feeding consumption estimates ----

fish_rate_estimates <- estimate_rate("Fish", 1, 1)


fish_biomass_estimates <- fish_rate_estimates %>% 
  group_by(SpeciesCode) %>% 
  mutate(Engulf_cap_m3 = sample(whale_lengths$Engulfment_m3[whale_lengths$SpeciesCode == SpeciesCode[1]],
                                size = n(),
                                replace = TRUE),
         Mass_est_t = sample(whale_lengths$Lockyer_mass_t[whale_lengths$SpeciesCode == SpeciesCode[1]],
                             size = n(),
                             replace = TRUE))





latitudes <- tribble(
  ~ region,    ~ latitude,
  "Polar",     65,
  "Temperate", 36
)

# Overall daily rate, Modified by MATT
estimate_daily <- function(rate_estimates, latitude, yday_center, season_len) {
  yday_start <- floor(yday_center - season_len / 2)
  yday_end <- floor(yday_center + season_len / 2)
  year_start <- as.POSIXct("2019-12-31", tz = "UTC")
  feeding_days <- year_start + lubridate::days(yday_start:yday_end)
  
  daylengths <- crossing(day = feeding_days, latitude) %>% 
    mutate(
      sunrise = maptools::sunriset(cbind(0, latitude), day, direction = "sunrise"),
      sunset = maptools::sunriset(cbind(0, latitude), day, direction = "sunset"),
      daylen = 24 * (sunset - sunrise),
      nightlen = 24 - daylen,
      yday = lubridate::yday(day)
    )
  
  #browser()
  # MATT: change rate_estimates to biomass_estimates
  daily_rates_fish <- fish_biomass_estimates %>% 
    left_join(daylengths, by = "region") %>% 
    # MATT: change daily_rate to daily_biomass
    group_by(SpeciesCode, region) %>% 
    mutate(
      # fish_biomass_best_high_m3 = rnorm(n(), 0.375, 0.243),
      #      fish_biomass_best_low_m3 = rnorm(n(), 0.375, 0.243)*0.29,
      daily_rate = day_rate * daylen + night_rate * nightlen) %>% 
    ungroup() %>% 
    mutate(
      daily_mean_low_kg = mean((rnorm(floor(daily_rate), 0.375, 0.243))*0.29)*7.8,
      daily_mean_high_kg = mean(rnorm(floor(daily_rate), 0.375, 0.243))*7.8,
      daily_consumption_low_kg = daily_mean_low_kg*Engulf_cap_m3*daily_rate, 
      daily_consumption_high_kg = daily_mean_high_kg*Engulf_cap_m3*daily_rate) 
  
}

# MATT: Change krill_rate_estimates to krill_biomass_estimates. Should have
# columns day_biomass and night_biomass.

aug1 <- 213
fish_daily <- estimate_daily(daily_rates_fish, latitudes, aug1, 120) %>% 
  mutate(Total_energy_intake_best_high_kJ = daily_consumption_high_kg*6000,
         Total_energy_intake_best_low_kJ = daily_consumption_low_kg*6000,
         Mass_specifc_energy_intake_best_high_kJ = Total_energy_intake_best_high_kJ/(Mass_est_t*1000),
         Mass_specifc_energy_intake_best_low_kJ = Total_energy_intake_best_low_kJ/(Mass_est_t*1000),
         Species = case_when(
           SpeciesCode == "mn" ~ "Megaptera novaeangliae",
           SpeciesCode == "be" ~ "Balaenoptera brydei")
  )


save(fish_daily, file = "daynights_fish.RData")
load("daynights_fish.RData")



fish_daily %>%  
  group_by(SpeciesCode) %>% 
  filter(SpeciesCode == "mn") %>%
  pull(daily_consumption_low_kg) %>%
  summary()



# FIGURES ----

##Figure 2----
#summary tables of lunges, water filtered, and prey consumed per day 
summ_stats_all <- krill_daily %>%    # toggle between fish and krill
  group_by(SpeciesCode, region) %>% 
  summarise(
    `Lunges per day` = median(daily_rate),
    Daily_rate_IQR25 = round(quantile(daily_rate, probs = 0.25, na.rm = TRUE), 2), 
    Daily_rate_IQR75 = round(quantile(daily_rate, probs = 0.75, na.rm = TRUE), 2),
    `Water filtered per day` = mean(Engulf_cap_m3*daily_rate),
    SE_filt = sd(Engulf_cap_m3*daily_rate)/sqrt(110),
    Water_filtered_IQR25 = round(quantile(Engulf_cap_m3*daily_rate, probs = 0.25, na.rm = TRUE), 2), 
    Water_filtered_IQR75 = round(quantile(Engulf_cap_m3*daily_rate, probs = 0.75, na.rm = TRUE), 2)) %>% 
  unite("Lunges per day IQR", c(Daily_rate_IQR25, Daily_rate_IQR75), sep = "-") %>% 
  unite("Water filtered per day IQR", c(Water_filtered_IQR25, Water_filtered_IQR75), sep = "-")

summ_prey_stats <- fish_daily %>%  # toggle between fish and krill
  #filter(!Region %in% c("North Atlantic", "Chile")) %>% 
  group_by(SpeciesCode, region) %>% 
  summarise(
    med_daily_consumpt_best = round(median(daily_consumption_low_kg, na.rm = TRUE), 2),
    med_daily_consumpt_best_IQR25 = round(quantile(daily_consumption_low_kg, probs = 0.25, na.rm = TRUE), 2), 
    med_daily_consumpt_best_IQR75 = round(quantile(daily_consumption_low_kg, probs = 0.75, na.rm = TRUE), 2), 
    med_daily_consumpt_top50 = round(median(daily_consumption_high_kg, na.rm = TRUE), 2),
    med_daily_consumpt_top50_IQR25 = round(quantile(daily_consumption_high_kg, probs = 0.25, na.rm = TRUE), 2), 
    med_daily_consumpt_top50_IQR75 = round(quantile(daily_consumption_high_kg, probs = 0.75, na.rm = TRUE), 2)) %>% 
  #unite("Daily consumption lower estimate IQR", c(med_daily_consumpt_hyp_low_IQR25, med_daily_consumpt_hyp_low_IQR75), sep = "-") %>% 
  unite("Daily consumption IQR", c(med_daily_consumpt_best_IQR25, med_daily_consumpt_best_IQR75), sep = "-") %>% 
  unite("Daily consumption Top 50% IQR", c(med_daily_consumpt_top50_IQR25, med_daily_consumpt_top50_IQR75), sep = "-")


summ_EnDens <- krill_daily %>% 
  #filter(!Region %in% c("North Atlantic", "Chile")) %>% 
  group_by(SpeciesCode, region) %>% 
  # Daily energy intake in gigajoules (GJ)
  summarise(
    med_daily_En_best = round(median(Total_energy_intake_best_low_kJ/1e6, na.rm = TRUE), 2),
    med_daily_En_best_IQR25 = round(quantile(Total_energy_intake_best_low_kJ/1e6, probs = 0.25, na.rm = TRUE), 2), 
    med_daily_En_best_IQR75 = round(quantile(Total_energy_intake_best_low_kJ/1e6, probs = 0.75, na.rm = TRUE), 2), 
    med_daily_mass_specific_En_best = round(median(Mass_specifc_energy_intake_best_low_kJ, na.rm = TRUE), 2),
    med_daily_En_top50 = round(median(Total_energy_intake_best_high_kJ/1e6, na.rm = TRUE), 2),
    med_daily_En_top50_IQR25 = round(quantile(Total_energy_intake_best_high_kJ/1e6, probs = 0.25, na.rm = TRUE), 2), 
    med_daily_En_top50_IQR75 = round(quantile(Total_energy_intake_best_high_kJ/1e6, probs = 0.75, na.rm = TRUE), 2),
    med_daily_mass_specific_En_high = round(median(Mass_specifc_energy_intake_best_high_kJ, na.rm = TRUE), 2)) %>% 
  #unite("Daily energy intake (GJ) lower estimate IQR", c(med_daily_En_hyp_low_IQR25, med_daily_En_hyp_low_IQR75), sep = "-") %>% 
  unite("Daily energy intake (GJ) IQR", c(med_daily_En_best_IQR25, med_daily_En_best_IQR75), sep = "-") %>% 
  unite("Daily energy intake (GJ) Top 50% IQR", c(med_daily_En_top50_IQR25, med_daily_En_top50_IQR75), sep = "-")


# Figure 2 Daily rorqual water filtration & krill consumption per individual ----

#ANTARCTIC
Daily_rate_Antarctic <- krill_daily %>% 
  filter(region == "Polar",
         daily_rate >75) %>% 
           
  ggplot(aes(x = fct_relevel(abbr_binom(Species), "B. bonaerensis"), y = daily_rate, 
             fill = abbr_binom(Species))) +
  geom_flat_violin(position = position_nudge(x = 0.1, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  #coord_flip() + 
  scale_fill_manual(values = pal) +
  labs(x = "",
       y = bquote(atop('Calculated feeding rate',
                       (lunges~ind^-1~d^-1)))) +
  #ylim(0, 2000) +
  scale_y_log10(labels = scales::comma, limits = c(25, 3000), 
                breaks = c(50,100,500,1000,2000)) +
  annotation_logticks(sides = "l") +
  
  theme_classic(base_size = 20) +
  theme(legend.position = "none",
        #axis.text.x = element_text(face = "italic", size = 16),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
        #axis.text.y = element_text(angle = 45, hjust = 1))
Daily_rate_Antarctic

dev.copy2pdf(file="Daily_rate_Antarctic_portrait.pdf", width=5.5, height=4.5)


Daily_filtration_Antarctic <- krill_daily %>% 
  filter(region == "Polar") %>%   
  ggplot(aes(x = fct_reorder(abbr_binom(Species), Engulf_cap_m3*daily_rate), y = Engulf_cap_m3*daily_rate, 
             fill = abbr_binom(Species))) +
  geom_flat_violin(position = position_nudge(x = 0.1, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  #coord_flip() + 
  scale_fill_manual(values = pal) +
  scale_y_log10(labels = scales::comma, 
                limits = c(100,75000), 
                breaks = c(100,500,1000,5000,10000,50000)) +
  annotation_logticks(sides = "l") +
  labs(x = "",
       y = bquote(atop('Calculated water filtered',
                       (m^3~ind^-1~d^-1)))) + 
  theme_classic(base_size = 20) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
        
Daily_filtration_Antarctic

dev.copy2pdf(file="Daily_filtration_Antarctic_portrait.pdf", width=6, height=4)



Daily_biomass_ingested_Antarctic <- krill_daily %>% 
  filter(region == "Polar") %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. bonaerensis")) %>% 
  
  # #best-low *our conservative estimate*
  ggplot(aes(x = Species, y = daily_consumption_low_kg,
             fill = abbr_binom(Species))) +
  geom_flat_violin(position = position_nudge(x = 0.1, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +

  scale_fill_manual(values = pal) +
  scale_y_log10(labels = scales::comma, limits = c(10, 40000), 
                breaks = c(0,100,250,1000,2500,10000,25000)) +
  annotation_logticks(sides = "l") +
  labs(x = "",
       y = bquote(atop('Calculated krill consumed',
                       (kg~ind^-1~d^-1)))) +
   theme_classic(base_size = 20) +
   theme(legend.position = "none",
         axis.text.x = element_text(face = "italic", size = 16))
        
Daily_biomass_ingested_Antarctic

dev.copy2pdf(file="Daily_biomass_ingested_Antarctic_portrait.pdf", width=6, height=4)


Fig_2_Antarctic <- ggarrange(Daily_rate_Antarctic, Daily_filtration_Antarctic, Daily_biomass_ingested_Antarctic, 
                             labels = c("A", "C", "E"), 
                             font.label = list(size = 18),
                             widths = c(1, 1), heights = c(2, 2, 2.2),
                             legend = "none",
                             ncol = 1, nrow = 3)
Fig_2_Antarctic

dev.copy2pdf(file="Fig_2_Antarctic_portrait.pdf", width=6, height=14)




#NON-ANTARCTIC
Daily_rate_nonAntarctic <- Daily_rate_Antarctic <- krill_daily %>% 
  filter(region == "Temperate") %>% 
  
  ggplot(aes(x = fct_relevel(abbr_binom(Species), "M. novaeangliae", "B. physalus"), 
             y = daily_rate, 
             fill = abbr_binom(Species))) +
  geom_flat_violin(position = position_nudge(x = 0.1, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  #coord_flip() + 
  scale_fill_manual(values = pal) +
  labs(x = "",
       y = "") +
  scale_y_log10(labels = scales::comma, limits = c(25, 3000), 
                breaks = c(50,100,500,1000,2000)) +
  annotation_logticks(sides = "l") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

Daily_rate_nonAntarctic

dev.copy2pdf(file="Daily_rate_nonAntarctic.pdf", width=7.25, height=4)


Daily_filtration_nonAntarctic <- krill_daily %>% 
  filter(region == "Temperate") %>%   
  ggplot(aes(x = fct_relevel(abbr_binom(Species), "M. novaeangliae", "B. physalus"), 
             y = Engulf_cap_m3*daily_rate, 
             fill = abbr_binom(Species))) +
  geom_flat_violin(position = position_nudge(x = 0.1, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  #coord_flip() + 
  scale_fill_manual(values = pal) +
  scale_y_log10(labels = scales::comma, 
                limits = c(100,75000), 
                breaks = c(100,500,1000,5000,10000,50000)) +
  annotation_logticks(sides = "l") +
  labs(x = "",
       y = "") + 
  theme_classic(base_size = 20) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
Daily_filtration_nonAntarctic

dev.copy2pdf(file="Daily_filtration_nonAntarctic.pdf", width=6, height=4)



Daily_biomass_ingested_nonAntarctic <- krill_daily %>% 
  filter(region == "Temperate") %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "M. novaeangliae", "B. physalus")) %>% 
  
  # #best-low *our conservative estimate*
  ggplot(aes(x = Species, y = daily_consumption_low_kg,
             fill = abbr_binom(Species))) +
  geom_flat_violin(position = position_nudge(x = 0.1, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  #coord_flip() + 
  
  scale_fill_manual(values = pal) +
  scale_y_log10(labels = scales::comma, limits = c(10, 40000), 
                breaks = c(0,100,250,1000,2500,10000,25000)) +
  annotation_logticks(sides = "l") +
  labs(x = "",
       y = "") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none",
        axis.text.x = element_text(face = "italic", size = 16))
Daily_biomass_ingested_nonAntarctic

dev.copy2pdf(file="Daily_biomass_ingested_nonAntarctic.pdf", width=6, height=4)


Fig_2_nonAntarctic <- ggarrange(Daily_rate_nonAntarctic, Daily_filtration_nonAntarctic, Daily_biomass_ingested_nonAntarctic, 
                                labels = c("B", "D", "F"), 
                                font.label = list(size = 18),
                                #widths = c(1.26, 1, 1),
                                legend = "none",
                                ncol = 1, nrow = 3)
Fig_2_nonAntarctic

dev.copy2pdf(file="Fig_2_nonAntarctic.pdf", width=6.5, height=13.5)


# making a full Figure 2
Fig_2_full <- ggarrange(Daily_rate_Antarctic, Daily_rate_nonAntarctic,  
                        Daily_filtration_Antarctic, Daily_filtration_nonAntarctic, 
                        Daily_biomass_ingested_Antarctic, Daily_biomass_ingested_nonAntarctic, 
                        labels = c("A","B","C","D","E","F"), 
                        font.label = list(size = 18),
                        widths = c(0.9, 1.1),
                        align = "hv",
                        legend = "none",
                        ncol = 2, nrow = 3)
Fig_2_full

 
dev.copy2pdf(file="Fig_2_full.pdf", width=12.5, height=14)



# Extended Figure 2 FISH-FEEDERS & hyp_low and Targeted feeding krill consumption ----
Daily_rate_fish <- fish_daily %>% 

  ggplot(aes(x = fct_relevel(abbr_binom(Species), "B. brydei"), y = daily_rate, 
             fill = abbr_binom(Species))) +
 
  geom_boxplot(position = position_nudge(x = 0, y = 0),
               width = .15, guides = FALSE, outlier.shape = NA, alpha = 0.8) +
  
  geom_flat_violin(position = position_nudge(x = 0.1, y = 0), alpha = 0.8) +
  
  scale_fill_manual(values = pal) +
  labs(x = "",
       y = bquote(atop('Estimated feeding rate',
                       ~(lunges~ind^-1~d^-1)))) + 
  scale_y_log10(labels = scales::comma, limits = c(5, 600),
                breaks = c(10,50,100,500)) +
  annotation_logticks(sides = "l") + 
  theme_classic(base_size = 20) +
  theme(legend.position = "none",        
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
  
Daily_rate_fish

dev.copy2pdf(file="Daily_rate_fish.pdf", width=9.5, height=4.75)


Daily_filtration_fish <- fish_daily %>% 
  filter(region == "Temperate") %>%   
  ggplot(aes(x = fct_reorder(abbr_binom(Species), Engulf_cap_m3*daily_rate), y = Engulf_cap_m3*daily_rate, 
             fill = abbr_binom(Species))) +
  geom_boxplot(position = position_nudge(x = 0, y = 0),
               width = .15, guides = FALSE, outlier.shape = NA, alpha = 0.8) +
  
  geom_flat_violin(position = position_nudge(x = 0.1, y = 0), alpha = 0.8) +

  scale_fill_manual(values = pal) +
  scale_y_log10(labels = scales::comma, 
                limits = c(10,20000), 
                breaks = c(50,100,500,1000,5000,10000)) +
  labs(x = "",
       y = bquote(atop('Estimated water filtered',
                       ~(m^3~ind^-1~d^-1)))) + 
  annotation_logticks(sides = "l") + 
  theme_classic(base_size = 20) +
  theme(legend.position = "none",        
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
Daily_filtration_fish

dev.copy2pdf(file="Daily_filtration_fish.pdf", width=9.5, height=4.75)



Daily_biomass_ingested_fish <- fish_daily %>% 
  ggplot() +
  
  # #best-low *our best estimate*
  geom_flat_violin(aes(x = abbr_binom(Species), y = daily_consumption_low_kg,
        fill = abbr_binom(Species)), 
    position = position_nudge(x = 0.1, y = 0), alpha = 0.6) +
  
  # #best-low *our best estimate*
  geom_boxplot(aes(x = abbr_binom(Species), y = daily_consumption_low_kg,
                   fill = abbr_binom(Species)), 
               position = position_nudge(x = -0.1, y = 0),
               width = .15, guides = FALSE, outlier.shape = NA, alpha = 0.8) +
  
  
  # best-upper, MOVED TO SUPPLEMENTAL
  geom_boxplot(aes(x = abbr_binom(Species), y = daily_consumption_high_kg,
                   fill = abbr_binom(Species)),
               position = position_nudge(x = -0.25, y = 0),
               width = .15, guides = FALSE, outlier.shape = NA, alpha = 0.8) +
  
# best-upper, MOVED TO SUPPLEMENTAL
  geom_flat_violin(aes(x = abbr_binom(Species), y = daily_consumption_high_kg,
        fill = abbr_binom(Species)),
  position = position_nudge(x = 0.1, y = 0), alpha = 0.6) +

  scale_fill_manual(values = pal) +
  scale_y_log10(labels = scales::comma, 
                limits = c(10,25000), 
                breaks = c(50,100,500,1000,5000,10000)) +
  labs(x = "Species",
       y = bquote(atop('Estimated fish consumed',
                       ~(kg~ind^-1~d^-1)))) +

  annotation_logticks(sides = "l") + 
  theme_classic(base_size = 20) +
  theme(legend.position = "none",
        axis.text.x = element_text(face = "italic", size = 16))

Daily_biomass_ingested_fish

dev.copy2pdf(file="Daily_biomass_ingested_fish.pdf", width=9.5, height=4.75)



# Extended Figure 2----
Fig_2_Extended_fish <- ggarrange(Daily_rate_fish, Daily_filtration_fish, Daily_biomass_ingested_fish,
                        labels = c("C","D","E"), 
                        font.label = list(size = 18),
                        align = "hv",
                        legend = "none",
                        ncol = 1, nrow = 3)
Fig_2_Extended_fish

dev.copy2pdf(file="Fig_2_Extended_fish.pdf", width=6.5, height=14)




Daily_biomass_ingested_AntarcticTOP50 <- krill_daily %>% 
  filter(region == "Polar") %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. bonaerensis")) %>% 
  
  ggplot(aes(fill = abbr_binom(Species))) +
  
  # best-upper, MOVED TO SUPPLEMENTAL
  geom_flat_violin(aes(x = Species, y = daily_consumption_high_kg),
                   position = position_nudge(x = 0.1, y = 0), alpha = 0.8) +
  geom_boxplot(aes(x = Species, y = daily_consumption_high_kg),
               width = .15, guides = FALSE, outlier.shape = NA, alpha = 0.8) +

  scale_fill_manual(values = pal) +
  scale_y_log10(labels = scales::comma, limits = c(50, 12000), 
                breaks = c(0,100,250,500,1000,2500,5000,10000)) +
  labs(x = "",
       y = bquote(atop('Calculated krill consumed',
                       ~(kg~ind^-1~d^-1)))) +
  annotation_logticks(sides = "l") + 
  theme_classic(base_size = 20) +
  theme(legend.position = "none",
        axis.text.x = element_text(face = "italic", size = 16))
Daily_biomass_ingested_AntarcticTOP50

dev.copy2pdf(file="Daily_biomass_ingested_AntarcticTOP50.pdf", width=9.5, height=4.75)


Daily_biomass_ingested_nonAntarcticTOP50 <- krill_daily %>% 
  filter(region == "Temperate") %>% 
  
  ggplot(aes(fill = abbr_binom(Species))) +
  
  # # Hyp_low, MOVED TO SUPPLEMENTAL
  geom_flat_violin(
    aes(x = fct_relevel(abbr_binom(Species), "M. novaeangliae", "B. physalus"),
        y = daily_consumption_hyp_low_kg),
    position = position_nudge(x = 0.1, y = 0), alpha = 0.4) +
  geom_boxplot(aes(x = fct_relevel(abbr_binom(Species), "M. novaeangliae", "B. physalus"),
                   y = daily_consumption_hyp_low_kg), 
               position = position_nudge(x = -0.1, y = 0),
               width = .15, guides = FALSE, outlier.shape = NA, alpha = 0.4) +
  
  # best-upper, MOVED TO SUPPLEMENTAL
  geom_flat_violin(
    aes(x = fct_relevel(abbr_binom(Species), "M. novaeangliae", "B. physalus"),
        y = daily_consumption_high_kg),
        position = position_nudge(x = 0.1, y = 0), alpha = 0.8) +
  geom_boxplot(aes(x = fct_relevel(abbr_binom(Species), "M. novaeangliae", "B. physalus"),
                   y = daily_consumption_high_kg), 
               position = position_nudge(x = -0.25, y = 0),
               width = .15, guides = FALSE, outlier.shape = NA, alpha = 0.8) +

  
  scale_fill_manual(values = pal) +
  scale_y_log10(labels = scales::comma, limits = c(100, 50000), 
                breaks = c(500,1000,2500,5000,10000,25000, 50000)) +
  labs(x = "",
       y = bquote(atop('Calculated krill consumed',
                        ~(kg~ind^-1~d^-1)))) +
  
  annotation_logticks(sides = "l") + 
  theme_classic(base_size = 20) +
  theme(legend.position = "none",
        axis.text.x = element_text(face = "italic", size = 16))

Daily_biomass_ingested_nonAntarcticTOP50

dev.copy2pdf(file="Daily_biomass_ingested_nonAntarcticTOP50.pdf", 
             width=10, height=5.5)


Fig_2_Extended_krill <- ggarrange(Daily_biomass_ingested_AntarcticTOP50, Daily_biomass_ingested_nonAntarcticTOP50,
                                 labels = c("A","B"), 
                                 font.label = list(size = 18),
                                 align = "hv",
                                 legend = "none",
                                 ncol = 1, nrow = 2)
Fig_2_Extended_krill

dev.copy2pdf(file="Fig_2_Extended_krill.pdf", width=6.5, height=10)


#summary tables of individual water filtered and prey consumed per year (Figure 3)----
krill_yearly_ind <- krill_daily %>% 
  mutate(filtration60 = Engulf_cap_m3*daily_rate*60,
         filtration90 = Engulf_cap_m3*daily_rate*90,
         filtration120 = Engulf_cap_m3*daily_rate*120,
         filtration182.5 = Engulf_cap_m3*daily_rate*182.5,
         prey60 = daily_consumption_low_kg*60,
         prey90 = daily_consumption_low_kg*90,
         prey120 = daily_consumption_low_kg*120,
         prey182.5 = daily_consumption_low_kg*182.5)

# Annual individual filtration, summary table
summ_filt_annual_ind_stats <- krill_yearly_ind %>%
  group_by(SpeciesCode, region) %>% 
  summarise(
    yr_filt_60days_IQR25 = round(quantile(filtration60, probs = 0.25, na.rm = TRUE), 0),
    yr_filt_60days_IQR75 = round(quantile(filtration60, probs = 0.75, na.rm = TRUE), 0),
    yr_filt_90days_IQR25 = round(quantile(filtration90, probs = 0.25, na.rm = TRUE), 0),
    yr_filt_90days_IQR75 = round(quantile(filtration90, probs = 0.75, na.rm = TRUE), 0),
    yr_filt_120days_IQR25 = round(quantile(filtration120, probs = 0.25, na.rm = TRUE), 0),
    yr_filt_120days_IQR75 = round(quantile(filtration120, probs = 0.75, na.rm = TRUE), 0),
    yr_filt_182.5days_IQR25 = round(quantile(filtration182.5, probs = 0.25, na.rm = TRUE), 0),
    yr_filt_182.5days_IQR75 = round(quantile(filtration182.5, probs = 0.75, na.rm = TRUE), 0)) %>% 
  unite("Filtration capacity (m3 ind yr), 60 days feeding, IQR", 
        c(yr_filt_60days_IQR25, yr_filt_60days_IQR75), sep = "-") %>% 
  unite("Filtration capacity (m3 ind yr), 90 days feeding, IQR", 
        c(yr_filt_90days_IQR25, yr_filt_90days_IQR75), sep = "-") %>% 
  unite("Filtration capacity (m3 ind yr), 120 days feeding, IQR", 
        c(yr_filt_120days_IQR25, yr_filt_120days_IQR75), sep = "-") %>% 
  unite("Filtration capacity (m3 ind yr), 182.5 days feeding, IQR", 
        c(yr_filt_182.5days_IQR25, yr_filt_182.5days_IQR75), sep = "-")




# Annual individual krill prey, summary table
summ_prey_annual_ind_stats <- krill_yearly_ind %>% 
  group_by(SpeciesCode, region) %>% 
  summarise(
    yr_prey_60days_IQR25 = round(quantile(prey60, probs = 0.25, na.rm = TRUE)/1000, 0),
    yr_prey_60days_IQR75 = round(quantile(prey60, probs = 0.75, na.rm = TRUE)/1000, 0),
    yr_prey_90days_IQR25 = round(quantile(prey90, probs = 0.25, na.rm = TRUE)/1000, 0),
    yr_prey_90days_IQR75 = round(quantile(prey90, probs = 0.75, na.rm = TRUE)/1000, 0),
    yr_prey_120days_IQR25 = round(quantile(prey120, probs = 0.25, na.rm = TRUE)/1000, 0),
    yr_prey_120days_IQR75 = round(quantile(prey120, probs = 0.75, na.rm = TRUE)/1000, 0),
    yr_prey_182.5days_IQR25 = round(quantile(prey182.5, probs = 0.25, na.rm = TRUE)/1000, 0),
    yr_prey_182.5days_IQR75 = round(quantile(prey182.5, probs = 0.75, na.rm = TRUE)/1000, 0)) %>% 
  unite("Krill consumption (tonnes ind yr), 60 days feeding, IQR", 
        c(yr_prey_60days_IQR25, yr_prey_60days_IQR75), sep = "-") %>% 
  unite("Krill consumption (tonnes ind yr), 90 days feeding, IQR", 
        c(yr_prey_90days_IQR25, yr_prey_90days_IQR75), sep = "-") %>% 
  unite("Krill consumption (tonnes ind yr), 120 days feeding, IQR", 
        c(yr_prey_120days_IQR25, yr_prey_120days_IQR75), sep = "-") %>% 
  unite("Krill consumption (tonnes ind yr), 182.5 days feeding, IQR", 
        c(yr_prey_182.5days_IQR25, yr_prey_182.5days_IQR75), sep = "-")



#summary tables of population water filtered and prey consumed per year (Figure 4)----
krill_yearly_pop <- krill_daily %>% 
  left_join(pop_data, by = "SpeciesCode") %>% 
  mutate(filtration60curr = Engulf_cap_m3*daily_rate*60*`Population estimate (IUCN 2019)`,
         filtration90curr = Engulf_cap_m3*daily_rate*90*`Population estimate (IUCN 2019)`,
         filtration120curr = Engulf_cap_m3*daily_rate*120*`Population estimate (IUCN 2019)`,
         filtration182.5curr = Engulf_cap_m3*daily_rate*182.5*`Population estimate (IUCN 2019)`,
         filtration60hist = Engulf_cap_m3*daily_rate*60*`Historical estimate`,
         filtration90hist = Engulf_cap_m3*daily_rate*90*`Historical estimate`,
         filtration120hist = Engulf_cap_m3*daily_rate*120*`Historical estimate`,
         filtration182.5hist = Engulf_cap_m3*daily_rate*182.5*`Historical estimate`,
         #prey numbers for NORTHERN HEMISPHERE ONLY, SEE BELOW FOR SOUTHERN HEMISPHERE
         prey_low60curr = daily_consumption_low_kg*(`Population estimate (Christensen 2006)`-`Southern hemisphere population estimate (Christensen 2006)`)*60,
         prey_low90curr = daily_consumption_low_kg*(`Population estimate (Christensen 2006)`-`Southern hemisphere population estimate (Christensen 2006)`)*90,
         prey_low120curr = daily_consumption_low_kg*(`Population estimate (Christensen 2006)`-`Southern hemisphere population estimate (Christensen 2006)`)*120,
         prey_low182.5curr = daily_consumption_low_kg*(`Population estimate (Christensen 2006)`-`Southern hemisphere population estimate (Christensen 2006)`)*182.5,
         prey_low60hist = daily_consumption_low_kg*(`Historical estimate`-`Southern hemisphere historic estimate (Christensen 2006)`)*60,
         prey_low90hist = daily_consumption_low_kg*(`Historical estimate`-`Southern hemisphere historic estimate (Christensen 2006)`)*90,
         prey_low120hist = daily_consumption_low_kg*(`Historical estimate`-`Southern hemisphere historic estimate (Christensen 2006)`)*120,
         prey_low182.5hist = daily_consumption_low_kg*(`Historical estimate`-`Southern hemisphere historic estimate (Christensen 2006)`)*182.5,
         
         prey_high60curr = daily_consumption_high_kg*(`Population estimate (Christensen 2006)`-`Southern hemisphere population estimate (Christensen 2006)`)*60,
         prey_high90curr = daily_consumption_high_kg*(`Population estimate (Christensen 2006)`-`Southern hemisphere population estimate (Christensen 2006)`)*90,
         prey_high120curr = daily_consumption_high_kg*(`Population estimate (Christensen 2006)`-`Southern hemisphere population estimate (Christensen 2006)`)*120,
         prey_high182.5curr = daily_consumption_high_kg*(`Population estimate (Christensen 2006)`-`Southern hemisphere population estimate (Christensen 2006)`)*182.5,
         prey_high60hist = daily_consumption_high_kg*(`Historical estimate`-`Southern hemisphere historic estimate (Christensen 2006)`)*60,
         prey_high90hist = daily_consumption_high_kg*(`Historical estimate`-`Southern hemisphere historic estimate (Christensen 2006)`)*90,
         prey_high120hist = daily_consumption_high_kg*(`Historical estimate`-`Southern hemisphere historic estimate (Christensen 2006)`)*120,
         prey_high182.5hist = daily_consumption_high_kg*(`Historical estimate`-`Southern hemisphere historic estimate (Christensen 2006)`)*182.5)%>% 
  
  mutate_at(vars(c("prey_low60curr":"prey_high182.5hist")), ~case_when(SpeciesCode == "bp" ~ .*0.8,     # correcting for proportion of diet that is fish; and proportion of individuals not in Southern Hemisphere
                                                                       SpeciesCode == "bw" ~ .,
                                                                       SpeciesCode == "mn" ~ .*0.55,
                                                                       SpeciesCode == "bb" ~ .))



summ_filt_annual_pop_stats_curr <- krill_yearly_pop %>% 
  filter(!Region %in% c("North Atlantic", "Chile")) %>% 
  mutate(Species = abbr_binom(Species)) %>%
  group_by(Species) %>% 
  summarise(
    yr_filt_60days_IQR25 = round(quantile(filtration60curr, probs = 0.25, na.rm = TRUE)/1e9, 0),
    yr_filt_60days_IQR75 = round(quantile(filtration60curr, probs = 0.75, na.rm = TRUE)/1e9, 0),
    yr_filt_90days_IQR25 = round(quantile(filtration90curr, probs = 0.25, na.rm = TRUE)/1e9, 0),
    yr_filt_90days_IQR75 = round(quantile(filtration90curr, probs = 0.75, na.rm = TRUE)/1e9, 0),
    yr_filt_120days_IQR25 = round(quantile(filtration120curr, probs = 0.25, na.rm = TRUE)/1e9, 0),
    yr_filt_120days_IQR75 = round(quantile(filtration120curr, probs = 0.75, na.rm = TRUE)/1e9, 0),
    yr_filt_182.5days_IQR25 = round(quantile(filtration182.5curr, probs = 0.25, na.rm = TRUE)/1e9, 0),
    yr_filt_182.5days_IQR75 = round(quantile(filtration182.5curr, probs = 0.75, na.rm = TRUE)/1e9, 0)) %>% 
  unite("Filtration capacity (km3 ind yr), 60 days feeding, IQR", 
        c(yr_filt_60days_IQR25, yr_filt_60days_IQR75), sep = "-") %>% 
  unite("Filtration capacity (km3 ind yr), 90 days feeding, IQR", 
        c(yr_filt_90days_IQR25, yr_filt_90days_IQR75), sep = "-") %>% 
  unite("Filtration capacity (km3 ind yr), 120 days feeding, IQR", 
        c(yr_filt_120days_IQR25, yr_filt_120days_IQR75), sep = "-") %>% 
  unite("Filtration capacity (km3 ind yr), 150 days feeding, IQR", 
        c(yr_filt_182.5days_IQR25, yr_filt_182.5days_IQR75), sep = "-")


summ_filt_annual_pop_stats_hist <- krill_yearly_pop %>% 
  filter(!Region %in% c("North Atlantic", "Chile")) %>% 
  mutate(Species = abbr_binom(Species)) %>%
  group_by(Species) %>% 
  summarise(
    yr_filt_60days_IQR25 = round(quantile(filtration60hist, probs = 0.25, na.rm = TRUE)/1e9, 0),
    yr_filt_60days_IQR75 = round(quantile(filtration60hist, probs = 0.75, na.rm = TRUE)/1e9, 0),
    yr_filt_90days_IQR25 = round(quantile(filtration90hist, probs = 0.25, na.rm = TRUE)/1e9, 0),
    yr_filt_90days_IQR75 = round(quantile(filtration90hist, probs = 0.75, na.rm = TRUE)/1e9, 0),
    yr_filt_120days_IQR25 = round(quantile(filtration120hist, probs = 0.25, na.rm = TRUE)/1e9, 0),
    yr_filt_120days_IQR75 = round(quantile(filtration120hist, probs = 0.75, na.rm = TRUE)/1e9, 0),
    yr_filt_182.5days_IQR25 = round(quantile(filtration182.5hist, probs = 0.25, na.rm = TRUE)/1e9, 0),
    yr_filt_182.5days_IQR75 = round(quantile(filtration182.5hist, probs = 0.75, na.rm = TRUE)/1e9, 0)) %>% 
  unite("Filtration capacity (km3 ind yr), 60 days feeding, IQR", 
        c(yr_filt_60days_IQR25, yr_filt_60days_IQR75), sep = "-") %>% 
  unite("Filtration capacity (km3 ind yr), 90 days feeding, IQR", 
        c(yr_filt_90days_IQR25, yr_filt_90days_IQR75), sep = "-") %>% 
  unite("Filtration capacity (km3 ind yr), 120 days feeding, IQR", 
        c(yr_filt_120days_IQR25, yr_filt_120days_IQR75), sep = "-") %>% 
  unite("Filtration capacity (km3 ind yr), 150 days feeding, IQR", 
        c(yr_filt_182.5days_IQR25, yr_filt_182.5days_IQR75), sep = "-")

# Figure 3----

pal2 <- c("B. bonaerensis" = "firebrick3", "B. borealis" = "goldenrod2", "B. brydei" = "darkorchid3",  
          "M. novaeangliae" = "gray30", "B. physalus" = "chocolate3", "B. musculus" = "dodgerblue2",
          "α = 0.035;  β = 1" = "black", "α = 0.06;  β = 0.75" = "dark gray", "α = 0.1;  β = 0.8" = "dimgray", 
          "α = 0.123;  β = 0.8" = "salmon2", "α = 0.17;  β = 0.773" = "skyblue3", "α = 0.42;  β = 0.67" = "thistle",
          "α = 1.66;  β = 0.559" = "palegreen3", "Barlow et al. 2008" = "blue")


# preparing comparative data (previous studies)
days_feeding = rep(60:183, each = 5)

whale_mass <- whale_lengths %>% 
  group_by(Species) %>% 
  summarize(med_wt_kg = median(Lockyer_mass_t*1000, na.rm = TRUE))
  

# read in data for direct predictions of daily ration (R) from literature
prey_predict <- read_excel("PreyIngestPredict.xlsx") %>% mutate(dummy = 1)
prey_predict_w_M <- tibble(M_kg = seq(5000,120000,5000), dummy =1) %>%
  full_join(prey_predict, by = "dummy") %>% 
  select(-"dummy") %>% 
  mutate(R = `Intercept (α)`*M_kg^`Exponent (β)`,
         R_compressed_90days = (R*365)/90,
         R_compressed_100days = (R*365)/100,
         R_compressed_120days = (R*365)/120)

prey_predict_from_BMR <- read_excel("PreyIngestPredict.xlsx", sheet = 2) %>% 
  mutate(KleiberBMR = 293.1*M_kg^0.75,
         dummy = 1)
BMRtoFMRprojection <- tibble(beta = seq(2,5,0.5), dummy = 1) %>% 
  full_join(prey_predict_from_BMR, by = "dummy") %>% 
  select(-dummy) %>% 
  mutate(ADMR = beta*KleiberBMR,
         R  = ADMR/(0.8*((3900*Z)+5400*(1-Z))),
         R_compressed_90days = (0.83*(R*365))/90,
         R_compressed_100days = (0.83*(R*365))/100,
         R_compressed_120days = (0.83*(R*365))/120)



days_feeding = rep(60:183, each = 5)

whale_prey_yr_Eq1 <- BMRtoFMRprojection %>% 
  dplyr::select(beta, Species, Z) %>% 
  left_join(whale_mass, by = "Species") %>% 
  filter(Species != "Balaenoptera acutorostrata") %>% 
  mutate(
    Species = abbr_binom(Species),
    KleiberBMR = 293.1*med_wt_kg^0.75,
    ADMR = beta*KleiberBMR,
    R  = ADMR/(0.8*((3900*Z)+5400*(1-Z)))) %>% 
  crossing(days_feeding) %>% 
  mutate(MDC = 0.83* R * 365 / days_feeding,
         MYC_t = MDC * days_feeding/1000)

whale_prey_yr_Eq1_summ <- whale_prey_yr_Eq1 %>% 
  group_by(Species) %>% 
  filter(beta == 2:3) %>% 
  summarise(avg_yr_cons = mean(MYC_t),
            SE_yr_cons = SE(MYC_t)) %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), 
                        "B. bonaerensis", "B. physalus", "M. novaeangliae", "B. musculus"))

         


whale_prey_yr_Eq2 <- whale_mass %>% 
  crossing(select(prey_predict, -dummy)) %>% 
  crossing(days_feeding) %>% 
  mutate(
    Species = abbr_binom(Species),
    MDC = 0.83 * (`Intercept (α)`*med_wt_kg^`Exponent (β)`)*365 / days_feeding,   # MDC consumption assuming that many days spent intensively feeding; actually decreases with increasing # of days_feeding because the yearly R is fixed
    MYC_t = MDC * days_feeding/1000) 


whale_prey_yr_Eq2_summ <- whale_prey_yr_Eq2 %>% 
  group_by(Species, Parameters) %>% 
  summarise(avg_yr_cons = mean(MYC_t),
            SE_yr_cons = SE(MYC_t))
  

whale_prey_yr_prev <- whale_prey_yr_Eq2 %>% 
  bind_rows(
    dplyr::select(
      filter(whale_prey_yr_Eq1, beta == 2.5:3), 
      Species, med_wt_kg, days_feeding, MDC, MYC_t, beta)) %>% 
  mutate(`Reference(s)` = ifelse(is.na(`Reference(s)`), "Barlow et al. 2008", `Reference(s)`),
         Parameters = ifelse(is.na(Parameters), "Barlow et al. 2008", Parameters),
         Species = fct_relevel(factor(abbr_binom(Species)), 
                               "B. bonaerensis", "B. physalus", 
                               "M. novaeangliae", "B. musculus")) 
  

whale_prey_yr_prev_summ <- whale_prey_yr_Eq2_summ %>% 
  bind_rows(whale_prey_yr_Eq1_summ) %>% 
  mutate(Parameters = ifelse(is.na(Parameters), "Barlow et al. 2008", Parameters)) 


# Species plots of previous estimates
prey_yr_prev_summ_bb_plot <- whale_prey_yr_prev_summ %>% 
  filter(Species == "B. bonaerensis") %>% 
  ggplot() +
  geom_pointrange(aes(x=fct_reorder(Parameters, avg_yr_cons), y=avg_yr_cons, 
                      ymin=(avg_yr_cons-SE_yr_cons), ymax=(avg_yr_cons+SE_yr_cons),
                      color = Parameters),
                  size=0.75) +
  scale_colour_manual(values = pal2) +
  scale_y_continuous(limits = c(0,550),
                     breaks = c(0, 100, 200, 300, 400, 500)) +
  # scale_y_log10(labels = scales::comma, limits = c(10, 1000),
  #               breaks = c(30,100,300,600)) +
  xlab(" \n") +
  theme_classic(base_size = 20) +
  
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

prey_yr_prev_summ_bb_plot


prey_yr_prev_summ_mn_plot <- whale_prey_yr_prev_summ %>% 
  filter(Species == "M. novaeangliae") %>% 
  ggplot() +
  geom_pointrange(aes(x=fct_reorder(Parameters, avg_yr_cons), y=avg_yr_cons, 
                      ymin=(avg_yr_cons-SE_yr_cons), ymax=(avg_yr_cons+SE_yr_cons),
                      color = Parameters),
                  size=0.75) +
  scale_colour_manual(values = pal2) +
  scale_y_continuous(limits = c(0,550),
                     breaks = c(0, 100, 200, 300, 400, 500)) +
  # scale_y_log10(labels = scales::comma, limits = c(10, 1000),
  #               breaks = c(30,100,300,600)) +
  xlab(" \n") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())
prey_yr_prev_summ_mn_plot


prey_yr_prev_summ_bp_plot <- whale_prey_yr_prev_summ %>% 
  filter(Species == "B. physalus") %>% 
  ggplot() +
  geom_pointrange(aes(x=fct_reorder(Parameters, avg_yr_cons), y=avg_yr_cons, 
                      ymin=(avg_yr_cons-SE_yr_cons), ymax=(avg_yr_cons+SE_yr_cons),
                      color = Parameters),
                  size=0.75) +
  scale_colour_manual(values = pal2) +
  scale_y_continuous(limits = c(20, 2500),
                     breaks = c(100,500,1000,1500,
                                2000,2500)) +
  xlab(" \n") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())
prey_yr_prev_summ_bp_plot


prey_yr_prev_summ_bw_plot <- whale_prey_yr_prev_summ %>% 
  filter(Species == "B. musculus") %>% 
  ggplot() +
  geom_pointrange(aes(x=fct_reorder(Parameters, avg_yr_cons), y=avg_yr_cons, 
                      ymin=(avg_yr_cons-SE_yr_cons), ymax=(avg_yr_cons+SE_yr_cons),
                      color = Parameters),
                  size=0.75) +
  scale_colour_manual(values = pal2) +
  scale_y_continuous(limits = c(20, 2500),
                     breaks = c(100,500,1000,1500,
                                2000,2500)) +
  xlab(" \n") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())
prey_yr_prev_summ_bw_plot



days_feeding = rep(60:183, each = 5)
Annual_filtfeed <- krill_daily %>% 
  group_by(Species, region) %>% 
  summarise(
    #Individual
    med_filt = median(Engulf_cap_m3*daily_rate),
    IQR25_filt = quantile(Engulf_cap_m3*daily_rate, probs = 0.25, na.rm = TRUE),
    IQR75_filt = quantile(Engulf_cap_m3*daily_rate, probs = 0.75, na.rm = TRUE),
    med_ingest_low_t = median(daily_consumption_low_kg, na.rm = TRUE)/1000,
    IQR25_ingest_low_t = quantile(daily_consumption_low_kg, probs = 0.25, na.rm = TRUE)/1000,
    IQR75_ingest_low_t = quantile(daily_consumption_low_kg, probs = 0.75, na.rm = TRUE)/1000,
    med_ingest_high_t = median(daily_consumption_high_kg, na.rm = TRUE)/1000,
    IQR25_ingest_high_t = quantile(daily_consumption_high_kg, probs = 0.25, na.rm = TRUE)/1000,
    IQR75_ingest_high_t = quantile(daily_consumption_high_kg, probs = 0.75, na.rm = TRUE)/1000
  ) %>%
  crossing(days_feeding) %>% 
  mutate_at(c("med_filt", "IQR25_filt", "IQR75_filt"), funs(yr = .*days_feeding)) %>% 
  mutate_at(c("med_ingest_low_t", "IQR25_ingest_low_t", "IQR75_ingest_low_t"), funs(yr = .*days_feeding)) %>% 
  mutate_at(c("med_ingest_high_t", "IQR25_ingest_high_t", "IQR75_ingest_high_t"), funs(yr = .*days_feeding)) %>% 
  left_join(pop_data, by = "Species") %>% 
  mutate(
    
    #Current population
    med_filt_curr = med_filt_yr*(`Population estimate (Christensen 2006)`-`Southern hemisphere population estimate (Christensen 2006)`),
    IQR25_filt_curr = IQR25_filt_yr*(`Population estimate (Christensen 2006)`-`Southern hemisphere population estimate (Christensen 2006)`),
    IQR75_filt_curr = IQR75_filt_yr*(`Population estimate (Christensen 2006)`-`Southern hemisphere population estimate (Christensen 2006)`),
    med_ingest_low_t_curr = med_ingest_low_t_yr*(`Population estimate (Christensen 2006)`-`Southern hemisphere population estimate (Christensen 2006)`),
    IQR25_ingest_low_t_curr = IQR25_ingest_low_t_yr*(`Population estimate (Christensen 2006)`-`Southern hemisphere population estimate (Christensen 2006)`),
    IQR75_ingest_low_t_curr = IQR75_ingest_low_t_yr*(`Population estimate (Christensen 2006)`-`Southern hemisphere population estimate (Christensen 2006)`),
    med_ingest_high_t_curr = med_ingest_high_t_yr*(`Population estimate (Christensen 2006)`-`Southern hemisphere population estimate (Christensen 2006)`),
    IQR25_ingest_high_t_curr = IQR25_ingest_high_t_yr*(`Population estimate (Christensen 2006)`-`Southern hemisphere population estimate (Christensen 2006)`),
    IQR75_ingest_high_t_curr = IQR75_ingest_high_t_yr*(`Population estimate (Christensen 2006)`-`Southern hemisphere population estimate (Christensen 2006)`),
    
    #Historic population, NORTHERN HEMISPHERE ONLY, NEED TO FIX!!!
    med_filt_hist = med_filt_yr*(`Historical estimate`-`Southern hemisphere historic estimate (Christensen 2006)`),
    IQR25_filt_hist = IQR25_filt_yr*(`Historical estimate`-`Southern hemisphere historic estimate (Christensen 2006)`),
    IQR75_filt_hist = IQR75_filt_yr*(`Historical estimate`-`Southern hemisphere historic estimate (Christensen 2006)`),
    med_ingest_low_t_hist = med_ingest_low_t_yr*(`Historical estimate`-`Southern hemisphere historic estimate (Christensen 2006)`),
    IQR25_ingest_low_t_hist = IQR25_ingest_low_t_yr*(`Historical estimate`-`Southern hemisphere historic estimate (Christensen 2006)`),
    IQR75_ingest_low_t_hist = IQR75_ingest_low_t_yr*(`Historical estimate`-`Southern hemisphere historic estimate (Christensen 2006)`),
    med_ingest_high_t_hist = med_ingest_high_t_yr*(`Historical estimate`-`Southern hemisphere historic estimate (Christensen 2006)`),
    IQR25_ingest_high_t_hist = IQR25_ingest_high_t_yr*(`Historical estimate`-`Southern hemisphere historic estimate (Christensen 2006)`),
    IQR75_ingest_high_t_hist = IQR75_ingest_high_t_yr*(`Historical estimate`-`Southern hemisphere historic estimate (Christensen 2006)`)
    
  ) %>% 
  
  mutate_at(vars(c("med_ingest_low_t_curr":"IQR75_ingest_high_t_curr",
                   "med_ingest_low_t_hist":"IQR75_ingest_high_t_hist")), ~case_when(SpeciesCode == "bp" ~ .*0.8,     # correcting for proportion of diet that is fish; and proportion of individuals not in Southern Hemisphere
                                                                                    SpeciesCode == "bw" ~ .,
                                                                                    SpeciesCode == "mn" | region == "Temperate" ~ .*0.55, # FILTER not working like I need
                                                                                    SpeciesCode == "bb" ~ .))




# And now for the prey
# YEARLY INDIVIDUAL PREY, ANTARCTIC
Annual_ingestion_ind_Antarctic_bb <- Annual_filtfeed %>% 
  filter(region == "Polar", abbr_binom(Species) == "B. bonaerensis") %>% 
  ggplot() +
  
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_yr, ymax = IQR75_ingest_low_t_yr,
                  x=days_feeding),
              fill = "grey80", alpha = 0.5) +
  geom_line(aes(days_feeding, med_ingest_low_t_yr, color = abbr_binom(Species)),
            size = 1.5) +
  
  scale_colour_manual(values = pal2) +
  labs(x = "Days feeding",
       y = bquote(atop('Projected krill consumed',
                       ~(tonnes~ind^-1~yr^-1)))) + 

  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_continuous(limits = c(0,550),
                     breaks = c(0, 100, 200, 300, 400, 500)) +
  
  # scale_y_log10(labels = scales::comma, limits = c(10, 1000),
  #               breaks = c(30,100,300,600)) +
  #annotation_logticks(sides = "l") +
  
  theme_classic(base_size = 35) +
  
  theme(strip.text = element_text(face = "italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")
Annual_ingestion_ind_Antarctic_bb

dev.copy2pdf(file="Annual_ingestion_ind_Antarctic_bb.pdf", width=14, height=8)



Annual_ingestion_ind_Antarctic_mn <- Annual_filtfeed %>% 
  filter(region == "Polar", abbr_binom(Species) == "M. novaeangliae") %>% 
  ggplot() +
  
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_yr, ymax = IQR75_ingest_low_t_yr,
                  x=days_feeding),
              fill = "grey80", alpha = 0.5) +
  geom_line(aes(days_feeding, med_ingest_low_t_yr, color = abbr_binom(Species)),
            size = 1.5) +
  
  scale_colour_manual(values = pal2) +
  labs(x = "Days feeding") + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_continuous(limits = c(0,550),
                     breaks = c(0, 100, 200, 300, 400, 500)) +
  # scale_y_log10(labels = scales::comma, limits = c(10, 1000),
  #               breaks = c(30,100,300,600)) +
  #annotation_logticks(sides = "l") +
  
  theme_classic(base_size = 35) +
  
  theme(strip.text = element_text(face = "italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")
Annual_ingestion_ind_Antarctic_mn


Fig_3_Antarctic <- ggarrange(Annual_ingestion_ind_Antarctic_bb, prey_yr_prev_summ_bb_plot,
                             Annual_ingestion_ind_Antarctic_mn, prey_yr_prev_summ_mn_plot, 
                                labels = c("D","","E",""), 
                                font.label = list(size = 24),
                                widths = c(3.25, 1, 3.25, 1),
                                legend = "none",
                                ncol = 4, nrow = 1)
Fig_3_Antarctic
  
dev.copy2pdf(file="Fig_3_Antarctic.pdf", width=20, height=6)





# YEARLY INDIVIDUAL PREY, NON-ANTARCTIC
Annual_ingestion_ind_nonAntarctic_mn <- Annual_filtfeed %>% 
  filter(region == "Temperate", abbr_binom(Species) == "M. novaeangliae") %>% 
  ggplot() +
  
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_yr, ymax = IQR75_ingest_low_t_yr,
                  x=days_feeding),
              fill = "grey80", alpha = 0.5) +
  geom_line(aes(days_feeding, med_ingest_low_t_yr, color = abbr_binom(Species)),
            size = 1.5) +
  
  scale_colour_manual(values = pal2) +
  # scale_y_log10(labels = scales::comma, limits = c(20, 3500),
  #               breaks = c(30,100,300,600,1500,3000)) +
  # annotation_logticks(sides = "l") +
  
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_continuous(limits = c(20, 2500),
                     breaks = c(100,500,1000,1500,
                                2000)) +
  theme_classic(base_size = 35) +

  theme(strip.text = element_text(face = "italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")
Annual_ingestion_ind_nonAntarctic_mn

dev.copy2pdf(file="Annual_ingestion_ind_nonAntarctic_mn.pdf", width=14, height=7)
   



Annual_ingestion_ind_nonAntarctic_bp <- Annual_filtfeed %>% 
  filter(region == "Temperate", abbr_binom(Species) == "B. physalus") %>% 
  ggplot() +
  
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_yr, ymax = IQR75_ingest_low_t_yr,
                  x=days_feeding),
              fill = "grey80", alpha = 0.5) +
  geom_line(aes(days_feeding, med_ingest_low_t_yr, color = abbr_binom(Species)),
            size = 1.5) +
  
  scale_colour_manual(values = pal2) +
  # scale_y_log10(labels = scales::comma, limits = c(20, 3500),
  #               breaks = c(30,100,300,600,1500,3000)) +
  # annotation_logticks(sides = "l") +
  
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_continuous(limits = c(20, 2500),
                     breaks = c(100,500,1000,1500,
                                2000)) +
  theme_classic(base_size = 35) +

  
  theme(strip.text = element_text(face = "italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")
Annual_ingestion_ind_nonAntarctic_bp

dev.copy2pdf(file="Annual_ingestion_ind_nonAntarctic_bp.pdf", width=14, height=7)



Annual_ingestion_ind_nonAntarctic_bw <- Annual_filtfeed %>% 
  filter(region == "Temperate", abbr_binom(Species) == "B. musculus") %>% 
  ggplot() +
  
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_yr, ymax = IQR75_ingest_low_t_yr,
                  x=days_feeding),
              fill = "grey80", alpha = 0.5) +
  geom_line(aes(days_feeding, med_ingest_low_t_yr, color = abbr_binom(Species)),
            size = 1.5) +
  
  scale_colour_manual(values = pal2) +
  
  # scale_y_log10(labels = scales::comma, limits = c(20, 3500),
  #               breaks = c(30,100,300,600,1500,3000)) +
  # annotation_logticks(sides = "l") +
  
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_continuous(limits = c(20, 2500),
                     breaks = c(100,500,1000,1500,
                                2000)) +
  theme_classic(base_size = 35) +
  
  theme(strip.text = element_text(face = "italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")
Annual_ingestion_ind_nonAntarctic_bw

dev.copy2pdf(file="Annual_ingestion_ind_nonAntarctic_bw.pdf", width=14, height=7)



Fig_3_nonAntarctic <- ggarrange(Annual_ingestion_ind_nonAntarctic_mn, prey_yr_prev_summ_mn_plot,
                                Annual_ingestion_ind_nonAntarctic_bp, prey_yr_prev_summ_bp_plot,
                                Annual_ingestion_ind_nonAntarctic_bw, prey_yr_prev_summ_bw_plot,
                             labels = c("A","","B","","C",""), 
                             font.label = list(size = 28),
                             widths = c(3, 1, 3, 1, 3, 1),
                             legend = "none",
                             ncol = 6, nrow = 1)
Fig_3_nonAntarctic

dev.copy2pdf(file="Fig_3_nonAntarctic.pdf", width=28, height=6)




# Figure 4----


# Southern Hemisphere calculations (Figure 4)----

krill_Ant_projection <- read_csv("krill_biomass_estimates_for_editing_10.21.20.csv")

#load("daynights_AntProj.RData")


latitudes <- tribble(
  ~ region,    ~ latitude,
  "Polar",     65,
  "Temperate", 36
)

# Overall daily rate, Modified by MSS
estimate_daily <- function(rate_estimates, latitude, yday_center, season_len) {
  yday_start <- floor(yday_center - season_len / 2)
  yday_end <- floor(yday_center + season_len / 2)
  year_start <- as.POSIXct("2019-12-31", tz = "UTC")
  feeding_days <- year_start + lubridate::days(yday_start:yday_end)
  
  daylengths <- crossing(day = feeding_days, latitude) %>% 
    mutate(
      sunrise = maptools::sunriset(cbind(0, latitude), day, direction = "sunrise"),
      sunset = maptools::sunriset(cbind(0, latitude), day, direction = "sunset"),
      daylen = 24 * (sunset - sunrise),
      nightlen = 24 - daylen,
      yday = lubridate::yday(day)
    )
  
  
  
  #browser()
  # MATT: change rate_estimates to biomass_estimates
  daily_rates <- krill_Ant_projection %>% 
    left_join(daylengths, by = "region") %>% 
    # MATT: change daily_rate to daily_biomass
    group_by(SpeciesCode, region) %>% 
    # mutate(biomass_best_low_m3 = rlnorm(n(), Biomass[1], `Biomass sd`[1]),
    #        biomass_best_high_m3 = rlnorm(n(), BiomassTop50[1], BiomassTop50sd[1]),
    mutate(daily_rate = day_rate * daylen + night_rate * nightlen) %>% 
    ungroup() %>% 
    mutate(
      daily_mean_low_kg = pmap_dbl(list(Biomass, `Biomass sd`, daily_rate), ~mean(rlnorm(floor(..3), ..1, ..2))),
      daily_mean_high_kg = pmap_dbl(list(BiomassTop50, BiomassTop50sd, daily_rate), ~mean(rlnorm(floor(..3), ..1, ..2))),
      daily_consumption_low_kg = daily_mean_low_kg*Engulf_cap_m3*daily_rate, 
      daily_consumption_high_kg = daily_mean_high_kg*Engulf_cap_m3*daily_rate) 
}

# MATT: Change krill_rate_estimates to krill_biomass_estimates. Should have
# columns day_biomass and night_biomass.

aug1 <- 213
krill_daily_Ant_projection <- estimate_daily(krill_Ant_projection, latitudes, aug1, 120)  %>% 
  mutate(Total_energy_intake_best_low_kJ = case_when(region == "Polar" ~ daily_consumption_low_kg*4575,
                                                     region == "Temperate" ~ daily_consumption_low_kg*3628),
         Total_energy_intake_best_high_kJ = case_when(region == "Polar" ~ daily_consumption_high_kg*4575,
                                                      region == "Temperate" ~ daily_consumption_high_kg*3628),
         Mass_specifc_energy_intake_best_high_kJ = Total_energy_intake_best_high_kJ/Mass_est_t*1000,
         Mass_specifc_energy_intake_best_low_kJ = Total_energy_intake_best_low_kJ/Mass_est_t*1000
  )


krill_daily_Ant_projection <- krill_daily_Ant_projection %>%  
  mutate(Species = case_when(
    SpeciesCode == "bw" ~ "Balaenoptera musculus",
    SpeciesCode == "bp" ~ "Balaenoptera physalus",
    SpeciesCode == "mn" ~ "Megaptera novaeangliae",
    SpeciesCode == "bb" ~ "Balaenoptera bonaerensis", 
    SpeciesCode == "be" ~ "Balaenoptera brydei",
    SpeciesCode == "bs" ~ "Balaenoptera borealis")) 

krill_daily_Ant_projection %>% 
  group_by(Species) %>% 
  filter(SpeciesCode == "bw") %>%
  pull(daily_consumption_high_kg) %>%
  summary()

#save(krill_daily_Ant_projection, file = "daynights_AntProj.RData")



krill_yearly_pop_AntProj <- krill_daily_Ant_projection %>% 
  left_join(pop_data, by = "SpeciesCode") %>% 
  mutate(
    prey_low60curr = daily_consumption_low_kg*`Southern hemisphere population estimate (Christensen 2006)`*60,
    prey_low90curr = daily_consumption_low_kg*`Southern hemisphere population estimate (Christensen 2006)`*90,
    prey_low120curr = daily_consumption_low_kg*`Southern hemisphere population estimate (Christensen 2006)`*120,
    prey_low182.5curr = daily_consumption_low_kg*`Southern hemisphere population estimate (Christensen 2006)`*182.5,
    prey_low60hist = daily_consumption_low_kg*`Southern hemisphere historic estimate (Christensen 2006)`*60,
    prey_low90hist = daily_consumption_low_kg*`Southern hemisphere historic estimate (Christensen 2006)`*90,
    prey_low120hist = daily_consumption_low_kg*`Southern hemisphere historic estimate (Christensen 2006)`*120,
    prey_low182.5hist = daily_consumption_low_kg*`Southern hemisphere historic estimate (Christensen 2006)`*182.5,
    prey_high60curr = daily_consumption_high_kg*`Southern hemisphere population estimate (Christensen 2006)`*60,
    prey_high90curr = daily_consumption_high_kg*`Southern hemisphere population estimate (Christensen 2006)`*90,
    prey_high120curr = daily_consumption_high_kg*`Southern hemisphere population estimate (Christensen 2006)`*120,
    prey_high182.5curr = daily_consumption_high_kg*`Southern hemisphere population estimate (Christensen 2006)`*182.5,
    prey_high60hist = daily_consumption_high_kg*`Southern hemisphere historic estimate (Christensen 2006)`*60,
    prey_high90hist = daily_consumption_high_kg*`Southern hemisphere historic estimate (Christensen 2006)`*90,
    prey_high120hist = daily_consumption_high_kg*`Southern hemisphere historic estimate (Christensen 2006)`*120,
    prey_high182.5hist = daily_consumption_high_kg*`Southern hemisphere historic estimate (Christensen 2006)`*182.5
  )


# Annual Southern Hemisphere population krill prey, summary table
# Southern Ocean current populations prey consumption (low and high)
summ_prey_annual_curr_pop_stats <- krill_yearly_pop_AntProj %>% 
  group_by(SpeciesCode, region) %>% 
  summarise(
    yr_prey_low_60days_IQR25 = round(quantile(prey_low60curr, probs = 0.25, na.rm = TRUE)/1000000000, 2),
    yr_prey_low_60days_IQR75 = round(quantile(prey_low60curr, probs = 0.75, na.rm = TRUE)/1000000000, 2),
    yr_prey_low_90days_IQR25 = round(quantile(prey_low90curr, probs = 0.25, na.rm = TRUE)/1000000000, 2),
    yr_prey_low_90days_IQR75 = round(quantile(prey_low90curr, probs = 0.75, na.rm = TRUE)/1000000000, 2),
    yr_prey_low_120days_IQR25 = round(quantile(prey_low120curr, probs = 0.25, na.rm = TRUE)/1000000000, 2),
    yr_prey_low_120days_IQR75 = round(quantile(prey_low120curr, probs = 0.75, na.rm = TRUE)/1000000000, 2),
    yr_prey_low_182.5days_IQR25 = round(quantile(prey_low182.5curr, probs = 0.25, na.rm = TRUE)/1000000000, 2),
    yr_prey_low_182.5days_IQR75 = round(quantile(prey_low182.5curr, probs = 0.75, na.rm = TRUE)/1000000000, 2),
    
    yr_prey_high_60days_IQR25 = round(quantile(prey_high60curr, probs = 0.25, na.rm = TRUE)/1000000000, 2),
    yr_prey_high_60days_IQR75 = round(quantile(prey_high60curr, probs = 0.75, na.rm = TRUE)/1000000000, 2),
    yr_prey_high_90days_IQR25 = round(quantile(prey_high90curr, probs = 0.25, na.rm = TRUE)/1000000000, 2),
    yr_prey_high_90days_IQR75 = round(quantile(prey_high90curr, probs = 0.75, na.rm = TRUE)/1000000000, 2),
    yr_prey_high_120days_IQR25 = round(quantile(prey_high120curr, probs = 0.25, na.rm = TRUE)/1000000000, 2),
    yr_prey_high_120days_IQR75 = round(quantile(prey_high120curr, probs = 0.75, na.rm = TRUE)/1000000000, 2),
    yr_prey_high_182.5days_IQR25 = round(quantile(prey_high182.5curr, probs = 0.25, na.rm = TRUE)/1000000000, 2),
    yr_prey_high_182.5days_IQR75 = round(quantile(prey_high182.5curr, probs = 0.75, na.rm = TRUE)/1000000000, 2)
  ) %>% 
  unite("Krill consumption (Mt ind yr), 'random', 60 days feeding, IQR", 
        c(yr_prey_low_60days_IQR25, yr_prey_low_60days_IQR75), sep = "-") %>% 
  unite("Krill consumption (Mt ind yr), 'Top50', 60 days feeding, IQR", 
        c(yr_prey_high_60days_IQR25, yr_prey_high_60days_IQR75), sep = "-") %>% 
  unite("Krill consumption (Mt ind yr), 'random', 90 days feeding, IQR", 
        c(yr_prey_low_90days_IQR25, yr_prey_low_90days_IQR75), sep = "-") %>% 
  unite("Krill consumption (Mt ind yr), 'Top50', 90 days feeding, IQR", 
        c(yr_prey_high_90days_IQR25, yr_prey_high_90days_IQR75), sep = "-") %>% 
  unite("Krill consumption (Mt ind yr), 'random', 120 days feeding, IQR", 
        c(yr_prey_low_120days_IQR25, yr_prey_low_120days_IQR75), sep = "-") %>% 
  unite("Krill consumption (Mt ind yr), 'Top50', 120 days feeding, IQR", 
        c(yr_prey_high_120days_IQR25, yr_prey_high_120days_IQR75), sep = "-") %>% 
  unite("Krill consumption (Mt ind yr), 'random', 182.5 days feeding, IQR", 
        c(yr_prey_low_182.5days_IQR25, yr_prey_low_182.5days_IQR75), sep = "-") %>% 
  unite("Krill consumption (Mt ind yr), 'Top50', 182.5 days feeding, IQR", 
        c(yr_prey_high_182.5days_IQR25, yr_prey_high_182.5days_IQR75), sep = "-")



# Southern Ocean historic populations prey consumption (low and high)
summ_prey_annual_hist_pop_stats <- krill_yearly_pop_AntProj %>% 
  group_by(SpeciesCode, region) %>% 
  summarise(
    yr_prey_low_60days_IQR25 = round(quantile(prey_low60hist, probs = 0.25, na.rm = TRUE)/1000000000, 2),
    yr_prey_low_60days_IQR75 = round(quantile(prey_low60hist, probs = 0.75, na.rm = TRUE)/1000000000, 2),
    yr_prey_low_90days_IQR25 = round(quantile(prey_low90hist, probs = 0.25, na.rm = TRUE)/1000000000, 2),
    yr_prey_low_90days_IQR75 = round(quantile(prey_low90hist, probs = 0.75, na.rm = TRUE)/1000000000, 2),
    yr_prey_low_120days_IQR25 = round(quantile(prey_low120hist, probs = 0.25, na.rm = TRUE)/1000000000, 2),
    yr_prey_low_120days_IQR75 = round(quantile(prey_low120hist, probs = 0.75, na.rm = TRUE)/1000000000, 2),
    yr_prey_low_182.5days_IQR25 = round(quantile(prey_low182.5hist, probs = 0.25, na.rm = TRUE)/1000000000, 2),
    yr_prey_low_182.5days_IQR75 = round(quantile(prey_low182.5hist, probs = 0.75, na.rm = TRUE)/1000000000, 2),
    
    yr_prey_high_60days_IQR25 = round(quantile(prey_high60hist, probs = 0.25, na.rm = TRUE)/1000000000, 2),
    yr_prey_high_60days_IQR75 = round(quantile(prey_high60hist, probs = 0.75, na.rm = TRUE)/1000000000, 2),
    yr_prey_high_90days_IQR25 = round(quantile(prey_high90hist, probs = 0.25, na.rm = TRUE)/1000000000, 2),
    yr_prey_high_90days_IQR75 = round(quantile(prey_high90hist, probs = 0.75, na.rm = TRUE)/1000000000, 2),
    yr_prey_high_120days_IQR25 = round(quantile(prey_high120hist, probs = 0.25, na.rm = TRUE)/1000000000, 2),
    yr_prey_high_120days_IQR75 = round(quantile(prey_high120hist, probs = 0.75, na.rm = TRUE)/1000000000, 2),
    yr_prey_high_182.5days_IQR25 = round(quantile(prey_high182.5hist, probs = 0.25, na.rm = TRUE)/1000000000, 2),
    yr_prey_high_182.5days_IQR75 = round(quantile(prey_high182.5hist, probs = 0.75, na.rm = TRUE)/1000000000, 2)
  ) %>% 
  unite("Krill consumption (Mt ind yr), 'random', 60 days feeding, IQR", 
        c(yr_prey_low_60days_IQR25, yr_prey_low_60days_IQR75), sep = "-") %>% 
  unite("Krill consumption (Mt ind yr), 'Top50', 60 days feeding, IQR", 
        c(yr_prey_high_60days_IQR25, yr_prey_high_60days_IQR75), sep = "-") %>% 
  unite("Krill consumption (Mt ind yr), 'random', 90 days feeding, IQR", 
        c(yr_prey_low_90days_IQR25, yr_prey_low_90days_IQR75), sep = "-") %>% 
  unite("Krill consumption (Mt ind yr), 'Top50', 90 days feeding, IQR", 
        c(yr_prey_high_90days_IQR25, yr_prey_high_90days_IQR75), sep = "-") %>% 
  unite("Krill consumption (Mt ind yr), 'random', 120 days feeding, IQR", 
        c(yr_prey_low_120days_IQR25, yr_prey_low_120days_IQR75), sep = "-") %>% 
  unite("Krill consumption (Mt ind yr), 'Top50', 120 days feeding, IQR", 
        c(yr_prey_high_120days_IQR25, yr_prey_high_120days_IQR75), sep = "-") %>% 
  unite("Krill consumption (Mt ind yr), 'random', 182.5 days feeding, IQR", 
        c(yr_prey_low_182.5days_IQR25, yr_prey_low_182.5days_IQR75), sep = "-") %>% 
  unite("Krill consumption (Mt ind yr), 'Top50', 182.5 days feeding, IQR", 
        c(yr_prey_high_182.5days_IQR25, yr_prey_high_182.5days_IQR75), sep = "-")





days_feeding = rep(60:183, each = 5)
Annual_filtfeed_Ant_projection <- krill_daily_Ant_projection %>% 
  group_by(Species, SpeciesCode, region) %>% 
  summarise(
    #Individual
    med_filt = median(Engulf_cap_m3*daily_rate),
    IQR25_filt = quantile(Engulf_cap_m3*daily_rate, probs = 0.25, na.rm = TRUE),
    IQR75_filt = quantile(Engulf_cap_m3*daily_rate, probs = 0.75, na.rm = TRUE),
    med_ingest_low_t = median(daily_consumption_low_kg, na.rm = TRUE)/1000,
    IQR25_ingest_low_t = quantile(daily_consumption_low_kg, probs = 0.25, na.rm = TRUE)/1000,
    IQR75_ingest_low_t = quantile(daily_consumption_low_kg, probs = 0.75, na.rm = TRUE)/1000,
    med_ingest_high_t = median(daily_consumption_high_kg, na.rm = TRUE)/1000,
    IQR25_ingest_high_t = quantile(daily_consumption_high_kg, probs = 0.25, na.rm = TRUE)/1000,
    IQR75_ingest_high_t = quantile(daily_consumption_high_kg, probs = 0.75, na.rm = TRUE)/1000
  ) %>%
  crossing(days_feeding) %>% 
  mutate_at(c("med_filt", "IQR25_filt", "IQR75_filt"), funs(yr = .*days_feeding)) %>% 
  mutate_at(c("med_ingest_low_t", "IQR25_ingest_low_t", "IQR75_ingest_low_t"), funs(yr = .*days_feeding)) %>% 
  mutate_at(c("med_ingest_high_t", "IQR25_ingest_high_t", "IQR75_ingest_high_t"), funs(yr = .*days_feeding)) %>% 
  left_join(pop_data, by = c("Species", "SpeciesCode")) %>% 
  mutate(
    
    #Current population
    med_filt_curr = med_filt_yr*`Southern hemisphere population estimate (Christensen 2006)`,
    IQR25_filt_curr = IQR25_filt_yr*`Southern hemisphere population estimate (Christensen 2006)`,
    IQR75_filt_curr = IQR75_filt_yr*`Southern hemisphere population estimate (Christensen 2006)`,
    med_ingest_low_t_curr = med_ingest_low_t_yr*`Southern hemisphere population estimate (Christensen 2006)`,
    IQR25_ingest_low_t_curr = IQR25_ingest_low_t_yr*`Southern hemisphere population estimate (Christensen 2006)`,
    IQR75_ingest_low_t_curr = IQR75_ingest_low_t_yr*`Southern hemisphere population estimate (Christensen 2006)`,
    med_ingest_high_t_curr = med_ingest_high_t_yr*`Southern hemisphere population estimate (Christensen 2006)`,
    IQR25_ingest_high_t_curr = IQR25_ingest_high_t_yr*`Southern hemisphere population estimate (Christensen 2006)`,
    IQR75_ingest_high_t_curr = IQR75_ingest_high_t_yr*`Southern hemisphere population estimate (Christensen 2006)`,
    
    #Historic population
    med_filt_hist = med_filt_yr*`Southern hemisphere historic estimate (Christensen 2006)`,
    IQR25_filt_hist = IQR25_filt_yr*`Southern hemisphere historic estimate (Christensen 2006)`,
    IQR75_filt_hist = IQR75_filt_yr*`Southern hemisphere historic estimate (Christensen 2006)`,
    med_ingest_low_t_hist = med_ingest_low_t_yr*`Southern hemisphere historic estimate (Christensen 2006)`,
    IQR25_ingest_low_t_hist = IQR25_ingest_low_t_yr*`Southern hemisphere historic estimate (Christensen 2006)`,
    IQR75_ingest_low_t_hist = IQR75_ingest_low_t_yr*`Southern hemisphere historic estimate (Christensen 2006)`,
    med_ingest_high_t_hist = med_ingest_high_t_yr*`Southern hemisphere historic estimate (Christensen 2006)`,
    IQR25_ingest_high_t_hist = IQR25_ingest_high_t_yr*`Southern hemisphere historic estimate (Christensen 2006)`,
    IQR75_ingest_high_t_hist = IQR75_ingest_high_t_yr*`Southern hemisphere historic estimate (Christensen 2006)`
    
  )



# YEARLY POPULATION FILTRATION, ANTARCTIC
Annual_filtration_Pop_Antarctic <- Annual_filtfeed_Ant_projection %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. bonaerensis", "M. novaeangliae", "B. physalus")) %>% 
  ggplot() +
  geom_ribbon(aes(ymin = IQR25_filt_curr/1e9, ymax = IQR75_filt_curr/1e9, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_ribbon(aes(ymin = IQR25_filt_hist/1e9, ymax = IQR75_filt_hist/1e9, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_line(aes(days_feeding, med_filt_curr/1e9, color = Species), 
            size = 1.5) +
  geom_line(aes(days_feeding, med_filt_hist/1e9, color = Species), 
            size = 1.5) +
  scale_colour_manual(values = pal) +
  facet_grid(~Species) +   # Can add region in here if desired
  labs(x = "Days feeding",
       y = bquote(atop('Projected water filtered',
                       ~(km^3~yr^-1)))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_continuous(limits = c(0.5, 1800),
                breaks = c(1,500,1000,1500)) +
  # scale_y_log10(labels = scales::comma, limits = c(0.5, 2500),
  #               breaks = c(1,5,10,25,100,250,1000,2500)) +
  theme_minimal(base_size = 32) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none",
        axis.text.x = element_text(size = 20)
        )
Annual_filtration_Pop_Antarctic

dev.copy2pdf(file="Annual_filtration_Pop_Antarctic.pdf", width=14, height=8)


# Southern Hemisphere prey, COMBINED
Annual_ingestion_PopComb_Antarctic <- Annual_filtfeed_Ant_projection %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. bonaerensis", "M. novaeangliae", "B. physalus")) %>% 
  ggplot() +
  
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_hist/1e6, ymax = IQR75_ingest_low_t_hist/1e6, 
                  x=days_feeding), 
              fill = "grey80", alpha = 0.5) +
  
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_curr/1e6, ymax = IQR75_ingest_low_t_curr/1e6, 
                  x=days_feeding), 
              fill = "grey80", alpha = 0.5) +
  
  geom_line(aes(days_feeding, med_ingest_low_t_curr/1e6, color = Species), 
            size = 1.5) +
  geom_line(aes(days_feeding, med_ingest_low_t_hist/1e6, color = Species), 
            size = 1.5) +

  scale_colour_manual(values = pal) +
  facet_grid(~Species) +   # Can add region in here if desired
  labs(x = "Days feeding",
       y = bquote(atop('Projected krill consumed',
                       ~(Mt~yr^-1)))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_continuous(limits = c(0.1, 275),
                breaks = c(1, 50, 100, 150, 200, 250)) +
  # scale_y_log10(labels = scales::comma, limits = c(0.1, 400),
  #               breaks = c(1, 5, 10, 50, 100, 200, 400)) +
  theme_minimal(base_size = 32) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none",
        axis.text.x = element_text(size = 20)
        )
Annual_ingestion_PopComb_Antarctic

dev.copy2pdf(file="Annual_ingestion_PopComb_Antarctic.pdf", width=14, height=8)




summ_SO_krill_surplus <- Annual_filtfeed_Ant_projection %>% 
  group_by(SpeciesCode) %>% 
  summarise(
    
    Med_low_curr = median(med_ingest_low_t_curr/1e6, na.rm = TRUE),
    IQR25_low_Mt_curr = round(quantile(IQR25_ingest_low_t_curr/1e6, probs = 0.25, na.rm = TRUE), 4),
    IQR75_low_Mt_curr = round(quantile(IQR75_ingest_low_t_curr/1e6, probs = 0.25, na.rm = TRUE), 4),
    
    Med_high_curr = median(med_ingest_high_t_curr/1e6, na.rm = TRUE),
    IQR25_high_Mt_curr = round(quantile(IQR25_ingest_high_t_curr/1e6, probs = 0.25, na.rm = TRUE), 4),
    IQR75_high_Mt_curr = round(quantile(IQR75_ingest_high_t_curr/1e6, probs = 0.25, na.rm = TRUE), 4),
    
    Med_low_hist = median(med_ingest_low_t_hist/1e6, na.rm = TRUE),
    IQR25_low_Mt_hist = round(quantile(IQR25_ingest_low_t_hist/1e6, probs = 0.25, na.rm = TRUE), 4),
    IQR75_low_Mt_hist = round(quantile(IQR75_ingest_low_t_hist/1e6, probs = 0.25, na.rm = TRUE), 4),
    
    Med_high_hist = median(med_ingest_high_t_hist/1e6, na.rm = TRUE),
    IQR25_high_Mt_hist = round(quantile(IQR25_ingest_high_t_hist/1e6, probs = 0.25, na.rm = TRUE), 4),
    IQR75_high_Mt_hist = round(quantile(IQR75_ingest_high_t_hist/1e6, probs = 0.25, na.rm = TRUE), 4)
    
  ) %>% 
  col_summ(sum)



#Maximum extent of krill surplus from these four species
sum(summ_SO_krill_surplus$IQR75_low_t_hist)-sum(summ_SO_krill_surplus$IQR25_low_t_curr)

#Minimum extent of krill surplus from these four species
sum(summ_SO_krill_surplus$IQR25_low_t_hist)-sum(summ_SO_krill_surplus$IQR75_low_t_curr)

#Most likely extent of krill surplus from these four species
sum(summ_SO_krill_surplus$Med_low_hist)-sum(summ_SO_krill_surplus$Med_low_curr)


#Maximum extent of krill surplus from these four species
sum(summ_SO_krill_surplus$IQR75_high_t_hist)-sum(summ_SO_krill_surplus$IQR25_high_t_curr)

#Minimum extent of krill surplus from these four species
sum(summ_SO_krill_surplus$IQR25_high_t_hist)-sum(summ_SO_krill_surplus$IQR75_high_t_curr)

#Most likely extent of krill surplus from these four species
sum(summ_SO_krill_surplus$Med_high_hist)-sum(summ_SO_krill_surplus$Med_high_curr)






 # PPR calculations ----
######################################################################################
# Estimating amount of primary production required to sustain global whale populations
######################################################################################

PPR_data <- Annual_filtfeed_Ant_projection %>% 
  group_by(Species) %>% 
  summarise(
    PPR_required_curr_low = median(med_ingest_low_t_curr*0.1111*((1/0.1)^(2.2-1))),
            PPR_required_hist = median(med_ingest_low_t_hist*0.1111*((1/0.1)^(2.2-1)))) %>% 
  ungroup %>% 
col_summ(sum)


# Annual Southern Ocean NPP 1825 Mt per yr (1997-2013). From Arrigo 2014



#prepare data for krill PPR in the CCE
PPR_data <- prey_master_varying_DperYr %>% 
  mutate(CCE_population = case_when(Species == "Balaenoptera musculus" ~ 1647,
                                    Species == "Balaenoptera physalus" ~ 9029,
                                    Species == "Megaptera novaeangliae" ~ 1918),
         PPR_scenario_CCE = TotalAnnualPreyConsumed_kg*0.1111*((1/0.1)^(2.2-1))*CCE_population,  # See Barlow et al. 2008 eq 6 for explanation
         PPR_scenario_global = TotalAnnualPreyConsumed_kg*0.1111*((1/0.1)^(2.2-1))*`Population estimate`,
         PPR_scenario_global_historical = TotalAnnualPreyConsumed_kg*0.1111*((1/0.1)^(2.2-1))*`Total removed`)    


#From Chavez and Messie 2009
npp_gm2 <- 479
CCEarea_km2 <- 225000
npp_kg <- npp_gm2/1000 * CCEarea_km2*1000^2

PPR_table <- PPR_data %>% 
  group_by(Species, scenario) %>% 
  summarise(PPR_mean_CCE = mean(PPR_scenario_CCE),
            PPR_se_CCE = SE(PPR_scenario_CCE),
            NPP_to_whales_CCE = PPR_mean_CCE/npp_kg,
            PPR_mean_global = mean(PPR_scenario_global),
            PPR_se_global = SE(PPR_scenario_global),
            NPP_to_whales_global = PPR_mean_global/(47500000000*1000),
            PPR_mean_global_historical = mean(PPR_scenario_global_historical),
            PPR_se_global_historical = SE(PPR_scenario_global_historical),
            NPP_to_whales_global_historical = PPR_mean_global_historical/(47500000000*1000)) # From Longhurst et al. 1995


#summary tables of population nutrients recycled per year (Figure 5) ----

#Southern Ocean iron recycling...

rnorm_trunc <- function(n, mean = 0, sd = 1, lower = 0, upper = Inf) {
  result <- rnorm(n, mean, sd)
  result[result < lower] <- lower
  result[result > upper] <- upper
  result
}

# THIS one needed for plots
Fe_total_dw <- function(x) {
  (x*0.8)*0.25*0.146 # Projection of the total weight of feces, multiplied by 0.25 to convert to dry weight, and the 0.146 is the amount of iron (in kg!!!!) per TONNE of whale feces (converted from 0.000146 kg Fe per kg of feces)
} # From Doughty et al. 2016, Roman and McCarthy 2010


# THIS one needed for numeric estimates
Fe_total_dw <- function(x) {
  (x*runif(length(x), 0.7, 0.9))*   # saying the amount of iron excreted ranges uniformly between 70-90% ingested, see Ratnarajah et al. 2016
    0.25*           # Projection of the total weight of feces, multiplied by 0.25 to convert to dry weight, and the 0.146 is the amount of iron (in kg!!!!) per TONNE of whale feces (converted from 0.000146 kg Fe per kg of feces)
    rnorm(length(x), 0.146, 0.135) # average and sd of Fe conc in whale feces from: Ratnarajah et al. 2014  
} 


AntKrill_FeResRey <- function(x) {
  (x *(1-0.7726) * # convert from wet weight (in tonnes) krill ingested to dry weightfrom Nicol S, Stolp M, Nordstrom O (1992) Change in the gross biochemistry and mineral content accompanying the moult cycle in the Antarctic krill Euphausia superba. Mar Biol 113(2).
    1000 * # convert from tonnes to kg 
    rnorm(length(x), 145, 133.7) / # proportion E. superba that is iron (whole krill) in mg kg-1
    1e9) *     #convert from mg to tonnes
    rnorm(length(x), 0.8, 0.05)      #saying that 80% (distribution) of ingested iron is excreted                     
  
}



# Fecal iron added to calculated carbon exported
Fe_t_to_C_export_t <- function(x) {
  (((((x*runif(length(x), 0.8, 0.9) * #assumes 85% of feces stays in photic zone long enough to be utilized by phytoplankton;  # can change runif to 0.85 for graphing
         1e6) *             #convert tonnes to grams
    0.017907) *           #convert grams iron to moles iron
    5e4) * 12.0107) /     #convert to moles carbon exported (Lavery et al. 2010), convert back to grams carbon
    1e6)              #convert to tonnes carbon exported
}




Annual_filtfeed_Ant_projection_Nutrients <- Annual_filtfeed_Ant_projection %>%
  
  #Using estimate of fecal production
  mutate_at(vars(c("med_ingest_low_t_curr":"IQR75_ingest_high_t_curr",
                   "med_ingest_low_t_hist":"IQR75_ingest_high_t_hist")),
            .funs = list(Fe = ~Fe_total_dw(.)/1000)) %>% # need to divide by 1000 to get the total iron recycled in tonnes

  # # Simply saying that 80% of ingested iron is excreted
  # mutate_at(vars(c("med_ingest_low_t_curr":"IQR75_ingest_high_t_curr",
  #                  "med_ingest_low_t_hist":"IQR75_ingest_high_t_hist")),
  #           .funs = list(Fe = ~AntKrill_FeResRey(.))) %>% # need to divide by 1000 to get the total iron recycled in tonnes
   
  
  mutate_at(vars(c("med_ingest_low_t_curr_Fe":"IQR75_ingest_high_t_curr_Fe",
                   "med_ingest_low_t_hist_Fe":"IQR75_ingest_high_t_hist_Fe")), 
            .funs = list(C_produced_Mt = ~(Fe_t_to_C_export_t(.)*3)/1e6)) %>% # Carbon production (~3x what is exported) stimulated by Fe defecation in Mt C
  
  mutate_at(vars(c("med_ingest_low_t_curr":"IQR75_ingest_high_t_curr",
                   "med_ingest_low_t_hist":"IQR75_ingest_high_t_hist")),
            .funs = list(C_respired_Mt = ~(.*(0.1*0.75))/1e6)) %>%   # Carbon respired by populations in Mt C
  
  mutate_at(vars(c("med_ingest_low_t_curr_Fe":"IQR75_ingest_high_t_curr_Fe",
                   "med_ingest_low_t_hist_Fe":"IQR75_ingest_high_t_hist_Fe")), 
            .funs = list(C_exported_Mt = ~Fe_t_to_C_export_t(.)/1e6))




# Using estimate of fecal production
summ_SO_pop_Fe <- Annual_filtfeed_Ant_projection_Nutrients %>% 
  group_by(SpeciesCode, region) %>% 
  summarise(
    
    Fe_produced_med_curr_low = median(med_ingest_low_t_curr_Fe, na.rm = TRUE),
    Fe_produced_IQR25_curr_low = median(IQR25_ingest_low_t_curr_Fe, na.rm = TRUE),
    Fe_produced_IQR75_curr_low = median(IQR75_ingest_low_t_curr_Fe, na.rm = TRUE),
    
    Fe_produced_med_curr_high = median(med_ingest_high_t_curr_Fe, na.rm = TRUE),
    Fe_produced_IQR25_curr_high = median(IQR25_ingest_high_t_curr_Fe, na.rm = TRUE),
    Fe_produced_IQR75_curr_high = median(IQR75_ingest_high_t_curr_Fe, na.rm = TRUE),
    
    Fe_produced_med_hist_low = median(med_ingest_low_t_hist_Fe, na.rm = TRUE),
    Fe_produced_IQR25_hist_low = median(IQR25_ingest_low_t_hist_Fe, na.rm = TRUE),
    Fe_produced_IQR75_hist_low = median(IQR75_ingest_low_t_hist_Fe, na.rm = TRUE),
    
    Fe_produced_med_hist_high = median(med_ingest_high_t_hist_Fe, na.rm = TRUE),
    Fe_produced_IQR25_hist_high = median(IQR25_ingest_high_t_hist_Fe, na.rm = TRUE),
    Fe_produced_IQR75_hist_high = median(IQR75_ingest_high_t_hist_Fe, na.rm = TRUE)
    
  ) %>% 
  col_summ(sum)



summ_SO_pop_C_export <- Annual_filtfeed_Ant_projection_Nutrients %>% 
  group_by(SpeciesCode, region) %>% 
  summarise(
    
    C_respired_med_curr_low = median(med_ingest_low_t_curr_C_respired_Mt, na.rm = TRUE),
    C_respired_IQR25_curr_low = median(IQR25_ingest_low_t_curr_C_respired_Mt, na.rm = TRUE),
    C_respired_IQR75_curr_low = median(IQR75_ingest_low_t_curr_C_respired_Mt, na.rm = TRUE),

    C_respired_med_curr_high = median(med_ingest_high_t_curr_C_respired_Mt, na.rm = TRUE),
    C_respired_IQR25_curr_high = median(IQR25_ingest_high_t_curr_C_respired_Mt, na.rm = TRUE),
    C_respired_IQR75_curr_high = median(IQR75_ingest_high_t_curr_C_respired_Mt, na.rm = TRUE),

    C_respired_med_hist_low = median(med_ingest_low_t_hist_C_respired_Mt, na.rm = TRUE),
    C_respired_IQR25_hist_low = median(IQR25_ingest_low_t_hist_C_respired_Mt, na.rm = TRUE),
    C_respired_IQR75_hist_low = median(IQR75_ingest_low_t_hist_C_respired_Mt, na.rm = TRUE),

    C_respired_med_hist_high = median(med_ingest_high_t_hist_C_respired_Mt, na.rm = TRUE),
    C_respired_IQR25_hist_high = median(IQR25_ingest_high_t_hist_C_respired_Mt, na.rm = TRUE),
    C_respired_IQR75_hist_high = median(IQR75_ingest_high_t_hist_C_respired_Mt, na.rm = TRUE),
    
    C_induced_med_curr_low = median(med_ingest_low_t_curr_Fe_C_produced_Mt, na.rm = TRUE),
    C_induced_IQR25_curr_low = median(IQR25_ingest_low_t_curr_Fe_C_produced_Mt, na.rm = TRUE),
    C_induced_IQR75_curr_low = median(IQR75_ingest_low_t_curr_Fe_C_produced_Mt, na.rm = TRUE),
    
    C_induced_med_curr_high = median(med_ingest_high_t_curr_Fe_C_produced_Mt, na.rm = TRUE),
    C_induced_IQR25_curr_high = median(IQR25_ingest_high_t_curr_Fe_C_produced_Mt, na.rm = TRUE),
    C_induced_IQR75_curr_high = median(IQR75_ingest_high_t_curr_Fe_C_produced_Mt, na.rm = TRUE),
    
    C_induced_med_hist_low = median(med_ingest_low_t_hist_Fe_C_produced_Mt, na.rm = TRUE),
    C_induced_IQR25_hist_low = median(IQR25_ingest_low_t_hist_Fe_C_produced_Mt, na.rm = TRUE),
    C_induced_IQR75_hist_low = median(IQR75_ingest_low_t_hist_Fe_C_produced_Mt, na.rm = TRUE),
    
    C_induced_med_hist_high = median(med_ingest_high_t_hist_Fe_C_produced_Mt, na.rm = TRUE),
    C_induced_IQR25_hist_high = median(IQR25_ingest_high_t_hist_Fe_C_produced_Mt, na.rm = TRUE),
    C_induced_IQR75_hist_high = median(IQR75_ingest_high_t_hist_Fe_C_produced_Mt, na.rm = TRUE),
    
    
    C_export_med_curr_low = median(med_ingest_low_t_curr_Fe_C_exported_Mt, na.rm = TRUE),
    C_export_IQR25_curr_low = median(IQR25_ingest_low_t_curr_Fe_C_exported_Mt, na.rm = TRUE),
    C_export_IQR75_curr_low = median(IQR75_ingest_low_t_curr_Fe_C_exported_Mt, na.rm = TRUE),
    
    C_export_med_curr_high = median(med_ingest_high_t_curr_Fe_C_exported_Mt, na.rm = TRUE),
    C_export_IQR25_curr_high = median(IQR25_ingest_high_t_curr_Fe_C_exported_Mt, na.rm = TRUE),
    C_export_IQR75_curr_high = median(IQR75_ingest_high_t_curr_Fe_C_exported_Mt, na.rm = TRUE),
    
    C_export_med_hist_low = median(med_ingest_low_t_hist_Fe_C_exported_Mt, na.rm = TRUE),
    C_export_IQR25_hist_low = median(IQR25_ingest_low_t_hist_Fe_C_exported_Mt, na.rm = TRUE),
    C_export_IQR75_hist_low = median(IQR75_ingest_low_t_hist_Fe_C_exported_Mt, na.rm = TRUE),
    
    C_export_med_hist_high = median(med_ingest_high_t_hist_Fe_C_exported_Mt, na.rm = TRUE),
    C_export_IQR25_hist_high = median(IQR25_ingest_high_t_hist_Fe_C_exported_Mt, na.rm = TRUE),
    C_export_IQR75_hist_high = median(IQR75_ingest_high_t_hist_Fe_C_exported_Mt, na.rm = TRUE),
    
  )  %>% 
  col_summ(sum)

# Total difference in C export in Mt C yr-1 
# In 1900 global CO2 emissions was ~1958 Mt; whale C export was 57-196Mt or 4-10% of global CO2 emissions at the time
# In 2000 global CO2 emissions was ~24,670 Mt; whale C export was 6-23Mt or 0.02-0.09% of global CO2 emissions at the time




# Total Fe recycled by E. superba population
133e12/0.486 * # number of krill in entire pop (first number is biomass in Mt); at 0.486g per krill (Atkinson et al. 2009)
  1.5917e-7 * # converting from 2.85 nmol iron recycled per day per krill in g
  182.5 /     # number of productive days in Austral summer
  1e12      # converting from grams to Mt



Fe_t_to_C_export_t(
  665e12/0.486 * # number of krill in entire pop (first number is biomass in Mt); at 0.486g per krill (Atkinson et al. 2009)
    1.5917e-7 * # converting from 2.85 nmol iron recycled per day per krill in g
    182.5 /     # number of productive days in Austral summer
    1e12      # converting from grams to Mt
)*
  0.2 # the porportion of iron in krill feces that persists in the photic zone for use by phytoplankton




# Southern Hemisphere Fe, COMBINED
Annual_ingestion_PopComb_Antarctic_Fe <- Annual_filtfeed_Ant_projection_Nutrients %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. bonaerensis", "M. novaeangliae", "B. physalus")) %>% 
  ggplot() +
  
  # geom_line(aes(days_feeding, med_ingest_high_t_curr_Fe, color = Species),
  #           size = 1.5, linetype = "dashed") +
  # geom_line(aes(days_feeding, med_ingest_high_t_hist_Fe, color = Species),
  #           size = 1.5, linetype = "dashed") +
  
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_hist_Fe, ymax = IQR75_ingest_low_t_hist_Fe, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  # geom_ribbon(aes(ymin = IQR25_ingest_high_t_hist_Fe, ymax = IQR75_ingest_high_t_hist_Fe,
  #                 x=days_feeding),
  #             fill = "grey70", alpha = 0.5) +
  
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_curr_Fe, ymax = IQR75_ingest_low_t_curr_Fe, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  # geom_ribbon(aes(ymin = IQR25_ingest_high_t_curr_Fe, ymax = IQR75_ingest_high_t_curr_Fe,
  #                 x=days_feeding),
  #             fill = "grey70", alpha = 0.5) +
  geom_line(aes(days_feeding, med_ingest_low_t_curr_Fe, color = Species), 
            size = 1.5) +
  
  geom_line(aes(days_feeding, med_ingest_low_t_hist_Fe, color = Species), 
            size = 1.5) +

  scale_colour_manual(values = pal) +
  facet_grid(~Species) +   # Can add region in here if desired
  labs(x = "Days feeding",
       y = bquote(atop('Projected Fe recycled',
                       ~(tonnes~yr^-1)))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_continuous(limits = c(1, 8000),
                breaks = c(1, 1000, 2500, 5000, 7500)) +
  # scale_y_log10(labels = scales::comma, limits = c(2, 12000),
  #               breaks = c(1, 10, 100, 500, 1000, 2500, 5000, 10000)) +
  theme_minimal(base_size = 32) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none",
        axis.text.x = element_text(size = 20)
        )
Annual_ingestion_PopComb_Antarctic_Fe

dev.copy2pdf(file="Annual_ingestion_PopComb_Antarctic_Fe.pdf", width=14, height=8)


Fig_4_full <- ggarrange(Annual_filtration_Pop_Antarctic, 
                        Annual_ingestion_PopComb_Antarctic,
                        Annual_ingestion_PopComb_Antarctic_Fe,
                            labels = c("A","B","C"), 
                            font.label = list(size = 20),
                            align = "hv",
                            legend = "none",
                            ncol = 1, nrow = 3)
Fig_4_full




# Figure 4 Extended----
Annual_filtration_Pop_nonAntarctic <- Annual_filtfeed %>% 
  filter(region == "Temperate") %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)),  "M. novaeangliae", "B. physalus")) %>% 
  ggplot() +
  geom_ribbon(aes(ymin = IQR25_filt_curr/1e9, ymax = IQR75_filt_curr/1e9, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_line(aes(days_feeding, med_filt_curr/1e9, color = Species), 
            size = 1.5) +
  geom_ribbon(aes(ymin = IQR25_filt_hist/1e9, ymax = IQR75_filt_hist/1e9, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_line(aes(days_feeding, med_filt_hist/1e9, color = Species), 
            size = 1.5) +
  scale_colour_manual(values = pal) +
  facet_grid(~Species) +   # Can add region in here if desired
  labs(x = "Days feeding",
       y = bquote(atop('Projected water filtered',
                       ~(km^3~yr^-1)))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_continuous(breaks = c(1, 100, 200, 300, 400), 
                     limits = c(1, 400)) +
  #               breaks = c(1,5,10,25,100,250)) +
  # scale_y_log10(labels = scales::comma, limits = c(1, 500),
  #               breaks = c(1,5,10,25,100,250)) +
  theme_minimal(base_size = 32) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none",
        axis.text.x = element_text(size = 20))
Annual_filtration_Pop_nonAntarctic

dev.copy2pdf(file="Annual_filtration_Pop_nonAntarctic.pdf",  width=14, height=8)


Total_filtered <- Annual_filtfeed %>%
  group_by(SpeciesCode, region) %>%
  filter(days_feeding == 100) %>% 
  summarise(sp_annual_avg_curr = med_filt_curr/1e9,
            sp_annual_avg_hist = med_filt_hist/1e9)



# Northern Hemisphere prey, COMBINED
Annual_ingestion_PopComb_nonAntarctic <- Annual_filtfeed %>% 
  filter(region == "Temperate") %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "M. novaeangliae", "B. physalus")) %>% 
  ggplot() +
  
  geom_line(aes(days_feeding, med_ingest_high_t_curr/1e6, color = Species), 
            size = 1.5, linetype = "dashed") +
  geom_line(aes(days_feeding, med_ingest_high_t_hist/1e6, color = Species), 
            size = 1.5, linetype = "dashed") +
  
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_hist/1e6, ymax = IQR75_ingest_low_t_hist/1e6, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_ribbon(aes(ymin = IQR25_ingest_high_t_hist/1e6, ymax = IQR75_ingest_high_t_hist/1e6, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_curr/1e6, ymax = IQR75_ingest_low_t_curr/1e6, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_ribbon(aes(ymin = IQR25_ingest_high_t_curr/1e6, ymax = IQR75_ingest_high_t_curr/1e6, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_line(aes(days_feeding, med_ingest_low_t_curr/1e6, color = Species), 
            size = 1.5) +
  
  
  
  geom_line(aes(days_feeding, med_ingest_low_t_hist/1e6, color = Species), 
            size = 1.5) +
  
  
  scale_colour_manual(values = pal) +
  facet_grid(~Species) +   # Can add region in here if desired
  labs(x = "Days feeding",
       y = bquote(atop('Projected krill consumed',
                       ~(Mt~yr^-1)))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_continuous(breaks = c(0, 50, 100, 150, 200), 
                     limits = c(0, 200)) +
  # scale_y_log10(labels = scales::comma, limits = c(0.5, 200),
  #               breaks = c(1, 5, 10, 50, 100, 200)) +
  theme_minimal(base_size = 32) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none",
        axis.text.x = element_text(size = 20))
Annual_ingestion_PopComb_nonAntarctic

dev.copy2pdf(file="Annual_ingestion_PopComb_nonAntarctic.pdf", width=14, height=8)


Fig_4_Extended <- ggarrange(Annual_filtration_Pop_nonAntarctic, 
                            Annual_ingestion_PopComb_nonAntarctic,
                            labels = c("A","B"), 
                            font.label = list(size = 24),
                            align = "hv",
                            legend = "none",
                            ncol = 2, nrow = 1)
Fig_4_Extended


dev.copy2pdf(file="Fig_4_Extended.pdf", width=22, height=8)








# Southern Hemisphere C export, COMBINED

Annual_ingestion_PopComb_Antarctic_C_export <- Annual_filtfeed_Ant_projection_Nutrients %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. bonaerensis", "M. novaeangliae", "B. physalus")) %>% 
  ggplot() +


  geom_ribbon(aes(ymin = IQR25_ingest_low_t_curr_Fe_C_exported_Mt, ymax = IQR75_ingest_low_t_curr_Fe_C_exported_Mt, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_hist_Fe_C_exported_Mt, ymax = IQR75_ingest_low_t_hist_Fe_C_exported_Mt, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +

  
  geom_line(aes(days_feeding, med_ingest_low_t_curr_Fe_C_exported_Mt, color = Species), 
            size = 1.5) +
  geom_line(aes(days_feeding, med_ingest_low_t_hist_Fe_C_exported_Mt, color = Species), 
            size = 1.5) +

  scale_colour_manual(values = pal) +
  facet_grid(~Species) +   # Can add region in here if desired
  labs(x = "Days feeding",
       y = bquote(atop('Projected C exported',
                       ~(Mt~yr^-1)))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_log10(labels = function(x) sprintf("%g", x),
                limits = c(0.025, 120),
                breaks = c(0.01,0.1, 1, 10, 100)) +
  theme_minimal(base_size = 32) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none",
        axis.text.x = element_text(size = 20)
        )
Annual_ingestion_PopComb_Antarctic_C_export

dev.copy2pdf(file="Annual_ingestion_PopComb_Antarctic_C_export.pdf", width=14, height=8)


