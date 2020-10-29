# Code to generate consumption estimates for baleen whales with prey estimates from the Scaling Paper----

library(tidyverse)
library(maptools)
library(readxl)
library(ggpubr)
library(purrr)
library(lubridate)

SE = function(x){sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x)))}



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



# summ_table <- All_rorqual_deployments %>% filter(Phase == "Total") %>% 
#   col_summ(sum)

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

# KRILL data from Scaling Paper (updated 10/26/20)----


krill_data_Scaling_paper <- read_csv("KrillData_Scaling_paper.csv") %>% 
  mutate(Biomass_hyp_low = ifelse(region == "Polar", NA, log(Biomass) + log(0.17))) %>% 
  rename(SpeciesCode = "Species") %>% 
  mutate_at(vars(Biomass, `Biomass sd`), log)



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
  ungroup() %>% 
  left_join(krill_data_Scaling_paper, by = c("SpeciesCode", "region"))



# quick look at species average rates; Assuming Sep 15, 12.5 hours of daylight
krill_rate_estimates %>% 
  filter(SpeciesCode == "bp") %>% 
  mutate(daily_rate = day_rate * 12.5 + night_rate * 11.5) %>% 
  pull(daily_rate) %>% 
  quantile(c(0.25, 0.5, 0.75, 0.95))

# Overall daily rate, Modified by MATT
estimate_daily <- function(rate_estimates, season_len) {
  # Peak of feeding season
  polar_peak <- as.POSIXct("2019-12-31", tz = "UTC") + days(48) # Feb 17
  temperate_peak <- as.POSIXct("2019-12-31", tz = "UTC") + days(231) # Aug 18
  # Latitudes
  polar_lat <- -65
  temperate_lat <- 36
  half_season <- floor(season_len / 2)
  # Day and night lengths
  feeding_days <- tibble(
    # Feeding days, polar and temperate
    polar_feeding = seq(polar_peak - days(half_season),
                        polar_peak + days(half_season),
                        by = "1 day"),
    temperate_feeding = seq(temperate_peak - days(half_season),
                            temperate_peak + days(half_season),
                            by = "1 day"),
    # Hours day/night, polar and temperate
    polar_sunrise = sunriset(cbind(0, polar_lat), 
                             polar_feeding, 
                             direction = "sunrise"),
    polar_sunset = sunriset(cbind(0, polar_lat), 
                            polar_feeding, 
                            direction = "sunset"),
    polar_daylen = 24 * (polar_sunset - polar_sunrise),
    polar_nightlen = 24 - polar_daylen,
    temperate_sunrise = sunriset(cbind(0, temperate_lat), 
                                 temperate_feeding, 
                                 direction = "sunrise"),
    temperate_sunset = sunriset(cbind(0, temperate_lat), 
                                temperate_feeding, 
                                direction = "sunset"),
    temperate_daylen = 24 * (temperate_sunset - temperate_sunrise),
    temperate_nightlen = 24 - temperate_daylen
  )
  
  # This function allows for NAs in logmean (i.e. Antarctic hypothetical low) and truncation
  # to the upper 50% (i.e. to account for selectivity)
  get_daily_mean <- function(logmean, logsd, lunges, trunc) {
    rand <- function(n, meanlog, sdlog) {
      if (n == 0) {
        0
      } else if (trunc) {
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
    result[is.na(logmean)] <- NA
    result
  }
  
  rate_estimates %>% 
    full_join(feeding_days, by = character()) %>%
    mutate(
      daylen = ifelse(region == "Polar", polar_daylen, temperate_daylen),
      nightlen = ifelse(region == "Polar", polar_nightlen, temperate_nightlen),
      daily_rate = day_rate * daylen + night_rate * nightlen,
      daily_mean_hyp_low_kg = get_daily_mean(Biomass_hyp_low, `Biomass sd`, daily_rate, trunc = FALSE),
      daily_mean_kg = get_daily_mean(Biomass, `Biomass sd`, daily_rate, trunc = FALSE),
      daily_mean_hyp_high_kg = get_daily_mean(Biomass, `Biomass sd`, daily_rate, trunc = TRUE),
      
      daily_consumption_hyp_low_kg = daily_mean_hyp_low_kg * Engulf_cap_m3 * daily_rate,
      daily_consumption_kg = daily_mean_kg * Engulf_cap_m3 * daily_rate, 
      daily_consumption_hyp_high_kg = daily_mean_hyp_high_kg * Engulf_cap_m3 * daily_rate)
}

# MATT: Change krill_rate_estimates to krill_biomass_estimates. Should have
# columns day_biomass and night_biomass.


krill_daily <- krill_biomass_estimates %>% 
  estimate_daily(120) %>% 
  mutate(
    Total_energy_intake_low_kJ = case_when(
      region == "Polar" ~ daily_consumption_hyp_low_kg * 4575,
      region == "Temperate" ~ daily_consumption_hyp_low_kg * 3628
    ),
    Total_energy_intake_kJ = case_when(
      region == "Polar" ~ daily_consumption_kg * 4575,
      region == "Temperate" ~ daily_consumption_kg * 3628
    ),
    Total_energy_intake_high_kJ = case_when(
      region == "Polar" ~ daily_consumption_hyp_high_kg * 4575,
      region == "Temperate" ~ daily_consumption_hyp_high_kg * 3628
    ),
    Mass_specifc_energy_intake_low_kJ = Total_energy_intake_low_kJ / (Mass_est_t * 1000),
    Mass_specifc_energy_intake_kJ = Total_energy_intake_kJ / (Mass_est_t * 1000),
    Mass_specifc_energy_intake_high_kJ = Total_energy_intake_high_kJ / (Mass_est_t*1000),
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
