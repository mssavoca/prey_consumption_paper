# ANALYSIS file for submission Savoca et al. ----


library(tidyverse)
library(maptools)
library(readxl)
library(ggpubr)

SE = function(x){sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x)))}

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
    map(tbl, ~ if(is.numeric(.x)) f(.x) else NA)
  )
}

# stuff from Max to generate data ----

load("rates_CATS_updated_5_14_20.RData") # CATS deployment summary

load("rates_MedTermTag.RData") # Medium-term deployment summary

load("rates[1].RData") # DTAG deployment summary

load("tag_guide4.RData")  #loads the tag guide 


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
  bind_rows(medterm_lunges_long)

summ_table <- All_rorqual_deployments %>% filter(Phase == "Total") %>% 
  col_summ(sum)


# droned whales---
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
whale_lengths <- read.csv("whale_masses.csv") %>% 
  rename(SpeciesCode = "species") %>% 
  full_join(select(KC_nontagged_lengths, name, length, SpeciesCode)) %>% 
  left_join(engulf_allo, by = "SpeciesCode") %>% 
  mutate(Engulfment_m3 = length ^ slope * 10 ^ intercept,
         #numbers from Lockyer 1976
         Lockyer_a = case_when(SpeciesCode == "bb" ~ 0.0496,   # From Table 1, the rest from Table 2
                               SpeciesCode == "be" ~ 0.0122,
                               SpeciesCode == "mn" ~ 0.0158,
                               SpeciesCode == "bs" ~ 0.0242,
                               SpeciesCode == "bp" ~ 0.01265,
                               SpeciesCode == "bw" ~ 0.0025
         ),
         Lockyer_b = case_when(SpeciesCode == "bb" ~ 2.31,    # From Table 1, the rest from Table 2
                               SpeciesCode == "be" ~ 2.74,
                               SpeciesCode == "mn" ~ 2.95,
                               SpeciesCode == "bs" ~ 2.43,
                               SpeciesCode == "bp" ~ 2.995,
                               SpeciesCode == "bw" ~ 3.25
         ),
         Lockyer_mass_t = (Lockyer_a*length ^Lockyer_b)*1.06,
         Species = case_when(
           SpeciesCode == "bw" ~ "Balaenoptera musculus",
           SpeciesCode == "bp" ~ "Balaenoptera physalus",
           SpeciesCode == "mn" ~ "Megaptera novaeangliae",
           SpeciesCode %in% c("bb","ba") ~ "Balaenoptera bonaerensis", 
           SpeciesCode == "be" ~ "Balaenoptera edeni",
           SpeciesCode == "bs" ~ "Balaenoptera borealis")) 



# Dave's KRILL data (from 11/21/19)----
#natural log transformed, which is what we were donig in the prior iteration

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


# read in, clean, combine data
MRY_krill_data <- read_csv("MontereyKrillData (1).csv") %>% 
  filter(!Species %in% c("bigBW", "bb")) %>% 
  select(Species:BiomassTop50sd) %>%    # remember: biomass is in kg/m3
  mutate(
    Biomass_hyp_low = exp(log(Biomass) + log(0.17)),
    Study_Area = "Monterey",
    Region = "Eastern North Pacific") %>% 
  rename(SpeciesCode = "Species") 

SoCal_krill_data <- read_csv("SoCalKrillData (1).csv") %>% 
  filter(!Species %in% c("bigBW", "bb")) %>% 
  select(Species:BiomassTop50sd) %>%    # remember: biomass is in kg/m3
  mutate(
    Biomass_hyp_low = exp(log(Biomass) + log(0.17)),
    Study_Area = "SoCal",
    Region = "Eastern North Pacific") %>% 
  rename(SpeciesCode = "Species") 

# combining Monterey and SoCal prey data
ENP_krill_data <- rbind(MRY_krill_data, SoCal_krill_data) %>% 
  #mutate(Study_Area2 = Study_Area) %>% 
  pivot_wider(names_from = Study_Area, values_from = c(`Num Days`:Biomass_hyp_low)) %>% 
  #rename(Study_Area = Study_Area2) %>% 
  mutate(
    `Num Days` = `Num Days_Monterey` + `Num Days_SoCal`, 
    Biomass_hyp_low = pooled_log_mean(Biomass_hyp_low_Monterey, Biomass_hyp_low_SoCal),
    Biomass = pooled_log_mean(Biomass_Monterey, Biomass_SoCal),
    `Biomass sd` = pooled_sd_mean(`Biomass sd_Monterey`, `Biomass sd_SoCal`, `Num Days_Monterey`, `Num Days_SoCal`),
    BiomassTop50 = pooled_log_mean(BiomassTop50_Monterey, BiomassTop50_SoCal),
    BiomassTop50sd = pooled_sd_mean(`BiomassTop50sd_Monterey`, `BiomassTop50sd_SoCal`, `Num Days_Monterey`, `Num Days_SoCal`),
    Study_Area = NA
    # 
    # `Num Days` = coalesce(`Num Days_Monterey`, `Num Days_SoCal`),
    # Biomass_hyp_low = coalesce(Biomass_hyp_low_Monterey, Biomass_hyp_low_SoCal),
    # Biomass = coalesce(Biomass_Monterey, Biomass_SoCal),
    # `Biomass sd` = coalesce(`Biomass sd_Monterey`, `Biomass sd_SoCal`),
    # BiomassTop50 = coalesce(BiomassTop50_Monterey, BiomassTop50_SoCal),
    # BiomassTop50sd = coalesce(`BiomassTop50sd_Monterey`, `BiomassTop50sd_SoCal`)
  ) %>% 
  select(-c(`Num Days_Monterey`:Biomass_hyp_low_SoCal)) %>% 
  mutate_at(vars(Biomass_hyp_low:BiomassTop50sd) , ~log(10^.x))


WAP_krill_data <- read_csv("AntarcticKrillData (1).csv") %>% 
  filter(!Species %in% c("bigBW")) %>% 
  select(Species:BiomassTop50sd) %>%    # remember: biomass is in kg/m3
  mutate(
    Biomass_hyp_low = NA,
    Study_Area = "Antarctic",
    Region = "Antarctic") %>% 
  rename(SpeciesCode = "Species")

SA_krill_data <- read_csv("SouthAfricaKrillData (1).csv") %>% 
  filter(!Species %in% c("bigBW", "bw", "bp", "bb")) %>% 
  select(Species:BiomassTop50sd) %>%    # remember: biomass is in kg/m3
  mutate(
    Biomass_hyp_low = exp(log(Biomass) + log(0.17)),
    Study_Area = "South Africa",
    Region = "South Africa") %>% 
  rename(SpeciesCode = "Species") 


# combining ENP and SA prey data
ENP_krill_data <- rbind(MRY_krill_data, SoCal_krill_data) %>% 
  #mutate(Study_Area2 = Study_Area) %>% 
  pivot_wider(names_from = Study_Area, values_from = c(`Num Days`:Biomass_hyp_low)) %>% 
  #rename(Study_Area = Study_Area2) %>% 
  mutate(
    `Num Days` = `Num Days_Monterey` + `Num Days_SoCal`, 
    Biomass_hyp_low = pooled_log_mean(Biomass_hyp_low_Monterey, Biomass_hyp_low_SoCal),
    Biomass = pooled_log_mean(Biomass_Monterey, Biomass_SoCal),
    `Biomass sd` = pooled_sd_mean(`Biomass sd_Monterey`, `Biomass sd_SoCal`, `Num Days_Monterey`, `Num Days_SoCal`),
    BiomassTop50 = pooled_log_mean(BiomassTop50_Monterey, BiomassTop50_SoCal),
    BiomassTop50sd = pooled_sd_mean(`BiomassTop50sd_Monterey`, `BiomassTop50sd_SoCal`, `Num Days_Monterey`, `Num Days_SoCal`),
    Study_Area = NA
    # 
    # `Num Days` = coalesce(`Num Days_Monterey`, `Num Days_SoCal`),
    # Biomass_hyp_low = coalesce(Biomass_hyp_low_Monterey, Biomass_hyp_low_SoCal),
    # Biomass = coalesce(Biomass_Monterey, Biomass_SoCal),
    # `Biomass sd` = coalesce(`Biomass sd_Monterey`, `Biomass sd_SoCal`),
    # BiomassTop50 = coalesce(BiomassTop50_Monterey, BiomassTop50_SoCal),
    # BiomassTop50sd = coalesce(`BiomassTop50sd_Monterey`, `BiomassTop50sd_SoCal`)
  ) %>% 
  select(-c(`Num Days_Monterey`:Biomass_hyp_low_SoCal)) 


All_krill_data <- rbind(MRY_krill_data, SoCal_krill_data, WAP_krill_data, SA_krill_data) 

All_krill_data_ENPcombined <- rbind(ENP_krill_data, WAP_krill_data, SA_krill_data)

# Whale population data---- 
#Current data from IUCN 2019 Redlist, whaling data compiled in Rocha et al. 2014
pop_data <- read_excel("Filtration project_Whale population and feeding information.xlsx", sheet = 1) %>%
  mutate(SpeciesCode = case_when(Species == "Balaenoptera musculus" ~ "bw",
                                 Species == "Balaenoptera physalus" ~ "bp",
                                 Species == "Megaptera novaeangliae" ~ "mn",
                                 Species == "Balaenoptera borealis" ~ "bs",
                                 Species == "Balaenoptera edeni" ~ "be",
                                 Species == "Balaenoptera acutorostrata" ~ "ba",
                                 Species == "Balaenoptera bonaerensis" ~ "bb",
                                 Species == "Eubalaena glacialis" ~ "eg",
                                 Species == "Eubalaena japonica" ~ "ej",
                                 Species == "Eubalaena australis" ~ "ea",
                                 Species == "Balaena mysticetus" ~ "bm"
  ))



# data projected for Antarctic blue and fin whales
Krill_data_Ant_projection <- read_excel("AntarcticKrillData (1).xlsx", sheet = 3) %>% 
  filter(Species %in% c("bp","bw")) %>% 
  rename(SpeciesCode = "Species") %>% 
  mutate(region = "Polar")

All_krill_data_ENPcombined_noSA <- rbind(ENP_krill_data, WAP_krill_data) %>% 
  mutate(region = ifelse(Region == "Antarctic",  "Polar","Temperate")) %>% 
  filter(!SpeciesCode %in% c("bp","bw") | region != "Polar") %>% 
  bind_rows(Krill_data_Ant_projection) %>%        # bringing in projected data for Antarctic blue and fin whales
  select(-c(Region, Sv:SvTop50sd)) %>%
  mutate(Study_Area = case_when(region == "Temperate" ~ "Eastern North Pacific",
                                region == "Polar" ~ "Antarctic")) %>% 
  mutate_at(vars(Biomass_hyp_low:BiomassTop50sd) , ~log(.x))  # natural log transformed, which is what we were donig in the prior iteration


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
  left_join(select(tag_guide4, ID, study_area), by = "ID") %>% 
  mutate(
    region = case_when(
      study_area == "Antarctic" ~ "Polar",
      TRUE ~ "Temperate"
    )
  )%>% 
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
  
  # Weight durations by the inverse of the predicted error
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
                             replace = TRUE
         )) %>% 
  ungroup %>% 
  left_join(All_krill_data_ENPcombined_noSA, by = c("SpeciesCode", "region"))


write_csv(krill_biomass_estimates, "krill_biomass_estimates_for_editing.csv")

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
  daily_rates <- krill_biomass_estimates %>% 
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
krill_daily <- estimate_daily(krill_biomass_estimates, latitudes, aug1, 120)  %>% 
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
           SpeciesCode == "be" ~ "Balaenoptera edeni",
           SpeciesCode == "bs" ~ "Balaenoptera borealis")
  )


# check to see prey consumption rates 
krill_daily %>%  
  group_by(SpeciesCode) %>% 
  filter(SpeciesCode == "bw", region == "Temperate") %>%
  pull(daily_consumption_high_kg) %>%
  summary()


save(krill_daily, file = "daynights_krill.RData")


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
           SpeciesCode == "be" ~ "Balaenoptera edeni")
  )


fish_daily %>%  
  group_by(SpeciesCode) %>% 
  filter(SpeciesCode == "mn") %>%
  pull(daily_consumption_high_kg) %>%
  summary()


save(fish_daily, file = "daynights_fish.RData")


# FIGURES ----



