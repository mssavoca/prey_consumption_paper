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

# stuff from Max to generate data ----

load("rates_CATS_updated_5_14_20.RData") # CATS deployment summary

load("rates_MedTermTag.RData") # Medium-term deployment summary

load("rates[1].RData") # DTAG deployment summary
  
dtag_lunges_long <- dtag_lunges_long %>% 
  mutate(prey_general = "Krill")

All_rorqual_deployments <- dtag_lunges_long %>% 
  bind_rows(lunge_rates)


load("tag_guide4.RData") 
#load("rates.RData")
# load("daynights_krill.RData") # This keeps its original name "krill_daily" here
# load("daynights_fish.RData") # This keeps its original name "fish_daily" here
# load("daynights_AntProj.RData")

load("balaenid_sampled.RData")


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




# Add in Bryde's data, MAX NEEDS TO CHECK THIS----
# Combining day and twilight lunge rates
# adding in Bryde's deployments from 2018
Brydes_lunges <- read_csv("Brydes_2018_lunge_rates.csv") %>% 
  select(ID, study_area, region, prey_general, duration_h, daylight_h, twilight_h, night_h, total_lunges) 


Brydes_lunges_long <- Brydes_lunges %>% 
  pivot_longer(duration_h:night_h, names_to = "Phase", values_to = "Time (hours)") %>% 
  rename(Study_Area = study_area,
         Region = region, 
         Lunges = total_lunges) %>% 
  mutate(prey_general = "Fish",
         tag_type = "CATS",
         Phase = case_when(Phase == "duration_h" ~ "Total",
                           Phase == "daylight_h" ~ "Day",
                           Phase == "twilight_h" ~ "Twilight",
                           Phase == "night_h" ~ "Night"),
         Rate = Lunges/`Time (hours)`) %>% 
  mutate_at(vars(Rate), ~replace(., is.infinite(.), 0))

Brydes_lunges_day <- Brydes_lunges_long %>% 
  filter(Phase %in% c("Day", "Twilight")) %>% 
  group_by(ID, Study_Area, Region, prey_general, tag_type) %>% 
  summarize(
    Rate = weighted.mean(Rate, by = `Time (hours)`),
    `Time (hours)` = sum(`Time (hours)`),
     Lunges = first(Lunges)) %>% 
  mutate(Phase = "Day")

Brydes_lunges_for_join <- Brydes_lunges_long %>%
  filter(Phase == "Night") %>% 
  bind_rows(Brydes_lunges_day) %>% 
  mutate(species_code = "be",
         region = "Temperate") %>% 
  pivot_wider(names_from = Phase, values_from = c(`Time (hours)`, Lunges, Rate)) %>% 
  rename(
    day_lunges = "Lunges_Day",
    night_lunges = "Lunges_Night",
    day_rate = "Rate_Day",
    night_rate = "Rate_Night",
    day_dur = "Time (hours)_Day",
    night_dur = "Time (hours)_Night"
    ) %>% 
  select(-c(Study_Area, tag_type, Region))


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


write.csv(x = pop)


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



  
# load("~/Documents/Research Data/Whale research/ShirelCh2/rates_DTAG only_9.17.19_long.RData")
# 
# lunge_rates_DTAG <- dtag_lunges_long %>% 
#   mutate_at(vars(Rate), ~replace(., is.nan(.), 0)) %>% 
#   rename(`Time (hours)` = Hours) %>% 
#   mutate(prey_general = "Krill",
#          Study_Area = "SoCal",
#          Region = "Eastern North Pacific",
#          tag_type = "DTAG")



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
daynight_rates <- lunge_rates %>%
  rbind(rename(dtag_lunges_long, `Time (hours)` = Hours)) %>% 
  left_join(select(tag_guide4, ID, study_area), by = "ID") %>% 
  mutate(
    region = case_when(
      study_area == "Antarctic" ~ "Polar",
      TRUE ~ "Temperate"
    )
  )%>% 
  group_by(ID, region, prey_general) %>% 
  group_modify(combine_rates) %>% 
  ungroup() %>% 
  mutate(species_code = substr(ID, 1, 2),
         prey_general = ifelse(is.na(prey_general), "Krill", prey_general)) %>% 
  #MODIFIED BY MATT TO ADD BRYDE'S
  bind_rows(Brydes_lunges_for_join)

ss_table_dn <- daynight_rates %>%
  group_by(species_code, region) %>%
  summarise(sample_size = n())


# %>% 
#   # Modified by Matt here
#   rename(SpeciesCode = "species_code") %>% 
#   left_join(whale_lengths, SpeciesCode, Engulfment_m3





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
    group_by(species_code, region) %>% 
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



# day_dist <- daynight_rates %>% 
#   filter(species_code == "bw") %>% 
#   mutate(dur_bucket = cut_width(day_dur, 1, boundary = 0)) %>% 
#   group_by(dur_bucket) %>% 
#   summarise(mean_rate = mean(day_rate),
#             sd_rate = sd(day_rate)) 
# 
# 
# daynight_rates %>% 
#   filter(species_code == "bw") %>% 
#   mutate(dur_bucket = cut_width(day_dur, 1, boundary = 0)) %>% 
#   ggplot(aes(day_dur, day_rate)) + geom_point() + geom_smooth(method = "lm")
# 
# ggplot(day_dist, aes(dur_bucket, mean_rate)) +
#   geom_pointrange(aes(ymin = mean_rate-sd_rate , ymax = mean_rate+sd_rate))




plot(dur, pred_err, type = "l")

dur_sd <- function(dur) {
  Asym <- -1.846
  R0 <- 12.67
  lrc <- -2.057
  
  pred_err <- Asym + (R0 - Asym) * exp(-exp(lrc) * dur)
  pred_err_10 <- Asym + (R0 - Asym) * exp(-exp(lrc) * 10)
  ifelse(dur > 10, pred_err_10, pred_err)/2 
}


estimate_rate_DEC <- function(prey, min_lunges, min_dur) {
  daynight_rates %>%
    filter(prey_general == prey,
           day_lunges + night_lunges >= min_lunges,
           day_dur + night_dur >= min_dur) %>% 
    group_by(species_code, region) %>% 
    group_modify(
      function(data, keys) {
        tibble(i = 1:iter) %>% 
          mutate(
            dep_row = sample(1:nrow(data), size = iter, replace = TRUE),
            mean_day = data$day_rate[dep_row],
            sd_day = dur_sd(data$day_dur[dep_row]),
            day_rate = rnorm(iter, mean_day, sd_day),
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


krill_rate_estimates_DEC <- estimate_rate_DEC("Krill", 1, 1)



krill_biomass_estimates <- krill_rate_estimates %>% 
  rename(SpeciesCode = "species_code") %>% 
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


###### write_csv(krill_biomass_estimates, "krill_biomass_estimates_for_editing.csv")



# Assuming Sep 15, 12.5 hours of daylight
krill_rate_estimates %>% 
  filter(species_code == "bw") %>% 
  mutate(daily_rate = day_rate * 12.5 + night_rate * 11.5) %>% 
  pull(daily_rate) %>% 
  quantile(c(0.25, 0.5, 0.75, 0.95))
# 25%      50%      75%      95% 
#   153.2969 208.8918 289.6198 471.0381


krill_rate_estimates_DEC %>% 
  filter(species_code == "bw") %>% 
  mutate(daily_rate = day_rate * 12.5 + night_rate * 11.5) %>% 
  pull(daily_rate) %>% 
  quantile(c(0.25, 0.5, 0.75, 0.95))
# 25%      50%      75%      95% 
#   140.1731 205.0902 287.1584 457.2846 

weighted_MFC <- krill_rate_estimates %>% 
  filter(species_code == "bw") %>% 
  mutate(daily_rate = day_rate * 12.5 + night_rate * 11.5)

Conf_int_DEC <- krill_rate_estimates_DEC %>% 
  filter(species_code == "bw") %>% 
  mutate(daily_rate = day_rate * 12.5 + night_rate * 11.5)


# two methods both used in meta analyses with dispararte SS. Uncer most condifiotns they are the same.
# we ran dave and Max's methods for the bw data and here are teh ecdf. as you can see the reuslts are the same.
  ggplot(weighted_MFC, aes(daily_rate)) +
    stat_ecdf(color = "red") +
    stat_ecdf(data = Conf_int_DEC, color = "blue") + 
    labs(y = "cumulative probability",
        x = "blue whale daily lunge rate") +
    xlim(0,500) +
    theme_minimal()


krill_rate_estimates %>% 
  filter(species_code == "bw") %>% 
  ggplot(aes(day_rate)) +
  geom_density() +
  theme_minimal()

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
      # daily_biomass_kg_best_low = biomass_best_low_m3*Engulf_cap_m3*daily_rate,
      #      daily_biomass_kg_best_high = biomass_best_high_m3*Engulf_cap_m3*daily_rate)
  
  
  
 # prey_biomass = rlnorm log mean 
  
  # daily_rate_summ <- daily_rates %>% 
  #   group_by(species_code, yday, region) %>% 
  #   summarize(
  #     q25 = quantile(daily_rate, 0.25),
  #     q50 = quantile(daily_rate, 0.50),
  #     q75 = quantile(daily_rate, 0.75),
  #     q95 = quantile(daily_rate, 0.95)
  #   ) %>% 
  #   ungroup()
  # 
  # filter(daily_rate_summ, species_code == "bw") %>% 
  #   ggplot(aes(yday)) +
  #   geom_ribbon(aes(ymin = q25, ymax = q75), fill = "grey80") +
  #   geom_line(aes(y = q50, group = 1)) +
  #   expand_limits(y = 0) +
  #   labs(x = "Julian day", y = "Lunges per day") +
  #   theme_classic()
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


krill_daily %>%  
  group_by(SpeciesCode) %>% 
  filter(SpeciesCode == "bw") %>%
  pull(daily_consumption_high_kg) %>%
  summary()


  ggplot(aes(daily_biomass_kg_best_low)) + 
  geom_density(aes(fill = SpeciesCode), alpha = 0.2) +
  facet_wrap(~SpeciesCode, scales = "free_x") +
  xlim(0,20000)
p

save(krill_daily, file = "daynights_krill.RData")




# Fish feeding information

fish_rate_estimates <- estimate_rate("Fish", 1, 1)


fish_biomass_estimates <- fish_rate_estimates %>% 
  rename(SpeciesCode = "species_code") %>% 
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
  

  # prey_biomass = rlnorm log mean 
  
  # daily_rate_summ <- daily_rates %>% 
  #   group_by(species_code, yday, region) %>% 
  #   summarize(
  #     q25 = quantile(daily_rate, 0.25),
  #     q50 = quantile(daily_rate, 0.50),
  #     q75 = quantile(daily_rate, 0.75),
  #     q95 = quantile(daily_rate, 0.95)
  #   ) %>% 
  #   ungroup()
  # 
  # filter(daily_rate_summ, species_code == "bw") %>% 
  #   ggplot(aes(yday)) +
  #   geom_ribbon(aes(ymin = q25, ymax = q75), fill = "grey80") +
  #   geom_line(aes(y = q50, group = 1)) +
  #   expand_limits(y = 0) +
  #   labs(x = "Julian day", y = "Lunges per day") +
  #   theme_classic()
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
  #filter(prey_general == "Fish", species_code != "bp") %>%
  pull(daily_consumption_high_kg) %>%
  summary()
  
ggplot(aes(daily_consumption_low_kg)) + 
  geom_density(aes(fill = SpeciesCode), alpha = 0.2) +
  facet_wrap(~SpeciesCode, scales = "free_x") +
  xlim(0,10000)
p

save(fish_daily, file = "daynights_fish.RData")

# generating plots and tables for paper ----

pal <- c("B. bonaerensis" = "firebrick3", "B. borealis" = "goldenrod2", "B. edeni" = "darkorchid3",  "M. novaeangliae" = "gray30", "B. physalus" = "chocolate3", "B. musculus" = "dodgerblue2")


#Supplemental figure of whale lengths
All_droned_lengths <- whale_lengths %>% 
  filter(!SpeciesCode %in% c("bs")) %>% 
  #filter(ID != "bw180904-44") %>% 
  group_by(name, SpeciesCode) %>% 
  ungroup %>% 
  #mutate(sp_lbl = factor(abbr_binom(Species)) %>% fct_rev) %>% 
  ggplot(aes(x = fct_relevel(factor(abbr_binom(Species)), "B. bonaerensis", "M. novaeangliae", "B. edeni","B. physalus", "B. musculus"), 
             y = length,
             fill = abbr_binom(Species))) +
  geom_flat_violin(position = position_nudge(x = 0.1, y = 0), alpha = .8) +
  geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
  geom_point(position = position_jitter(width = .05), alpha = 0.6) +
  coord_flip() +
  scale_fill_manual(values = pal) +
  labs(x = "Species",
       y = "Measured length (m)") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"))
All_droned_lengths

dev.copy2pdf(file="All_droned_lengths.pdf", width=12, height=8)


Engulfment_capacity <- whale_lengths %>% 
  filter(!SpeciesCode %in% c("bs")) %>% 
  #filter(ID != "bw180904-44") %>% 
  group_by(name, SpeciesCode) %>% 
  ungroup %>% 
  #mutate(sp_lbl = factor(abbr_binom(Species)) %>% fct_rev) %>% 
  ggplot(aes(x = fct_relevel(factor(abbr_binom(Species)), "B. bonaerensis", "B. edeni","M. novaeangliae", "B. physalus", "B. musculus"), 
             y = Engulfment_m3,
             fill = abbr_binom(Species))) +
  geom_flat_violin(position = position_nudge(x = 0.1, y = 0), alpha = .8) +
  geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
  geom_point(position = position_jitter(width = .05), alpha = 0.6) +
  coord_flip() +
  scale_fill_manual(values = pal) +
  labs(x = "Species",
       y = bquote('Engulfment capacity'~(m^3))) +
  scale_y_log10(labels = scales::comma,
                breaks = c(1,5,10,50,100)) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"))
Engulfment_capacity

dev.copy2pdf(file="Engulfment_capacity.pdf", width=14, height=6.5)

 ##Figure 2----
#summary tables of lunges, water filtered, and prey consumed per day 
summ_stats_all <- fish_daily %>%    # toggle between fish and krill
  group_by(SpeciesCode, region) %>% 
  summarise(
    `Lunges per day` = median(daily_rate),
    Daily_rate_IQR25 = round(quantile(daily_rate, probs = 0.25, na.rm = TRUE), 2), 
    Daily_rate_IQR75 = round(quantile(daily_rate, probs = 0.75, na.rm = TRUE), 2),
    `Water filtered per day` = median(Engulf_cap_m3*daily_rate),
    Water_filtered_IQR25 = round(quantile(Engulf_cap_m3*daily_rate, probs = 0.25, na.rm = TRUE), 2), 
    Water_filtered_IQR75 = round(quantile(Engulf_cap_m3*daily_rate, probs = 0.75, na.rm = TRUE), 2)) %>% 
  unite("Lunges per day IQR", c(Daily_rate_IQR25, Daily_rate_IQR75), sep = "-") %>% 
  unite("Water filtered per day IQR", c(Water_filtered_IQR25, Water_filtered_IQR75), sep = "-")

summ_prey_stats <- fish_daily %>%  # toggle between fish and krill
  #filter(!Region %in% c("North Atlantic", "Chile")) %>% 
  group_by(SpeciesCode, region) %>% 
  summarise(
    # med_daily_consumpt_hyp_low = round(median(prey_mass_per_day_hyp_low_kg), 2),
    # med_daily_consumpt_hyp_low_IQR25 = round(quantile(prey_mass_per_day_hyp_low_kg, probs = 0.25, na.rm = TRUE), 2),
    # med_daily_consumpt_hyp_low_IQR75 = round(quantile(prey_mass_per_day_hyp_low_kg, probs = 0.75, na.rm = TRUE), 2),
    med_daily_consumpt_best = round(median(daily_consumption_low_kg, na.rm = TRUE), 2),
    med_daily_consumpt_best_IQR25 = round(quantile(daily_consumption_low_kg, probs = 0.25, na.rm = TRUE), 2), 
    med_daily_consumpt_best_IQR75 = round(quantile(daily_consumption_low_kg, probs = 0.75, na.rm = TRUE), 2), 
    med_daily_consumpt_top50 = round(median(daily_consumption_high_kg, na.rm = TRUE), 2),
    med_daily_consumpt_top50_IQR25 = round(quantile(daily_consumption_high_kg, probs = 0.25, na.rm = TRUE), 2), 
    med_daily_consumpt_top50_IQR75 = round(quantile(daily_consumption_high_kg, probs = 0.75, na.rm = TRUE), 2)) %>% 
  #unite("Daily consumption lower estimate IQR", c(med_daily_consumpt_hyp_low_IQR25, med_daily_consumpt_hyp_low_IQR75), sep = "-") %>% 
  unite("Daily consumption IQR", c(med_daily_consumpt_best_IQR25, med_daily_consumpt_best_IQR75), sep = "-") %>% 
  unite("Daily consumption Top 50% IQR", c(med_daily_consumpt_top50_IQR25, med_daily_consumpt_top50_IQR75), sep = "-")


summ_EnDens <- fish_daily %>% 
  #filter(!Region %in% c("North Atlantic", "Chile")) %>% 
  group_by(SpeciesCode, region) %>% 
  # Daily energy intake in gigajoules (GJ)
  summarise(
    # med_daily_En_hyp_low = round(median(EnDens_hyp_low/1e6),2),
    # med_daily_En_hyp_low_IQR25 = round(quantile(EnDens_hyp_low/1e6, probs = 0.25, na.rm = TRUE), 2),
    # med_daily_En_hyp_low_IQR75 = round(quantile(EnDens_hyp_low/1e6, probs = 0.75, na.rm = TRUE), 2),
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
  filter(region == "Polar") %>%  
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
       y = bquote('Estimated feeding rate'~(lunges~ind^-1~d^-1))) + 
  #ylim(0, 2000) +
  scale_y_log10(labels = scales::comma, limits = c(50, 3000), 
                breaks = c(10,100,250,500,1000,2500)) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic"))
Daily_rate_Antarctic

dev.copy2pdf(file="Daily_rate_Antarctic.pdf", width=16, height=6)


Daily_filtration_Antarctic <- krill_daily %>% 
  filter(region == "Polar") %>%   
  ggplot(aes(x = fct_reorder(abbr_binom(Species), Engulf_cap_m3*daily_rate), y = Engulf_cap_m3*daily_rate, 
             fill = abbr_binom(Species))) +
  geom_flat_violin(position = position_nudge(x = 0.1, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  coord_flip() + 
  scale_fill_manual(values = pal) +
  scale_y_log10(labels = scales::comma, breaks = c(100,1000,5000,10000,50000,100000)) +
  #ylim(0,60000) +
  labs(x = "Species",
       y = bquote('Estimated water filtered'~(m^3~ind^-1~d^-1))) + 
  theme_classic(base_size = 24) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic"))
Daily_filtration_Antarctic

dev.copy2pdf(file="Daily_filtration_Antarctic.pdf", width=16, height=6)



Daily_biomass_ingested_Antarctic <- krill_daily %>% 
  filter(region == "Polar") %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. bonaerensis")) %>% 
  
  ggplot() +
  # #best-low *our best estimate*
geom_flat_violin(aes(x = Species, y = daily_consumption_low_kg/1000,
                     fill = Species), 
                 position = position_nudge(x = 0.2, y = 0), alpha = 0.5, adjust = 2) +
geom_boxplot(aes(x = Species, y = daily_consumption_low_kg/1000,
                 fill = abbr_binom(Species)), width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5,
             position = position_nudge(x = 0.12, y = 0)) +
# best-upper
geom_flat_violin(aes(x = Species, y = daily_consumption_high_kg/1000,
                     fill = Species), 
                 position = position_nudge(x = 0.2, y = 0), alpha = 0.5, adjust = 2) +
geom_boxplot(aes(x = abbr_binom(Species), y = daily_consumption_high_kg/1000,
                 fill = abbr_binom(Species)), width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +

  coord_flip() +
  scale_fill_manual(values = pal) +
  scale_y_log10(labels = scales::comma, limits = c(0.05, 15), breaks = c(0,1,5,10,15)) +
  labs(x = "Species",
       y = bquote('Estimated prey consumed'~(tonnes~ind^-1~d^-1))) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic"))
Daily_biomass_ingested_Antarctic

dev.copy2pdf(file="Daily_biomass_ingested_Antarctic.pdf", width=16, height=6)


Fig_2_Antarctic <- ggarrange(Daily_rate_Antarctic, Daily_filtration_Antarctic, Daily_biomass_ingested_Antarctic, 
                             #labels = c("A", "B", "C"), # THIS IS SO COOL!!
                             font.label = list(size = 18),
                             legend = "none",
                             ncol = 1, nrow = 3)
Fig_2_Antarctic

dev.copy2pdf(file="Fig_2_Antarctic.pdf", width=12, height=18)




#NON-ANTARCTIC
Daily_rate_nonAntarctic <- krill_daily %>% 
  filter(region == "Temperate") %>%  
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
       y = bquote('Estimated feeding rate'~(lunges~ind^-1~d^-1))) + 
  #ylim(0, 2000) +
  scale_y_log10(labels = scales::comma,
                limits = c(10, 1500), 
                breaks = c(10,100,250,500,1000)) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic"))
Daily_rate_nonAntarctic

dev.copy2pdf(file="Daily_rate_nonAntarctic.pdf", width=16, height=6)


Daily_filtration_nonAntarctic <- krill_daily %>% 
  filter(region == "Temperate") %>%   
  ggplot(aes(x = fct_reorder(abbr_binom(Species), Engulf_cap_m3*daily_rate), y = Engulf_cap_m3*daily_rate, 
             fill = abbr_binom(Species))) +
  geom_flat_violin(position = position_nudge(x = 0.1, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  coord_flip() + 
  scale_fill_manual(values = pal) +
  scale_y_log10(labels = scales::comma,
                limits = c(250, 60000),
                breaks = c(500,1000,5000,10000,50000)) +
  labs(x = "Species",
       y = bquote('Estimated water filtered'~(m^3~ind^-1~d^-1))) + 
  theme_classic(base_size = 24) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic"))
Daily_filtration_nonAntarctic

dev.copy2pdf(file="Daily_filtration_nonAntarctic.pdf", width=16, height=6)



Daily_biomass_ingested_nonAntarctic <- krill_daily %>% 
  filter(region == "Temperate") %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. physalus", "M. novaeangliae")) %>% 
  
  ggplot() +
  # #best-low *our best estimate*
  geom_flat_violin(aes(x = Species, y = daily_consumption_low_kg/1000,
                       fill = Species), 
                   position = position_nudge(x = 0.2, y = 0), alpha = 0.5, adjust = 2) +
  geom_boxplot(aes(x = Species, y = daily_consumption_low_kg/1000,
                   fill = abbr_binom(Species)), width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5,
               position = position_nudge(x = 0.12, y = 0)) +
  # best-upper
  geom_flat_violin(aes(x = Species, y = daily_consumption_high_kg/1000,
                       fill = Species), 
                   position = position_nudge(x = 0.2, y = 0), alpha = 0.5, adjust = 2) +
  geom_boxplot(aes(x = abbr_binom(Species), y = daily_consumption_high_kg/1000,
                   fill = abbr_binom(Species)), width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  
  coord_flip() +
  scale_fill_manual(values = pal) +
  scale_y_log10(labels = scales::comma, 
                limits = c(0.1, 60), 
                breaks = c(0,1,5,10,25,50)) +
  labs(x = "Species",
       y = bquote('Estimated prey consumed'~(tonnes~ind^-1~d^-1))) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic"))
Daily_biomass_ingested_nonAntarctic

dev.copy2pdf(file="Daily_biomass_ingested_nonAntarctic.pdf", width=16, height=6)


Fig_2_nonAntarctic <- ggarrange(Daily_rate_nonAntarctic, Daily_filtration_nonAntarctic, Daily_biomass_ingested_nonAntarctic, 
                                #labels = c("D", "E", "F"), # THIS IS SO COOL!!
                                font.label = list(size = 18),
                                legend = "none",
                                ncol = 1, nrow = 3)
Fig_2_nonAntarctic

dev.copy2pdf(file="Fig_2_nonAntarctic.pdf", width=12, height=18)



# FISH-FEEDERS, Individual/daily
Daily_rate_fish <- fish_daily %>% 
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
       y = bquote('Estimated feeding rate'~(lunges~ind^-1~d^-1))) + 
  #ylim(0, 2000) +
  scale_y_log10(labels = scales::comma,
                limits = c(1, 500), 
                breaks = c(1,10,100,250,500)) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic"))
Daily_rate_fish

dev.copy2pdf(file="Daily_rate_fish.pdf", width=16, height=6)


Daily_filtration_fish <- fish_daily %>% 
  filter(region == "Temperate") %>%   
  ggplot(aes(x = fct_reorder(abbr_binom(Species), Engulf_cap_m3*daily_rate), y = Engulf_cap_m3*daily_rate, 
             fill = abbr_binom(Species))) +
  geom_flat_violin(position = position_nudge(x = 0.1, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  coord_flip() + 
  scale_fill_manual(values = pal) +
  scale_y_log10(labels = scales::comma,
                limits = c(25, 15000),
                breaks = c(100,500,1000,5000,10000)) +
  labs(x = "Species",
       y = bquote('Estimated water filtered'~(m^3~ind^-1~d^-1))) + 
  theme_classic(base_size = 24) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic"))
Daily_filtration_fish

dev.copy2pdf(file="Daily_filtration_fish.pdf", width=16, height=6)



Daily_biomass_ingested_fish <- fish_daily %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "M. novaeangliae")) %>% 
  
  ggplot() +
  # #best-low *our best estimate*
  geom_flat_violin(aes(x = Species, y = daily_consumption_low_kg/1000,
                       fill = Species), 
                   position = position_nudge(x = 0.2, y = 0), alpha = 0.5, adjust = 2) +
  geom_boxplot(aes(x = Species, y = daily_consumption_low_kg/1000,
                   fill = abbr_binom(Species)), width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5,
               position = position_nudge(x = 0.12, y = 0)) +
  # best-upper
  geom_flat_violin(aes(x = Species, y = daily_consumption_high_kg/1000,
                       fill = Species), 
                   position = position_nudge(x = 0.2, y = 0), alpha = 0.5, adjust = 2) +
  geom_boxplot(aes(x = abbr_binom(Species), y = daily_consumption_high_kg/1000,
                   fill = abbr_binom(Species)), width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  
  coord_flip() +
  scale_fill_manual(values = pal) +
  scale_y_log10(labels = scales::comma, 
                limits = c(0.1, 60), 
                breaks = c(0,1,5,10,25,50)) +
  labs(x = "Species",
       y = bquote('Estimated prey consumed'~(tonnes~ind^-1~d^-1))) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic"))
Daily_biomass_ingested_fish

dev.copy2pdf(file="Daily_biomass_ingested_fish.pdf", width=16, height=6)


Fig_2_fish <- ggarrange(Daily_rate_fish, Daily_filtration_fish, Daily_biomass_ingested_fish, 
                                #labels = c("D", "E", "F"), # THIS IS SO COOL!!
                                font.label = list(size = 18),
                                legend = "none",
                                ncol = 1, nrow = 3)
Fig_2_fish

dev.copy2pdf(file="Fig_2_fish.pdf", width=12, height=18)






#summary tables of individual water filtered and prey consumed per year (Figure 3)----
krill_yearly_ind <- krill_daily %>% 
  mutate(filtration60 = Engulf_cap_m3*daily_rate*60,
         filtration90 = Engulf_cap_m3*daily_rate*90,
         filtration120 = Engulf_cap_m3*daily_rate*120,
         filtration182.5 = Engulf_cap_m3*daily_rate*182.5,
         prey60 = daily_biomass_kg_best_high*60,
         prey90 = daily_biomass_kg_best_high*90,
         prey120 = daily_biomass_kg_best_high*120,
         prey182.5 = daily_biomass_kg_best_high*182.5)

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







# YEARLY INDIVIDUAL FILTRATION, ANTARCTIC
Annual_filtration_ind_Antarctic <- Annual_filtfeed %>% 
  filter(region == "Polar") %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. bonaerensis")) %>% 
  ggplot() +
  geom_ribbon(aes(ymin = IQR25_filt_yr, ymax = IQR75_filt_yr, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_line(aes(days_feeding, med_filt_yr, color = Species), 
            size = 1.5) +
  scale_colour_manual(values = pal) +
  facet_grid(~Species) +   # Can add region in here if desired
  labs(x = "Days feeding",
       y = bquote('Projected water filtered'~(m^3~ind^-1~yr^-1))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_log10(labels = scales::comma, limits = c(8e4, 8e6),
                breaks = c(100000,500000,1e6,5e6)) +
  theme_minimal(base_size = 24) +
    theme(strip.text = element_text(face = "italic"),
          legend.position = "none")
Annual_filtration_ind_Antarctic

dev.copy2pdf(file="Annual_filtration_ind_Antarctic.pdf", width=14, height=8)



# YEARLY INDIVIDUAL FILTRATION, NON-ANTARCTIC
Annual_filtration_ind_nonAntarctic <- Annual_filtfeed %>% 
  filter(region == "Temperate") %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. physalus", "M. novaeangliae")) %>% 
  ggplot() +
  geom_ribbon(aes(ymin = IQR25_filt_yr, ymax = IQR75_filt_yr, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_line(aes(days_feeding, med_filt_yr, color = Species), 
            size = 1.5) +
  scale_colour_manual(values = pal) +
  facet_grid(~Species) +   # Can add region in here if desired
  labs(x = "Days feeding",
       y = bquote('Projected water filtered'~(m^3~ind^-1~yr^-1))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_log10(labels = scales::comma, limits = c(8e4, 8e6),
                breaks = c(100000,500000,1e6,5e6)) +
  theme_minimal(base_size = 24) + 
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none")
Annual_filtration_ind_nonAntarctic

dev.copy2pdf(file="Annual_filtration_ind_nonAntarctic.pdf", width=14, height=8)

# And now for the prey
# YEARLY INDIVIDUAL PREY, ANTARCTIC
Annual_ingestion_ind_Antarctic <- Annual_filtfeed %>% 
  filter(region == "Polar") %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. bonaerensis")) %>% 
  ggplot() +
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_yr, ymax = IQR75_ingest_low_t_yr, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_ribbon(aes(ymin = IQR25_ingest_high_t_yr, ymax = IQR75_ingest_high_t_yr, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_line(aes(days_feeding, med_ingest_low_t_yr, color = Species), 
            size = 1.5) +
  geom_line(aes(days_feeding, med_ingest_high_t_yr, color = Species), 
            size = 1.5, linetype = "dashed") +
  scale_colour_manual(values = pal) +
  facet_grid(~Species) +   # Can add region in here if desired
  labs(x = "Days feeding",
       y = bquote('Projected prey consumption'~(tonnes~ind^-1~yr^-1))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_log10(labels = scales::comma, limits = c(10, 4000),
                breaks = c(30,100,300,600,1500,3000)) +
  theme_minimal(base_size = 24) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none")
Annual_ingestion_ind_Antarctic

dev.copy2pdf(file="Annual_ingestion_ind_Antarctic.pdf", width=14, height=8)



# YEARLY INDIVIDUAL PREY, NON-ANTARCTIC
Annual_ingestion_ind_nonAntarctic <- Annual_filtfeed %>% 
  filter(region == "Temperate") %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. physalus", "M. novaeangliae")) %>% 
  ggplot() +
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_yr, ymax = IQR75_ingest_low_t_yr, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_ribbon(aes(ymin = IQR25_ingest_high_t_yr, ymax = IQR75_ingest_high_t_yr, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_line(aes(days_feeding, med_ingest_low_t_yr, color = Species), 
            size = 1.5) +
  geom_line(aes(days_feeding, med_ingest_high_t_yr, color = Species), 
            size = 1.5, linetype = "dashed") +
  scale_colour_manual(values = pal) +
  facet_grid(~Species) +   # Can add region in here if desired
  labs(x = "Days feeding",
       y = bquote('Projected prey consumption'~(tonnes~ind^-1~yr^-1))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_log10(labels = scales::comma, limits = c(10, 4000),
                breaks = c(30,100,300,600,1500,3000)) +
  theme_minimal(base_size = 24) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none")
Annual_ingestion_ind_nonAntarctic

dev.copy2pdf(file="Annual_ingestion_ind_nonAntarctic.pdf", width=14, height=8)


# Figure 4----
# YEARLY POPULATION FILTRATION, NON-ANTARCTIC
Annual_filtration_Pop_nonAntarctic <- Annual_filtfeed %>% 
  filter(region == "Temperate") %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. musculus", "M. novaeangliae")) %>% 
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
       y = bquote('Projected water filtered'~(km^3~yr^-1))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_log10(labels = scales::comma, limits = c(1, 2500),
                breaks = c(1,5,10,25,100,250,1000,2500)) +
  theme_minimal(base_size = 24) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none")
Annual_filtration_Pop_nonAntarctic

dev.copy2pdf(file="Annual_filtration_Pop_nonAntarctic.pdf",  width=14, height=8)


Total_filtered <- Annual_filtfeed %>%
  group_by(SpeciesCode, region) %>%
  filter(days_feeding == 100) %>% 
  summarise(sp_annual_avg_curr = med_filt_curr/1e9,
            sp_annual_avg_hist = med_filt_hist/1e9)


# Northern Hemisphere prey, current
Annual_ingestion_PopCurr_nonAntarctic <- Annual_filtfeed %>% 
  filter(region == "Temperate") %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. musculus", "M. novaeangliae")) %>% 
  ggplot() +
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_curr/1e6, ymax = IQR75_ingest_low_t_curr/1e6, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_ribbon(aes(ymin = IQR25_ingest_high_t_curr/1e6, ymax = IQR75_ingest_high_t_curr/1e6, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_line(aes(days_feeding, med_ingest_low_t_curr/1e6, color = Species), 
            size = 1.5) +
  geom_line(aes(days_feeding, med_ingest_high_t_curr/1e6, color = Species), 
            size = 1.5, linetype = "dashed") +
  scale_colour_manual(values = pal) +
  facet_grid(~Species) +   # Can add region in here if desired
  labs(x = "Days feeding",
       y = bquote('Projected prey consumption'~(Mt~yr^-1))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_log10(limits = c(1.5, 175),
                breaks = c(5, 10, 25, 50, 100)) +
  theme_minimal(base_size = 24) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none")
Annual_ingestion_PopCurr_nonAntarctic

dev.copy2pdf(file="Annual_ingestion_PopCurr_nonAntarctic.pdf", width=14, height=9)


# Northern Hemisphere prey, historic
Annual_ingestion_PopHist_nonAntarctic <- Annual_filtfeed %>% 
  filter(region == "Temperate") %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. musculus", "M. novaeangliae")) %>% 
  ggplot() +
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_hist/1e6, ymax = IQR75_ingest_low_t_hist/1e6, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_ribbon(aes(ymin = IQR25_ingest_high_t_hist/1e6, ymax = IQR75_ingest_high_t_hist/1e6, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_line(aes(days_feeding, med_ingest_low_t_hist/1e6, color = Species), 
            size = 1.5) +
  geom_line(aes(days_feeding, med_ingest_high_t_hist/1e6, color = Species), 
            size = 1.5, linetype = "dashed") +
  scale_colour_manual(values = pal) +
  facet_grid(~Species) +   # Can add region in here if desired
  labs(x = "Days feeding",
       y = bquote('Projected prey consumption'~(Mt~yr^-1))) +  
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_log10(limits = c(1.5, 175),
                breaks = c(5, 10, 25, 50, 100)) +
  theme_minimal(base_size = 24) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none")
Annual_ingestion_PopHist_nonAntarctic

dev.copy2pdf(file="Annual_ingestion_PopHist_nonAntarctic.pdf", width=14, height=9)


# Northern Hemisphere prey, COMBINED
Annual_ingestion_PopComb_nonAntarctic <- Annual_filtfeed %>% 
  filter(region == "Temperate") %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. musculus", "M. novaeangliae")) %>% 
  ggplot() +
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
  geom_line(aes(days_feeding, med_ingest_high_t_curr/1e6, color = Species), 
            size = 1.5, linetype = "dashed") +
  
  
  geom_line(aes(days_feeding, med_ingest_low_t_hist/1e6, color = Species), 
            size = 1.5) +
  geom_line(aes(days_feeding, med_ingest_high_t_hist/1e6, color = Species), 
            size = 1.5, linetype = "dashed") +
  
  scale_colour_manual(values = pal) +
  facet_grid(~Species) +   # Can add region in here if desired
  labs(x = "Days feeding",
       y = bquote('Projected prey consumption'~(Mt~yr^-1))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_log10(labels = scales::comma, limits = c(0.1, 400),
                breaks = c(1, 5, 10, 50, 100, 200, 400)) +
  theme_minimal(base_size = 24) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none")
Annual_ingestion_PopComb_nonAntarctic

dev.copy2pdf(file="Annual_ingestion_PopComb_nonAntarctic.pdf", width=14, height=8)



# Annual Northern Hemisphere population krill prey, summary table
# Northern Hemisphere current populations prey consumption (low and high)
summ_prey_annual_curr_pop_stats <- krill_yearly_pop %>% 
  filter(region != "Polar") %>% 
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


# Northern Hemisphere historical populations prey consumption (low and high)
summ_prey_annual_hist_pop_stats <- krill_yearly_pop %>% 
  filter(region != "Polar") %>% 
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



# Southern Hemisphere calculations (Figure 4)----



krill_Ant_projection <- read_csv("krill_biomass_estimates_for_editing.csv")


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
  # daily_biomass_kg_best_low = biomass_best_low_m3*Engulf_cap_m3*daily_rate,
  #      daily_biomass_kg_best_high = biomass_best_high_m3*Engulf_cap_m3*daily_rate)
  
  
  
  
  # prey_biomass = rlnorm log mean 
  
  # daily_rate_summ <- daily_rates %>% 
  #   group_by(species_code, yday, region) %>% 
  #   summarize(
  #     q25 = quantile(daily_rate, 0.25),
  #     q50 = quantile(daily_rate, 0.50),
  #     q75 = quantile(daily_rate, 0.75),
  #     q95 = quantile(daily_rate, 0.95)
  #   ) %>% 
  #   ungroup()
  # 
  # filter(daily_rate_summ, species_code == "bw") %>% 
  #   ggplot(aes(yday)) +
  #   geom_ribbon(aes(ymin = q25, ymax = q75), fill = "grey80") +
  #   geom_line(aes(y = q50, group = 1)) +
  #   expand_limits(y = 0) +
  #   labs(x = "Julian day", y = "Lunges per day") +
  #   theme_classic()
}

# MATT: Change krill_rate_estimates to krill_biomass_estimates. Should have
# columns day_biomass and night_biomass.

aug1 <- 213
krill_daily_Ant_projection <- estimate_daily(krill_Ant_projection, latitudes, aug1, 120)  %>% 
  mutate(Total_energy_intake_best_low_kJ = case_when(region == "Polar" ~ daily_consumption_low_kg*4575,
                                                     region == "Temperate" ~ daily_consumption_low_kg*3628),
         Total_energy_intake_best_high_kJ = case_when(region == "Polar" ~ daily_consumption_high_kg*4575,
                                                      region == "Temperate" ~ daily_consumption_high_kg*3628),
         Mass_specifc_energy_intake_best_high_kJ = Total_energy_intake_best_high_kJ/Mass_est_kg,
         Mass_specifc_energy_intake_best_low_kJ = Total_energy_intake_best_low_kJ/Mass_est_kg
  )


krill_daily_Ant_projection <- krill_daily_Ant_projection %>%  
  mutate(Species = case_when(
    SpeciesCode == "bw" ~ "Balaenoptera musculus",
    SpeciesCode == "bp" ~ "Balaenoptera physalus",
    SpeciesCode == "mn" ~ "Megaptera novaeangliae",
    SpeciesCode == "bb" ~ "Balaenoptera bonaerensis", 
    SpeciesCode == "be" ~ "Balaenoptera edeni",
    SpeciesCode == "bs" ~ "Balaenoptera borealis")) 

krill_daily_Ant_projection %>% 
  group_by(Species) %>% 
  filter(Species == "Balaenoptera musculus") %>%
  pull(daily_consumption_high_kg) %>%
  summary()

save(krill_daily_Ant_projection, file = "daynights_AntProj.RData")



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
       y = bquote('Projected water filtered'~(km^3~yr^-1))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_log10(labels = scales::comma, limits = c(1, 2500),
                breaks = c(1,5,10,25,100,250,1000,2500)) +
  theme_minimal(base_size = 24) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none")
Annual_filtration_Pop_Antarctic

dev.copy2pdf(file="Annual_filtration_Pop_Antarctic.pdf", width=14, height=8)


# Southern Hemisphere prey, current
Annual_ingestion_PopCurr_Antarctic <- Annual_filtfeed_Ant_projection %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. bonaerensis", "M. novaeangliae", "B. physalus")) %>% 
  ggplot() +
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_curr/1e6, ymax = IQR75_ingest_low_t_curr/1e6, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_ribbon(aes(ymin = IQR25_ingest_high_t_curr/1e6, ymax = IQR75_ingest_high_t_curr/1e6, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_line(aes(days_feeding, med_ingest_low_t_curr/1e6, color = Species), 
            size = 1.5) +
  geom_line(aes(days_feeding, med_ingest_high_t_curr/1e6, color = Species), 
            size = 1.5, linetype = "dashed") +
  scale_colour_manual(values = pal) +
  facet_grid(~Species) +   # Can add region in here if desired
  labs(x = "Days feeding",
       y = bquote('Projected prey consumption'~(Mt~yr^-1))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_log10(labels = scales::comma, limits = c(0.1, 400),
                breaks = c(1, 5, 10, 50, 100, 200, 400)) +
  theme_minimal(base_size = 24) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none")
Annual_ingestion_PopCurr_Antarctic

dev.copy2pdf(file="Annual_ingestion_PopCurr_Antarctic.pdf", width=14, height=9)



# Southern Hemisphere prey, historic
Annual_ingestion_PopHist_Antarctic <- Annual_filtfeed_Ant_projection %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. bonaerensis", "M. novaeangliae", "B. physalus")) %>% 
  ggplot() +
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_hist/1e6, ymax = IQR75_ingest_low_t_hist/1e6, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_ribbon(aes(ymin = IQR25_ingest_high_t_hist/1e6, ymax = IQR75_ingest_high_t_hist/1e6, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_line(aes(days_feeding, med_ingest_low_t_hist/1e6, color = Species), 
            size = 1.5) +
  geom_line(aes(days_feeding, med_ingest_high_t_hist/1e6, color = Species), 
            size = 1.5, linetype = "dashed") +
  scale_colour_manual(values = pal) +
  facet_grid(~Species) +   # Can add region in here if desired
  labs(x = "Days feeding",
       y = bquote('Projected prey consumption'~(Mt~yr^-1))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_log10(labels = scales::comma, limits = c(0.1, 400),
                breaks = c(1, 5, 10, 50, 100, 200, 400)) +
  theme_minimal(base_size = 24) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none")
Annual_ingestion_PopHist_Antarctic

dev.copy2pdf(file="Annual_ingestion_PopHist_Antarctic.pdf", width=14, height=8)



# Southern Hemisphere prey, COMBINED
Annual_ingestion_PopComb_Antarctic <- Annual_filtfeed_Ant_projection %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. bonaerensis", "M. novaeangliae", "B. physalus")) %>% 
  ggplot() +
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
  geom_line(aes(days_feeding, med_ingest_high_t_curr/1e6, color = Species), 
            size = 1.5, linetype = "dashed") +
  
  
  geom_line(aes(days_feeding, med_ingest_low_t_hist/1e6, color = Species), 
            size = 1.5) +
  geom_line(aes(days_feeding, med_ingest_high_t_hist/1e6, color = Species), 
            size = 1.5, linetype = "dashed") +
  scale_colour_manual(values = pal) +
  facet_grid(~Species) +   # Can add region in here if desired
  labs(x = "Days feeding",
       y = bquote('Projected prey consumption'~(Mt~yr^-1))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_log10(labels = scales::comma, limits = c(0.1, 400),
                breaks = c(1, 5, 10, 50, 100, 200, 400)) +
  theme_minimal(base_size = 24) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none")
Annual_ingestion_PopComb_Antarctic

dev.copy2pdf(file="Annual_ingestion_PopComb_Antarctic.pdf", width=14, height=8)


col_summ <- function(tbl, f) {
  rbind(
    tbl,
    map(tbl, ~ if(is.numeric(.x)) f(.x) else NA)
  )
}



summ_SO_krill_surplus <- Annual_filtfeed_Ant_projection %>% 
  group_by(SpeciesCode, region) %>% 
  summarise(
    
    Med_low_curr = median(med_ingest_low_t_curr, na.rm = TRUE),
    IQR25_low_t_curr = round(quantile(IQR25_ingest_low_t_curr, probs = 0.25, na.rm = TRUE), 4),
    IQR75_low_t_curr = round(quantile(IQR75_ingest_low_t_curr, probs = 0.25, na.rm = TRUE), 4),
    
    Med_high_curr = median(med_ingest_high_t_curr, na.rm = TRUE),
    IQR25_high_t_curr = round(quantile(IQR25_ingest_high_t_curr, probs = 0.25, na.rm = TRUE), 4),
    IQR75_high_t_curr = round(quantile(IQR75_ingest_high_t_curr, probs = 0.25, na.rm = TRUE), 4),
    
    Med_low_hist = median(med_ingest_low_t_hist, na.rm = TRUE),
    IQR25_low_t_hist = round(quantile(IQR25_ingest_low_t_hist, probs = 0.25, na.rm = TRUE), 4),
    IQR75_low_t_hist = round(quantile(IQR75_ingest_low_t_hist, probs = 0.25, na.rm = TRUE), 4),
    
    Med_high_hist = median(med_ingest_high_t_hist, na.rm = TRUE),
    IQR25_high_t_hist = round(quantile(IQR25_ingest_high_t_hist, probs = 0.25, na.rm = TRUE), 4),
    IQR75_high_t_hist = round(quantile(IQR75_ingest_high_t_hist, probs = 0.25, na.rm = TRUE), 4)
    
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



# Krill surplus plot----
krill_surplus <- krill_daily_Ant_projection %>% 
  left_join(pop_data, by = "SpeciesCode") %>% 
  mutate(
    yearly_consumption_curr = (daily_consumption_low_kg/1000)*`Southern hemisphere population estimate (Christensen 2006)`*120,
    yearly_consumption_hist = (daily_consumption_low_kg/1000)*`Southern hemisphere historic estimate (Christensen 2006)`*120
    
  ) %>%
 
   ggplot() +
  geom_flat_violin(aes(region, yearly_consumption_curr), alpha = .5) +
  geom_flat_violin(aes(region, yearly_consumption_hist), alpha = .5) +
  coord_flip() +
  scale_y_log10(labels = scales::comma) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"))

krill_surplus
  
  
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_hist/1e6, ymax = IQR75_ingest_low_t_hist/1e6, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_ribbon(aes(ymin = IQR25_ingest_high_t_hist/1e6, ymax = IQR75_ingest_high_t_hist/1e6, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_line(aes(days_feeding, med_ingest_low_t_hist/1e6, color = Species), 
            size = 1.5) +
  geom_line(aes(days_feeding, med_ingest_high_t_hist/1e6, color = Species), 
            size = 1.5, linetype = "dashed") +
  scale_colour_manual(values = pal) +
  facet_grid(~Species) +   # Can add region in here if desired
  labs(x = "Days feeding",
       y = bquote('Estimated prey consumed'~(Mt~yr^-1))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_log10(labels = scales::comma, limits = c(0.1, 400),
                breaks = c(1, 5, 10, 50, 100, 200, 400)) +
  theme_minimal(base_size = 24) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none")
Annual_ingestion_PopHist_Antarctic




#summary tables of population nutrients recycled per year (Figure 5) ----

# Just doing Southern Ocean iron for now...

rnorm_trunc <- function(n, mean = 0, sd = 1, lower = 0, upper = Inf) {
  result <- rnorm(n, mean, sd)
  result[result < lower] <- lower
  result[result > upper] <- upper
  result
}


# Fe_total_dw <- function(x) {
#   (x*0.8)*0.25*0.146 # Projection of the total weight of feces, multiplied by 0.25 to convert to dry weight, and the 0.146 is the amount of iron (in kg!!!!) per TONNE of whale feces (converted from 0.000146 kg Fe per kg of feces)
# } # From Doughty et al. 2016, Roman and McCarthy 2010 

Fe_total_dw <- function(x) {
  (x*runif(length(x), 0.7, 0.9))*   # saying the amount of iron excreted ranges uniformly betweek 70-90% ingested, see Ratnarajah et al. 2016
    0.25*
    rnorm(length(x), 0.146, 0.135) # Projection of the total weight of feces, multiplied by 0.25 to convert to dry weight, and the 0.146 is the amount of iron (in kg!!!!) per TONNE of whale feces (converted from 0.000146 kg Fe per kg of feces)
} # average and sd of Fe conc in whale feces from: Ratnarajah et al. 2014  


## THESE ARE NOT CORRECT ANYMORE
# N_calc_kg <- function(x) {
#   x*0.8*0.25*0.0056 # THINK ABOUT THE 0.25
# } # see Roman et al. 2016, converted moles PON to g to kg
# 
# P_calc_kg <- function(x) {
#   x*0.8*0.0089
# } # whale fecal average from Ratnarajah et al. 2014 NEED TO CHECK




Fe_t_to_C_export_t <- function(x) {
  x * 1e6 *             #convert tonnes to grams
    0.01791 *           #convert grams iron to moles iron
    5e4 * 12.0107 /     #convert to moles carbon exported (Lavery et al. 2010), convert back to grams carbon
    1e6                 #convert to tonnes carbon
}



Annual_filtfeed_Ant_projection_Nutrients <- Annual_filtfeed_Ant_projection %>%
  
  mutate_at(vars(c("med_ingest_low_t_curr":"IQR75_ingest_high_t_curr",
                   "med_ingest_low_t_hist":"IQR75_ingest_high_t_hist")), 
            .funs = list(Fe = ~Fe_total_dw(.)/1000)) %>% # need to divide by 1000 to get the total iron recycled in tonnes
  
  mutate_at(vars(c("med_ingest_low_t_curr_Fe":"IQR75_ingest_high_t_curr_Fe",
                   "med_ingest_low_t_hist_Fe":"IQR75_ingest_high_t_hist_Fe")), 
            .funs = list(C_produced_Mt = ~(Fe_t_to_C_export_t(.)*3)/1e6)) %>% # Carbon production (~3x what is exported) stimulated by Fe defecation in Mt C
  
  mutate_at(vars(c("med_ingest_low_t_curr_Fe":"IQR75_ingest_high_t_curr_Fe",
                   "med_ingest_low_t_hist_Fe":"IQR75_ingest_high_t_hist_Fe")), 
            .funs = list(C_exported_Mt = ~Fe_t_to_C_export_t(.)*
                           runif(length(.), 0.6, 1)/1e6)) 
  
## estimated C respired is as simple as:
# pre-whaling
quantile(runif(10000, 0.003, 0.006)*4760) # amount respired is 0.03-0.06% of that, see van Franeker et al 1997; PPR are my calculations

# post-whaling
quantile(runif(10000, 0.003, 0.006)*580) # amount respired is 0.03-0.06% of that, see van Franeker et al 1997; PPR are my calculations


  # mutate_at(vars(c("med_ingest_low_t_curr":"IQR75_ingest_high_t_curr",
  #                  "med_ingest_low_t_hist":"IQR75_ingest_high_t_hist")), 
  #           .funs = list(C_respired_Mt = ~(.*(0.1*
  #                                               runif(length(.), 0.25, 0.75))/1e6)))  # Carbon respired by populations in Mt C





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
665e12/0.486 * # number of krill in entire pop (first number is biomass in Mt); at 0.486g per krill (Atkinson et al. 2009)
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


# Southern Hemisphere Fe, current
Annual_ingestion_PopCurr_Antarctic_Fe <- Annual_filtfeed_Ant_projection_Nutrients %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. bonaerensis", "M. novaeangliae", "B. physalus")) %>% 
  ggplot() +
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_curr_Fe, ymax = IQR75_ingest_low_t_curr_Fe, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_ribbon(aes(ymin = IQR25_ingest_high_t_curr_Fe, ymax = IQR75_ingest_high_t_curr_Fe, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_line(aes(days_feeding, med_ingest_low_t_curr_Fe, color = Species), 
            size = 1.5) +
  geom_line(aes(days_feeding, med_ingest_high_t_curr_Fe, color = Species), 
            size = 1.5, linetype = "dashed") +
  scale_colour_manual(values = pal) +
  facet_grid(~Species) +   # Can add region in here if desired
  labs(x = "Days feeding",
       y = bquote('Estimated Fe recycled'~(tonnes~yr^-1))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_log10(labels = scales::comma, limits = c(1, 5000),
                breaks = c(1,10, 100, 500, 1000, 2500, 5000)) +
  theme_minimal(base_size = 24) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none")
Annual_ingestion_PopCurr_Antarctic_Fe

dev.copy2pdf(file="Annual_ingestion_PopCurr_Antarctic_Fe.pdf", width=14, height=9)




# Southern Hemisphere Fe, historic
Annual_ingestion_PopHist_Antarctic_Fe <- Annual_filtfeed_Ant_projection_Nutrients %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. bonaerensis", "M. novaeangliae", "B. physalus")) %>% 
  ggplot() +
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_hist_Fe, ymax = IQR75_ingest_low_t_hist_Fe, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_ribbon(aes(ymin = IQR25_ingest_high_t_hist_Fe, ymax = IQR75_ingest_high_t_hist_Fe, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_line(aes(days_feeding, med_ingest_low_t_hist_Fe, color = Species), 
            size = 1.5) +
  geom_line(aes(days_feeding, med_ingest_high_t_hist_Fe, color = Species), 
            size = 1.5, linetype = "dashed") +
  scale_colour_manual(values = pal) +
  facet_grid(~Species) +   # Can add region in here if desired
  labs(x = "Days feeding",
       y = bquote('Projected Fe recycled'~(tonnes~yr^-1))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_log10(labels = scales::comma, limits = c(1, 12000),
                breaks = c(1,10, 100, 500, 1000, 2500, 5000, 10000)) +
  theme_minimal(base_size = 24) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none")
Annual_ingestion_PopHist_Antarctic_Fe

dev.copy2pdf(file="Annual_ingestion_PopHist_Antarctic_Fe.pdf", width=14, height=9)


# Southern Hemisphere Fe, COMBINED
Annual_ingestion_PopComb_Antarctic_Fe <- Annual_filtfeed_Ant_projection_Nutrients %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. bonaerensis", "M. novaeangliae", "B. physalus")) %>% 
  ggplot() +
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_hist_Fe, ymax = IQR75_ingest_low_t_hist_Fe, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_ribbon(aes(ymin = IQR25_ingest_high_t_hist_Fe, ymax = IQR75_ingest_high_t_hist_Fe, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_curr_Fe, ymax = IQR75_ingest_low_t_curr_Fe, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_ribbon(aes(ymin = IQR25_ingest_high_t_curr_Fe, ymax = IQR75_ingest_high_t_curr_Fe, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_line(aes(days_feeding, med_ingest_low_t_curr_Fe, color = Species), 
            size = 1.5) +
  geom_line(aes(days_feeding, med_ingest_high_t_curr_Fe, color = Species), 
            size = 1.5, linetype = "dashed") +
  
  
  
  geom_line(aes(days_feeding, med_ingest_low_t_hist_Fe, color = Species), 
            size = 1.5) +
  geom_line(aes(days_feeding, med_ingest_high_t_hist_Fe, color = Species), 
            size = 1.5, linetype = "dashed") +
  scale_colour_manual(values = pal) +
  facet_grid(~Species) +   # Can add region in here if desired
  labs(x = "Days feeding",
       y = bquote('Projected Fe recycled'~(tonnes~yr^-1))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_log10(labels = scales::comma, limits = c(2, 12000),
                breaks = c(1, 10, 100, 500, 1000, 2500, 5000, 10000)) +
  theme_minimal(base_size = 24) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none")
Annual_ingestion_PopComb_Antarctic_Fe

dev.copy2pdf(file="Annual_ingestion_PopComb_Antarctic_Fe.pdf", width=14, height=8)




# Southern Hemisphere C net, COMBINED

Annual_ingestion_PopComb_Antarctic_C_export <- Annual_filtfeed_Ant_projection_Nutrients %>% 
  mutate(Species = fct_relevel(factor(abbr_binom(Species)), "B. bonaerensis", "M. novaeangliae", "B. physalus")) %>% 
  ggplot() +
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_curr_Fe_C_produced_Mt_C_exported_Mt, ymax = IQR75_ingest_low_t_curr_Fe_C_produced_Mt_C_exported_Mt, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_ribbon(aes(ymin = IQR25_ingest_high_t_curr_Fe_C_produced_Mt_C_exported_Mt, ymax = IQR75_ingest_high_t_curr_Fe_C_produced_Mt_C_exported_Mt, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  
  geom_ribbon(aes(ymin = IQR25_ingest_low_t_hist_Fe_C_produced_Mt_C_exported_Mt, ymax = IQR75_ingest_low_t_hist_Fe_C_produced_Mt_C_exported_Mt, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_ribbon(aes(ymin = IQR25_ingest_high_t_hist_Fe_C_produced_Mt_C_exported_Mt, ymax = IQR75_ingest_high_t_hist_Fe_C_produced_Mt_C_exported_Mt, 
                  x=days_feeding), 
              fill = "grey70", alpha = 0.5) +
  geom_line(aes(days_feeding, med_ingest_low_t_curr_Fe_C_produced_Mt_C_exported_Mt, color = Species), 
            size = 1.5) +
  geom_line(aes(days_feeding, med_ingest_high_t_curr_Fe_C_produced_Mt_C_exported_Mt, color = Species), 
            size = 1.5, linetype = "dashed") +
  
  
  
  geom_line(aes(days_feeding, med_ingest_low_t_hist_Fe_C_produced_Mt_C_exported_Mt, color = Species), 
            size = 1.5) +
  geom_line(aes(days_feeding, med_ingest_high_t_hist_Fe_C_produced_Mt_C_exported_Mt, color = Species), 
            size = 1.5, linetype = "dashed") +
  scale_colour_manual(values = pal) +
  facet_grid(~Species) +   # Can add region in here if desired
  labs(x = "Days feeding",
       y = bquote('Projected C exported'~(Mt~yr^-1))) + 
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_log10(labels = function(x) sprintf("%g", x),
                limits = c(0.01, 120),
                breaks = c(0.01, 0.1, 1, 10, 100)) +
  theme_minimal(base_size = 24) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none")
Annual_ingestion_PopComb_Antarctic_C_export

dev.copy2pdf(file="Annual_ingestion_PopComb_Antarctic_C_export.pdf", width=14, height=8)


# Old code below here

Annual_filtfeed_Ant_projection_Nutrients <- Annual_filtfeed_Ant_projection %>% 
  mutate_at(vars(c("med_ingest_low_t_curr":"IQR75_ingest_high_t_curr",
                   "med_ingest_low_t_hist":"IQR75_ingest_high_t_hist")), .funs = list(Fe = ~Fe_total_dw(.)/1000)) %>% # need to divide by 1000 to get the total iron recycled in tonnes
  mutate_at(vars(c("med_ingest_low_t_curr_Fe":"IQR75_ingest_high_t_curr_Fe",
                   "med_ingest_low_t_hist_Fe":"IQR75_ingest_high_t_hist_Fe")), .funs = list(C_produced_Mt = ~((.*0.75)/1e6)*11111.11)) %>% # Carbon production stimulated by Fe defecation in Mt C
  mutate_at(vars(c("med_ingest_low_t_curr":"IQR75_ingest_high_t_curr",
                   "med_ingest_low_t_hist":"IQR75_ingest_high_t_hist")), .funs = list(C_respired_Mt = ~(.*(0.1*0.75))/1e6)) %>%   # CHECK THIS, SHOULD BE 0.1, not 0.45? Carbon respired by populations in Mt C
  mutate_at(vars(c("med_ingest_low_t_curr_Fe_C_produced_Mt":"IQR75_ingest_high_t_curr_Fe_C_produced_Mt",
                   "med_ingest_low_t_hist_Fe_C_produced_Mt":"IQR75_ingest_high_t_hist_Fe_C_produced_Mt")), 
            .funs = list(C_exported_Mt = ~.*0.92))   #This is in fact C export; Check Lavery et al. 2010 ref, pg 3


# C_respired_med_curr_low = median(med_ingest_low_t_curr_C_respired_Mt, na.rm = TRUE),
# C_respired_IQR25_curr_low = median(IQR25_ingest_low_t_curr_C_respired_Mt, na.rm = TRUE),
# C_respired_IQR75_curr_low = median(IQR75_ingest_low_t_curr_C_respired_Mt, na.rm = TRUE),
# 
# C_respired_med_curr_high = median(med_ingest_high_t_curr_C_respired_Mt, na.rm = TRUE),
# C_respired_IQR25_curr_high = median(IQR25_ingest_high_t_curr_C_respired_Mt, na.rm = TRUE),
# C_respired_IQR75_curr_high = median(IQR75_ingest_high_t_curr_C_respired_Mt, na.rm = TRUE),     
# 
# C_respired_med_hist_low = median(med_ingest_low_t_hist_C_respired_Mt, na.rm = TRUE),
# C_respired_IQR25_hist_low = median(IQR25_ingest_low_t_hist_C_respired_Mt, na.rm = TRUE),
# C_respired_IQR75_hist_low = median(IQR75_ingest_low_t_hist_C_respired_Mt, na.rm = TRUE), 
# 
# C_respired_med_hist_high = median(med_ingest_high_t_hist_C_respired_Mt, na.rm = TRUE),
# C_respired_IQR25_hist_high = median(IQR25_ingest_high_t_hist_C_respired_Mt, na.rm = TRUE),
# C_respired_IQR75_hist_high = median(IQR75_ingest_high_t_hist_C_respired_Mt, na.rm = TRUE),                            


