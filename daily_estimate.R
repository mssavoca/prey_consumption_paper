library(tidyverse)
library(maptools)
library(readxl)


# stuff from Max to generate data ----

load("rates.RData")
load("tag_guide4.RData")
load("rates.RData")
load("rates[1].RData")

# bringing in whale length data from Paolo; using allometric equations from Shirel to convert to engulfment capacity
whale_lengths <- read.csv("whale_masses.csv") %>% 
  rename(SpeciesCode = "species") %>% 
  left_join(engulf_allo, by = "SpeciesCode") %>% 
  mutate(Engulfment_m3 = length ^ slope * 10 ^ intercept)

# Add in Bryde's data, MAX NEEDS TO CHECK THIS
# Combining day and twilight lunge rates
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


# Dave's KRILL data (from 11/21/19), natural log transformed, which is what we were donig in the prior iteration

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
     xlab = "dur",
     ylab = "weight")
 
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




krill_biomass_estimates <- krill_rate_estimates %>% 
  rename(SpeciesCode = "species_code") %>% 
  group_by(SpeciesCode) %>% 
  mutate(Engulf_cap_m3 = sample(whale_lengths$Engulfment_m3[whale_lengths$SpeciesCode == SpeciesCode[1]],
                                size = n(),
                                replace = TRUE),
         Mass_est_kg = sample(whale_lengths$mass[whale_lengths$SpeciesCode == SpeciesCode[1]],
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
#        25%      50%      75%      95% 
#   148.2142 206.9413 289.0931 454.3588 
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
                                                    region == "Temperate" ~ daily_consumption_high_kg*3628),
         Total_energy_intake_best_high_kJ = case_when(region == "Polar" ~ daily_consumption_low_kg*4575,
                                                     region == "Temperate" ~ daily_consumption_high_kg*3628),
         Mass_specifc_energy_intake_best_high_kJ = Total_energy_intake_best_high_kJ/Mass_est_kg,
         Mass_specifc_energy_intake_best_low_kJ = Total_energy_intake_best_low_kJ/Mass_est_kg,
  )


krill_daily %>%  
  group_by(SpeciesCode) %>% 
  filter(SpeciesCode == "bw") %>%
  pull(daily_consumption_low_kg) %>%
  summary()


  ggplot(aes(daily_biomass_kg_best_low)) + 
  geom_density(aes(fill = SpeciesCode), alpha = 0.2) +
  facet_wrap(~SpeciesCode, scales = "free_x") +
  xlim(0,20000)
p

save(daynight_rates, krill_daily, file = "data/outputs/daynights.RData")




# Fish feeding information

fish_rate_estimates <- estimate_rate("Fish", 1, 1)


fish_biomass_estimates <- fish_rate_estimates %>% 
  rename(SpeciesCode = "species_code") %>% 
  group_by(SpeciesCode) %>% 
  mutate(Engulf_cap_m3 = sample(whale_lengths$Engulfment_m3[whale_lengths$SpeciesCode == SpeciesCode[1]],
                                size = n(),
                                replace = TRUE),
         Mass_est_kg = sample(whale_lengths$mass[whale_lengths$SpeciesCode == SpeciesCode[1]],
                              size = n(),
                              replace = TRUE))
  




latitudes <- tribble(
  ~ region,    ~ latitude,
  "Polar",     65,
  "Temperate", 36
)

# Overall daily rate, Modified by MA
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
      daily_mean_low_kg = mean((rnorm(floor(daily_rate), 0.375, 0.243))*0.29),
      daily_mean_high_kg = mean(rnorm(floor(daily_rate), 0.375, 0.243)),
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
         Mass_specifc_energy_intake_best_high_kJ = Total_energy_intake_best_high_kJ/Mass_est_kg,
         Mass_specifc_energy_intake_best_low_kJ = Total_energy_intake_best_low_kJ/Mass_est_kg,
  )


fish_daily %>%  
  group_by(SpeciesCode) %>% 
  # filter(SpeciesCode == "mn") %>%
  # pull(Mass_specifc_energy_intake_best_high_kJ) %>%
  # summary()
  
ggplot(aes(daily_consumption_high_kg)) + 
  geom_density(aes(fill = SpeciesCode), alpha = 0.2) +
  facet_wrap(~SpeciesCode, scales = "free_x") +
  xlim(0,5000)
p


# generating plots and tables for paper ----

#summary tables of lunges, water filtered, and prey consumed per day (Figure 2)----
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
    med_daily_En_best = round(median(Total_energy_intake_best_low_kJ/1e6), 2),
    med_daily_En_best_IQR25 = round(quantile(Total_energy_intake_best_low_kJ/1e6, probs = 0.25, na.rm = TRUE), 2), 
    med_daily_En_best_IQR75 = round(quantile(Total_energy_intake_best_low_kJ/1e6, probs = 0.75, na.rm = TRUE), 2), 
    med_daily_mass_specific_En_best = round(median(Mass_specifc_energy_intake_best_low_kJ), 2),
    med_daily_En_top50 = round(median(Total_energy_intake_best_high_kJ/1e6), 2),
    med_daily_En_top50_IQR25 = round(quantile(Total_energy_intake_best_high_kJ/1e6, probs = 0.25, na.rm = TRUE), 2), 
    med_daily_En_top50_IQR75 = round(quantile(Total_energy_intake_best_high_kJ/1e6, probs = 0.75, na.rm = TRUE), 2),
    med_daily_mass_specific_En_high = round(median(Mass_specifc_energy_intake_best_high_kJ), 2)) %>% 
  #unite("Daily energy intake (GJ) lower estimate IQR", c(med_daily_En_hyp_low_IQR25, med_daily_En_hyp_low_IQR75), sep = "-") %>% 
  unite("Daily energy intake (GJ) IQR", c(med_daily_En_best_IQR25, med_daily_En_best_IQR75), sep = "-") %>% 
  unite("Daily energy intake (GJ) Top 50% IQR", c(med_daily_En_top50_IQR25, med_daily_En_top50_IQR75), sep = "-")





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
         filtration182.5hist = Engulf_cap_m3*daily_rate*182.5*`Historical estimate`
    )


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

# And now for the prey

# Southern Hemisphere current----



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
    mutate(biomass_best_low_m3 = rlnorm(n(), Biomass[1], `Biomass sd`[1]),
           biomass_best_high_m3 = rlnorm(n(), BiomassTop50[1], BiomassTop50sd[1]),
           daily_rate = day_rate * daylen + night_rate * nightlen) %>% 
    ungroup() %>% 
    mutate(daily_biomass_kg_best_low = biomass_best_low_m3*Engulf_cap_m3*daily_rate,
           daily_biomass_kg_best_high = biomass_best_high_m3*Engulf_cap_m3*daily_rate)
  
  
  
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
  mutate(Total_energy_intake_best_low_kJ = case_when(region == "Polar" ~ daily_biomass_kg_best_low*4575,
                                                     region == "Temperate" ~ daily_biomass_kg_best_low*3628),
         Total_energy_intake_best_high_kJ = case_when(region == "Polar" ~ daily_biomass_kg_best_high*4575,
                                                      region == "Temperate" ~ daily_biomass_kg_best_high*3628),
         Mass_specifc_energy_intake_best_high_kJ = Total_energy_intake_best_high_kJ/Mass_est_kg,
         Mass_specifc_energy_intake_best_low_kJ = Total_energy_intake_best_low_kJ/Mass_est_kg,
  )


krill_daily_Ant_projection %>%  
  group_by(SpeciesCode) %>% 
  filter(SpeciesCode == "bw") %>%
  pull(daily_biomass_kg_best_high) %>%
  summary()
  






