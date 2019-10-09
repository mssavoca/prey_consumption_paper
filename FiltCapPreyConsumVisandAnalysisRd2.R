## Analyses and visualizations for Filtration capacity and prey consumption paper ----


# Load functions and packages ----

library(ggplot2)
library(dplyr)
library(data.table)
library(readxl)
library(forcats)
library(tidyverse)
library(mgcv)
library(lme4)
library(lmerTest)
library(MuMIn)
library(ggpubr)


# formula for standard error
SE = function(x){sd(x)/sqrt(sum(!is.na(x)))}

# Abbreviate a binomial e.g. Balaenoptera musculus -> B. musculus
abbr_binom = function(binom) {
  paste(str_sub(binom, 1, 1), 
        str_extract(binom, " .*"), 
        sep = ".")
}


# Data importing, cleaning, combining ----

## SpeciesCode will be "key" variable


# Species-specific average length data (from Shirel's paper), to recalculate average engulfment capacity
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

# Allometric equations from Shirel's paper
# creating fucntions from Shirel's paper for MW (in kg) for engulfment capacity in liters for each species where we have a known length
engulf_allo <- tribble(
  ~SpeciesCode, ~slope,   ~intercept,
  "bb",     3.10910,  0.69969,
  "be",     3.1453,   0.5787,
  "bp",     3.54883,  0.15604,
  "bw",     3.667316, -0.014078,
  "mn",     3.24603,  0.85934
)


# Skips first two rows
tag_guide <- read_excel("TAG GUIDE_2.18.19.xlsx", skip = 2) %>% 
  rename(Study_Area = `Study_Area     _`,
         SpeciesCode = `Spec      _`)

# read in whale population data. Current data from IUCN 2019 Redlist, whaling data compiled in Rocha et al. 2014
pop_data <- read_excel("Filtration project_Whale population and feeding information.xlsx", sheet = 1)


# Load prey distributions from Jeremy's paper
d_full_NULL <- read_csv("Cetacea model output NULL_EXTANT.csv") %>% 
  mutate(Species = ifelse(Species == "bonarensis", "bonaerensis", Species),
         SpeciesCode = case_when(Species ==  "musculus" ~ "bw",
                                 Species == "physalus"~ "bp",
                                 Species ==  "novaeangliae" ~ "mn",
                                 Species == "bonaerensis" ~ "bb"))
d_full_BOUT <- read_csv("Cetacea model output BOUT_EXTANT.csv") %>% 
  mutate(Species = ifelse(Species == "bonarensis", "bonaerensis", Species),
         SpeciesCode = case_when(Species ==  "musculus" ~ "bw",
                                 Species == "physalus"~ "bp",
                                 Species ==  "novaeangliae" ~ "mn",
                                 Species == "bonaerensis" ~ "bb"))

# summarize prey distributions from Jeremy's paper
d_sum_NULL <- d_full_NULL %>%
  group_by(SpeciesCode) %>%
  dplyr::summarize(wgtMeanNULL_wt_g = weighted.mean(`Prey W (g)`, Percent),
                   medNULL_wt_g = median(`Prey W (g)`),
                   wgtMeanNULL_E = weighted.mean(`Energy (kJ)`, Percent),
                   medNULL_E = median(`Energy (kJ)`))
d_sum_NULL <- na.omit(d_sum_NULL)

d_sum_BOUT = d_full_BOUT %>% 
  group_by(SpeciesCode) %>% 
  dplyr::summarize(wgtMeanBOUT_wt_g = weighted.mean(`Prey W (g)`, Percent), 
                   medBOUT_wt_g = median(`Prey W (g)`), 
                   wgtMeanBOUT_E = weighted.mean(`Energy (kJ)`, Percent),
                   medBOUT_E = median(`Energy (kJ)`))
d_sum_NULL <- na.omit(d_sum_NULL)



# Rorqual feeding data from Jeremy's scaling paper
RorqualData <- read_csv("lunge_rates_from_Paolo.csv") %>% 
  mutate(`deployment-time_h` = `deployment-time_secs`/60/60) %>% 
  select(-c(species)) 

# From Danuta, some of the deployments differ (for example see 2014 blue whales)
OdontoceteData <- read_csv("foragestats_combined_ko2.csv") %>%
  filter(Species %in% c("Balaenoptera_physalus", "Balaenoptera_musculus", "Balaenoptera_bonaerensis", "Megaptera_novaeangliae"))

cetacean_data <- left_join(OdontoceteData, RorqualData, by = c("ID")) %>%
  # feeding events are lunges, buzzes for rorquals, odontocetes
  mutate(TotalFeedingEvents = coalesce(total_lunges, total_buzz_count),
         TotalTagTime_h = coalesce(`deployment-time_h`, total_duration_h))   # fixing the ba/bb confusion


#joining dataframes into the master to be used for all analyses
master_filtprey_data <- read_excel("ALLPRHS 2.5.2019.xls") %>% 
  full_join(select(cetacean_data, c(ID, lunge_quality, sonar_exp, TotalFeedingEvents, TotalTagTime_h)), by = "ID") %>%   # WHICH JOIN SHOULD THIS BE?!?!?!?
  left_join(select(tag_guide, ID, Study_Area), by = "ID") %>% 
  mutate(         
    # fixing the ba/bb confusion, BUT THE LOCATION IS SOCAL?? SHOULD I JUST OMIT THESE WEIRD ONES???
    SpeciesCode = substr(ID,1,2),
    SpeciesCode = if_else(SpeciesCode == "ba", "bb", SpeciesCode),
    # Replace NAs in location with "SoCal"
    Study_Area = replace_na(Study_Area, "SoCal"),
    # Make whale length a number
    whaleLength = parse_number(whaleLength),
    CommonName = case_when(
      SpeciesCode == "bw" ~ "Blue Whale",
      SpeciesCode == "bp" ~ "Fin Whale",
      SpeciesCode == "mn" ~ "Humpback Whale",
      SpeciesCode %in% c("ba", "bb") ~ "Minke Whale", 
      TRUE ~ "Bryde's Whale"),
    Species = case_when(
      SpeciesCode == "bw" ~ "Balaenoptera musculus",
      SpeciesCode == "bp" ~ "Balaenoptera physalus",
      SpeciesCode == "mn" ~ "Megaptera novaeangliae",
      SpeciesCode %in% c("bb","ba") ~ "Balaenoptera bonaerensis", 
      SpeciesCode == "be" ~ "Balaenoptera edeni")) %>% 
  select(-notes) %>% 
  # Join population numbers
  left_join(select(pop_data, -c(SpeciesCode, `Current Range`)), by = "Species") %>% 
  # Join species-median size data
  left_join(select(v_data, -CommonName), by = "SpeciesCode") %>% 
  # Join allometric equations
  left_join(engulf_allo, by = "SpeciesCode") %>% 
  # Calculate engulfment volumes
  mutate(
    BestLengthEst = coalesce(whaleLength, med_TLm),
    Engulfment_L = BestLengthEst ^ slope * 10 ^ intercept * 0.9766,
    # Feeding rates
    TotalLunges = dayhourslunges + nighthourslunges + twilightzonelunges,
    TotalHours = dayhours + nighthours + twilightzone,
    TotalTagHours = coalesce(TotalHours, TotalTagTime_h),
    TotalLungeEvents = coalesce(TotalLunges, as.double(TotalFeedingEvents)),
    LungesPerHour = TotalLungeEvents/TotalTagHours,
    # These will definitely change once Shirel, Dave, and I get updated daytime numbers
    LungesPerDayHour = dayhourslunges/dayhours,
    LungesPerNightHour = nighthourslunges/nighthours,
    LungesPerTwHour = twilightzonelunges/twilightzone,
    EngulfVolPerHr = Engulfment_L*LungesPerHour) %>%
  # Replace NaNs with NAs in LungesPer* columns
  mutate_at(c("LungesPerDayHour", "LungesPerNightHour", "LungesPerTwHour"), 
            function(col) if_else(!is.nan(col), col, NA_real_)) %>% 
  # remove old columns that were coalesced
  select(-c(TotalHours, TotalFeedingEvents, TotalLunges)) %>% 
  # Prey data
  separate(prey, into = c("Prey", "Prey notes"), sep = " ") %>% 
  filter(!Prey %in% c("Milk", "N")) %>% 
  mutate(
    PreyClean = case_when(
      Species == "Balaenoptera bonaerensis" ~ "Krill-feeding",
      SpeciesCode == "bw" ~ "Krill-feeding",     # ASK DAVE ABOUT FIN AND HUMPBACK WHALE PREY
      SpeciesCode == "be" ~ "Fish-feeding",
      Study_Area == "Antarctic" ~ "Krill-feeding",
      Prey %in% c("Anchovies", "Fish", "Herring", "SandLance", "Sardines")  ~ "Fish-feeding",
      Prey %in% c("Inverts", "Krill") ~ "Krill-feeding"), 
    Region = case_when(
      Study_Area %in% c("Monterey", "SoCal", "Cordell Bank", "San Diego", "WA Coast")  ~ "Eastern North Pacific",
      Study_Area %in% c("Stellwagen", "Norway", "Azores", "Greenland") ~ "North Atlantic",
      Study_Area == "South Africa" ~ "South Africa",
      Study_Area == "Antarctic" ~ "Antarctic",
      Study_Area == "Chile" ~ "Chile"))



write_csv(master_filtprey_data, "Deployments and data for filtprey analyses.csv")


# preliminary plots for filtration capacity ----

pal <- c("B. bonaerensis" = "firebrick3", "B. edeni" = "darkorchid3",  "M. novaeangliae" = "gray30", "B. physalus" = "chocolate3", "B. musculus" = "dodgerblue2")
Shape <- c("B. bonaerensis" = 15, "B. edeni" = 8, "M. novaeangliae" = 17, "B. physalus" = 18, "B. musculus" = 19)

EmpEngulfCapPLot <- ggplot(master_filtprey_data, aes(whaleLength, Engulfment_L)) +
  geom_point(inherit.aes = TRUE, aes(shape = abbr_binom(Species), color = abbr_binom(Species)), size = 2) + 
  #labs(col="Species", shape = "Species") +
  geom_smooth(method = lm, color = "black", size  = 0.5) +
  xlab("Length (m)") + ylab("Engulfment capacity (L)") +
  ggtitle("Engulfment capacity by length for whales measured by drone") +
  scale_colour_manual(values = pal) +
  scale_shape_manual(values= Shape) +
  theme_bw() +
  theme(legend.text = element_text(face="italic"))
EmpEngulfCapPLot + guides(color = guide_legend(title = "Species"), 
                          shape = guide_legend(title = "Species"),
                          face = "italic") +
  # the specifics of the linear equation will likely change and need to be updated
  annotate("text", x = 10, y = 110000,
           label = c("y = 7516.47x - 76066.55")) +
  annotate("text", x = 10, y = 100000, 
           label = "italic(R)^2 == 0.88",
           parse = TRUE) +
  annotate("text", x = 10, y = 90000, 
           label = "italic(p) <2e-16",
           parse = TRUE)
  
EmpEngulfCapMod <- lmer(Engulfment_L ~ whaleLength + (1|Species), data = master_filtprey_data)
summary(EmpEngulfCapMod)
r.squaredGLMM(EmpEngulfCapMod)  # intereestin in R2m, the R2 on the marginal (i.e., fixed) effects, see: https://ecologyforacrowdedplanet.wordpress.com/2013/08/27/r-squared-in-mixed-models-the-easy-way/


master_filtprey_data$Species <- fct_relevel(master_filtprey_data$Species, 
                                            "Balaenoptera edeni", "Balaenoptera bonaerensis","Megaptera novaeangliae","Balaenoptera physalus","Balaenoptera musculus")

v_all_deploy <- ggplot(filter(master_filtprey_data, TotalTagHours > 2 & sonar_exp %in% c("none", NA) & PreyClean %in% c("Krill-feeding", "Fish-feeding")),
                       aes(x = fct_reorder(abbr_binom(Species), EngulfVolPerHr, fun = median), y = EngulfVolPerHr, 
                           color = abbr_binom(Species), shape = abbr_binom(Species))) + 
  geom_point(inherit.aes = T, alpha = 0.8, position = position_jitter(width = .25)) + 
  geom_boxplot(inherit.aes = T, guides = FALSE, outlier.shape = NA, alpha = 0) +
  facet_grid(.~PreyClean, scales = "free_x") +
  #scale_y_log10(labels = scales::comma) +
  scale_colour_manual(values = pal) +
  scale_shape_manual(values = Shape) +
  xlab("Species") + ylab("Filtration capacity (liters per hour)") + 
  ggtitle("Water filtered per individual per hour (all delpoyments >2 hours)") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1, face="italic"),
        axis.text.y = element_text(size=12),
        axis.title=element_text(size=13, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        strip.text.x = element_text(size = 12))
v_all_deploy + theme(legend.position="none") 

# summary results table for water filtration per individual per hour
vol_summ <- master_filtprey_data %>% 
  group_by(SpeciesCode) %>% 
  filter(TotalHours > 2 & TotalLunges > 0 & sonar_exp =="none", PreyClean == "Krill-feeding") %>% 
  dplyr::summarize(MeanWaterFiltered = mean(EngulfVolPerHr),
            SEWaterFiltered = SE(EngulfVolPerHr))



v_24_deploy <- ggplot(filter(master_filtprey_data, TotalTagHours > 24 & sonar_exp %in% c("none", NA)),
                      aes(x = reorder(abbr_binom(Species), EngulfVolPerHr, FUN = median), y = EngulfVolPerHr*24, 
                                      color = abbr_binom(Species), shape = abbr_binom(Species))) + 
  geom_point(inherit.aes = T, alpha = 0.8, position = position_jitter(width = .25)) + 
  geom_boxplot(inherit.aes = T, guides = FALSE, outlier.shape = NA, alpha = 0) +
  facet_grid(.~PreyClean, scales = "free_x") +
  scale_colour_manual(values = pal) +
  scale_y_log10(labels = scales::comma) + 
  scale_shape_manual(values = Shape) +
  xlab("Species") + ylab("Total water filtered (liters per day)") + 
  ggtitle("Water filtered per individual per day (tags on >24 hours)") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1, face="italic"),
        axis.text.y = element_text(size=12),       
        axis.title=element_text(size=13, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        strip.text.x = element_text(size = 12))
v_24_deploy + theme(legend.position="none") 


# now for varying hours feeding within a day

vol_master_for_join <- master_filtprey_data %>% 
  mutate(dummy = 1)
vol_master_varying_HrperD <- tibble(hours_feeding = seq(1,12,1), dummy = 1) %>% 
  full_join(vol_master_for_join, by = "dummy") %>% 
  select(-dummy) %>% 
  mutate(TotalWaterFiltered = hours_feeding*EngulfVolPerHr)


v_HrperDfish <- ggplot(filter(vol_master_varying_HrperD, TotalTagHours > 2 & sonar_exp %in% c("none", NA) & PreyClean == "Fish-feeding"),
                       aes(x = hours_feeding, y = TotalWaterFiltered, color = abbr_binom(Species), shape = abbr_binom(Species))) + 
  geom_point(inherit.aes = T, aes(group = abbr_binom(Species)), alpha = 0.3) + 
  geom_smooth(inherit.aes = T, aes(group = abbr_binom(Species)), color = "black", size = 0.5, method = "lm") +
  facet_grid(.~abbr_binom(Species)) + #HERE
  scale_colour_manual(values = pal) +
  scale_shape_manual(values = Shape) +
  scale_x_continuous(breaks=seq(0, 12, 3)) +
  #scale_y_log10(labels = scales::comma) + 
  ylim(0, 5000000) +
  xlab("Hours feeding") + ylab("Total water filtered (L)") + 
  ggtitle("Water filtered per individual per day (all deployments >2 hours)") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        strip.text.x = element_text(size = 12, face = "italic"))
v_HrperDfish + theme(legend.position="none")

vol_master_varying_HrperD$Species <- fct_relevel(vol_master_varying_HrperD$Species, 
                                            "B. edeni", "B. bonaerensis","M. novaeangliae","B. physalus","B. musculus")


v_HrperDkrill <- ggplot(filter(vol_master_varying_HrperD, TotalTagHours > 2 & sonar_exp %in% c("none", NA) & PreyClean == "Krill-feeding"),
                        aes(x = hours_feeding, y = TotalWaterFiltered, color = abbr_binom(Species), shape = abbr_binom(Species))) + 
  geom_point(inherit.aes = T, aes(group = abbr_binom(Species), alpha = 0.3)) + 
  geom_smooth(inherit.aes = T, aes(group = abbr_binom(Species)), color = "black", size = 0.5, method = "lm") +
  facet_grid(.~abbr_binom(Species)) +
  scale_colour_manual(values = pal) +
  scale_shape_manual(values = Shape) +
  scale_x_continuous(breaks=seq(0, 12, 3)) +
  #scale_y_log10(labels = scales::comma) + 
  ylim(0, 25000000) +
  xlab("Hours feeding") + ylab("Total water filtered (L)") + 
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        strip.text.x = element_text(size = 12, face = "italic"))
v_HrperDkrill + theme(legend.position="none")

ggarrange(v_HrperDfish, v_HrperDkrill, 
          labels = c("A", "B"), # THIS IS SO COOL!!
          legend = "none",
          ncol = 1, nrow = 2)


# annual amount of water filtered per individual, varying days
vol_master_for_year_join <- vol_master_varying_HrperD %>% 
  mutate(dummy = 1)
vol_master_varying_DperYr <- tibble(days_feeding = seq(60,182.5,10), dummy = 1) %>% 
  full_join(vol_master_for_year_join, by = "dummy") %>% 
  select(-dummy) %>% 
  mutate(TotalAnnualWaterFiltered = days_feeding*TotalWaterFiltered)


v_DperYrfish <- ggplot(filter(vol_master_varying_DperYr, hours_feeding %in% (6:12) & TotalTagHours > 2 & sonar_exp %in% c("none", NA) & PreyClean == "Fish-feeding"),
                       aes(x = days_feeding, y = TotalAnnualWaterFiltered, color = abbr_binom(Species), shape = abbr_binom(Species))) +
  geom_point(inherit.aes = T, aes(group = abbr_binom(Species)), alpha = 0.2) +
  geom_smooth(inherit.aes = T, aes(group = abbr_binom(Species)), color = "black", size = 0.5, method = "lm") +
  facet_grid(.~abbr_binom(Species)) +
  scale_colour_manual(values = pal) +
  scale_shape_manual(values = Shape) +
  #scale_y_log10(labels = scales::comma) + 
  #ylim(0, 25000000) +
  xlab("Days feeding") + ylab("Total water filtered (L)") +
  ggtitle("Water filtered per individual per year (all deployments >2 hours)") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        strip.text.x = element_text(size = 12, face = "italic"))
v_DperYrfish + theme(legend.position="none")

vol_master_varying_DperYr$SpeciesCode <- fct_relevel(vol_master_varying_DperYr$SpeciesCode, "be", "bb","mn","bp","bw")
v_DperYrkrill <- ggplot(filter(vol_master_varying_DperYr, hours_feeding %in% (6:12) & TotalTagHours > 2 & sonar_exp %in% c("none", NA) & PreyClean == "Krill-feeding"),
                        aes(x = days_feeding, y = TotalAnnualWaterFiltered, color = abbr_binom(Species), shape = abbr_binom(Species))) +
  geom_point(inherit.aes = T, aes(group = abbr_binom(Species)), alpha = 0.2) +
  geom_smooth(inherit.aes = T, aes(group = abbr_binom(Species)), color = "black", size = 0.5, method = "lm") +
  facet_grid(.~abbr_binom(Species)) +
  scale_colour_manual(values = pal) +
  scale_shape_manual(values = Shape) +
  #scale_y_log10(labels = scales::comma) + 
  #ylim(0, 25000000) +
  xlab("Days feeding") + ylab("Total water filtered (L)") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        strip.text.x = element_text(size = 12, face = "italic"))
v_DperYrkrill + theme(legend.position="none")

ggarrange(v_DperYrfish, v_DperYrkrill, 
          labels = c("A", "B"), # THIS IS SO COOL!!
          legend = "none",
          ncol = 1, nrow = 2)

# annual amount of water filtered current global populartion, varying days
v_DperYrfish_currentpop <- ggplot(filter(vol_master_varying_DperYr, hours_feeding %in% (6:12) & TotalTagHours > 2 & sonar_exp %in% c("none", NA) & PreyClean == "Fish-feeding"),
                                  aes(x = days_feeding, y = TotalAnnualWaterFiltered*`Population estimate`, color = abbr_binom(Species), shape = abbr_binom(Species))) +
  geom_point(inherit.aes = T, aes(group = abbr_binom(Species)), alpha = 0.2) +
  geom_smooth(inherit.aes = T, aes(group = abbr_binom(Species)), color = "black", size = 0.5, method = "lm") +
  facet_grid(.~abbr_binom(Species)) +
  scale_colour_manual(values = pal) +
  scale_shape_manual(values = Shape) +
  #scale_y_log10(labels = scales::comma) + 
  xlab("Days feeding") + ylab("Total water filtered (L)") +
  ggtitle("Water filtered annually by current global population (all deployments >2 hours)") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        strip.text.x = element_text(size = 12, face = "italic"))
v_DperYrfish_currentpop + theme(legend.position="none")


v_DperYrkrill_currentpop <- ggplot(filter(vol_master_varying_DperYr, hours_feeding %in% (6:12) & TotalTagHours > 2 & sonar_exp %in% c("none", NA) & PreyClean == "Krill-feeding"),
                                   aes(x = days_feeding, y = TotalAnnualWaterFiltered*`Population estimate`, color = abbr_binom(Species), shape = abbr_binom(Species))) +
  geom_point(inherit.aes = T, aes(group = abbr_binom(Species)), alpha = 0.2) +
  geom_smooth(inherit.aes = T, aes(group = abbr_binom(Species)), color = "black", size = 0.5, method = "lm") +
  facet_grid(.~abbr_binom(Species)) +
  scale_colour_manual(values = pal) +
  scale_shape_manual(values = Shape) +
  #scale_y_log10(labels = scales::comma) + 
  ylim(0, 4e14) +
  xlab("Days feeding") + ylab("Total water filtered (liters)") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        strip.text.x = element_text(size = 12, face = "italic"))
v_DperYrkrill_currentpop + theme(legend.position="none")


ggarrange(v_DperYrfish_currentpop, v_DperYrkrill_currentpop, 
          labels = c("A", "B"), # THIS IS SO COOL!!
          legend = "none",
          ncol = 1, nrow = 2)











# Excess code and plots below here ----

# annual amount of water filtered historical global population, varying days
v_DperYrfish_histpop <- ggplot(filter(vol_master_varying_DperYr, hours_feeding %in% (6:12) & TotalTagTime_h > 2 & TotalLunges > 0 & sonar_exp =="none", PreyClean == "Fish-feeding"),
                               aes(x = days_feeding, y = TotalAnnualWaterFiltered*`Total removed`, color = SpeciesCode, shape = SpeciesCode)) +
  geom_point(inherit.aes = T, aes(group = SpeciesCode), alpha = 0.2) +
  geom_smooth(inherit.aes = T, aes(group = SpeciesCode), color = "black", size = 0.5) +
  facet_grid(.~SpeciesCode, labeller = as_labeller(SpCodetoName)) +
  scale_colour_manual(values = pal,
                      labels = c("B. acutorostrata", "B. bonaerensis", "B. edeni", "B. musculus", "B. physalus", "M. novaeangliae")) +
  scale_shape_manual(values = Shape,
                     labels = c("B. acutorostrata", "B. bonaerensis", "B. edeni", "B. musculus", "B. physalus", "M. novaeangliae")) +
  scale_y_log10(labels = scales::comma) + 
  xlab("Days feeding") + ylab("Total water filtered (liters)") +
  ggtitle("Annual filtration capacity removed by 20th century whaling (all deployments >2 hours)") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        strip.text.x = element_text(size = 12, face = "italic"))
v_DperYrfish_currentpop + theme(legend.position="none")

v_DperYrkrill_histpop <- ggplot(filter(vol_master_varying_DperYr, hours_feeding %in% (6:12) & TotalTagTime_h > 2 & TotalLunges > 0 & sonar_exp =="none", PreyClean == "Krill-feeding"),
                                aes(x = days_feeding, y = TotalAnnualWaterFiltered*`Total removed`, color = SpeciesCode, shape = SpeciesCode)) +
  geom_point(inherit.aes = T, aes(group = SpeciesCode), alpha = 0.2) +
  geom_smooth(inherit.aes = T, aes(group = SpeciesCode), color = "black", size = 0.5) +
  facet_grid(.~SpeciesCode, labeller = as_labeller(SpCodetoName)) +
  scale_colour_manual(values = pal,
                      labels = c("B. acutorostrata", "B. bonaerensis", "B. edeni", "B. musculus", "B. physalus", "M. novaeangliae")) +
  scale_shape_manual(values = Shape,
                     labels = c("B. acutorostrata", "B. bonaerensis", "B. edeni", "B. musculus", "B. physalus", "M. novaeangliae")) +
  scale_y_log10(labels = scales::comma) + 
  xlab("Days feeding") + ylab("Total water filtered (liters)") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        strip.text.x = element_text(size = 12, face = "italic"))
v_DperYrkrill_histpop + theme(legend.position="none")


ggarrange(v_DperYrfish_histpop, v_DperYrkrill_histpop, 
          labels = c("A", "B"), # THIS IS SO COOL!!
          legend = "none",
          ncol = 1, nrow = 2)