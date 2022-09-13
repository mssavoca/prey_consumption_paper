####################################################################
# Southern Ocean whale prey consumption analysis for Anaëlle Durfort 
####################################################################

#Scaling relationships

engulfallo_FR_biomass <- tribble(
  ~Species, ~slope,   ~intercept,   ~per25th_FR,   ~med_FR,   ~per75th_FR,    ~Biomass_perm3, 
  "Minke",    3.11151,  -2.31328,    849,           1207,     1615,           0.15,   
  "Fin",      3.54883,  -2.85446,    159,           222,      302,            0.17, 
  "Blue",     3.667115, -3.024580,   138,           204,      270,            0.17,
  "Humpback", 3.24603,  -2.15150,    402,           622,      914,            0.19
)


# Engulfment Capacity data from Kahane-Rapport and Goldbogen 2018. J Morph
# Feeding rate data from Savoca et al. 2021 Nature
# Krill biomass per m3 from Savoca et al. 2021 Nature and Goldbogen et al. 2019 Science



AD_whale_f <- read_xlsx("whale_size_at_age_durfort.xlsx", sheet = 2) %>% 
  rename("Age" = "...1") %>% 
  pivot_longer(cols = Blue:Minke, names_to = "Species", values_to = "length_m") %>% 
  left_join(engulfallo_FR_biomass, by = "Species") %>% 
  mutate(Engulfment_m3 = length_m ^ slope * 10 ^ intercept,
         DailyPreyCons_kg_25th = per25th_FR * Engulfment_m3 * Biomass_perm3,
         DailyPreyCons_kg_Med = med_FR * Engulfment_m3 * Biomass_perm3,
         DailyPreyCons_kg_75th = per75th_FR * Engulfment_m3 * Biomass_perm3,
         )

write_csv(AD_whale_f, "AD_whale_f.csv")


AD_whale_m <- read_xlsx("whale_size_at_age_durfort.xlsx", sheet = 3) %>% 
  rename("Age" = "...1") %>% 
  pivot_longer(cols = Blue:Minke, names_to = "Species", values_to = "length_m") %>% 
  left_join(engulfallo_FR_biomass, by = "Species") %>% 
  mutate(Engulfment_m3 = length_m ^ slope * 10 ^ intercept,
         DailyPreyCons_kg_25th = per25th_FR * Engulfment_m3 * Biomass_perm3,
         DailyPreyCons_kg_Med = med_FR * Engulfment_m3 * Biomass_perm3,
         DailyPreyCons_kg_75th = per75th_FR * Engulfment_m3 * Biomass_perm3,
  )


write_csv(AD_whale_m, "AD_whale_m.csv")



