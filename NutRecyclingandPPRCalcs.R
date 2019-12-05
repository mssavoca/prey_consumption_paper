##############################################################
# Calculation of nutrient recycling potential of baleen whales 
##############################################################


#############################
# Load functions and packages
#############################
library(ggplot2)
library(dplyr)
library(data.table)
library(readxl)
library(forcats)
library(tidyverse)
library(mgcv)
library(lme4)
library(lmerTest)
library(ggpubr)


# formula for standard error
SE = function(x){sd(x)/sqrt(sum(!is.na(x)))}

# Abbreviate a binomial e.g. Balaenoptera musculus -> B. musculus
abbr_binom = function(binom) {
  paste(str_sub(binom, 1, 1), 
        str_extract(binom, " .*"), 
        sep = ".")
}








#trying some stuff----

Fe_calc_kg <- function(x) {
  x*0.8*0.000146
} # From Doughty et al. 2016, Roman and McCarthy 2010
N_calc_kg <- function(x) {
  x*0.8*0.0056
} # see Roman et al. 2016, converted moles PON to g to kg
P_calc_kg <- function(x) {
  x*0.8*0.0089
} # whale fecal average from Ratnarajah et al. 2014

d_strapped

# d_strapped_Ant_fecal_projection <- d_strapped_Ant_projection %>% 
#   mutate(
#     total_feces_best_low_kg = prey_mass_per_day_best_low_kg*0.8,    # From Doughty et al. 2016, Roman and McCarthy 2010
#             total_Fe_recycled_best_low_kg = total_feces_best_low_kg*0.000146,    # whale fecal average from Ratnarajah et al. 2014
#             total_N_recycled_best_low_kg = total_feces_best_low_kg*0.0056,     # see Roman et al. 2016, converted moles PON to g to kg  
#             total_P_recycled_best_low_kg = total_feces_best_low_kg*0.0089,
#     total_feces_best_upper_kg = prey_mass_per_day_best_upper_kg*0.8,    # From Doughty et al. 2016, Roman and McCarthy 2010
#             total_Fe_recycled_best_upper_kg = total_feces_best_upper_kg*0.000146,    # whale fecal average from Ratnarajah et al. 2014
#             total_N_recycled_best_upper_kg = total_feces_best_upper_kg*0.0056,     # see Roman et al. 2016, converted moles PON to g to kg  
#             total_P_recycled_best_low_kg = total_feces_best_upper_kg*0.0089) %>%     # whale fecal average from Ratnarajah et al. 2014
# 
#   


d_Ant_nutrients <- d_strapped_Ant_projection %>% 
  pivot_longer(cols = c(prey_mass_per_day_best_low_kg, prey_mass_per_day_best_upper_kg),
               names_to = "best_prey_est",
               values_to = "prey_mass_per_day_best_kg") %>% 
  
  # Taken from fig 3 pre-code
  filter(daily_rate >5) %>%  
  #taking the average of Dave's two distributions
  mutate(prey60_currpop = prey_mass_per_day_best_kg*`Southern hemisphere population estimate (Christensen 2006)`*60,
         prey90_currpop = prey_mass_per_day_best_kg*`Southern hemisphere population estimate (Christensen 2006)`*90,
         prey120_currpop = prey_mass_per_day_best_kg*`Southern hemisphere population estimate (Christensen 2006)`*120,
         prey150_currpop = prey_mass_per_day_best_kg*`Southern hemisphere population estimate (Christensen 2006)`*150,
         prey60_histpop = prey_mass_per_day_best_kg*`Southern hemisphere historic estimate (Christensen 2006)`*60,
         prey90_histpop = prey_mass_per_day_best_kg*`Southern hemisphere historic estimate (Christensen 2006)`*90,
         prey120_histpop = prey_mass_per_day_best_kg*`Southern hemisphere historic estimate (Christensen 2006)`*120,
         prey150_histpop = prey_mass_per_day_best_kg*`Southern hemisphere historic estimate (Christensen 2006)`*150) %>%
  pivot_longer(cols = c("prey60_currpop", "prey90_currpop", "prey120_currpop", "prey150_currpop",
                        "prey60_histpop", "prey90_histpop", "prey120_histpop", "prey150_histpop"),
               names_to = "feeding_days",
               values_to = "krill_consumed_yr") %>% 
  
  mutate(
    time_rec = ifelse(str_detect(feeding_days, "hist"), "historical", "current"),
    feeding_days = case_when(
      feeding_days %in% c("prey60_currpop", "prey60_histpop") ~ "60 days feeding",
      feeding_days %in% c("prey90_currpop", "prey90_histpop") ~ "90 days feeding",
      feeding_days %in% c("prey120_currpop", "prey120_histpop") ~ "120 days feeding",
      feeding_days %in% c("prey150_currpop", "prey150_histpop") ~ "150 days feeding")) %>% 
  
  dplyr::select(Species, feeding_days, krill_consumed_yr, time_rec) %>% 
  mutate(Fe_best_est_kg = Fe_calc_kg(krill_consumed_yr),
         N_best_est_kg = N_calc_kg(krill_consumed_yr),
         P_best_est_kg = P_calc_kg(krill_consumed_yr),
         C_respired = krill_consumed_yr*(0.45*0.75))



d_Ant_nutrients %>% filter(time_rec == "historical") %>% group_by(Species) %>% summarise(Fe_med = median(Fe_best_est_kg))

(20252140*0.75)/1000*5e4/1e6

- C_respired/1e9

pal <- c("B. bonaerensis" = "firebrick3", "B. borealis" = "goldenrod2", "B. edeni" = "darkorchid3",  "M. novaeangliae" = "gray30", "B. physalus" = "chocolate3", "B. musculus" = "dodgerblue2")

Fe_recycled_Antarctic <- ggplot(data = d_Ant_nutrients) +
  geom_boxplot(aes(x = Species, 
                   y = Fe_best_est_kg/1000,
                   fill = time_rec),
               width = .5, outlier.shape = NA, alpha = 0.5) +
  #facet_grid(time_rec~feeding_days, scales = "free", space = "free") +
  #coord_flip() + 
  #scale_fill_manual(values = pal) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "Species",
       y = bquote('Estimated Fe recycled'~(tonnes~yr^-1))) + 
  theme_classic(base_size = 20) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))
Fe_recycled_Antarctic

Carbon_production_Antarctic <- ggplot(data = d_Ant_nutrients) +
  geom_boxplot(aes(x = fct_relevel(Species, "B. bonaerensis", "M. novaeangliae", "B. physalus"), 
                   y = (((Fe_best_est_kg*0.75)/1000)*5e4/1e6),   #From Lavery et al 2010 Proc Roy Soc B
                   fill = time_rec),
               width = .5, guide = TRUE, outlier.shape = NA, alpha = 0.5) +
  #facet_grid(time_rec~feeding_days, scales = "free", space = "free") +
  #coord_flip() + 
  #scale_fill_manual(values = pal) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "Species",
       y = bquote('Fe-induced C producton'~(Mt~yr^-1))) + 
  theme_classic(base_size = 20) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) 
Carbon_production_Antarctic

Carbon_respired_Antarctic <- ggplot(data = d_Ant_nutrients) +
  geom_boxplot(aes(x = fct_relevel(Species, "B. bonaerensis", "M. novaeangliae", "B. physalus"), 
                   y = (C_respired/1e9),   #From Lavery et al 2010 Proc Roy Soc B
                   fill = time_rec),
               width = .5, outlier.shape = NA, alpha = 0.5) +
  #facet_grid(time_rec~feeding_days, scales = "free", space = "free") +
  #coord_flip() + 
  #scale_fill_manual(values = pal) +
  scale_y_log10(labels = scales::comma, limits = c(0.1, 100), breaks = c(0, 1, 10, 100)) +
  labs(x = "Species",
       y = bquote('Estimated C respired'~(Mt~yr^-1))) + 
  theme_classic(base_size = 20) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))
Carbon_respired_Antarctic


Carbon_flux_Antarctic <- ggplot(data = d_Ant_nutrients) +
  geom_boxplot(aes(x = fct_relevel(Species, "B. bonaerensis", "M. novaeangliae", "B. physalus"), 
                   y = (((Fe_best_est_kg*0.75)/1000)*5e4/1e6) - C_respired/1e9 ,   #From Lavery et al 2010 Proc Roy Soc B
                   fill = time_rec),
               width = .5, guide = TRUE, outlier.shape = NA, alpha = 0.5) +
  #facet_grid(time_rec~feeding_days, scales = "free", space = "free") +
  #coord_flip() + 
  #scale_fill_manual(values = pal) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "Species",
       y = bquote('Estimated C exported'~(Mt~yr^-1)),
       fill = "Population") + 
  theme_classic(base_size = 20) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))
Carbon_flux_Antarctic

#Doesnt work out, will do manually
C_plot <- ggarrange(Carbon_production_Antarctic, Carbon_respired_Antarctic, Carbon_flux_Antarctic,
                      labels = c("A", "B", "C"), # THIS IS SO COOL!!
                      font.label = list(size = 18),
                      legend = "none",
                      widths = c(1,2), 
                      ncol = 1, nrow = 3)
C_plot


d_summ_Ant_nutrients <- d_Ant_nutrients %>% 
  group_by(Species, time_rec) %>%
  summarise(
    Fe_best_est_t = median(Fe_best_est_kg)/1000,
    Fe_IQR_25 = round(quantile(Fe_best_est_kg, probs = 0.25, na.rm = TRUE)/1000, 2),
    Fe_IQR_75 = round(quantile(Fe_best_est_kg, probs = 0.75, na.rm = TRUE)/1000, 2),
    N_best_est_t = median(N_best_est_kg)/1000,
    N_IQR_25 = round(quantile(N_best_est_kg, probs = 0.25, na.rm = TRUE)/1000, 2),
    N_IQR_75 = round(quantile(N_best_est_kg, probs = 0.75, na.rm = TRUE)/1000, 2),
    P_best_est_t = median(P_best_est_kg)/1000,
    P_IQR_25 = round(quantile(P_best_est_kg, probs = 0.25, na.rm = TRUE)/1000, 2),
    P_IQR_75 = round(quantile(P_best_est_kg, probs = 0.75, na.rm = TRUE)/1000, 2)
    ) %>% 
  unite("Fe recycled (t pop yr) IQR", 
        c(Fe_IQR_25, Fe_IQR_75), sep = "-") %>%
  unite("N recycled (t pop yr) IQR", 
        c(N_IQR_25, N_IQR_75), sep = "-") %>%
  unite("P recycled (t pop yr) IQR", 
        c(P_IQR_25, P_IQR_75), sep = "-") %>% 
  mutate(Region = "Antarctic") 
  
  #pivot_wider(names_from = time_rec, values_from = c(Fe_best_est_t, `Fe recycled (t pop yr) IQR`))






# Creates nutrient recycling and PPR estimates for NON-ANTARCTIC data and projection ----
d_nonAnt_nutrients <- d_strapped %>% 
  filter(Region != "Antarctic") %>% 
  pivot_longer(cols = c(prey_mass_per_day_best_low_kg, prey_mass_per_day_best_upper_kg),
               names_to = "best_prey_est",
               values_to = "prey_mass_per_day_best_kg") %>% 
  
  # Taken from fig 3 pre-code
  filter(daily_rate >5) %>%  
  #taking the average of Dave's two distributions
  mutate(prey60_currpop = prey_mass_per_day_best_kg*(`Population estimate (Christensen 2006)`-`Southern hemisphere population estimate (Christensen 2006)`)*60,
         prey90_currpop = prey_mass_per_day_best_kg*(`Population estimate (Christensen 2006)`-`Southern hemisphere population estimate (Christensen 2006)`)*90,
         prey120_currpop = prey_mass_per_day_best_kg*(`Population estimate (Christensen 2006)`-`Southern hemisphere population estimate (Christensen 2006)`)*120,
         prey150_currpop = prey_mass_per_day_best_kg*(`Population estimate (Christensen 2006)`-`Southern hemisphere population estimate (Christensen 2006)`)*150,
         prey60_histpop = prey_mass_per_day_best_kg*(`Historical estimate`-`Southern hemisphere historic estimate (Christensen 2006)`)*60,
         prey90_histpop = prey_mass_per_day_best_kg*(`Historical estimate`-`Southern hemisphere historic estimate (Christensen 2006)`)*90,
         prey120_histpop = prey_mass_per_day_best_kg*(`Historical estimate`-`Southern hemisphere historic estimate (Christensen 2006)`)*120,
         prey150_histpop = prey_mass_per_day_best_kg*(`Historical estimate`-`Southern hemisphere historic estimate (Christensen 2006)`)*150) %>%
  pivot_longer(cols = c("prey60_currpop", "prey90_currpop", "prey120_currpop", "prey150_currpop",
                        "prey60_histpop", "prey90_histpop", "prey120_histpop", "prey150_histpop"),
               names_to = "feeding_days",
               values_to = "prey_consumed_yr") %>% 
  
  mutate(
    time_rec = ifelse(str_detect(feeding_days, "hist"), "historical", "current"),
    feeding_days = case_when(
      feeding_days %in% c("prey60_currpop", "prey60_histpop") ~ "60 days feeding",
      feeding_days %in% c("prey90_currpop", "prey90_histpop") ~ "90 days feeding",
      feeding_days %in% c("prey120_currpop", "prey120_histpop") ~ "120 days feeding",
      feeding_days %in% c("prey150_currpop", "prey150_histpop") ~ "150 days feeding"),
    krill_consumed_yr = case_when(Species == "B. physalus" ~ prey_consumed_yr*0.8,     # correcting for proportion of diet that is fish; and proportion of individuals not in Southern Hemisphere
                                  Species == "B. musculus" ~ prey_consumed_yr,
                                  Species == "M. novaeangliae" ~ prey_consumed_yr*0.55)) %>% 
  
  dplyr::select(Species, feeding_days, krill_consumed_yr, time_rec) %>% 
  mutate(Fe_best_est_kg = Fe_calc_kg(krill_consumed_yr),
         N_best_est_kg = N_calc_kg(krill_consumed_yr),
         P_best_est_kg = P_calc_kg(krill_consumed_yr),
         C_respired = krill_consumed_yr*(0.45*0.75)) %>% 
  filter(!is.na(krill_consumed_yr)) 


d_summ_nonAnt_nutrients <- d_nonAnt_nutrients %>% 
  group_by(Species, time_rec) %>%
  summarise(
    Fe_best_est_t = median(Fe_best_est_kg)/1000,
    Fe_IQR_25 = round(quantile(Fe_best_est_kg, probs = 0.25, na.rm = TRUE)/1000, 2),
    Fe_IQR_75 = round(quantile(Fe_best_est_kg, probs = 0.75, na.rm = TRUE)/1000, 2),
    N_best_est_t = median(N_best_est_kg)/1000,
    N_IQR_25 = round(quantile(N_best_est_kg, probs = 0.25, na.rm = TRUE)/1000, 2),
    N_IQR_75 = round(quantile(N_best_est_kg, probs = 0.75, na.rm = TRUE)/1000, 2),
    P_best_est_t = median(P_best_est_kg)/1000,
    P_IQR_25 = round(quantile(P_best_est_kg, probs = 0.25, na.rm = TRUE)/1000, 2),
    P_IQR_75 = round(quantile(P_best_est_kg, probs = 0.75, na.rm = TRUE)/1000, 2)
  ) %>% 
  unite("Fe recycled (t pop yr) IQR", 
        c(Fe_IQR_25, Fe_IQR_75), sep = "-") %>%
  unite("N recycled (t pop yr) IQR", 
        c(N_IQR_25, N_IQR_75), sep = "-") %>%
  unite("P recycled (t pop yr) IQR", 
        c(P_IQR_25, P_IQR_75), sep = "-") %>% 
  mutate(Region = "non-Antarctic")



# Combining summary tables for figure

d_nut_summ_combined <- bind_rows(d_summ_Ant_nutrients, d_summ_nonAnt_nutrients) %>% 
  pivot_longer(cols = c(Fe_best_est_t, N_best_est_t, P_best_est_t),
               names_to = "element",
               values_to = "tons excreted") %>% 
  mutate(element = ifelse(str_detect(element, "Fe"), "Total iron", 
                          ifelse(str_detect(element, "N"), "Total nitrogen", "Total phosphorus"))
         ) %>% 
  select(c(Species, time_rec, Region, element, `tons excreted`)) 
  
  

Fe_plot <-  d_nut_summ_combined %>% 
  filter(element == "Total iron") %>% 
  ggplot() +
  geom_bar(aes(x = fct_relevel(Species, "B. bonaerensis", "M. novaeangliae", "B. musculus"), 
               y = `tons excreted`, 
               fill = Region), stat = "identity") +
  facet_grid(time_rec~element, scales = "free", space = "free") +
  #scale_y_log10(labels = scales::comma) +
  scale_fill_manual(values = c("navy", "olivedrab")) +
  labs(x = "Species",
       y = bquote('Estimated quantity excreted'~(t~yr^-1))) + 
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))
Fe_plot


NandP_plot <-  d_nut_summ_combined %>% 
  filter(element != "Total iron") %>% 
  ggplot() +
  geom_bar(aes(x = fct_relevel(Species, "B. bonaerensis", "M. novaeangliae", "B. musculus"), 
               y = `tons excreted`, 
               fill = Region), stat = "identity") +
  facet_grid(time_rec~element, scales = "free", space = "free") +
  #scale_y_log10(labels = scales::comma) +
  scale_fill_manual(values = c("navy", "olivedrab")) +
  labs(x = "Species",
       y = "") + 
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))
NandP_plot

nut_plot <- ggarrange(Fe_plot, NandP_plot,
                      labels = c("A", "B"), # THIS IS SO COOL!!
                      font.label = list(size = 18),
                      legend = "none",
                      widths = c(1,2), 
                      ncol = 2, nrow = 1)
nut_plot



######################################################################################
# Estimating amount of primary production required to sustain global whale populations
######################################################################################

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


# now to plot
pal <- c("Balaenoptera bonaerensis" = "firebrick3", "Megaptera novaeangliae" = "gray30", "Balaenoptera physalus" = "chocolate3", "Balaenoptera musculus" = "dodgerblue2")



PPR_plot_CCE <- PPR_table %>% 
  ungroup %>% 
  filter(scenario %in% c("geoMeanNULL_wt_g_DC", 
                         "MedNULL_wt_g_DC", 
                         "geoMeanBOUT_wt_g_DC", 
                         "MedBOUT_wt_g_DC")) %>% 
  drop_na(NPP_to_whales_CCE) %>% 
  mutate(scenario = fct_reorder(scenario, NPP_to_whales_CCE),
         Species = fct_reorder(factor(Species), NPP_to_whales_CCE)) %>% 
  ggplot(aes(x = scenario, y = NPP_to_whales_CCE, fill = Species)) +
  geom_col(position = "stack") + 
  ylab("Proportion of NPP to whales in the CCE") + 
  theme_bw() +
  scale_fill_manual(values = pal) +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1, face="italic"),
        axis.text.y = element_text(size=12),       
        axis.title=element_text(size=13, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        strip.text.x = element_text(size = 12))
PPR_plot_CCE + theme(legend.position="top")



PPR_plot_global <- PPR_table %>% 
  ungroup %>% 
  filter(scenario %in% c("geoMeanNULL_wt_g_DC", 
                         "MedNULL_wt_g_DC", 
                         "geoMeanBOUT_wt_g_DC", 
                         "MedBOUT_wt_g_DC")) %>% 
  mutate(scenario = fct_reorder(scenario, NPP_to_whales_global),
         Species = fct_reorder(factor(Species), NPP_to_whales_global)) %>% 
  ggplot(aes(x = scenario, y = NPP_to_whales_global, fill = Species)) +
  geom_col(position = "stack") + 
  ylab("Current NPP to whales") + 
  theme_bw() +
  scale_fill_manual(values = pal) +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1, face="italic"),
        axis.text.y = element_text(size=12),       
        axis.title=element_text(size=13, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        strip.text.x = element_text(size = 12))
PPR_plot_global + theme(legend.position="none")


PPR_plot_global_historical <- PPR_table %>% 
  ungroup %>% 
  filter(scenario %in% c("geoMeanNULL_wt_g_DC", 
                         "MedNULL_wt_g_DC", 
                         "geoMeanBOUT_wt_g_DC", 
                         "MedBOUT_wt_g_DC")) %>% 
  drop_na(NPP_to_whales_global_historical) %>% 
  mutate(scenario = fct_reorder(scenario, NPP_to_whales_global_historical),
         Species = fct_reorder(factor(Species), NPP_to_whales_global_historical)) %>% 
  ggplot(aes(x = scenario, y = NPP_to_whales_global_historical, fill = Species)) +
  geom_col(position = "stack") + 
  ylab("Historical NPP to whales") + 
  theme_bw() +
  scale_fill_manual(values = pal) +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1, face="italic"),
        axis.text.y = element_text(size=12),       
        axis.title=element_text(size=13, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        strip.text.x = element_text(size = 12))
PPR_plot_global_historical + theme(legend.position="none")

NPPwhales_combined_plot <- ggarrange(PPR_plot_global, PPR_plot_global_historical, 
                                     ncol = 2, nrow = 1,
                                     common.legend = TRUE)
annotate_figure(NPPwhales_combined_plot,
                top = text_grob("Global net primary productivity (NPP) to whales", 
                                face = "bold", size = 14))



# estimating amount of food consumed over the course of a feeding season (see Lockyer 2007).
# using the equation Nr_ind = (MDC*prop_excreted)*nut_conc
nutrient_sims <- crossing(TotalPreyConsumed_kg = 18500,
                          E_blubber = seq(27e3, 37e3, by = 1e3),
                          A = 0.8,
                          beta = seq(2, 5, by = 0.5),
                          M = 42240,
                          D_feeding = seq(60, 182.5, by = 10)) %>% 
  mutate(E_preyseason = ((W_blubber*E_blubber)/A) + (((beta*293.1*M^0.75)*D_feeding)/A))

sims_plot <- sims %>% 
  ggplot(aes(E_preyseason)) +
  geom_histogram() +
  theme_classic()

sims_summ <- summarize_at(sims, vars(E_preyseason), list(mean, sd, median, min, max))













# old code below here ----

# read in data for direct predictions of daily ration (R) from literature
prey_predict <- read_excel("PreyIngestPredict.xlsx") %>% mutate(dummy = 1)
prey_predict_w_M <- tibble(M_kg = seq(5000,120000,5000), dummy =1) %>%
  full_join(prey_predict, by = "dummy") %>% 
  select(-"dummy") %>% 
  mutate(R = `Intercept (A)`*M_kg^`Exponent (B)`,
         R_compressed_90days = (R*365)/90,
         R_compressed_120days = (R*365)/120)


# read in data for predictions from BMR-->FMR of daily ration (R) from literature
#formulas: 

# krill diet percentages (from Pauly et al. 1998)
# Z-values (proportion krill in diet)
# bw = 1
# bp = 0.8
# bb - 0.9
# ba = 0.65
# bbor = 0.8
# be = 0.4
# mn = 0.55

# average weight (in kg) from Jeremy's data or Lockyer 1976
# bw = 93000
# bp = 53000
# bb/ba = 6700
# bbor = mean(c(15,18.5, 15.5,17, 17,18.5)*1000)  # From Wiki, from Lockyer 1976
# OR  mean(c(8.53, 10.25, 11.38, 11.28, 16.08, 15.56, 8.58, 13.32, 9.90, 10.61, 13.80, 15.76, 12.94, 
# #             8.89, 16.36, 12.34, 12.22, 21.62)*1000)
# be = mean(c(11.98,11.77,13.87,13.46,11.75,13.76,11.49,11.64,11.93,13.61,8.39,13.91,15.43,14.12,13.78,
#             13.91,12.55, 12.46, 11.89, 12.79, 11.32, 9.99, 16.15, 15.96, 14.85, 12.78,15.47)*1000) 
# mn = 36000


# estimate of BMR extrapolated from Kleiber 1975: 
# BMR = 293.1*(M^0.75)
# estimate of ADMR ~ FMR from Lavigne, Leaper and Lavigne 2007, Barlow et al. 2008: 
# ADMR = beta*BMR, where beat is a value between 2-5 (see Barlow et al. 2008 p.287 for more)
# estimate of R from Leaper and Lavigne 2007, Barlow et al. 2008:
# R = ADMR/(0.8*[3900*Z + 5450*(1-Z)])
# to compute prey intake per day when feeding if ALL feeding is compressed into 120 days
# (R*365)/120

prey_predict_from_BMR <- read_excel("PreyIngestPredict.xlsx", sheet = 2) %>% 
  mutate(KleiberBMR = 293.1*M_kg^0.75,
         dummy = 1)
BMRtoFMRprojection <- tibble(beta = seq(2,5,0.5), dummy = 1) %>% 
  full_join(prey_predict_from_BMR, by = "dummy") %>% 
  select(-dummy) %>% 
  mutate(ADMR = beta*KleiberBMR,
         R  = ADMR/(0.8*((3900*Z)+5400*(1-Z))),
         R_compressed_90days = (R*365)/90,
         R_compressed_120days = (R*365)/120)


# read in whale population data. Current data from IUCN 2019 Redlist, whaleing data compiled in Rocha et al. 2014
pop_data <- read_excel("Filtration project_Whale population and feeding information.xlsx", sheet = 1)

# load data
d_full_NULL <- read_csv("Cetacea model output NULL_EXTANT.csv") %>% 
  mutate(Species = ifelse(Species == "bonarensis", "bonaerensis", Species))
d_full_BOUT <- read_csv("Cetacea model output BOUT_EXTANT.csv") %>% 
  mutate(Species = ifelse(Species == "bonarensis", "bonaerensis", Species))

f_data <- read_excel("ALLPRHS 2.5.2019.xls")

tag_guide <- read_excel("TAG GUIDE_2.18.19.xlsx")
tag_guide=tag_guide[c(2:nrow(tag_guide)),] ##get rid of first two rows
colnames(tag_guide)=as.character(unlist(tag_guide[1,])) # makes first row the header
tag_guide=tag_guide[-1,] # remove first row, leaving header

v_data <- read_excel("mwmrmeasures.xlsx")
v_data <- v_data[,c(1:3)] #keeps only columns 1-3
names(v_data)[names(v_data) == 'Species'] <- 'CommonName'
v_data$L <- v_data$MW*0.9766  #creates column that converts MW (kg) to liters

# creating fucntions from Shirel's paper for MW (in kg) for engulfment capacity in liters for each species where we have a known length
bw_L <- function(x) {
  (x^3.667316*10^-0.014078)*0.9766}   # change to 0.9737098 CHANGE!!!

bp_L <- function(x) {
  (x^3.54883*10^0.15604)*0.9766}

mn_L <- function(x) {
  (x^3.24603*10^0.85934)*0.9766}

ba_L <- function(x) {
  (x^3.10910*10^0.69969)*0.9766}

be_L <- function(x) {
  (x^3.1453*10^0.5787)*0.9766} 

# creating a new column where I recalcuate all the engulfment capacities using the functions above
v_data$Recalc_L <- ifelse(v_data$CommonName == "Blue Whale", bw_L(v_data$TLm), 
                          ifelse(v_data$CommonName == "Fin Whale", bp_L(v_data$TLm),
                                 ifelse(v_data$CommonName == "Humpback Whale", mn_L(v_data$TLm),
                                        ifelse(v_data$CommonName == "Minke Whale", ba_L(v_data$TLm), be_L(v_data$TLm)))))

# create table looking at averages of MW 
v_data_species <- v_data %>% group_by(CommonName) %>% 
  summarize(Mean_L = mean(L), Mean_Recalc_L = mean(Recalc_L), 
            Med_L = median(L), Med_Recalc_L = median(Recalc_L), 
            Mean_TL = mean(TLm), Med_TL = median(TLm))

#joining f_data and tag guide to have location data
f_data <- left_join(f_data, tag_guide[ , c("ID", "Study_Area     _")], by = "ID")
names(f_data)[names(f_data) == "Study_Area     _"] <- "Study_Area"
f_data$Study_Area <- replace_na(f_data$Study_Area, "SoCal")   # Replace NAs in location with "SoCal"
f_data <- within(f_data, rm(notes)) #removing unwanted columns and rows
# cleaning data and adding columns
f_data = f_data %>% 
  mutate(Species = substr(f_data$ID,1,2),
         CommonName = ifelse(Species == "bw", "Blue Whale",
                             ifelse(Species == "bp", "Fin Whale",
                                    ifelse(Species == "mn", "Humpback Whale",
                                           ifelse(Species == "bb", "Minke Whale", "Bryde's Whale")))),
         TotalLunges = dayhourslunges + nighthourslunges + twilightzonelunges,
         TotalHours = dayhours + nighthours + twilightzone,
         LungesPerHour = TotalLunges/TotalHours,
         LungesPerDayHour = dayhourslunges/dayhours,
         LungesPerNightHour = nighthourslunges/nighthours,
         LungesPerTwHour = twilightzonelunges/twilightzone,
         Length = as.numeric(gsub(" m", "", f_data$whaleLength)),
         EmpEngulfCap = ifelse(CommonName == "Blue Whale", bw_L(Length), 
                               ifelse(CommonName == "Fin Whale", bp_L(Length),
                                      ifelse(CommonName == "Humpback Whale", mn_L(Length),
                                             ifelse(CommonName == "Minke Whale", ba_L(Length), be_L(Length)))))) %>% 
  select(-whaleLength) %>% 
  separate(prey, into = c("Prey", "Prey notes"), sep = " ") %>% 
  drop_na(Prey) %>% 
  filter(!Prey %in% c("Milk", "N")) %>% 
  mutate(PreyClean = case_when(
    Prey %in% c("Anchovies", "Fish", "Herring", "SandLance", "Sardines")  ~ "Fish-feeding",
    Prey %in% c("Inverts", "Krill") ~ "Krill-feeding"), 
    Region = case_when(
      Study_Area %in% c("Monterey", "SoCal", "Cordell Bank", "San Diego", "WA Coast")  ~ "Eastern North Pacific",
      Study_Area %in% c("Stellwagen", "Norway", "Azores", "Greenland") ~ "North Atlantic",
      Study_Area == "South Africa" ~ "South Africa",
      Study_Area == "Antarctic" ~ "Antarctic",
      Study_Area == "Chile" ~ "Chile"))
f_data$LungesPerDayHour[is.nan(f_data$LungesPerDayHour)] <- NA
f_data$LungesPerTwHour[is.nan(f_data$LungesPerTwHour)] <- NA
f_data$LungesPerNightHour[is.nan(f_data$LungesPerNightHour)] <- NA

RorqualData <- read_csv("lunge_rates_from_Paolo.csv") %>% 
  left_join(f_data, by= "ID") %>% 
  mutate(`deployment-time_h` = `deployment-time_secs`/60/60,
         SpeciesCode = substr(ID,1,2),)  %>% 
  select(-c(species)) %>% 
  arrange(dayhours, ID)
RorqualData$SpeciesCode[RorqualData$SpeciesCode == "ba"] <- "bb"

# # sweet tidy code from Max
d_sum_NULL <- d_full_NULL %>%
  group_by(Genus, Species) %>%
  summarize(wgtMeanNULL_wt_g = weighted.mean(`Prey W (g)`, Percent),
            medNULL_wt_g = median(`Prey W (g)`),
            wgtMeanNULL_E = weighted.mean(`Energy (kJ)`, Percent),
            medNULL_E = median(`Energy (kJ)`))

d_sum_BOUT = d_full_BOUT %>% 
  group_by(Genus, Species) %>% 
  summarize(wgtMeanBOUT_wt_g = weighted.mean(`Prey W (g)`, Percent), 
            medBOUT_wt_g = median(`Prey W (g)`), 
            wgtMeanBOUT_E = weighted.mean(`Energy (kJ)`, Percent),
            medBOUT_E = median(`Energy (kJ)`))

OdontoceteData <- read_csv("foragestats_combined_ko2.csv") %>% 
  separate(Species, into = c("Genus", "Species"), sep = "_") %>% 
  left_join(d_sum_NULL, by = c("Genus", "Species"))
##FIXME, I think this works though
OdontoceteData<- OdontoceteData %>% 
  left_join(d_sum_BOUT, by = c("Genus", "Species"))


# REMEMBER: OdontoceteData is all of the data we need. 
# Add rorqual data to odontocete data
cetacean_data <- left_join(OdontoceteData, RorqualData, by = c("ID")) %>%
  # feeding events are lunges, buzzes for rorquals, odontocetes
  mutate(TotalFeedingEvents = coalesce(total_lunges, total_buzz_count),
         TotalTagTime_h = coalesce(`deployment-time_h`, total_duration_h)) %>% 
  left_join(v_data_species, by = "CommonName") %>% 
  mutate(EngulfVolPerHr = LungesPerHour*Med_Recalc_L,
         EngulfVolPerDayHr = LungesPerDayHour*Med_Recalc_L) %>% 
  select(-c(Med_L, Mean_L)) %>% 
  left_join(pop_data, by = "SpeciesCode") %>% 
  #drop_na(Species) %>% 
  mutate(feeding_rate = TotalFeedingEvents / TotalTagTime_h,    # FEEDING RATE
         SpeciesCode = substr(ID,1,2),
         geoMeanNULL_wt_g_DC = case_when(
           Species =="Balaenoptera musculus" ~ 0.629646079*79862.436,        # kg to g multiply by 1000, m3 to l divide by 1000, so they cancel each other out; The second number is the Med_recalc_L from vol_data_species
           Species =="Balaenoptera physalus" ~ 0.620866221*60010.441,        # data from BaleenWhaleForagingDistBigKrill100Bins.xlsx
           Species =="Megaptera novaeangliae" ~ 0.608799363*25683.815,
           Species =="Balaenoptera bonaerensis" ~ 0.498822479*3069.295),
         MedNULL_wt_g_DC = case_when(
           Species =="Balaenoptera musculus" ~ 0.705480231*79862.436,        #kg to g multiply by 1000, m3 to l divide by 1000, so they cancel each other out; The second number is the Med_recalc_L from vol_data_species
           Species =="Balaenoptera physalus" ~ 0.657933225*60010.441,        # data from BaleenWhaleForagingDistBigKrill100Bins.xlsx
           Species =="Megaptera novaeangliae" ~ 0.657933225*25683.815,
           Species =="Balaenoptera bonaerensis" ~ 0.497702356*3069.295), 
         geoMeanBOUT_wt_g_DC = case_when(
           Species =="Balaenoptera musculus" ~ 1.665454025*79862.436,        #kg to g multiply by 1000, m3 to l divide by 1000, so they cancel each other out; The second number is the Med_recalc_L from vol_data_species
           Species =="Balaenoptera physalus" ~ 1.895174907*60010.441,        # data from BaleenWhaleForagingDistBigKrillBout100Bins.xlsx
           Species =="Megaptera novaeangliae" ~ 1.365715958*25683.815,
           Species =="Balaenoptera bonaerensis" ~ 0.555868207*3069.295),
         MedBOUT_wt_g_DC = case_when(
           Species =="Balaenoptera musculus" ~ 1.873817423*79862.436,        #kg to g multiply by 1000, m3 to l divide by 1000, so they cancel each other out; The second number is the Med_recalc_L from vol_data_species
           Species =="Balaenoptera physalus" ~ 2.15443469*60010.441,        # data from BaleenWhaleForagingDistBigKrillBout100Bins.xlsx
           Species =="Megaptera novaeangliae" ~ 1.519911083*25683.815,
           Species =="Balaenoptera bonaerensis" ~ 0.572236766*3069.295)
  )      
#rename(Species = Species.y)
cetacean_data$SpeciesCode <- sub("ba", "bb", cetacean_data$SpeciesCode)
cetacean_data$Species <- sub("Balaenoptera acutorostrata", "Balaenoptera bonaerensis", cetacean_data$Species)

Prey_consumpt_hr <- cetacean_data %>% 
  filter(TotalTagTime_h > 2, Species %in% c("Balaenoptera musculus", "Balaenoptera physalus", "Megaptera novaeangliae", "Balaenoptera bonaerensis") & sonar_exp %in% c("none", NA)) %>%  
  select(ID, Species, feeding_rate, wgtMeanNULL_wt_g:medBOUT_E, TotalTagTime_h, `Population estimate`, `Total removed`,
         geoMeanNULL_wt_g_DC, MedNULL_wt_g_DC, geoMeanBOUT_wt_g_DC, MedBOUT_wt_g_DC) %>%
  gather(scenario, prey_wgt_g, c(wgtMeanNULL_wt_g, medNULL_wt_g, wgtMeanBOUT_wt_g, medBOUT_wt_g,
                                 geoMeanNULL_wt_g_DC, MedNULL_wt_g_DC, geoMeanBOUT_wt_g_DC, MedBOUT_wt_g_DC)) %>% 
  #  gather(scenario_E, prey_E, c(wgtMeanNULL_E, medNULL_E, wgtMeanBOUT_E, medBOUT_E)) %>%  Switch these as necessary
  mutate(SpeciesCode = substr(ID,1,2), 
         hourly_prey_in_g = prey_wgt_g * feeding_rate,
         hourly_prey_in_kg = hourly_prey_in_g/1000,
         scenario_type = ifelse(scenario %in% c("wgtMeanNULL_wt_g", "medNULL_wt_g", "geoMeanNULL_wt_g_DC", "MedNULL_wt_g_DC"), "NULL", "BOUT"),
         calc_type = ifelse(scenario %in% c("wgtMeanNULL_wt_g", "wgtMeanBOUT_wt_g", "geoMeanNULL_wt_g_DC", "geoMeanBOUT_wt_g_DC"), "mean", "med")) %>% 
  mutate_if(is.character, as.factor) %>% 
  unite("scenario_calc", c("scenario_type", "calc_type"), sep = "_", remove = FALSE)
Prey_consumpt_hr$SpeciesCode <- sub("ba", "bb", Prey_consumpt_hr$SpeciesCode)
Prey_consumpt_hr$Species <- sub("Balaenoptera acutorostrata", "Balaenoptera bonaerensis", Prey_consumpt_hr$Species)
Prey_consumpt_hr <- mutate(Prey_consumpt_hr, Binomial = abbr_binom(Species))

# now for varying hours feeding within a day
prey_master_for_join <- Prey_consumpt_hr %>% 
  mutate(dummy = 1)
prey_master_varying_HrperD <- tibble(hours_feeding = seq(6,12,1), dummy = 1) %>%  # assumed a feeding day was feeding for 6-12 hours
  full_join(prey_master_for_join, by = "dummy") %>% 
  select(-dummy) %>% 
  mutate(TotalPreyConsumed_kg = hours_feeding*hourly_prey_in_kg)

excreta_day_master <- prey_master_varying_HrperD %>% 
  filter(scenario %in% c("geoMeanNULL_wt_g_DC", "MedNULL_wt_g_DC", "geoMeanBOUT_wt_g_DC", "MedBOUT_wt_g_DC")) %>% 
  mutate(total_feces_kg = TotalPreyConsumed_kg*0.8,    # From Doughty et al. 2016, Roman and McCarthy 2010
         total_Fe_recycled_kg = total_feces_kg*0.000146,    # whale fecal average from Ratnarajah et al. 2014
         total_N_recycled_kg = total_feces_kg*0.0056,     # see Roman et al. 2016, converted moles PON to g to kg  
         total_P_recycled_kg = total_feces_kg*0.0089)     # whale fecal average from Ratnarajah et al. 2014

# different scenarios
## consider that 90% of ingested iron is excreted** My 2014 paper. So consider iron reservoir in krill consumed by whales
## consider total feces excreted, likely around 75-80% wet weight of what's consumed (see Doughty et al. 2016)
### multiply total weight of feces by concentration of different elements in Fe

pal <- c("ba" = "gold3", "bb" = "firebrick3", "be" = "darkorchid3",  "mn" = "gray30", "bp" = "chocolate3", "bw" = "dodgerblue2" )
Shape <- c("ba" = 10, "bb" = 15, "be" = 8, "mn" = 17, "bp" = 18, "bw" = 19)

ExcretedFe_HrperD <- ggplot(filter(excreta_day_master, TotalTagTime_h > 2 & scenario %in% c("geoMeanNULL_wt_g_DC", "MedNULL_wt_g_DC", "geoMeanBOUT_wt_g_DC", "MedBOUT_wt_g_DC")),
                            aes(x = hours_feeding, y = total_Fe_recycled_kg, color = SpeciesCode, shape = SpeciesCode)) + 
  geom_point(inherit.aes = T, aes(group = SpeciesCode), alpha = 0.3) + 
  geom_smooth(inherit.aes = T, aes(group = SpeciesCode), size = 0.5) +
  facet_grid(Binomial~scenario_calc) +
  scale_colour_manual(name = "Species",
                      values = pal, 
                      labels = c( "B. bonaerensis", "M. novaeangliae", "B. physalus", "B. musculus")) +
  scale_shape_manual(name = "Species",
                     values = Shape,
                     labels = c( "B. bonaerensis", "M. novaeangliae", "B. physalus", "B. musculus")) +
  scale_x_continuous(breaks=seq(0, 12, 3)) +
  xlab("Hours feeding") + ylab("Total Fe excreted (kg)") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 10, face="italic"))
ExcretedFe_HrperD + theme(legend.position = "none")


# now for varying days feeding within a year
excreta_day_master_year_join <- excreta_day_master %>% 
  mutate(dummy = 1)
excreta_year_master <- tibble(days_feeding = seq(60,182.5,10), dummy = 1) %>% 
  full_join(excreta_day_master_year_join, by = "dummy") %>% 
  select(-dummy) %>% 
  mutate(TotalAnnualPreyConsumed_kg = days_feeding*TotalPreyConsumed_kg,
         total_annual_feces_kg = TotalAnnualPreyConsumed_kg*0.8,    # From Doughty et al. 2016, Roman and McCarthy 2010
         total_annual_Fe_recycled_kg = total_annual_feces_kg*0.000146,    # whale fecal average from Ratnarajah et al. 2014
         total_annual_N_recycled_kg = total_annual_feces_kg*0.0056,     # see Roman et al. 2016, converted moles PON to g to kg  
         total_annual_P_recycled_kg = total_annual_feces_kg*0.0089)    # whale fecal average from Ratnarajah et al. 2014

NutRecTableDisc <- excreta_year_master %>% filter(days_feeding == 120) %>% 
  group_by(Binomial) %>% 
  summarise(Individual_Fe = median(total_annual_Fe_recycled_kg),
            Individual_Fe_SE = SE(total_annual_Fe_recycled_kg),
            Population_Fe = median(total_annual_Fe_recycled_kg*`Population estimate`),
            Population_Fe_SE = SE(total_annual_Fe_recycled_kg*`Population estimate`),
            Individual_N = median(total_annual_N_recycled_kg),
            Individual_N_SE = SE(total_annual_N_recycled_kg),
            Population_N = median(total_annual_N_recycled_kg*`Population estimate`),
            Population_N_SE = SE(total_annual_N_recycled_kg*`Population estimate`),
            Individual_P = median(total_annual_P_recycled_kg),
            Individual_P_SE = SE(total_annual_P_recycled_kg),
            Population_P = median(total_annual_P_recycled_kg*`Population estimate`),
            Population_P_SE = SE(total_annual_P_recycled_kg*`Population estimate`))
View(NutRecTableDisc)

ExcretedFe_DperYr <- ggplot(filter(excreta_year_master, TotalTagTime_h > 2 & scenario %in% c("geoMeanNULL_wt_g_DC", "MedNULL_wt_g_DC", "geoMeanBOUT_wt_g_DC", "MedBOUT_wt_g_DC")),
                            aes(x = days_feeding, y = total_annual_Fe_recycled_kg, color = SpeciesCode, shape = SpeciesCode)) + 
  geom_point(inherit.aes = T, aes(group = SpeciesCode), alpha = 0.3) + 
  geom_smooth(inherit.aes = T, aes(group = SpeciesCode), size = 0.5) +
  facet_grid(Binomial~scenario_calc) +
  scale_colour_manual(name = "Species",
                      values = pal, 
                      labels = c( "B. bonaerensis", "M. novaeangliae", "B. physalus", "B. musculus")) +
  scale_shape_manual(name = "Species",
                     values = Shape,
                     labels = c( "B. bonaerensis", "M. novaeangliae", "B. physalus", "B. musculus")) +
  scale_x_continuous(breaks=seq(0, 12, 3)) +
  xlab("days feeding") + ylab("Total Fe excreted (kg)") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 10, face="italic"))
ExcretedFe_DperYr + theme(legend.position = "none")

