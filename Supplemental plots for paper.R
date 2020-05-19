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
library(ggmap)
library(maps)


# formula for standard error
SE = function(x){sd(x)/sqrt(sum(!is.na(x)))}

# Abbreviate a binomial e.g. Balaenoptera musculus -> B. musculus
abbr_binom = function(binom) {
  paste(str_sub(binom, 1, 1), 
        str_extract(binom, " .*"), 
        sep = ".")
}


# Supplemental map figure of tag deployments ----

pal <- c("B. bonaerensis" = "firebrick3", "B. borealis" = "goldenrod2", 
         "B. edeni" = "darkorchid3",  "M. novaeangliae" = "gray30", 
         "B. physalus" = "chocolate3", "B. musculus" = "dodgerblue2",
         "B. mysticetus" = "#4304ff", "E. glacialis" = "#ff9f04")


bm <- data.frame(
  lon = rep(seq(-62, -57, length.out = 3), each = 2), 
  lat = rep(seq(68, 70, length.out = 2)),
  Species = "Balaena mysticetus"
)

eg <- data.frame(
  lon = rep(seq(-50, -40, length.out = 5), each = 4), 
  lat = rep(seq(38, 44, length.out = 4)),
  Species = "Eubalaena glacialis"
)
eg2 <- data.frame(
  lon = seq(-47.5, -42.5, length.out = 3), 
  lat = 46,
  Species = "Eubalaena glacialis"
)

bw_az <- data.frame(
  lon = -25, 
  lat = 40,
  Species = "Balaenoptera musculus"
)

bw_cce <- data.frame(
  lon = rep(seq(-150, -125, length.out = 14), each = 8), 
  lat = rep(seq(26, 40, length.out = 8)),
  Species = "Balaenoptera musculus"
)

bp_cce <- data.frame(
  lon = rep(seq(-145, -130, length.out = 7), each = 3), 
  lat = rep(seq(42, 45, length.out = 3)),
  Species = "Balaenoptera physalus"
)

bp_gr <- data.frame(
  lon = -59.5, 
  lat = 65,
  Species = "Balaenoptera physalus"
)

bp_sw <- data.frame(
  lon = seq(-62, -57, length.out = 3), 
  lat = 38,
  Species = "Balaenoptera physalus"
)

bp_az <- data.frame(
  lon = rep(seq(-22, -28, length.out = 3), each = 2), 
  lat = rep(seq(36, 38, length.out = 2)),
  Species = "Balaenoptera physalus"
)


be <- data.frame(
  lon = rep(seq(3, 15, length.out = 7)), 
  lat = -34,
  Species = "Balaenoptera edeni"
)

bb <- data.frame(
  lon = rep(seq(-80, -65, length.out = 6), each = 3), 
  lat = rep(seq(-60, -63, length.out = 3)),
  Species = "Balaenoptera bonaerensis"
)
bb2 <- data.frame(
  lon = -62.25, 
  lat = -60,
  Species = "Balaenoptera bonaerensis"
)

mn_wap <- data.frame(
  lon = rep(seq(-90, -74, length.out = 7), each = 3), 
  lat = rep(seq(-65, -69, length.out = 3)),
  Species = "Megaptera novaeangliae"
)

mn_cce <- data.frame(
  lon = rep(seq(-160, -154, length.out = 4), each = 20), 
  lat = rep(seq(25, 50, length.out = 20)),
  Species = "Megaptera novaeangliae"
)
mn_cce2 <- data.frame(
  lon = rep(seq(-159, -155, length.out = 3)), 
  lat = 51.32,
  Species = "Megaptera novaeangliae"
)

mn_sa <- data.frame(
  lon = rep(seq(6, 12, length.out = 4), each = 2), 
  lat = rep(seq(-30, -32, length.out = 2)),
  Species = "Megaptera novaeangliae"
)

mn_sw <- data.frame(
  lon = rep(seq(-72, -56, length.out = 7), each = 3), 
  lat = rep(seq(32, 36, length.out = 3)),
  Species = "Megaptera novaeangliae"
)

mn_nor <- data.frame(
  lon = 20, 
  lat = 68,
  Species = "Megaptera novaeangliae"
)

whale_points <- bind_rows(bm, eg, eg2, bw_az, bw_cce, bp_cce,
                          bp_gr, bp_sw, bp_az, be, bb, bb2,
                          mn_wap, mn_sw, mn_sa, mn_cce, mn_cce2)


map.world = map_data(map="world")

map <- ggplot() + 
  geom_map(data=map.world,
           map=map.world,aes(map_id=region,x=long,y=lat), 
           fill = "gray20") +
  #coord_cartesian() +
  geom_point(aes(x = lon, y = lat, color = abbr_binom(Species)), 
             data = whale_points, size = 1.2) +
  geom_segment(aes(x=-52, y=43, xend=-69, yend=43), size = 0.25) +
  geom_segment(aes(x=-60, y=39, xend=-69, yend=43), size = 0.25) +
  geom_segment(aes(x=-69, y=38, xend=-69, yend=43), size = 0.25) +
  labs(color = "Species") +
  guides(color = guide_legend(override.aes = list(size=2))) +
  scale_colour_manual(values = pal) +
  theme(legend.justification=c(0,0), 
        legend.position=c(0.1,0.2),
        legend.text = element_text(face = "italic",
                                   size = 12),
        legend.title = element_text(size = 14),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.box.background = element_rect(colour = "black")) +
  theme_nothing(legend = TRUE)
map 


dev.copy2pdf(file="Map of tag deployments.pdf", width=15, height=8)



##########################################
# load, clean, combine and summarize data----
##########################################

# read in data for direct predictions of daily ration (R) from literature
prey_predict <- read_excel("PreyIngestPredict.xlsx") %>% mutate(dummy = 1)
prey_predict_w_M <- tibble(M_kg = seq(5000,120000,5000), dummy =1) %>%
  full_join(prey_predict, by = "dummy") %>% 
  select(-"dummy") %>% 
  mutate(R = `Intercept (A)`*M_kg^`Exponent (B)`,
         R_compressed_90days = (R*365)/90,
         R_compressed_100days = (R*365)/100,
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
         R_compressed_90days = (0.83*(R*365))/90,
         R_compressed_100days = (0.83*(R*365))/100,
         R_compressed_120days = (0.83*(R*365))/120)


# estimating amount of food consumed over the course of a feeding season (see Lockyer 2007).
# using the equation I developed with Dave Cade:E_preyseason = ((W_blubber*E_blubber)/A) + (((beta*293.1*M^0.75)*D_feeding)/A)
sims <- crossing(W_blubber = 18500,
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







###########################################
# preliminary plots for prey consumption
########################################### 
# plot predictions from literature, with ours
pal <- c("B. acutorostrata" = "gold3", "B. bonaerensis" = "firebrick3", "B. borealis" = "navy", "B. edeni" = "darkorchid3",  "M. novaeangliae" = "gray30", "B. physalus" = "chocolate3", "B. musculus" = "dodgerblue2")
Shape <- c("B. acutorostrata" = 10, "B. bonaerensis" = 15, "B. edeni" = 8, "B. borealis" = 7, "M. novaeangliae" = 17, "B. physalus" = 18, "B. musculus" = 19)


Eq1_ingestion_prediction <- BMRtoFMRprojection %>% 
  group_by(Species, M_kg) %>% 
  summarise_at(vars(KleiberBMR:R_compressed_120days), 
               list(~mean(.), ~min(.), ~max(.)), na.rm = TRUE) 


ingest_MSpredict_plot_R <- ggplot(Eq1_ingestion_prediction, 
                                aes(color = abbr_binom(Species))) +
  geom_pointrange(aes(x = M_kg, y = R_mean, 
                      ymin = R_min, ymax = R_max)) + 
  scale_x_log10(labels = scales::comma,
                limits = c(6000,105000)) +
  scale_y_log10(labels = scales::comma,
                limits = c(100,10000),
                breaks = c(100,500,1000,5000,10000)) +
  scale_colour_manual(values = pal) +
  #scale_shape_manual(values = Shape) +
  labs(title = "Daily Ration (R)",
       x = "Body mass (kg)",
       y = bquote('Prey consumed'~(kg~d^-1)), 
       color = "Species") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5))
ingest_MSpredict_plot_R


ingest_MSpredict_plot_MDC120 <- ggplot(Eq1_ingestion_prediction, 
                                       aes(color = abbr_binom(Species))) +
  geom_pointrange(aes(x = M_kg, y = R_compressed_120days_mean, 
                      ymin = R_compressed_120days_min, ymax = R_compressed_120days_max), shape = 17) + 
  scale_x_log10(labels = scales::comma,
                limits = c(6000,105000)) +
  scale_y_log10(labels = scales::comma,
                limits = c(100,10000),
                breaks = c(100,500,1000,5000,10000)) +
  scale_colour_manual(values = pal) +
  #scale_shape_manual(values = Shape) +
  labs(title = "MDC - 120 days feeding",
       x = "Body mass (kg)",
       y = "",
       color = "Species") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5)) 
ingest_MSpredict_plot_MDC120


ingest_MSpredict_plot_MDC90 <- ggplot(Eq1_ingestion_prediction, 
                                  aes(color = abbr_binom(Species))) +
  geom_pointrange(aes(x = M_kg, y = R_compressed_90days_mean, 
                      ymin = R_compressed_90days_min, ymax = R_compressed_90days_max), shape = 15) + 
  scale_x_log10(labels = scales::comma,
                limits = c(6000,105000)) +
  scale_y_log10(labels = scales::comma,
                limits = c(100,10000),
                breaks = c(100,500,1000,5000,10000)) +
  scale_colour_manual(values = pal) +
  #scale_shape_manual(values = Shape) +
  labs(title = "MDC - 90 days feeding",
       x = "Body mass (kg)",
       y = "",
       color = "Species") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5)) 
ingest_MSpredict_plot_MDC90



PreyConsumpt_Eq1 <- ggarrange(ingest_MSpredict_plot_R, ingest_MSpredict_plot_MDC120, ingest_MSpredict_plot_MDC90, 
                             labels = c("A", "B", "C"), # THIS IS SO COOL!!
                             font.label = list(size = 18),
                             common.legend = TRUE, legend="bottom",
                             ncol = 3, nrow = 1)
PreyConsumpt_Eq1

dev.copy2pdf(file="PreyConsumpt_Eq1.pdf", width=18, height=6)






ingest_directpredict_plot_R <- ggplot(prey_predict_w_M, aes(M_kg, R)) +
  geom_line(aes(color = str_wrap(`Reference(s)`, 75)), size = 1.15) +
  
  scale_x_log10(labels = scales::comma)+
  scale_y_log10(labels = scales::comma,
                limits = c(20,18000),
                breaks = c(100,500,1000,5000,10000,15000)) +
  scale_colour_manual(values = c("black", "dark gray", "dimgray", "salmon2", 
                                 "slategray1", "thistle", "palegreen3", "skyblue3")) +
  labs(title = "Daily Ration (R)",
       x = "Body mass (kg)",
       y = bquote('Prey consumed'~(kg~d^-1)), 
       color = "Reference(s)") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        plot.title = element_text(hjust = 0.5)) 
ingest_directpredict_plot_R


ingest_directpredict_plot_R120 <- ggplot(prey_predict_w_M, aes(M_kg, R_compressed_120days)) +
  geom_line(aes(color = str_wrap(`Reference(s)`, 20)), size = 1.15) +
  
  scale_x_log10(labels = scales::comma)+
  scale_y_log10(labels = scales::comma,
                limits = c(20,18000),
                breaks = c(100,500,1000,5000,10000,15000)) +
  scale_colour_manual(values = c("black", "dark gray", "dimgray", "salmon2", 
                                 "slategray1", "thistle", "palegreen3", "skyblue3")) +
  labs(title = "MDC - 120 days feeding",
       x = "Body mass (kg)",
       y = "",
       color = "Reference(s)") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        plot.title = element_text(hjust = 0.5)) 
ingest_directpredict_plot_R120


ingest_directpredict_plot_R90 <- ggplot(prey_predict_w_M, aes(M_kg, R_compressed_90days)) +
  geom_line(aes(color = str_wrap(`Reference(s)`, 20)), size = 1.15) +
  
  scale_x_log10(labels = scales::comma)+
  scale_y_log10(labels = scales::comma,
                limits = c(20,18000),
                breaks = c(100,500,1000,5000,10000,15000)) +
  scale_colour_manual(values = c("black", "dark gray", "dimgray", "salmon2", 
                                 "slategray1", "thistle", "palegreen3", "skyblue3")) +
  labs(title = "MDC - 90 days feeding",
       x = "Body mass (kg)",
       y = "",
       color = "Reference(s)") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        plot.title = element_text(hjust = 0.5)) 
ingest_directpredict_plot_R90


PreyConsumpt_Eq2 <- ggarrange(ingest_directpredict_plot_R, ingest_directpredict_plot_R120, ingest_directpredict_plot_R90, 
                              labels = c("A", "B", "C"), # THIS IS SO COOL!!
                              font.label = list(size = 18),
                              common.legend = TRUE, legend="bottom",
                              ncol = 3, nrow = 1)
PreyConsumpt_Eq2

#Save pdf of plot
dev.copy2pdf(file="PreyConsumpt_Eq2.pdf", width=18, height=6)
       


# Mass Specific Energy Intake ----
## of all species and all prey types!

load("daynights_krill.RData") # This keeps its original name "krill_daily" here
load("daynights_fish.RData") # This keeps its original name "fish_daily" here
load("daynights_AntProj.RData") # This keeps its original name "krill_daily_Ant_proj" here
load("balaenid_sampled.RData")

krill_daily_Ant_proj_forjoin <- krill_daily_Ant_projection %>% 
  filter(SpeciesCode %in% c("bw", "bp")) %>% 
  mutate(prey_type = "Euphausia superba")
fish_daily <- fish_daily %>% 
  mutate(prey_type = "forage fish")
krill_daily <- krill_daily %>% 
  mutate(prey_type = case_when(region == "Polar" ~ "E. superba",
                               region == "Temperate" ~ "non-E. superba krill"))

balaenid_sampled <- balaenid_sampled %>% 
  mutate(Mass_est_kg = case_when(Species == "Balaena mysticetus" ~ 70000,
                                 Species == "Eubalaena glacialis" ~ 55000),
         Mass_specifc_energy_intake_best_high_kJ = prey_t_longday*1000*6620/Mass_est_kg,
         Mass_specifc_energy_intake_best_low_kJ = prey_t_shortday*1000*6620/Mass_est_kg,
         prey_type = "copepods")

combined_en <- balaenid_sampled  %>% 
  bind_rows(
    select(krill_daily,Species, Mass_specifc_energy_intake_best_low_kJ, Mass_specifc_energy_intake_best_high_kJ, prey_type),
    select(fish_daily, prey_type, Species, Mass_specifc_energy_intake_best_low_kJ, Mass_specifc_energy_intake_best_high_kJ),
    #select(krill_daily_Ant_proj_forjoin, prey_type, Species, Mass_specifc_energy_intake_best_low_kJ, Mass_specifc_energy_intake_best_high_kJ)
  ) %>% 
  rename(`Lower best estimate` = Mass_specifc_energy_intake_best_low_kJ,
         `Upper best estimate` = Mass_specifc_energy_intake_best_high_kJ) %>% 
  pivot_longer(`Lower best estimate`:`Upper best estimate`,
               names_to = "est_type", values_to = "Mass_specifc_energy_intake_kJ") %>% 
  mutate(est_type = fct_relevel(factor(est_type), "Low"))


MS_En_in <- ggplot(combined_en, 
              aes(fill = abbr_binom(Species))) +
  geom_boxplot(aes(x = abbr_binom(Species), 
                   y = Mass_specifc_energy_intake_kJ),
               outlier.shape = NA) +
  facet_grid(.~prey_type, scales = "free") +
  geom_hline(yintercept = 80, linetype = "dashed") +
  scale_fill_manual(values = pal) +
  labs(x = "Species",
       y = bquote('Mass specific energy intake'~(kJ~kg^-1~d^-1)),
       fill = "Species") + 
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                   face = "italic"),
        legend.text = element_text(face = "italic")) +
  scale_y_log10(labels = scales::comma, limits = c(50, 5000),
                breaks = c(50,100,250,500,1000,2500,5000))
MS_En_in

  #Save pdf of plot
dev.copy2pdf(file="Mass_specific_En_in.pdf", width=14, height=6)
 
