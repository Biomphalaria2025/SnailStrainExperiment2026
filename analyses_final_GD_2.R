#### Code associated with the snail dataset to test for 
#### temporal differences between four B. glabrata snail strains infected with S. mansoni 
#### Guillaume Delhaye. Last edited  05/08/2026

#Article title: Optimising the Biomphalaria glabrata - Schistosoma mansoni 
#laboratory life-cycle: improving efficiency via enhancing snail-schistosome compatibility

####RStudio version used: 2026.05.1
####R version used:R version 4.6.0 
####RTools version used:4.5 

# install.packages(c("tidyverse", "readxl", "janitor", "survival", "survminer", 
#                    "Rcpp", "brms", "tidybayes","here"))

library(tidyverse)
library(readxl)
library(janitor)
library(survival)
library(survminer)
library(Rcpp)
library(brms)
library(tidybayes)
library(here)
library(ggsurvfit)

## options to use STAN (not always necessary)
pkgbuild::find_rtools()
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
options(mc.cores = 4)

# Confirm project location
here()

## Create an output folder for figures 
dir.create(here("output_2026"), showWarnings = FALSE)

##### 1. load and clean the data ------------------------------------------------
infected_raw <- read_xlsx(here("strains_data.xlsx"), sheet = 1)
non_infected_raw <- read_xlsx(here("strains_data.xlsx"), sheet = 2)

# clean the data 
data <- data.frame(infected_raw) |>
  clean_names() |>
  # change name of BB022 into BB19
  mutate(strain=ifelse(strain == "BB022", "BB19", strain)) |>
  # transform format of wrong columns
  mutate(weeks = as.numeric(weeks), 
         infected = as.factor("I"),
         snail_id = as.factor(snail_id), 
         strain = as.factor(strain), 
         alive_y_n = as.factor(alive_y_n), 
         shell_size_mm = as.numeric(shell_size_mm), 
         shedding_y_n = as.factor(shedding_y_n), 
         no_cerc_1 = as.numeric(no_cerc_1), 
         no_cerc_2 = as.numeric(no_cerc_2),
         no_cerc_3 = as.numeric(no_cerc_3), 
         notes = as.character(notes))|> 
  select(date, weeks, infected, snail_id, strain, alive_y_n, 
         shell_size_mm, shedding_y_n, no_cerc_1, 
         no_cerc_2, no_cerc_3, notes) |>
  # for weeks 1 and 2, the dilution is lower (6.67) than weeks afte (44.44)
  mutate(cerc1_conv = case_when(weeks %in% c(1) ~ no_cerc_1 * 6.67, 
                                weeks == 2 & strain == "BB19" ~ no_cerc_1 * 6.67,
                                weeks == 2 & strain != "BB19" ~ no_cerc_1 * 44.44,
                                weeks > 2 ~ no_cerc_1 * 44.44),
         
         cerc2_conv = case_when(weeks %in% c(1) ~ no_cerc_2 * 6.67, 
                                weeks == 2 & strain == "BB19" ~ no_cerc_2 * 6.67,
                                weeks == 2 & strain != "BB19" ~ no_cerc_2 * 44.44,
                                weeks > 2 ~ no_cerc_2 * 44.44),
         
         cerc3_conv = case_when(weeks %in% c(1) ~ no_cerc_3 * 6.67, 
                                weeks == 2 & strain == "BB19" ~ no_cerc_3 * 6.67,
                                weeks == 2 & strain != "BB19" ~ no_cerc_3 * 44.44,
                                weeks > 2 ~ no_cerc_3 * 44.44)) |>
  rowwise() |>
  mutate(cerc_aver = mean(c(cerc1_conv, cerc2_conv,cerc3_conv)),
         cerc_sd = sd(c(cerc1_conv, cerc2_conv,cerc3_conv)),
         shell_size_mm = ceiling(shell_size_mm)) |> 
  ungroup()|>
  # Attach the data for non-infected
  bind_rows(non_infected_raw |>
              clean_names()  |>
              # change name of BB022 into BB19
              mutate(strain=ifelse(strain == "BB022", "BB19", strain)) |>
             # filter(weeks != 3)|> # remove week 3 beause of measurement problems
              mutate(date = as.Date(date, format = "%d.%m.%Y"),
                     infected = as.factor("NI"),
                     snail_id = as.factor(snail_id), 
                     strain = as.factor(strain), 
                     alive_y_n = as.factor(alive_y_n), 
                     shell_size_mm = as.numeric(shell_size_mm), 
                     notes = as.character(notes))) |>
  filter(weeks != 7) |>
  mutate(shell_size_mm = ifelse(weeks %in% c(1, 6), shell_size_mm, NA)) |>
  mutate(dead = abs(as.numeric(alive_y_n)-2))
str(data)
summary(data)

##############################################################################
# 2. Infection prevalence ---------------------
###########################################################################
# Some individuals died before the week one of the experiment, 
# First, check how many individuals are alive at week 1
alive_w_1 <- data |>
  filter(infected == "I")|>
  filter(weeks == 1) |> 
  summarise(alive_at_beginning = sum(alive_y_n == "Y"), 
            .by = c(strain, infected)) 

## Number of infected individuals
number_infected <- data |>
  filter(infected == "I")|> 
  group_by(snail_id) |>
  fill(shedding_y_n, .direction = "updown") |>
  #calculate if an individual has been infected (i.e. has 
  #shed cercariae at least once before a given date)
  mutate(cs = (cumsum(as.numeric(shedding_y_n)-1)>0)*1) |> 
  ungroup() |>
  summarise(cs = sum(cs, na.rm = T), 
            .by = c(strain, weeks)) |>
  left_join(alive_w_1) |>
  # Calculate the proportion of infecte as the number of individuals
  # that have produced cercariae / number of individuals alive at week 1
  mutate(prop_infected = (cs/alive_at_beginning)*100)

#########################################################
#### SUPPLEMENTARY MATERIAL S1 - ####FIGURE#####
#################################################################
# plot figure of infection rate. 
ggplot(number_infected, 
       aes(y = prop_infected, x = weeks, colour = strain)) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 1.2) +
  scale_y_continuous(
    name = "Percentage of infected individuals",
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = c(0, 0)) +
  xlab("Weeks") +
  labs(colour = "Strain") +
  theme_light(base_size = 20) +
  theme(
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 18)
  )

#SAVING AS TIFF
ggsave("output_2026/percentage_infected.tiff",
  dpi = 600, compression = "lzw")

##SAVING AS Jpeg
ggsave("output_2026/percentage_infected.jpeg")


##################################################################
## Test if infection prevalence is different between strains at week 6
number_infected_6 <- number_infected |>
  filter(weeks == 6) #|>

# Binomial model testing if proportion of infected individuals at week 6 
# is different between strain
mod_infected <- brm(cs|trials(alive_at_beginning) ~ strain, 
                    data = number_infected_6,
                    family = binomial(link = "logit"))

summary(mod_infected)

# Test for pairwise differences (a > b), with real probability
hypothesis(mod_infected, "Intercept + strainBB19 > Intercept + strainNo2")
hypothesis(mod_infected, "Intercept + strainBB19 > Intercept + strainLSHTM")
hypothesis(mod_infected, "Intercept + strainBB19 > Intercept")
hypothesis(mod_infected, "Intercept > Intercept + strainNo2")
hypothesis(mod_infected, "Intercept > Intercept + strainLSHTM")
hypothesis(mod_infected, "Intercept + strainLSHTM > Intercept + strainNo2")

############################################################################
## 3. Cercarial output ------------------------------------------------
##########################################################################
# many 0 so use a hurdle negative binomial distribution
data2 = data.frame(data) |> 
  filter(infected == "I" & alive_y_n == "Y")|> # only for infected individuals
  group_by(snail_id) |>
  mutate(cum_cerc = ceiling(cumsum(cerc_aver)), 
         strain = as.character(strain), 
         snail_id = as.character(snail_id)) |>
  ungroup()
str(data2)

## Graph 
data2 |>
  ggplot(aes(y = cum_cerc, x = weeks, colour = snail_id)) +
  geom_point() +
  geom_line() +
  facet_wrap(~strain) + 
  theme(legend.position = "none")

## Run model with hurdle negative bionmial family
mod_cerc_cum_hu_negbin <- brm(
  bf(cum_cerc ~ weeks * strain,
     hu ~ strain),
  data = data2,
  family = hurdle_negbinomial(),
  chains = 2,
  iter = 4000, 
  control = list(max_treedepth = 10, 
                 adapt_delta = 0.80))

mod_cerc_cum_hu_negbin

saveRDS(mod_cerc_cum_hu_negbin, "output_2026/model_cerc_hurdle_negbin.rds")

#Plotting Fig 2 - cercarial output
plot(conditional_effects(mod_cerc_cum_hu_negbin), 
     points = TRUE, 
     effects = "weeks*strain",
     ask = FALSE)[[3]] +
  facet_wrap(~strain) +
  theme_light() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 12)
  ) +
  ylab("Cumulative cercariae production") + 
  scale_x_continuous(breaks = 1:6) +
  xlab("Weeks")

##SAVING AS LOW QUALITY JPEG
ggsave("output_2026/cumulative_cerc_production_poisson.jpeg")

##SAVING AS HIGH QUALITY TIFF
ggsave("output_2026/cumulative_cerc_production_poisson.tiff",dpi = 600,
  compression = "lzw")

### extract fixed effects
(fe <- fixef(mod_cerc_cum_hu_negbin))
data.frame(t(fe)) |>
  transmute(BB02 = exp(weeks), 
            BB19 = exp(weeks + weeks.strainBB19), 
            LSHTM = exp(weeks + weeks.strainLSHTM), 
            No2 = exp(weeks + weeks.strainNo2))

### Test for differences in slopes between populations
# BB19 - BB02
hypothesis(mod_cerc_cum_hu_negbin, "weeks  + weeks:strainBB19 > weeks")
# BB19 - LSHTM
hypothesis(mod_cerc_cum_hu_negbin, "weeks  + weeks:strainBB19 > weeks  + weeks:strainLSHTM")
# BB19 - No2
hypothesis(mod_cerc_cum_hu_negbin, "weeks  + weeks:strainBB19 > weeks  + weeks:strainNo2")
# BB02 - LSHTM
hypothesis(mod_cerc_cum_hu_negbin, "weeks > weeks  + weeks:strainLSHTM")
# BB02 - No2
hypothesis(mod_cerc_cum_hu_negbin, "weeks > weeks  + weeks:strainNo2")
# LSHTM - No2
hypothesis(mod_cerc_cum_hu_negbin, "weeks  + weeks:strainLSHTM > weeks  + weeks:strainNo2")

## Extract the non zero part of the model -> cannot be used because does not 
## include the hurdle part of the model
# Number of cerc produced by BB19 at week 6
exp(fe[1,1] + fe[4,1] + 6*(fe[3,1] + fe[7,1]))
# Number of cerc produced by BB02 at week 6
exp(fe[1,1] + 6*(fe[3,1]))
# Number of cerc produced by LSHTM at week 6
exp(fe[1,1] + fe[5,1] + 6*(fe[3,1] + fe[8,1]))
# Number of cerc produced by No2 at week 6
exp(fe[1,1] + fe[5,1] + 6*(fe[2,1] + fe[8,1]))

## Extract the mean predicted values at week 6 using the posterior distributions
# create a data frame with the specific values for week 6
pred_val <- expand_grid(strain = unique(data2$strain), 
                        weeks = c(1:6))

cerc_week_6 <- mod_cerc_cum_hu_negbin |> 
  predicted_draws(newdata = pred_val) |>
  filter(weeks == 6) |>
  ungroup()

# Extract and calculate the mean
cerc_week_6 |>
  group_by(strain) |>
  summarise(mean(.prediction))

# Test if there are differences at week 6 between cercariae produced by pairs of populations
# BB19 > LSHTM 
hypothesis(mod_cerc_cum_hu_negbin, 
           "Intercept + strainBB19 + 6*(weeks  + weeks:strainBB19) > 
            Intercept + strainLSHTM + 6*(weeks  + weeks:strainLSHTM)" )

# BB19 > No2 
hypothesis(mod_cerc_cum_hu_negbin, 
           "Intercept + strainBB19 + 6*(weeks  + weeks:strainBB19) > 
            Intercept + strainNo2 + 6*(weeks  + weeks:strainNo2)" )

# BB19 > BB02 
hypothesis(mod_cerc_cum_hu_negbin, 
           "Intercept + strainBB19 + 6*(weeks  + weeks:strainBB19) > 
            Intercept +  6*(weeks)" )

# BB02 > LSHTM 
hypothesis(mod_cerc_cum_hu_negbin, 
           "Intercept +  6*(weeks) > 
            Intercept + strainLSHTM + 6*(weeks  + weeks:strainLSHTM)" )

# BB02 > No2 
hypothesis(mod_cerc_cum_hu_negbin, 
           "Intercept +  6*(weeks) > 
           Intercept + strainNo2 + 6*(weeks  + weeks:strainNo2)")


# LSHTM > No2 
hypothesis(mod_cerc_cum_hu_negbin, 
           "Intercept + strainLSHTM + 6*(weeks  + weeks:strainLSHTM) > 
            Intercept + strainNo2 + 6*(weeks  + weeks:strainNo2)" )


#### Graphs ---
## Cercariae production
# Raw data
(cerc_production <- 
    ggplot(data |> 
             filter(infected == "I" & alive_y_n == "Y"), 
           aes(y = cerc_aver, x = weeks)) +
    geom_point(alpha = 0.2) + 
    facet_wrap(~ strain) +
    geom_line(aes(colour =  snail_id), alpha = 0.2) +
    geom_smooth(data = data |> filter(infected == "I" & cerc_aver > 0),
                aes(y = cerc_aver, x = weeks),
                method = "lm") +
    scale_x_continuous(name="Time (weeks)", limits=c(1, 6), n.breaks = 6) +
    scale_y_continuous(name="Cercariae production per individual") +
    theme_minimal() +
    theme(legend.position = 'none'))
ggsave("output/cerc_prod_time.jpeg", height = 4, width = 4)


#################################################################
# 4. MORTALITY  ------------------------------------------------
######################################################################
# Test of Average survival between infected and non infected over all strains
(survival_analysis_I_NI <- survdiff(Surv(weeks, dead) ~ infected, 
                                    data = data))

# figure of mortality risk for all strains split to infected and uninfected - NOT USED IN PAPER
survfit2(Surv(weeks, dead) ~ infected, 
         data = data) |>
  ggsurvfit(type = "risk") +
  scale_ggsurvfit()  +
  add_confidence_interval() +
  ylab("Mortality risk") +
  xlab("Time (weeks)")
ggsave("output_2026/mean_mortality_risk.jpeg", height = 4, width = 4)

################################################################
############### FIGURE 4 - MORTALITY OVER TIME ==============
################################################################
# BB19
(survival_analysis_BB19 <- survdiff(Surv(weeks, dead) ~ infected, 
                                    data = data |> filter(strain == "BB19")))

mort_BB19 <- survfit2(Surv(weeks, dead) ~ infected, 
                      data = data |> filter(strain == "BB19")) |>
  ggsurvfit(type = "risk") +
  scale_ggsurvfit(y_scales = list(limits = c(0,1), 
                                  breaks = seq(0, 1, by = 0.2))) +
  add_confidence_interval() +
  ylab("") +
  xlab("") +
  ggtitle("BB19") + 
  theme(plot.title = element_text(hjust = 0.5))


# BB02
(survival_analysis_BB02 <- survdiff(Surv(weeks, dead) ~ infected, 
                                    data = data |> filter(strain == "BB02")))

mort_BB02 <- survfit2(Surv(weeks, dead) ~ infected, 
                      data = data |> filter(strain == "BB02")) |>
  ggsurvfit(type = "risk") +
  scale_ggsurvfit(y_scales = list(limits = c(0,1), 
                                  breaks = seq(0, 1, by = 0.2))) +
  add_confidence_interval() +
  ylab("Mortality risk") +
  xlab("") +
  ggtitle("BB02") + 
  theme(plot.title = element_text(hjust = 0.5))

# LSHTM
(survival_analysis_LSHTM <- survdiff(Surv(weeks, dead) ~ infected, 
                                     data = data |> filter(strain == "LSHTM")))
mort_LSHTM <- survfit2(Surv(weeks, dead) ~ infected, 
                       data = data |> filter(strain == "LSHTM")) |>
  ggsurvfit(type = "risk") +
  scale_ggsurvfit(y_scales = list(limits = c(0,1), 
                                  breaks = seq(0, 1, by = 0.2))) +
  add_confidence_interval() +
  ylab("Mortality risk") +
  xlab("Time (weeks)") +
  ggtitle("LSHTM") + 
  theme(plot.title = element_text(hjust = 0.5))

# No2
(survival_analysis_No2 <- survdiff(Surv(weeks, dead) ~ infected, 
                                   data = data |> filter(strain == "No2")))

mort_No2 <- survfit2(Surv(weeks, dead) ~ infected, 
                     data = data |> filter(strain == "No2")) |>
  ggsurvfit(type = "risk") +
  scale_ggsurvfit(y_scales = list(limits = c(0,1), 
                                  breaks = seq(0, 1, by = 0.2))) +
  add_confidence_interval() +
  ylab("") +
  xlab("Time (weeks)") +
  ggtitle("No2") + 
  theme(plot.title = element_text(hjust = 0.5))

#generate one fraph with the 4 panels for each strain
ggpubr::ggarrange(mort_BB02, mort_BB19, mort_LSHTM, mort_No2, 
                  legend = "none")

#### SAVE AS JPEG
ggsave("output_2026/strain_mortality_risk.jpeg", height = 5, width = 5)

###SAVE AS TIFF
ggsave("output_2026/strain_mortality_risk.tiff",
       dpi = 600, compression = "lzw")

#######################################################################

# Test of survival between infected strains
(survival_analysis_all_I <- survdiff(Surv(weeks, dead) ~ strain, 
                                     data = data |> filter(infected == "I")))

## Try different pairs of strains here
(survival_analysis_all_I <- survdiff(Surv(weeks, dead) ~ strain, 
                                     data = data |> 
                                       filter(infected == "I" & 
                                                strain %in% c("No2", "LSHTM"))))
(survival_analysis_all_I <- survdiff(Surv(weeks, dead) ~ strain, 
                                     data = data |> 
                                       filter(infected == "I" & 
                                                strain %in% c("No2", "BB02"))))
(survival_analysis_all_I <- survdiff(Surv(weeks, dead) ~ strain, 
                                     data = data |> 
                                       filter(infected == "I" & 
                                                strain %in% c("No2", "BB19"))))
(survival_analysis_all_I <- survdiff(Surv(weeks, dead) ~ strain, 
                                     data = data |> 
                                       filter(infected == "I" & 
                                                strain %in% c("LSHTM", "BB02"))))
(survival_analysis_all_I <- survdiff(Surv(weeks, dead) ~ strain, 
                                     data = data |> 
                                       filter(infected == "I" & 
                                                strain %in% c("LSHTM", "BB19"))))
(survival_analysis_all_I <- survdiff(Surv(weeks, dead) ~ strain, 
                                     data = data |> 
                                       filter(infected == "I" & 
                                                strain %in% c("BB02", "BB19"))))

###############################################################################

### ADDITIONAL FIGURES
### Mortality rates figure - NOT USED
data |>
  select(weeks, snail_id, strain, alive_y_n, infected) |>
  group_by(weeks, strain, infected) |>
  summarise(n_dead = sum(alive_y_n == "N", na.rm = TRUE)) |>
  ggplot(aes(y = n_dead, x = weeks, colour = infected)) +
  geom_line(linewidth = 1) +
  facet_wrap(~strain) +
  scale_colour_discrete(name = "", 
                        labels = c("Infected", "Non infected")) +
  scale_x_continuous(name="Time (weeks)", limits=c(1, 6), n.breaks = 6) +
  scale_y_continuous(name="Dead individuals", limits=c(0, 18), n.breaks = 6)+
  theme_minimal() +
  theme(legend.position="none")
ggsave("output_2026/number_dead_time.jpeg")

##########################################################################
#######FIGURE SUPLEMENTARY 2##############################################
###########################################################################

### Does cumulative cercariae production predict mortality probability ?
## do individuals producing more cerc die more? and is it true for all populations? 
## cumulative cercariae production
## modify data -> replace NA in cumulative cerc production by the last value while alive
## then take the sqrt to decrease variance in predictor for model convergence
data2 <- data |> 
  group_by(snail_id) |>
  mutate(cum_cerc_prod = sqrt(cumsum(replace_na(cerc_aver, 0)))) |>
  ungroup()|>
  filter(weeks == 6) |>
  mutate(alive = as.numeric(as.factor(alive_y_n))-1) |>
  select(weeks, infected, snail_id, strain, alive, cerc_aver, cum_cerc_prod)
str(data2)

# check for a few snails
data2 |> filter(snail_id == "INo2_8")
data2 |> filter(snail_id == "IBB02_16")
data2 |> filter(snail_id == "IBB02_20")

## Model
mod_cerc_mort <- brm(alive ~ 1 + (1 + cum_cerc_prod | strain), 
                     family = bernoulli(), data = data2, 
                     iter = 6000, 
                     control = list(max_treedepth = 19, 
                                    adapt_delta = 0.99))
pairs(mod_cerc_mort)

summary(mod_cerc_mort) 
pp_check(mod_cerc_mort, ndraws = 100)

### formatted figure
get_variables(mod_cerc_mort)

##SUPPLEMENTARY FIGURE 2
mod_cerc_mort |>
  spread_draws(r_strain[condition, te]) |>
  median_qi(.width = c(.95, .66)) |>
  filter(te == "cum_cerc_prod") |>
  ggplot(aes(y = condition, x = r_strain, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  geom_vline(xintercept = 0, linewidth = 1.5, colour = "red") +
  theme_light() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)) +
  ylab("") +
  xlab("Influence of cercariae production on survival")

##Saving as JPEG
ggsave("output_2026/cercariae_production.jpeg", height = 5, width = 5)

###SAVE AS TIFF
ggsave("output_2026/cercariae_production.tiff",
       dpi = 600, compression = "lzw")

### No clear effect of Cerc production on mortality.


############################################################################
## 4. SHELL SIZE ------------------------------------------------
##########################################################################
## model testing if there is a difference in growth rate between infected and non-infected
mod_growth_diff <- brm(shell_size_mm ~ weeks * infected + (1 + weeks | snail_id), 
                       data = data, family = gaussian,
                       prior = c(prior(lognormal(0, 0.2), class = sigma)), 
                       iter = 8000, 
                       control = list(max_treedepth = 15, 
                                      adapt_delta = 0.99))
pp_check(mod_growth_diff)

get_variables(mod_growth_diff)

mod_growth_diff |>
  spread_draws(b_weeks, `b_weeks:infectedNI`) |>
  transmute(slope_infected = b_weeks, 
            slope_non_infected = b_weeks + `b_weeks:infectedNI`) |>
  median_qi(.width = c(.95, .8, .5)) |> view()

## create a model of shell size through time while controlling for 
## non independence of individual observations (random effect per individual snail)

## Test if there are differences in growth rates between the different strains 
## of infected snails (what we are actually interested in)
mod_shell_size_I <- brm(shell_size_mm ~ weeks * strain + (1 + weeks | snail_id), 
                        data = data |> 
                          filter(infected == "I"), 
                        family = gaussian,
                        prior = c(prior(lognormal(0, 0.2), class = sigma)), 
                        iter = 8000, control = list(max_treedepth = 15, adapt_delta = 0.99))

## we see that BB19 grows faster than BB02 (positive estimate, CI does not contain 0)
## and No2 grows slower than BB02 (negative estimate, CI does not contain 0)
pp_check(mod_shell_size_I)
fe_growth <- fixef(mod_shell_size_I)

# Average growth rate per week (in mm)
growth_rate_infected <- data.frame(strain = c("BB02", "BB19", "LSTHM", "No2"), 
                                   growth_rate = round(c(fe_growth[2, 1], 
                                                         fe_growth[2, 1] + fe_growth[6, 1], 
                                                         fe_growth[2, 1] + fe_growth[7, 1],
                                                         fe_growth[2, 1] + fe_growth[8, 1]),2))
growth_rate_infected ##average Growth results for each strains

# Test if there are significant differences in the average slope (contrasts)
hypothesis(mod_shell_size_I, "weeks + weeks:strainBB19 > weeks")
hypothesis(mod_shell_size_I, "weeks + weeks:strainBB19 > weeks + weeks:strainLSHTM")
hypothesis(mod_shell_size_I, "weeks + weeks:strainBB19 > weeks + weeks:strainNo2")
hypothesis(mod_shell_size_I, "weeks + weeks:strainLSHTM > weeks + weeks:strainNo2")
hypothesis(mod_shell_size_I, "weeks + weeks:strainLSHTM > weeks")
hypothesis(mod_shell_size_I, "weeks + weeks:strainNo2 > weeks")

## Test if there is a difference in slope for the non-infected ones (negative 
## control, important to see if there are some genetic differences impacting 
## growth independantly from infection
mod_shell_size_NI <- brm(shell_size_mm ~ weeks * strain + (1 + weeks | snail_id), 
                         data = data |> 
                           filter(infected == "NI"),# select non-infected individuals
                         family = gaussian,
                         prior = c(prior(lognormal(0, 0.2), class = sigma)), 
                         iter = 8000, control = list(max_treedepth = 15, 
                                                     adapt_delta = 0.99))

pp_check(mod_shell_size_NI)
fixef(mod_shell_size_NI)

hypothesis(mod_shell_size_NI, "weeks + weeks:strainBB19 > weeks + weeks:strainNo2")
hypothesis(mod_shell_size_NI, "weeks + weeks:strainLSHTM > weeks + weeks:strainNo2 ")
hypothesis(mod_shell_size_NI, "weeks > weeks + weeks:strainNo2 ")
hypothesis(mod_shell_size_NI, "weeks + weeks:strainBB19 > weeks + weeks:strainLSHTM")
hypothesis(mod_shell_size_NI, "weeks + weeks:strainBB19 > weeks")
hypothesis(mod_shell_size_NI, "weeks + weeks:strainLSHTM > weeks")

################################################################
############### FIGURE 4 - SHELL GROWTH CHANGE ==============
################################################################

## Growth rate 
# infected 
(growth_infected <- data |> 
   #  filter(weeks %in% c(1,6)) |>
   filter(infected == "I") |>
   ggplot(aes(y = shell_size_mm, x = weeks)) +
   geom_point() + 
   facet_wrap(~ strain) +
   geom_line(aes(colour =  snail_id), alpha = 0.3) +
   geom_smooth(method = "lm", se = FALSE, linewidth = 2) +
   theme_minimal(base_size = 22) +
   theme(legend.position = 'none'))
ggsave("output_2026/shell_size_time.jpeg")

# non-infected
(growth_non_infected <- data |> 
    # filter(weeks %in% c(1,6)) |>
    filter(infected == "NI") |>
    ggplot(aes(y = shell_size_mm, x = weeks)) +
    geom_point() + 
    facet_wrap(~ strain) +
    geom_line(aes(colour =  snail_id)) +
    ### confidence intervals (the gray around lines) are off (se = FALSE)
    geom_smooth(method = "lm", se = FALSE, linewidth = 2) +
    theme_minimal(base_size = 22) +
    theme(legend.position = 'none'))

# Growth infected vs non-infected
(growth_i_ni <- data |> 
    #  filter(weeks %in% c(1,6)) |>
    ggplot(aes(y = shell_size_mm, x = weeks, 
               colour = infected, alpha = 0.2)) + 
    geom_jitter(width = 0.2, alpha = 0.2) +
    facet_wrap(~ strain) +
    ### confidence intervals (the gray around lines) are off (se = FALSE)
    geom_smooth(method = "lm", se = FALSE, linewidth = 2) +
    scale_colour_discrete(name = "", 
                          labels = c("Infected", "Non infected")) +
    scale_x_continuous(name="Time (weeks)", limits=c(1, 6), n.breaks = 6) +
    scale_y_continuous(name="Shell size (mm)", limits=c(5, 20))+
    theme_minimal(base_size = 22) +
    theme(legend.position="none"))

##SAVING AS LOW QUALITY JPEG
ggsave("output_2026/shell_size_time.jpeg")

##SAVING AS HIGH QUALITY TIFF
ggsave("output_2026/shell_size_time.tiff",dpi = 600,
       compression = "lzw")
############### END ################################################