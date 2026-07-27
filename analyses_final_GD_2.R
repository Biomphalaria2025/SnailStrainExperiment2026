#### Code associated with the snail dataset to test for 
#### temporal differences between 4 B. glabrata snail strains infected with S. mansoni 
#### Guillaume Delhaye 24-09-2023

#Article title: Optimising the Biomphalaria glabrata - Schistosoma mansoni 
#laboratory life-cycle: improving efficiency via enhancing snail-schistosome compatibility

##Before analysis install these programs:

####RStudio version used: 2026.05.1
####R version used:R version 4.6.0 
####RTools version used:4.5 

install.packages(c("tidyverse", "readxl", "janitor", "survival", "survminer", 
                   "Rcpp", "brms", "tidybayes","here")) #just once

library(tidyverse)
library(readxl)
library(janitor)
library(survival)
library(survminer)
library(Rcpp)
library(brms)
library(tidybayes)
library(here)

# Confirm project location
here()

# Source custom functions
source(here("hurdle_gaussian.R"))
options(mc.cores = 4)

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
  group_by(strain, infected) |>
  summarise(alive_at_beginning = sum(alive_y_n == "Y"))

## Number of infected individuals
number_infected <- data |>
  filter(infected == "I")|> 
  group_by(snail_id) |>
  fill(shedding_y_n, .direction = "updown") |>
  #complicated way to calculate if an individual has been infected (i.e. has 
  #shed cercariae at least once before a given date)
  mutate(cs = (cumsum(as.numeric(shedding_y_n)-1)>0)*1) |> 
  ungroup() |>
  group_by(strain, weeks) |>
  summarise(cs = sum(cs, na.rm = T)) |>
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
    expand = c(0, 0)
  ) +
  
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

##TEACHING R where to find these STAN tools (not always necessary)
pkgbuild::find_rtools()
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")

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

# many 0 and otherwise a lognormal distribution (right skewed 
# with few very high values). 
# Use a zero-inflated (or hurdle) lognormal model -> NOT USED because the model
# does not fit the data. However, I would leave it because it is what would 
# be used traditionally. I used a hurdle gaussian on sqrt transformed data 

# ## Try with custom hurdle gaussian on sqrt response data (because lognormal overestimate 
# ## low values and underestimate high). Function in "hurdle_gaussian.R"
# ## Thanks to Andrew Heiss DOI: 10.59350/ety2j-09566
 model_cerc_gauss <- bf(cerc_aver  ~ weeks*strain + (1 | snail_id), # production > 0 depends on strain and time
                        hu ~ weeks) #(1 | snail_id)) # number of 0 depends on strain and time
# 
 priors_gauss = c(prior(normal(0, 10), class = Intercept, lb = 0),
                  prior(normal(0, 1000), class = b, coef = "weeks"),
                  prior(lognormal(0, 5), class = sigma)
)

 mod_cerc_gauss <- brm(model_cerc_gauss,
                       data = data |> 
                         mutate(weeks = weeks - 1, 
                                cerc_aver = sqrt(cerc_aver)) |> # set week 1 as the intercept for the model
                         filter(infected == "I" & alive_y_n == "Y"), # only for infected individuals
                       family = hurdle_gaussian,
                       prior = priors_gauss,                
                       iter = 10000, 
                       control = list(max_treedepth = 15, 
                                      adapt_delta = 0.99),
                       stanvars = stanvars,
                       sample_prior = "yes"
 )
 
 ### RUN THIS BEFORE PROGRESSING FURTHER cum_cerc data set generated
data2 = data.frame(data) |> 
  #mutate(weeks = weeks - 1) |> # set week 1 as the intercept for the model
  filter(infected == "I" & alive_y_n == "Y")|> # only for infected individuals
  group_by(snail_id) |>
  mutate(cum_cerc = ceiling(cumsum(cerc_aver))) |>
  ungroup()
str(data2)

##############################################################################
##################F I G U R E 2 - Change in cercarial output over time. #######
###############################################################################

#this Figure version not used at the end
data2 |>
  ggplot(aes(y = cum_cerc, x = weeks, colour = snail_id)) +
  geom_point() +
  geom_line() +
  facet_wrap(~strain) + 
  theme(legend.position = "none")
  

mod_cerc_cum_poisson <- brm(
  bf(cum_cerc ~ weeks*strain),
                      data = data2,
                      family = poisson(),
                      chains = 2,
                      iter = 4000, 
                      control = list(max_treedepth = 15, 
                                     adapt_delta = 0.90),
                      stanvars = stanvars,
                      sample_prior = "yes")
saveRDS(mod_cerc_cum_poisson, "output_2026/model_poisson.rds")

#Plotting Fig 2 - cercarial output
pp_check(mod_cerc_cum_poisson) + xlim (0, 50000)
plot(conditional_effects(mod_cerc_cum_poisson), 
     points = TRUE, 
     effects = "weeks*strain",
     facet_args = list("strain"), 
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

#############################################################################
 mod_cerc_lognorm <- brm(cerc_aver  ~ weeks*strain + (1 + weeks | snail_id),
                       data = data |> 
                         mutate(weeks = weeks) |> # set week 1 as the intercept for the model
                         filter(infected == "I" & alive_y_n == "Y"), # only for infected individuals
                       family = hurdle_lognormal(),
                       #prior = priors_gauss,                
                       iter = 2000, 
                       control = list(max_treedepth = 15, 
                                      adapt_delta = 0.95),
                       stanvars = stanvars,
                       sample_prior = "yes"
 )
#### GRAPH NOT USED
 pp_check(mod_cerc_lognorm, ndraws = 100) + 
   xlim(0, 50)
 plot(conditional_effects(mod_cerc_lognorm), 
      points = TRUE, 
      effects = "weeks*strain",
      facet_args = list("strain"), 
      ask = F)[[3]] 

### extract fixed effects
(fe <- fixef(mod_cerc_cum_poisson))
data.frame(t(fe)) |>
  transmute(BB02 = exp(weeks), 
         BB19 = exp(weeks + weeks.strainBB19), 
         LSHTM = exp(weeks + weeks.strainLSHTM), 
         No2 = exp(weeks + weeks.strainNo2))
  
### Test for differences in slopes between populations
# BB19 - BB02
hypothesis(mod_cerc_cum_poisson, "weeks  + weeks:strainBB19 > weeks")
# BB19 - LSHTM
hypothesis(mod_cerc_cum_poisson, "weeks  + weeks:strainBB19 > weeks  + weeks:strainLSHTM")
# BB19 - Mo2
hypothesis(mod_cerc_cum_poisson, "weeks  + weeks:strainBB19 > weeks  + weeks:strainNo2")
# BB02 - LSHTM
hypothesis(mod_cerc_cum_poisson, "weeks > weeks  + weeks:strainLSHTM")
# BB02 - No2
hypothesis(mod_cerc_cum_poisson, "weeks > weeks  + weeks:strainNo2")
# LSHTM - No2
hypothesis(mod_cerc_cum_poisson, "weeks  + weeks:strainLSHTM > weeks  + weeks:strainNo2")

# Number of cerc produced by BB02 at week 6
exp(fe[1,1] + 6*(fe[2,1]))
# Number of cerc produced by BB19 at week 6
exp(fe[1,1] + fe[3,1] + 6*(fe[2,1] + fe[6,1]))
# Number of cerc produced by LSHTM at week 6
exp(fe[1,1] + fe[4,1] + 6*(fe[2,1] + fe[7,1]))
# Number of cerc produced by No2 at week 6
exp(fe[1,1] + fe[5,1] + 6*(fe[2,1] + fe[8,1]))

#TESTS INCOMPLETE?
# Test if there are differences at week 6 between cercariae produced by pairs of populations
# BB19 > LSHTM 
hypothesis(mod_cerc_cum_poisson, 
           "exp(Intercept + strainBB19 + 6*(weeks  + weeks:strainBB19)) > 
            exp(Intercept + strainLSHTM + 6*(weeks  + weeks:strainLSHTM))" )
#-> 8128.25, p > 0.99

# BB19 > No2 
hypothesis(mod_cerc_cum_poisson, 
           "exp(Intercept + strainBB19 + 6*(weeks  + weeks:strainBB19)) > 
            exp(Intercept + strainNo2 + 6*(weeks  + weeks:strainNo2))" )


# BB19 > BB02 
hypothesis(mod_cerc_cum_poisson, 
           "exp(Intercept + strainBB19 + 6*(weeks  + weeks:strainBB19)) > 
            exp(Intercept +  6*(weeks))" )


# BB02 > LSHTM 
hypothesis(mod_cerc_gauss, 
           "(Intercept +  6*(weeks))^2 > 
            (Intercept + strainLSHTM + 6*(weeks  + weeks:strainLSHTM))^2" )
# -> 1495.39, p = 0.86

# BB02 > No2 
hypothesis(mod_cerc_gauss, 
           "(Intercept +  6*(weeks))^2 > 
           (Intercept + strainNo2 + 6*(weeks  + weeks:strainNo2))^2")
# -> 2318.55, p = 0.96

# LSHTM > No2 
hypothesis(mod_cerc_gauss, 
           "(Intercept + strainLSHTM + 6*(weeks  + weeks:strainLSHTM))^2 > 
            (Intercept + strainNo2 + 6*(weeks  + weeks:strainNo2))^2" )
# -> 823.16, p = 0.72

#### Graphs --- NOT USED - testing for the best visualisation of cercarial production
plot(conditional_effects(mod_cerc_gauss))
## Cercariae production GRAPH NOT USED
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

####GRAPH NOT USED
(cercaria_prod_model <- data |>
    filter(infected == "I" & alive_y_n == "Y") |> 
    ggplot(aes(y = cerc_aver, x = weeks)) + 
    geom_jitter(width = 0.1, alpha = 0.1) + 
    geom_smooth(method = "lm") +
    facet_wrap(~strain)+
    scale_x_continuous(name="Weeks", limits=c(1, 6), n.breaks = 7) +
    scale_y_continuous(name="Cercariae production", limits=c(0, 25000), n.breaks = 6)+
    theme_minimal())
ggsave("output/cercariae_production.jpeg")

############################################################################
## 4. SHELL SIZE ------------------------------------------------
##########################################################################

## model testing if there is a diference in growth rate between infected and non-infected
mod_growth_diff <- brm(shell_size_mm ~ weeks * infected + (1 + weeks | snail_id), 
                       data = data, family = gaussian,
                       prior = c(prior(lognormal(0, 0.2), class = b, coef = "weeks"),
                                 prior(lognormal(0, 0.2), class = sigma)), 
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
                        prior = c(prior(lognormal(0, 0.2), class = b, coef = "weeks"),
                                  prior(lognormal(0, 0.2), class = sigma)), 
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
hypothesis(mod_shell_size_I, "weeks + weeks:strainBB19 = weeks + weeks:strainLSHTM")
hypothesis(mod_shell_size_I, "weeks + weeks:strainBB19 = weeks + weeks:strainNo2")
hypothesis(mod_shell_size_I, "weeks + weeks:strainLSHTM = weeks + weeks:strainNo2")

## Test if there is a difference in slope for the non-infected ones (negative 
## control, important to see if there are some genetic differences impacting 
## growth independantly from infection
mod_shell_size_NI <- brm(shell_size_mm ~ weeks * strain + (1 + weeks | snail_id), 
                         data = data |> 
                           filter(infected == "NI"),# select non-infected individuals
                         family = gaussian,
                         prior = c(prior(lognormal(0, 0.2), class = b, coef = "weeks"),
                                   prior(lognormal(0, 0.2), class = sigma)), 
                         iter = 8000, control = list(max_treedepth = 15, 
                                                     adapt_delta = 0.99))

pp_check(mod_shell_size_NI)
fixef(mod_shell_size_NI)

hypothesis(mod_shell_size_NI, "weeks + weeks:strainBB19 = weeks + weeks:strainLSHTM")
hypothesis(mod_shell_size_NI, "weeks + weeks:strainNo2 = weeks + weeks:strainBB19 ")
hypothesis(mod_shell_size_NI, "weeks + weeks:strainNo2 = weeks + weeks:strainLSHTM")

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
ggsave("output/shell_size_time.jpeg")

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

#################################################################
#5. MORTALITY - SURVIVAL - DEATHS ------------------------------------------------
######################################################################

library(ggsurvfit)
# library("survival")
# library("survminer")
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

#now lets do it for each strain
################################################################
############### FIGURE 4 - MORTALITY OVER TIME ==============
################################################################
# and for each strain (little statistical power because no replicate)
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
ggsave("output/number_dead_time.jpeg")

##########################################################################
#######FIGURE SUPLEMENTARY 2##############################################
###########################################################################

### Does cumulative cercariae production predict mortality probability ?
## do individuals producing more cerc die more? and is it true for all populations? 
## cumulative cercariae production
data |> 
  group_by(snail_id) |>
  mutate(cum_cerc_prod = cumsum(cerc_aver)) |>
  ggplot(aes(y = cum_cerc_prod, x = weeks)) +
  geom_point() + 
  facet_wrap(~ strain) +
  geom_line(aes(colour =  snail_id)) +
  theme_minimal() +
  theme(legend.position = 'none')

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
data2 |> filter(snail_id == "IBB19_1") # OK

###sometimes there are problems in communication between R and RTools. 
##follow the steps below to fix it but don't use if all works well
#make sure you use RTools 4.5 first (otherwise the model below wont run)

#add this to force R to use Rtools 4.5 (only if probem appears)
# writeLines('PATH="${RTOOLS45_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")

#restart R - important
#test if RTools is installed by using - should return OK
# pkgbuild::has_build_tools(debug = TRUE)

mod_cerc_mort <- brm(alive ~ 1 + (1 + cum_cerc_prod | strain), 
                     family = bernoulli(), data = data2, 
                     iter = 6000, 
                     control = list(max_treedepth = 19, 
                                    adapt_delta = 0.99))
pairs(mod_cerc_mort)

summary(mod_cerc_mort) ####it reports one divergent. but the model is good enough
pp_check(mod_cerc_mort)

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
    axis.text = element_text(size = 12)
  ) +
  ylab("") +
  xlab("Influence of cercariae production on survival")

##Saving as JPEG
ggsave("output/cercariae_production.jpeg", height = 5, width = 5)

###SAVE AS TIFF
ggsave("output/cercariae_production.tiff",
       dpi = 600, compression = "lzw")

### No clear effect of Cerc production on mortality.


# ######################## Exploratory figures --------------------------------
# # check distribution of the data
# 
# ## shell size
# # infected
# data |> 
#   filter(infected == "I") |>
#   ggplot(aes(shell_size_mm)) +
#   geom_histogram(stat = "count") + 
#   facet_wrap(~strain)
# ggsave("output/hist_shell_I.jpeg")
# 
# # non-infected
# data |> 
#   filter(infected == "NI") |>
#   ggplot(aes(shell_size_mm)) +
#   geom_histogram(stat = "count") + 
#   facet_wrap(~strain)
# ggsave("output/hist_shell_NI.jpeg")
# 
# ## cercaria production (only infected)
# data |> 
#   filter(infected == "I")|>
#   ggplot(aes(((cerc_aver)))) +
#   geom_histogram(bins = 10) + 
#   facet_wrap(~strain)
# ggsave("output/hist_log_cerc_I.jpeg")
