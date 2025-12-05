################################################################################
#########  Discrete-event simulation of a sex-specific 3-state model   #########
################################################################################

# Fernando Alarid-Escudero, PhD (1) 
# Jeremy Goldhaber-Fiebert, PhD (1)

# Affiliations: 		
# 1 Department of Health Policy, School of Medicine, and Stanford Health Policy, 
#   and Freeman-Spogli Institute for International Studies, Stanford University, 
#   Stanford, California, USA

################################################################################

rm(list = ls()) # clear memory (removes all the variables from the workspace)

# 01 Load packages & functions ----
library(data.table)
library(reshape2)
library(dplyr)

source(file = "code/functions_des.R")

# 02 Input Model Parameters ----

set.seed(1)                      # set the seed 

# Model structure
n_i <- 100000 # number of individuals
n_t <- 60     # number of cycles
v_names_cycles <- paste("cycle", 0:n_t, sep = " ") # names of cycles
v_names_states  <- c("H", "S", "D") # state names
n_s  <- length(v_names_states)  # number of states
d_e <- 0.03 # discount rate for QALYs at 3%
d_c <- 0.03 # discount rate for costs at 3%

# Event rates by sex
r_HS_female <- 0.05  # rate to become Sick when Healthy for females
hr_HS_male  <- 1.10  # hazard ratio for males
r_SD_female <- 0.10  # rate to die when Sick for females
hr_SD_male  <- 0.90  # hazard ratio for males

## 02.1 Background mortality ----
# Probabilities of death from healthy state is age-dependent
# load probability of death by age from Human Mortalitiy Database
dt_hmd_usa_2015 <- read.csv("data/HMD_USA_Mx_2015.csv")[, -c(1, 7)] 
max_age <- max(dt_hmd_usa_2015$Age)
dt_hmd_usa_2015_long <- reshape2::melt(dt_hmd_usa_2015, 
                                       id.vars = c("Year", "Age"), 
                                       variable.name = "Sex", 
                                       value.name = "death_rate")

# Generate a wide/long block-diagonal data.frame with age-specific mortality 
# rates conditional on being alive at each year of age
df_p_HD <- obtain_probs_des(dt_hmd_usa_2015_long)

# Convert data frame to a data table for efficiency
dt_p_HD <- data.table(df_p_HD) # convert data frame to a data table for efficiency
dt_p_HD
setkey(dt_p_HD, Year, Age, Sex)      # set the data table to be indexed by Age and Sex

# Population characteristics
p_female <- 0.51                # proportion population who are female

# Costs and utilities inputs 
c_H  <- 400                     # cost of remaining one cycle Healthy
c_S  <- 1000                    # cost of remaining one cycle Sick
c_D  <- 0                       # cost of remaining one cycle Dead
u_H  <- 0.8                     # utility when Healthy 
u_S  <- 0.5                     # utility when Sick
u_D  <- 0                       # utility when Dead

# 03 Sample individual level characteristics ----
set.seed(2024025)
## 03.1 Static characteristics ----
# Store these characteristics in a single data frame for convenient referencing

# Initial Age
age_dist <- read.csv("data/MyPopulation_AgeDistribution.csv")   # load age distribution
plot(prop ~ age, age_dist)
v_age0 <- sample(x    = age_dist$age,
                 prob = age_dist$prop,
                 size = n_i, 
                 replace = TRUE)         # sample from age distribution
# Sex
v_sex <- sample(x    = c("Female", "Male"),
                prob = c(p_female, 1 - p_female),
                size = n_i,
                replace = TRUE)                 # randomly sample sex

# Store these static individual characteristics in one data frame
dt_Pop <- data.table(Ind  = 1:n_i, 
                     Year = 2015,
                     Age = v_age0, 
                     Sex  = v_sex)
dt_Pop[, SexNum := fifelse(Sex == "Female", 
                           yes = 0, no = 1)]
dt_Pop
t_init <- Sys.time()
# Compute individual- and age-specific probabilities of dying from other causes
m_p_HD <- as.matrix(dt_p_HD[.(dt_Pop$Year, dt_Pop$Age, dt_Pop$Sex), DP0:DP110])

# 04 Discrete-event simulation (DES) ----
## 04.1 Sample age of death as a function of Yea, Age, and Sex ----
dt_Pop[, Age_death_oc := nps_nhppp(m_probs = m_p_HD, correction = "uniform")]

## 04.2 Sample time and age of getting sick ----
dt_Pop[, Time_to_sick := rexp(n = .N, 
                              rate = exp(log(r_HS_female) + 
                                           log(hr_HS_male)*SexNum))]
dt_Pop[, Age_to_sick := Age + Time_to_sick]

## 04.3 Sample time and age of dying from Sick ----
dt_Pop[, Time_to_die_from_sick := rexp(n = .N, 
                                       rate = exp(log(r_SD_female) + 
                                                    log(hr_SD_male)*SexNum))]
dt_Pop[, Age_to_die_from_sick := Age_to_sick + Time_to_die_from_sick]
dt_Pop[, Got_sick := (Age_to_sick < Age_death_oc)]

## 04.4 Age and cause of death ----
dt_Pop[, Age_death := matrixStats::rowMins(cbind(Age_death_oc, 
                                                 Age_to_die_from_sick))]
dt_Pop[, Cause_of_death := fifelse(test = (Age_death == Age_death_oc),
                                   yes  = "oc",
                                   no   = "dc")]
## 04.5 Years in healthy ----
# Undiscounted
dt_Pop[, ly_healthy := matrixStats::rowMins(cbind(Age_to_sick, 
                                                  Age_death)) - Age]
# Discounted for QALYs
dt_Pop[, ly_healthy_disc := calc_discount(x = 1, 
                                          disc_factor = d_e,
                                          time_init = 0,
                                          time_fin = ly_healthy)]
# Discounted for costs
dt_Pop[, ly_healthy_disc_costs := calc_discount(x = 1, 
                                                disc_factor = d_c,
                                                time_init = 0,
                                                time_fin = ly_healthy)]
## 04.5 Years in Sick ----
# Undiscounted
dt_Pop[, ly_Sick := matrixStats::rowMaxs(cbind(Age_death - Age_to_sick, 
                                               0))]
# Discounted for QALYs
dt_Pop[, ly_Sick_disc := calc_discount(x = 1, 
                                       disc_factor = d_e,
                                       time_init = Age_to_sick,
                                       time_fin = Age_to_sick + ly_Sick)]
# Discounted for costs
dt_Pop[, ly_Sick_disc_costs := calc_discount(x = 1, 
                                             disc_factor = d_e,
                                             time_init = Age_to_sick,
                                             time_fin = Age_to_sick + ly_Sick)]
## 04.6 Years alive ----
dt_Pop[, `:=`(ly_total = ly_healthy + ly_Sick,
              ly_total_disc = ly_healthy_disc + ly_Sick_disc)]
if (sum((dt_Pop$Age_death - dt_Pop$Age) != dt_Pop$ly_total) > 1) { 
  "Time alive is not calculated correctly"
}

# 05 Compute QALYs ----
### Healthy
dt_Pop[, `:=`(qaly_healthy = ly_healthy * u_H,
              qaly_healthy_disc = ly_healthy_disc * u_H)]
### Sick
dt_Pop[, `:=`(qaly_Sick = ly_Sick * u_S,
              qaly_Sick_disc = ly_Sick_disc * u_S)]

### Total
dt_Pop[, `:=`(qaly_total = qaly_healthy + qaly_Sick,
              qaly_total_disc = qaly_healthy_disc + qaly_Sick_disc)]

# 06 Compute costs ----
### Healthy
dt_Pop[, `:=`(cost_healthy = ly_healthy * c_H,
              cost_healthy_disc = ly_healthy_disc_costs * c_H)]
### Disease
dt_Pop[, `:=`(cost_Sick = ly_Sick * c_S,
              cost_Sick_disc = ly_Sick_disc * c_S)]

### Total
dt_Pop[, `:=`(cost_total = cost_healthy + cost_Sick,
              cost_total_disc = cost_healthy_disc + cost_Sick_disc)]
comp_time = Sys.time() - t_init
comp_time
# 07 Summary CE outcomes ----
dt_Pop_ce <- dt_Pop[, .(Exp_LY = mean(ly_total),
                        Exp_LY_disc = mean(ly_total_disc),
                        Exp_QALY = mean(qaly_total),
                        Exp_QALY_disc = mean(qaly_total_disc),
                        Exp_COST = mean(cost_total),
                        Exp_COST_disc = mean(cost_total_disc))]
dt_Pop_ce

## By sex ----
dt_Pop_ce_sex <- dt_Pop[, .(Exp_LY = mean(ly_total),
                            Exp_LY_disc = mean(ly_total_disc),
                            Exp_QALY = mean(qaly_total),
                            Exp_QALY_disc = mean(qaly_total_disc),
                            Exp_COST = mean(cost_total),
                            Exp_COST_disc = mean(cost_total_disc)),
                        by = .(Sex)]
dt_Pop_ce_sex
