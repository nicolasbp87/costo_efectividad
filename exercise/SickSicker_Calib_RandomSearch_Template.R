###########################################################################
# This code was created by the DARTH workgroup (www.darthworkgroup.com). 
# When using or modifying this code, please do so with attribution and 
# cite our publications:

# - Alarid-Escudero F, Maclehose RF, Peralta Y, Kuntz KM, Enns EA. 
#   Non-identifiability in model calibration and implications for 
#   medical decision making. Med Decis Making. 2018; 38(7):810-821.

# - Jalal H, Pechlivanoglou P, Krijkamp E, Alarid-Escudero F, Enns E, 
#   Hunink MG. An Overview of R in Health Decision Sciences. 
#   Med Decis Making. 2017; 37(3): 735-746. 
###########################################################################


###################  Calibration Specifications  ###################

# Model: 4-State Sick-Sicker Markov Model
# Inputs to be calibrated: p.S1S2, hr.S1, hr.S2
# Targets: Surv, Prev, PropSick

# Search method: Random search using Latin-Hypercube Sampling
# Goodness-of-fit measure: Sum of Log-Likelihood

####################################################################

# ******************************************************************************
# 01 Calibration Overview ------------------------------------------------------
# ******************************************************************************

### 01.01 Model description  ---------------------------------------------------
# Model: Sick-Sicker 4-state Markov Model
# Inputs to be calibrated: p_S1S2, hr_S1, hr_S2
# Targets: Surv, Prev, PropSick

### 01.02 Calibration method  --------------------------------------------------
# Search method: Random search using Latin-Hypercube Sampling
# Goodness-of-fit measure: Sum of Log-Likelihood

# ******************************************************************************
# 02 Setup ---------------------------------------------------------------------
# ******************************************************************************

### 02.01 Clear environment  ---------------------------------------------------
rm(list = ls())

### 02.02 Load packages  -------------------------------------------------------
# Install pacman if not present
# if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
# Load pacman
# library(pacman)
# Load (install if needed) CRAN packages
# p_load(
#   lhs,         # Latin Hypercube Sampling
#   plotrix,     # Plotting with confidence intervals
#   psych,       # Pairs panels
#   scatterplot3d# 3D scatter plots
# )

# calibration functionality
library(lhs)

# visualization
library(plotrix)
library(psych)
library(scatterplot3d) # now that we have three inputs to estimate, we'll need higher dimension visualization

# ******************************************************************************
# 03 Load calibration targets --------------------------------------------------
# ******************************************************************************
### 03.01 Load target data  ----------------------------------------------------
load("data/SickSicker_CalibTargets.RData")
lst_targets <- SickSicker_targets

### 03.02 Visualize calibration targets  ---------------------------------------
# TARGET 1: Survival ("Surv")
plotrix::plotCI(x = lst_targets$Surv$time, y = lst_targets$Surv$value,
                ui = lst_targets$Surv$ub,
                li = lst_targets$Surv$lb,
                ylim = c(0, 1),
                xlab = "Time", ylab = "Pr Survive")

# TARGET 2: Prevalence ("Prev")
plotrix::plotCI(x = lst_targets$Prev$time, y = lst_targets$Prev$value,
                ui = lst_targets$Prev$ub,
                li = lst_targets$Prev$lb,
                ylim = c(0, 1),
                xlab = "Time", ylab = "Prev")

# TARGET 3: Proportion who are Sick ("PropSick"), among all those afflicted (Sick+Sicker)
plotrix::plotCI(x = lst_targets$PropSick$time, y = lst_targets$PropSick$value,
                ui = lst_targets$PropSick$ub,
                li = lst_targets$PropSick$lb,
                ylim = c(0, 1),
                xlab = "Time", ylab = "PropSick")

# ******************************************************************************
# 04 Load model as a function --------------------------------------------------
# ******************************************************************************

### 04.01 Source model function  -----------------------------------------------
# Function inputs: parameters to be estimated through calibration
# Function outputs: model predictions corresponding to target data
source("code/Functions/SickSicker_MarkovModel_Function.R") # creates the function run_sick_sicker_markov()

### 04.02 Test model function  -------------------------------------------------
v_params_test <- c(p_S1S2 = 0.105, hr_S1 = 3, hr_S2 = 10)
run_sick_sicker_markov(v_params_test) # Test: function works correctly

# ******************************************************************************
# 05 Calibration specifications ------------------------------------------------
# ******************************************************************************

### 05.01 Set random seed  -----------------------------------------------------
set.seed(072218) # For reproducible sequence of random numbers

### 05.02 Define calibration parameters  ---------------------------------------
# number of random samples
n_samp <- 1000

# names and number of input parameters to be calibrated
v_param_names <- c("p_S1S2", "hr_S1", "hr_S2")
n_param <- length(v_param_names)

# Search space bounds
lb <- c(p_S1S2 = 0.01, hr_S1 = 1.0, hr_S2 = 5) # lower bound
ub <- c(p_S1S2 = 0.50, hr_S1 = 4.5, hr_S2 = 15) # upper bound

### 05.03 Define calibration targets  ------------------------------------------
v_target_names <- c("Surv", "Prev", "PropSick")
n_target <- length(v_target_names)

# ******************************************************************************
# 06 Run calibration using Random Search (LHS) ---------------------------------
# ******************************************************************************

### 06.01 Record start time  ---------------------------------------------------
t_init <- Sys.time()

### 06.02 Sample input values using Latin Hypercube Sampling  ------------------

# Sample unit Latin Hypercube
m_lhs_unit <- randomLHS(n_samp, n_param)

# Rescale to min/max of each parameter
m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
for (i in 1:n_param) {
  m_param_samp[, i] <- qunif(m_lhs_unit[, i],
                             min = lb[i],
                             max = ub[i])
}
colnames(m_param_samp) <- v_param_names

### 06.03 Visualize parameter samples  -----------------------------------------
pairs.panels(m_param_samp)

### 06.04 Run the model and calculate goodness-of-fit for each sample  ---------

# initialize goodness-of-fit matrix
m_GOF <- matrix(nrow = n_samp, ncol = n_target)
colnames(m_GOF) <- paste0(v_target_names, "_fit")

# loop through sampled sets of input values
for (j in 1:n_samp) {
  
  # Run model for a given parameter set
  # model_res <- # ADD YOUR CODE HERE
  
  # Calculate goodness-of-fit of model outputs to targets
  
  # TARGET 1: Survival ("Surv") - log likelihood
  # m_GOF[j, 1] <- # ADD YOUR CODE HERE
  
  # TARGET 2: "Prev" - log likelihood
  # m_GOF[j, 2] <- # ADD YOUR CODE HERE
  
  # TARGET 3: "PropSick" - log likelihood
  # m_GOF[j, 3] <- # ADD YOUR CODE HERE
  
  
} # End loop over sampled parameter sets


### 06.05 Calculate overall GOF  -----------------------------------------------
# can give different targets different weights
v_weights <- matrix(1, nrow = n_target, ncol = 1)
# matrix multiplication to calculate weight sum of each GOF matrix row
v_GOF_overall <- c(m_GOF %*% v_weights)
# Store in GOF matrix with column name "Overall"
m_GOF <- cbind(m_GOF, Overall_fit = v_GOF_overall)

### 06.06 Calculate computation time  ------------------------------------------
comp_time <- Sys.time() - t_init
comp_time

# ******************************************************************************
# 07 Explore calibration results -----------------------------------------------
# ******************************************************************************

### 07.01 Identify best-fitting parameter sets  --------------------------------
# Arrange parameter sets in order of fit
m_calib_res <- cbind(m_param_samp, m_GOF)
m_calib_res <- m_calib_res[order(-m_calib_res[, "Overall_fit"]), ]

# Examine the top 10 best-fitting sets
m_calib_res[1:10, ]

### 07.02 Visualize best-fitting parameter sets in 3D  -------------------------
# Plot the top 100 (top 10%)
# Note: When visualizing the top 100 (top 10%) input sets, use scatterplot3d() 
# instead of plot().

# ADD YOUR CODE HERE

### 07.03 Pairwise plots of top parameter sets  --------------------------------
pairs.panels(m_calib_res[1:100, v_param_names])

### 07.04 Compare model predictions to targets  --------------------------------
# Compute output from best and worst parameter set
# v_out_best  <- # ADD YOUR CODE HERE
# v_out_worst <- # ADD YOUR CODE HERE

# Plot: TARGET 1: Survival ("Surv")
# ADD YOUR CODE HERE

# Plot: TARGET 2: Prevalence ("Prev")
# ADD YOUR CODE HERE

# Plot: TARGET 3: PropSick
# ADD YOUR CODE HERE
