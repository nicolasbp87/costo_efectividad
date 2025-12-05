# title: "Sensitivity Analysis using Linear Regression Metamodeling"
# authors: "Fernando Alarid-Escudero & Jeremy Goldhaber-Fiebert"
# date: 2024/02/16

# Initial setup ----

## Install Packages (If you haven't done so yet!)
# install.packages("stargazer")
# install.packages("xtable")
# install.packages("gdata")
# install.packages("xlsx")
# install.packages("ggplot2")
# install.packages("reshape2")
# install.packages("reshape")
# install.packages("VGAM")
# install.packages("matrixStats")
# install.packages("ellipse")
library(dplyr) # For data.frame wrangling and manipulation
library(dampack) # For CEA and CEA plotting
library(ggplot2) # For plotting

## Clean everything from workspace
rm(list = ls())

#Load the plotting functions
source("code/MM_functions.R")

## Load PSA dataset
# Read the `.csv` simulation file into `R`.
#Sim <- read.csv("PSA_100k.csv",header=TRUE) # Feel free to try a larger PSA!
Sim <- read.csv("data/PSA.csv", header = TRUE)

#Display first five observations of the data fram using the command `head`
head(Sim)[, 1:7] #Only display the 6 columns of the simulation file

# Format PSA objects ----
# Label our strategies
v_names_str <- c("Chemo", "Radio", "Surgery")
v_names_str

# Obtain the number of strategies
n_str <- length(v_names_str)
n_str

#* Create data.frames with PSA output 
#* data.frame of costs
df_c <- Sim %>%
  select(Chemo_Cost, Radio_Cost, Surg_Cost)
colnames(df_c) <- v_names_str
#* data.frame of effectiveness
df_e <- Sim %>%
  select(Chemo_Eff, Radio_Eff, Surg_Eff)
colnames(df_e) <- v_names_str
#* data.frame of input parameters
df_psa_input <- Sim %>%
  select(pFailChemo, pFailRadio, pFailSurg, pDieSurg, 
         muDieCancer, 
         cChemo, cRadio, cSurg)
v_names_params <- colnames(df_psa_input)
v_names_params
n_params <- length(v_names_params)
n_params

## Create PSA object ----
#* `make_psa_obj` function included in `dampack` package
l_psa <- make_psa_obj(cost          = df_c, 
                      effectiveness = df_e, 
                      parameters    = df_psa_input, 
                      strategies    = v_names_str)

## Create matrices for outcomes and parameters
Outcomes <- Sim[, 2:(n_str*2 + 1)]
head(Outcomes)
Parms <- Sim[, (n_str*2 + 2):(ncol(Sim))]
head(Parms)


## Compute Net Health Benefit (NHB) ----

#* Define a willingness to pay (WTP), $\lambda=\$50,000$
lambda <- 50000

#* Create data.frame with NHB per strategy
df_nhb <- df_e - df_c/lambda
head(df_nhb)

#* Combine data.frames for metamodeling analysis
df_metamodel <- bind_cols(df_nhb, df_psa_input)
# Calculate NHB for each strategy and add it to the original data frame
# Sim$Chemo_NHB <- Sim$Chemo_Eff - Sim$Chemo_Cost/lambda
# Sim$Radio_NHB <- Sim$Radio_Eff - Sim$Radio_Cost/lambda
# Sim$Surg_NHB  <- Sim$Surg_Eff  - Sim$Surg_Cost/lambda

#Create a NHB matrix
#NHB <- Sim[,((dim(Sim)[2]-(dep-1)):dim(Sim)[2])]
# NHB <- Sim[,((ncol(Sim)-(n_str-1)):ncol(Sim))]
# head(NHB)

# Linear Regression Metamodel ----
## Linear regression metamodel on NHB
#Linear reression model for each strategy's NHB using the function `lm`
ChemoNHB.fit <- lm(Chemo ~ pFailChemo + pFailRadio + pFailSurg + pDieSurg + 
                     muDieCancer + cChemo + cRadio + cSurg, 
                   data = df_metamodel)
summary(ChemoNHB.fit) #Display results
RadioNHB.fit <- lm(Radio ~ pFailChemo + pFailRadio + pFailSurg + pDieSurg + 
                     muDieCancer + cChemo + cRadio + cSurg, 
                   data = df_metamodel)
summary(RadioNHB.fit) #Display results
SurgNHB.fit <- lm(Surgery ~ pFailChemo + pFailRadio + pFailSurg + pDieSurg + 
                      muDieCancer + cChemo + cRadio + cSurg, 
                    data = df_metamodel)
summary(SurgNHB.fit) #Display results

## Display results with publication quality
# We will use the `stargazer` package
library(stargazer) #More info: http://www.jakeruss.com/a-stargazer-cheatsheet.html

#Linear regresion metamodel results
  stargazer(ChemoNHB.fit, type = 'text', #Change type="latex" for latex output or "html" for html
            dep.var.caption = " NHB ", dep.var.labels = c(" Chemo "), 
            notes = c('NHB: Net health benefit  '),
            digits = 2, digits.extra = 1, keep.stat = c("rsq", "n"))

#Lets interpret the reuslts!
# Magnitude
# Significance (i.e., p-value)

## Standardized parameters
# We can standardized the simulation parameters within the `lm` function using the command `scale`.
ChemoNHB.fit.std <- lm(Chemo ~ scale(pFailChemo) + scale(pFailRadio) + 
                         scale(pFailSurg) + scale(pDieSurg) + 
                         scale(muDieCancer) + scale(cChemo) + scale(cRadio) + 
                         scale(cSurg), data = df_metamodel)

## Standardized parameters {.smaller .vcenter .flexbox}
stargazer(ChemoNHB.fit.std, type = 'text', 
          dep.var.caption = " NHB ", dep.var.labels = c(" Chemo "),
          covariate.labels = v_names_params,
          notes = c('NHB: Net health benefit  '),
          digits = 2, digits.extra = 1, keep.stat = c("rsq","n")) 

## Linear regression metamodel for all strategies
RadioNHB.fit.std <- lm(Radio ~ scale(pFailChemo) + scale(pFailRadio) + 
                         scale(pFailSurg) + scale(pDieSurg) + 
                         scale(muDieCancer) + scale(cChemo) + scale(cRadio) + 
                         scale(cSurg), data = df_metamodel)
SurgNHB.fit.std <- lm(Surgery ~ scale(pFailChemo) + scale(pFailRadio) + 
                        scale(pFailSurg) + scale(pDieSurg) + 
                        scale(muDieCancer) + scale(cChemo) + scale(cRadio) + 
                        scale(cSurg), data = df_metamodel)

## Display results for all strategies {.smaller .vcenter .flexbox}
#Regression coefficients from metamodels on the NHB of all strategies
stargazer(ChemoNHB.fit.std, RadioNHB.fit.std, SurgNHB.fit.std, type = 'text',
          covariate.labels = v_names_params, dep.var.caption = " NHB", 
          dep.var.labels = v_names_str,
          digits = 2, digits.extra = 1, keep.stat = c("rsq","n"),
          notes = c('NHB: Net health benefit  '), notes.align = 'l',
          omit.table.layout = "#")

# One-Way Sensitivity Analysis (OWSA) ---

## Define parameter and range
# Parameter and its range to do one-way sensitivity analysis. 
parm <- 'muDieCancer'

## OWSA using `base` R ----
v_ranges <- c(0.036, 0.367)

# Create a vector with 400 values between the min and max `v_ranges` values
parm.range <- seq(v_ranges[1], v_ranges[2], length.out = 400)

## Multiple Multivariate Regression (MMR) Metamodel
# Run a linear regression for the NHB of each strategy
# We will estimate a MMR metamodel using `lm` again, but with multiple dependent variables to estimate all linear regression simultaneously
# Tell `lm` which are the multiple dependent variables by creating a matrix of dependent variables using `cbind` command
# To account the nonlinearities we will include a second order polynomial using the `poly`command.
# Note: Second order polynomial makes sense based on Talyor series expansion.
Oway.mlm <- lm(cbind(Chemo, Radio, Surgery) ~ poly(muDieCancer, degree = 2)+
                 pFailChemo + pFailRadio + pFailSurg + pDieSurg + 
                 cChemo + cRadio + cSurg,  data = df_metamodel)
summary(Oway.mlm)

## Generate predicted values
# Generate matrix as a base for prediction and name the columns as the parameters
Oway.mat <- matrix(rep(colMeans(df_psa_input)), 
                   nrow = length(parm.range), ncol = n_params, 
                   byrow = T)
colnames(Oway.mat) <- v_names_params

# Substitute the column corresponding to the parameter of interest with the 
# vector defined by the ranges
Oway.mat[, parm] <- parm.range

#Transform to data frame, the format required for `predict` command
Oway.mat <- data.frame(Oway.mat)

## Generate predicted values
# Predict outcomes (i.e., NHB) for each strategy on the values of the parameters 
# defined previously using the MMMR Metamodel 
Oway.pred <- data.frame(predict(Oway.mlm, newdata = Oway.mat))

#Name the predicted outcome columns with the names of the strategies
colnames(Oway.pred) <- v_names_str
head(Oway.pred)

## Plot the results
#We will plot the graphs of this course using the package `ggplot2`
# More info on the package: (http://link.springer.com/book/10.1007%2F978-0-387-98141-3). 
# To use `ggplot` our data frame needs to be modified into a more suitable form. 
# Use the command `stack` to concatenate multiple vectors into a single vector 
# (i.e., column)
Oway.pred <- stack(Oway.pred, select = v_names_str) #Strategies define which variables to select from the data frame
colnames(Oway.pred) <- c("NHB", "Strategy")
head(Oway.pred)

## Ploting the results 
# Combine the predicted outcomes with the parameter values
Oway.plot <- cbind(Oway.pred, Oway.mat)
head(Oway.plot)

# Use `ggplot` function from the `ggplot2` package to plot the results
ggplot(data = Oway.plot, aes(x = muDieCancer, y = NHB, color = Strategy)) +
  geom_line()

## Ploting the results {.smaller}
# Lets do some few editing
txtsize <- 16 #Text size for the graphs
ggplot(data = Oway.plot, aes(x = muDieCancer, y = NHB, color = Strategy)) +
  geom_line() +
  ggtitle("One-way sensitivity analysis \n Net Health Benefit") + 
  xlab(parm) +
  scale_colour_hue("Strategy", l = 50) +
  scale_x_continuous(n.breaks = 7) + #Adjust for number of ticks in x axis
  scale_y_continuous("E[NHB]", n.breaks = 7) +
  theme_bw(base_size = txtsize) +
  theme(legend.position = "bottom",
        legend.key = element_rect(colour = "black"))
# OneWaySA(v_names_str, Parms, df_nhb, parm, v_ranges) 

## OWSA using `dampack` ----
#* Create `owsa` object with function included in `dampack` package
# For one parameter 
l_owsa_one_param <- owsa(sa_obj = l_psa, params = parm, 
                         outcome = "nhb", wtp = lambda)
gg_owsa_one_param <- plot(l_owsa_one_param) + 
  scale_y_continuous("E[NHB]", n.breaks = 8) +
  theme(legend.position = "bottom")
gg_owsa_one_param

# For all parameters
l_owsa <- owsa(sa_obj = l_psa, outcome = "nhb", wtp = lambda)
gg_owsa <- plot(l_owsa) + 
  scale_y_continuous("E[NHB]", n.breaks = 8) +
  theme(legend.position = "bottom")
gg_owsa

# Threshold Analysis ----

## Threshold Analysis using base R ----
## Break even values
# Mean INHB. E.i., difference between intercepts
theta0ChemoRad = ChemoNHB.fit.std$coef[1] - RadioNHB.fit.std$coef[1]
theta0ChemoSurg = ChemoNHB.fit.std$coef[1] - SurgNHB.fit.std$coef[1]
theta0RadSurg = RadioNHB.fit.std$coef[1] - SurgNHB.fit.std$coef[1]
# Define vector for thresholds (useful when looping)
theta.ChemoRad  <- rep(0, n_params)
theta.ChemoSurg <- rep(0, n_params)
theta.RadSurg   <- rep(0, n_params)
# Coefficients corresponding to each parameters
for (i in 1:(n_params)) {
  theta.ChemoRad[i] <- ChemoNHB.fit.std$coef[i + 1] - RadioNHB.fit.std$coef[i + 1]
  theta.ChemoSurg[i] <- ChemoNHB.fit.std$coef[i + 1] - SurgNHB.fit.std$coef[i + 1]
  theta.RadSurg[i] <- RadioNHB.fit.std$coef[i + 1] - SurgNHB.fit.std$coef[i + 1]
}

z.ChemoRad <- -theta0ChemoRad/theta.ChemoRad
z.ChemoRad[c(3, 4, 8)] <- NA
z.ChemoSurg <- -theta0ChemoSurg/theta.ChemoSurg
z.ChemoSurg[c(2, 7)] <- NA
z.RadSurg <- -theta0RadSurg/theta.RadSurg
z.RadSurg[c(1, 7)] <- NA

## Table of Thresholds `stargazer` 
#Matrix with parameter names, and threshold for each INHB
threshold <- cbind(v_names_params,
                   round(z.ChemoRad, 2), 
                   round(z.ChemoSurg, 2), 
                   round(z.RadSurg, 2))
colnames(threshold) <- c('Parameter', '[Chemo>Radio]', '[Chemo>Surg]', 
                         '[Radio>Surgery]')

#use stargazer to display results
stargazer(threshold, type = 'text', title = 'One-Way Threshold Analysis',
          digits = 1, digits.extra = 2)

# Two-Way Sensitivity Analysis (TWSA) ----

# Parameters
parm1 <- 'muDieCancer'#5
parm2 <- 'pFailChemo'#1

## OWSA using `base` R ----
## Define parameter ranges
range1 <- c(0.036, 0.367)
range2 <- c(0.36, 0.558)

# Create a couple of vectors with 301 values between the min and max `range` values
parm1.range <- seq(from = range1[1], to = range1[2], length.out = 301)
parm2.range <- seq(from = range2[1], to = range2[2], length.out = 301)

## Multiple Multivariate Regression (MMR) Metamodel
#Run a linear regression for the NHB of each strategy
#We will estimate a MMR metamodel using `lm` again, but with multiple dependent variables to estimate all linear regression simultaneously
#Tell `lm` which are the multiple dependent variables by creating a matrix of dependent variables using `cbind` command
#To account the nonlinearities we will include a second order polynomial using the `poly`command.
  #Note: Second order polynomial makes sense based on Talyor series expansion.
Tway.mlm <- lm(cbind(Chemo, Radio, Surgery) ~ 
                 poly(muDieCancer,2) + poly(pFailChemo,2) +
                 muDieCancer:pFailChemo + pFailRadio + pFailSurg + pDieSurg +
                 cChemo + cRadio + cSurg, data = df_metamodel)
summary(Tway.mlm)

## Generate predicted values
#Generate an expanded grid with all teh combination of both parameters using `expand.grid`
TWSA <- expand.grid(parm1 = parm1.range, parm2 = parm2.range)
#Generate matrix as a base for prediction and name the columns as the parameters
Tway.mat <- matrix(rep(colMeans(Parms)), 
                   nrow = nrow(TWSA), ncol = ncol(Parms), 
                   byrow = T)
colnames(Tway.mat) <- colnames(Parms)
#Substitute the column corresponding to the parameter of interest with the vector defined by the ranges
Tway.mat[, parm1] <- TWSA[, 1]
Tway.mat[, parm2] <- TWSA[, 2]

## Generate predicted values
#Transform to data frame, the format required for `predict` command
Tway.mat <- data.frame(Tway.mat)
head(Tway.mat)

## Generate predicted values
#Predict outcomes (i.e., NHB) for each strategy on the values of the parameters defined previously using the MMMR Metamodel 
Tway.pred <- data.frame(predict(Tway.mlm, newdata = Tway.mat))

#Name the predicted outcome columns with the names of the strategies
colnames(Tway.pred) <- v_names_str
head(Tway.pred)

## Determine optimal strategy
#Find optimal strategy in terms of maximum Outcome
Optimal <- max.col(Tway.pred)

#Get outcome's name of optimal strategy
OptimalOut <- apply(Tway.pred, 1, max)

## Plotting the results
#Add the the previous vectors to our data frame
Tway.pred$Strategy <- factor(Optimal, labels = v_names_str)
Tway.pred$value <- OptimalOut

## Plotting the results 
#Combine the predicted outcomes with the parameter values into a new data frame.
Tway.plot <- cbind(Tway.pred, Tway.mat)
head(Tway.plot)

## Plotting the results {.smaller}
#Use the `geom_tile` option of `ggplot` to create a surface plot and fill it accordingly to different `Strategy`.
ggplot(Tway.plot, aes(x = muDieCancer, y = pFailChemo)) + 
  geom_tile(aes(fill = Strategy))

## Plotting the results {.smaller}
#After a few editing
ggplot(Tway.plot, aes(x = muDieCancer, y = pFailChemo)) + 
  geom_tile(aes(fill = Strategy)) +
  theme_bw() +
  ggtitle(expression(atop("Two-way sensitivity analysis", 
  atop("Net Health Benefit")))) + 
  scale_fill_discrete("Strategy: ", l = 50) +
  xlab(parm1) +
  ylab(parm2) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = txtsize),
  legend.key = element_rect(colour = "black"),
  legend.text = element_text(size = txtsize),
  title = element_text(face = "bold", size = 15),
  axis.title.x = element_text(face = "bold", size = txtsize),
  axis.title.y = element_text(face = "bold", size = txtsize),
  axis.text.y = element_text(size = txtsize),
  axis.text.x = element_text(size = txtsize))

#The two-way sensitivity analysis of the NHB on **`r parm1`** and **`r parm2`** looks like
TwoWaySA(Strategies = v_names_str, Parms = Parms, Outcomes = df_nhb,
         parm1 = parm1, parm2 = parm2, range1 = range1, range2 = range2)

## TWSA using `dampack` ----
#* Create `twsa` object with function included in `dampack` package
# For one parameter 
l_twsa <- twsa(sa_obj = l_psa, param1 = parm1, param2 = parm2, poly.order = 2,
                         outcome = "nhb", wtp = lambda, nsamp = 300)
gg_twsa <- plot(l_twsa) + 
  theme(legend.position = "bottom")
gg_twsa

### Decision-Sensitivity Analysis | Multinomial Logistic (MNL) Metamodel ####

## Define parameter and range to do MNL decision-sensitivity analysis. 
#Name the parameter
parm.mnl <- 'muDieCancer'
#The range will be determined by a 95% percentile coverage of the PSA sample
range.mnl <- seq(2.5, 97.5, length = 400) #vector to define 400 samples between the 2.5th and 97.5th percentiles
j <- round(range.mnl*(length(Parms[, parm.mnl])/100)) #indexing vector;j=round(y*n/100) where n is the size of vector of interest
parm.range.mnl <- sort(Parms[j, parm.mnl])

#Create a vector with 400 values between the min and max `range` values
#parm.range.mnl<-seq(range[1],range[2],length.out=400)

#Calculate the preferred strategy (i.e., optimal startegy) in terms of its cost-efectiveness for each of the simulations.
Optimal.mnl <- data.frame(max.col(df_nhb))
names(Optimal.mnl) <- "Strategy"

## Multinomial Logistic (MNL) Metamodel
#Merge the optimal strategy indicator with the rest of the parameters
df_mnl <- data.frame(Optimal.mnl, Parms)
#Estimate the MNL model over the Strategy variable using the function `vglm` from the `VGAM`package
library(VGAM)
MNL.fit <- vglm(Strategy ~ ., 
                family = multinomial,
                data = df_mnl)
summary(MNL.fit)

## Predicted values
# Generate matrix as a base for prediction and name the columns as the parameters
m_mnl <- matrix(rep(colMeans(Parms)), 
                nrow = length(parm.range.mnl), ncol = ncol(Parms), 
                byrow = T)
colnames(m_mnl) <- colnames(Parms)

#Substitute the column corresponding to the parameter of interest with the vector defined by the ranges
m_mnl[, parm] <- parm.range.mnl

#Transform to data frame, the format required for `predict` command
df_mnl <- data.frame(m_mnl)

## Generate predicted values
# Predict outcomes (i.e., probability of each strategy being optimal)on the values of the parameters defined previously using the MNL Metamodel 
MNL.pred <- data.frame(predict(MNL.fit, newdata = df_mnl, type = "response"))
#Name the predicted outcome columns with the names of the strategies  
colnames(MNL.pred) <- v_names_str 
head(MNL.pred)

# Concatenate multiple vectors into a single vector along with a factor indicating where each observation originated. 
MNL.pred <- stack(MNL.pred, select = v_names_str) #Strategies define which variables to select from the data frame
colnames(MNL.pred) <- c("Probability", "Strategy")
head(MNL.pred)
MNL.plot <- cbind(MNL.pred, df_mnl) 
head(MNL.plot)


## Plotting the results {.smaller}
#Use `ggplot` to plot the results
ggplot(data = MNL.plot, 
       aes(x = muDieCancer, y = Probability, color = Strategy)) +
  geom_point(size = 2) +
  geom_line()

## Plotting the results {.smaller}
# Lets do some few editing
ggplot(data = MNL.plot, 
       aes(x = muDieCancer, y = Probability, color = Strategy)) +
  geom_point(size = 2) +
  geom_line() +
  ggtitle("One-way sensitivity analysis \n Net Health Benefit") + 
  xlab(parm) +
  ylab("E[NHB]") +
  scale_colour_hue("Strategy", l=50) +
  scale_x_continuous() + #Adjust for number of ticks in x axis
  scale_y_continuous() +
  theme_bw() +
  theme(legend.position="bottom",legend.title=element_text(size = txtsize),
        legend.key = element_rect(colour = "black"),
        legend.text = element_text(size = txtsize),
        title = element_text(face="bold", size=15),
        axis.title.x = element_text(face="bold", size=txtsize),
        axis.title.y = element_text(face="bold", size=txtsize),
        axis.text.y = element_text(size=txtsize),
        axis.text.x = element_text(size=txtsize))

## Plotting the results {.vcenter .flexbox}
# Using the user-built function
MNL.SA(parm, range, v_names_str, Parms, df_nhb)

# Generate CE figures ----
## Visualize PSA results and CEA ----
#* Vector with willingness-to-pay (WTP) thresholds.
v_wtp <- seq(0, 150000, by = 5000)

### Cost-Effectiveness Scatter plot ----
#* `plot.psa` function included in `dampack` packages.
#* Plot with colorblind-friendly colors
gg_scattter <- plot(l_psa, txtsize = txtsize) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  scale_y_continuous("Cost (Thousand $)", 
                     n.breaks = 10,
                     labels = function(x) x/1000) +
  scale_x_continuous("Effectiveness (QALYs)",
                     n.breaks = 10) +
  guides(col = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom")
gg_scattter

### Incremental cost-effectiveness ratios (ICERs) with probabilistic output ----
#* Compute expected costs and effects for each strategy from the PSA
#* `summary.psa` function included in `dampack` package
df_out_ce_psa <- summary(l_psa)

#* Function included in "R/Functions.R"; depends on the `dplyr` package
#* The latest version can be found in `dampack` package
df_cea_psa <- calculate_icers(cost       = df_out_ce_psa$meanCost, 
                              effect     = df_out_ce_psa$meanEffect,
                              strategies = df_out_ce_psa$Strategy)
df_cea_psa

### Plot cost-effectiveness frontier with probabilistic output ----
#* `plot.icers` function included in `dampack` package
plot(df_cea_psa, label = "all", txtsize = txtsize) +
  # expand_limits(x = max(table_cea$QALYs) + 0.1) +
  theme(legend.position = c(0.7, 0.2))

## Cost-effectiveness acceptability curves (CEACs) and frontier (CEAF) ---
#* `ceac`, and `plot.ceac` functions included in `dampack` package
ceac_obj <- ceac(wtp = v_wtp, psa = l_psa)
#* Regions of highest probability of cost-effectiveness for each strategy
summary(ceac_obj)
#* CEAC & CEAF plot
gg_ceac <- plot(ceac_obj, txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  scale_y_continuous("Pr Most Cost-Effective", 
                     limits = c(0, 1),
                     n.breaks = 10) +
  theme(legend.position = c(0.8, 0.8))
gg_ceac

## Expected Loss Curves (ELCs) ----
#* `calc_exp_loss` and `plot.exp_loss` functions included in `dampack` package
elc_obj <- calc_exp_loss(wtp = v_wtp, psa = l_psa)
elc_obj

#* ELC plot
gg_elc <- plot(elc_obj, log_y = FALSE, 
               txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14,
               col = "full") +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  # geom_point(aes(shape = as.name("Strategy"))) +
  scale_y_continuous("Expected Loss (Thousand $)", 
                     n.breaks = 20,
                     labels = function(x) x/1000) +
  theme(legend.position = c(0.3, 0.7),)
gg_elc

## Expected value of perfect information (EVPI) ----
#* `calc_exp_loss` function included in `dampack` package
evpi <- calc_evpi(wtp = v_wtp, psa = l_psa)
#* EVPI plot
gg_evpi <- plot(evpi, effect_units = "QALY", 
                txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14) +
  scale_y_continuous("EVPI (Thousand $)", 
                     n.breaks = 20,
                     labels = function(x) x/1000)
gg_evpi
