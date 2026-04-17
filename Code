#####################################################################################################################

# Analysis of CARTaGENE’s Representativeness - A Case Study on Smoking and Lung Cancer Risk

#####################################################################################################################

# Karim Maher

# This code is part of the study: Analysis of CARTaGENE’s Representativeness - A Case Study on Smoking and Lung Cancer Risk

# by Karim Maher, Vikki Ho, Rado Ramasy, Laura Pelland-St-Pierre, Bouchra Nasri

############ Libraries ############
library(tidyverse)
library(mice)
library(survey)
library(mitools)
library(survival)

############################################################
# 1. Comprehensive Smoking Index (CSI)
############################################################

# Parameters
delta <- 1     # latency parameter (years)
tau   <- 25    # half-life of smoking effect (years)

# Smoking components

# intensity:
# mean number of cigarettes smoked per day
intensity <- cigarettes_per_day

# duration:
# number of years between smoking initiation and cessation
duration <- age_at_stop - age_at_start

# time since cessation:
# number of years between cessation and reference date
time_since_cessation <- age_at_reference - age_at_stop

dur_star <- duration + delta
tsc_star <- time_since_cessation + delta

CSI <- log(1 + intensity * dur_star * exp(-tsc_star / tau))

############################################################
# 2. Absolute Bias per Category (ABc)
############################################################

# p_sample = proportion of category in cohort
# p_population = proportion in reference population

ABc <- abs(p_sample - p_population)

# calculated separately for each category of each variable

############################################################
# 3. Conditional Model before inputation and weights 
############################################################

model_unweighted <- clogit(
  outcome ~ CSI + covariates + strata(risk_set_id),
  data = CaG
)

############################################################
# 4. Multiple Imputation
############################################################

method <- make.method(data)

method["continuous_var"]  <- "pmm"
method["binary_var"]      <- "logreg"
method["categorical_var"] <- "polyreg"

pred <- make.predictorMatrix(data)

pred[,] <- 0
pred["var_to_impute", c("predictor1","predictor2","predictor3")] <- 1

imp <- mice(
  data,
  m = 50,
  method = method,
  predictorMatrix = pred,
  seed = 123
)

############################################################
# 5. Conditional Model after Imputation (Unweighted)
############################################################
#Creation of Imputed Dataset
imputed_datasets <- lapply(1:50, function(i) {
  CaG_i <- complete(imp, i)
})

# Fit the model separately in each imputed dataset
models <- lapply(imputed_datasets, function(data) {
  data$outcome <- as.numeric(as.character(data$outcome))
  clogit(
    outcome ~ CSI + covariates + strata(risk_set_id),
    data = data
  )
})

mira_obj <- as.mira(models)
pooled <- pool(mira_obj)
summary(pooled)

############################################################
# 6. Raking Weights
############################################################

# create survey design for each imputed dataset
design_imp <- lapply(complete(imp, "all"), function(d){
  svydesign(
    ids = ~1,
    data = d
  )
})

# apply raking calibration
svy_rake <- lapply(design_imp, function(des){
  rake(
    design = des,
    sample.margins = list(~varA, ~varB, ~varC),
    population.margins = list(popA, popB, popC)
  )
})

# varA, varB, varC: sociodemographic variables used to calibrate the sample

# popA, popB, popC: marginal distributions of these variables in the target population

# different raking specifications can be tested depending on which variables are included in the calibration

############################################################
# 7. Conditional Model after Raking (Weighted)
############################################################

res_raked <- lapply(svy_rake, function(des){
  svycoxph(
    Surv(time_variable, outcome) ~
      CSI + covariates + strata(risk_set_id),
    design = des
  )
})
