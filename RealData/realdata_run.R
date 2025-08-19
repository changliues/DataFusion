source("DataFusion/minimax_core.R")
source("DataFusion/datafusion_dataset.R")
source("DataFusion/datafusion_minimax.R")
source("DataFusion/compare_run_para_functions.R")
source("DataFusion/create_df.R")
source("DataFusion/compare_model_specification.R")

library(caret)
library(plyr)

alpha <- 0.05
critical_t <- qnorm(1 - alpha/2)

# functions to be used
expit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}

###############################################read data###############################################
dat <- readRDS(file = "dat.rds")
sample.size <- nrow(dat)
###############################################proxci###############################################
# number of folds
n.flds <- 4
###############################################BSIV###############################################
#################################### test robustness ####################################
proximal_model_spec_urban <- readRDS(file = "proximal_model_spec_urban.rds")
Athey_model_spec <- readRDS(file = "Athey_model_spec.rds")
naive_model_spec <- readRDS(file = "naive_model_spec.rds")
##########################################################################################
# number of folds
k <- 4
##################################################################################################

# Create data
dataset <- ProxCIData$new(X = cbind(dat$X1,dat$X2,dat$B), Z = dat$Z, M = dat$M, M.fac = dat$M.fac, Y = dat$Y, A = dat$A, G = dat$G, Xs = cbind(dat$X1,dat$X2,dat$B))
################--proxci estimation(location proxcy)--################
source("DataFusion/datafusion_proximal_inference_discrete_urban.R")
proxci_estimator <- ProximalInference$new(
  dataset,
  crossfit_folds =n.flds,
  mod_spec = proximal_model_spec_urban)
  
est <- proxci_estimator$est()

result.proxci <- list(est1 = est$estimates1,
                      est2 = est$estimates2,
                      sample.size = sample.size)
################--Athey estimation--################
cf_inds <- dataset$create_crossfit_split(k)
Athey_estimator <- run_para_Athey(dataset,Athey_model_spec,k,cf_inds)
################--naive & equi_confounding estimation--################
cf_inds <- dataset$create_crossfit_split(k)
naive_estimator <- run_para_naive(dataset,naive_model_spec,k,cf_inds)
