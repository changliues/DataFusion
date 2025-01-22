source("DataFusion/minimax_core.R")
source("DataFusion/datafusion_dataset.R")
source("DataFusion/datafusion_minimax.R")
source("DataFusion/datafusion_proximal_inference.R")
source("DataFusion/datafusion_datagen.R")
source("DataFusion/BSIV_run_para_functions.R")
source("DataFusion/BSIV_model_specification.R")

library(caret)
library(plyr)

# sample sizes
sample.size.for.sim = 1000

alpha <- 0.05
critical_t <- qnorm(1 - alpha/2)

# functions to be used
expit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}

###############################################proxci###############################################
# number of folds
split_ratio <- 0.25
n.flds <- 4
# gammas and lambdas sequences
gammash <- c(0.5,0.6,0.7)
lambdash <- 2^seq(-1,1)*log(sample.size.for.sim)/sample.size.for.sim
lambdasfh <- 2^seq(-1,1)*log(sample.size.for.sim)/sample.size.for.sim
gammasq <- c(0.25,0.3,0.35)
lambdasq <- 2^seq(-2,0)*log(sample.size.for.sim)/sample.size.for.sim
lambdasfq <- 2^seq(-2,0)*log(sample.size.for.sim)/sample.size.for.sim
gammash_nt <- 1.5
lambdash_nt <- 4*log(sample.size.for.sim)/sample.size.for.sim
lambdasfh_nt <- 4*log(sample.size.for.sim)/sample.size.for.sim
gammasq_nt <- 1
lambdasq_nt <- 2*log(sample.size.for.sim)/sample.size.for.sim
lambdasfq_nt <- 2*log(sample.size.for.sim)/sample.size.for.sim

params0 <- list()
for(j in 1:n.flds){
  params0[[j]] <- list()
  params0[[j]]$h_params <- as.data.frame(matrix(NA, nrow = 1, ncol = 4))
  names(params0[[j]]$h_params) <- c("lambda_h", "lambda_f", "gamma", "score")
  params0[[j]]$q_params <- as.data.frame(matrix(NA, nrow = 1, ncol = 4))
  names(params0[[j]]$q_params) <- c("lambda_q", "lambda_f", "gamma", "score")
}

params1 <- list()
for(j in 1:n.flds){
  params1[[j]] <- list()
  params1[[j]]$h_params <- as.data.frame(matrix(NA, nrow = 1, ncol = 4))
  names(params1[[j]]$h_params) <- c("lambda_h", "lambda_f", "gamma", "score")
  params1[[j]]$q_params <- as.data.frame(matrix(NA, nrow = 1, ncol = 4))
  names(params1[[j]]$q_params) <- c("lambda_q", "lambda_f", "gamma", "score")
}

params2 <- list()
for(j in 1:n.flds){
  params2[[j]] <- list()
  params2[[j]]$h_params <- as.data.frame(matrix(NA, nrow = 1, ncol = 4))
  names(params2[[j]]$h_params) <- c("lambda_h", "lambda_f", "gamma", "score")
  params2[[j]]$q_params <- as.data.frame(matrix(NA, nrow = 1, ncol = 4))
  names(params2[[j]]$q_params) <- c("lambda_q", "lambda_f", "gamma", "score")
}

params3 <- list()
for(j in 1:n.flds){
  params3[[j]] <- list()
  params3[[j]]$h_params <- as.data.frame(matrix(NA, nrow = 1, ncol = 4))
  names(params3[[j]]$h_params) <- c("lambda_h", "lambda_f", "gamma", "score")
  params3[[j]]$q_params <- as.data.frame(matrix(NA, nrow = 1, ncol = 4))
  names(params3[[j]]$q_params) <- c("lambda_q", "lambda_f", "gamma", "score")
}

params4 <- list()
for(j in 1:n.flds){
  params4[[j]] <- list()
  params4[[j]]$h_params <- as.data.frame(matrix(NA, nrow = 1, ncol = 4))
  names(params4[[j]]$h_params) <- c("lambda_h", "lambda_f", "gamma", "score")
  params4[[j]]$q_params <- as.data.frame(matrix(NA, nrow = 1, ncol = 4))
  names(params4[[j]]$q_params) <- c("lambda_q", "lambda_f", "gamma", "score")
}

###############################################BSIV###############################################
#################################### test robustness ####################################
BSIV_model_spec <- readRDS(file = "BSIV_model_spec.rds")
proximal_model_spec <- readRDS(file = "proximal_model_spec.rds")
##########################################################################################
# number of folds
k <- 4
##################################################################################################

run_sim <- function(n.sim){
  
  # Generate data
  data <- generate_data(sample.size.for.sim)
  dataset <- data$dataset
  target <- data$target
  
  ################--proxci estimation--################
  proxci_estimator <- ProximalInference$new(
    dataset,
    crossfit_folds =n.flds,
    gammash = gammash,
    lambdash = lambdash,
    lambdasfh = lambdasfh,
    gammasq = gammasq,
    lambdasq = lambdasq,
    lambdasfq = lambdasfq,
    gammash_nt = gammash_nt,
    lambdash_nt = lambdash_nt,
    lambdasfh_nt = lambdasfh_nt,
    gammasq_nt = gammasq_nt,
    lambdasq_nt = lambdasq_nt,
    lambdasfq_nt = lambdasfq_nt,
    split_ratio = split_ratio,
    print_best_params = T,
    mod_spec = proximal_model_spec)
  
  ################--proxci estimation (all true)--################
  est0 <- proxci_estimator$est0()
  
  ifb.lower0 <- est0$estimates4 - critical_t * est0$estimates4.sd/sqrt(sample.size.for.sim)
  ifb.upper0 <- est0$estimates4 + critical_t * est0$estimates4.sd/sqrt(sample.size.for.sim)
  
  h0 <- proxci_estimator$h
  q0 <- proxci_estimator$q
  # dens_m <- estimator$dens_m
  
  for(j in 1:n.flds){
    params0[[j]]$h_params <- c(h0[[j]]$lambda_h,h0[[j]]$lambda_f,h0[[j]]$gamma, h0[[j]]$score)
    params0[[j]]$q_params <- c(q0[[j]]$lambda_q,q0[[j]]$lambda_f,q0[[j]]$gamma, q0[[j]]$score)
    # params[[j]]$dens_m_params <- c(dens_m[[j]]$lambda_kde,dens_m[[j]]$mse)
    
  }
  
  result.proxci0 <- list(true.ETT = target$ett.true,
                         est1 = est0$estimates1,
                         est2 = est0$estimates2,
                         est3 = est0$estimates3,
                         ifb.est = est0$estimates4,
                         ifb.sd = est0$estimates4.sd,
                         ifb.lower = ifb.lower0,
                         ifb.upper = ifb.upper0,
                         ifb.cov = target$ett.true >= ifb.lower0 & target$ett.true <= ifb.upper0,
                         sample.size = sample.size.for.sim,
                         params = params0
  )
  
  ################--proxci estimation (case1)--################
  est1 <- proxci_estimator$est1()
  
  ifb.lower1 <- est1$estimates4 - critical_t * est1$estimates4.sd/sqrt(sample.size.for.sim)
  ifb.upper1 <- est1$estimates4 + critical_t * est1$estimates4.sd/sqrt(sample.size.for.sim)
  
  h1 <- proxci_estimator$h
  q1 <- proxci_estimator$q_nt
  # dens_m <- estimator$dens_m
  
  for(j in 1:n.flds){
    params1[[j]]$h_params <- c(h1[[j]]$lambda_h,h1[[j]]$lambda_f,h1[[j]]$gamma, h1[[j]]$score)
    params1[[j]]$q_params <- c(q1[[j]]$lambda_q,q1[[j]]$lambda_f,q1[[j]]$gamma, q1[[j]]$score)
  }
  
  result.proxci1 <- list(true.ETT = target$ett.true,
                         est1 = est1$estimates1,
                         est2 = est1$estimates2,
                         est3 = est1$estimates3,
                         ifb.est = est1$estimates4,
                         ifb.sd = est1$estimates4.sd,
                         ifb.lower = ifb.lower1,
                         ifb.upper = ifb.upper1,
                         ifb.cov = target$ett.true >= ifb.lower1 & target$ett.true <= ifb.upper1,
                         sample.size = sample.size.for.sim,
                         params = params1
  )
  
  ################--proxci estimation (case2)--################
  est2 <- proxci_estimator$est2()
  
  ifb.lower2 <- est2$estimates4 - critical_t * est2$estimates4.sd/sqrt(sample.size.for.sim)
  ifb.upper2 <- est2$estimates4 + critical_t * est2$estimates4.sd/sqrt(sample.size.for.sim)
  
  h2 <- proxci_estimator$h
  q2 <- proxci_estimator$q_nt
  # dens_m <- estimator$dens_m
  
  for(j in 1:n.flds){
    params2[[j]]$h_params <- c(h2[[j]]$lambda_h,h2[[j]]$lambda_f,h2[[j]]$gamma, h2[[j]]$score)
    params2[[j]]$q_params <- c(q2[[j]]$lambda_q,q2[[j]]$lambda_f,q2[[j]]$gamma, q2[[j]]$score)
  }
  
  result.proxci2 <- list(true.ETT = target$ett.true,
                         est1 = est2$estimates1,
                         est2 = est2$estimates2,
                         est3 = est2$estimates3,
                         ifb.est = est2$estimates4,
                         ifb.sd = est2$estimates4.sd,
                         ifb.lower = ifb.lower2,
                         ifb.upper = ifb.upper2,
                         ifb.cov = target$ett.true >= ifb.lower2 & target$ett.true <= ifb.upper2,
                         sample.size = sample.size.for.sim,
                         params = params2
  )
  
  ################--proxci estimation (case3)--################
  est3 <- proxci_estimator$est3()
  
  ifb.lower3 <- est3$estimates4 - critical_t * est3$estimates4.sd/sqrt(sample.size.for.sim)
  ifb.upper3 <- est3$estimates4 + critical_t * est3$estimates4.sd/sqrt(sample.size.for.sim)
  
  h3 <- proxci_estimator$h_nt
  q3 <- proxci_estimator$q
  # dens_m <- estimator$dens_m
  
  for(j in 1:n.flds){
    params3[[j]]$h_params <- c(h3[[j]]$lambda_h,h3[[j]]$lambda_f,h3[[j]]$gamma, h3[[j]]$score)
    params3[[j]]$q_params <- c(q3[[j]]$lambda_q,q3[[j]]$lambda_f,q3[[j]]$gamma, q3[[j]]$score)
  }
  
  result.proxci3 <- list(true.ETT = target$ett.true,
                         est1 = est3$estimates1,
                         est2 = est3$estimates2,
                         est3 = est3$estimates3,
                         ifb.est = est3$estimates4,
                         ifb.sd = est3$estimates4.sd,
                         ifb.lower = ifb.lower3,
                         ifb.upper = ifb.upper3,
                         ifb.cov = target$ett.true >= ifb.lower3 & target$ett.true <= ifb.upper3,
                         sample.size = sample.size.for.sim,
                         params = params3
  )
  
  ################--proxci estimation (all FALSE)--################
  est4 <- proxci_estimator$est4()
  
  ifb.lower4 <- est4$estimates4 - critical_t * est4$estimates4.sd/sqrt(sample.size.for.sim)
  ifb.upper4 <- est4$estimates4 + critical_t * est4$estimates4.sd/sqrt(sample.size.for.sim)
  
  h4 <- proxci_estimator$h_nt
  q4 <- proxci_estimator$q_nt
  # dens_m <- estimator$dens_m
  
  for(j in 1:n.flds){
    params4[[j]]$h_params <- c(h4[[j]]$lambda_h,h4[[j]]$lambda_f,h4[[j]]$gamma, h4[[j]]$score)
    params4[[j]]$q_params <- c(q4[[j]]$lambda_q,q4[[j]]$lambda_f,q4[[j]]$gamma, q4[[j]]$score)
  }
  
  result.proxci4 <- list(true.ETT = target$ett.true,
                         est1 = est4$estimates1,
                         est2 = est4$estimates2,
                         est3 = est4$estimates3,
                         ifb.est = est4$estimates4,
                         ifb.sd = est4$estimates4.sd,
                         ifb.lower = ifb.lower4,
                         ifb.upper = ifb.upper4,
                         ifb.cov = target$ett.true >= ifb.lower4 & target$ett.true <= ifb.upper4,
                         sample.size = sample.size.for.sim,
                         params = params4
  )
  
  ################--BSIV estimation (all TRUE)--################
  cf_inds <- proxci_estimator$cf_inds
  BSIV_estimator0 <- run_para_BSIV(dataset,BSIV_model_spec$all_true,k,cf_inds)
  BSIV_estimator1 <- run_para_BSIV(dataset,BSIV_model_spec$case1,k,cf_inds)
  BSIV_estimator2 <- run_para_BSIV(dataset,BSIV_model_spec$case2,k,cf_inds)
  BSIV_estimator3 <- run_para_BSIV(dataset,BSIV_model_spec$case3,k,cf_inds)
  BSIV_estimator4 <- run_para_BSIV(dataset,BSIV_model_spec$case4,k,cf_inds)
  BSIV_estimator5 <- run_para_BSIV(dataset,BSIV_model_spec$all_false,k,cf_inds)
  
  BSIV.pi.est0 <- BSIV_estimator0$pi.est.ETTb
  BSIV.if.est0 <- BSIV_estimator0$if.est.ETTb
  BSIV.if.sd0 <- BSIV_estimator0$if.sd.ETTb
  BSIV.if.lower0 <- BSIV_estimator0$if.lower.ETTb
  BSIV.if.upper0 <- BSIV_estimator0$if.upper.ETTb
  BSIV.if.cov0 <- BSIV_estimator0$if.cov.ETTb
  BSIV.n.rows.removed0 <- BSIV_estimator0$n.rows.removed
  
  BSIV.pi.est1 <- BSIV_estimator1$pi.est.ETTb
  BSIV.if.est1 <- BSIV_estimator1$if.est.ETTb
  BSIV.if.sd1 <- BSIV_estimator1$if.sd.ETTb
  BSIV.if.lower1 <- BSIV_estimator1$if.lower.ETTb
  BSIV.if.upper1 <- BSIV_estimator1$if.upper.ETTb
  BSIV.if.cov1 <- BSIV_estimator1$if.cov.ETTb
  BSIV.n.rows.removed1 <- BSIV_estimator1$n.rows.removed
  
  BSIV.pi.est2 <- BSIV_estimator2$pi.est.ETTb
  BSIV.if.est2 <- BSIV_estimator2$if.est.ETTb
  BSIV.if.sd2 <- BSIV_estimator2$if.sd.ETTb
  BSIV.if.lower2 <- BSIV_estimator2$if.lower.ETTb
  BSIV.if.upper2 <- BSIV_estimator2$if.upper.ETTb
  BSIV.if.cov2 <- BSIV_estimator2$if.cov.ETTb
  BSIV.n.rows.removed2 <- BSIV_estimator2$n.rows.removed
  
  BSIV.pi.est3 <- BSIV_estimator3$pi.est.ETTb
  BSIV.if.est3 <- BSIV_estimator3$if.est.ETTb
  BSIV.if.sd3 <- BSIV_estimator3$if.sd.ETTb
  BSIV.if.lower3 <- BSIV_estimator3$if.lower.ETTb
  BSIV.if.upper3 <- BSIV_estimator3$if.upper.ETTb
  BSIV.if.cov3 <- BSIV_estimator3$if.cov.ETTb
  BSIV.n.rows.removed3 <- BSIV_estimator3$n.rows.removed
  
  BSIV.pi.est4 <- BSIV_estimator4$pi.est.ETTb
  BSIV.if.est4 <- BSIV_estimator4$if.est.ETTb
  BSIV.if.sd4 <- BSIV_estimator4$if.sd.ETTb
  BSIV.if.lower4 <- BSIV_estimator4$if.lower.ETTb
  BSIV.if.upper4 <- BSIV_estimator4$if.upper.ETTb
  BSIV.if.cov4 <- BSIV_estimator4$if.cov.ETTb
  BSIV.n.rows.removed4 <- BSIV_estimator4$n.rows.removed
  
  BSIV.pi.est5 <- BSIV_estimator5$pi.est.ETTb
  BSIV.if.est5 <- BSIV_estimator5$if.est.ETTb
  BSIV.if.sd5 <- BSIV_estimator5$if.sd.ETTb
  BSIV.if.lower5 <- BSIV_estimator5$if.lower.ETTb
  BSIV.if.upper5 <- BSIV_estimator5$if.upper.ETTb
  BSIV.if.cov5 <- BSIV_estimator5$if.cov.ETTb
  BSIV.n.rows.removed5 <- BSIV_estimator5$n.rows.removed
  
  result.BSIV0 <- list(true.ETT = target$ett.true,
                       BSIV.pi.est = BSIV.pi.est0,
                       BSIV.if.est = BSIV.if.est0,
                       BSIV.if.sd = BSIV.if.sd0,
                       BSIV.if.lower = BSIV.if.lower0,
                       BSIV.if.upper = BSIV.if.upper0,
                       BSIV.if.cov = BSIV.if.lower0 <= target$ett.true & BSIV.if.upper0 >= target$ett.true,
                       BSIV.n.rows.removed = BSIV.n.rows.removed0,
                       sample.size = sample.size.for.sim)
  
  result.BSIV1 <- list(true.ETT = target$ett.true,
                       BSIV.pi.est = BSIV.pi.est1,
                       BSIV.if.est = BSIV.if.est1,
                       BSIV.if.sd = BSIV.if.sd1,
                       BSIV.if.lower = BSIV.if.lower1,
                       BSIV.if.upper = BSIV.if.upper1,
                       BSIV.if.cov = BSIV.if.lower1 <= target$ett.true & BSIV.if.upper1 >= target$ett.true,
                       BSIV.n.rows.removed = BSIV.n.rows.removed1,
                       sample.size = sample.size.for.sim)
  
  result.BSIV2 <- list(true.ETT = target$ett.true,
                       BSIV.pi.est = BSIV.pi.est2,
                       BSIV.if.est = BSIV.if.est2,
                       BSIV.if.sd = BSIV.if.sd2,
                       BSIV.if.lower = BSIV.if.lower2,
                       BSIV.if.upper = BSIV.if.upper2,
                       BSIV.if.cov = BSIV.if.lower2 <= target$ett.true & BSIV.if.upper2 >= target$ett.true,
                       BSIV.n.rows.removed = BSIV.n.rows.removed2,
                       sample.size = sample.size.for.sim)
  
  result.BSIV3 <- list(true.ETT = target$ett.true,
                       BSIV.pi.est = BSIV.pi.est3,
                       BSIV.if.est = BSIV.if.est3,
                       BSIV.if.sd = BSIV.if.sd3,
                       BSIV.if.lower = BSIV.if.lower3,
                       BSIV.if.upper = BSIV.if.upper3,
                       BSIV.if.cov = BSIV.if.lower3 <= target$ett.true & BSIV.if.upper3 >= target$ett.true,
                       BSIV.n.rows.removed = BSIV.n.rows.removed3,
                       sample.size = sample.size.for.sim)
  
  result.BSIV4 <- list(true.ETT = target$ett.true,
                       BSIV.pi.est = BSIV.pi.est4,
                       BSIV.if.est = BSIV.if.est4,
                       BSIV.if.sd = BSIV.if.sd4,
                       BSIV.if.lower = BSIV.if.lower4,
                       BSIV.if.upper = BSIV.if.upper4,
                       BSIV.if.cov = BSIV.if.lower4 <= target$ett.true & BSIV.if.upper4 >= target$ett.true,
                       BSIV.n.rows.removed = BSIV.n.rows.removed4,
                       sample.size = sample.size.for.sim)
  
  result.BSIV5 <- list(true.ETT = target$ett.true,
                       BSIV.pi.est = BSIV.pi.est5,
                       BSIV.if.est = BSIV.if.est5,
                       BSIV.if.sd = BSIV.if.sd5,
                       BSIV.if.lower = BSIV.if.lower5,
                       BSIV.if.upper = BSIV.if.upper5,
                       BSIV.if.cov = BSIV.if.lower5 <= target$ett.true & BSIV.if.upper5 >= target$ett.true,
                       BSIV.n.rows.removed = BSIV.n.rows.removed5,
                       sample.size = sample.size.for.sim)
  
  return(list(result.proxci0 = result.proxci0,
              result.proxci1 = result.proxci1,
              result.proxci2 = result.proxci2,
              result.proxci3 = result.proxci3,
              result.proxci4 = result.proxci4,
              result.BSIV0 = result.BSIV0,
              result.BSIV1 = result.BSIV1,
              result.BSIV2 = result.BSIV2,
              result.BSIV3 = result.BSIV3,
              result.BSIV4 = result.BSIV4,
              result.BSIV5 = result.BSIV5))
}



sim.result <- run_sim(l)


