library(caret)
library(kernlab)
library(pracma)


source("DataFusion/datafusion_dataset.R")
source("DataFusion/minimax_core.R")

# Define the MinimaxRKHS class
MinimaxRKHS <- setRefClass(
  "MinimaxRKHS",
  fields = list(
    hdim = "numeric",
    fdim = "numeric",
    lambda_h = "numeric",
    lambda_f = "numeric",
    lambda_score = "numeric",
    gamma = "numeric",
    gamma_score = "numeric",
    Kh_ = "matrix",
    Kf_ = "matrix",
    alpha_ = "numeric",
    beta_ = "numeric",
    r1_ = "matrix",
    r2_ = "matrix"
  ),
  
  methods = list(
    initialize = function(hdim, fdim, lambda_h = 1, lambda_f = 1, lambda_score = 1, gamma = 1, gamma_score = 1) {
      .self$hdim <- hdim
      .self$fdim <- fdim
      .self$lambda_h <- lambda_h
      .self$lambda_f <- lambda_f
      .self$lambda_score <- lambda_score
      .self$gamma <- gamma
      .self$gamma_score <- gamma_score
    },
    
    .create_kernel = function(X, Y = NULL) {
      kernelMatrix(kernel = rbfdot(sigma = .self$gamma), 
                   x = X, 
                   y = Y)
    },
    
    .unpack_data = function(data) {
      r1 <- data[, 1:.self$hdim]
      r2 <- data[, (.self$hdim + 1):(.self$hdim + .self$fdim)]
      g1 <- data[, .self$hdim + .self$fdim + 1]
      g2 <- data[, .self$hdim + .self$fdim + 2]
      list(r1 = r1, r2 = r2, g1 = g1, g2 = g2)
    },
    
    fit = function(data) {
      unpacked <- .self$.unpack_data(data)
      r1 <- unpacked$r1
      r2 <- unpacked$r2
      g1 <- unpacked$g1
      g2 <- unpacked$g2
      
      .self$Kh_ <- .self$.create_kernel(r1)
      .self$Kf_ <- .self$.create_kernel(r2)
      
      # Solve minimax problem
      solution <- minimax_solve(.self$Kh_, .self$Kf_, g1, g2, .self$lambda_h, .self$lambda_f)
      .self$alpha_ <- solution$alpha
      .self$beta_ <- solution$beta
      
      .self$r1_ <- r1
      .self$r2_ <- r2
      
      .self
    },
    
    predict = function(data) {
      if (is.null(.self$alpha_) | is.null(.self$beta_)) {
        stop("Model is not fitted.")
      }
      r1 <- .self$.unpack_data(data)$r1
      .self$h_(r1)
    },
    
    h_ = function(r){
      .self$alpha_ %*% .self$.create_kernel(.self$r1_, r)
    },
    
    f_ = function(r){
      .self$beta_ %*% .self$.create_kernel(.self$r2_, r)
    },
    
    score = function(data) {
      
      unpacked <- .self$.unpack_data(data)
      r1 <- unpacked$r1
      r2 <- unpacked$r2
      g1 <- unpacked$g1
      g2 <- unpacked$g2
      gammas <- .self$gamma_score
      lambdafs <- .self$lambda_score
      
      h_val <- .self$h_(r1)
      
      scrs <- rep(NA,length(gammas)*length(lambdafs))
      idx <- 1
      for(i in 1:length(gammas)){
        Kf_score <- kernelMatrix(kernel = rbfdot(sigma = gammas[i]),
                                 x = r2,
                                 y = NULL)
        for(j in 1:length(lambdafs)){
          lambda2 <- lambdafs[j]
          
          scrs[idx] <- score_nuisance_function(h_val, Kf_score, g1, g2, lambda2)
          
          idx <- idx+1
        }
      }
      
      mean(scrs)
      
    }
  )
)


library(caret)

# Define the custom model

MinimaxRKHSCV <- function(data, hdim, fdim, lambdash = NULL, lambdasf = NULL, gammas = NULL, split_ratio = 0.25) {
  if (is.null(lambdash)) {
    lambdash <- 10^seq(-5, 0)
  }
  if (is.null(lambdasf)) {
    lambdasf <- 10^seq(-5, 0)
  }
  if (is.null(gammas)) {
    gammas <- 1
  }
  
  tune_grid <- expand.grid(lambda_h = lambdash, lambda_f = lambdasf, gamma = gammas)
  
  minimax_models <- list()
  scores <- rep(NA,nrow(tune_grid))
  
  # initialize models
  for(i in 1:nrow(tune_grid)){
    
    train_model <- MinimaxRKHS$new(hdim = hdim, fdim = fdim, 
                                   lambda_h = tune_grid$lambda_h[i],
                                   lambda_f = tune_grid$lambda_f[i],
                                   lambda_score = lambdasf,
                                   gamma = tune_grid$gamma[i],
                                   gamma_score = gammas)
    minimax_models[[i]] <- train_model
  }
  
  split_idx <- createDataPartition(1:nrow(data), times = 1, p = split_ratio, list = TRUE)
  scores <- rep(NA, nrow(tune_grid))
  if(split_ratio==1){
    for(i in 1:nrow(tune_grid)){
      minimax_models[[i]]$fit(data[split_idx[[1]],])
      # train_model$predict(data[flds[[k]],])
      scores[i] <- minimax_models[[i]]$score(data[split_idx[[1]],])
    }
  }
  else{
    for(i in 1:nrow(tune_grid)){
        minimax_models[[i]]$fit(data[-(split_idx[[1]]),])
        # train_model$predict(data[flds[[k]],])
        scores[i] <- minimax_models[[i]]$score(data[split_idx[[1]],])
    }
  }

  best_idx <- which.max(scores)
  
  best_model <- minimax_models[[best_idx]]
  best_params <- list(lambda_h = tune_grid$lambda_h[best_idx], lambda_f = tune_grid$lambda_f[best_idx], gamma = tune_grid$gamma[best_idx])
  best_score <- scores[best_idx]
  
  return(list(best_model_ = best_model, best_params_ = best_params, best_score_ = best_score))
}
