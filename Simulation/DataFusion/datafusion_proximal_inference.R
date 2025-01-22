source("DataFusion/datafusion_dataset.R")
source("DataFusion/datafusion_minimax.R")
source("DataFusion/proximal_model_specification.R")
proximal_model_spec <- readRDS(file = "proximal_model_spec.rds")

library(dplyr)
library(caret)
library(kernlab)
library(pracma)
library(evmix)

ProximalInference <- setRefClass(
  "ProximalInference",
  fields = list(
    data = "ProxCIData",
    crossfit_folds = "numeric",
    lambdash = "numeric",
    lambdasfh = "numeric",
    gammash = "numeric",
    lambdasq = "numeric",
    lambdasfq = "numeric",
    gammasq = "numeric",
    lambdash_nt = "numeric",
    lambdasfh_nt = "numeric",
    gammash_nt = "numeric",
    lambdasq_nt = "numeric",
    lambdasfq_nt = "numeric",
    gammasq_nt = "numeric",
    split_ratio = "numeric",
    print_best_params = "logical",
    cf_inds = "list",
    h = "list",
    q = "list",
    hdim = "numeric",
    fdim = "numeric",
    h_nt = "list",
    q_nt = "list",
    hdim_nt = "numeric",
    fdim_nt = "numeric",
    mod_spec = "list"
  ),
  methods = list(
    initialize = function(
    proxci_dataset,
    crossfit_folds = 1,
    lambdash = NULL,
    lambdasfh = NULL,
    gammash = NULL,
    lambdasq = NULL,
    lambdasfq = NULL,
    gammasq = NULL,
    lambdash_nt = NULL,
    lambdasfh_nt = NULL,
    gammash_nt = NULL,
    lambdasq_nt = NULL,
    lambdasfq_nt = NULL,
    gammasq_nt = NULL,
    split_ratio = 0.25,
    print_best_params = FALSE,
    mod_spec = proximal_model_spec
    ) {
      stopifnot(crossfit_folds >= 1)
      .self$data <- proxci_dataset
      .self$crossfit_folds <- crossfit_folds
      .self$lambdash <- lambdash
      .self$lambdasfh <- lambdasfh
      .self$gammash <- gammash
      .self$lambdasq <- lambdasq
      .self$lambdasfq <- lambdasfq
      .self$gammasq <- gammasq
      .self$lambdash_nt <- lambdash_nt
      .self$lambdasfh_nt <- lambdasfh_nt
      .self$gammash_nt <- gammash_nt
      .self$lambdasq_nt <- lambdasq_nt
      .self$lambdasfq_nt <- lambdasfq_nt
      .self$gammasq_nt <- gammasq_nt
      .self$split_ratio <- split_ratio
      .self$print_best_params <- print_best_params
      
      .self$hdim <- ncol(.self$data$r1())
      .self$fdim <- ncol(.self$data$r2())
      .self$hdim_nt <- ncol(.self$data$r1())
      .self$fdim_nt <- ncol(.self$data$r2())
      
      .self$cf_inds <- .self$data$create_crossfit_split(crossfit_folds)
      .self$mod_spec <- mod_spec
      
      .self$h <- lapply(1:crossfit_folds, function(i) {
        .self$estimate_h(fold = i)
      })
      
      .self$q <- lapply(1:crossfit_folds, function(i) {
        .self$estimate_q(fold = i)
      })
      
      .self$h_nt <- lapply(1:crossfit_folds, function(i) {
        .self$estimate_h_nt(fold = i)
      })
      
      .self$q_nt <- lapply(1:crossfit_folds, function(i) {
        .self$estimate_q_nt(fold = i)
      })
      
    },
    
    estimate_h = function(fold = 0) {
      pgo <- mean(.self$data$G == 0)
      g1 <- as.matrix((.self$data$G == 0) / pgo)
      g2 <- as.matrix(-1 * .self$data$Y * (.self$data$G == 0)/pgo)
      dat <- join_data(.self$data$r1(), .self$data$r2(), g1, g2)[.self$cf_inds[[fold]]$train, ]
      
      search <- MinimaxRKHSCV(
        dat,
        .self$hdim,
        .self$fdim,
        lambdash = .self$lambdash,
        lambdasf = .self$lambdasfh,
        gammas = .self$gammash,
        split_ratio = .self$split_ratio
      )
      
      if (.self$print_best_params) {
        cat("h: gamma = ", search$best_params_$gamma, ", lambda_f = ", search$best_params_$lambda_f, ", lambda_h = ", search$best_params_$lambda_h, ", best score: ", search$best_score_, "\n")
      }
      
      return(list(h_mod = search$best_model_$h_,
                  lambda_h = search$best_params_$lambda_h,
                  lambda_f = search$best_params_$lambda_f,
                  gamma = search$best_params_$gamma,
                  score = search$best_score_))
    },
    
    estimate_q = function(fold = 0) {
      fold_train_idx <- .self$cf_inds[[fold]]$train
      x1 <- .self$data$X[fold_train_idx,1]
      x2 <- .self$data$X[fold_train_idx,2]
      x3 <- .self$data$X[fold_train_idx,3]
      g <- .self$data$G[fold_train_idx]
      a <- .self$data$A[fold_train_idx]
      pge.x.fit <- glm(.self$mod_spec$proximal_correct_model$g.x,
                       data =  data.frame(g=g,x1=x1,x2=x2,x3=x3),
                       family = "binomial")
      pge.x.hat <- expit(predict(pge.x.fit, newdata = data.frame(x1 = .self$data$X[fold_train_idx,1],
                                                                 x2 = .self$data$X[fold_train_idx,2],
                                                                 x3 = .self$data$X[fold_train_idx,3])))
      
      x1e <- x1[g==1]
      x2e <- x2[g==1]
      x3e <- x3[g==1]
      ae <- a[g==1]
      pa.ex.fit <- glm(.self$mod_spec$proximal_correct_model$a.ex,
                       data =  data.frame(ae=ae,x1e=x1e,x2e=x2e,x3e=x3e),
                       family = "binomial")
      pa.ex.hat <- expit(predict(pa.ex.fit, newdata = data.frame(x1e = .self$data$X[fold_train_idx,1],
                                                                 x2e = .self$data$X[fold_train_idx,2],
                                                                 x3e = .self$data$X[fold_train_idx,3])))
      pae.x.hat <- (a*pa.ex.hat+(1-a)*(1-pa.ex.hat))*pge.x.hat
      
      g1 <- as.matrix((g == 0)/(1-pge.x.hat))
      g2 <- as.matrix(-1 * (g == 1)/pae.x.hat)
      dat <- join_data(.self$data$r2()[fold_train_idx,], .self$data$r1()[fold_train_idx,], g1, g2)
      
      search <- MinimaxRKHSCV(
        dat,
        .self$hdim,
        .self$fdim,
        lambdash = .self$lambdasq,
        lambdasf = .self$lambdasfq,
        gammas = .self$gammasq,
        split_ratio = .self$split_ratio
      )
      
      if (.self$print_best_params) {
        cat("q: gamma = ", search$best_params_$gamma, ", lambda_f = ", search$best_params_$lambda_f, ", lambda_q = ", search$best_params_$lambda_h, ", best score: ", search$best_score_, "\n")
      }
      
      return(list(q_mod = search$best_model_$h_,
                  lambda_q = search$best_params_$lambda_h,
                  lambda_f = search$best_params_$lambda_f,
                  gamma = search$best_params_$gamma,
                  score = search$best_score_))
    },
    
    estimate_h_nt = function(fold = 0) {
      pgo <- mean(.self$data$G == 0)
      g1 <- as.matrix((.self$data$G == 0) / pgo)
      g2 <- as.matrix(-1 * .self$data$Y * (.self$data$G == 0)/pgo)
      dat <- join_data(.self$data$r1(), .self$data$r2(), g1, g2)[.self$cf_inds[[fold]]$train, ]
      
      search <- MinimaxRKHSCV(
        dat,
        .self$hdim_nt,
        .self$fdim_nt,
        lambdash = .self$lambdash_nt,
        lambdasf = .self$lambdasfh_nt,
        gammas = .self$gammash_nt,
        split_ratio = .self$split_ratio
      )
      
      if (.self$print_best_params) {
        cat("h: gamma = ", search$best_params_$gamma, ", lambda_f = ", search$best_params_$lambda_f, ", lambda_h = ", search$best_params_$lambda_h, ", best score: ", search$best_score_, "\n")
      }
      
      return(list(h_mod = search$best_model_$h_,
                  lambda_h = search$best_params_$lambda_h,
                  lambda_f = search$best_params_$lambda_f,
                  gamma = search$best_params_$gamma,
                  score = search$best_score_))
    },
    
    estimate_q_nt = function(fold = 0) {
      fold_train_idx <- .self$cf_inds[[fold]]$train
      x1s <- .self$data$Xs[fold_train_idx,1]
      x2s <- .self$data$Xs[fold_train_idx,2]
      x3s <- .self$data$Xs[fold_train_idx,3]
      g <- .self$data$G[fold_train_idx]
      a <- .self$data$A[fold_train_idx]
      pge.x.fit <- glm(.self$mod_spec$proximal_incorrect_model$g.x,
                       data =  data.frame(g=g,x1s=x1s,x2s=x2s,x3s=x3s),
                       family = "binomial")
      pge.x.hat <- expit(predict(pge.x.fit, newdata = data.frame(x1s = .self$data$X[fold_train_idx,1],
                                                                 x2s = .self$data$X[fold_train_idx,2],
                                                                 x3s = .self$data$X[fold_train_idx,3])))
      
      x1se <- x1s[g==1]
      x2se <- x2s[g==1]
      x3se <- x3s[g==1]
      ae <- a[g==1]
      pa.ex.fit <- glm(.self$mod_spec$proximal_incorrect_model$a.ex,
                       data =  data.frame(ae=ae,x1se=x1se,x2se=x2se,x3se=x3se),
                       family = "binomial")
      pa.ex.hat <- expit(predict(pa.ex.fit, newdata = data.frame(x1se = .self$data$X[fold_train_idx,1],
                                                                 x2se = .self$data$X[fold_train_idx,2],
                                                                 x3se = .self$data$X[fold_train_idx,3])))
      pae.x.hat <- (a*pa.ex.hat+(1-a)*(1-pa.ex.hat))*pge.x.hat
      
      g1 <- as.matrix((g == 0)/(1-pge.x.hat))
      g2 <- as.matrix(-1 * (g == 1)/pae.x.hat)
      dat <- join_data(.self$data$r2()[fold_train_idx,], .self$data$r1()[fold_train_idx,], g1, g2)
      
      search <- MinimaxRKHSCV(
        dat,
        .self$hdim_nt,
        .self$fdim_nt,
        lambdash = .self$lambdasq_nt,
        lambdasf = .self$lambdasfq_nt,
        gammas = .self$gammasq_nt,
        split_ratio = .self$split_ratio
      )
      
      if (.self$print_best_params) {
        cat("q: gamma = ", search$best_params_$gamma, ", lambda_f = ", search$best_params_$lambda_f, ", lambda_q = ", search$best_params_$lambda_h, ", best score: ", search$best_score_, "\n")
      }
      
      return(list(q_mod = search$best_model_$h_,
                  lambda_q = search$best_params_$lambda_h,
                  lambda_f = search$best_params_$lambda_f,
                  gamma = search$best_params_$gamma,
                  score = search$best_score_))
    },
    
    est0 = function(reduction = mean) {
      estimates1 <- c()
      estimates2 <- c()
      estimates3 <- c()
      estimates4 <- c()
      estimates4.sd <- c()
      for (fold in 1:.self$crossfit_folds) {
        ### regression
        fold_train_idx <- .self$cf_inds[[fold]]$train
        x1 <- .self$data$X[fold_train_idx,1]
        x2 <- .self$data$X[fold_train_idx,2]
        x3 <- .self$data$X[fold_train_idx,3]
        g <- .self$data$G[fold_train_idx]
        pge.x.fit <- glm(.self$mod_spec$all_true$g.x,
                         data =  data.frame(g=g,x1=x1,x2=x2,x3=x3),
                         family = "binomial")
        
        x1e <- x1[g==1]
        x2e <- x2[g==1]
        x3e <- x3[g==1]
        a <- .self$data$A[fold_train_idx]
        ae <- a[g==1]
        pa.ex.fit <- glm(.self$mod_spec$all_true$a.ex,
                         data =  data.frame(ae=ae,x1e=x1e,x2e=x2e,x3e=x3e),
                         family = "binomial")
        
        m <- .self$data$M[fold_train_idx]
        me <- m[g==1]
        m0e <- m[g==1 & a==0]
        x10e <- x1[g==1 & a==0]
        x20e <- x2[g==1 & a==0]
        x30e <- x3[g==1 & a==0]
        
        rax <- matrix(cbind(m0e,0,x10e,x20e,x30e),ncol = 5)
        w <- as.vector(.self$h[[fold]]$h_mod(rax))
        
        eta.fit <-  lm(.self$mod_spec$all_true$eta,
                       data =  data.frame(w=w,x10e=x10e,x20e=x20e,x30e=x30e),)
        
        ### evaluation
        fold_idx <- .self$cf_inds[[fold]]$eval
        Igo <- .self$data$G[fold_idx] == 0
        Ige <- .self$data$G[fold_idx] == 1
        pgo <- mean(.self$data$G[fold_idx] == 0)
        pa.o <- mean(.self$data$A[fold_idx][.self$data$G[fold_idx] == 0])
        Ia0 <- .self$data$A[fold_idx] == 0
        Ia1 <- .self$data$A[fold_idx] == 1
        rh <- .self$data$r1()[fold_idx, ]
        rq <- .self$data$r2()[fold_idx, ]
        y <- .self$data$Y[fold_idx]
        
        pge.x.hat <- expit(predict(pge.x.fit, newdata = data.frame(x1 = .self$data$X[fold_idx,1],
                                                                   x2 = .self$data$X[fold_idx,2],
                                                                   x3 = .self$data$X[fold_idx,3])))
        pa.ex.hat <- expit(predict(pa.ex.fit, newdata = data.frame(x1e = .self$data$X[fold_idx,1],
                                                                   x2e = .self$data$X[fold_idx,2],
                                                                   x3e = .self$data$X[fold_idx,3])))
        eta0.hat <- predict(eta.fit,newdata = data.frame(x10e = .self$data$X[fold_idx,1],
                                                         x20e = .self$data$X[fold_idx,2],
                                                         x30e = .self$data$X[fold_idx,3]))
        
        est1 <- mean((Igo/(pgo*pa.o)) * (y-eta0.hat))
        
        est2 <- mean((1/(pgo*pa.o)) *(Igo*y - (Ige*Ia0)/(1-pa.ex.hat)*.self$h[[fold]]$h_mod(rh)*(1/pge.x.hat-1)))
        
        est3 <- mean(Igo/(pgo*pa.o)) * y*(1 - Ia0*.self$q[[fold]]$q_mod(rq))
        
        ests4 <- (1/(pgo*pa.o)) * (Igo * (y-Ia0*.self$q[[fold]]$q_mod(rq)*(y-.self$h[[fold]]$h_mod(rh))-eta0.hat)-
                                     Ige*Ia0/(1-pa.ex.hat)*(.self$h[[fold]]$h_mod(rh)-eta0.hat)*(1/pge.x.hat-1))
        est4 <- mean(ests4)
        est4.sd <- sd(ests4)
        
        estimates1 <- c(estimates1, est1)
        estimates2 <- c(estimates2, est2)
        estimates3 <- c(estimates3, est3)
        estimates4 <- c(estimates4, est4)
        estimates4.sd <- c(estimates4.sd, est4.sd)
      }
      if (is.null(reduction)) {
        return(list(estimates1 = estimates1,
                    estimates2 = estimates2,
                    estimates3 = estimates3,
                    estimates4 = estimates4,
                    estimates4.sd = estimates4.sd
        ))
      } else {
        return(list(estimates1 = reduction(estimates1),
                    estimates2 = reduction(estimates2),
                    estimates3 = reduction(estimates3),
                    estimates4 = reduction(estimates4),
                    estimates4.sd = sqrt(reduction(estimates4.sd^2))
        ))
      }
    },
    
    est1 = function(reduction = mean) {
      estimates1 <- c()
      estimates2 <- c()
      estimates3 <- c()
      estimates4 <- c()
      estimates4.sd <- c()
      for (fold in 1:.self$crossfit_folds) {
        ### regression
        fold_train_idx <- .self$cf_inds[[fold]]$train
        x1 <- .self$data$X[fold_train_idx,1]
        x2 <- .self$data$X[fold_train_idx,2]
        x3 <- .self$data$X[fold_train_idx,3]
        x1s <- .self$data$Xs[fold_train_idx,1]
        x2s <- .self$data$Xs[fold_train_idx,2]
        x3s <- .self$data$Xs[fold_train_idx,3]
        g <- .self$data$G[fold_train_idx]
        pge.x.fit <- glm(.self$mod_spec$case1$g.x,
                         data =  data.frame(g=g,x1s=x1s,x2s=x2s,x3s=x3s),
                         family = "binomial")
        
        x1e <- x1[g==1]
        x2e <- x2[g==1]
        x3e <- x3[g==1]
        x1se <- x1s[g==1]
        x2se <- x2s[g==1]
        x3se <- x3s[g==1]
        a <- .self$data$A[fold_train_idx]
        ae <- a[g==1]
        pa.ex.fit <- glm(.self$mod_spec$case1$a.ex,
                         data =  data.frame(ae=ae,x1se=x1se,x2se=x2se,x3se=x3se),
                         family = "binomial")
        
        m <- .self$data$M[fold_train_idx]
        me <- m[g==1]
        m0e <- m[g==1 & a==0]
        x10e <- x1[g==1 & a==0]
        x20e <- x2[g==1 & a==0]
        x30e <- x3[g==1 & a==0]
        x1s0e <- x1s[g==1 & a==0]
        x2s0e <- x2s[g==1 & a==0]
        x3s0e <- x3s[g==1 & a==0]
        
        rax <- matrix(cbind(m0e,0,x10e,x20e,x30e),ncol = 5)
        w <- as.vector(.self$h[[fold]]$h_mod(rax))
        
        eta.fit <-  lm(.self$mod_spec$case1$eta,
                       data =  data.frame(w=w,x10e=x10e,x20e=x20e,x30e=x30e),)
        
        ### evaluation
        fold_idx <- .self$cf_inds[[fold]]$eval
        Igo <- .self$data$G[fold_idx] == 0
        Ige <- .self$data$G[fold_idx] == 1
        pgo <- mean(.self$data$G[fold_idx] == 0)
        pa.o <- mean(.self$data$A[fold_idx][.self$data$G[fold_idx] == 0])
        Ia0 <- .self$data$A[fold_idx] == 0
        Ia1 <- .self$data$A[fold_idx] == 1
        rh <- .self$data$r1()[fold_idx, ]
        rq <- .self$data$r2()[fold_idx, ]
        y <- .self$data$Y[fold_idx]
        
        pge.x.hat <- expit(predict(pge.x.fit, newdata = data.frame(x1s = .self$data$X[fold_idx,1],
                                                                   x2s = .self$data$X[fold_idx,2],
                                                                   x3s = .self$data$X[fold_idx,3])))
        pa.ex.hat <- expit(predict(pa.ex.fit, newdata = data.frame(x1se = .self$data$X[fold_idx,1],
                                                                   x2se = .self$data$X[fold_idx,2],
                                                                   x3se = .self$data$X[fold_idx,3])))
        eta0.hat <- predict(eta.fit,newdata = data.frame(x10e = .self$data$X[fold_idx,1],
                                                         x20e = .self$data$X[fold_idx,2],
                                                         x30e = .self$data$X[fold_idx,3]))
        
        est1 <- mean((Igo/(pgo*pa.o)) * (y-eta0.hat))
        
        est2 <- mean((1/(pgo*pa.o)) *(Igo*y - (Ige*Ia0)/(1-pa.ex.hat)*.self$h[[fold]]$h_mod(rh)*(1/pge.x.hat-1)))
        
        est3 <- mean((Igo/(pgo*pa.o)) * y*(1 - Ia0*.self$q_nt[[fold]]$q_mod(rq)))
        
        ests4 <- (1/(pgo*pa.o)) * (Igo * (y-Ia0*.self$q_nt[[fold]]$q_mod(rq)*(y-.self$h[[fold]]$h_mod(rh))-eta0.hat)-
                                     Ige*Ia0/(1-pa.ex.hat)*(.self$h[[fold]]$h_mod(rh)-eta0.hat)*(1/pge.x.hat-1))
        est4 <- mean(ests4)
        est4.sd <- sd(ests4)
        
        estimates1 <- c(estimates1, est1)
        estimates2 <- c(estimates2, est2)
        estimates3 <- c(estimates3, est3)
        estimates4 <- c(estimates4, est4)
        estimates4.sd <- c(estimates4.sd, est4.sd)
      }
      if (is.null(reduction)) {
        return(list(estimates1 = estimates1,
                    estimates2 = estimates2,
                    estimates3 = estimates3,
                    estimates4 = estimates4,
                    estimates4.sd = estimates4.sd
        ))
      } else {
        return(list(estimates1 = reduction(estimates1),
                    estimates2 = reduction(estimates2),
                    estimates3 = reduction(estimates3),
                    estimates4 = reduction(estimates4),
                    estimates4.sd = sqrt(reduction(estimates4.sd^2))
        ))
      }
    },
    
    est2 = function(reduction = mean) {
      estimates1 <- c()
      estimates2 <- c()
      estimates3 <- c()
      estimates4 <- c()
      estimates4.sd <- c()
      for (fold in 1:.self$crossfit_folds) {
        ### regression
        fold_train_idx <- .self$cf_inds[[fold]]$train
        x1 <- .self$data$X[fold_train_idx,1]
        x2 <- .self$data$X[fold_train_idx,2]
        x3 <- .self$data$X[fold_train_idx,3]
        x1s <- .self$data$Xs[fold_train_idx,1]
        x2s <- .self$data$Xs[fold_train_idx,2]
        x3s <- .self$data$Xs[fold_train_idx,3]
        g <- .self$data$G[fold_train_idx]
        pge.x.fit <- glm(.self$mod_spec$case2$g.x,
                         data =  data.frame(g=g,x1=x1,x2=x2,x3=x3),
                         family = "binomial")
        
        x1e <- x1[g==1]
        x2e <- x2[g==1]
        x3e <- x3[g==1]
        x1se <- x1s[g==1]
        x2se <- x2s[g==1]
        x3se <- x3s[g==1]
        a <- .self$data$A[fold_train_idx]
        ae <- a[g==1]
        pa.ex.fit <- glm(.self$mod_spec$case2$a.ex,
                         data =  data.frame(ae=ae,x1e=x1e,x2e=x2e,x3e=x3e),
                         family = "binomial")
        
        m <- .self$data$M[fold_train_idx]
        me <- m[g==1]
        m0e <- m[g==1 & a==0]
        x10e <- x1[g==1 & a==0]
        x20e <- x2[g==1 & a==0]
        x30e <- x3[g==1 & a==0]
        x1s0e <- x1s[g==1 & a==0]
        x2s0e <- x2s[g==1 & a==0]
        x3s0e <- x3s[g==1 & a==0]
        
        rax <- matrix(cbind(m0e,0,x10e,x20e,x30e),ncol = 5)
        w <- as.vector(.self$h[[fold]]$h_mod(rax))
        
        eta.fit <-  lm(.self$mod_spec$case2$eta,
                       data =  data.frame(w=w,x1s0e=x1s0e,x2s0e=x2s0e,x3s0e=x3s0e),)
        
        ### evaluation
        fold_idx <- .self$cf_inds[[fold]]$eval
        Igo <- .self$data$G[fold_idx] == 0
        Ige <- .self$data$G[fold_idx] == 1
        pgo <- mean(.self$data$G[fold_idx] == 0)
        pa.o <- mean(.self$data$A[fold_idx][.self$data$G[fold_idx] == 0])
        Ia0 <- .self$data$A[fold_idx] == 0
        Ia1 <- .self$data$A[fold_idx] == 1
        rh <- .self$data$r1()[fold_idx, ]
        rq <- .self$data$r2()[fold_idx, ]
        y <- .self$data$Y[fold_idx]
        
        pge.x.hat <- expit(predict(pge.x.fit, newdata = data.frame(x1 = .self$data$X[fold_idx,1],
                                                                   x2 = .self$data$X[fold_idx,2],
                                                                   x3 = .self$data$X[fold_idx,3])))
        pa.ex.hat <- expit(predict(pa.ex.fit, newdata = data.frame(x1e = .self$data$X[fold_idx,1],
                                                                   x2e = .self$data$X[fold_idx,2],
                                                                   x3e = .self$data$X[fold_idx,3])))
        eta0.hat <- predict(eta.fit,newdata = data.frame(x1s0e = .self$data$X[fold_idx,1],
                                                         x2s0e = .self$data$X[fold_idx,2],
                                                         x3s0e = .self$data$X[fold_idx,3]))
        
        est1 <- mean((Igo/(pgo*pa.o)) * (y-eta0.hat))
        
        est2 <- mean((1/(pgo*pa.o)) *(Igo*y - (Ige*Ia0)/(1-pa.ex.hat)*.self$h[[fold]]$h_mod(rh)*(1/pge.x.hat-1)))
        
        est3 <- mean((Igo/(pgo*pa.o)) * y*(1 - Ia0*.self$q_nt[[fold]]$q_mod(rq)))
        
        ests4 <- (1/(pgo*pa.o)) * (Igo * (y-Ia0*.self$q_nt[[fold]]$q_mod(rq)*(y-.self$h[[fold]]$h_mod(rh))-eta0.hat)-
                                     Ige*Ia0/(1-pa.ex.hat)*(.self$h[[fold]]$h_mod(rh)-eta0.hat)*(1/pge.x.hat-1))
        est4 <- mean(ests4)
        est4.sd <- sd(ests4)
        
        estimates1 <- c(estimates1, est1)
        estimates2 <- c(estimates2, est2)
        estimates3 <- c(estimates3, est3)
        estimates4 <- c(estimates4, est4)
        estimates4.sd <- c(estimates4.sd, est4.sd)
      }
      if (is.null(reduction)) {
        return(list(estimates1 = estimates1,
                    estimates2 = estimates2,
                    estimates3 = estimates3,
                    estimates4 = estimates4,
                    estimates4.sd = estimates4.sd
        ))
      } else {
        return(list(estimates1 = reduction(estimates1),
                    estimates2 = reduction(estimates2),
                    estimates3 = reduction(estimates3),
                    estimates4 = reduction(estimates4),
                    estimates4.sd = sqrt(reduction(estimates4.sd^2))
        ))
      }
    },
    
    est3 = function(reduction = mean) {
      estimates1 <- c()
      estimates2 <- c()
      estimates3 <- c()
      estimates4 <- c()
      estimates4.sd <- c()
      for (fold in 1:.self$crossfit_folds) {
        ### regression
        fold_train_idx <- .self$cf_inds[[fold]]$train
        x1 <- .self$data$X[fold_train_idx,1]
        x2 <- .self$data$X[fold_train_idx,2]
        x3 <- .self$data$X[fold_train_idx,3]
        x1s <- .self$data$Xs[fold_train_idx,1]
        x2s <- .self$data$Xs[fold_train_idx,2]
        x3s <- .self$data$Xs[fold_train_idx,3]
        g <- .self$data$G[fold_train_idx]
        pge.x.fit <- glm(.self$mod_spec$case3$g.x,
                         data =  data.frame(g=g,x1=x1,x2=x2,x3=x3),
                         family = "binomial")
        
        x1e <- x1[g==1]
        x2e <- x2[g==1]
        x3e <- x3[g==1]
        x1se <- x1s[g==1]
        x2se <- x2s[g==1]
        x3se <- x3s[g==1]
        a <- .self$data$A[fold_train_idx]
        ae <- a[g==1]
        pa.ex.fit <- glm(.self$mod_spec$case3$a.ex,
                         data =  data.frame(ae=ae,x1e=x1e,x2e=x2e,x3e=x3e),
                         family = "binomial")
        
        m <- .self$data$M[fold_train_idx]
        me <- m[g==1]
        m0e <- m[g==1 & a==0]
        x10e <- x1[g==1 & a==0]
        x20e <- x2[g==1 & a==0]
        x30e <- x3[g==1 & a==0]
        x1s0e <- x1s[g==1 & a==0]
        x2s0e <- x2s[g==1 & a==0]
        x3s0e <- x3s[g==1 & a==0]
        
        rax <- matrix(cbind(m0e,0,x10e,x20e,x30e),ncol = 5)
        w <- as.vector(.self$h_nt[[fold]]$h_mod(rax))
        
        eta.fit <-  lm(.self$mod_spec$case3$eta,
                       data =  data.frame(w=w,x1s0e=x1s0e,x2s0e=x2s0e,x3s0e=x3s0e),)
        
        ### evaluation
        fold_idx <- .self$cf_inds[[fold]]$eval
        Igo <- .self$data$G[fold_idx] == 0
        Ige <- .self$data$G[fold_idx] == 1
        pgo <- mean(.self$data$G[fold_idx] == 0)
        pa.o <- mean(.self$data$A[fold_idx][.self$data$G[fold_idx] == 0])
        Ia0 <- .self$data$A[fold_idx] == 0
        Ia1 <- .self$data$A[fold_idx] == 1
        rh <- .self$data$r1()[fold_idx, ]
        rq <- .self$data$r2()[fold_idx, ]
        y <- .self$data$Y[fold_idx]
        
        pge.x.hat <- expit(predict(pge.x.fit, newdata = data.frame(x1 = .self$data$X[fold_idx,1],
                                                                   x2 = .self$data$X[fold_idx,2],
                                                                   x3 = .self$data$X[fold_idx,3])))
        pa.ex.hat <- expit(predict(pa.ex.fit, newdata = data.frame(x1e = .self$data$X[fold_idx,1],
                                                                   x2e = .self$data$X[fold_idx,2],
                                                                   x3e = .self$data$X[fold_idx,3])))
        eta0.hat <- predict(eta.fit,newdata = data.frame(x1s0e = .self$data$X[fold_idx,1],
                                                         x2s0e = .self$data$X[fold_idx,2],
                                                         x3s0e = .self$data$X[fold_idx,3]))
        
        est1 <- mean((Igo/(pgo*pa.o)) * (y-eta0.hat))
        
        est2 <- mean((1/(pgo*pa.o)) *(Igo*y - (Ige*Ia0)/(1-pa.ex.hat)*.self$h_nt[[fold]]$h_mod(rh)*(1/pge.x.hat-1)))
        
        est3 <- mean((Igo/(pgo*pa.o)) * y*(1 - Ia0*.self$q[[fold]]$q_mod(rq)))
        
        ests4 <- (1/(pgo*pa.o)) * (Igo * (y-Ia0*.self$q[[fold]]$q_mod(rq)*(y-.self$h_nt[[fold]]$h_mod(rh))-eta0.hat)-
                                     Ige*Ia0/(1-pa.ex.hat)*(.self$h_nt[[fold]]$h_mod(rh)-eta0.hat)*(1/pge.x.hat-1))
        est4 <- mean(ests4)
        est4.sd <- sd(ests4)
        
        estimates1 <- c(estimates1, est1)
        estimates2 <- c(estimates2, est2)
        estimates3 <- c(estimates3, est3)
        estimates4 <- c(estimates4, est4)
        estimates4.sd <- c(estimates4.sd, est4.sd)
      }
      if (is.null(reduction)) {
        return(list(estimates1 = estimates1,
                    estimates2 = estimates2,
                    estimates3 = estimates3,
                    estimates4 = estimates4,
                    estimates4.sd = estimates4.sd
        ))
      } else {
        return(list(estimates1 = reduction(estimates1),
                    estimates2 = reduction(estimates2),
                    estimates3 = reduction(estimates3),
                    estimates4 = reduction(estimates4),
                    estimates4.sd = sqrt(reduction(estimates4.sd^2))
        ))
      }
    },
    
    est4 = function(reduction = mean) {
      estimates1 <- c()
      estimates2 <- c()
      estimates3 <- c()
      estimates4 <- c()
      estimates4.sd <- c()
      for (fold in 1:.self$crossfit_folds) {
        ### regression
        fold_train_idx <- .self$cf_inds[[fold]]$train
        x1 <- .self$data$X[fold_train_idx,1]
        x2 <- .self$data$X[fold_train_idx,2]
        x3 <- .self$data$X[fold_train_idx,3]
        x1s <- .self$data$Xs[fold_train_idx,1]
        x2s <- .self$data$Xs[fold_train_idx,2]
        x3s <- .self$data$Xs[fold_train_idx,3]
        g <- .self$data$G[fold_train_idx]
        pge.x.fit <- glm(.self$mod_spec$all_false$g.x,
                         data =  data.frame(g=g,x1s=x1s,x2s=x2s,x3s=x3s),
                         family = "binomial")
        
        x1e <- x1[g==1]
        x2e <- x2[g==1]
        x3e <- x3[g==1]
        x1se <- x1s[g==1]
        x2se <- x2s[g==1]
        x3se <- x3s[g==1]
        a <- .self$data$A[fold_train_idx]
        ae <- a[g==1]
        pa.ex.fit <- glm(.self$mod_spec$all_false$a.ex,
                         data =  data.frame(ae=ae,x1se=x1se,x2se=x2se,x3se=x3se),
                         family = "binomial")
        
        m <- .self$data$M[fold_train_idx]
        me <- m[g==1]
        m0e <- m[g==1 & a==0]
        x10e <- x1[g==1 & a==0]
        x20e <- x2[g==1 & a==0]
        x30e <- x3[g==1 & a==0]
        x1s0e <- x1s[g==1 & a==0]
        x2s0e <- x2s[g==1 & a==0]
        x3s0e <- x3s[g==1 & a==0]
        
        rax <- matrix(cbind(m0e,0,x10e,x20e,x30e),ncol = 5)
        w <- as.vector(.self$h_nt[[fold]]$h_mod(rax))
        
        eta.fit <-  lm(.self$mod_spec$all_false$eta,
                       data =  data.frame(w=w,x1s0e=x1s0e,x2s0e=x2s0e,x3s0e=x3s0e),)
        
        ### evaluation
        fold_idx <- .self$cf_inds[[fold]]$eval
        Igo <- .self$data$G[fold_idx] == 0
        Ige <- .self$data$G[fold_idx] == 1
        pgo <- mean(.self$data$G[fold_idx] == 0)
        pa.o <- mean(.self$data$A[fold_idx][.self$data$G[fold_idx] == 0])
        Ia0 <- .self$data$A[fold_idx] == 0
        Ia1 <- .self$data$A[fold_idx] == 1
        rh <- .self$data$r1()[fold_idx, ]
        rq <- .self$data$r2()[fold_idx, ]
        y <- .self$data$Y[fold_idx]
        
        pge.x.hat <- expit(predict(pge.x.fit, newdata = data.frame(x1s = .self$data$X[fold_idx,1],
                                                                   x2s = .self$data$X[fold_idx,2],
                                                                   x3s = .self$data$X[fold_idx,3])))
        pa.ex.hat <- expit(predict(pa.ex.fit, newdata = data.frame(x1se = .self$data$X[fold_idx,1],
                                                                   x2se = .self$data$X[fold_idx,2],
                                                                   x3se = .self$data$X[fold_idx,3])))
        eta0.hat <- predict(eta.fit,newdata = data.frame(x1s0e = .self$data$X[fold_idx,1],
                                                         x2s0e = .self$data$X[fold_idx,2],
                                                         x3s0e = .self$data$X[fold_idx,3]))
        
        est1 <- mean((Igo/(pgo*pa.o)) * (y-eta0.hat))
        
        est2 <- mean((1/(pgo*pa.o)) *(Igo*y - (Ige*Ia0)/(1-pa.ex.hat)*.self$h_nt[[fold]]$h_mod(rh)*(1/pge.x.hat-1)))
        
        est3 <- mean((Igo/(pgo*pa.o)) * y*(1 - Ia0*.self$q_nt[[fold]]$q_mod(rq)))
        
        ests4 <- (1/(pgo*pa.o)) * (Igo * (y-Ia0*.self$q_nt[[fold]]$q_mod(rq)*(y-.self$h_nt[[fold]]$h_mod(rh))-eta0.hat)-
                                     Ige*Ia0/(1-pa.ex.hat)*(.self$h_nt[[fold]]$h_mod(rh)-eta0.hat)*(1/pge.x.hat-1))
        est4 <- mean(ests4)
        est4.sd <- sd(ests4)
        
        estimates1 <- c(estimates1, est1)
        estimates2 <- c(estimates2, est2)
        estimates3 <- c(estimates3, est3)
        estimates4 <- c(estimates4, est4)
        estimates4.sd <- c(estimates4.sd, est4.sd)
      }
      if (is.null(reduction)) {
        return(list(estimates1 = estimates1,
                    estimates2 = estimates2,
                    estimates3 = estimates3,
                    estimates4 = estimates4,
                    estimates4.sd = estimates4.sd
        ))
      } else {
        return(list(estimates1 = reduction(estimates1),
                    estimates2 = reduction(estimates2),
                    estimates3 = reduction(estimates3),
                    estimates4 = reduction(estimates4),
                    estimates4.sd = sqrt(reduction(estimates4.sd^2))
        ))
      }
    }
    
    
  )
)
