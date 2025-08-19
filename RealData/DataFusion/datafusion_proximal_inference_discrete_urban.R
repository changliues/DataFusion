source("DataFusion/datafusion_dataset.R")
source("DataFusion/proximal_model_specification_urban.R")
proximal_model_spec <- readRDS(file = "proximal_model_spec_urban.rds")

library(dplyr)
library(caret)
library(kernlab)
library(pracma)
library(evmix)
library(nnet)

ProximalInference <- setRefClass(
  "ProximalInference",
  fields = list(
    data = "ProxCIData",
    crossfit_folds = "numeric",
    cf_inds = "list",
    h = "list",
    q = "list",
    mod_spec = "list",
    lbd = "numeric",
    Z_cat = "numeric",
    Z_cat_n = "integer",
    M_cat = "numeric",
    M_cat_n = "integer"
  ),
  methods = list(
    
    initialize = function(
    proxci_dataset,
    crossfit_folds = 1,
    mod_spec = proximal_model_spec,
    lbd = 0
    ) {
      stopifnot(crossfit_folds >= 1)
      .self$data <- proxci_dataset
      .self$crossfit_folds <- crossfit_folds
      .self$cf_inds <- .self$data$create_crossfit_split(crossfit_folds)
      .self$mod_spec <- mod_spec
      .self$lbd <- lbd
      
      .self$h <- lapply(1:crossfit_folds, function(i) {
        .self$estimate_h(fold = i)
      })
      
      .self$q <- lapply(1:crossfit_folds, function(i) {
        .self$estimate_q(fold = i)
      })
      
      .self$Z_cat <- sort(unique(proxci_dataset$Z))
      .self$Z_cat_n <- length(unique(proxci_dataset$Z))
      
      .self$M_cat <- sort(unique(proxci_dataset$M.fac))
      .self$M_cat_n <- length(unique(proxci_dataset$M.fac))
      
    },
    
    estimate_h = function(fold = 0) {
      dat <- data.frame(.self$data$Y, .self$data$M.fac, .self$data$Z, .self$data$A, .self$data$X, .self$data$G)[.self$cf_inds[[fold]]$train, ]
      colnames(dat) <- c("Y", "M", "Z", "A", "X1", "X2", "X3", "G")
      
      y.azo.fit <- glm(.self$mod_spec$y.azo,
                   data =  data.frame(y=dat$Y,a=dat$A, z=dat$Z,x1=dat$X1,x2=dat$X2,x3=dat$X3)[dat$G==0,],
                   family = "gaussian")
      
      y.aZo <- function(a,x1,x2,x3){
        predict(y.azo.fit, newdata=data.frame(a=a,
                                              z=.self$Z_cat,
                                              x1=x1,
                                              x2=x2,
                                              x3=x3))
      }
      
      # pm.azo.fit <- function(m,a,z,x1,x2,x3){
      #   nrow(dat[dat$M==m & dat$A==a & dat$Z==z & dat$X1==x1 & dat$X2==x2 & dat$X3==x3 & dat$G==0,])/nrow(dat[dat$A==a & dat$Z==z & dat$X1==x1 & dat$X2==x2 & dat$X3==x3 & dat$G==0,])
      # }
      
      pm.azo.fit <- multinom(.self$mod_spec$m.azo,
                             data =  data.frame(m=as.factor(dat$M),a=dat$A, z=dat$Z,x1=dat$X1,x2=dat$X2,x3=dat$X3)[dat$G==0,],
                             trace = F)
      
      pM.aZo <- function(a,x1,x2,x3){
        
        t(predict(pm.azo.fit, newdata = data.frame(a=a,z=.self$Z_cat,x1=x1,x2=x2,x3=x3), type = "probs"))
      }
      
      h_mod <- function(m,a,x1,x2,x3){
        
        n <- length(x1)
        
        if(length(m)==1){
          m.vec <- rep(m,n)
        }else{
          m.vec <- m
        }
        
        if(length(a)==1){
          a.vec <- rep(a,n)
        }else{
          a.vec <- a
        }
        
        result <- c()
        for(i in 1:n){
          mat <- pM.aZo(a.vec[i],x1[i],x2[i],x3[i])
          #+matrix(rnorm(12,0.01,0.1),nrow=4)
          inv.mat <- solve(t(mat) %*% mat + .self$lbd*diag(3)) %*% t(mat)
          # inv.mat <- pinv(mat)
          result <- c(result,c(y.aZo(a.vec[i],x1[i],x2[i],x3[i]) %*% inv.mat)[m.vec[i]])
        }
        
        return(result)
      }
      
      return(list(h_mod = h_mod))
    },
    
    estimate_q = function(fold = 0) {
      dat <- data.frame(.self$data$M.fac, .self$data$Z, .self$data$A, .self$data$X, .self$data$G)[.self$cf_inds[[fold]]$train, ]
      colnames(dat) <- c("M", "Z", "A", "X1", "X2", "X3", "G")
      
      # pm.ao.fit <- function(m,a,x1,x2,x3){
      #   nrow(dat[dat$M==m & dat$A==a & dat$X1==x1 & dat$X2==x2 & dat$X3==x3  & dat$G==0,])/nrow(dat[dat$A==a & dat$X1==x1 & dat$X2==x2 & dat$X3==x3  & dat$G==0,])
      # }
      # 
      # pm.ae.fit <- function(m,a,x1,x2,x3){
      #   nrow(dat[dat$M==m & dat$A==a & dat$X1==x1 & dat$X2==x2 & dat$X3==x3  & dat$G==1,])/nrow(dat[dat$A==a & dat$X1==x1 & dat$X2==x2 & dat$X3==x3  & dat$G==1,])
      # }
      # 
      # pa.o.fit <- function(a,x1,x2,x3){
      #   nrow(dat[dat$A==a & dat$X1==x1 & dat$X2==x2 & dat$X3==x3  & dat$G==0,])/nrow(dat[dat$X1==x1 & dat$X2==x2 & dat$X3==x3  & dat$G==0,])
      # }
      
      pm.ao.fit <- multinom(.self$mod_spec$m.ao,
                            data =  data.frame(m=as.factor(dat$M),a=dat$A,x1=dat$X1,x2=dat$X2,x3=dat$X3)[dat$G==0,],
                            trace = F)
      
      pm.ae.fit <- multinom(.self$mod_spec$m.ae,
                            data =  data.frame(m=as.factor(dat$M),a=dat$A,x1=dat$X1,x2=dat$X2,x3=dat$X3)[dat$G==1,],
                            trace = F)
      
      pa.o.fit <- glm(.self$mod_spec$a.ox,
                      data =  data.frame(a=dat$A, x1=dat$X1,x2=dat$X2,x3=dat$X3)[dat$G==0,],
                      family = "binomial")
      
      pa.Mo <- function(a,x1,x2,x3){
        pa.o.hat <- expit(predict(pa.o.fit, newdata = data.frame(a=a,x1=x1,x2=x2,x3=x3)))
        pm.ao.hat <- predict(pm.ao.fit, newdata = data.frame(a=a,x1=x1,x2=x2,x3=x3), type = "probs")
        pm.ae.hat <- predict(pm.ae.fit, newdata = data.frame(a=a,x1=x1,x2=x2,x3=x3), type = "probs")
        pm.ae.hat/(pm.ao.hat*pa.o.hat)
        
        # pa.1o <- pm.ae.fit(1,a,x1,x2,x3)/(pm.ao.fit(1,a,x1,x2,x3)*pa.o.fit(a,x1,x2,x3))
        # pa.2o <- pm.ae.fit(2,a,x1,x2,x3)/(pm.ao.fit(2,a,x1,x2,x3)*pa.o.fit(a,x1,x2,x3))
        # pa.3o <- pm.ae.fit(3,a,x1,x2,x3)/(pm.ao.fit(3,a,x1,x2,x3)*pa.o.fit(a,x1,x2,x3))
        # c(pa.1o,pa.2o,pa.3o)
      }
      
      # pz.amo.fit <- function(z,a,m,x1,x2,x3){
      #   nrow(dat[dat$Z==z & dat$A==a & dat$M==m & dat$X1==x1 & dat$X2==x2 & dat$X3==x3 & dat$G==0,])/nrow(dat[dat$A==a & dat$M==m & dat$X1==x1 & dat$X2==x2 & dat$X3==x3 & dat$G==0,])
      # }
      
      pz.amo.fit <- multinom(.self$mod_spec$z.amo,
                             data =  data.frame(z=as.factor(dat$Z),a=dat$A, m=dat$M,x1=dat$X1,x2=dat$X2,x3=dat$X3)[dat$G==0,],
                             trace = F)
      
      
      pZ.aMo <- function(a,x1,x2,x3){
        t(predict(pz.amo.fit, newdata = data.frame(a=a,m=.self$M_cat,x1=x1,x2=x2,x3=x3), type = "probs"))
      }
      
      q_mod <- function(z,a,x1,x2,x3){
        
        n <- length(z)
        
        if(length(a)==1){
          a.vec <- rep(a,n)
        }else{
          a.vec <- a
        }
        
        result <- c()
        for(i in 1:n){
          result <- c(result,(pa.Mo(a.vec[i],x1[i],x2[i],x3[i]) %*% pinv(pZ.aMo(a.vec[i],x1[i],x2[i],x3[i]),tol = 1e-30))[z[i]])
        }
        
        return(result)
      }
      
      return(list(q_mod = q_mod))
    },
    
    est = function(reduction = mean) {
      estimates1 <- c()
      estimates2 <- c()
      estimates3 <- c()
      estimates4 <- c()
      estimates4.if <- list()
      for (fold in 1:.self$crossfit_folds) {
        ### regression
        fold_train_idx <- .self$cf_inds[[fold]]$train
        x1 <- .self$data$X[fold_train_idx,1]
        x2 <- .self$data$X[fold_train_idx,2]
        x3 <- .self$data$X[fold_train_idx,3]
        g <- .self$data$G[fold_train_idx]
        pge.x.fit <- glm(.self$mod_spec$g.x,
                         data =  data.frame(g=g,x1=x1,x2=x2,x3=x3),
                         family = "binomial")
        
        x1e <- x1[g==1]
        x2e <- x2[g==1]
        x3e <- x3[g==1]
        a <- .self$data$A[fold_train_idx]
        ae <- a[g==1]
        pa.ex.fit <- glm(.self$mod_spec$a.ex,
                         data =  data.frame(ae=ae,x1e=x1e,x2e=x2e,x3e=x3e),
                         family = "binomial")
        
        # m <- .self$data$M[fold_train_idx]
        # me <- m[g==1]
        # m0e <- m[g==1 & a==0]
        m <- .self$data$M.fac[fold_train_idx]
        me <- m[g==1]
        m0e <- m[g==1 & a==0]
        x10e <- x1[g==1 & a==0]
        x20e <- x2[g==1 & a==0]
        x30e <- x3[g==1 & a==0]
        
        pm.ae.fit <-  multinom(.self$mod_spec$m.ae,
                               data =  data.frame(m=as.factor(me),a=ae,x1=x1e,x2=x2e,x3=x3e),
                               trace = F)
        
        ### evaluation
        fold_idx <- .self$cf_inds[[fold]]$eval
        Igo <- .self$data$G[fold_idx] == 0
        Ige <- .self$data$G[fold_idx] == 1
        pgo <- mean(.self$data$G[fold_idx] == 0)
        pa.o <- mean(.self$data$A[fold_idx][.self$data$G[fold_idx] == 0])
        Ia0 <- .self$data$A[fold_idx] == 0
        Ia1 <- .self$data$A[fold_idx] == 1
        y <- .self$data$Y[fold_idx]
        
        pge.x.hat <- expit(predict(pge.x.fit, newdata = data.frame(x1 = .self$data$X[fold_idx,1],
                                                                   x2 = .self$data$X[fold_idx,2],
                                                                   x3 = .self$data$X[fold_idx,3])))
        pa.ex.hat <- expit(predict(pa.ex.fit, newdata = data.frame(x1e = .self$data$X[fold_idx,1],
                                                                   x2e = .self$data$X[fold_idx,2],
                                                                   x3e = .self$data$X[fold_idx,3])))
        pM.0e.hat <- predict(pm.ae.fit, newdata = data.frame(a=0,x1=.self$data$X[fold_idx,1],x2=.self$data$X[fold_idx,2],x3=.self$data$X[fold_idx,3]), type = "probs")
        
        hm0.hat <- matrix(NA, nrow = length(fold_idx), ncol = .self$M_cat_n)
        for(l in 1:.self$M_cat_n){
          hm0.hat[,l] <- .self$h[[fold]]$h_mod(.self$M_cat[l],0,.self$data$X[fold_idx,1],.self$data$X[fold_idx,2],.self$data$X[fold_idx,3])
        }
        eta0.hat <- c(rowSums(pM.0e.hat * hm0.hat))
        
        h0.hat <- .self$h[[fold]]$h_mod(.self$data$M.fac[fold_idx],0,.self$data$X[fold_idx,1],.self$data$X[fold_idx,2],.self$data$X[fold_idx,3])
        q0.hat <- .self$q[[fold]]$q_mod(.self$data$Z[fold_idx],0,.self$data$X[fold_idx,1],.self$data$X[fold_idx,2],.self$data$X[fold_idx,3])
        
        est1 <- mean((Igo/(pgo*pa.o)) * (y-eta0.hat))
        
        est2 <- mean((1/(pgo*pa.o)) *(Igo*y - (Ige*Ia0)/(1-pa.ex.hat)*h0.hat*(1/pge.x.hat-1)))
        
        est3 <- mean((Igo/(pgo*pa.o)) * (y - Ia0*y*q0.hat))
        
        ests4 <- (1/(pgo*pa.o)) * (Igo * (y-Ia0*q0.hat*(y-h0.hat)-eta0.hat)-Ige*Ia0/(1-pa.ex.hat)*(h0.hat-eta0.hat)*(1/pge.x.hat-1))
        est4 <- mean(ests4)
        estimates4.if[[fold]] <- ests4
        
        estimates1 <- c(estimates1, est1)
        estimates2 <- c(estimates2, est2)
        estimates3 <- c(estimates3, est3)
        estimates4 <- c(estimates4, est4)
      }
      
      estimates4=mean(estimates4)
      estimates4.sd <- c()
      
      for(fold in 1:(.self$crossfit_folds)) {
        estimates4.sd <- c(estimates4.sd, sd(estimates4.if[[fold]]-estimates4))
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
