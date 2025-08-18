run_para_BSIV <- function(dataset,used_model,crossfit_folds=4,BSIV_cf_inds){
  
  n <- dataset$n
  
  y <- dataset$Y
  m <- dataset$M
  x1 <- dataset$X[,1]
  x2 <- dataset$X[,2]
  z <- dataset$X[,3]
  a <- dataset$A
  gn <- dataset$G
  x1s <- dataset$Xs[,1]
  x2s <- dataset$Xs[,2]
  
  # Create dataframe
  df <- data.frame(y = y, m = m, x1 = x1, x2 = x2, z = z, a = a, gn = gn, x1s = x1s, x2s = x2s)
  colnames(df) <- c("y", "m", "x1", "x2", "z", "a", "gn", "x1s", "x2s")
  
  # Estimation
  ## P(A=1,G=O)
  p.1o <- nrow(df[a == 1 & gn == 0,])/nrow(df)
  
  ## P(A=1|G=O)
  pa.o <- nrow(df[a == 1 & gn == 0,])/nrow(df[gn == 0,])
  
  ## E(M|G=O)
  m.o <- mean(df$m[gn == 0])
  
  cf_inds <- BSIV_cf_inds
  g.xz.hat <- rep(NA,n)
  z.hat <- rep(NA,n)
  z.o.hat <- rep(NA,n)
  z.e.hat <- rep(NA,n)
  a.ze.hat <- rep(NA,n)
  a.zo.hat <- rep(NA,n)
  a.1o.hat <- rep(NA,n)
  a.0o.hat <- rep(NA,n)
  m.11o.hat <- rep(NA,n)
  m.01o.hat <- rep(NA,n)
  m.10o.hat <- rep(NA,n)
  m.00o.hat <- rep(NA,n)
  m.0ze.hat <- rep(NA,n)
  y.11o.hat <- rep(NA,n)
  y.01o.hat <- rep(NA,n)
  y.10o.hat <- rep(NA,n)
  y.00o.hat <- rep(NA,n)
  for(fold in 1:crossfit_folds){
    
    fold_train_idx <- cf_inds[[fold]]$train
    fold_idx <- cf_inds[[fold]]$eval
    
    df.train <- df[fold_train_idx,]
    
    #Models are correctly specified(or at least closely approximated, 3rd power)
    
    ## P(G=E|X)
    g.xz.fit <- glm(used_model$g.xz, data = df[fold_train_idx,], family = "binomial")
    g.xz.hat[fold_idx] <- expit(predict(g.xz.fit, 
                                        newdata = data.frame(x1=x1[fold_idx],
                                                             x1s=x1[fold_idx],
                                                             x2=x2[fold_idx],
                                                             x2s=x2[fold_idx],
                                                             z=z[fold_idx])))
    
    ## P(Z=1|X,G=O)
    z.o.fit <- glm(used_model$z.o, data = df.train[df.train$gn == 0,], family = "binomial")
    z.o.hat[fold_idx] <- expit(predict(z.o.fit,
                                       newdata = data.frame(x1=x1[fold_idx],
                                                            x1s=x1[fold_idx],
                                                            x2=x2[fold_idx],
                                                            x2s=x2[fold_idx])))
    
    ## P(Z=1|X,G=E)
    z.e.fit <- glm(used_model$z.e, data = df.train[df.train$gn == 1,], family = "binomial")
    z.e.hat[fold_idx] <- expit(predict(z.e.fit,
                                       newdata = data.frame(x1=x1[fold_idx],
                                                            x1s=x1[fold_idx],
                                                            x2=x2[fold_idx],
                                                            x2s=x2[fold_idx])))
    
    ## P(A=1|X,Z,G=E)
    a.ze.fit <- glm(used_model$a.ze, data = df.train[df.train$gn == 1,], family = "binomial")
    a.ze.hat[fold_idx] <- expit(predict(a.ze.fit,
                                        newdata = data.frame(x1=x1[fold_idx],
                                                             x1s=x1[fold_idx],
                                                             x2=x2[fold_idx],
                                                             x2s=x2[fold_idx],
                                                             z=z[fold_idx])))
    
    ## P(A=1|X,Z,G=O)
    a.zo.fit <- glm(used_model$a.zo, data = df.train[df.train$gn == 0,], family = "binomial")
    a.zo.hat[fold_idx] <- expit(predict(a.zo.fit, 
                                        newdata = data.frame(x1=x1[fold_idx],
                                                             x1s=x1[fold_idx],
                                                             x2=x2[fold_idx],
                                                             x2s=x2[fold_idx],
                                                             z=z[fold_idx])))
    
    ## P(A=1|Z=1,X,G=O)
    a.1o.hat[fold_idx] <- expit(predict(a.zo.fit, 
                                        newdata = data.frame(x1=x1[fold_idx],
                                                             x1s=x1[fold_idx],
                                                             x2=x2[fold_idx],
                                                             x2s=x2[fold_idx],
                                                             z=1)))
    
    ## P(A=1|Z=0,X,G=O)
    a.0o.hat[fold_idx] <- expit(predict(a.zo.fit,
                                        newdata = data.frame(x1=x1[fold_idx],
                                                             x1s=x1[fold_idx],
                                                             x2=x2[fold_idx],
                                                             x2s=x2[fold_idx],
                                                             z=0)))
    
    m.azo.fit <- glm(used_model$m.azo, data = df.train[df.train$gn == 0,])
    y.azo.fit <- glm(used_model$y.azo, data = df.train[df.train$gn == 0,])
    ## E(M|A=1,Z=1,X,G=O)
    m.11o.hat[fold_idx] <- predict(m.azo.fit,
                                   newdata = data.frame(x1=x1[fold_idx],
                                                        x1s=x1[fold_idx],
                                                        x2=x2[fold_idx],
                                                        x2s=x2[fold_idx],
                                                        z=1,a=1))
    
    ## E(Y|A=1,Z=1,X,G=O)
    y.11o.hat[fold_idx] <- predict(y.azo.fit,
                                   newdata = data.frame(x1=x1[fold_idx],
                                                        x1s=x1[fold_idx],
                                                        x2=x2[fold_idx],
                                                        x2s=x2[fold_idx],
                                                        z=1,a=1))
    
    ## E(M|A=1,Z=0,X,G=O)
    m.10o.hat[fold_idx] <- predict(m.azo.fit,
                                   newdata = data.frame(x1=x1[fold_idx],
                                                        x1s=x1[fold_idx],
                                                        x2=x2[fold_idx],
                                                        x2s=x2[fold_idx],
                                                        z=0,a=1))
    
    ## E(Y|A=1,Z=0,X,G=O)
    y.10o.hat[fold_idx] <- predict(y.azo.fit,
                                   newdata = data.frame(x1=x1[fold_idx],
                                                        x1s=x1[fold_idx],
                                                        x2=x2[fold_idx],
                                                        x2s=x2[fold_idx],
                                                        z=0,a=1))
    
    ## E(M|A=0,Z=1,X,G=O)
    m.01o.hat[fold_idx] <- predict(m.azo.fit,
                                   newdata = data.frame(x1=x1[fold_idx],
                                                        x1s=x1[fold_idx],
                                                        x2=x2[fold_idx],
                                                        x2s=x2[fold_idx],
                                                        z=1,a=0))
    
    ## E(Y|A=0,Z=1,X,G=O)
    y.01o.hat[fold_idx] <- predict(y.azo.fit,
                                   newdata = data.frame(x1=x1[fold_idx],
                                                        x1s=x1[fold_idx],
                                                        x2=x2[fold_idx],
                                                        x2s=x2[fold_idx],
                                                        z=1,a=0))
    
    ## E(M|A=0,Z=0,X,G=O)
    m.00o.hat[fold_idx] <- predict(m.azo.fit,
                                   newdata = data.frame(x1=x1[fold_idx],
                                                        x1s=x1[fold_idx],
                                                        x2=x2[fold_idx],
                                                        x2s=x2[fold_idx],
                                                        z=0,a=0))
    
    ## E(Y|A=0,Z=0,X,G=O)
    y.00o.hat[fold_idx] <- predict(y.azo.fit,
                                   newdata = data.frame(x1=x1[fold_idx],
                                                        x1s=x1[fold_idx],
                                                        x2=x2[fold_idx],
                                                        x2s=x2[fold_idx],
                                                        z=0,a=0))
    
    ## E(M|A=0,X,Z,G=E)
    m.aze.fit <- glm(used_model$m.aze, data = df.train[df.train$gn == 1,])
    m.0ze.hat[fold_idx] <- predict(m.aze.fit,
                                   newdata = data.frame(x1=x1[fold_idx],
                                                        x1s=x1[fold_idx],
                                                        x2=x2[fold_idx],
                                                        x2s=x2[fold_idx],
                                                        z=z[fold_idx],
                                                        a=0))
  }
  
  idx.rm <- which(abs(a.1o.hat-a.0o.hat) <= sqrt(1/n))
  if(length(idx.rm) > 0){
    
    df[idx.rm,] <- NA
    # y[idx.rm] <- NA
    # m[idx.rm] <- NA
    # x1[idx.rm] <- NA
    # x2[idx.rm] <- NA
    z[idx.rm] <- NA
    a[idx.rm] <- NA
    gn[idx.rm] <- NA
    # x1s[idx.rm] <- NA
    # x2s[idx.rm] <- NA
    
    g.xz.hat[idx.rm] <- NA
    z.o.hat[idx.rm] <- NA
    z.e.hat[idx.rm] <- NA
    a.ze.hat[idx.rm] <- NA
    a.zo.hat[idx.rm] <- NA
    a.1o.hat[idx.rm] <- NA
    a.0o.hat[idx.rm] <- NA
    m.11o.hat[idx.rm] <- NA
    y.11o.hat[idx.rm] <- NA
    m.10o.hat[idx.rm] <- NA
    y.10o.hat[idx.rm] <- NA
    m.01o.hat[idx.rm] <- NA
    y.01o.hat[idx.rm] <- NA
    m.00o.hat[idx.rm] <- NA
    y.00o.hat[idx.rm] <- NA
    m.0ze.hat[idx.rm] <- NA
  }
  
  # P(A=1|X,G=O)=P(A=1|Z=1,X,G=O)P(Z=1|X,G=O)+P(A=1|Z=0,X,G=O)(1-P(Z=1|X,G=O))
  a.o.hat <- a.1o.hat*z.o.hat + a.0o.hat*(1-z.o.hat)
  
  ## influence-function of psi_1^az
  z.1o <- df$z[a==1 & gn == 0]
  ## influence-function of psi_1^11
  y.11o.hat.1o <- y.11o.hat[a==1 & gn == 0]
  m.11o.hat.1o <- m.11o.hat[a==1 & gn == 0]
  
  psi1.11.hat <- mean((y.11o.hat.1o-m.11o.hat.1o)*z.1o, na.rm = T)
  
  if.psi1.11 <- ((gn == 0)/p.1o) * ((a==1)*(z==1)*(y-m-(y.11o.hat-m.11o.hat))+
                                      (a==1)*((y.11o.hat-m.11o.hat)*z-psi1.11.hat))
  
  ## influence-function of psi_1^01
  y.01o.hat.1o <- y.01o.hat[a==1 & gn == 0]
  m.01o.hat.1o <- m.01o.hat[a==1 & gn == 0]
  
  psi1.01.hat <- mean((y.01o.hat.1o-m.01o.hat.1o)*z.1o, na.rm = T)
  
  if.psi1.01 <- ((gn == 0)/p.1o)*((a==0)*(z==1)*(a.1o.hat/(1-a.1o.hat))*(y-m-(y.01o.hat-m.01o.hat))+
                                    (a==1)*((y.01o.hat-m.01o.hat)*z-psi1.01.hat))
  
  ## influence-function of psi_1^10
  y.10o.hat.1o <- y.10o.hat[a==1 & gn == 0]
  m.10o.hat.1o <- m.10o.hat[a==1 & gn == 0]
  
  psi1.10.hat <- mean((y.10o.hat.1o-m.10o.hat.1o)*z.1o, na.rm = T)
  
  if.psi1.10 <- ((gn == 0)/p.1o)*((a==1)*(z==0)*((z.o.hat*a.1o.hat)/((1-z.o.hat)*a.0o.hat))*(y-m-(y.10o.hat-m.10o.hat))+
                                    (a==1)*((y.10o.hat-m.10o.hat)*z-psi1.10.hat))
  
  ## influence-function of psi_1^00
  y.00o.hat.1o <- y.00o.hat[a==1 & gn == 0]
  m.00o.hat.1o <- m.00o.hat[a==1 & gn == 0]
  
  psi1.00.hat <- mean((y.00o.hat.1o-m.00o.hat.1o)*z.1o, na.rm = T)
  
  if.psi1.00 <- ((gn == 0)/p.1o)*((a==0)*(z==0)*((z.o.hat*a.1o.hat)/((1-z.o.hat)*(1-a.0o.hat)))*(y-m-(y.00o.hat-m.00o.hat))+
                                    (a==1)*((y.00o.hat-m.00o.hat)*z-psi1.00.hat))
  
  ## influence-function of psi_2^az
  ## influence-function of psi_2^10
  psi2.10.hat <- mean((y.10o.hat.1o-m.10o.hat.1o), na.rm = T)
  
  if.psi2.10 <- ((gn == 0)/p.1o)*(((a==1)*(z==0))/a.0o.hat*(a.o.hat/(1-z.o.hat))*(y-m-(y.10o.hat-m.10o.hat))+
                                    (a==1)*(y.10o.hat-m.10o.hat-psi2.10.hat))
  
  ## influence-function of psi_2^00
  psi2.00.hat <- mean((y.00o.hat.1o-m.00o.hat.1o), na.rm = T)
  
  if.psi2.00 <- ((gn == 0)/p.1o)*(((a==0)*(z==0))/(1-a.0o.hat)*(a.o.hat/(1-z.o.hat))*(y-m-(y.00o.hat-m.00o.hat))+
                                    (a==1)*(y.00o.hat-m.00o.hat-psi2.00.hat))
  
  ## influence-function of psi_3
  a.1o.hat.1o <- a.1o.hat[a==1 & gn == 0]
  a.0o.hat.1o <- a.0o.hat[a==1 & gn == 0]
  
  psi3.hat <- mean(((y.01o.hat.1o-m.01o.hat.1o)-(y.00o.hat.1o-m.00o.hat.1o))/((1-a.1o.hat.1o)-(1-a.0o.hat.1o)), na.rm = T)
  
  if.psi3 <- ((gn == 0)/((1-a.1o.hat)-(1-a.0o.hat))) * (a.o.hat/p.1o) *
    ( (a==0)*((z==1)/((1-a.1o.hat)*z.o.hat)*(y-m-(y.01o.hat-m.01o.hat)) -
                (z==0)/((1-a.0o.hat)*(1-z.o.hat))*(y-m-(y.00o.hat-m.00o.hat))) +
        ((y.01o.hat-m.01o.hat)-(y.00o.hat-m.00o.hat))/((1-a.1o.hat)-(1-a.0o.hat)) *
        (-(z==1)/z.o.hat*((a==0)-(1-a.1o.hat))+(z==0)/(1-z.o.hat)*((a==0)-(1-a.0o.hat)))) +
    ((a==1)*(gn == 0)/p.1o)*(((y.01o.hat-m.01o.hat)-(y.00o.hat-m.00o.hat))/((1-a.1o.hat)-(1-a.0o.hat))-psi3.hat)
  
  ## influence-function of psi_5
  m.0ze.hat.1o <- m.0ze.hat[a==1 & gn == 0]
  a.zo.hat.1o <- a.zo.hat[a==1 & gn == 0]
  
  psi5.hat <- mean(m.0ze.hat.1o/a.zo.hat.1o, na.rm = T)
  
  if.psi5 <- 1/p.1o * (((a==0)*(gn == 1)/(1-a.ze.hat))*(1/g.xz.hat-1)*(m-m.0ze.hat) +
                         (gn == 0)*m.0ze.hat - (a==1)*(gn == 0)*psi5.hat)
  
  
  ## influence-function of psi_8 = psi_6 + psi_7
  psi8.hat <- m.o/pa.o
  
  if.psi8 <- (gn == 0)/p.1o * (m-(a==1)*psi8.hat)
  
  ## influence-function of ETT
  if.ett <- if.psi1.11-if.psi1.01-if.psi1.10+if.psi1.00+
    if.psi2.10-if.psi2.00-
    if.psi3-if.psi5+if.psi8
  
  ## plug-in estimator of ETT
  pi.est.ETTb <- psi1.11.hat-psi1.01.hat-psi1.10.hat+psi1.00.hat+
    psi2.10.hat-psi2.00.hat-
    psi3.hat-psi5.hat+psi8.hat
  
  ## influence-function-based estimator of ETT
  if.est.ETTb  <- pi.est.ETTb + mean(if.ett, na.rm = T)
  
  # 95% CI
  ### get standard derivation of each fold
  var.vec <- rep(NA,k)
  for(fold in 1:crossfit_folds){
    fold_idx <- cf_inds[[fold]]$eval
    var.vec[fold] <- var(if.ett[fold_idx]-if.est.ETTb, na.rm = T)
  }
  if.sd.ETTb <- sqrt(mean(var.vec, na.rm = T))/sqrt(n-length(idx.rm))
  lower <- if.est.ETTb  - critical_t * if.sd.ETTb
  upper <- if.est.ETTb  + critical_t * if.sd.ETTb
  
  result <- list(pi.est.ETTb = pi.est.ETTb,
                 if.est.ETTb = if.est.ETTb,
                 if.sd.ETTb = if.sd.ETTb,
                 if.lower.ETTb = lower,
                 if.upper.ETTb = upper,
                 n.rows.removed = length(idx.rm))
  
  return(result)
}
