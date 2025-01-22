run_para_BSIV <- function(dataset,used_model,crossfit_folds=4,BSIV_cf_inds){
  
  n <- dataset$n
  
  y <- dataset$Y
  m <- dataset$M
  x1 <- dataset$X[,1]
  x2 <- dataset$X[,2]
  b <- dataset$X[,3]
  a <- dataset$A
  gn <- dataset$G
  x1s <- dataset$Xs[,1]
  x2s <- dataset$Xs[,2]
  
  # Create dataframe
  df <- data.frame(y = y, m = m, x1 = x1, x2 = x2, b = b, a = a, gn = gn, x1s = x1s, x2s = x2s)
  colnames(df) <- c("y", "m", "x1", "x2", "b", "a", "gn", "x1s", "x2s")
  
  # Estimation
  ## P(A=1,G=O)
  p.1o <- nrow(df[a == 1 & gn == 0,])/nrow(df)
  
  ## P(A=1|G=O)
  pa.o <- nrow(df[a == 1 & gn == 0,])/nrow(df[gn == 0,])
  
  ## E(M|G=O)
  m.o <- mean(df$m[gn == 0])
  
  cf_inds <- BSIV_cf_inds
  g.xb.hat <- rep(NA,n)
  b.hat <- rep(NA,n)
  b.o.hat <- rep(NA,n)
  b.e.hat <- rep(NA,n)
  a.be.hat <- rep(NA,n)
  a.bo.hat <- rep(NA,n)
  a.1o.hat <- rep(NA,n)
  a.0o.hat <- rep(NA,n)
  m.11o.hat <- rep(NA,n)
  m.01o.hat <- rep(NA,n)
  m.10o.hat <- rep(NA,n)
  m.00o.hat <- rep(NA,n)
  m.0be.hat <- rep(NA,n)
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
    g.xb.fit <- glm(used_model$g.xb, data = df[fold_train_idx,], family = "binomial")
    g.xb.hat[fold_idx] <- expit(predict(g.xb.fit, 
                                        newdata = data.frame(x1=x1[fold_idx],
                                                             x1s=x1[fold_idx],
                                                             x2=x2[fold_idx],
                                                             x2s=x2[fold_idx],
                                                             b=b[fold_idx])))
    
    ## P(B=1|X,G=O)
    b.o.fit <- glm(used_model$b.o, data = df.train[df.train$gn == 0,], family = "binomial")
    b.o.hat[fold_idx] <- expit(predict(b.o.fit,
                                       newdata = data.frame(x1=x1[fold_idx],
                                                            x1s=x1[fold_idx],
                                                            x2=x2[fold_idx],
                                                            x2s=x2[fold_idx])))
    
    ## P(B=1|X,G=E)
    b.e.fit <- glm(used_model$b.e, data = df.train[df.train$gn == 1,], family = "binomial")
    b.e.hat[fold_idx] <- expit(predict(b.e.fit,
                                       newdata = data.frame(x1=x1[fold_idx],
                                                            x1s=x1[fold_idx],
                                                            x2=x2[fold_idx],
                                                            x2s=x2[fold_idx])))
    
    ## P(A=1|X,B,G=E)
    a.be.fit <- glm(used_model$a.be, data = df.train[df.train$gn == 1,], family = "binomial")
    a.be.hat[fold_idx] <- expit(predict(a.be.fit,
                                        newdata = data.frame(x1=x1[fold_idx],
                                                             x1s=x1[fold_idx],
                                                             x2=x2[fold_idx],
                                                             x2s=x2[fold_idx],
                                                             b=b[fold_idx])))
    
    ## P(A=1|X,B,G=O)
    a.bo.fit <- glm(used_model$a.bo, data = df.train[df.train$gn == 0,], family = "binomial")
    a.bo.hat[fold_idx] <- expit(predict(a.bo.fit, 
                                        newdata = data.frame(x1=x1[fold_idx],
                                                             x1s=x1[fold_idx],
                                                             x2=x2[fold_idx],
                                                             x2s=x2[fold_idx],
                                                             b=b[fold_idx])))
    
    ## P(A=1|B=1,X,G=O)
    a.1o.hat[fold_idx] <- expit(predict(a.bo.fit, 
                                        newdata = data.frame(x1=x1[fold_idx],
                                                             x1s=x1[fold_idx],
                                                             x2=x2[fold_idx],
                                                             x2s=x2[fold_idx],
                                                             b=1)))
    
    ## P(A=1|B=0,X,G=O)
    a.0o.hat[fold_idx] <- expit(predict(a.bo.fit,
                                        newdata = data.frame(x1=x1[fold_idx],
                                                             x1s=x1[fold_idx],
                                                             x2=x2[fold_idx],
                                                             x2s=x2[fold_idx],
                                                             b=0)))
    
    m.abo.fit <- glm(used_model$m.abo, data = df.train[df.train$gn == 0,])
    y.abo.fit <- glm(used_model$y.abo, data = df.train[df.train$gn == 0,])
    ## E(M|A=1,B=1,X,G=O)
    m.11o.hat[fold_idx] <- predict(m.abo.fit,
                                   newdata = data.frame(x1=x1[fold_idx],
                                                        x1s=x1[fold_idx],
                                                        x2=x2[fold_idx],
                                                        x2s=x2[fold_idx],
                                                        b=1,a=1))
    
    ## E(Y|A=1,B=1,X,G=O)
    y.11o.hat[fold_idx] <- predict(y.abo.fit,
                                   newdata = data.frame(x1=x1[fold_idx],
                                                        x1s=x1[fold_idx],
                                                        x2=x2[fold_idx],
                                                        x2s=x2[fold_idx],
                                                        b=1,a=1))
    
    ## E(M|A=1,B=0,X,G=O)
    m.10o.hat[fold_idx] <- predict(m.abo.fit,
                                   newdata = data.frame(x1=x1[fold_idx],
                                                        x1s=x1[fold_idx],
                                                        x2=x2[fold_idx],
                                                        x2s=x2[fold_idx],
                                                        b=0,a=1))
    
    ## E(Y|A=1,B=0,X,G=O)
    y.10o.hat[fold_idx] <- predict(y.abo.fit,
                                   newdata = data.frame(x1=x1[fold_idx],
                                                        x1s=x1[fold_idx],
                                                        x2=x2[fold_idx],
                                                        x2s=x2[fold_idx],
                                                        b=0,a=1))
    
    ## E(M|A=0,B=1,X,G=O)
    m.01o.hat[fold_idx] <- predict(m.abo.fit,
                                   newdata = data.frame(x1=x1[fold_idx],
                                                        x1s=x1[fold_idx],
                                                        x2=x2[fold_idx],
                                                        x2s=x2[fold_idx],
                                                        b=1,a=0))
    
    ## E(Y|A=0,B=1,X,G=O)
    y.01o.hat[fold_idx] <- predict(y.abo.fit,
                                   newdata = data.frame(x1=x1[fold_idx],
                                                        x1s=x1[fold_idx],
                                                        x2=x2[fold_idx],
                                                        x2s=x2[fold_idx],
                                                        b=1,a=0))
    
    ## E(M|A=0,B=0,X,G=O)
    m.00o.hat[fold_idx] <- predict(m.abo.fit,
                                   newdata = data.frame(x1=x1[fold_idx],
                                                        x1s=x1[fold_idx],
                                                        x2=x2[fold_idx],
                                                        x2s=x2[fold_idx],
                                                        b=0,a=0))
    
    ## E(Y|A=0,B=0,X,G=O)
    y.00o.hat[fold_idx] <- predict(y.abo.fit,
                                   newdata = data.frame(x1=x1[fold_idx],
                                                        x1s=x1[fold_idx],
                                                        x2=x2[fold_idx],
                                                        x2s=x2[fold_idx],
                                                        b=0,a=0))
    
    ## E(M|A=0,X,B,G=E)
    m.abe.fit <- glm(used_model$m.abe, data = df.train[df.train$gn == 1,])
    m.0be.hat[fold_idx] <- predict(m.abe.fit,
                                   newdata = data.frame(x1=x1[fold_idx],
                                                        x1s=x1[fold_idx],
                                                        x2=x2[fold_idx],
                                                        x2s=x2[fold_idx],
                                                        b=b[fold_idx],
                                                        a=0))
  }
  
  idx.rm <- which(abs(a.1o.hat-a.0o.hat) <= sqrt(1/n))
  if(length(idx.rm) > 0){
    
    df[idx.rm,] <- NA
    # y[idx.rm] <- NA
    # m[idx.rm] <- NA
    # x1[idx.rm] <- NA
    # x2[idx.rm] <- NA
    b[idx.rm] <- NA
    a[idx.rm] <- NA
    gn[idx.rm] <- NA
    # x1s[idx.rm] <- NA
    # x2s[idx.rm] <- NA
    
    g.xb.hat[idx.rm] <- NA
    b.o.hat[idx.rm] <- NA
    b.e.hat[idx.rm] <- NA
    a.be.hat[idx.rm] <- NA
    a.bo.hat[idx.rm] <- NA
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
    m.0be.hat[idx.rm] <- NA
  }
  
  # P(A=1|X,G=O)=P(A=1|B=1,X,G=O)P(B=1|X,G=O)+P(A=1|B=0,X,G=O)(1-P(B=1|X,G=O))
  a.o.hat <- a.1o.hat*b.o.hat + a.0o.hat*(1-b.o.hat)
  
  ## influence-function of psi_1^ab
  b.1o <- df$b[a==1 & gn == 0]
  ## influence-function of psi_1^11
  y.11o.hat.1o <- y.11o.hat[a==1 & gn == 0]
  m.11o.hat.1o <- m.11o.hat[a==1 & gn == 0]
  
  psi1.11.hat <- mean((y.11o.hat.1o-m.11o.hat.1o)*b.1o, na.rm = T)
  
  if.psi1.11 <- ((gn == 0)/p.1o) * ((a==1)*(b==1)*(y-m-(y.11o.hat-m.11o.hat))+
                                      (a==1)*((y.11o.hat-m.11o.hat)*b-psi1.11.hat))
  
  ## influence-function of psi_1^01
  y.01o.hat.1o <- y.01o.hat[a==1 & gn == 0]
  m.01o.hat.1o <- m.01o.hat[a==1 & gn == 0]
  
  psi1.01.hat <- mean((y.01o.hat.1o-m.01o.hat.1o)*b.1o, na.rm = T)
  
  if.psi1.01 <- ((gn == 0)/p.1o)*((a==0)*(b==1)*(a.1o.hat/(1-a.1o.hat))*(y-m-(y.01o.hat-m.01o.hat))+
                                    (a==1)*((y.01o.hat-m.01o.hat)*b-psi1.01.hat))
  
  ## influence-function of psi_1^10
  y.10o.hat.1o <- y.10o.hat[a==1 & gn == 0]
  m.10o.hat.1o <- m.10o.hat[a==1 & gn == 0]
  
  psi1.10.hat <- mean((y.10o.hat.1o-m.10o.hat.1o)*b.1o, na.rm = T)
  
  if.psi1.10 <- ((gn == 0)/p.1o)*((a==1)*(b==0)*((b.o.hat*a.1o.hat)/((1-b.o.hat)*a.0o.hat))*(y-m-(y.10o.hat-m.10o.hat))+
                                    (a==1)*((y.10o.hat-m.10o.hat)*b-psi1.10.hat))
  
  ## influence-function of psi_1^00
  y.00o.hat.1o <- y.00o.hat[a==1 & gn == 0]
  m.00o.hat.1o <- m.00o.hat[a==1 & gn == 0]
  
  psi1.00.hat <- mean((y.00o.hat.1o-m.00o.hat.1o)*b.1o, na.rm = T)
  
  if.psi1.00 <- ((gn == 0)/p.1o)*((a==0)*(b==0)*((b.o.hat*a.1o.hat)/((1-b.o.hat)*(1-a.0o.hat)))*(y-m-(y.00o.hat-m.00o.hat))+
                                    (a==1)*((y.00o.hat-m.00o.hat)*b-psi1.00.hat))
  
  ## influence-function of psi_2^ab
  ## influence-function of psi_2^10
  psi2.10.hat <- mean((y.10o.hat.1o-m.10o.hat.1o), na.rm = T)
  
  if.psi2.10 <- ((gn == 0)/p.1o)*(((a==1)*(b==0))/a.0o.hat*(a.o.hat/(1-b.o.hat))*(y-m-(y.10o.hat-m.10o.hat))+
                                    (a==1)*(y.10o.hat-m.10o.hat-psi2.10.hat))
  
  ## influence-function of psi_2^00
  psi2.00.hat <- mean((y.00o.hat.1o-m.00o.hat.1o), na.rm = T)
  
  if.psi2.00 <- ((gn == 0)/p.1o)*(((a==0)*(b==0))/(1-a.0o.hat)*(a.o.hat/(1-b.o.hat))*(y-m-(y.00o.hat-m.00o.hat))+
                                    (a==1)*(y.00o.hat-m.00o.hat-psi2.00.hat))
  
  ## influence-function of psi_3
  a.1o.hat.1o <- a.1o.hat[a==1 & gn == 0]
  a.0o.hat.1o <- a.0o.hat[a==1 & gn == 0]
  
  psi3.hat <- mean(((y.01o.hat.1o-m.01o.hat.1o)-(y.00o.hat.1o-m.00o.hat.1o))/((1-a.1o.hat.1o)-(1-a.0o.hat.1o)), na.rm = T)
  
  if.psi3 <- ((gn == 0)/((1-a.1o.hat)-(1-a.0o.hat))) * (a.o.hat/p.1o) *
    ( (a==0)*((b==1)/((1-a.1o.hat)*b.o.hat)*(y-m-(y.01o.hat-m.01o.hat)) -
                (b==0)/((1-a.0o.hat)*(1-b.o.hat))*(y-m-(y.00o.hat-m.00o.hat))) +
        ((y.01o.hat-m.01o.hat)-(y.00o.hat-m.00o.hat))/((1-a.1o.hat)-(1-a.0o.hat)) *
        (-(b==1)/b.o.hat*((a==0)-(1-a.1o.hat))+(b==0)/(1-b.o.hat)*((a==0)-(1-a.0o.hat)))) +
    ((a==1)*(gn == 0)/p.1o)*(((y.01o.hat-m.01o.hat)-(y.00o.hat-m.00o.hat))/((1-a.1o.hat)-(1-a.0o.hat))-psi3.hat)
  
  ## influence-function of psi_5
  m.0be.hat.1o <- m.0be.hat[a==1 & gn == 0]
  a.bo.hat.1o <- a.bo.hat[a==1 & gn == 0]
  
  psi5.hat <- mean(m.0be.hat.1o/a.bo.hat.1o, na.rm = T)
  
  if.psi5 <- 1/p.1o * (((a==0)*(gn == 1)/(1-a.be.hat))*(1/g.xb.hat-1)*(m-m.0be.hat) +
                         (gn == 0)*m.0be.hat - (a==1)*(gn == 0)*psi5.hat)
  
  
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
    var.vec[fold] <- var(if.ett[fold_idx], na.rm = T)
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
