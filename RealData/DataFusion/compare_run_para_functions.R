run_para_Athey <- function(dataset,used_model,crossfit_folds=4,Athey_cf_inds){
  
  n <- dataset$n
  
  y <- dataset$Y
  m <- dataset$M
  x1 <- dataset$X[,1]
  x2 <- dataset$X[,2]
  x3 <- dataset$Z
  x4 <- dataset$X[,3]
  a <- dataset$A
  gn <- dataset$G
  
  # Create dataframe
  df <- data.frame(y = y, m = m, x1 = x1, x2 = x2, x3 = x3, x4 = x4, a = a, gn = gn)
  colnames(df) <- c("y", "m", "x1", "x2", "x3", "x4", "a", "gn")
  
  # Estimation
  ## P(A=1|G=O)
  pa.o <- nrow(df[a == 1 & gn == 0,])/nrow(df[gn == 0,])
  
  ## E(Y|A=1,G=O)
  y.1o <- mean(df$y[a==1 & gn == 0])
  
  ## E(Y|A=0,G=O)
  y.0o <- mean(df$y[a==0 & gn == 0])
  
  cf_inds <- Athey_cf_inds
  y.m0o.e.o.hat <- rep(NA,n)
  for(fold in 1:crossfit_folds){
    
    fold_train_idx <- cf_inds[[fold]]$train
    fold_idx <- cf_inds[[fold]]$eval
    
    df.train <- df[fold_train_idx,]
    
    ## E(Y|M,A=0,X,G=O)
    y.m0o.fit <- glm(used_model$y.m0o, data = df.train[df.train$gn == 0 & df.train$a==0,], family = "gaussian")
    y.m0o.hat <- predict(y.m0o.fit,newdata = df.train[df.train$gn == 1,])
    
    y.m0o.e.fit <- glm(y.m0o.hat ~ x1*x2*x3*x4*I(x3^2)*I(x4^2), df.train[df.train$gn == 1,], family = "gaussian")
    y.m0o.e.hat <- predict(y.m0o.e.fit,newdata = df.train[df.train$gn == 0,])
    
    y.m0o.e.o.fit <- glm(y.m0o.e.hat ~ x1*x2*x3*x4*I(x3^2)*I(x4^2), df.train[df.train$gn == 0,], family = "gaussian")
    y.m0o.e.o.hat[fold_idx] <- predict(y.m0o.e.o.fit, newdata=df[fold_idx,])
    }
  
  y.m0o.e.o <- y.m0o.e.o.hat[gn==0]
  Athey.est.ETT  <- y.1o-mean(y.m0o.e.o)/pa.o+(1-pa.o)/pa.o*y.0o
  
  result <- list(Athey.est.ETT = Athey.est.ETT)
  
  return(result)
}

run_para_naive <- function(dataset,used_model,crossfit_folds=4,naive_cf_inds){
  
  n <- dataset$n
  
  y <- dataset$Y
  m <- dataset$M
  x1 <- dataset$X[,1]
  x2 <- dataset$X[,2]
  x3 <- dataset$Z
  x4 <- dataset$X[,3]
  a <- dataset$A
  gn <- dataset$G
  
  # Create dataframe
  df <- data.frame(y = y, m = m, x1 = x1, x2 = x2, x3 = x3, x4 = x4, a = a, gn = gn)
  colnames(df) <- c("y", "m", "x1", "x2", "x3", "x4", "a", "gn")
  
  # Estimation
  ## E(Y|A=1,G=O)
  y.1o <- mean(df$y[a==1 & gn == 0])
  y.0o <- mean(df$y[a==0 & gn == 0])
  m.0o <- mean(df$m[a==0 & gn == 0])
  p.1o <- nrow(df[a == 1 & gn == 0,])/nrow(df)
  
  cf_inds <- naive_cf_inds
  y.0o.hat <- rep(NA,n)
  m.0e.hat <- rep(NA,n)
  m.0o.hat <- rep(NA,n)
  a.e.hat <- rep(NA,n)
  a.o.hat <- rep(NA,n)
  g.hat <- rep(NA,n)
  for(fold in 1:crossfit_folds){
    
    fold_train_idx <- cf_inds[[fold]]$train
    fold_idx <- cf_inds[[fold]]$eval
    
    df.train <- df[fold_train_idx,]
    
    ## E(Y|A=0,X,G=O)
    y.0o.fit <- glm(used_model$y.0o, data = df.train[df.train$gn == 0 & df.train$a==0,], family = "gaussian")
    y.0o.hat[fold_idx] <- predict(y.0o.fit, newdata=df[fold_idx,])
    
    ## E(M|A=0,X,G=E)
    m.0e.fit <- glm(used_model$m.0e, data = df.train[df.train$gn == 1 & df.train$a==0,], family = "gaussian")
    m.0e.hat[fold_idx] <- predict(m.0e.fit, newdata=df[fold_idx,])
    
    ## E(M|A=0,X,G=O)
    m.0o.fit <- glm(used_model$m.0o, data = df.train[df.train$gn == 0 & df.train$a==0,], family = "gaussian")
    m.0o.hat[fold_idx] <- predict(m.0o.fit, newdata=df[fold_idx,])
    
    ## P(A=1|X,G=O)
    a.o.fit <- glm(used_model$a.o, data = df.train[df.train$gn == 0,], family = "binomial")
    a.o.hat[fold_idx] <- expit(predict(a.o.fit, newdata=df[fold_idx,]))
    
    ## P(A=1|X,G=E)
    a.e.fit <- glm(used_model$a.e, data = df.train[df.train$gn == 1,], family = "binomial")
    a.e.hat[fold_idx] <- expit(predict(a.e.fit, newdata=df[fold_idx,]))
    
    ## P(G=E|X)
    g.fit <- glm(used_model$g, data = df.train, family = "binomial")
    g.hat[fold_idx] <- expit(predict(g.fit, newdata=df[fold_idx,]))
  }
  
  y0o.1o <- y.0o.hat[a==1 & gn==0]
  m0e.1o <- m.0e.hat[a==1 & gn==0]
  
  equi.est.ETT <- y.1o-y.0o+m.0o-mean(m0e.1o)
  
  naive.est.ETT  <- y.1o-mean(y0o.1o)
  
  psi1 <- ((a==1)*(gn==0)/p.1o)*(y-y.1o)
  psi2 <- ((a==0)*(gn==0)/p.1o)*(a.o.hat/(1-a.o.hat))*(y-y.0o.hat)+((a==1)*(gn==0)/p.1o)*(y.0o.hat-mean(y0o.1o))
  
  psi3 <- ((a==0)*(gn==0)/p.1o)*(1/(1-a.o.hat))*(m-m.0o.hat)
  psi4 <- ((a==0)*(gn==1)/p.1o)*(1/g.hat-1)*(m-m.0e.hat)
  psi5 <- ((a==0)*(gn==0)/p.1o)*(a.o.hat/(1-a.o.hat))*(y-y.0o.hat)
  psi6 <- ((gn==0)/p.1o)*(m.0o.hat-m.0e.hat+(a==1)*(y-y.0o.hat))
  
  naive.if.ETT <- psi1-psi2
  naive.if.est.ETT <- mean(naive.if.ETT)+y.1o-mean(y0o.1o)
  
  equi.if.ETT <- psi3+psi4+psi5+psi6
  equi.if.est.ETT <- mean(equi.if.ETT)
  
  result <- list(naive.est.ETT = naive.est.ETT,
                 naive.if.est.ETT = naive.if.est.ETT,
                 equi.est.ETT = equi.est.ETT,
                 equi.if.est.ETT = equi.if.est.ETT)
  
  return(result)
}
