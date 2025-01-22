# Load necessary library
library(MASS)

source("DataFusion/datafusion_dataset.R")
source("DataFusion/datafusion_dgp_fun.R")
source("DataFusion/datafusion_true.R")
datafusion_true <- readRDS(file = "datafusion_true.rds")


# Define the generate_data function
generate_data <- function(n) {
  X1 <- rnorm(n,0.1,1)
  X2 <- rnorm(n,-0.1,1)
  B <- rbinom(n,1,p.b.fun(X1,X2))
  G <- rbinom(n,1,p.g.fun(X1,X2,B)) #G=1 -- exp; G=0 -- obs
  U <- rnorm(n, mu.u.fun(X1,X2,B), 1)
  p.A <- G*p.a.e.fun(X1,X2,B)+(1-G)*p.a.o.fun(U,X1,X2,B)
  A <- rbinom(n,1,p.A)
  Z <- rnorm(n,mu.z.fun(U,X1,X2,B,A),1)
  
  # generate fake covariates to test robustness
  X1s <- ifelse(X1 <= -0.2, -1, ifelse(X1>=0.3,1,0))
  X2s <- ifelse(X2 <= -0.3, 3, ifelse(X2>=0.2,1,2))
  Bs <- B+rnorm(n,0,1)
  
  # Generate outcomes
  M0 <- rnorm(n,mu.m.fun(U,X1,X2,0,B),1)
  M1 <- rnorm(n,mu.m.fun(U,X1,X2,1,B),1)
  M <- (1-A)*M0 + A * M1
  
  a.xb.exp <- rep(NA,n)
  u.abx.exp <- rep(NA,n)
  m.abx.exp <- rep(NA,n)
  for(i in 1:n){
    
    if(G[i]==0){
      # get P(A=1|X,B) (or E[A|X,B])
      p.a.xb <- function(u0){
        p.a.o.fun(u0,X1[i],X2[i],B[i])*dnorm(u0,mu.u.fun(X1[i],X2[i],B[i]),1)
      }
      a.xb.exp[i] <- sum(p.a.xb(seq(-8,8,by=0.04)))*0.04
      
      # get E[U|A,B,X]
      p.u.abx <- function(u0){
        u0*((p.a.o.fun(u0,X1[i],X2[i],B[i])*A[i]+(1-p.a.o.fun(u0,X1[i],X2[i],B[i]))*(1-A[i]))*dnorm(u0,mu.u.fun(X1[i],X2[i],B[i]),1))/(A[i]*a.xb.exp[i]+(1-A[i])*(1-a.xb.exp[i]))
      }
      u.abx.exp[i] <- sum(p.u.abx(seq(-8,8,by=0.04)))*0.04
    }else{
      # get P(A=1|X,B) (or E[A|X,B])
      a.xb.exp[i] <- p.a.e.fun(X1[i],X2[i],B[i])
      # get E[U|A,B,X] (which is just get E[U|B,X])
      u.abx.exp[i] <- mu.u.fun(X1[i],X2[i],B[i])
    }
    
    m.abx.exp[i] <- mu.m.fun(u.abx.exp[i],X1[i],X2[i],A[i],B[i])
  }
  
  
  W <- beta1.fun(X1,X2,B)*A+gamma0.fun(X1,X2)*(A-a.xb.exp)+ym0.fun(X1,X2)
  Y <- rnorm(n,m.abx.exp+W,1)
  
  # mean(Y[G == 0 & A == 1])-mean(Y[G == 0 & A == 0])
  
  X <- cbind(X1,X2,B)
  colnames(X) <- c("X1","X2","X3")
  
  Xs <- cbind(X1s,X2s,Bs)
  colnames(Xs) <- c("X1s","X2s","X3s")
  
  dataset <- ProxCIData$new(X = X, Z = Z, M = M, Y = Y, A = A, G = G, Xs = Xs)
  return(list(dataset = dataset, target = datafusion_true))
}
