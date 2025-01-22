# Load necessary library
library(MASS)

source("DataFusion/datafusion_dgp_fun.R")

set.seed(12345)
n <- 1e8
# Generate U, X
X1 <- rnorm(n,0.1,1)
X2 <- rnorm(n,-0.1,1)
# X <- cbind(X1,X2)
# colnames(X) <- c("X1","X2")
B <- rbinom(n,1,p.b.fun(X1,X2))
G <- rbinom(n,1,p.g.fun(X1,X2,B)) #G=1 -- exp; G=0 -- obs
U <- rnorm(n, mu.u.fun(X1,X2,B), 1)
p.A <- G*p.a.e.fun(X1,X2,B)+(1-G)*p.a.o.fun(U,X1,X2,B)
A <- rbinom(n,1,p.A)

# Generate outcomes
M0 <- rnorm(n,mu.m.fun(U,X1,X2,0,B),1)
M1 <- rnorm(n,mu.m.fun(U,X1,X2,1,B),1)
M <- (1-A)*M0 + A * M1

ett.true <- mean(beta1.fun(X1,X2,B)[G == 0 & A == 1]) + mean(M1[G == 0 & A == 1]-M0[G == 0 & A == 1])
# ett.true
# [1] 0.7005945


datafusion_true <- list(ett.true = ett.true)
saveRDS(datafusion_true, file = "datafusion_true.rds")
