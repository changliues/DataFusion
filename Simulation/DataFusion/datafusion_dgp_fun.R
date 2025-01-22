expit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}


p.b.fun <- function(x1,x2){
  expit(-0.43+0.15*x1+0.18*x2)
}

mu.u.fun <- function(x1,x2,b){
  0.15-0.35*x1+0.8*x2+0.6*b
}

p.g.fun <- function(x1,x2,b){
  expit(0.2+0.15*x1+0.1*x2-0.35*b)
}

p.a.o.fun <- function(u,x1,x2,b){
  expit(-0.3-0.1*u+1.3*b+0.1*x1+0.15*x2)
}

p.a.e.fun <- function(x1,x2,b){
  expit(-0.23+0.68*b-0.13*x1)
}

mu.z.fun <- function(u,x1,x2,b,a){
  0.2+1.4*u+1.5*a+0.1*x1-0.5*x2+1.3*b
}

mu.m.fun <- function(u,x1,x2,a,b){
  0.75*u+0.2*x1-0.5*x2+0.4*a+0.12*b
}

beta1.fun <- function(x1,x2,b){
  0.25+0.3*x1-0.6*x2-0.1*b
}

gamma0.fun <- function(x1,x2){
  0.12-0.5*x1+0.45*x2
}

ym0.fun <- function(x1,x2){
  0.54-0.28*x1+0.35*x2
}