ECLSK8_complete <- readRDS("DataFusion/ECLSK8_complete.rds")
STAR_complete <- readRDS("DataFusion/STAR_complete.rds")

#A assigned does not necessarily match A receive (96 unmatched records)

### STAR
# A_assign.e <- ifelse(STAR_complete$g3classtype==unique(STAR_complete$g3classtype)[3],1,0)
A.e <- ifelse(STAR_complete$g3classsize<=19,1,0)
X1.e <- ifelse(STAR_complete$gender=="MALE",1,0)
X2.e <- ifelse(STAR_complete$race=="WHITE",1,0)
B.e <- ifelse(STAR_complete$g3freelunch=="FREE LUNCH",1,0)
Z.e <- ifelse(STAR_complete$g3surban=="RURAL",3,ifelse(STAR_complete$g3surban=="SUBURBAN",2,1)) #(NOT NEEDED)

M_math.e <- STAR_complete$g3tmathss

M_math.fac.e <- ifelse(M_math.e <= 580,1,ifelse(M_math.e <= 590,2,ifelse(M_math.e <= 600,3,ifelse(M_math.e <= 610, 4,ifelse(M_math.e <= 620 , 5,6)))))

Y_math.e <- STAR_complete$g8tmathss #NULL

### ECLS
A.o <- ifelse(ECLSK8_complete$A5TOTRA<=19,1,0)
X1.o <- ifelse(ECLSK8_complete$GENDER==1,1,0)
X2.o <- ifelse(ECLSK8_complete$RACE==1,1,0)
B.o <- ifelse(ECLSK8_complete$P5LUNCHS==1,1,0)
Z.o <- ECLSK8_complete$R5URBAN

M_math.o <- ECLSK8_complete$C5R4MTHT_R

M_math.fac.o <- ifelse(M_math.o <= 580,1,ifelse(M_math.o <= 590,2,ifelse(M_math.o <= 600,3,ifelse(M_math.o <= 610, 4,ifelse(M_math.o <= 620 , 5, 6)))))

Y_math.o <- ECLSK8_complete$C7R4MTHT_R

### Pooled
G <- c(rep(1,length(A.e)),rep(0,length(A.o)))
A <- c(A.e,A.o)
X1 <- c(X1.e,X1.o)
X2 <- c(X2.e,X2.o)
Z <- c(Z.e,Z.o)
B <- c(B.e,B.o)
M_math <- c(M_math.e,M_math.o)
M_math.fac <- c(M_math.fac.e,M_math.fac.o)
Y_math <- c(rep(0,length(A.e)),Y_math.o)

dat <- data.frame(G=G,A=A,X1=X1,X2=X2,Z=Z,B=B,M=M_math,M.fac=M_math.fac,Y=Y_math)
dat <- dat[dat$M >= 570 & dat$M <= 630,]

saveRDS(dat,"dat.rds")

