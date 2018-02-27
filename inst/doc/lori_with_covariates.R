## ----load-data-----------------------------------------------------------
library(ade4)
data(aravo)
library(lori)
library(FactoMineR)


# Reduce size of data so that code runs quicker (remove lines to perform full experiment -> 45')
Y <- aravo$spe[1:35, 1:40]
R <- aravo$env[1:35, ]
C <- aravo$traits[1:40, ]
data(aravo)
m1 <- 35
m2 <- 40

# Code categorical variables as indicators
R$ZoogDno <- rep(0,m1)
R$ZoogDsome <- rep(0,m1)
R$ZoogDhigh <- rep(0,m1)
R$ZoogDno[R$ZoogD=="no"] <- 1
R$ZoogDsome[R$ZoogD=="some"] <- 1
R$ZoogDhigh[R$ZoogD=="high"] <- 1

R$Form1 <- rep(0,m1)
R$Form2 <- rep(0,m1)
R$Form3 <- rep(0,m1)
R$Form4 <- rep(0,m1)
R$Form5 <- rep(0,m1)
R$Form1[R$Form==1] <- 1
R$Form2[R$Form==2] <- 1
R$Form3[R$Form==3] <- 1
R$Form4[R$Form==4] <- 1
R$Form5[R$Form==5] <- 1

R <- R[, c(1,2,4,6:14)]

## ----qut-rc--------------------------------------------------------------
# Run LORI algorithm without using covariates 
qut_no_covariates = lori(Y)

## ----plots-qut-rc--------------------------------------------------------

plot_interaction(qut_no_covariates$Theta, selectRow = "contrib 10", selectCol = "contrib 10")
plot_interaction(cbind(qut_no_covariates$Theta, R), quanti.sup = 41:52, choix = "quanti.sup")


## ----qut-covariates------------------------------------------------------
# Define projection on aravo covariates
Xr <- R
Xr <- qr(Xr)
Xr <- qr.Q(Xr)
Xc <- C
Xc <- qr(Xc)
Xc <- qr.Q(Xc)

aravo_projection=function(X){
  X = as.matrix(X)
  p = function(X, Xr, Xc){
  return(X - Xr%*%t(Xr)%*%X - X%*%Xc%*%t(Xc) + Xr%*%t(Xr)%*%X%*%Xc%*%t(Xc))
}
  return(p(X, aravo$R, aravo$C))
}


# Run algorithm with QUT using the known covariates 
#------------------------------------------------------------------------------------
# this takes a while with entire data set (45 minutes)
qut_covariates = lori(Y, R, C)

plot_interaction(qut_no_covariates$Theta, selectRow = "contrib 10", selectCol = "contrib 10")
plot_interaction(cbind(qut_no_covariates$Theta, R), quanti.sup = 41:52, choix = "quanti.sup")


