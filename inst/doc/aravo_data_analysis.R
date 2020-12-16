## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(lori)
set.seed(123)

## ----load the data------------------------------------------------------------
data(aravo)
Y <- aravo$spe
R <- aravo$env
R <- R[, c(1,2,4,6)]
C <- aravo$traits
d <- dim(Y)
n <- d[1]
p <- d[2]
U <- covmat(n,p,R,C)
U <- scale(U)

## -----------------------------------------------------------------------------
# Tune regularization parameter
res_cv <- cv.lori(Y, U, reff=F, ceff=F, trace.it=F, len=5) 
res_lori <- lori(Y, U, lambda1 = res_cv$lambda1, lambda2=res_cv$lambda2, reff=F, ceff=F)

## -----------------------------------------------------------------------------
# multiple imputation
res_mi <- mi.lori(Y, U, lambda1 = res_cv$lambda1, lambda2=res_cv$lambda2, reff=F, ceff=F, M=20)
boxplot(res_mi$mi.epsilon)

