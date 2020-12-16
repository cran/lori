## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(lori)
set.seed(123)

## ----simulate covariates------------------------------------------------------
## covariates
n <- 10 # number of rows
p <- 6 # number of columns
K1 <- 2 # number of row covariates
K2 <- 2 # number of column covariates
K3 <- 3 # number of (rowxcolumn) covariates
R <- matrix(rnorm(n*K1), nrow=n) # matrix of row covariates
C <- matrix(rnorm(p*K2), nrow=p) # matrix of column covariates
E <- matrix(rnorm(n*p*K3), nrow=n*p) # matrix of  (rowxcolumn) covariates
U <- covmat(n, p, R, C, E)
U <- scale(U)


## ----eval=FALSE---------------------------------------------------------------
#  U <- covmat(n, p, R)
#  U <- covmat(n, p, C)
#  U <- covmat(n, p, R, C)
#  U <- covmat(n, p, R=R, E=E)
#  U <- covmat(n, p, C=C, E=E)

## ----simulate parameters------------------------------------------------------
## parameters
alpha0 <- rep(0, n)
alpha0[1:6] <- 1
beta0 <- rep(0, p)
beta0[1:4] <- 1
epsilon0 <- rep(0, K1+K2+K3)
epsilon0[5:6] <- 0.2
r <- 2 #rank of interaction matrix theta0
theta0 <- 0.1*matrix(rnorm(n*r), nrow=n)%*%diag(r:1)%*%matrix(rnorm(p*r), nrow=r)
theta0 <- sweep(theta0, 2, colMeans(theta0))
theta0 <- sweep(theta0, 1, rowMeans(theta0))


## ----intensities and data-----------------------------------------------------
## construct x0
x0 <- matrix(rep(alpha0, p), nrow=n) #row effects
x0 <- x0 + matrix(rep(beta0, each=n), nrow=n) #add column effects
x0 <- x0 + matrix(U%*% epsilon0, nrow=n) #add cov effects
x0 <- x0 + theta0 #add interactions

## sample count data y
y0 <- matrix(rpois(n*p, lambda = c(exp(x0))), nrow = n)

## add missing values
p_miss <- 0.2
y <- y0
y[sample(1:(n*p), round(p_miss*n*p))] <- NA

## ----estimation---------------------------------------------------------------
## lori estimation
lambda1 <- 0.1
lambda2 <- 0.1
m <- sum(!is.na(y))
t <- Sys.time()
res.lori <- lori(y, U, 0.1, 0.1)
t <- Sys.time()-t

## ----cross-validation---------------------------------------------------------
res.cv <- cv.lori(y, U, trace.it = F, len=5)
res.lori <- lori(y, U, res.cv$lambda1, res.cv$lambda2)

## -----------------------------------------------------------------------------
res.lori$alpha
res.lori$beta
res.lori$espilon
res.lori$theta


## -----------------------------------------------------------------------------
res.lori$imputed

## ----multiple imputation------------------------------------------------------
res.mi <- mi.lori(y, U, res.cv$lambda1, res.cv$lambda2)
res.pool <- pool.lori(res.mi)
boxplot(res.mi$mi.beta, pch="", names=paste("col", 1:p))

