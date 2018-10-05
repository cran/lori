#' covmat
#'
#' @param R nxK1 matrix of row covariates
#' @param C nxK2 matrix of column covariates
#' @param n number of rows
#' @param p number ofcolumns
#'
#' @return the joint product of R and C, a (np)x(K1+K1) matrix in order row1col1,row2col1,...,rowncol1, row1col2, row2col2,...,rowncolp
#' @export
#'
#' @examples
#' R <- matrix(rnorm(10), 5)
#' C <- matrix(rnorm(9), 3)
#' covs <- covmat(R,C,5,3)
covmat <- function(R,C,n,p){
  if(is.null(R)) {
    covs <- covmatC(C,n)
  } else if(is.null(C)) {
    covs <- covmatR(R,p)
  } else {
    R <- as.matrix(R)
    C <- as.matrix(C)
    dR <- dim(R)
    dC <- dim(C)
    K1 <- dR[2]
    K2 <- dC[2]
    covs <- cbind(do.call(rbind, replicate(nrow(C), R, simplify=FALSE)),
                  C[rep(seq_len(nrow(C)), each=nrow(R)),])
  }
  return(covs)
}
#' covmatR
#'
#' @param R nxK1 matrix of row covariates
#' @param p number ofcolumns
#'
#' @return repeats every row of R p times
#' @export
#'
#' @examples
#' R <- matrix(rnorm(10), 5)
#' cov <- covmatR(R,3)
covmatR <- function(R, p){
  R <- as.matrix(R)
  dR <- dim(R)
  n <- dR[1]
  K1 <- dR[2]
  covs <- do.call(rbind, replicate(p, R, simplify=FALSE))
  return(covs)
}

#' covmatC
#'
#' @param C nxK2 matrix of column covariates
#' @param n number of rows
#'
#' @return repeats C n times
#' @export
#'
#' @examples
#' C <- matrix(rnorm(10), 5)
#' cov <- covmatC(C,3)
covmatC <- function(C, n){
  C <- as.matrix(C)
  dC <- dim(C)
  p <- dC[1]
  K2 <- dC[2]
  covs <- C[rep(seq_len(p), each=n),]
  return(covs)
}

#' Automatic selection of nuclear norm regularization parameter
#' @param Y A matrix of counts (contingency table).
#' @param covmat A (np)xK matrix of K covariates about rows and columns
#' @param lambda2 A positive number, the regularization parameter for covariates main effects
#' @param q A number between \code{0} and \code{1}. The quantile of the distribution of $lambda_{QUT}$ to take.
#' @param N An integer. The number of parametric bootstrap samples to draw.
#' @return the value of $lambda_{QUT}$ to use in LoRI.
#' @export
#' @examples
#' X = matrix(rnorm(30), 15)
#' Y = matrix(rpois(15, 1:15), 5)
#' lambda = qut(Y,X, 10, N=10)
qut <- function(Y, covmat, lambda2 = 0, q = 0.95, N = 100){
  d <- dim(Y)
  n <- d[1]
  p <- d[2]
  nullest <- null_model(Y, lambda2, covmat, trace.it = F)
  X0 <- nullest$X
  lambdas <- rep(0, N)
  m <- sum(!is.na(Y))
  for(i in 1:N)
  {
    Ysimul <- matrix(stats::rpois(n*p, exp(c(X0))), nrow = n)
    nullest2 <- null_model(Ysimul, lambda2, covmat)
    X0simul <- nullest2$X
    dat <- sweep(Ysimul-exp(X0simul), 2, colMeans(Ysimul-exp(X0simul)))
    dat <- sweep(Ysimul-exp(X0simul), 1, rowMeans(Ysimul-exp(X0simul)))
    lambdas[i]<- (1 /m)*svd::propack.svd(dat, neig = 1, opts = list(maxiter = 1e5))$d
    cat('\r',i,"/",N)
  }
  return(stats::quantile(lambdas, q)[[1]])
}
