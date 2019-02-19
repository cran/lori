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
covmat <- function(R, C, n, p) {
  if (is.null(R)) {
    covs <- covmatC(C, n)
  } else if (is.null(C)) {
    covs <- covmatR(R, p)
  } else {
    R <- as.matrix(R)
    C <- as.matrix(C)
    dR <- dim(R)
    dC <- dim(C)
    K1 <- dR[2]
    K2 <- dC[2]
    covs <-
      cbind(do.call(rbind, replicate(nrow(C), R, simplify = FALSE)),
            C[rep(seq_len(nrow(C)), each = nrow(R)), ])
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
covmatR <- function(R, p) {
  R <- as.matrix(R)
  dR <- dim(R)
  n <- dR[1]
  K1 <- dR[2]
  covs <- do.call(rbind, replicate(p, R, simplify = FALSE))
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
covmatC <- function(C, n) {
  C <- as.matrix(C)
  dC <- dim(C)
  p <- dC[1]
  K2 <- dC[2]
  covs <- C[rep(seq_len(p), each = n), ]
  return(covs)
}

#' Automatic selection of nuclear norm regularization parameter
#' @param Y A matrix of counts (contingency table).
#' @param cov A (np)xK matrix of K covariates about rows and columns
#' @param lambda2 A positive number, the regularization parameter for covariates main effects
#' @param q A number between \code{0} and \code{1}. The quantile of the distribution of $lambda_{QUT}$ to take.
#' @param N An integer. The number of parametric bootstrap samples to draw.
#' @return the value of $lambda_{QUT}$ to use in LoRI.
#' @export
#' @examples
#' X = matrix(rnorm(30), 15)
#' Y = matrix(rpois(15, 1:15), 5)
#' lambda = qut(Y,X, 10, N=10)
qut <- function(Y,
                cov,
                lambda2 = 0,
                q = 0.95,
                N = 100) {
  d <- dim(Y)
  n <- d[1]
  p <- d[2]
  nullest <-
    null_model(Y, as.matrix(cov), lambda2 = lambda2, trace.it = F)
  X0 <- nullest$X
  lambdas <- rep(0, N)
  m <- sum(!is.na(Y))
  for (i in 1:N)
  {
    Ysimul <- matrix(stats::rpois(n * p, exp(c(X0))), nrow = n)
    nullest2 <- null_model(Ysimul, cov, lambda2 = lambda2)
    X0simul <- nullest2$X
    dat <-
      sweep(Ysimul - exp(X0simul), 2, colMeans(Ysimul - exp(X0simul)))
    dat <-
      sweep(Ysimul - exp(X0simul), 1, rowMeans(Ysimul - exp(X0simul)))
    lambdas[i] <-
      (1 / m) * svd::propack.svd(dat, neig = 1, opts = list(maxiter = 1e5))$d
    cat('\r', i, "/", N)
  }
  return(stats::quantile(lambdas, q)[[1]])
}


#' cv.lori
#'
#' @param Y [matrix, data.frame] abundance table (nxp)
#' @param cov [matrix, data.frame] design matris (npxq)
#' @param N [integer] number of cross-validation folds
#' @param thresh [positive number] convergence threshold, default is 1e-5
#' @param maxit [integer] maximum number of iterations, default is 100
#' @param rank.max [integer] maximum rank of interaction matrix, default is 2
#' @param trace.it [boolean] whether information about convergence should be printed
#' @param parallel [boolean] whether the N-fold cross-validation should be parallelized, default value is TRUE
#' @param len [integer] the size of the grid
#'
#' @return A list with the following elements
#' \item{lambda1}{regularization parameter estimated by cross-validation for nuclear norm penalty (interaction matrix)}
#' \item{lambda2}{regularization parameter estimated by cross-validation for l1 norm penalty (main effects)}
#' \item{errors}{a table containing the prediction errors for all pairs of parameters}

#' @export
#' @import data.table doParallel parallel corpcor foreach
#'
#' @examples
#' X <- matrix(rnorm(50), 10)
#' Y <- matrix(rpois(10, 1:10), 5)
#' res <- cv.lori(Y, X, N=2, len=2)
cv.lori <- function(Y,
                    cov = NULL,
                    N = 10,
                    thresh = 1e-5,
                    maxit = 100,
                    rank.max = 5,
                    trace.it = F,
                    parallel = F,
                    len = 20) {
  Y <- as.matrix(Y)
  Y2 <- Y
  Y2[is.na(Y2)] <- 0
  m <- sum(!is.na(Y))
  d <- dim(Y)
  n <- d[1]
  p <- d[2]
  lambda2.max <-
    max(c(rowSums(Y, na.rm = T) / m, colSums(Y, na.rm = T) / m))
  lambda1.max <- max(svd(Y2 / m)$d)
  lambda2.min <- 1e-4 * lambda2.max
  lambda1.min <- svd(Y2 / m)$d[min(p,rank.max+1)]
  grid.lambda1 <-
    exp(seq(log(lambda1.min), log(lambda1.max), length.out = len))
  grid.lambda2 <-
    exp(seq(log(lambda2.min), log(lambda2.max), length.out = len))
  grid <- as.matrix(data.table::CJ(grid.lambda1, grid.lambda2))
  grid <- grid[nrow(grid):1,]
  na_func <- function(x, prob = 0.1) {
    x <- as.matrix(x)
    omega <- !is.na(x)
    obs.idx <- which(omega)
    yp <- x
    yp[sample(obs.idx, round(prob * sum(omega)))] <- NA
    return(yp)
  }
  if (parallel) {
    nbco <- detectCores()
    cl <- makeCluster(nbco)
    registerDoParallel(cl)
    res.cv <-
      foreach(k = 1:N,
              .packages = c("lori", "corpcor", "parallel", "glmnet")) %dopar% {
                sapply(1:nrow(grid),
                       function(i) {
                         yy <- na_func(as.matrix(Y), prob = 0.1)
                         if (trace.it)
                           print(paste("lambda", i))
                         res <-
                           lori(as.matrix(yy),
                                cov,
                                lambda1 = grid[i, 1],
                                lambda2 = grid[i, 2])$imputed
                         return(sqrt(sum((res - Y) ^ 2, na.rm = T)))
                       })
              }
  } else{
    ylist <-
      lapply(1:N, function(k)
        na_func(as.matrix(Y), prob = 0.05))
    res.cv <-   lapply(1:N, function(k) {
      if (trace.it)
        print(paste("boot", k))
      sapply(1:nrow(grid),
             function(i) {
               if (trace.it)
                 print(paste("lambda", i))
               res <-
                 lori(ylist[[k]], cov, lambda1 = grid[i, 1], lambda2 = grid[i, 2])$imputed
               return(sqrt(sum((res - Y) ^ 2, na.rm = T)))
             })
    })
  }
  res.cv <- colMeans(do.call(rbind, res.cv))
  l <- which.min(res.cv)
  lambda1 <- grid[l, 1]
  lambda2 <- grid[l, 2]
  dat <-
    data.frame(errors = res.cv,
               lambda1 = grid[, 1],
               lambda2 = grid[, 2])
  return(list(
    lambda1 = lambda1,
    lambda2 = lambda2,
    errors = dat
  ))
}
