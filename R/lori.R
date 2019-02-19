#' Estimation algorithm for the LORI parameters.
#'
#' @param Y [matrix, data.frame] count table (nxp).
#' @param lambda1 [positive number] the regularization parameter for the interaction matrix.
#' @param lambda2 [positive number] the regularization parameter for the covariate effects.
#' @param cov [matrix, data.frame] design matrix (np*q) in order row1xcol1,row2xcol2,..,rownxcol1,row1xcol2,row2xcol2,...,...,rownxcolp
#' @param rank.max [integer] maximum rank of interaction matrix (smaller than min(n-1,p-1))
#' @param thresh [positive number] convergence tolerance of algorithm, by default \code{1e-6}.
#' @param maxit [integer] maximum allowed number of iterations.
#' @param trace.it [boolean] whether convergence information should be printed
#' @param reff [boolean] whether row effects should be fitted, default value is TRUE
#' @param ceff [boolean] whether column effects should be fitted, default value is TRUE
#' @import glmnet corpcor
#'
#' @examples
#' \dontrun{
#' X = matrix(rnorm(50), 25)
#' Y <- matrix(rpois(25, 1:25), 5)
#' res <- altmin(Y, 100, 100, X)
#' }
altmin <-
  function(Y,
           lambda1,
           lambda2,
           cov = NULL,
           rank.max = 5,
           thresh = 1e-6,
           maxit = 1e3,
           trace.it = T,
           reff = T,
           ceff = T)
  {
    Y = as.matrix(Y)
    n = nrow(Y)
    p = ncol(Y)
    rank.max = min(n - 1, p - 1, rank.max)
    Omega = 1 - 1 * is.na(Y)
    m = sum(!is.na(Y))
    theta = matrix(0, nrow = n, ncol = p)
    alpha <- rep(0, n)
    beta <- rep(0, p)
    if (!is.null(cov))
      epsilon <- rep(0, ncol(cov))
    else
      epsilon <- 0
    alpmat <- matrix(rep(alpha, p), nrow = n)
    betmat <- matrix(rep(beta, each = n), nrow = n)
    epsmat <- matrix(cov %*% epsilon, nrow = n)
    error = 1
    objective = rep(0, maxit)
    objective[1] <- 1
    count = 1
    Y2 <- Y
    Y2[is.na(Y)] <- 0
    X <- alpmat + betmat + epsmat + theta
    Ytalp <- c(Y[!is.na(Y)])
    while (error > thresh && count < maxit) {
      d.tmp <- fast.svd(theta)$d
      alpha <- rep(0, n)
      if (reff) {
        vec1 <-
          (rowSums(Y, na.rm = T) - lambda2) / rowSums(Omega * exp(betmat + epsmat +
                                                                    theta))
        alpha[vec1 > 1] <-
          log(((rowSums(Y, na.rm = T) - lambda2) / rowSums(Omega * exp(
            betmat + epsmat +
              theta
          )))[vec1 > 1])
        vec2 <-
          (rowSums(Y, na.rm = T) + lambda2) / rowSums(Omega * exp(betmat + epsmat +
                                                                    theta))
        alpha[vec2 < 1] <-
          log(((rowSums(Y, na.rm = T) + lambda2) / rowSums(Omega * exp(
            betmat + epsmat +
              theta
          )))[vec2 < 1])
        alpmat <- matrix(rep(alpha, p), nrow = n)
      }
      beta <- rep(0, p)
      if (ceff) {
        vec1 <-
          (colSums(Y, na.rm = T) - lambda2) / colSums(Omega * exp(alpmat + epsmat +
                                                                    theta))
        beta[vec1 > 1] <-
          log(((colSums(Y, na.rm = T) - lambda2) / colSums(Omega * exp(
            alpmat + epsmat + theta
          )))[vec1 > 1])
        vec2 <-
          (colSums(Y, na.rm = T) + lambda2) / colSums(Omega * exp(alpmat + epsmat +
                                                                    theta))
        beta[vec2 < 1] <-
          log(((colSums(Y, na.rm = T) + lambda2) / colSums(Omega * exp(
            alpmat + epsmat + theta
          )))[vec2 < 1])
        betmat <- matrix(rep(beta, each = n), nrow = n)
      }
      off <- alpmat + betmat + theta
      epsilon <-
        glmnet(
          cov[!is.na(Y),],
          Ytalp,
          family = "poisson",
          offset = c(off[!is.na(Y)]),
          lambda = lambda2,
          intercept = F
        )
      epsilon <- epsilon$beta
      epsmat <- matrix(cov %*% epsilon, nrow = n)
      flag <- TRUE
      grad_theta <-
        (-Y2 + exp(alpmat + betmat + epsmat + theta)) / m
      grad_theta <- sweep(grad_theta, 2, colMeans(grad_theta))
      grad_theta <- sweep(grad_theta, 1, rowMeans(grad_theta))
      step <- 1
      armijo.step <- 1e-2
      ref.obj <-
        sum(-Y * (alpmat + betmat + epsmat + theta) + exp(alpmat + betmat +
                                                            epsmat + theta),
            na.rm = T) / m + lambda1 * sum(d.tmp) + lambda2 * sum(abs(epsilon))
      while (flag) {
        svd_theta <-
          fast.svd(theta - step * grad_theta)
        d <- svd_theta$d[1:rank.max]
        u <- svd_theta$u[, 1:rank.max]
        v <- svd_theta$v[, 1:rank.max]
        d <- d[d > (lambda1 * step)] - lambda1 * step
        if (length(d) == 0) {
          theta2 <- matrix(rep(0, n * p), nrow = n)
        } else if (length(d) == 1) {
          u <- svd_theta$u[, 1:length(d)]
          v <- svd_theta$v[, 1:length(d)]
          theta2 <- d * u %*% t(v)
        } else {
          u <- svd_theta$u[, 1:length(d)]
          v <- svd_theta$v[, 1:length(d)]
          theta2 <- u %*% diag(d) %*% t(v)
        }
        diff <-
          sum(-Y * (alpmat + betmat + epsmat + theta2) + exp(alpmat + betmat +
                                                               epsmat + theta2),
              na.rm = T) / m + lambda1 * sum(d) + lambda2 * sum(abs(epsilon)) - ref.obj
        if (((diff > step * armijo.step * sum(grad_theta * (theta2 - theta)))&&step>1e-5))
          step <- 0.5 * step
        else
          flag <- FALSE
      }
      theta <- theta2
      X <- alpmat + betmat + epsmat + theta
      count = count + 1
      objective[count] = -sum(Y * X - Omega * exp(X), na.rm = T) / m + lambda1 *
        sum(d) + lambda2 * sum(abs(epsilon))
      error = abs(objective[count] - objective[count - 1]) / abs(objective[count])
      epsilon <- as.vector(epsilon)
      names(epsilon) <- colnames(cov)
      rownames(X) <- rownames(Y)
      colnames(X) <- colnames(Y)
      rownames(theta) <- rownames(Y)
      colnames(theta) <- colnames(Y)
    }
    return(structure(
      list(
        X = X,
        alpha = alpha,
        beta = beta,
        epsilon = epsilon,
        theta = theta,
        objective = unlist(objective),
        iter = count,
        rank = rank,
        convergence = count < maxit
      )
    ))
  }


#' LORI method.
#'
#' @param Y [matrix, data.frame] count table (nxp).
#' @param cov [matrix, data.frame] design matrix (np*q) in order row1xcol1,row2xcol2,..,rownxcol1,row1xcol2,row2xcol2,...,...,rownxcolp
#' @param lambda1 [positive number] the regularization parameter for the interaction matrix.
#' @param lambda2 [positive number] the regularization parameter for the covariate effects.
#' @param reff [boolean] whether row effects should be fitted, default value is TRUE
#' @param ceff [boolean] whether column effects should be fitted, default value is TRUE
#' @param rank.max [integer] maximum rank of interaction matrix (smaller than min(n-1,p-1))
#' @param thresh [positive number] convergence tolerance of algorithm, by default \code{1e-6}.
#' @param maxit [integer] maximum allowed number of iterations.
#' @param trace.it [boolean] whether convergence information should be printed
#' @import glmnet
#' @export
#' @return A list with the following elements
#' \item{X}{nxp matrix of log of expected counts}
#' \item{alpha}{row effects}
#' \item{beta}{column effects}
#' \item{epsilon}{covariate effects}
#' \item{theta}{nxp matrix of row-column interactions}
#' \item{imputed}{nxp matrix of imputed counts}
#' \item{means}{nxp matrix of expected counts (exp(X))}
#' \item{cov}{npxK matrix of covariates}
#' @examples
#' \dontshow{
#' X <- matrix(rnorm(50), 25)
#' Y <- matrix(rpois(25, 1:25), 5)
#' res <- lori(Y, X, 10, 10)
#' }
lori <-
  function(Y,
           cov = NULL,
           lambda1 = NULL,
           lambda2 = NULL,
           reff = T,
           ceff = T,
           rank.max = 5,
           thresh = 1e-5,
           maxit = 1e3,
           trace.it = F)
  {
    Y <- as.matrix(Y)
    d <- dim(Y)
    n <- d[1]
    p <- d[2]
    if (!is.null(cov))
      q <- ncol(cov)
    else
      q <- 0
    rank.max = min(n - 1, p - 1, rank.max)
    res <- altmin(
      Y,
      lambda1,
      lambda2,
      cov = cov,
      thresh = thresh,
      trace.it = trace.it,
      rank.max = rank.max,
      reff = reff,
      ceff = ceff
    )
    X <- res$X
    alpha <- res$alpha
    beta <- res$beta
    epsilon <- res$epsilon
    theta <- res$theta
    imputed <- exp(X)
    imputed[!is.na(Y)] <- Y[!is.na(Y)]
    means <- exp(X)
    res <-
      structure(
        list(
          X = X,
          alpha = alpha,
          beta = beta,
          epsilon = epsilon,
          theta = theta,
          imputed = imputed,
          means = means,
          cov = cov
        )
      )
    return(res)
  }



#' Estimation algorithm for the LORI parameters without interactions.
#'
#' @param Y [matrix, data.frame] count table (nxp).
#' @param lambda2 [positive number] the regularization parameter for the covariate effects.
#' @param cov [matrix, data.frame] design matrix (np*q) in order row1xcol1,row2xcol2,..,rownxcol1,row1xcol2,row2xcol2,...,...,rownxcolp
#' @param thresh [positive number] convergence tolerance of algorithm, by default \code{1e-6}.
#' @param maxit [integer] maximum allowed number of iterations.
#' @param trace.it [boolean] whether convergence information should be printed
#' @param reff [boolean] whether row effects should be fitted, default value is TRUE
#' @param ceff [boolean] whether column effects should be fitted, default value is TRUE
#' @import glmnet
#' @examples
#' \dontrun{
#' X = matrix(rnorm(50), 25)
#' Y <- matrix(rpois(25, 1:25), 5)
#' res <- null_model(Y, 100, 100, X)
#' }
null_model <-
  function(Y,
           lambda2,
           cov = NULL,
           thresh = 1e-6,
           maxit = 1e3,
           trace.it = T,
           reff = T,
           ceff = T)
  {
    Y = as.matrix(Y)
    n = nrow(Y)
    p = ncol(Y)
    Omega = 1 - 1 * is.na(Y)
    m = sum(!is.na(Y))
    alpha <- rep(0, n)
    beta <- rep(0, p)
    if (!is.null(cov))
      epsilon <- rep(0, ncol(cov))
    else
      epsilon <- 0
    alpmat <- matrix(rep(alpha, p), nrow = n)
    betmat <- matrix(rep(beta, each = n), nrow = n)
    epsmat <- matrix(cov %*% epsilon, nrow = n)
    error = 1
    objective = rep(0, maxit)
    objective[1] <- 1
    count = 1
    Y2 <- Y
    Y2[is.na(Y)] <- 0
    X <- alpmat + betmat + epsmat
    Ytalp <- c(Y[!is.na(Y)])
    while (error > thresh && count < maxit) {
      alpha <- rep(0, n)
      if (reff) {
        vec1 <-
          (rowSums(Y, na.rm = T) - lambda2) / rowSums(Omega * exp(betmat + epsmat))
        alpha[vec1 > 1] <-
          log(((rowSums(Y, na.rm = T) - lambda2) / rowSums(Omega * exp(
            betmat + epsmat
          )))[vec1 > 1])
        vec2 <-
          (rowSums(Y, na.rm = T) + lambda2) / rowSums(Omega * exp(betmat + epsmat))
        alpha[vec2 < 1] <-
          log(((rowSums(Y, na.rm = T) + lambda2) / rowSums(Omega * exp(
            betmat + epsmat
          )))[vec2 < 1])
        alpmat <- matrix(rep(alpha, p), nrow = n)
      }
      beta <- rep(0, p)
      if (ceff) {
        vec1 <-
          (colSums(Y, na.rm = T) - lambda2) / colSums(Omega * exp(alpmat + epsmat))
        beta[vec1 > 1] <-
          log(((colSums(Y, na.rm = T) - lambda2) / colSums(Omega * exp(
            alpmat + epsmat
          )))[vec1 > 1])
        vec2 <-
          (colSums(Y, na.rm = T) + lambda2) / colSums(Omega * exp(alpmat + epsmat))
        beta[vec2 < 1] <-
          log(((colSums(Y, na.rm = T) + lambda2) / colSums(Omega * exp(
            alpmat + epsmat
          )))[vec2 < 1])
        betmat <- matrix(rep(beta, each = n), nrow = n)
      }
      off <- alpmat + betmat
      epsilon <-
        glmnet(
          cov[!is.na(Y),],
          Ytalp,
          family = "poisson",
          offset = c(off[!is.na(Y)]),
          lambda = lambda2,
          intercept = F
        )
      epsilon <- epsilon$beta
      epsmat <- matrix(cov %*% epsilon, nrow = n)
      X <- alpmat + betmat + epsmat
      count = count + 1
      objective[count] = -sum(Y * X - Omega * exp(X), na.rm = T) / m + lambda2 * sum(abs(epsilon))
      error = abs(objective[count] - objective[count - 1]) / abs(objective[count])
      epsilon <- as.vector(epsilon)
      names(epsilon) <- colnames(cov)
      rownames(X) <- rownames(Y)
      colnames(X) <- colnames(Y)
    }
    return(structure(
      list(
        X = X,
        alpha = alpha,
        beta = beta,
        epsilon = epsilon,
        objective = unlist(objective),
        iter = count,
        convergence = count < maxit
      )
    ))
  }

