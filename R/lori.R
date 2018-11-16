#' Estimation algorithm for the LoRI parameters.
#'
#' @param Y A matrix of counts (nxp).
#' @param lambda1 A number, the regularization parameter for the interaction matrix.
#' @param lambda2 A number, the regularization parameter for the covariate effects.
#' @param rank.max An integer, maximum rank of interaction matrix (smaller than min(n-1,p-1))
#' @param cov A matrix of row and column covariates (np*K1) in order row1xcol1,row2xcol2,..,rownxcol1,row1xcol2,row2xcol2,...,...,rownxcolp
#' @param thresh A number, convergence tolerance of algorithm, by default \code{1e-6}.
#' @param maxit An integer, maximum allowed number of iterations.
#' @param mu_init A number, initial intercept value
#' @param alpha_init A vector of length K1, initial value of covariate effects
#' @param theta_init An nxp matrix, initial value of the interactions.
#' @param trace.it A boolean, whether convergence information should be printed
#' @import glmnet
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
           rank.max = 2,
           thresh = 1e-6,
           maxit = 1e3,
           mu_init = 0,
           alpha_init = NULL,
           theta_init = NULL,
           trace.it = T)
  {
    Y = as.matrix(Y)
    n = nrow(Y)
    p = ncol(Y)
    rank.max = min(n - 1, p - 1, rank.max)
    Omega = 1 - 1 * is.na(Y)
    m = sum(!is.na(Y))
    if (is.null(theta_init)) {
      theta = matrix(0, nrow = n, ncol = p)
    } else{
      theta = theta_init
    }
    mu <- mu_init
    if (is.null(alpha_init)) {
      if (!is.null(cov))
        alpha_init <- rep(0, ncol(cov))
      else
        alpha_init <- NULL
    }
    alpha <- alpha_init
    alpmat <- matrix(cov %*% alpha, nrow = n)
    error = 1
    objective = rep(0, maxit)
    objective[1] <- 1
    count = 1
    Y2 <- Y
    Y2[is.na(Y)] <- 0
    X <- mu + alpmat + theta
    Ytalp <- c(Y[!is.na(Y)])
    while (error > thresh && count < maxit) {
      X.tmp = X
      alpha.tmp <- alpha
      mu.tmp <- mu
      theta.tmp <- theta
      d.tmp <- svd(theta.tmp)$d
      alpha <-
        glmnet(
          cov[!is.na(Y),],
          Ytalp,
          family = "poisson",
          offset = c(theta[!is.na(Y)]),
          lambda = lambda2
        )
      mu <- alpha$a0
      alpha <- alpha$beta
      alpmat <- matrix(cov %*% alpha, nrow = n)
      flag <- TRUE
      grad_theta <- (-Y2 + exp(mu + alpmat + theta)) / m
      grad_theta <- sweep(grad_theta, 2, colMeans(grad_theta))
      grad_theta <- sweep(grad_theta, 1, rowMeans(grad_theta))
      step <- 1
      number <-
        sum(-Y * (mu + alpmat + theta) + exp(mu + alpmat + theta), na.rm = T) /
        m + lambda1 * sum(d.tmp) + lambda2 * sum(abs(alpha))
      while (flag) {
        svd_theta <-
          svd(theta - step * grad_theta, nu = rank.max, nv = rank.max)
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
          number - sum(-Y * (mu + alpmat + theta2) + exp(mu + alpmat + theta2), na.rm = T) /
          m - lambda1 * sum(d) - lambda2 * sum(abs(alpha))
        if (diff < -abs(number) * thresh)
          step <- 0.5 * step
        else
          flag <- FALSE
      }
      theta <- theta2
      X <- mu + alpmat + theta
      count = count + 1
      objective[count] = -sum(Y * X - Omega * exp(X), na.rm = T) / m + lambda1 * sum(d) + lambda2 *
        sum(abs(alpha))
      error = abs(objective[count] - objective[count - 1]) / abs(objective[count])
      alpha <- as.vector(alpha)
      names(alpha) <- colnames(cov)
      rownames(X) <- rownames(Y)
      colnames(X) <- colnames(Y)
      rownames(theta) <- rownames(Y)
      colnames(theta) <- colnames(Y)
    }
    rank = sum(svd(theta)$d > 5e-06)
    return(structure(
      list(
        X = X,
        theta = theta,
        mu = mu,
        alpha = alpha,
        objective = unlist(objective),
        iter = count,
        rank = rank,
        convergence = count < maxit
      )
    ))
  }


#' LORI method.
#'
#' @param Y A matrix of counts (nxp).
#' @param cov A matrix of row and column covariates (np*K1) in order row1xcol1,row2xcol1,..,rownxcol1,row1xcol2,row2xcol2,...,...,rownxcolp
#' @param lambda1 A positive number, the regularization parameter for the interaction matrix.
#' @param lambda2 A positive number, the regularization parameter for the covariate effects.
#' @param thresh A number, convergence tolerance of algorithm, by default \code{1e-6}.
#' @param maxit An integer, maximum allowed number of iterations.
#' @param trace.it A boolean, whether convergence information should be printed
#' @param rank.max An integer, maximum rank of interaction matrix (smaller than min(n-1,p-1))
#' @param plots A boolean indicating whether graphical displays should be returned
#' @param r.cov a vector of indices indicating the indices of the row covariates
#' @param c.cov a vector of indices indicating the indices of the column covariates
#' @param rc.cov a vector of indices indicating the indices of the row-column covariates
#' @param spe_time a boolean indicating whether columns of Y are temporal (months, years, etc.)
#' @import glmnet
#' @export
#' @return A list with the following elements
#' \item{X}{nxp matrix of log of expected counts}
#' \item{mu}{intercept}
#' \item{alpha}{coefficients of covariates}
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
           thresh = 1e-5,
           maxit = 1e3,
           trace.it = F,
           rank.max = 2,
           plots = F,
           r.cov = NULL,
           c.cov = NULL,
           rc.cov = NULL,
           spe_time = T)
  {
    Y <- as.matrix(Y)
    d <- dim(Y)
    n <- d[1]
    p <- d[2]
    rank.max = min(n - 1, p - 1, rank.max)
    m <- sum(!is.na(Y))
    ymis <- c(Y)
    if (is.null(lambda2)) {
      print("computing regularization parameter 1")
      lambda2 <-
        cv.glmnet(cov[!is.na(ymis),], ymis[!is.na(ymis)], family = "poisson")$lambda.min
      lambda1tmp <- qut(Y, cov, lambda2)
      th <-
        lori(
          Y,
          cov,
          lambda1tmp,
          1e5,
          trace.it = F,
          plots = F,
          maxit = maxit,
          rank.max = rank.max
        )$theta
      lambda2 <-
        cv.glmnet(cov[!is.na(ymis),], ymis[!is.na(ymis)], family = "poisson",
                  offset = c(th)[!is.na(ymis)])$lambda.min
    }
    if (is.null(lambda1)) {
      print("computing regularization parameter 2")
      lambda1 <- qut(Y, cov, lambda2)
    }
    mu <- 0
    X = matrix(0, n, p)
    theta = matrix(0, nrow = n, ncol = p)
    if (!is.null(cov))
      alpha <- rep(0, ncol(cov))
    else
      alpha <- NULL
    alpmat <- matrix(cov %*% alpha, nrow = n)
    iter <- 1
    if (trace.it)
      print("fitting model...")
    alpha_init <- if (is.null(cov))
      NULL
    else
      as.numeric(alpha)
    res <- altmin(
      Y,
      lambda1,
      lambda2,
      cov = cov,
      mu_init = mu,
      alpha_init = alpha_init,
      theta_init = theta,
      thresh = thresh,
      trace.it = trace.it,
      rank.max = rank.max
    )
    X <- res$X
    mu <- res$mu
    alpha <- res$alpha
    theta <- res$theta
    imputed <- exp(X)
    imputed[!is.na(Y)] <- Y[!is.na(Y)]
    means <- exp(X)
    res <-
      structure(
        list(
          X = X,
          mu = mu,
          alpha = alpha,
          theta = theta,
          imputed = imputed,
          means = means,
          cov = cov
        )
      )
    if (plots) {
      opar <- par()
      plot_cov(res)
      plot_counts(res, r.cov, c.cov, rc.cov)
      if (spe_time) {
        plot(1:p,
             colSums(res$imputed),
             ylab = "Estimated total counts",
             xlab = "Time")
        mod <-
          lm(count ~ time, data = cbind.data.frame(count = colSums(res$imputed), time =
                                                     1:p))
        lines(predict(mod), col = 2, lwd = 1)
      }
      par <- opar
    }
    return(res)
  }


#' null_model
#'
#' @param Y nxp count matrix
#' @param lambda2 regularization parameter
#' @param cov npxK matrix of covariates
#' @param thresh convergence threshold
#' @param maxit maximum number of iterations
#' @param trace.it boolean indicating whether information about convergence should be displayed
#' @return A list with the following elements
#' \item{X}{nxp matrix of log of expected counts}
#' \item{mu}{intercept}
#' \item{alpha}{coefficients of covariates}
#'
#' @import glmnet
#' @export
#'
#' @examples
#' \dontshow{
#' X <- matrix(rnorm(50), 25)
#' Y <- matrix(rpois(25, 1:25), 5)
#' res <- null_model(Y, 10, X)
#' }
null_model <-
  function(Y,
           lambda2,
           cov = NULL,
           thresh = 1e-6,
           maxit = 1e3,
           trace.it = F)
  {
    Y = as.matrix(Y)
    n = nrow(Y)
    p = ncol(Y)
    Omega = 1 - 1 * is.na(Y)
    m = sum(!is.na(Y))
    X = matrix(0, n, p)
    mu <- 0
    if (!is.null(cov))
      alpha <- rep(0, ncol(cov))
    else
      alpha <- NULL
    alpmat <- matrix(cov %*% alpha, nrow = n)
    error = 1
    objective = rep(0, maxit)
    objective[1] <- 1
    count = 1
    Ytalp <- c(Y[!is.na(Y)])
    while (error > thresh && count < maxit) {
      if ((count == 1) &&
          (trace.it))
        cat("\r", "starting optimization...")
      X.tmp = X
      alpha.tmp <- alpha
      mu.tmp <- mu
      alpha <-
        glmnet(cov[!is.na(Y),], Ytalp, family = "poisson", lambda = lambda2)
      mu <- alpha$a0
      alpha <- alpha$beta
      alpmat <- matrix(cov %*% alpha, nrow = n)
      X <- mu + alpmat
      count = count + 1
      objective[count] = -sum(Y * X - Omega * exp(X), na.rm = T) / m + lambda2 *
        sum(abs(alpha))
      error = abs(objective[count] - objective[count - 1]) / abs(objective[count])
      if (trace.it)
        cat("\r", "objective:", objective[count])
    }
    if (trace.it) {
      if (count < maxit)
        cat("\r", "model converged !")
      if (count == maxit)
        cat("\r",
            "model did not converge, stopped after",
            maxit,
            "iterations.")
    }

    return(structure(list(
      X = X, mu = mu, alpha = alpha
    )))
  }
