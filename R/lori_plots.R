#' plot_cov
#'
#' @param res_lori lori output
#' @importFrom grDevices grey
#' @importFrom graphics abline axis image lines par plot plot.new polygon
#' @importFrom stats cov lm lowess predict quantile
#' @import gridExtra
#' @import lattice
#' @import grid
#' @export
#' @examples
#' Y <- matrix(rpois(9, lambda=1:9),nrow=3)
#' Y[sample(1:9, 2)] <- NA
#' cov <- matrix(rnorm(18), nrow=9)
#' colnames(cov) <- c("cov1","cov2")
#' res <- lori(Y,cov,1,1,rank.max=2,maxit=1)
#' plot_cov(res)
plot_cov <- function(res_lori){
  plot.new()
  alpha <- res_lori$alpha
  if(sum(abs(alpha)>0)>0){
    sig <- min(abs(alpha)[abs(alpha)>0])
    dig <- TRUE
    t <- 1
    while(dig){
      dig <- round(sig,digits=1)>0
      t <- t+1
      sig <- 10*sig
      if(t>10) dig <- F
    }
    sgn <- sign(alpha)
    ord <- sort(abs(alpha),decreasing=T, index.return=T)
    alpha <- ord$x
    idx <- ord$ix
    nms <- names(alpha)
    alpha <- sgn[idx]*alpha
    names(alpha) <- nms
    grid.table(round(alpha,digits=t+1), rows=nms)
  } else grid.table(alpha, rows=rownames(alpha))
}

#' plot_counts
#'
#' @param res_lori lori output
#' @param r.cov a vector of indices indicating the indices of the row covariates
#' @param c.cov a vector of indices indicating the indices of the column covariates
#' @param rc.cov a vector of indices indicating the indices of the row-column covariates
#' @importFrom grDevices grey
#' @importFrom graphics abline axis image lines par plot plot.new polygon
#' @importFrom stats cov lm lowess predict quantile
#' @export
#'
#' @examples
#' Y <- matrix(rpois(9, lambda=1:9),nrow=3)
#' Y[sample(1:9, 2)] <- NA
#' cov <- matrix(rnorm(18), nrow=9)
#' colnames(cov) <- c("cov1","cov2")
#' res <- lori(Y,cov,1,1,rank.max=2,maxit=1)
#' plot_counts(res, r.cov=1:2)
plot_counts <- function(res_lori, r.cov =NULL, c.cov=NULL, rc.cov=NULL){
  cov <- as.matrix(res_lori$cov)
  theta <- res_lori$theta
  d <- dim(theta)
  n <- d[1]
  p <- d[2]
  if(!is.null(r.cov)){
    par(mai=c(1.02, 0.82, 0.82, 0.42))
    par(xpd = FALSE)
    alpha <- res_lori$alpha[r.cov]
    qq <- quantile(rowMeans(exp(matrix(cov[, r.cov]%*%alpha, nrow = n))), 0.99)
    colo <- sapply(rowMeans(exp(matrix(cov[, r.cov]%*%alpha, nrow = n))),
                   function(t) if(t<qq) "black" else "red")
    plot(1:n, rowMeans(exp(matrix(cov[, r.cov]%*%alpha, nrow = n))),
          ylab = "Row covariate effects (count scale)", xlab = "Row", col=colo)
    abline(h=qq, col="darkblue", lwd=0.5)
  }
  if(!is.null(c.cov)){
    par(mai=c(1.02, 0.82, 0.82, 0.42))
    beta <- res_lori$alpha[c.cov]
    plot(1:p, colMeans(exp(matrix(cov[, c.cov]%*%beta, nrow = n))), ylab = "Column effects (count scale)", xlab = "Column")
    lines(lowess(1:p, colMeans(exp(matrix(cov[, c.cov]%*%beta, nrow = n)))), col = 2, lwd = 1)
  }
  if(!is.null(rc.cov)){
    ab <- res_lori$alpha[rc.cov]
    par(mai=c(1,1,1,2))
    data <- exp(matrix(cov[,rc.cov]%*%ab,nrow=n))
    x <- (1:nrow(data))
    y <- (1:ncol(data))
    image(y, x, t(data), col = grey(seq(1, 0, length = 100)),
          axes=FALSE,xlab="",ylab="",srt=45)
    axis(3, at = 1:ncol(data), labels=colnames(data),srt=45,tick=T)
    axis(2, at = 1:nrow(data), labels=rownames(data),srt=45,tick=T)
    legend.col(grey(seq(0, 1, length = 256)), lev = seq(0,max(data),length.out = 3))
  }
  par(mai=c(1,1,1,2))
  data <- exp(theta)
  x <- (1:nrow(data))
  y <- (1:ncol(data))
  image(y, x, t(data), col = grey(seq(1, 0, length = 100)),
        axes=FALSE,xlab="",ylab="",srt=45)
  axis(3, at = 1:ncol(data), labels=colnames(data),srt=45,tick=T)
  axis(2, at = 1:nrow(data), labels=rownames(data),srt=45,tick=T)
  legend.col(grey(seq(0, 1, length = 256)), lev = seq(0,max(data),length.out = 3))
}


#' legend.col
#'
#' @param col color
#' @param lev vector of levels
#' @importFrom grDevices grey
#' @importFrom graphics abline axis image lines par plot plot.new polygon
#' @importFrom stats cov lm lowess predict quantile
#' @export
#' @examples
#' image(matrix(rnorm(18), nrow=9), col = grey(seq(1, 0, length = 100)))
#' legend.col(grey(seq(0, 1, length = 256)), lev = seq(0,100,length.out = 3))
legend.col <- function(col, lev){

  opar <- par

  n <- length(col)

  bx <- par("usr")

  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 100,
              bx[2] + (bx[2] - bx[1]) / 100 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n

  xx <- rep(box.cx, each = 2)

  par(xpd = TRUE)
  for(i in 1:n){

    yy <- c(box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[n-i+1], border = col[n-i+1])

  }
  par(new = TRUE)
  plot(0, 0, type = "n",
       ylim = c(lev[1], lev[length(lev)]),
       yaxt = "n", ylab = "",
       xaxt = "n", xlab = "",
       frame.plot = FALSE)
  axis(side = 4, las = 2, tick = FALSE, line = .2)
  par <- opar
}
