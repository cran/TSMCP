#' Two stage multiple change point detection in linear models.
#'
#' This function provides a two‐stage procedure for simultaneously detecting multiple
#' change‐points in linear models. In the cutting stage, the change‐point problem is
#' converted into a model selection problem so that a modern model selection method
#' can be applied. In the refining stage, the change‐points obtained in the cutting
#' stage are finalized via a refining method. The tuning parameter lambda is chosen
#' by BIC.
#'
#'
#' @param Y an response vector.
#' @param X the n-by-p design matrix.
#' @param method the method to be used by lasso, adaptive lasso, mcp or scad.
#' See \code{\link[plus]{plus}} in R packages \pkg{plus} for details.
#' @param c ceiling(c*sqrt(length(Y))) is the length of each segments in spliting stage.
#' @return tsmcplm returns an object of class "tsmcplm".
#' An object of class "tsmcplm" is a list containing the following components:
#' @return \item{change.points}{estimators of change points.}
#' @export tsmcplm
#' @import plus lars
#' @seealso plus lars
#' @references {Jin, B., Wu, Y., & Shi, X. (2016). Consistent two‐stage multiple
#' change‐point detection in linear models. \emph{Canadian Journal of Statistics},
#' 44(2), 161-179.}
#' @examples ## example 1 : mean shift model
#' ## true change point location:
#' ## 100, 130, 150, 230, 250, 400, 440, 650, 760, 780, 810
#' Y <- rnorm(1000, 0, 0.5) +
#'   c(rep(0,100), rep(4,30),rep(-1,20), rep(2,80), rep(-2,20),
#'     rep(3,150), rep(-1, 40), rep(1,210), rep(5,110), rep(2,20),
#'     rep(7,30), rep(3,190))
#' ts.plot(Y)
#'
#' ##estimate change points
#' tsmcplm(Y = Y, X = NULL, method = "adapt", c = 0.3)
#'
#'
#' ## example 2: linear model:
#' ## a periodic auto correlation series with period 122 and
#' ## order of auto correlation 1
#'
#' ###
#' ## true change point location:
#' ## 200, 350, 450, 550, 700, and 850
#'
#'
#' n=1000
#' y <- rnorm(n)
#' for (t in 2:n) {
#'   y[t] <- cos(t*pi/61) + 3*sin(t*pi/61) + 0.5*y[t-1] + (2*sin(t*pi/61) +
#'                                                           0.1 * y[t-1])*(200 < t)+ (2* cos(t*pi/61) - 4 *sin(t*pi/61) -
#'                                                                                       0.6* y[t-1] )*(350 < t) + (2* sin(t*pi/61) +
#'                                                                                                                    0.7* y[t-1] )*(450 < t) + (-3* sin(t*pi/61) -
#'                                                                                                                                                 0.3* y[t-1] )*(550 < t) + (-3* cos(t*pi/61) +
#'                                                                                                                                                                              5* sin(t*pi/61))* (700 < t) + (3* cos(t*pi/61) -
#'                                                                                                                                                                                                               5* sin(t*pi/61) - 0.4*y[t-1] )* (850 < t) + rnorm(1)
#'
#' }
#' ts.plot(y)
#'
#' x <- sapply(2:n, function(t){cbind(cos(t*pi/61), sin(t*pi/61), y[t-1])}, simplify = FALSE)
#' x <- do.call(rbind, x)
#' tsmcplm(Y = y[-1], X = x, method = "adapt", c = 2)
#'
#'



tsmcplm <- function(Y, X, method = c("lasso", "adapt", "mcp", "scad"), c) {

  # library("lars")
  n <- length(Y)
  m <- ceiling(c * sqrt(n))
  q <- floor(n/m)
  K_temp <- matrix(0, nrow = q, ncol = q, byrow = TRUE)
  p <- dim(X)[2]
  if (length(p) == 0) {
    p <- 0
  }

  ########### transform###################


  X_temp <- cbind(rep(1, n), X)

  Y_temp <- c(Y)
  for (i in 1:q) K_temp[i, 1:i] <- rep(1, i)




  x <- NULL
  y <- NULL
  x[[1]] <- as.matrix(X_temp[1:((n - (q - 1) * m)), ])
  y[[1]] <- Y_temp[1:((n - (q - 1) * m))]

  for (i in 2:q) {
    x[[i]] <- as.matrix(X_temp[(n - (q - i + 1) * m + 1):((n - (q -
                                                                  i) * m)), ])
    y[[i]] <- Y_temp[(n - (q - i + 1) * m + 1):((n - (q - i) * m))]
  }
  X_temp1 <- lapply(1:length(x), function(j, mat, list) kronecker(K_temp[j,
                                                                         , drop = FALSE], x[[j]]), mat = K_temp, list = x)
  Xn <- do.call("rbind", X_temp1)


  w1 <- stats::lm(Y_temp ~ Xn - 1)$coefficients
  w <- w1

  for (i in 1:q) {
    w[c(((i - 1) * (p + 1) + 1):(i * (p + 1)))] <- sum(w1[c(((i - 1) *
                                                               (p + 1) + 1):(i * (p + 1)))]^2)

  }

  #################### adlasso#########################
  if (method == "adapt") {
    X_temp <- scale(Xn, center = FALSE, scale = 1/w)  # x1 times the weights with group
    object <- lars(X_temp, Y_temp, type = "lasso", intercept = FALSE,
                   normalize = FALSE)
    # get min BIC
    bic <- log(n) * object$df + n * log(as.vector(object$RSS)/n)  # rss/n version
    step.bic2 <- which.min(bic)  # step with min BIC
    coeff <- coef.lars(object, s = step.bic2, mode = "step")




    adcp.coef.s <- sum(abs(coeff))

    adcp.coef.v.m <- abs(matrix(c(coeff), q, (p + 1), byrow = T))

    adcp.coef.m <- c(apply(adcp.coef.v.m, 1, max))

    adcp.cp <- which(adcp.coef.m != 0)

    if (length(adcp.cp) > 1) {
      for (i in 2:length(adcp.cp)) {
        if (adcp.cp[i] - adcp.cp[i - 1] == 1)
          adcp.cp[i] <- 0
      }
    }


    adcp.cp1 <- adcp.cp[adcp.cp > 1 & adcp.cp < q]

    d1 <- length(adcp.cp1)

    if (d1 == 0) {
      adcpcss.cp <- 0

    }

    if (d1 >= 1) {
      # step 4 find change points
      adcpcss.cp <- NULL


      adcp.cp1 <- c(0, adcp.cp1, q + 1)
      for (i in 1:d1) {


        y1 <- NULL
        x1 <- NULL
        for (k in (adcp.cp1[i + 1] - 1):(adcp.cp1[i + 1] + 1)) # for(k in (adcp.cp1[i]+1):(adcp.cp1[i+2]-1))
        {
          y1 <- c(y1, y[[k]])
          x1 <- rbind(x1, as.matrix(x[[k]]))
        }


        # using cusum find change point.
        cp <- css(y1, x1)
        if (cp == 0)
          next
        if (cp != 0) {
          if (adcp.cp1[i + 1] == 0)
            adcpcss.cp <- c(adcpcss.cp, cp)
          if (adcp.cp1[i + 1] > 0)
            adcpcss.cp <- c(adcpcss.cp, cp + n - (q - adcp.cp1[i +
                                                                 1] + 2) * m)
          # mcpcss.cp=c(mcpcss.cp,cp+n-(q-mcp.cp1[i])*m)
        }


      }

      if (length(adcpcss.cp) == 0)
        adcpcss.cp <- 0 else adcpcss.cp <- adcpcss.cp

    }

    tt <- which(abs(diff(adcpcss.cp)) < 10)
    if (length(tt) > 0)
      adcpcss.cp <- adcpcss.cp[-tt]

    if (length(adcpcss.cp) >= 1 & min(adcpcss.cp) > 0) {

      cp1 <- c(0, adcpcss.cp, n)
      q <- length(cp1) - 1
      K_temp <- matrix(0, nrow = q, ncol = q, byrow = TRUE)

      for (i in 1:q) K_temp[i, 1:i] <- rep(1, i)

      X_temp <- cbind(rep(1, n), X)
      Y_temp <- c(Y)
      x <- NULL
      y <- NULL
      for (t in 2:length(cp1)) {
        x[[t - 1]] <- as.matrix(X_temp[(cp1[t - 1] + 1):cp1[t],
                                       ])
        y[[t - 1]] <- Y_temp[(cp1[t - 1] + 1):cp1[t]]
      }



      cp0 <- NULL
      for (t in 2:(length(cp1) - 1)) {
        x1 <- rbind(x[[t - 1]], x[[t]])
        y1 <- c(y[[t - 1]], y[[t]])
        tt <- css(y1, x1)
        if (tt != 0)
          cp0 <- c(cp0, tt + cp1[t - 1])
      }



      adcpcss.cp <- cp0
      tt <- which(abs(diff(adcpcss.cp)) < 10)
      if (length(tt) > 0)
        adcpcss.cp <- adcpcss.cp[-tt]
    }
  }
  ################################# lasso#######################

  if (method == "lasso") {
    X_temp <- Xn  # x1 times the weights with group
    object <- lars(X_temp, Y_temp, type = "lasso", intercept = FALSE,
                   normalize = FALSE)
    # get min BIC
    bic <- log(n) * object$df + n * log(as.vector(object$RSS)/n)  # rss/n version
    step.bic2 <- which.min(bic)  # step with min BIC
    coeff <- coef.lars(object, s = step.bic2, mode = "step")




    adcp.coef.s <- sum(abs(coeff))


    adcp.coef.v.m <- abs(matrix(c(coeff), q, (p + 1), byrow = T))
    adcp.coef.m <- c(apply(adcp.coef.v.m, 1, max))




    adcp.cp <- which(adcp.coef.m != 0)

    if (length(adcp.cp) > 1) {
      for (i in 2:length(adcp.cp)) {
        if (adcp.cp[i] - adcp.cp[i - 1] == 1)
          adcp.cp[i] <- 0
      }
    }


    adcp.cp1 <- adcp.cp[adcp.cp > 1 & adcp.cp < q]

    d1 <- length(adcp.cp1)

    if (d1 == 0) {
      adcpcss.cp <- 0

    }

    if (d1 >= 1) {
      # step 4 find change points
      adcpcss.cp <- NULL


      adcp.cp1 <- c(0, adcp.cp1, q + 1)
      for (i in 1:d1) {


        y1 <- NULL
        x1 <- NULL
        for (k in (adcp.cp1[i + 1] - 1):(adcp.cp1[i + 1] + 1)) # for(k in (adcp.cp1[i]+1):(adcp.cp1[i+2]-1))
        {
          y1 <- c(y1, y[[k]])
          x1 <- rbind(x1, as.matrix(x[[k]]))
        }


        # using cusum find change point.
        cp <- css(y1, x1)
        if (cp == 0)
          next
        if (cp != 0) {
          if (adcp.cp1[i + 1] == 0)
            adcpcss.cp <- c(adcpcss.cp, cp)
          if (adcp.cp1[i + 1] > 0)
            adcpcss.cp <- c(adcpcss.cp, cp + n - (q - adcp.cp1[i +
                                                                 1] + 2) * m)
          # mcpcss.cp=c(mcpcss.cp,cp+n-(q-mcp.cp1[i])*m)
        }


      }
      if (length(adcpcss.cp) == 0)
        adcpcss.cp <- 0 else adcpcss.cp <- adcpcss.cp
    }



    tt <- which(abs(diff(adcpcss.cp)) < 10)
    if (length(tt) > 0)
      adcpcss.cp <- adcpcss.cp[-tt]

    if (length(adcpcss.cp) >= 1 & min(adcpcss.cp) > 0) {
      cp1 <- c(0, adcpcss.cp, n)
      q <- length(cp1) - 1
      K_temp <- matrix(0, nrow = q, ncol = q, byrow = TRUE)

      for (i in 1:q) K_temp[i, 1:i] <- rep(1, i)

      X_temp <- cbind(rep(1, n), X)
      Y_temp <- c(Y)
      x <- NULL
      y <- NULL
      for (t in 2:length(cp1)) {
        x[[t - 1]] <- as.matrix(X_temp[(cp1[t - 1] + 1):cp1[t],
                                       ])
        y[[t - 1]] <- Y_temp[(cp1[t - 1] + 1):cp1[t]]
      }



      cp0 <- NULL
      for (t in 2:(length(cp1) - 1)) {
        x1 <- rbind(x[[t - 1]], x[[t]])
        y1 <- c(y[[t - 1]], y[[t]])
        tt <- css(y1, x1)
        if (tt != 0)
          cp0 <- c(cp0, tt + cp1[t - 1])
      }



      adcpcss.cp <- cp0

      tt <- which(abs(diff(adcpcss.cp)) < 10)
      if (length(tt) > 0)
        adcpcss.cp <- adcpcss.cp[-tt]


    }
  }
  ############################### mcp
  if (method == "mcp") {
    object <- plus(Xn, Y_temp, method = "mc+", gamma = 2.4, intercept = F,
                   normalize = F, eps = 1e-30)
    # step 1 estimate coef using BIC.
    bic <- log(dim(Xn)[1]) * object$dim + dim(Xn)[1] * log(as.vector((1 -
                                                                        object$r.square) * sum(Y_temp^2))/length(Y_temp))
    step.bic <- which.min(bic)
    mcp.coef <- coef.plus(object, lam = object$lam[step.bic])

    mcp.coef.s <- sum(abs(mcp.coef))


    mcp.coef.v.m <- abs(matrix(c(mcp.coef), q, (p + 1), byrow = T))
    mcp.coef.m <- c(apply(mcp.coef.v.m, 1, max))




    mcp.cp <- which(mcp.coef.m != 0)

    if (length(mcp.cp) > 1) {
      for (i in 2:length(mcp.cp)) {
        if (mcp.cp[i] - mcp.cp[i - 1] == 1)
          mcp.cp[i] <- 0
      }
    }


    mcp.cp1 <- mcp.cp[mcp.cp > 1 & mcp.cp < q]

    d1 <- length(mcp.cp1)

    if (d1 == 0) {
      mcpcss.cp <- 0
      adcpcss.cp <- mcpcss.cp
    }

    if (d1 >= 1) {
      # step 4 find change points
      mcpcss.cp <- NULL

      mcp.cp1 <- c(0, mcp.cp1, q + 1)
      for (i in 1:d1) {
        y1 <- NULL
        x1 <- NULL
        for (k in (mcp.cp1[i + 1] - 1):(mcp.cp1[i + 1] + 1)) # for(k in (mcp.cp1[i]+1):(mcp.cp1[i+2]-1))
        {
          y1 <- c(y1, y[[k]])
          x1 <- rbind(x1, as.matrix(x[[k]]))
        }


        # using cusum find change point.

        cp <- css(y1, x1)
        if (cp == 0)
          next
        if (cp != 0) {
          if (mcp.cp1[i + 1] == 0)
            mcpcss.cp <- c(mcpcss.cp, cp)
          if (mcp.cp1[i + 1] > 0)
            mcpcss.cp <- c(mcpcss.cp, cp + n - (q - mcp.cp1[i +
                                                              1] + 2) * m)
          # mcpcss.cp=c(mcpcss.cp,cp+n-(q-mcp.cp1[i])*m)

        }

      }
      if (length(mcpcss.cp) == 0)
        adcpcss.cp <- 0 else adcpcss.cp <- mcpcss.cp
    }


    tt <- which(abs(diff(adcpcss.cp)) < 10)
    if (length(tt) > 0)
      adcpcss.cp <- adcpcss.cp[-tt]

    if (length(adcpcss.cp) >= 1 & min(adcpcss.cp) > 0) {
      n <- length(Y)
      cp1 <- c(0, adcpcss.cp, n)
      q <- length(cp1) - 1
      K_temp <- matrix(0, nrow = q, ncol = q, byrow = TRUE)

      for (i in 1:q) K_temp[i, 1:i] <- rep(1, i)

      X_temp <- cbind(rep(1, n), X)
      Y_temp <- c(Y)
      x <- NULL
      y <- NULL
      for (t in 2:length(cp1)) {
        x[[t - 1]] <- as.matrix(X_temp[(cp1[t - 1] + 1):cp1[t],
                                       ])
        y[[t - 1]] <- Y_temp[(cp1[t - 1] + 1):cp1[t]]
      }


      cp0 <- NULL
      for (t in 2:(length(cp1) - 1)) {
        x1 <- rbind(x[[t - 1]], x[[t]])
        y1 <- c(y[[t - 1]], y[[t]])
        tt <- css(y1, x1)
        if (tt != 0)
          cp0 <- c(cp0, tt + cp1[t - 1])
      }


      adcpcss.cp <- cp0
      tt <- which(abs(diff(adcpcss.cp)) < 10)
      if (length(tt) > 0)
        adcpcss.cp <- adcpcss.cp[-tt]

    }


  }

  ############################### scad

  if (method == "scad") {
    object <- plus(Xn, Y_temp, method = "scad", gamma = 2.4, intercept = F,
                   normalize = F, eps = 1e-30)
    # step 1 estimate coef using BIC.
    bic <- log(dim(Xn)[1]) * object$dim + dim(Xn)[1] * log(as.vector((1 -
                                                                        object$r.square) * sum(Y_temp^2))/length(Y_temp))
    step.bic <- which.min(bic)
    mcp.coef <- coef.plus(object, lam = object$lam[step.bic])
    mcp.coef.s <- sum(abs(mcp.coef))
    mcp.coef.v.m <- abs(matrix(c(mcp.coef), q, (p + 1), byrow = T))
    mcp.coef.m <- c(apply(mcp.coef.v.m, 1, max))
    mcp.cp <- which(mcp.coef.m != 0)

    if (length(mcp.cp) > 1) {
      for (i in 2:length(mcp.cp)) {
        if (mcp.cp[i] - mcp.cp[i - 1] == 1)
          mcp.cp[i] <- 0
      }
    }


    mcp.cp1 <- mcp.cp[mcp.cp > 1 & mcp.cp < q]

    d1 <- length(mcp.cp1)

    if (d1 == 0) {
      mcpcss.cp <- 0
      adcpcss.cp <- mcpcss.cp
    }

    if (d1 >= 1) {
      # step 4 find change points
      mcpcss.cp <- NULL

      mcp.cp1 <- c(0, mcp.cp1, q + 1)
      for (i in 1:d1) {
        y1 <- NULL
        x1 <- NULL
        for (k in (mcp.cp1[i + 1] - 1):(mcp.cp1[i + 1] + 1)) # for(k in (mcp.cp1[i]+1):(mcp.cp1[i+2]-1))
        {
          y1 <- c(y1, y[[k]])
          x1 <- rbind(x1, as.matrix(x[[k]]))
        }


        # using cusum find change point.

        cp <- css(y1, x1)
        if (cp == 0)
          next
        if (cp != 0) {
          if (mcp.cp1[i + 1] == 0)
            mcpcss.cp <- c(mcpcss.cp, cp)
          if (mcp.cp1[i + 1] > 0)
            mcpcss.cp <- c(mcpcss.cp, cp + n - (q - mcp.cp1[i +
                                                              1] + 2) * m)
          # mcpcss.cp=c(mcpcss.cp,cp+n-(q-mcp.cp1[i])*m)

        }

      }
      if (length(mcpcss.cp) == 0)
        adcpcss.cp <- 0 else adcpcss.cp <- mcpcss.cp
    }


    tt <- which(abs(diff(adcpcss.cp)) < 10)
    if (length(tt) > 0)
      adcpcss.cp <- adcpcss.cp[-tt]

    if (length(adcpcss.cp) >= 1 & min(adcpcss.cp) > 0) {
      n <- length(Y)
      cp1 <- c(0, adcpcss.cp, n)
      q <- length(cp1) - 1
      K_temp <- matrix(0, nrow = q, ncol = q, byrow = TRUE)

      for (i in 1:q) K_temp[i, 1:i] <- rep(1, i)

      X_temp <- cbind(rep(1, n), X)
      Y_temp <- c(Y)
      x <- NULL
      y <- NULL
      for (t in 2:length(cp1)) {
        x[[t - 1]] <- as.matrix(X_temp[(cp1[t - 1] + 1):cp1[t],
                                       ])
        y[[t - 1]] <- Y_temp[(cp1[t - 1] + 1):cp1[t]]
      }

      cp0 <- NULL
      for (t in 2:(length(cp1) - 1)) {
        x1 <- rbind(x[[t - 1]], x[[t]])
        y1 <- c(y[[t - 1]], y[[t]])
        tt <- css(y1, x1)
        if (tt != 0)
          cp0 <- c(cp0, tt + cp1[t - 1])
      }



      adcpcss.cp <- cp0
      tt <- which(abs(diff(adcpcss.cp)) < 10)
      if (length(tt) > 0)
        adcpcss.cp <- adcpcss.cp[-tt]

    }


  }

  if (length(adcpcss.cp) == 0)
    return(0) else return(change.points = adcpcss.cp)

}





