#' @name expected_berg
#'
#' @title Expected Frequencies by the BerG Regression Model
#'
#' @description Provides the expected frequencies after fit of
#' the BerG regression model according to the definition given in
#' Kleiber and Zeileis (2016).
#'
#' @param y vector of the observed values of the response variable.
#' @param mu vector of the fitted values.
#' @param phi vector of the fitted values for dispersion index.
#'
#' @return Expected frequencies after fitting the BerG regression model
#' to a dataset.
#'
# #' @details For a count y with possible outcomes y = 0, 1, 2, \ldots, the
# #' expected frequencies were defined in Kleiber and Zeileis (2016).
#'
#' @references Bourguignon, M. & Medeiros, R. (2020). A simple and useful regression model for fitting count data
#'
#' @references Kleiber, C., & Zeileis, A. (2016). Visualizing count data regressions using rootograms. The American Statistician, 70, 296-303
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
#'
#' @examples
#'
#' # Sample size
#' n <- 100
#
#' # Covariates
#' X <- cbind(rep(1,n), runif(n, 0,1), runif(n,0,1))
#' Z <- cbind(rep(1,n), runif(n, 0,1), runif(n,0,1))
#'
#' # Parameters and relations
#' beta <- c(1, 1.2, 0.2)
#' gama <- c(2, 1.5, 1.2)
#'
#' mu <- exp(X%*%beta)
#' phi <- exp(Z%*%gama)
#'
#' y <- rberg(n, mu, phi)
#'
#' expect_berg(y, mu, phi)
#'
expect_berg <- function(y, mu, phi){

  x <- sort(unique(y))
  n <- length(y)
  s <- rep(0, length(x))

  if (length(phi) == 1) phi = as.matrix(rep(phi, n))

  for(i in 1:n){
    s = s + dberg(x, mu[i], phi[i])
  }

  return(s)
}

skewness <- function (x, na.rm = FALSE, type = 3)
{
  if (any(ina <- is.na(x))) {
    if (na.rm)
      x <- x[!ina]
    else return(NA)
  }
  if (!(type %in% (1:3)))
    stop("Invalid 'type' argument.")
  n <- length(x)
  x <- x - mean(x)
  y <- sqrt(n) * sum(x^3)/(sum(x^2)^(3/2))
  if (type == 2) {
    if (n < 3)
      stop("Need at least 3 complete observations.")
    y <- y * sqrt(n * (n - 1))/(n - 2)
  }
  else if (type == 3)
    y <- y * ((1 - 1/n))^(3/2)
  y
}

kurtosis <- function (x, na.rm = FALSE, type = 3)
{
  if (any(ina <- is.na(x))) {
    if (na.rm)
      x <- x[!ina]
    else return(NA)
  }
  if (!(type %in% (1:3)))
    stop("Invalid 'type' argument.")
  n <- length(x)
  x <- x - mean(x)
  r <- n * sum(x^4)/(sum(x^2)^2)
  y <- if (type == 1)
    r - 3
  else if (type == 2) {
    if (n < 4)
      stop("Need at least 4 complete observations.")
    ((n + 1) * (r - 3) + 6) * (n - 1)/((n - 2) * (n - 3))
  }
  else r * (1 - 1/n)^2 - 3
  y
}

#' @name rqr_berg
#'
#' @title Randomized quantile residuals in the BerG regression model
#'
#' @description Randomized quantile residuals resulting from the fit of
#' the BerG regression model.
#'
#' @param y vector; observed values of the response variable.
#' @param mu fitted values.
#' @param phi fitted values for dispersion index.
#'
#' @return Randomized quantile residuals resulting from the fit of
#' the BerG regression model.
#'
#' @export
#'
#' @references
#' Bourguignon, M. and Medeiros, R. (2019). A flexibe, simple and useful regression model for fitting count data.
#' Dunn, P. K., and Smyth, G. K. (1996). Randomized quantile residuals. Journal of Computational and Graphical Statistics, 5, 236-244.
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
#'
rqr_berg <- function(y, mu, phi){

  n <- length(y)
  a <- vector()
  b <- vector()
  u <- vector()

  for(i in 1:n){
    a[i] <- pberg(y[i] - 1, mu[i], phi[i])
    b[i] <- pberg(y[i], mu[i], phi[i])
    u[i] <- stats::runif(1, a[i], b[i])
  }
  return(stats::qnorm(u))
}

#' @name rootogram_berg
#'
#' @title Rootograms for the BerG Regression Model
#'
#' @description Provides a rootograms that graphically compare
#'  (square roots) of empirical frequencies with fitted frequencies.
#'  For details see Kleiber and Zeileis (2016).
#'
#' @param object object of class 'bergrm'.
#'
#' @return The rootogram for a fitted BerG regression.
#'
#' @references Bourguignon, M. & Medeiros, R. (2020). A simple and useful regression model for fitting count data.
#'
#' @references Kleiber, C., & Zeileis, A. (2016). Visualizing count data regressions using rootograms. The American Statistician, 70, 296-303
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
#'
#' @examples
#'
#' data(grazing)
#'
#' # BerG fit
#' fit <- glm.bg(birds ~ ., data = grazing)
#'
#' # Rootogram
#' root_berg(fit)
#'
root_berg <- function(object){

  y <- object$response
  mu.h <- object$fitted.values
  phi.h <- object$phi.h
  n <- object$n.obs

  ob <- sort(unique(y))
  obs <- table(y)
  esp <-expect_berg(y, mu.h, phi.h)

  # Rootogram
  op <- graphics::par("mar")
  graphics::par(mar = c(5, 4.5, 4, 2) + 0.1)
  x.axis <- graphics::barplot(sqrt(obs), col = "lightgray",
                              xlab = "y", ylab = expression(sqrt("Frequency")),
                              ylim = c(0, max(sqrt(obs), sqrt(esp)) + 0.5))
  graphics::points(x.axis, sqrt(esp), col = "red4", type = "b", pch = 16)
  par(mar = op)
}


#' @name envel_berg
#'
#' @title Envelope Graph of the BerG Regression Model
#'
#' @description Provides the normal envelope simulated
#' probability plot of Pearson's residuals, and of the randomized
#' quantile residuals resulting from the BerG regression model fit.
#'
#' @param object object of class 'bergrm'.
#' @param residual character; specifies which residual should be produced
#'     in the normal probability plot. The available arguments are "all" (default)
#'     for both randomized quantile and Pearson residuals; "quantile" for randomized
#'     quantile residuals; and "pearson" for Pearson residuals.
#' @param R number of replicates.
#' @param control a list of control arguments specified via \code{berg_control}.
#' @param ... arguments passed to \code{berg_control}.
#'
#' @references
#' Bourguignon, M. and Medeiros, R. (2019). A flexibe, simple and useful regression model for fitting count data.
#' Kleiber, C., & Zeileis, A. (2016). Visualizing count data regressions using rootograms. The American Statistician, 70, 296-303
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @return Returns the normal envelope simulated probability plot of
#' Pearson residuals and random quantile residuals, respectively.
#'
#' @export
envel_berg <- function(object, residual = c("all", "quantile", "pearson"),
                       R = 99, control = berg_control(...), ...)
{
  # Model specifications
  y <- object$response
  X <- object$X
  Z <- object$Z
  p <- NCOL(X)
  k <- NCOL(Z)
  n <- length(y)
  delta <- as.numeric(y == 0)

  link <- object$link
  link.phi <- object$link.phi

  # Link functions
  g1 = g(link)$fun
  g1.inv = g(link)$inv
  g1. = g(link)$deriv.
  g2 = g(link.phi)$fun
  g2.inv = g(link.phi)$inv
  g2. = g(link.phi)$deriv.

  # Estimatives
  mu = stats::fitted.values(object)
  phi = object$phi.hat

  ######################
  #     Residuals      #
  #####################

  # Randomized quantile residuals
  rq <- object$residuals

  # Pearson residuals
  rp <- object$pearson.residuals

  ###################
  #     Envelope    #
  ###################
  rs_1 <- matrix(0, n, R)
  rs_2 <- matrix(0, n, R)

  for(i in 1:R){

    # Simulated sample
    y.tilde <- rberg(n, mu, phi)

    # Estimatives
    est.tilde <- mle_berg(y.tilde, X, Z, link = link, link.phi = link.phi,
                          control = berg_control(start = c(stats::coef(object, what = "mean"),
                                                           stats::coef(object, what = "dispersion"))))$est

    mu.tilde  <-  g1.inv(X%*%est.tilde[1:p])
    phi.tilde <- g2.inv(Z%*%est.tilde[(p + 1):(p + k)])

    # Empirical residuals
    rs_1[,i] <- as.numeric((y.tilde - mu.tilde) /
                             sqrt(phi.tilde * mu.tilde))

    rs_2[,i] <- rqr_berg(y.tilde, mu.tilde, phi.tilde)
  }

  # Sort
  rs_1 = apply(rs_1, 2, sort)
  rs_2 = apply(rs_2, 2, sort)

  # Min and max
  mint = apply(rs_1, 1, min);  minq = apply(rs_2, 1, min)
  Maxt = apply(rs_1, 1, max);  Maxq = apply(rs_2, 1, max)

  # 0.5 and 99.5 quantiles for the envelope
  mmt = apply(rs_1, 1, stats::quantile, probs = 0.005);   mmq = apply(rs_2, 1, stats::quantile, probs = 0.005)
  MMt = apply(rs_1, 1, stats::quantile, probs = 0.995); MMq = apply(rs_2, 1, stats::quantile, probs = 0.995)

  # 2.5, 5 e 97.5 quantiles for the envelope
  mt = apply(rs_1, 1, stats::quantile, probs = 0.025);   mq = apply(rs_2, 1, stats::quantile, probs = 0.025)
  Mt = apply(rs_1, 1, stats::quantile, probs = 0.975);  Mq = apply(rs_2, 1, stats::quantile, probs = 0.975)

  # Median
  at = apply(rs_1, 1, stats::quantile, probs = 0.5);aq = apply(rs_2, 1, stats::quantile, probs = 0.5)

  # Theoretical quantiles for the normal distribution
  qq = stats::qqnorm(1:n, axes = FALSE, xlab = " ", ylab = " ",
                     type = "l", lty = 1, plot.it = FALSE)$x

  ######################################
  # Envelope                           #
  ######################################

  residual <- match.arg(residual)

  if (residual != "all"){
    ask <- FALSE
  }else{
    ask <- TRUE
  }

  # Pearson residual
  if (residual != "quantile") {
  graphics::par(ask = ask)
  stats::qqnorm(rp, main = " ", xlab = "Theoretical quantile",
                ylab = "Pearson residual", type = "n")
  graphics::polygon (c(qq, sort(qq, decreasing = T)),
                     c(mint, sort(Maxt, decreasing = T)), col = "lightgray", border=NA)
  graphics::polygon(c(qq, sort(qq, decreasing = T)),
                    c(mmt, sort(MMt, decreasing = T)), col = "gray", border = NA)
  graphics::polygon(c(qq, sort(qq, decreasing = T)),
                    c(mt, sort(Mt, decreasing = T)), col = "darkgray", border = NA)

  graphics::points(qq, sort(rp), pch = "+")
  graphics::points(qq, at, type="l", lty=2)
  graphics::box()
  graphics::par(ask = FALSE)
  }


  # Randomized quantile residual
  if (residual != "pearson" ) {
  graphics::par(ask = ask)
  stats::qqnorm(sort(rq), main = " ", xlab = "Theoretical quantile",
                ylab = "Randomized quantile residual", type = "n")
  graphics::polygon(c(qq, sort(qq, decreasing = T)),
                    c(minq, sort(Maxq, decreasing = T)), col = "lightgray", border = NA)
  graphics::polygon(c(qq, sort(qq, decreasing = T)),
                    c(mmq, sort(MMq, decreasing = T)),col = "gray", border = NA)
  graphics::polygon(c(qq, sort(qq, decreasing = T)),
                    c(mq, sort(Mq, decreasing = T)), col = "darkgray", border = NA)

  graphics::points(qq, sort(rq), pch = "+")
  graphics::points(qq, aq, type = "l", lty = 2)
  graphics::box()
  graphics::par(mfrow=c(1,1), ask = FALSE)
  }
}
