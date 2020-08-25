#' @name bergrm-methods
#' @title Methods for 'bergrm' objects
#' @param object an object of class \code{bergrm}.
#'
#' @return .
#'
#' @references Bourguignon, M. and Medeiros, R. M. R. (2020). A simple
#'     and useful regression model for fitting count data.
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
NULL

bergrm <- function(x) UseMethod("bergrm")

# Print
#' @export
print.bergrm <- function(object)
{
  cat("Call:\n")
  print(object$call)
  cat("\nmean Coefficients:\n")
  print((object$coefficients)$mean)
  cat("\ndispersion Coefficients:\n")
  print((object$coefficients)$dispersion)
}

# Summary
#' @export
summary.bergrm <- function(object)
{
  n <- object$n.obs
  p <- object$p
  k <- object$k

  # Summary for residuals
  rq <- object$residuals
  TAB.residuals <- round(cbind(mean(rq), stats::sd(rq), skewness(rq), kurtosis(rq) + 3), 6)
  colnames(TAB.residuals) <- c("Mean", "Sd", "Skewness", "Kurtosis")
  rownames(TAB.residuals) <- " "

  # Summary for mu
  est.beta <- stats::coef(object, what = "mean")
  se.beta <- sqrt(diag(stats::vcov(object))[1:p])
  zval.beta <- est.beta/se.beta
  pval.beta <- 2*stats::pnorm(abs(zval.beta), lower.tail = FALSE)

  TAB.mu <- cbind(Estimate = est.beta,
                  `Std. Error` = se.beta,
                  `z value` = zval.beta,
                  `Pr(>|z|)` = pval.beta)
  rownames(TAB.mu) <- object$names.mean

  # Summary for phi
  est.gamma <- stats::coef(object, what = "dispersion")
  se.gamma <- sqrt(diag(stats::vcov(object))[(p + 1):(p + k)])
  zval.gamma <- est.gamma/se.gamma
  pval.gamma <- 2*stats::pnorm(abs(zval.gamma), lower.tail = FALSE)

  TAB.phi <- cbind(Estimate = est.gamma,
                   `Std. Error` = se.gamma,
                   `z value` = zval.gamma,
                   `Pr(>|z|)` = pval.gamma)
  rownames(TAB.phi) <- object$names.dispersion

  Log.lik <- stats::logLik(object)#; names(Log.lik) <- " "
  AIC <- as.numeric(object$AIC); names(AIC) <- " "
  BIC <- as.numeric(object$BIC); names(BIC) <- " "

  if (!is.null(object$test)){
    p.val <- pchisq(object$test, k-1, lower.tail = FALSE)
    TAB.test <- round(rbind(object$test,p.val), 4)
    colnames(TAB.test) <- c("S","W","LR","G")
    rownames(TAB.test) <- c("Value", "P value")
  }else{
    TAB.test <- NULL
  }

  out <- list(call = object$call, residuals = TAB.residuals,link = object$link,
              link.phi = object$link.phi, mean = TAB.mu, dispersion = TAB.phi, test = TAB.test,
              Log.lik = Log.lik, AIC = AIC, BIC = BIC)

  class(out) <- "summary.bergrm"
  out
}

# Print summary
#' @export
print.summary.bergrm <- function(object)
{
  cat("Call:\n")
  print(object$call)
  cat("\nSummary for residuals:\n")
  print(object$residuals)
  cat("\n----------------------------------------------------------------\n")
  cat("Mean:\n")
  cat("Coefficients:\n")
  stats::printCoefmat(object$mean)
  cat("\n----------------------------------------------------------------\n")
  cat("Dispersion:\n")
  cat("\nLink function:",object$link,"\n")
  cat("Coefficients:\n")
  stats::printCoefmat(object$dispersion)
  cat("\n----------------------------------------------------------------")
  if (!is.null(object$test)){
    cat("\n\nTest for constant dispersion:\n")
    print(object$test)
  }
  cat("\nIn addition, Log-lik value:",object$Log.lik,
      "\nAIC:",object$AIC,"and BIC:",object$BIC)
}


# Plot
#' @export
#' @rdname  bergrm-methods
plot.bergrm <- function(object)
{
  y <- object$response
  rq <- object$residuals
  mu.h <- object$fitted.values
  phi.h <- object$phi.h
  n <- object$n.obs

  graphics::par(mfrow=c(2,2))
  graphics::plot(mu.h, rq, xlab = "Fitted values", ylab = "Residuals", pch = "+")
  graphics::plot(1:n, rq, xlab = "Index", ylab = "Residuals", pch = "+")
  graphics::plot(stats::density(rq), xlab = "Residuals", ylab = "Density", main = " ", ylim = c(0, stats::dnorm(0)))
  graphics::curve(stats::dnorm(x), lty = 2, col = 2, add = T)
  stats::qqnorm(rq, xlab = "Theoretical quantile", ylab = "Residuals", pch = "+", main = " ")
  graphics::abline(0, 1, lty = 2)
  graphics::par(mfrow=c(1,1))

  ob <- sort(unique(y))
  obs <- table(y)
  esp <-expect_berg(y, mu.h, phi.h)

  # Observed vs expected
  graphics::barplot(t(cbind(sqrt(obs), sqrt(esp))), beside = T, col = c("gray40", "lightgray"), ylab = "Frequency", xlab = "y")
  graphics::legend("topright", c("Observed", "Expected"), col = c("gray40", "lightgray"),bty = "n", pch = 15)

  # Rootogram
  op <- par()
  par(mar = c(5, 4.5, 4, 2) + 0.1)
  x.axis <- graphics::barplot(sqrt(obs), col = "lightgray",
                    xlab = "y", ylab = expression(sqrt("Frequency")),
                    ylim = c(0, max(sqrt(obs), sqrt(esp)) + 0.5))
  points(x.axis, sqrt(esp), col = "red4", type = "b", pch = 16)
  par(op)
}

# Log-likelihood
#' @export
#' @rdname bergrm-methods
logLik.bergrm <- function(object) {
  ll <- object$logLik
  class(ll) <- "logLik"
  return(ll)
}

# AIC
#' @export
#' @rdname bergrm-methods
#' @param numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.
AIC.bergrm <- function(object, k = 2) {
  AIC <- - 2 * object$logLik + k * (object$p + object$k)
  class(AIC) <- "AIC"
  return(AIC)
}

# BIC
#' @export
#' @rdname bergrm-methods
BIC.bergrm <- function(object) {
  n <- object$n.obs
  BIC <- - 2 * object$logLik + log(n) * (object$p + object$k)
  return(BIC)
}

# Parameter estimates
#' @rdname bergrm-methods
#' @export
#' @param what a character indicating which parameter coefficients are
#'   required, parameters for the \code{"mean"} or for the
#'   \code{"dispersion"} model. If \code{"all"} (default), a list with
#'   coefficients for the \code{mean} and for the \code{dispersion}
#'   model is returned.
coef.bergrm <- function(object,
                       what = c("all", "mean", "dispersion")) {
  what <- match.arg(what)
  out <- switch(what,
                "all"        = list(
                  mean       = (object$coef)$mean,
                  dispersion = (object$coef)$disp),
                "mean"       = (object$coef)$mean,
                "dispersion" = (object$coef)$disp)
  return(out)
}

#  Variance-covariance matrix
#' @rdname bergrm-methods
#' @export
vcov.bergrm <- function(object) {
  return(object$vcov)
}

# Design matrices
#' @rdname bergrm-methods
#' @export
#' @param matrix a character indicating which model matrix is
#'   required, the model matrix for the mean (\code{"mean"}) or for the
#'   dispersion parameter (\code{"dispersion"}). If \code{"all"} (default), a list with
#'   with both matrices are returned.
model.matrix.bergrm <- function (object, matrix = c("all", "mean", "dispersion")) {

  what <- match.arg(what)
  out <- switch(what,
                "all" = list(mean = object$X,
                             dispersion = object$Z),
                "mean" = object$X,
                "dispersion" = object$Z)
  return(out)
}
