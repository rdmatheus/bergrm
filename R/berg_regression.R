#' @name bergrm
#'
#' @title BerG Regression for Count Data
#'
#' @description Fit of the BerG regression model via maximum
#' likelihood for a new parameterization of this distribution
#' that is indexed by the mean and the dispersion index.
#'
#' @param formula a symbolic description of the model, of type y ~ x or
#'  y ~ x | z.See details below.
#' @param data an optional data frame containing the variables in the formula.
#'  By default the variables are taken from environment(formula).
#' @param link,link.phi character; specification of the link function in the
#'  mean model. The links \code{"log"}, \code{"sqrt"}, and \code{"identity"}
#'  can be used. For the mean model, the default link function is "log".
#'  For the dispersion index model, the default link function is also "log"
#'  unless formula is of type y ~ x where the default is "identity".
#' @param y a numeric vector of the response variable, in counts. This argument
#'  is used only in the \code{berg_fit} function.
#' @param X,Z model matrices associated with the mean and the dispersion index
#'  parameters, respectively, which are used only in the \code{berg_fit} function.
#' @param disp.test logical; if TRUE, the function \code{glm.bg} returns the
#'  test for constant dispersion.
#' @param control a list of control arguments specified via \code{berg_control}
#' (under development).
#' @param ... arguments passed to \code{berg_control} (under development).
#'
#' @return  The \code{glm.bg} function returns an object of class "bergrm",
#'  which consists of a list with the following components. \code{berg_fit} is
#'  an auxiliary function that returns an unclassed list formed by some of the
#'  components below. This function can be used to fit the BerG regression
#'  in terms of the y, X, and Z objects instead of the formula.
#' \describe{
#'   \item{coefficients}{a list containing the elements "mean" and
#'         "dispersion" that consist of the estimates of the
#'          coefficients associated with the mean and the dispersion index,
#'          respectively.}
#'   \item{link, link.phi}{the link function used for the mean model, and for
#'    the dispersion index model, respectively.}
#'   \item{vcov}{asymptotic covariance matrix of the maximum likelihood
#'    estimator of the model parameters vector. Specifically, the inverse of
#'    the Fisher information matrix.}
#'   \item{logLik}{log-likelihood of the fitted model.}
#'   \item{AIC}{model Akaike information criteria.}
#'   \item{BIC}{model bayesian information criteria.}
#'   \item{n.obs, p, k}{Sample size, number of coefficients in the mean model,
#'    and number of coefficients in the dispersion index model, respectively.}
#'   \item{feasible}{Logical. If \code{TRUE}, the estimates obtained belong to
#'    the parametric space.}
#'   \item{pearson.residuals}{a vector with the Pearson residuals.}
#'   \item{residuals}{a vector with the randomized quantile residuals.}
#'   \item{fitted.values}{a vector with the fitted means.}
#'   \item{phi.hat}{a vector with the fitted dispersion indexes.}
#'   \item{linear.predictors}{a list containing the elements "eta1" and
#'         "eta2" that consist of the fitted linear predictor for the mean and
#'          the dispersion index.}
#'   \item{response}{the vector of the response.}
#'   \item{X, Z}{model matrices associated with the mean and the dispersion
#'    parameter, respectively.}
#'   \item{call}{the function call.}
#'   \item{formula}{the formula used to specify the model in \code{glm.bg}.}
#'  }
#'
#' @details The basic formula is of type \code{y ~ x1 + x2 + ... + xp} which
#'  specifies the model for the mean response only. Following the syntax of the
#'  \code{betareg} package (Cribari-Neto and Zeileis, 2010), the model for the dispersion index, say in terms of
#'  z1, z2, ..., zk, is specified as \code{y ~ x1 + x2 + ... + xp | z1 + z2 +
#'  ... +zk} using functionalities inherited from package \code{Formula}
#'  (Zeileis and Croissant, 2010).
#'
#' @references Bourguignon, M. & Medeiros, R. (2020). A simple and useful regression model for fitting count data.
#' @references Cribari-Neto F, Zeileis A (2010). Beta Regression in R. Journal
#'  of Statistical Software, 34, 1–24
#' @references Zeileis A, Croissant Y (2010). Extended Model Formulas in
#'  R: Multiple Parts and Multiple Responses. Journal of Statistical Software,
#'  34, 1–13.
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @examples
#' # Dataset: grazing
#'
#' data(grazing)
#' attach(grazing); head(grazing)
#'
#' # Response variable (Number of understorey birds)
#' plot(table(birds), xlab = "Number of understorey birds", ylab = "Frequency")
#'
#' # Explanatory variables
#' boxplot(split(birds, when), ylab = "Number of understorey birds",
#'         xlab = "When the bird count was conduct", pch = 16)
#' boxplot(split(birds, grazed), ylab = "Number of understorey birds",
#'         xlab = " Which side of the stockproof fence", pch = 16)
#'
#'
#' # Fit of the BerG regression model with varying dispersion
#' fit_disp <- glm.bg(birds ~ when + grazed | when + grazed, data = grazing)
#' summary(fit_disp)
#'
#' # Fit of the BerG regression model with constant dispersion
#' fit <- glm.bg(birds ~ when + grazed, data = grazing, link.phi = "identity")
#' summary(fit)
#'
#' # Diagnostic
#' plot(fit)
#' envel_berg(fit)
#'
#' detach(grazing)
NULL

###################################################################
# Link function                                                   #
##################################################################
g <- function(link){

  switch(link,

         identity = {
           fun <- function(theta) theta
           inv <- function(eta) eta
           deriv. <- function(theta) rep.int(1, length(theta))
           deriv.. <- function(theta) rep.int(0, length(theta))
           valideta <- function(eta) TRUE
         },

         log = {
           fun <- function(theta) log(theta)
           inv <- function(eta) pmax(exp(eta), .Machine$double.eps)
           deriv. <- function(theta) 1 / theta
           deriv.. <- function(theta) -1 / (theta ^ 2)
           valideta <- function(eta) TRUE
         },

         sqrt = {
           fun <- function(theta) sqrt(theta)
           inv <- function(eta) eta^2
           deriv. <- function(theta) 1 / (2 * sqrt(theta))
           deriv.. <- function(theta) -1 / (4 * (theta ^ (3 / 2)))
           valideta <- function(eta) all(is.finite(eta)) && all(eta > 0)
         },

         stop(gettextf("link %s not available", sQuote(link)), domain = NA))

  environment(fun) <- environment(inv) <- environment(deriv.) <-
    environment(deriv..) <- environment(valideta) <- asNamespace("stats")

  structure(list(fun = fun, inv = inv, deriv. = deriv.,
                 deriv.. = deriv.., valideta = valideta,
                 name = link), class = "link-sdlrm")
}


###################################################################
# Log-lokelihood function                                         #
##################################################################
ll_berg <- function(par, y, X, Z, link = "log", link.phi = "log"){

  # Link functions
  g1.inv <- g(link)$inv
  g2.inv <- g(link.phi)$inv

  # Design matrices and necessary quantities
  X <- as.matrix(X); Z <- as.matrix(Z)
  p <- NCOL(X); k <- NCOL(Z)

  # Relations
   mu <- g1.inv(X%*%par[1 : p])
  phi <- g2.inv(Z%*%par[(p + 1):(p + k)])

  l0 <- as.numeric(log(1 - mu[y == 0] + phi[y == 0]) - log(1 + mu[y == 0]+phi[y == 0]))
  l  <- as.numeric(log(4 * mu[y > 0]) + (y[y > 0] - 1)*log(mu[y > 0]+phi[y > 0] - 1) -
                     (y[y > 0] + 1)*log(mu[y > 0]+phi[y > 0] + 1))

  return(sum(c(l0,l)))
}


###################################################################
# Score function                                                  #
##################################################################
U_berg <- function(par, y, X, Z, link = "log", link.phi = "log"){

  # Link functions
  g1.inv <- g(link)$inv
     g1. <- g(link)$deriv.
    g1.. <- g(link)$deriv..

  g2.inv <- g(link.phi)$inv
     g2. <- g(link.phi)$deriv.
    g2.. <- g(link.phi)$deriv..

  # Design matrices and necessary quantities
  X <- as.matrix(X); Z <- as.matrix(Z)
  p <- NCOL(X); k <- NCOL(Z)
  delta <- (as.numeric(y == 0))

  # Relations
   mu <- g1.inv(X%*%par[1 : p])
  phi <- g2.inv(Z%*%par[(p + 1):(p + k)])

  u1 <- - 2 * delta * (1 + phi)/((1 - mu + phi) * (1 + mu + phi)) +
    (1 - delta)*( 1 / mu + 2 * (y - mu - phi) / ((mu + phi - 1)*(mu + phi + 1)) )

  u2 <- 2 * delta * mu/((1 - mu + phi)*(1 + mu + phi)) +
    2 * (1 - delta) * (y - mu - phi) / ((mu + phi - 1) * (mu + phi + 1))

  D1 <- diag(as.numeric(1/g1.(mu)))
  D2 <- diag(as.numeric(1/g2.(phi)))

  Ub <- t(X)%*%D1%*%u1; Ug = t(Z)%*%D2%*%u2

  U <- c(Ub, Ug)

  return(U)
}

###################################################################
# Fisher information matrix                                       #
##################################################################
K_berg <- function(par, X, Z, link = "log", link.phi = "log"){

  # Link functions
  g1.inv <- g(link)$inv
     g1. <- g(link)$deriv.
    g1.. <- g(link)$deriv..

  g2.inv <- g(link.phi)$inv
     g2. <- g(link.phi)$deriv.
    g2.. <- g(link.phi)$deriv..

  # Design matrices and necessary quantities
  X = as.matrix(X); Z = as.matrix(Z)
  p = NCOL(X); k = NCOL(Z)

  # Relations
  mu =  g1.inv(X%*%par[1 : p])
  phi = g2.inv(Z%*%par[(p + 1):(p + k)])

  D1 <- diag(as.numeric(1 / g1.(mu)))
  D2 <- diag(as.numeric(1 / g2.(phi)))

  K <- matrix(rep(NA, (p + k) * (p + k)), p + k, p + k)

  W1 <- diag(-as.numeric(
    (- 4 * mu * (1 + phi) / ((1 - mu + phi) * ((mu + phi + 1)^3)) - 2 /
       (mu * (1 + mu + phi)) + 4 * mu * ((mu + phi)^2 + 1) /
         (((mu + phi - 1)^2) * ((mu + phi + 1)^3)) - 2 * mu * (mu + phi) /
            (((mu + phi - 1)^2) * ((mu + phi + 1)^2))) * (1 / g1.(mu)) -
              (- 2 * (1 + phi) / ((mu + phi + 1)^2) + 2 / (1 + mu + phi) -
                 4 * mu * (mu + phi) / ((mu + phi - 1) * ((mu + phi + 1)^2)) +
                   2 * mu / ((mu + phi - 1) * (mu + phi + 1))) *
                     (g1..(mu) / (g1.(mu)^2))
  ))

  W2 <- diag(-as.numeric(
    (2 * (mu^2) + 2*((1 + phi)^2)) / ((1 - mu + phi)*((1 + mu + phi)^3)) -
      4 * mu * (mu + phi) / (((mu + phi - 1)^2) * ((mu + phi + 1)^2)) +
      4 * mu * ((mu + phi)^2 + 1) / (((mu + phi - 1)^2) * ((mu + phi + 1)^3))
  ))

  W3 <- diag(-as.numeric(
    (- 4 * mu * (1 + phi) / ((1 - mu + phi) * ((1 + mu + phi)^3)) -
       4 * mu * (mu + phi)/(((mu + phi-1)^2) * ((mu + phi + 1)^2))  +
          4 * mu * ((mu + phi)^2  +  1) / (((mu + phi - 1)^2) *
            ((mu + phi + 1)^3))) * (1 / g2.(phi)) - (2 * mu /
               ((1 + mu + phi)^2) + 2 * mu / ((mu + phi - 1) *
                  (mu + phi + 1)) - 4 * mu * (mu + phi) / ((mu + phi - 1) *
                    ((mu + phi + 1)^2))) *(g2..(phi) / (g2.(phi)^2))
  ))

  Kbb <- t(X)%*%W1%*%D1%*%X
  Kbg <- t(X)%*%D1%*%W2%*%D2%*%Z
  Kgb <- t(Kbg)
  Kgg <- t(Z)%*%W3%*%D2%*%Z

  K[1:p,1:p] <- Kbb; K[1:p,(p + 1):(p + k)] <- Kbg; K[(p + 1):(p + k),(1:p)] <- Kgb
  K[(p + 1):(p + k),(p + 1):(p + k)] <- Kgg

  return(K)

}

###################################################################
# Control optimization function in nloptr                         #
##################################################################
berg_control <- function(start = NULL,
                         start2 = NULL,
                         constant = 1e-7,
                         error = 1e-8,
                         optimizer = "nloptr",
                         algorithm = "NLOPT_LD_SLSQP", ...){

  rval <- list(start = start,
               start2 = start,
               constant = constant,
               error = error,
               optimizer = optimizer,
               algorithm = algorithm)

  rval <- c(rval, list(...))
  return(rval)
}

#' @rdname bergrm
#' @export
glm.bg <- function(formula, data, link=c("log", "sqrt", "identity"),
                   link.phi = NULL, disp.test = FALSE,
                   control = berg_control(...), ...)
{
  cl <- match.call()
  if (missing(data))  data <- environment(formula)

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  oformula <- stats::as.formula(formula)
  formula <- Formula::as.Formula(formula)

  if (length(formula)[2L] < 2L) {
    formula <- Formula::as.Formula(formula(formula), ~ 1)
    simple_formula <- TRUE
  }else {
    if (length(formula)[2L] > 2L) {
      formula <- Formula::Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
    simple_formula <- FALSE
  }

  mf$formula <- formula
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- stats::terms(formula, data = data)
  mtX <- stats::terms(formula, data = data, rhs = 1L)
  mtZ <- stats::delete.response(stats::terms(formula, data = data, rhs = 2L))
  y <- stats::model.response(mf, "numeric")
  X <- stats::model.matrix(mtX, mf)
  Z <- stats::model.matrix(mtZ, mf)
  p <- NCOL(X); k <- NCOL(Z)

  if (length(y) < 1)
    stop("empty model")
  if (any(y < 0))
    stop("invalid dependent variable, all observations must be non-negatives integers")

  if (is.character(link)){
    link <- match.arg(link)
  }

  if (is.null(link.phi)){
    link.phi <- if (k == 1) "identity" else "log"
  }

  out <- berg_fit(y, X, Z, link = link, link.phi = link.phi, control = control)
  out$call <- cl
  out$formula <- formula
  out$names.mean <- c("(Intercept)", colnames(X)[2:p])

  if (k > 1) {

    out$names.dispersion <- c("(Intercept)", colnames(Z)[2:k])

    if (disp.test == TRUE){
      start2 <- control$start2
      out$test <- round(disp_test(y, X, Z, cols = 2:k, link = link,
                                  link.phi = link.phi, start = start2), 5)
    }

  }else{

    if (out$link.phi == "identity"){
      out$names.dispersion <- "phi"
    }else{
      out$names.dispersion <- "g(phi)"
    }

    out$test <- NULL
  }

  est <- c((out$coe)$mean, (out$coe)$disp)
  out$vcov <- solve(K_berg(est, X, Z, link = link, link.phi = link.phi))

  class(out) <- "bergrm"
  out
}

###################################################################
# Fit function                                                    #
##################################################################
#' @rdname bergrm
#' @export
berg_fit <- function(y, X, Z = NULL, link = "log", link.phi = "log",
                     control = berg_control(...), ...)
{

  n <- length(y); delta <- (as.numeric(y == 0))
  if (is.null(Z)) Z <- matrix(1, nrow = n, ncol = 1)

  # Vector lenght
  p <- NCOL(X); k <- NCOL(Z)

  # Link functions
  g1     <- g(link)$fun
  g1.inv <- g(link)$inv
  g1.    <- g(link)$deriv.

  g2     <- g(link.phi)$fun
  g2.inv <- g(link.phi)$inv
  g2.    <- g(link.phi)$deriv.

  mle <- mle_berg(y, X, Z, link = link, link.phi = link.phi,
                  control = control)

  # Estimates
  est <- mle$est

  # Coefficients
  beta <- est[1:p]; gama <- est[(p + 1):(p + k)]

  # Feasible?
  feasible <- mle$feasible

  # Fit
  mu  <- g1.inv(X%*%beta)
  phi <- g2.inv(Z%*%gama)

  # Loglikelihood and informations criteria
  logLik <- mle$logLik
  AIC <- - 2 * logLik + 2 * (p + k)
  BIC <- - 2 * logLik + log(n) * (p + k)

  # Pearson residuals
  rp = as.numeric((y - mu) / sqrt(phi * mu))

  # Randomized quantile residuals
  rq = rqr_berg(y, mu, phi)

  rval <- list(coefficients = list(mean = beta, dispersion = gama),
               link = link, link.phi = link.phi,
               logLik = logLik, AIC = AIC,
               BIC = BIC, n.obs = n, p = p, k = k, feasible = feasible, pearson.residuals = rp, residuals = rq,
               fitted.values = structure(mu, .Names = names(y)),
               linear.predictors=list(eta1 = g(link)$fun(mu), eta2 = g(link.phi)$fun(phi)),
               phi.hat = phi, response=y, Z=Z, X=X)
  rval
}
