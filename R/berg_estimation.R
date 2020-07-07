#' @name bergrm_estimation
#'
#' @title Maximum likelihood estimates for the BerG regression model
#'
#' @description Maximum likelihood estimates of the coefficients of a fit
#' of the BerG regression model.
#'
#' @param y a numeric vector of the response variable, in counts.
#' @param X,Z model matrices associated with the mean and the dispersion index, respectively.
#' @param link,link.phi character specification of the link function in the mean and in the dispersion index. The links \code{"log"} (default) \code{"sqrt"} and \code{"identity"} can be used.
#' @param control a list of control arguments specified via \code{berg_control}.
#' @param eq_constraint a function to evaluate equality constraints under the null hypothesis, in hypothesis testing contexts.
#' @param eq_constraint_jac a function to evaluate the jacobian of the equality constraints, which were passed by the \code{eq_constraint} argument.
#' @param ... arguments passed to \code{berg_control}.
#'
#' @return A list with the values of the maximum likelihood estimates (\code{est});
#' the value of the log likelihood at estimates (\code{logLik}); and a logical value,
#' where \code{TRUE} indicates that the point obtained is feasible (\code{feasible}).
#'
#' @references Bourguignon, M. & Medeiros, R. (2020). A simple and useful regression model for fitting count data.
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
#'
mle_berg <- function(y, X, Z, link = "log", link.phi = "log",
                     control = berg_control(...), eq_constraint = NULL,
                     eq_constraint_jac = NULL,...){

  # Control list
  start     <- control$start
  constant  <- control$constant
  error     <- control$error
  optimizer <- control$optimizer
  algorithm <- control$algorithm

  # Link funtions

  # Dispersion link function
      g1 <- g(link)$fun
  g1.inv <- g(link)$inv
     g1. <- g(link)$deriv.

  # Dispersion link function
      g2 <- g(link.phi)$fun
  g2.inv <- g(link.phi)$inv
     g2. <- g(link.phi)$deriv.


  # Design matrices and necessary quantities
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  p <- NCOL(X)
  k <- NCOL(Z)
  n <- as.numeric(length(y))

  if(is.null(start)){
        bs <- solve(t(X)%*%X)%*%t(X)%*%g1(y + 0.1)
       mu. <- g1.inv(X%*%bs)
    sigma. <- stats::var(g1(y + 0.1))/(g1.(mu.)^2)
      phi. <- sigma. / mu. + abs(mu. - 1)
        gs <- solve(t(Z)%*%Z)%*%t(Z)%*%g2(phi.)
     start <- c(bs, gs)
  }

  # Log-likelihood
  ll <- function(par) -ll_berg(par, y, X, Z, link, link.phi)

  # Score function
  U <- function(par) -U_berg(par, y, X, Z, link, link.phi)

  # Constraints and its Jacobian
  hj <- function(par){


    # Relations
    mu  <- g1.inv(X%*%par[1:p])
    phi <- g2.inv(Z%*%par[(p + 1):(p + k)])

    # Function h (h (theta) < 0 if theta is admissible)
    hj <- rbind(- phi + mu - 1 + constant, - phi - mu + 1 + constant)

    return(hj)

  }

  Jh <- function(par){

    # Relations
    mu  <- g1.inv(X%*%par[1:p])
    phi <- g2.inv(Z%*%par[(p + 1):(p + k)])

    # Diagonal matrix
    D1 <- diag(as.numeric(1 / g1.(mu)))
    D2 <- diag(as.numeric(1 / g2.(phi)))

    # Jacobian matrix
    Jh <- matrix(NA, 2 * n, p + k)
    Jh[1:n, 1:p] <- D1%*%X
    Jh[1:n, (p + 1):(p + k)] <- - D2%*%Z
    Jh[(n + 1):(2 * n), 1:p] <- - D1%*%X
    Jh[(n + 1):(2 * n),(p + 1):(p + k)] <- - D2%*%Z


    return(Jh)
  }


  if (optimizer == "nloptr"){

    if (!is.null(eq_constraint)){
      eq_constraint <- match.fun(eq_constraint)
      eq_constraint_jac <- match.fun(eq_constraint_jac)
    }

    est <- nloptr::nloptr(x0 = start,
                          eval_f = ll,
                          eval_grad_f = U,
                          eval_g_ineq = hj,
                          eval_jac_g_ineq = Jh,
                          eval_g_eq = eq_constraint,
                          eval_jac_g_eq = eq_constraint_jac,
                          opts = list("algorithm" = algorithm,
                                      "xtol_rel"= error))$solution

    logLik <- -ll(est)

    feasible <- all(hj(est) < 0)
  }

  if (optimizer == "optim"){

    est <- suppressWarnings(stats::optim(par = start,
                                         fn = ll,
                                         gr = U,
                                         method = "BFGS")$par)
    logLik <- -ll(est)

    feasible <- all(hj(est) < 0)
  }

  return(list(est = est, logLik = logLik, feasible = feasible))
}
