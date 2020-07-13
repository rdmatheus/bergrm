mean_test <- function(y, X, Z, cols, beta0 = NULL,
                      link = "log", link.phi = "log",
                      control = berg_control(...), ...){

  # Necessary quantities
  X <- as.matrix(X); Z <- as.matrix(Z)
  p <- NCOL(X); k <- NCOL(Z);
  n <- length(y); r <- length(cols)

  if (is.null(beta0)){
    beta0 <- as.matrix(rep(0, r))
  }else{
    beta0 <- as.matrix(beta0)
  }

  # Partition in X mean model matrix
  X1 <- as.matrix(X[, -cols]); X2 <- as.matrix(X[, cols])
  X_aux <- cbind(X1, X2)

  # Unrestricted mle
  theta_hat <- mle_berg(y, X_aux, Z, link, link.phi, control)$est

  # Restricted mle
  eq_constraint <- function(theta){
    return(c(rep(0, p - r), rep(1, r), rep(0, k))%*%theta -
             beta0)
  }

  eq_constraint_jac <- function(theta){
    return(c(rep(0, p - r), rep(1, r), rep(0, k)))
  }

  theta_tilde <- mle_berg(y, X_aux, Z, link, link.phi, control,
                         eq_constraint = eq_constraint,
                         eq_constraint_jac = eq_constraint_jac)$est

  # Tests statistics
  Ub2 <- U_berg(theta_tilde, y, X_aux, Z, link, link.phi)[(p - r + 1):p]
  Kb2 <- solve(K_berg(theta_tilde, X_aux, Z, link, link.phi))[(p - r + 1):p, (p - r + 1):p]

  # Score
  S <- t(Ub2)%*%Kb2%*%Ub2

  # Wald
  W <- t(theta_hat[(p - r + 1):p] - beta0)%*%solve(
    solve(K_berg(theta_hat, X_aux, Z, link, link.phi))[(p - r + 1):p, (p - r + 1):p]
  )%*%(theta_hat[(p - r + 1):p] - beta0)

  # Likelihood ratio
  LR <- 2 * (ll_berg(theta_hat, y, X_aux, Z, link, link.phi) -
               ll_berg(theta_tilde, y, X_aux, Z, link, link.phi))

  # Gradient
  G <- t(Ub2)%*%(theta_hat[(p - r + 1):p] - beta0)

  #t(U_berg(theta_tilde, y, X_aux, Z, link, link.phi))%*%
  #            (theta_hat - theta_tilde)

  out <- as.matrix(c(S, W, LR, G))
  rownames(out) <- c("Score","Wald","Lik. Ratio","Gradient")
  colnames(out) <- "Value"

  return(t(out))
}


disp_test <- function(y, X, Z, cols, gamma0 = NULL,
                      link = "log", link.phi = "log", control = berg_control(...), ...){

  # Necessary quantities
  X <- as.matrix(X); Z <- as.matrix(Z)
  p <- NCOL(X); k <- NCOL(Z);
  n <- length(y); q <- length(cols)

  if (is.null(gamma0)){
    gamma0 <- as.matrix(rep(0, q))
  }else{
    gamma0 <- as.matrix(gamma0)
  }

  # Partition in X mean model matrix
  Z1 <- as.matrix(Z[, -cols]); Z2 = as.matrix(Z[, cols])
  Z_aux <- cbind(Z1, Z2)

  # Unrestricted mle
  theta_hat <- mle_berg(y, X, Z_aux, link, link.phi, control)$est

  # Restricted mle
  eq_constraint <- function(theta){
    return(cbind(matrix(0, q, p + k - q), diag(rep(1, q)))%*%theta -
             gamma0)
  }

  eq_constraint_jac <- function(theta){
    return(cbind(matrix(0, q, p + k - q), diag(rep(1, q))))
  }

  theta_tilde <- mle_berg(y, X, Z_aux, link, link.phi, control,
                         eq_constraint = eq_constraint,
                         eq_constraint_jac = eq_constraint_jac)$est


  #Ug2 <- U_berg(theta_tilde, y, X, Z_aux, link, link.phi)[(p + k - q + 1):(p + k)]
  #Kg2 <- solve(K_berg(theta_tilde, X, Z_aux, link, link.phi))[(p + k - q + 1):(p + k), (p + k - q + 1):(p + k)]

  # Tests statistics

  # Score
  S <- t(U_berg(theta_tilde, y, X, Z_aux, link, link.phi))%*%
         solve(K_berg(theta_tilde, X, Z_aux, link, link.phi))%*%
         U_berg(theta_tilde, y, X, Z_aux, link, link.phi)

  # Wald
  W <- t(theta_hat[(p + k - q + 1):(p + k)] - gamma0)%*%solve(
    solve(K_berg(theta_hat, X, Z_aux, link, link.phi))[(p + k - q + 1):(p + k), (p + k - q + 1):(p + k)]
  )%*%(theta_hat[(p + k - q + 1):(p + k)] - gamma0)

  # Likelihood ratio
  LR <- 2 * (ll_berg(theta_hat, y, X, Z_aux, link, link.phi) -
               ll_berg(theta_tilde, y, X, Z_aux, link, link.phi))

  # Gradient
  G <- t(U_berg(theta_tilde, y, X, Z_aux, link, link.phi))%*%
              (theta_hat - theta_tilde)

  #t(Ug2)%*%(theta_hat[(p + k - q + 1):(p + k)] - gamma0)

  out <- as.matrix(c(S, W, LR, G))
  rownames(out) <- c("Score","Wald","Lik. Ratio","Gradient")
  colnames(out) <- "Value"

  return(t(out))
}
