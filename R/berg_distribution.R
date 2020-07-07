#' @name berg
#'
#' @title The BerG Distribution
#'
#' @description Probability mass function, distribution function,
#' quantile function and random generation for the BerG distribution
#' with parameters \code{mu} and \code{phi}.
#'
#' @param x vector of non-negative integer quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param mu numeric; non-negative mean.
#' @param phi numeric; the dispersion index (greather than \code{abs(mu - 1)}),
#' @param lower.tail logical; if TRUE (default), probabilities are \code{P(X <= x)}, otherwise, \code{P(X > x)}.
#'
#' @return \code{dberg} returns the probability function, \code{pberg}
#' gives the distribution function, \code{qberg} gives the quantile function,
#' and \code{rberg} generates random observations.
#'
#' @details This set of functions represents the probability function, the cumulative distribution
#' function, quantile function and a random number generator for the BerG distribution
#' parameterized in terms of the mean and the dispersion index. This new
#' parameterization was proposed by Bourguignon, M. & Medeiros, R. (2020).
#'
#' @references Bourguignon, M. & Weiss, C. (2017). An INAR(1) process for modeling count
#' time series with equidispersion, underdispersion and overdispersion. Test, 26, 847--868.
#'
#' @references Bourguignon, M. & Medeiros, R. (2020). A simple and useful regression model for fitting count data
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @seealso For the populational skewness and kurtosis of the BerG distribution, use the functions \code{\link{skberg}} and \code{\link{ktberg}}, respectively.
#'
#' @examples
#' # BerG observations as categorical data:
#' Ni <- rberg(50, mu = 4, phi = 3.2); table(factor(Ni, 0:max(Ni)))
#' plot(prop.table(table(Ni)))
#' points(sort(unique(Ni)), dberg(sort(unique(Ni)), mu = 4, phi = 3.2), pch=16, col="red")
#'
#' # Probability function:
#'
#' # Parameters
#' mu = c(2.3, 3, 3.9); phi = 3
#'
#' plot(0:17-0.25,dberg(0:17,phi,mu[1]),type="h",xlab=" ",ylab="Pmf",
#' col="gray70",lwd=5);mtext("y",side=1,line=2.5)
#' segments(0:17,rep(0,length(0:17)),0:17,dberg(0:17,phi,mu[2]),col="gray40",lwd=5)
#' segments(0:17+0.25,rep(0,length(0:17)),0:17+0.25,dberg(0:17,phi,mu[3]),col="gray20",lwd=5)
#'
#' legend(13,0.32,legend =c(expression(mu==2.3),
#'           expression(mu==3), expression(mu==3.9)),
#'           lty=1,lwd=5,col=c("gray70","gray40","gray20"),bty="n")
#'
NULL

#' @rdname berg
#' @export
dberg <- function(x, mu, phi){
  if (any(x < 0))
    stop("x must be a nonnegative integer value")
  if ((any(mu<=0)) || (any(phi<=0)))
    stop("The parameters must be positives")
  if (any(phi < abs(mu - 1)))
    warning("Constraints are not satisfied")

  p0 <- (1 - mu + phi) / (1 + mu + phi)
  p  <- 4 * mu * ((mu + phi - 1)^(x[x > 0] - 1)) / ((mu + phi + 1)^(x[x > 0] + 1))

  prob <- c(rep(p0, sum(x == 0)),p)
  index <- c(which(x == 0), which(x > 0))

  return(prob[sort(index, index.return = T)$ix])
}

#' @rdname berg
#' @export
pberg <- function(q, mu, phi, lower.tail=TRUE){
  if ((any(mu <= 0)) || (any(phi <= 0)))
    stop("The parameters must be positives")
  if (any(phi < abs(mu - 1)))
    warning("Constraints are not satisfied")

  q <- floor(q)
  p <- (1 - mu + phi) / (1 + mu + phi) +
    (2 * mu / (1 + mu + phi)) * (1 - ((mu + phi - 1) /
                                        (mu + phi + 1))^q[q >= 0])

  prob <- c(rep(0, sum(q < 0)), p)
  index <- c(which(q < 0), which(q >= 0))

  prob <- prob[sort(index, index.return = TRUE)$ix]

  ifelse(lower.tail == TRUE, prob, 1 - prob)

}

#' @rdname berg
#' @export
qberg <- function(p, mu, phi, lower.tail = TRUE){
  if ((any(p < 0)) || (any(p > 1)))
    stop("p must be in the unit interval: (0, 1)")
  if ((any(mu <= 0)) || (any(phi <= 0)))
    stop("The parameters must be positives")
  if (any(phi < abs(mu - 1)))
    warning("Constraints are not satisfied")

  if(lower.tail == FALSE)
    p <- 1 - p

  p0 <- (1 - mu + phi)/(1 + mu + phi)

  ifelse(length(p) > 1, p.star <- p[p > p0], p.star <- p)
  ifelse(length(mu) > 1, mu.star <- mu[p > p0], mu.star <- mu)
  ifelse(length(phi) > 1, phi.star <- phi[p > p0], phi.star <- phi)

  q <- ceiling(
    round(log((1 - p.star) * (1 + mu.star + phi.star) / (2 * mu.star)) /
            log((mu.star + phi.star - 1) / (mu.star + phi.star + 1)), 2)
  )

  quanti <- c(rep(0, sum(p <= p0)), q)
  index <- c(which(p <= p0), which(p > p0))

  return(quanti[sort(index, index.return = TRUE)$ix])
}

#' @rdname berg
#' @export
rberg <- function(n, mu, phi){
  if ((any(mu<=0)) || (any(phi<=0)))
    stop("The parameters must be positives")
  if (any(phi < abs(mu - 1)))
    warning("Constraints are not satisfied")

  u <- stats::runif(n)
  return(qberg(u, mu, phi))
}


#' @name berg_sk_kt
#'
#' @title The BerG Skewness and Kurtosis
#'
#' @description Skewness and kurtosis for the BerG distribution
#' with parameters \code{mu} and \code{phi}.
#'
#' @param mu numeric; non-negative mean.
#' @param phi numeric; the dispersion index (greather than \code{abs(mu - 1)}).
#'
#' @return The functions \code{skberg} and \code{ktberg} returns the skewness and kurtosis of the BerG distribution
#' respectively.
#'
#' @references Bourguignon, M. & Medeiros, R. (2020). A simple and useful regression model for fitting count data
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @examples
#' # For mu = 2 and phi = 2.5, the skewness and kurtosis of the BerG
#' # distribution are, respectively,
#' skberg(mu = 2, phi = 2.5)
#' ktberg(mu = 2, phi = 2.5)
NULL

#' @rdname berg_sk_kt
#' @export
skberg <- function(mu, phi){
  if((any(mu<=0)) || (any(phi<=0)))
    warning("The parameters must be positives")
  if(any(phi < abs(mu - 1)))
    warning("Constraints are not satisfied")

  m3 <- mu*(3*(mu^2)+6*mu*phi+3*(phi^2)-1)/2
  return((m3 - 3*mu*(mu*phi) - mu^3)/((phi*mu)^(3/2)))
}
#' @rdname berg_sk_kt
#' @export
ktberg <- function(mu, phi){
  if((any(mu<=0)) || (any(phi<=0)))
    warning("The parameters must be positives")
  if(any(phi < abs(mu - 1)))
    warning("Constraints are not satisfied")

  m2 <- mu*(mu+phi)
  m3 <- mu*(3*(mu^2)+6*mu*phi+3*(phi^2)-1)/2
  m4 <- mu*(mu+phi)*(3*(mu^2)+6*mu*phi+3*(phi^2)-2)
  (m4 - 4*mu*m3 + 6*(mu^2)*m2 - 3*(mu^4))/((phi*mu)^2)
}
