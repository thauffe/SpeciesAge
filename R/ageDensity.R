#' @title Probability Distribution of Species Age
#'
#' @description This function provides the probability for a species age
#' given speciation and extinction rate, and the sampling fraction
#'
#' @param t branching time
#' @param lambda speciation rate
#' @param mu extinction rate
#' @param rho sampling fraction (Default: 1.0)
#'
#' @details Some additional details
#'
#' @seealso \code{\link{simSpAge}} for stochastic age simulation
#'
#' @author Daniele Silvestro
#'
#' @examples
#' ageDensity(lambda = 0.2, mu = 0.1, t = 3.0, rho = 0.7)
#'
#' @export ageDensity

ageDensity <- function(t, lambda, mu, rho = 1.0) {
  # get the integrated rate
  integrated_rate <- integrateRate(lambda, mu, 0, t, rho)
  # compute the density of the age
  density <- 2 * lambda * p0t(t, lambda, mu, rho) * dpois(0, integrated_rate)
  return(density)
}


p0t <- function(t, lambda, mu, rho) {
  1 - (rho * (lambda - mu) / (rho * lambda + (lambda * (1 - rho) - mu) * exp((mu - lambda) * t)))
}


integrateRate <- function(lambda, mu, t1, t2, rho) {
  # the definite integral, Lambda(t1, t2)
  2 * (mu * t2 - log(lambda * rho - (lambda * (rho - 1) + mu) * exp((mu - lambda) * t2))) -
    2 * (mu * t1 - log(lambda * rho - (lambda * (rho - 1) + mu) * exp((mu - lambda) * t1)))
}


