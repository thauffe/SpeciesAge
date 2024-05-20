#' @title Expected Mean Species Age
#'
#' @description This function calculates the expected mean species age.
#'
#' @param v Branching time
#' @param lambda Speciation rate
#' @param mu Extinction rate
#' @param rho Sampling fraction (Default: 1.0)
#'
#' @author Daniele Silvestro
#'
#' @examples
#' meanAge(v = 3.0, lambda = 0.2, mu = 0.1, rho = 0.7)

#' @export meanAge

meanAge <- function(v, lambda, mu, rho = 1.0) {
  if (rho > 1) {
    stop("Sampling fraction should not be greater than 1")
  }
  # get the integrated rate
  integrated_rate <- integrateRate(lambda, mu, 0, v, rho)
  # compute the probability of zero events
  p_no_events <- dpois(0, integrated_rate)
  # compute the mean, given that the age is not v
  mean_with_events <- integrate(tAgeDensity,
                                lambda, mu, rho,
                                lower = 0, upper = v)$value
  # combine the probabilities
  mean_age <- p_no_events * v + mean_with_events
  return(mean_age)
}

tAgeDensity <- function(t, lambda, mu, rho) {
  t * ageDensity(t, lambda, mu, rho)
}
