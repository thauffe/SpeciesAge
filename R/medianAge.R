#' @title Expected Median Species Age
#'
#' @description This function calculates the expected median species age.
#'
#' @param v branching time
#' @param lambda speciation rate
#' @param mu extinction rate
#' @param rho sampling fraction (Default: 1.0)
#'
#' @author Daniele Silvestro
#'
#' @examples
#' medianAge(v = 5.0, lambda = 0.3, mu = 0.1, rho = 0.7)

#' @export medianAge

medianAge <- function(lambda, mu, v, rho) {
  if (rho > 1) {
    stop("Sampling fraction should not be greater than 1")
  }
  # get the integrated rate
  integrated_rate <- integrateRate(lambda, mu, 0, v, rho)
  # compute the probability of zero events
  p_no_events <- dpois(0, integrated_rate)
  # compute the median by optimization
  median_age <- optim(par = v / 2, fn = function(t) {
    # compute the probability older than t
    p <- integrate(ageDensity,
                   lambda, mu, rho,
                   lower = t, upper = v)$value + p_no_events
    # compute the distance from the median
    return(abs(p - 0.5))
  },
  v,
  lower = 0, upper = v, method = "Brent")$par
  return(median_age)
}
