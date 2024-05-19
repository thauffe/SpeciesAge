#' @title Esitmate Speciation and Extinction rate from phylogeny
#'
#' @description This function estimates speciation and extinction rate from a phylogeny under a given extinction and sampling fraction
#'
#' @param phy An ultrametric bifurcating phylogenetic tree, in ape "phylo" format.
#' @param epsilon Extinction fraction
#' @param rho Sampling fraction
#' @param ml_optim Method to use for optimisation (Default: "subplex"). May be one of "optim", "subplex", "nlminb", "nlm" (partial unambigious string is allowed).
#'
#' @author Torsten Hauffe
#'
#' @return A named vector of two parameters, lambda and mu.
#'
#' @examples
#' estimateBD(phy = , epsilon = 0.5)


estimateBD <- function(phy, epsilon = 0, rho = 1, ml_optim = "subplex") {
  if(epsilon >= 1) {
    stop("Extinction fraction should not be greater than 1")
  }
  # Birth-death likelihood
  bd_lik <- make.bd(tree = phy, sampling.f = rho)
  # Constrain extinction fraction
  con <- paste0("mu ~ ", epsilon, " * lambda")
  bd_lik <- constrain(bd_lik, con)
  # Initial birth rate for ML search with yule rate
  b_init <- phy$Nnode / sum(phy$edge.length)
  if (epsilon > 0) {
    b_init <- b_init / epsilon
  }
  # Finde maximum likelihood rates
  bd_fit <- find.mle(bd_lik, x.init = b_init, method = ml_optim)
  # Format birth and death rate for output
  lambda <- coef(bd_fit)
  mu <- epsilon * lambda
  bd_rates <- c(lambda, mu)
  names(bd_rates) <- c("lambda", "mu")
  return(bd_rates)
}
