#' @title Esitmate Speciation and Extinction rate from phylogeny
#'
#' @description This function estimates speciation and extinction rate from a phylogeny under a given extinction and sampling fraction
#'
#' @param tree An ultrametric bifurcating phylogenetic tree, in ape "phylo" format.
#' @param epsilon Extinction fraction (Default: 0.0)
#' @param rho Sampling fraction (Default: 1.0)
#' @param ml_optim Method to use for optimisation (Default: "subplex"). May be one of "optim", "subplex", "nlminb", "nlm" (partial unambigious string is allowed).
#'
#' @details
#' This function estimates speciation and extinction rate given a fixed
#' extinction fraction (i.e. extinction rate divided by speciation rate).
#' Calculation is based on the reconstructed birth-death process (Nee et al., 1994),
#' implemented in the make.bd function of the diversitree package (FitzJohn, 2012).
#'
#' @author Torsten Hauffe
#'
#' @references
#' FitzJohn R.G. 2012.
#' diversitree: comparative phylogenetic analyses of diversification in R.
#' Methods Ecol. Evol. 3(6): 1084-1092.
#'
#' Nee S., May R.M., and Harvey P.H. 1994.
#' The reconstructed evolutionary process.
#' Philos. Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
#'
#' @return A named vector of two parameters, lambda and mu.
#'
#' @examples
#' estimateBD(tree = , epsilon = 0.5)


estimateBD <- function(tree, epsilon = 0.0, rho = 1.0, ml_optim = "subplex") {
  if(epsilon >= 1) {
    stop("Extinction fraction should not be greater than 1")
  }
  # Birth-death likelihood
  bd_lik <- make.bd(tree = tree, sampling.f = rho)
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
