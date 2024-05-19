#' @title Simulation of Species Age
#'
#' @description This function stochastically simulates the distribution of species age
#' given its branching time, speciation and extinction rate, and the sampling fraction
#'
#' @param branching_time time over which to run the simulation
#' @param lambda speciation rate
#' @param mu extinction rate
#' @param rho sampling fraction (Default: 1.0)
#' @param reps number of stochastic ages
#' @param batch batch size of how many phylogenies are simualted at once
#'
#' @details Some additional details
#'
#' @return A vector of stochastic species ages
#'
#' @seealso \code{\link{ageDensity}} for probabilistic age estimation
#'
#' @author Daniele Silvestro
#'
#' @references
#' Nee S., May R.M., and Harvey P.H. 1994.
#' The reconstructed evolutionary process.
#' Philos. Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
#'
#' @examples
#' simSpAge(branching_time = 3.0, lambda = 0.2, mu = 0.1, rho = 0.7, reps = 10000)


simSpAge <-function(branching_time, lambda, mu, rho = 1.0, reps, batch = 100){
  if(rho > 1) {
    stop("Sampling fraction should not be greater than 1")
  }
  sp_ages <- c()
  repeat {
    trees <- TreeSim::sim.bd.age(age = branching_time,
                                 lambda, mu,
                                 frac = rho,
                                 numbsim = batch,
                                 complete = TRUE,
                                 mrca = FALSE)
    for (i in 1:batch){
      tree_i <- trees[[i]]
      if (class(tree_i) == "numeric"){
        if (tree_i == 1 & runif(1) < rho){
          sp_ages <- c(sp_ages, branching_time)
        }
      }
      else{
        extant <- getExtant(tree_i)
        N <- sum(rbinom(length(extant), size = 1, prob = rho))
        if (N == 1){
          age_all <- calculateTipAges(tree_i)
          a <- age_all[age_all$tip == sample(extant, 1), 2]
          sp_ages <- c(sp_ages, a)
        }
      }
      if (length(sp_ages) >= reps){
        break
      }
    }
    if (length(sp_ages) >= reps){
      break
    }
  }
  return(sp_ages)
}
