#' @title Get Phylogenetic Ages for the Tips of a Phylogeny

calculateTipAges <- function(tree) {
  phy.age <- picante::node.age(tree)
  BL.position <- cbind(phy.age$edge, phy.age$age, tree$edge.length)
  dist.tip <- max(phy.age$age) - BL.position[, 3]
  BL.positions <- cbind(BL.position, dist.tip)
  ages <- BL.positions[, 5] + BL.positions[, 4]
  BL.positions <- cbind(BL.positions,ages)
  node.ages <- as.data.frame(BL.positions)
  colnames(node.ages) <- c("parental.node","daughter.node",
                           "dist.root","BL",
                           "dist.tip","mrca.age")
  ## node.ages is a data frame listing as variables the identity of parental and
  # daughter nodes, the distance from the root and from the present of each node,
  # the branch length and the age of the most recent common ancestor
  species.ages <- picante::node.ages[node.ages[,2] < length(tree$tip) + 1, ]
  rownames(species.ages) <- tree$tip[species.ages$daughter.node]
  # species ages is node.ages data frame reduced to the tips (species)
  species.ages <- species.ages[order(row.names(species.ages)), ]
  output.table <- as.data.frame(cbind(row.names(species.ages),
                                      species.ages$mrca.age))
  colnames(output.table) <- c('tip','tip.age')
  output.table$tip.age <- as.numeric(output.table$tip.age)
  return(output.table)
}