#' @title retrodictNodes
#' @name retrodictNodes
#' @description Given two trees with support values, creates a dataframe containing support values of tree 1 and the occurrence of the clade in tree 2. This dataframe can be used for logistic regressions testing whether support values of tree 1 retrodict the occurrence of clades in tree 2.
#' @author Daniel YM Nakamura
#'
#' @param tree1 A .phylo tree that can be loaded using ape::read.tree for NEWICK files or TreeTools::ReadTntTree for TNT files. The .phylo must contain $node.label.
#' @param tree2 Another .phylo tree.
#' @param outgroup Optional. Specify outgroup taxa to remove (by default, the function assumes that the user does not want to remove outgroup taxa)
#' @param root Optional. Specify the same root for both trees, which is recommended to facilitate tree comparisons (by default, the function assumes that trees share the same root)
#' @param dataframe Optional. Write a TSV file in current directory containing the output dataframe (by default, no .TSV is written).
#' @examples
#' # Example 1 (simplest case)
#' tree1 = read.tree (text="(t1,(t2,(t3,(t4,t5)75)32)45);")
#' tree2 = read.tree (text="(t1,(t6,(t3,(t4,t5)47)53)94);")
#' retrodictNodes (tree1, tree2)
#'
#' @export
retrodictNodes = function (tree1,tree2,
                        outgroup=NULL,
                        root=NULL,
                        dataframe=F){
  # Initial warnings
  missing_params <- c()
  if (is.null(tree1)) missing_params <- c(missing_params, "tree1")
  if (is.null(tree2)) missing_params <- c(missing_params, "tree2")
  if (length(missing_params) > 0) {   # Check if there are any missing parameters and print a message
    message("The following parameters are missing: ", paste(missing_params, collapse = ", "))
  } else {
    message("All required parameters provided.") # Proceed with the main functionality if all parameters are provided
  }

  #################
  # PREPROCESSING #
  #################

  # If specified, prune outgroup terminals
  if (!is.null(outgroup)) {
    tree1 <- drop.tip(tree1, outgroup)
    tree2 <- drop.tip(tree2, outgroup)
  }

  # Prune terminals not shared by both trees
  shared_terminals <- intersect(tree1$tip.label, tree2$tip.label)
  tree1_pruned <- drop.tip(tree1, tree1$tip.label[!(tree1$tip.label %in% shared_terminals)])
  tree2_pruned <- drop.tip(tree2, tree2$tip.label[!(tree2$tip.label %in% shared_terminals)])

  tree1_pruned = ladderize(tree1_pruned, right = TRUE) # Sort nodes in the tree according to clade size
  tree2_pruned = ladderize(tree2_pruned, right = TRUE) # Sort nodes in the tree according to clade size

  # If specified, reroot both trees using the same terminal
  if (!is.null(root)) {
    tree1_pruned <- root(tree1_pruned, outgroup = root)
    tree2_pruned <- root(tree2_pruned, outgroup = root)
  }

  #######################
  # DATASET PREPARATION #
  #######################

  # Create a dataframe with shared nodes of both trees
  s = RNODE::sharedNodes(tree1, tree2, spearman=F)
  # Create a list with unique nodes of each tree
  u = RNODE::uniqueNodes(tree1, tree2, composition = F)

  # Create vector with MOL support
  s1 = as.numeric(s$Support_Tree_1) # shared clades
  u1 = as.numeric(u[[1]]$Support) # unique clades
  # Create vector with TE occurrence
  s2 <- rep(1, length(s1))  # Create s2 with the same length as s1, filled with 1
  u2 <- rep(0, length(u1))  # Create u2 with the same length as u1, filled with 0

  # Combine the vectors into a data frame
  df <- data.frame(
    support_tree1 = c(s1, u1),        # Append u1 below s1
    occurrence_tree2 = c(s2, u2)     # Append u2 below s2
  )

  # Remove NA value from the root
  df = na.omit(df)

  # If requested, write a TSV file
  if (dataframe){write.table(df, "logistic.nodes.tsv", sep = "\t", row.names = FALSE)}

  return(df)
}
