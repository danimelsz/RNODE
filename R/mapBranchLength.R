#' @title mapBranchLength
#' @name mapBranchLength
#' @description Given a tree 1 without branch lengths (e.g. strict consensus of optimal trees) and a tree 2 with branch lengths (e.g. each MPT before reconciliation into a strict consensus), return tree 1 with branch lengths from shared clades with tree 2.
#' @author Daniel YM Nakamura
#'
#' @param tree1 A phylo object without branch lengths
#' @param tree2 A phylo object with branch lengths
#' @param write Optional. Specify the name of tree file to be written locally (nothing is written if this parameter is not specified)
#' @param outgroup Optional. Specify outgroup taxa to remove (by default, outgroup = F assumes that the user does not want to remove outgroup taxa)
#' @param root Optional. Specify the same root for both trees, which is recommended to facilitate tree comparisons (by default, root = F assumes that trees share the same root)
#' @examples
#' # Example 1 (identify unique nodes)
#' tree1 = read.tree (text="(t1,(t3,(t2,(t4,t5))));")
#' tree2 = read.tree (text="(t1,(t2,(t3,(t4,t5)47)53)94);")
#' mapBranchLength (tree1, tree2)
#'
#' @export
mapBranchLength = function (tree1,tree2,
                       write=NULL,
                       outgroup=NULL,
                       root=NULL,
                       plotTrees=F, node.numbers=T, tree.width=10, tree.height=10, tree.fsize=0.5, tree.adj=c(-1.5,0.5), tree.cex=2, tree.output="trees_unique_nodes.pdf"
){
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
  
  # If specified, reroot both trees using the same terminal
  if (!is.null(root)) {
    tree1_pruned <- root(tree1_pruned, outgroup = root)
    tree2_pruned <- root(tree2_pruned, outgroup = root)
  }
  
  # Ladderize trees
  tree1_pruned = ladderize(tree1_pruned, right = TRUE) # Sort nodes in the tree according to clade size
  tree2_pruned = ladderize(tree2_pruned, right = TRUE) # Sort nodes in the tree according to clade size
  
}