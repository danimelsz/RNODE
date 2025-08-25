#' @title mapSupport
#' @name pasteSupport
#' @description Given a tree 1 without support values (e.g. strict consensus of optimal trees) and tree 2 with support values (e.g. majority consensus from suboptimal bootstrap pseudo-replicates), return tree 1 with support values from shared clades with tree 2.
#' @author Daniel YM Nakamura
#'
#' @param tree1 A phylo object without support values
#' @param tree2 A phylo object with support values
#' @param write Optional. Specify the name of tree file to be written locally (nothing is written if this parameter is not specified)
#' @param outgroup Optional. Specify outgroup taxa to remove (by default, outgroup = F assumes that the user does not want to remove outgroup taxa)
#' @param root Optional. Specify the same root for both trees, which is recommended to facilitate tree comparisons (by default, root = F assumes that trees share the same root)
#' @examples
#' # Example 1 (identify unique nodes)
#' tree1 = read.tree (text="(t1,(t3,(t2,(t4,t5))));")
#' tree2 = read.tree (text="(t1,(t2,(t3,(t4,t5)47)53)94);")
#' pasteSupport (tree1, tree2)
#'
#' @export
pasteSupport = function (tree1,tree2,
                         write=NULL,
                         outgroup=NULL,
                         root=NULL
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

  ########################
  # PASTE SUPPORT VALUES #
  ########################

  # Generate a matrix of clades in tree1 (no support e.g. consensus) present and absent in tree2 (with support)
  m = matchNodes(tree1_pruned, tree2_pruned, method="descendants")
  # Ensure first two columns are tr1/tr2 (adjust if already named)
  colnames(m)[1:2] <- c("tr1", "tr2")
  # Offset: node.labels are for internal nodes  (Ntip+1 ... Ntip+Nnode)
  ntip2 <- ape::Ntip(tree2_pruned)
  # Clean labels to numeric (handles "100%" or "")
  lab2 <- as.numeric(gsub("%", "", tree2_pruned$node.label))
  # Index of each tr2 node into node.label
  idx <- as.numeric(m[, "tr2"]) - ntip2
  # Build support with NA where tr2 is NA or out of range
  tr2_support <- rep(NA_real_, nrow(m))
  ok <- !is.na(idx) & idx >= 1 & idx <= length(lab2)
  tr2_support[ok] <- lab2[idx[ok]]
  # Bind as new column
  m <- cbind(m, tr2_support = tr2_support)
  # Paste support values from tree 2 to tree 1
  tree1_pruned$node.label <- m[, "tr2_support"]

  # Write tree 1 (originally lacking support values, now with support from tree 2)
  if (!is.null(write)) {
    file_name = paste0(write, ".nwk")
    write.tree(ladderize(tree1_pruned, right=T), file=file_name)
  }

  return (list(tree1_pruned, m))
}
