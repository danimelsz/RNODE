#' @title mapSupport
#' @name mapSupport
#' @description Given a tree 1 without support values (e.g. strict consensus of optimal trees) and tree 2 with support values (e.g. majority consensus from suboptimal bootstrap pseudo-replicates), return tree 1 with support values from shared clades with tree 2.
#' @author Daniel YM Nakamura
#'
#' @param tree1 A phylo object without support values
#' @param tree2 A phylo object with support values
#' @param write Optional. Specify the name of tree file to be written locally (nothing is written if this parameter is not specified)
#' @param outgroup Optional. Specify outgroup taxa to remove (by default, outgroup = F assumes that the user does not want to remove outgroup taxa)
#' @param root Optional. Specify the same root for both trees, which is recommended to facilitate tree comparisons (by default, root = F assumes that trees share the same root)
#' @param plotTrees Optional. Plot the two trees (tree1 with support values on left and tree2 on right) in \code{PDF} format. If \code{plot = T}, the user should also adjust \code{PDF} dimensions (e.g. \code{width = 8}, \code{height = 8}), label size (e.g. \code{fsize = 4}), and position and size of support values (e.g. \code{adj = c(-1.5,0.5)}, \code{cex = 0.6}).
#' @param tree.output Optional. Name of the output figure.
#' @param tree.width Optional. Width of trees in PDF if plotTrees = T.
#' @param tree.height Optional. Height of trees in PDF if plotTrees = T.
#' @param tree.fsize Optional. Font size in PDF if plotTrees = T.
#' @param tree.adj Optional. Adjust horizontal and vertical position if plotTrees = T.
#' @param tree.cex Optional. Adjust support size in nodes if plotTrees = T.
#' @param node.numbers Optional. If plotTrees = T, show node index (do not confuse with support values'by default, True).
#' @examples
#' # Example 1 (identify unique nodes)
#' tree1 = read.tree (text="(t1,(t3,(t2,(t4,t5))));")
#' tree2 = read.tree (text="(t1,(t2,(t3,(t4,t5)47)53)94);")
#' mapSupport (tree1, tree2)
#'
#' @export
mapSupport = function (tree1,tree2,
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

  #######################################
  # PLOTTING HIGHLIGHTING UNIQUE CLADES #
  #######################################

  # If specified, plot trees with node index (inside squares) and support values
  if (plotTrees) {
    pdf(file=tree.output, width = tree.width, height = tree.height)  # Save plotted tree to PDF, adjust width and height as needed
    par(mfrow = c(1, 2), oma=c(1,0.5,1,0.5))
    # Tree 1
    plotTree(tree1_pruned, fsize = tree.fsize, ftype="i", node.numbers=node.numbers, color="black") # Adjust font size as needed
    df_unique = uniqueNodes(tree1_pruned, tree2_pruned)
    nodelabels(node=df_unique[[1]]$Node,
               cex=tree.cex, # Adjust the size of circles
               pch=21, bg="blue")
    nodelabels(tree1_pruned$node.label,
               cex = tree.cex,       # font size
               frame = "none",  # no box
               adj = tree.adj) # adjust position
    plotTree(tree2_pruned, fsize = tree.fsize, ftype="i", node.numbers=node.numbers, color="black", direction="leftwards") # Adjust font size as needed
    nodelabels(node=df_unique[[2]]$Node,
               cex=tree.cex, # Adjust the size of circles
               pch=21, bg="red")
    nodelabels(tree2_pruned$node.label,
               cex = tree.cex,       # font size
               frame = "none",  # no box
               adj = tree.adj) # adjust position
    dev.off()
  }

  return (list(tree1_pruned, m))
}
