#' @title compareBranchLength
#' @name compareBranchLength
#' @description Given two trees with branch lengths, return a dataframe containing lengths from matching branches.
#' @author Daniel YM Nakamura
#'
#' @param tree1 A phylo object without branch lengths
#' @param tree2 A phylo object with branch lengths
#' @param write Optional. Specify the name of the dataframe file to be written locally (nothing is written if this parameter is not specified)
#' @param outgroup Optional. Specify outgroup taxa to remove (by default, outgroup = F assumes that the user does not want to remove outgroup taxa)
#' @param root Optional. Specify the same root for both trees, which is recommended to facilitate tree comparisons (by default, root = F assumes that trees share the same root)
#' @examples
#' # Example 1 (identify unique nodes)
#' tree1 = read.tree (text="(t1,(t3,(t2,(t4,t5))));")
#' tree2 = read.tree (text="(t1,(t2,(t3,(t4,t5)47)53)94);")
#' compareBranchLength (tree1, tree2)
#'
#' @export
mapBranchLength = function (tree1,tree2,
                            write=NULL,
                            outgroup=NULL,
                            root=NULL){
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
  
  #############
  # DATAFRAME #
  #############
  
  get_clades <- function(tree) {
    n_tips <- length(tree$tip.label)
    clades <- vector("list", nrow(tree$edge))
    edge_type <- character(nrow(tree$edge))
    
    for (i in seq_len(nrow(tree$edge))) {
      desc_node <- tree$edge[i, 2]
      if (desc_node <= n_tips) {
        # Descendant is a tip -> terminal edge
        edge_type[i] <- "terminal"
        desc <- desc_node
      } else {
        # Internal edge -> clade of descendant tips
        edge_type[i] <- "internal"
        desc <- Descendants(tree, desc_node, "tips")[[1]]
      }
      clades[[i]] <- sort(tree$tip.label[desc])
    }
    
    data.frame(
      ParentNode = tree$edge[, 1],
      ChildNode = tree$edge[, 2],
      EdgeType = edge_type,
      Clade = sapply(clades, paste, collapse = ","),
      EdgeLength = tree$edge.length,
      stringsAsFactors = FALSE
    )
  }
  
  
  # Extract clade info for both trees
  df1 <- get_clades(tree1_pruned)
  df2 <- get_clades(tree2_pruned)
  
  # Merge by matching clade composition
  df_merged <- merge(df1, df2, by = "Clade", all = TRUE, suffixes = c("_tree1", "_tree2"))
  
  return(df_merged)
  }