#' @title mapBranchLength
#' @name mapBranchLength
#' @description Given a tree 1 without branch lengths (e.g. strict consensus of optimal trees) and a tree 2 with branch lengths (e.g. each MPT before reconciliation into a strict consensus), return tree 1 with branch lengths from shared clades with tree 2.
#' @author Daniel YM Nakamura
#'
#' @param tree1 A phylo object without branch lengths
#' @param trees2 A phylo object with branch lengths
#' @param method. Optional. Method to sample branch lengths. If "random", randomly select one of the trees from trees2. If "minimum" (default), map the minimum values from the pool of MPTs (trees2) to each edge present in tree1
#' @examples
#' # Example 1 (identify unique nodes)
#' tree1 = read.tree (text="(t1,(t3,(t2,(t4,t5))));")
#' tree2 = read.tree (text="(t1,(t2,(t3,(t4,t5)47)53)94);")
#' mapBranchLength (tree1, tree2)
#'
#' @export
mapBranchLength <- function(tree1, trees2, method = "minimum") {
  # Dependencies
  if (!requireNamespace("ape", quietly = TRUE))
    stop("Package 'ape' is required.")

  method <- match.arg(method)

  # Ensure trees2 is a list
  if (inherits(trees2, "phylo")) trees2 <- list(trees2)

  # Check taxon sets
  tips1 <- tree1$tip.label
  for (tr in trees2) {
    if (!setequal(tr$tip.label, tips1))
      stop("All trees must have the same tip set as tree1.")
  }

  # Helper to get all clades from a tree
  get_clades <- function(tree) {
    n_tips <- length(tree$tip.label)
    lapply((n_tips + 1):(n_tips + tree$Nnode), function(node)
      sort(ape::extract.clade(tree, node)$tip.label))
  }

  # Store all clades
  clades1 <- get_clades(tree1)
  clades2_list <- lapply(trees2, get_clades)

  # Prepare result
  n_edges <- nrow(tree1$edge)
  mapped_lengths <- numeric(n_edges)

  # --- Step 1: Map terminal (leaf) edges ---
  for (i in seq_along(tips1)) {
    tip <- tips1[i]
    edge_idx <- which(tree1$edge[, 2] == i)
    lengths <- numeric()

    for (tr in trees2) {
      idx2 <- which(tr$edge[, 2] == which(tr$tip.label == tip))
      lengths <- c(lengths, tr$edge.length[idx2])
    }

    if (length(lengths) == 0) {
      mapped_lengths[edge_idx] <- 0  # Should not happen normally
    } else if (method == "random") {
      mapped_lengths[edge_idx] <- sample(lengths, 1)
    } else if (method == "minimum") {
      mapped_lengths[edge_idx] <- min(lengths)
    }
  }

  # --- Step 2: Map internal edges ---
  n_tips <- length(tree1$tip.label)
  for (i in seq_along(clades1)) {
    clade1 <- clades1[[i]]
    edge_idx <- which(tree1$edge[, 2] == (n_tips + i))
    lengths <- numeric()

    # Collect branch lengths from trees2 for shared clades
    for (j in seq_along(trees2)) {
      clades2 <- clades2_list[[j]]
      match_idx <- which(sapply(clades2, function(cl) identical(cl, clade1)))
      if (length(match_idx) > 0) {
        edge2 <- which(trees2[[j]]$edge[, 2] == (length(tips1) + match_idx))
        lengths <- c(lengths, trees2[[j]]$edge.length[edge2])
      }
    }

    # Assign branch length or 0 if no match (collapsed edge)
    if (length(lengths) == 0) {
      mapped_lengths[edge_idx] <- 0
    } else if (method == "random") {
      mapped_lengths[edge_idx] <- sample(lengths, 1)
    } else if (method == "minimum") {
      mapped_lengths[edge_idx] <- min(lengths)
    }
  }

  # Attach lengths
  tree1$edge.length <- mapped_lengths
  return(tree1)
}


