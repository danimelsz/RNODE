#' @title uniqueNodes
#' @name uniqueNodes
#' @description Creates two separate dataframes containing the list of descendants and support values of unique clades between two trees.
#' @author Daniel YM Nakamura, Taran Grant
#'
#' @param tree1 A .phylo tree that can be loaded using ape::read.tree for NEWICK files or TreeTools::ReadTntTree for TNT files
#' @param tree2 Another .phylo tree
#' @param composition Optional. Specify if composition of corresponding clades should be present in the dataframe (by default, composition = T)
#' @param outgroup Optional. Specify outgroup taxa to remove (by default, outgroup = F assumes that the user does not want to remove outgroup taxa)
#' @param root Optional. Specify the same root for both trees, which is recommended to facilitate tree comparisons (by default, root = F assumes that trees share the same root)
#' @param dataframe Optional. Write a TSV file in current directory containing the output dataframe (by default, dataset = T).
#'
#' @examples
#' # Example 1 (identify unique nodes)
#' tree1 = read.tree (text="(t1,(t2,(t3,(t4,t5)75)32)45);")
#' tree2 = read.tree (text="(t1,(t6,(t3,(t4,t5)47)53)94);")
#' uniqueNodes (tree1, tree2)
#'
#' # Example 2 (count unique nodes)
#' tree1 = read.tree (text="(t1,(t2,(t3,(t4,t5)75)32)45);")
#' tree2 = read.tree (text="(t1,(t6,(t3,(t4,t5)47)53)94);")
#' nrow (uniqueNodes (tree1, tree2, composition = F))
#'
#' @export
uniqueNodes = function(tree1, tree2,
                       composition=T,
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

  # Check if the input trees contains $node.label (support)
  if (is.null(tree1$node.label) && is.null(tree2$node.label)) {
    stop("Input trees with no support values stored at node.labels")
  } else if (is.null(tree1$node.label) && is.null(tree2$node.label)==F){
    stop("Input tree 1 with no support values stored at node.labels")
  } else if (is.null(tree1$node.label)==F && is.null(tree2$node.label)){
    stop("Input tree 2 with no support values stored at node.labels")
  } else if (is.null(tree1$node.label)==F && is.null(tree2$node.label)==F) {
    print("Both trees with support values.")
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

  ################################################
  #                                              #
  #    Clades unique to tree1 (after pruning)    #
  #                                              #
  ################################################

  # Select clades unique to tree1_pruned
  m <- matchNodes(tree1_pruned, tree2_pruned, method="descendants")
  m_NA <- which(is.na(m[, 2]))

  # Extract the support values for the unique nodes in tree1_pruned
  m_NA_support_values <- tree1_pruned$node.label[m_NA]
  m_NA_support_values <- gsub("=", "", m_NA_support_values) # Clean the support values (remove "=" symbol if present)
  print(m_NA_support_values)

  # Looping to store descendants from unique clades in tree1
  c1_NA = c() # empty vector
  for (i in m_NA + length(tree1_pruned$tip.label)) {
    c1_NA.temp = extract.clade(tree1_pruned, i)
    c1_NA.temp.nwk = writeTree(c1_NA.temp)
    c1_NA.temp.nwk.noBranchLength = str_replace_all(c1_NA.temp.nwk, ":\\d+(\\.\\d+)?", ":") # clean branch lengths
    c1_NA.temp.nwk.noColon = gsub(":", "", c1_NA.temp.nwk.noBranchLength) # clean colons :
    c1_NA = c(c1_NA, c1_NA.temp.nwk.noColon)
  }

  # Dataframe with unique clades in tree 1
  if (composition==T && !is.null(tree1_pruned$node.label)){
    df_unique_Tree1 <- data.frame(Node = m_NA + length(tree1_pruned$tip.label),
                                  Support = m_NA_support_values,
                                  Descendants = c1_NA)}
  else if (composition==F && !is.null(tree1_pruned$node.label)){
    df_unique_Tree1 <- data.frame(Node = m_NA + length(tree1_pruned$tip.label),
                                  Support = m_NA_support_values)}
  else if (composition==F && is.null(tree1_pruned$node.label)){
    df_unique_Tree1 <- data.frame(Node = m_NA + length(tree1_pruned$tip.label))
  }
  if (dataframe){write.table(df_unique_Tree1, "Tree1_unique.clades.tsv", sep = "\t", row.names = FALSE)}

  ################################################
  #                                              #
  #   Clades unique to tree2 (after pruning)     #
  #                                              #
  ################################################

  # Select clades unique to tree2_pruned
  n <- matchNodes(tree2_pruned, tree1_pruned, method="descendants")
  n_NA <- which(is.na(n[, 2]))

  # Extract the support values for the unique nodes in tree2_pruned
  n_NA_support_values <- tree2_pruned$node.label[n_NA]
  n_NA_support_values <- gsub("=", "", n_NA_support_values) # Clean the support values (remove "=" symbol if present)
  print(n_NA_support_values)

  # Looping to store descendants from unique clades in tree1
  c2_NA = c() # empty vector
  for (i in n_NA + length(tree2_pruned$tip.label)) {
    c2_NA.temp = extract.clade(tree2_pruned, i)
    c2_NA.temp.nwk = writeTree(c2_NA.temp)
    c2_NA.temp.nwk.noBranchLength = str_replace_all(c2_NA.temp.nwk, ":\\d+(\\.\\d+)?", ":") # clean branch lengths
    c2_NA.temp.nwk.noColon = gsub(":", "", c2_NA.temp.nwk.noBranchLength) # clean colons :
    c2_NA = c(c2_NA, c2_NA.temp.nwk.noColon)
  }

  # Dataframe with unique clades in tree 2
  if (composition==T && !is.null(tree2_pruned$node.label)){
    df_unique_Tree2 <- data.frame(Node = n_NA + length(tree2_pruned$tip.label),
                                  Support = n_NA_support_values,
                                  Descendants = c2_NA)}
  else if (composition==F && !is.null(tree2_pruned$node.label)){
    df_unique_Tree2 <- data.frame(Node = n_NA + length(tree2_pruned$tip.label),
                                  Support = n_NA_support_values)}
  else if (composition==F && is.null(tree2_pruned$node.label)){
    df_unique_Tree2 <- data.frame(Node = n_NA + length(tree2_pruned$tip.label))
  }

  if (dataframe){write.table(df_unique_Tree2, "Tree2_unique.clades.tsv", sep = "\t", row.names = FALSE)}

  return(list(df_unique_Tree1, df_unique_Tree2))
}
