#' @title sharedNodes
#' @name sharedNodes
#' @description Compare support values of shared clades between two trees. The outputs are (1) basic statistics about number of shared clades and support values; (2) a dataframe with node labels, descendants, and support values of shared clades, which facilitates descriptive and statistical comparisons of clade composition and support between corresponding nodes.
#' @author Daniel YM Nakamura, Taran Grant
#'
#' @param tree1 A \code{phylo} tree that can be loaded using \code{ape::read.tree} for \code{NEWICK} files or \code{TreeTools::ReadTntTree} for \code{TNT} files. The \code{phylo} must contain the element \code{$node.label}.
#' @param tree2 Another \code{phylo} tree
#' @param composition Optional. Specify if composition of corresponding clades should be present in the dataframe (by default, \code{composition = F})
#' @param outgroup Optional. Specify outgroup taxa to remove (by default, the function assumes that the user does not want to remove outgroup taxa; i.e. \code{outgroup = NULL})
#' @param root Optional. Specify the same root for both trees, which is recommended to facilitate tree comparisons (by default, the function assumes that trees share the same root; i.e. \code{root = NULL})
#' @param plotTrees Optional. Plot the two trees after taxa pruning in \code{PDF} format. If \code{plot = T}, the user should also adjust \code{PDF} dimensions (e.g. \code{width = 8}, \code{height = 8}), label size (e.g. \code{fsize = 4}), and position and size of support values (e.g. \code{adj = c(-1.5,0.5)}, \code{cex = 0.6}).
#' @param node.numbers Optional. If plotTrees = T, show node index (do not confuse with support values'by default, True).
#' @param tanglegram Optional. Plot a tanglegram minimizing the number of crosses of lines linking two trees in \code{PDF} format.
#' @param dataframe Optional. Write a \code{TSV} file in current directory containing the output dataframe (by default, no \code{TSV} is written).
#' @param spearman Optional. Test the correlation between support values using a Spearman test (by default, \code{spearman = T}).
#'
#' @examples
#' # Example 1 (simplest case)
#' tree1 = read.tree (text="(t1,(t2,(t3,(t4,t5)75)32)45);")
#' tree2 = read.tree (text="(t1,(t6,(t3,(t4,t5)47)53)94);")
#' sharedNodes (tree1, tree2)
#'
#' # Example 2 (show internal topology of each node, remove outgroup taxa t9 and t8, and reroot in t1 in both trees)
#' sharedNodes (tree1, tree2, composition=T, outgroup=c("t9", "t8"), root="t1",)
#'
#' # Example 3 (plot  two trees)
#' sharedNodes (tree1, tree2, plotTrees=T, width=8, height=8, fsize=3, adj=c(-1.5,0.5), cex=0.6)
#'
#' @export
sharedNodes = function (tree1,tree2,
                        composition=F,
                        outgroup=NULL,
                        root=NULL,
                        plotTrees=F, node.numbers=T, width=NULL, height=NULL, fsize=NULL, adj=NULL, cex=NULL,
                        tanglegram=F,
                        dataframe=F, messages=T,
                        spearman=F){
  # Initial warnings
  missing_params <- c()
  if (is.null(tree1)) missing_params <- c(missing_params, "tree1")
  if (is.null(tree2)) missing_params <- c(missing_params, "tree2")
  if (length(missing_params) > 0) {   # Check if there are any missing parameters and print a message
    message("The following parameters are missing: ", paste(missing_params, collapse = ", "))
  } else {
    message("All required parameters provided.") # Proceed with the main functionality if all parameters are provided
  }

  '
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
  '

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

  # If specified, plot trees with node index (inside squares) and support values
  tree1_pruned = ladderize(tree1_pruned, right = TRUE) # Sort nodes in the tree according to clade size
  tree2_pruned = ladderize(tree2_pruned, right = TRUE) # Sort nodes in the tree according to clade size
  if (plotTrees) {
    pdf("tree1_pruned.pdf", width = width, height = height)  # Save plotted tree to PDF, adjust width and height as needed
    plotTree(tree1_pruned, fsize = fsize, ftype="i", node.numbers=node.numbers, color="blue") # Adjust font size as needed
    nodelabels(tree1_pruned$node.label,
               adj=adj, # Adjust horizontal and vertical position
               frame="none", # Specify the borders of support values
               cex=cex) # Adjust the size of support values
    dev.off()
    pdf("tree2_pruned.pdf", width = width, height = height)  # Save plotted tree to PDF, adjust width and height as needed
    plotTree(tree2_pruned, fsize = fsize, ftype="i", node.numbers=node.numbers, color="red") # Adjust font size as needed
    nodelabels(tree2_pruned$node.label,
               adj=adj, # Adjust horizontal and vertical position
               frame="none", # Specify borders of support values
               cex=cex) # Adjust support size
    dev.off()
  }

  # If specified, reroot both trees using the same terminal
  if (!is.null(root)) {
    tree1_pruned <- root(tree1_pruned, outgroup = root)
    tree2_pruned <- root(tree2_pruned, outgroup = root)
  }

  ##########################
  # ASSESSING SHARED NODES #
  ##########################

  # Generate a matrix of clades in tree1 present and absent in tree2
  m = matchNodes(tree1_pruned, tree2_pruned, method="descendants")
  # Show only shared clades in the matrix
  m_noNA <- m[complete.cases(m), ]
  # Vector with only shared clades in pruned tree 1
  m1 = m_noNA[,1] - length(tree1_pruned$tip.label)
  # Vector with only shared clades in pruned tree 2
  m2 = m_noNA[,2] - length(tree2_pruned$tip.label)

  # Create a vector of support in shared clades (values from pruned tree 1)
  s1 = c() # empty vector
    for (x in m1) {
      s1.temp = tree1_pruned$node.label[x]
      s1.temp <- gsub("=", "", s1.temp) # Clean the symbol "=" before support values from TNT output
      s1 = c(s1, s1.temp)
    }

  # Create a vector of support in shared clades (values from pruned tree 2)
  s2 = c() # empty vector
    for (y in m2) {
      s2.temp = tree2_pruned$node.label[y]
      s2.temp <- gsub("=", "", s2.temp)  # Clean the symbol "=" before support values from TNT output
      s2 = c(s2, s2.temp)
    }

  # Create a vector of terminal composition in each shared node (clades from pruned tree 1)
  c1 = c() # empty vector
  for (i in m_noNA[,1]) {
    i
    c1.temp = extract.clade(tree1_pruned, i)
    c1.temp.nwk = writeTree(c1.temp)
    c1.temp.nwk.noBranchLength = str_replace_all(c1.temp.nwk, ":\\d+(\\.\\d+)?", ":") # clean branch lengths
    c1.temp.nwk.noColon = gsub(":", "", c1.temp.nwk.noBranchLength) # clean colons :
    c1 = c(c1, c1.temp.nwk.noColon)
  }

  # Create a vector of terminal composition in each shared node (clades from pruned tree 2)
  c2 = c() # empty vector
  for (i in m_noNA[,2]) {
    i
    c2.temp = extract.clade(tree2_pruned, i)
    c2.temp.nwk = writeTree(c2.temp)
    c2.temp.nwk.noBranchLength = str_replace_all(c2.temp.nwk, ":\\d+(\\.\\d+)?", ":") # clean branch lengths
    c2.temp.nwk.noColon = gsub(":", "", c2.temp.nwk.noBranchLength) # clean colons :
    c2 = c(c2, c2.temp.nwk.noColon)
  }

  # Create dataframe
  if (composition==TRUE && !is.null(tree1_pruned$node.label)) { # if composition should be given and support values are present in the tree 1
    df <- data.frame(Node_Tree_1 = m_noNA[,1],
                   Node_Tree_2 = m_noNA[,2],
                   Support_Tree_1 = s1,
                   Support_Tree_2 = s2,
                   Composition_Tree_1 = c1,
                   Composition_Tree_2 = c2)
  } else if (composition==F && !is.null(tree1_pruned$node.label)) { # if composition should be ignored and support values are present in the tree 1
    df <- data.frame(Node_Tree_1 = m_noNA[,1],
                     Node_Tree_2 = m_noNA[,2],
                     Support_Tree_1 = s1,
                     Support_Tree_2 = s2)
  } else if (composition==F && is.null(tree1_pruned$node.label)) { # if composition should be ignored and support values are present in the tree 1
    df <- data.frame(Node_Tree_1 = m_noNA[,1],
                     Node_Tree_2 = m_noNA[,2])
  }

  # Basic statistics
  nClades1 = tree1_pruned$Nnode
  nClades2 = tree2_pruned$Nnode
  meanSupport1 = mean(na.omit(as.numeric(tree1_pruned$node.label)))
  sup1 = as.numeric(df$Support_Tree_1[df$Support_Tree_1 != ""]) # delete empty value of the root
  meanSupport2 = mean(na.omit(as.numeric(tree2_pruned$node.label)))
  sup2 = as.numeric(df$Support_Tree_2[df$Support_Tree_2 != ""]) # delete empty value of the root
  if (spearman) {spearman = cor.test (sup1, sup2, method="spearman")}

  # Print
  if (messages) {
  print ("")
  print ("Tree comparisons done!")
  print (paste("Number of shared clades: ", nrow(df)))
  print (paste("Tree 1: Total number of clades =", nClades1, "; Mean support =", meanSupport1))
  print (paste("Support of shared clades in tree 1: ", min(sup1), "–", max(sup1), " (", round(mean(sup1),2), ")", sep=""))
  print (paste("Tree 2: Total number of clades =", nClades2, "; Mean support =", meanSupport2))
  print (paste("Support of shared clades in tree 2: ", min(sup2), "–", max(sup2), " (", round(mean(sup2),2), ")", sep=""))
  if (spearman) {print (paste("Spearman's test: RHO = ", spearman$estimate, "; P-value = ", spearman$p.value))}
  }

  # Output
  if (dataframe) {
    write.table(df, "shared.clades.tsv", sep = "\t", row.names = FALSE)
  }
  return (df)

  ############################
  # TANGLEGRAM VISUALIZATION #
  ############################

  if (tanglegram==T) {

  # Force the trees to be ultrametric
  tree1_ultrametric <- force.ultrametric(tree1_pruned, method = "extend")
  tree2_ultrametric <- force.ultrametric(tree2_pruned, method = "extend")

  # Resolve polytomies in the trees; multi2di() resolves polytomies by adding zero-length branches
  tree1_binary <- multi2di(tree1_ultrametric)
  tree2_binary <- multi2di(tree2_ultrametric)

  # Convert binary ultrametric trees to hclust objects
  hclust1 <- ape::as.hclust.phylo(tree1_binary)
  hclust2 <- ape::as.hclust.phylo(tree2_binary)

  # Convert hclust objects to dendrograms
  dend1 <- as.dendrogram(hclust1)
  dend2 <- as.dendrogram(hclust2)

  # Save tanglegram as a PDF
  pdf("tanglegram_comparison.pdf", width = 5, height = 5) # Specify dimensions if needed
  tanglegram(a, b,
             highlight_distinct_edges = TRUE,   # Highlights differences
             common_subtrees_color_lines = TRUE, # Colors common subtrees
             edge.lwd = c(1, 1),  # Adjust edge thickness (left, right)
             lab.cex = 0.5,  # Adjust font size of terminal names
             main_left = "Tree 1", # Title Tree 1
             main_right = "Tree 2", # Title Tree 2
             margin_inner = 1) # Distance between trees

  dev.off()  # Close the PDF device
  }

}
