#' @title splitOrdFromUnord
#' @name splitOrdFromUnord
#' @description Splits a morphological matrix containing ordered and unordered characters into two matrices (ordered and unordered). This is useful to run phylogenetic analyses with MK and ORDERED models in IQTREE, especially if a concatenated matrix containing invariant characters is given as input.
#' @author Daniel YM Nakamura
#' @param input Input file (concatenated morphological matrix in Nexus format).
#' @param output_index Output index (e.g. if the user specify it as "Desktop/Index", the output files will be "Desktop/Index_ORDERED.nexus" and "Desktop/Index_UNORDERED.nexus")
#' @param list_ordered List of ordered characters e.g. c(1, 3, 9, 13). Character numbering starts with 1 (even if input data is .TNT).
#' @examples
#' # Example
#' splitOrdFromUnord
#'
#' @export
splitOrdFromUnord = function(input,
                             output_index,
                             list_ordered) {
  # Input
  morphANDunordered = TreeTools::ReadCharacters(input)

  # Create list of unordered characters
  list_unordered = setdiff(1:ncol(morphANDunordered), list_ordered)

  # Ordered matrix
  ordered = morphANDunordered[, list_ordered]

  # Output file name
  ordered_name = paste0(output_index, "_ORDERED.nexus")
  # Write a temporary file
  write.nexus.data(ordered, file=ordered_name, format="standard", interleaved=F)
  # Remove the interleave section to avoid errors in IQTREE
  ordered_temp = gsub("INTERLEAVE=NO", "", readLines(ordered_name))
  # Write a permanent file
  writeLines(ordered_temp, ordered_name)

  # Unordered matrix
  unordered = morphANDunordered[, list_unordered]

  # Output file name
  unordered_name = paste0(output_index, "_UNORDERED.nexus")
  # Write a temporary file
  write.nexus.data(unordered, file=unordered_name, format="standard", interleaved=F)
  # Remove the interleave section to avoid errors in IQTREE
  unordered_temp = gsub("INTERLEAVE=NO", "", readLines(unordered_name))
  # Write a permanent file
  writeLines(unordered_temp, unordered_name)
}
