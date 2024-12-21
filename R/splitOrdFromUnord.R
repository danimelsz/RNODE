#' @title splitOrdFromUnord
#' @name splitOrdFromUnord
#' @description Splits a morphological matrix containing ordered and unordered characters into two matrices (ordered and unordered). This is useful to run phylogenetic analyses with MK and ORDERED models in IQTREE, especially if a concatenated matrix containing invariant characters is given as input.
#' @author Daniel YM Nakamura
#' @param input Input file (concatenated morphological matrix in Nexus or TNT format).
#' @param output_index Output index (e.g. if the user specify it as "Desktop/Index", the output files will be "Desktop/Index_ORDERED.nexus" and "Desktop/Index_UNORDERED.nexus")
#' @param list_ordered List of ordered characters e.g. c(1, 3, 9, 13). Character numbering starts with 1 (even if input data is .TNT).
#' @param invariant Logical value specifying if invariant characters should be removed. Following IQTREE, invariants are columns with (1) a unique state (i.e. all cells with the same value) or (2) a unique state or ? or - (i.e. some cells are the same number, other cells are ?, and/or other cells are -)
#' @examples
#' # Example
#' splitOrdFromUnord
#'
#' @export
splitOrdFromUnord = function(input,
                             output_index,
                             list_ordered,
                             invariant = F) {
  # Input
  morphANDunordered = TreeTools::ReadCharacters(input)

  # Create list of unordered characters
  list_unordered = setdiff(1:ncol(morphANDunordered), list_ordered)

  # Ordered matrix
  ordered = morphANDunordered[, list_ordered]
  # Remove invariant characters
  if (invariant) {
    ordered_ <- ordered[, apply(ordered, 2, function(col) {
      unique_values <- unique(col)
      valid_values <- unique_values[!unique_values %in% c("?", "-")]
      !(length(valid_values) == 1 && all(unique_values %in% c("?", "-", valid_values)))
    })]
    print(paste("Number of removed characters in the ordered matrix: ", ncol(ordered)-ncol(ordered_)))
    ordered = ordered_
  }
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
  # Remove invariant characters
  if (invariant) {
    unordered_ <- unordered[, apply(unordered, 2, function(col) {
      unique_values <- unique(col)
      valid_values <- unique_values[!unique_values %in% c("?", "-")]
      !(length(valid_values) == 1 && all(unique_values %in% c("?", "-", valid_values)))
    })]
    print(paste("Number of removed characters in the unordered matrix: ", ncol(unordered)-ncol(unordered_)))
    unordered = unordered_
  }
  # Output file name
  unordered_name = paste0(output_index, "_UNORDERED.nexus")
  # Write a temporary file
  write.nexus.data(unordered, file=unordered_name, format="standard", interleaved=F)
  # Remove the interleave section to avoid errors in IQTREE
  unordered_temp = gsub("INTERLEAVE=NO", "", readLines(unordered_name))
  # Write a permanent file
  writeLines(unordered_temp, unordered_name)
}
