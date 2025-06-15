#' @title filterMissing
#' @name filterMissing
#' @description Filter taxa (rows) and/or characters (columns) in a matrix containing only missing data.
#' @author Daniel YM Nakamura
#' @param input Input file (morphological matrix in Nexus or TNT format).
#' @param output_path Output path (e.g. if the user specify it as "Desktop/Index", the output files will be "Desktop/Index_ORDERED.nexus" and "Desktop/Index_UNORDERED.nexus")
#' @param missing Parameter specifying if rows and/or columns in which all cells are missing data (?) should be removed. Options: "row" (default), "column", "both".
#' @examples
#' # Example
#' filterMissing (input="testdata/test_filterMissing.nexus", output_path="testdata/test", missing="row")
#'
#' @export
filterMissing = function(input,
                       output_path,
                       missing="row") {

  # Input
  data = TreeTools::ReadCharacters(input)

  # If missing is 'row', delete rows containing only ?
  if (missing == 'row' || missing == "both") {
    rows_to_delete <- apply(data, 1, function(row) all(row == "?"))
    cat("Rows deleted:", which(rows_to_delete), "\n")
    data <- data[!rows_to_delete, ]
  }

  # If missing is 'row', delete rows containing only ?
  if (missing == 'column' || missing == "both") {
    cols_to_delete <- apply(data, 2, function(col) all(col == "?"))
    cat("Columns deleted:", which(cols_to_delete), "\n")
    data <- data[!cols_to_delete, ]
  }

  # Output file name
  name = paste0(output_path, "_FILTERED.nexus")
  # Write a temporary file
  write.nexus.data(data, file=name, format="standard", interleaved=F)
  # Remove the interleave section to avoid errors in IQTREE
  temp = gsub("INTERLEAVE=NO", "", readLines(name))
  # Write a permanent file
  writeLines(temp, name)
}
