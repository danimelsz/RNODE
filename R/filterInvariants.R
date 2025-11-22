#' @title filterInvariants
#' @name filterInvariants
#' @description Delete invariant characters in a morphological matrix. We follow the definition of invariant site from IQ-Tree: (1) constant sites containing only a single character state in all sequences, (2) partially constant sites (N and/or -), and (3) ambiguously constant sites (e.g. C, Y and -).
#' @author Daniel YM Nakamura
#' @param input Input file (molecular or morphological matrix in Nexus or already loaded as a matrix object in R).
#' @param output_index Output index (e.g. if the user specify it as "Desktop/Index", the output files will be "Desktop/Index_onlyVARIANTS.nexus")
#' @examples
#' # Example
#' filterInvariants(input="../testdata/015_MORPH_data.nexus", output_index="../testdata/015_MORPH_data")
#'
#' @export
filterInvariants <- function(input, output_index) {

  # --- Load input ---
  if (is.character(input) && file.exists(input)) {
    mat <- TreeTools::ReadCharacters(input)
  } else {
    mat <- input
    if (!is.matrix(mat) && !is.data.frame(mat)) {
      stop("Input must be either a matrix/data.frame or a valid file path.")
    }
  }

  mat <- as.matrix(mat)

  # --- Helper: parse states ---
  parse_state <- function(x) {
    x <- toupper(x)

    # Nucleotide ambiguity codes
    nuc_map <- list(
      A="A", C="C", G="G", T="T", U="T",
      R=c("A","G"), Y=c("C","T"), S=c("G","C"), W=c("A","T"),
      K=c("G","T"), M=c("A","C"),
      B=c("C","G","T"), D=c("A","G","T"),
      H=c("A","C","T"), V=c("A","C","G"),
      N=c("A","C","G","T")
    )

    # Missing data
    if (x %in% c("?", "-")) return(NA)

    # Nucleotide?
    if (x %in% names(nuc_map)) return(nuc_map[[x]])

    # Morphological polymorphisms [01] {01} (01)
    if (grepl("^\\[[0-9]+\\]$", x) ||
        grepl("^\\{[0-9]+\\}$", x) ||
        grepl("^\\([0-9]+\\)$", x)) {
      inner <- gsub("\\[|\\]|\\{|\\}|\\(|\\)", "", x)
      return(strsplit(inner, "")[[1]])
    }

    # Single morphological digit
    if (grepl("^[0-9]$", x)) return(x)

    stop(paste("Unrecognized symbol:", x))
  }

  morph_states <- as.character(0:9)

  # Determine possible states in a column
  get_possible_states <- function(chars) {
    sets <- lapply(chars, parse_state)

    # Missing data (NA) expands to all morph states
    sets <- lapply(sets, function(s) {
      if (length(s) == 1 && is.na(s)) morph_states else s
    })

    Reduce(intersect, sets)
  }

  # Determine variable vs invariant columns
  is_variable <- sapply(seq_len(ncol(mat)), function(i) {
    col <- mat[, i]
    col_no_na <- col[!is.na(col)]

    possible <- get_possible_states(col_no_na)

    # invariant only if intersection == 1 state
    !(length(possible) == 1)
  })

  deleted_columns <- which(!is_variable)

  # --- Report ---
  cat("Deleted", length(deleted_columns), "invariant columns.\n")
  if (length(deleted_columns) > 0)
    cat("Columns deleted:", paste(deleted_columns, collapse = ", "), "\n")

  # Filter to variable columns
  filtered <- mat[, is_variable, drop = FALSE]

  # --- Output file ---
  output_name <- paste0(output_index, "_onlyVARIANTS.nexus")

  # Prepare data for TreeTools output
  df <- as.data.frame(filtered)
  write.nexus.data(df, file = output_name, format = "standard", interleaved = FALSE)

  # Fix IQ-TREE incompatibility
  lines <- readLines(output_name)
  lines <- gsub("INTERLEAVE=NO", "", lines)
  writeLines(lines, output_name)

  return(filtered)
}
