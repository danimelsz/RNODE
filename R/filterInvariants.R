#' @title filterInvariants
#' @name filterInvariants
#' @description Delete invariant characters in a morphological matrix. We follow the definition of invariant site from IQ-Tree: (1) constant sites containing only a single character state in all sequences, (2) partially constant sites (N and/or -), and (3) ambiguously constant sites (e.g. C, Y and -).
#' @author Daniel YM Nakamura
#' @param input Input file (molecular or morphological matrix in Nexus).
#' @param output_index Output index (e.g. if the user specify it as "Desktop/Index", the output files will be "Desktop/Index_onlyVARIANTS.nexus")
#' @examples
#' # Example
#' filterInvariants
#'
#' @export
filterInvariants = function(input,
                             output_index) {
  # Input
  mat = TreeTools::ReadCharacters(input)

  # --- Helper: parse morphological and nucleotide states ---
  parse_state <- function(x) {
    x <- toupper(x)

    # Nucleotide ambiguity mapping
    nuc_map <- list(
      A="A", C="C", G="G", T="T", U="T",
      R=c("A","G"), Y=c("C","T"), S=c("G","C"), W=c("A","T"),
      K=c("G","T"), M=c("A","C"), B=c("C","G","T"), D=c("A","G","T"),
      H=c("A","C","T"), V=c("A","C","G"), N=c("A","C","G","T")
    )

    # Missing data → full ambiguity
    if (x %in% c("?", "-")) return(NA)

    # Nucleotide?
    if (x %in% names(nuc_map)) return(nuc_map[[x]])

    # Morphological polymorphism formats: [01], {01}, (01)
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

  # Full morphological state space
  morph_states <- as.character(0:9)

  get_possible_states <- function(chars) {
    sets <- lapply(chars, parse_state)

    # Corrected: only treat NA as missing when length == 1
    sets <- lapply(sets, function(s) {
      if (length(s) == 1 && is.na(s)) morph_states else s
    })

    # Intersection of all sets
    Reduce(intersect, sets)
  }

  is_variable <- sapply(seq_len(ncol(mat)), function(i) {
    col <- mat[, i]
    col <- col[!is.na(col)]
    u <- unique(col)

    possible <- get_possible_states(u)

    # invariant if intersection has exactly one element
    !(length(possible) == 1)
  })

  deleted_columns <- which(!is_variable)

  cat("Deleted", length(deleted_columns), "invariant columns.\n")
  if (length(deleted_columns) > 0)
    cat("Columns deleted:", paste(deleted_columns, collapse = ", "), "\n")

  mat[, is_variable, drop = FALSE]

  # Output file name
  name = paste0(output_path, "_onlyVARIANTS.nexus")
  # Write a temporary file
  write.nexus.data(data, file=name, format="standard", interleaved=F)
  # Remove the interleave section to avoid errors in IQTREE
  temp = gsub("INTERLEAVE=NO", "", readLines(name))
  # Write a permanent file
  writeLines(temp, name)
}
