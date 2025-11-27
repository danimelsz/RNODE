#' @title splitNoStates
#' @name splitNoStates
#' @description Splits a morphological matrix according to the number of character-states. This procedure has been recommended to run phylogenetic analyses with the MK and MKv models (the 'K' refers to the number of states). Khakurel et al. (2024) demonstrated that MK models with high K values can understimate the branch lengths, whereas MK models with small K values can overstimate them. As such, some recent studies have partitioned morphological characters according to their number of states (e.g. Černý & Simonoff 2023).
#' @author Daniel YM Nakamura
#' @param input Input file (morphological matrix, either loaded previously in R or loaded from a local file).
#' @param inpu_format To load from a local file: 'nexus' or 'tnt'.
#' @param ambiguity_addState If T, ambiguities are counted as additional character-states (default: ambiguity_addState = F).
#' @param inapplicable_addState If T, inapplicable states are counted toward the sum of unique character states.
#' @param output_index Output index (e.g. if the user specify it as "Desktop/Index", the output files will be "Desktop/Index_ORDERED.nexus" and "Desktop/Index_UNORDERED.nexus")
#' @examples
#' # Example
#' splitNoStates(input="../testdata/048_MORPH_data.nex", output_index = "../testdata/048_MORPH", list_ordered=list_ordered)
#'
#' @references Černý, D., & Simonoff, A. L. (2023). Statistical evaluation of character support reveals the instability of higher-level dinosaur phylogeny. Scientific Reports, 13(1), 9273.
#' @references Khakurel, B., Grigsby, C., Tran, T. D., Zariwala, J., Höhna, S., & Wright, A. M. (2024). The fundamental role of character coding in Bayesian morphological phylogenetics. Systematic biology, 73(5), 861-871.
#' @export
#'
