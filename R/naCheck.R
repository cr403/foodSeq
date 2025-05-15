#' @title NA Checker
#'
#' @description Checks a phyloseq object for any taxa that have no superkingdom/kingdom assignment
#'
#' @param physeq phyloseq object of interest
#' @param trnl trnl or 12sv5
#' @param group grouping variable (e.g., study, location, etc.)
#'
#' @return for each study in the phyloseq object, returns the number of NA's present
#'
#' @export
naCheck <- function(physeq,
                    trnl = TRUE,
                    group = NULL){

  if (!is.null(group)) {
    groups <- physeq@sam_data %>%
      data.frame() %>%
      pull(.data[[group]]) %>%
      unique()
  } else {
    groups <- 1
  }

  if (trnl) {
    taxlevel <- "superkingdom"
  } else {
    taxlevel <- "kingdom"
  }

  for(i in groups) {

    df <- physeq %>%
      prune_samples(get_variable(physeq, .data[[group]]) %in% i, .) %>%
      prune_taxa(taxa_sums(.) > 0, .) %>%
      tax_table() %>%
      data.frame()
      filter(is.na(.data[[taxlevel]]))

    cat(i, " has ", nrow(df), " NA's", "\n")
  }
}
