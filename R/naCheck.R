#' @title NA Checker
#'
#' @description Checks a phyloseq object for any taxa that have no superkingdom assignment
#'
#' @param physeq phyloseq object of interest
#'
#' @return for each study in the phyloseq object, returns the number of NA's present
#'
#' @export
naCheck <- function(physeq){
  studies <- physeq@sam_data %>% data.frame() %>% pull(study) %>% unique()

  for(i in studies) {

    df <- physeq %>%
      prune_samples(get_variable(physeq, "study") %in% i, .) %>%
      prune_taxa(taxa_sums(.) > 0, .) %>%
      tax_table() %>%
      data.frame() %>%
      filter(is.na(superkingdom))

    cat(i, " has ", nrow(df), " NA's", "\n")
  }
}
