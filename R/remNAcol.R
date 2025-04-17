#' @title Remove NA Columns
#'
#' @description removes columns with all NA rows from the sample data table of a phyloseq object
#'
#' @param ps phyloseq object with sample data table
#'
#' @return ps phyloseq object with NA columns removed
#'
#' @export
remNAcol <- function(ps # phyloseq object
) {
  object_name <- deparse(substitute(ps)) # Captures the name of the phyloseq object

  before_removal <- dim(ps@sam_data)

  df <- ps@sam_data %>%
    data.frame() %>%
    select(where(~ !all(is.na(.))))

  ps@sam_data <- sample_data(df)

  after_removal <- dim(ps@sam_data)

  print(paste0("Dimensions of ", object_name, "@sam_data ", "before removal: ", before_removal[1], " ", before_removal[2]))
  print(paste0("Dimensions of ", object_name,"@sam_data ", "after removal: ", after_removal[1], " ", after_removal[2]))
  print(paste0("Columns removed: ", before_removal[2] - after_removal[2]))

  return(ps)
}
