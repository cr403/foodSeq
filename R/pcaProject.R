#' @title PCA Projection
#'
#' @description projects new data into existing PCA space
#'
#' @param new_ps new phyloseq object; with NA's, non-clr transformed
#' @param old_ps old phyloseq object; with NA's, non-clr transformed
#' @param pca_output pcaPlot()$pca_output
#'
#' @return pseudocount = number used when running clr transformation
#' @return ps.filt.clr = new_ps after going through clr trasnformation exactly the same as old_ps
#' @return projection.scores = data frame with PC scores for each sample in new_ps
#'
#' @export
pcaProject <- function(new_ps, # new phyloseq object; with NA's, NON-clr transformed
                       old_ps, # old phyloseq object; with NA's, NON-clr trasnformed
                       pca_output # pcaPlot()$pca_output
){
  # extract names of taxa from original phyloseq
  desired_taxa <- taxa_names(old_ps)

  # pull out the new counts & ensure taxa are columns
  mat_new <- as(otu_table(new_ps), "matrix")
  if (taxa_are_rows(otu_table(new_ps))) {
    mat_new <- t(mat_new)
  }

  # add taxa missing from new phyloseq in as zeros
  missing_taxa <- setdiff(desired_taxa, colnames(mat_new))
  if (length(missing_taxa)) {
    zero_mat <- matrix(
      0,
      nrow = nrow(mat_new),
      ncol = length(missing_taxa),
      dimnames = list(rownames(mat_new), missing_taxa)
    )
    mat_new <- cbind(mat_new, zero_mat)
  }

  # ensure that order of taxa is the same as in original phyloseq and drop taxa with no match
  mat_new <- mat_new[, desired_taxa, drop = FALSE]

  # build the combined tax_table by subsetting the OLD one
  taxtab_old <- as(tax_table(old_ps), "matrix")
  combined_taxtab <- taxtab_old[desired_taxa, , drop = FALSE]

  # regenerate phyloseq that's been cleaned
  new_ps_clean <- phyloseq(
    otu_table(mat_new, taxa_are_rows = FALSE),
    tax_table(combined_taxtab),
    sample_data(new_ps)
  )

  # extract psuedocount used in clr transformation of old phyloseq
  comp_mat <- microbiome::transform(old_ps, transform = "compositional") %>%
    otu_table() %>%
    as.matrix() %>%
    as.vector()

  pseudocount <- 0.65 * min(comp_mat[comp_mat > 0], na.rm = TRUE)

  # manually perform clr transformation using psuedocount from old phyloseq transformation
  mat_count <- as.matrix(otu_table(new_ps_clean))
  comp <- sweep(mat_count, 1, rowSums(mat_count), "/")
  comp <- comp + pseudocount
  gm <- exp(rowMeans(log(comp)))
  clr_mat <- log(comp / gm)

  # filter out non-foods
  taxtab <- as.matrix(tax_table(new_ps_clean))
  food_taxa <- rownames(taxtab)[!is.na(taxtab[, "superkingdom"])]
  clr_mat <- clr_mat[, food_taxa, drop = FALSE]

  # Project into old pca space
  clr_mat <- clr_mat[, rownames(pca_fit$rotation), drop = FALSE]
  scores <- predict(pca_fit, newdata = clr_mat)

  return(list(pseudocount = pseudocount,
              ps.filt.clr = ps.filt.clr,
              projection.scores = scores))

}
