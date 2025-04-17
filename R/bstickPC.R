#' @title Broken stick analysis
#'
#' @description This code runs Broken Stick method for Principal Component selection
#'
#' @param pcaOutput object assigned to output of pcaPlot()
#'
#' @return scree plot from pcaPlot() with dotted red line indicating suggested importance cutoff
#'
#' @export
bstickPC <- function(pcaOutput # object assigned to output of pcaPlot()
){
  eigenvalues <- pcaOutput$scree.table %>%
    pull(Eigenvalue)

  # Compute broken stick distribution
  bs_values <- vegan::bstick(n = length(eigenvalues))

  # Determine number of components to keep
  bs_k <- sum(eigenvalues > bs_values)

  # Add broken stick line to scree plot
  pcaOutput$scree.plot +
    geom_vline(xintercept = bs_k, color = "red", linetype = "dashed") +
    annotate("text", x = 100, y = 8, label = paste0("Retain PC's before PC", bs_k), hjust = 1, vjust = -0.5, color = "red")
}
