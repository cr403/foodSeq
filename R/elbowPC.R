#' --------------------
#' Elbow Method
#' Code written by Caroline
#' --------------------
#'
#' This function performs Elbow Method for Principal Component importance
#'
#' @param pcaOutput object assigned to the output of pcaPlot()
#'
#' @return scree plot from pcaPlot() with suggested cutoff marked by a dotted red line
#'
#' @export
elbowPC <- function(pcaOutput # object assigned to the output of pcaPlot()
){
  eigenvalues <- pcaOutput$scree.table %>%
    pull(Eigenvalue)

  n <- length(eigenvalues)
  points <- cbind(1:n, eigenvalues)

  # line from first to last point
  line_vec <- points[n, ] - points[1, ]
  line_vec <- line_vec / sqrt(sum(line_vec^2))

  # distance from each point to the line
  distances <- apply(points, 1, function(point){
    vec <- point - points[1, ]
    proj <- sum(vec * line_vec) * line_vec
    orth_vec <- vec - proj
    sqrt(sum(orth_vec^2))
  })

  elbow <- which.max(distances)

  # Add line to scree plot and print elbow point
  pcaOutput$scree.plot +
    geom_vline(xintercept = elbow, color = "red", linetype = "dashed") +
    annotate("text", x = 100, y = 8, label = paste0("Elbow = PC", elbow), hjust = 1, vjust = -0.5, color = "red")
}
