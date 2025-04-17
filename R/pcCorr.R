#' @title PC correlations
#'
#' @description makes scatterplot with PC and variable of interest
#'
#' @param data pcaPlot()$pca.df
#' @param xVar variable name of interest (e..g, PC1)
#' @param yVar variable name of interest (e.., age)
#' @param colorVar variable by which to color points
#' @param continuous T/F
#' @param mid determines midpoint of gradient color scale
#'
#' @export
pcCorr <- function(data,
         xVar,
         yVar,
         colorVar,
         continuous = TRUE,
         mid = 50) {
  plot <- data %>%
    ggplot(aes(x = .data[[xVar]], y = .data[[yVar]], color = .data[[colorVar]])) +
    geom_point() +
    geom_smooth(method = "lm") +
    ggpmisc::stat_poly_eq(
      formula = y ~ x,
      aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
      color = "blue",
      size = 6,
      parse = TRUE
    ) +
    theme(axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12),
          legend.position = "none")

  if(continuous) {
    plot <- plot + scale_color_gradient2(low = "navy", mid = "gold", high = "red", midpoint = mid)
  }

  return(plot)
}
