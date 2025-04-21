#' @title Density Plot
#'
#' @description generates density plot for two variables of interest
#'
#' @param data pcaPlot() output
#' @param x "PC" of interest
#' @param y variable of interest (e.g., country)
#' @param title optional input of title
#' @param scale = 2
#' @param alpha = 0.7
#' @param title.size = 14
#' @param text.size = 12
#' @param customColors optional color palette inclusion
#'
#' @return density plot
#'
#' @export
pcDens <- function(data, # pcaPlot() output
                   x = "PC1", # PC of interest
                   y = "country", # grouping variable
                   title = NULL,
                   scale = 2,
                   alpha = 0.7,
                   title.size = 14,
                   text.size = 12,
                   customColors = NULL
){
  plot <- data$pca.df %>%
    data.frame() %>%
    ggplot(aes(x = .data[[x]], y = fct_reorder(.data[[y]], -.data[[x]]), fill = .data[[y]])) +
    ggridges::geom_density_ridges(scale = 2, alpha = 0.7) +
    labs(title = title) +
    theme(axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12),
          axis.title.y = element_blank(),
          legend.position = "none")

  if (!is.null(customColors)) {
    plot <- plot + ggplot2::scale_color_manual(values = customColors)
  } else {
    plot <- plot
  }

  return(plot)
}
