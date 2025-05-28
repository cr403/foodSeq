#' @title Factor Loading Analysis
#'
#' @description This function generates factor loading plots.
#'
#' @param physeq clr filtered phyloseq object on which pcaPlot() was run
#' @param pca_output output object from pcaPlot()
#' @param PCx principal component for x-axis of graph
#' @param nTaxa number of taxa to include on plot
#' @param name name of common name column
#' @param labWidth character length breakpoint for wrapping labels
#' @param title optional title
#' @param titleSize optional title size
#' @param axisTitleSize optional axis title size
#' @param textSize optional text size
#'
#' @return load.df = loadings data frame in descending order by PCx
#' @return plot = plot of top nTaxa taxa for PCx
#'
#' @export
fctLoad <- function(physeq,
                    pca_output,
                    PCx = 1,
                    nTaxa = 10,
                    name = "taxlabel",
                    labWidth = 60,
                    title = NULL,
                    titleSize = 18,
                    axisTitleSize = 16,
                    textSize = 14){

  # Extract loadings
  load.df <- pca_output$loadings %>%
    data.frame() %>%
    rownames_to_column(var = "asv")

  # Extract taxa information
  taxtab <- physeq@tax_table %>%
    data.frame() %>%
    rownames_to_column(var = "asv") %>%
    left_join(load.df, by = "asv") %>%
    mutate(label = .data[[name]],
           label = wrapLabels(label, labWidth))

  # Order loadings by absolute value
  PC_col <- sym(paste0("PC", PCx))

  result <- taxtab %>%
    slice_max(order_by = abs(-!!PC_col), n = nTaxa) %>%
    relocate(label, .after = asv) %>%
    relocate(PC_col, .after = label)

  # Plot top nTaxa
  plot <- result %>%
    ggplot(aes(x = fct_reorder(label, !!PC_col), y = !!PC_col)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme(
      plot.title = element_text(size = titleSize, face = "bold"),
      axis.title = element_text(size = axisTitleSize, face = "bold"),
      axis.title.y = element_blank(),
      axis.text = element_text(size = textSize))

  if(!is.null(title)) {
    plot <- plot + labs(title = title)
  }

  return(list(load.df = result, plot = plot))

}
