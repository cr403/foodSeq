#'
#' @title PCA plot and biplot
#'
#' @description This function runs a Principal Component Analysis.
#'
#' @param ps clr transformed and filtered data
#' @param colorVar variable from samdf to color samples by
#' @param colorName what to display variable name as in legend
#' @param nTaxa number of taxa to display
#' @param customColors optional named vector of colors
#' @param customGradient optional df of low/high colors for gradient
#' @param mid for customGradient -- midpoint to median, mean, middle of range
#' @param xPC Principal Component for x-axis
#' @param yPC Principal Component for y-axis
#' @param ellipse optional add centroid ellipses
#' @param bplab name of column that you want to use for labeling the biplot
#'
#' @return pca biplot with nTaxa factor loadings displayed
#' @return samdf with all pc's added as columns
#' @return data frame with loadings for pca space
#' @return scree plot
#' @return table with scree plot numbers
#' @return output of prcomp()
#' @export
pcaPlot <- function(ps, # clr transformed and filtered data
                    colorVar = NULL, # variable from samdf to color samples by
                    colorName = NULL, # what to display variable name as in legend
                    nTaxa = 10, # number of taxa to display
                    customColors = NULL, # optional named vector of colors
                    customGradient = NULL, # optional df of low/high colors for gradient
                    mid = "mean", # for customGradient -- midpoint to median, mean, middle of range
                    xPC = 1, # Principal Component for x-axis
                    yPC = 2,  # Principal Component for y-axis
                    ellipse = FALSE, # optional add centroid ellipses
                    bplab = NULL # optional variable name for biplot arrow labels
) {

  if ("name" %in% colnames(data.frame(ps@sam_data))) {
    # Prevent conflict with 'name' column
    sample_data(ps) <- ps@sam_data %>%
      data.frame() %>%
      dplyr::rename(name.x = name)
  }

  if (!is.null(bplab)) {
    # Relocate bplab to the end of tax table so it will be used for labeling
    tax_table(ps) <- ps@tax_table %>%
      data.frame() %>%
      relocate(.data[[bplab]], .after = everything()) %>%
      as.matrix() %>%
      tax_table()
  }

  samdf <- data.frame(ps@sam_data) %>%
    rownames_to_column(var = 'name')

  # PCA
  pca <- prcomp(ps@otu_table, center = TRUE, scale = FALSE)

  # % variance explained
  eigs <- pca$sdev^2
  varExplained <- 100 * eigs / sum(eigs)
  names(varExplained) <- paste0('PC', seq_along(varExplained))

  # Create a scree table with eigenvalues, variance explained, and cumulative variance
  scree.table <- data.frame(
    PC = paste0("PC", seq_along(varExplained)),
    Eigenvalue = eigs,
    VarianceExplained = varExplained,
    CumulativeVariance = cumsum(varExplained)
  )

  # Generate a scree plot using ggplot2
  scree.plot <- ggplot(scree.table, aes(x = as.numeric(gsub("PC", "", PC)), y = VarianceExplained)) +
    geom_line() +
    geom_point() +
    labs(title = "Scree Plot", x = "Principal Component", y = "Variance Explained (%)") +
    theme_classic()

  # Extract variance explained for specified PCs
  ve.xPC <- as.character(round(varExplained[paste0('PC', xPC)], 3))
  ve.yPC <- as.character(round(varExplained[paste0('PC', yPC)], 3))

  # PCA scores
  pca.df <- data.frame(pca$x) %>%
    rownames_to_column(var = 'name')

  # Add back sample data
  pca.df <- left_join(pca.df, samdf, by = "name")

  # Calculate plotting limits based on specified PCs
  limit <- max(abs(pca.df[, c(paste0('PC', xPC), paste0('PC', yPC))])) +
    0.05 * max(abs(pca.df[, c(paste0('PC', xPC), paste0('PC', yPC))]))

  # Initialize PCA plot
  pca.plot <- ggplot(pca.df, aes_string(x = paste0('PC', xPC), y = paste0('PC', yPC), color = colorVar)) +
    geom_point(size = 2, alpha = 0.5) +
    coord_equal() +
    labs(x = paste0('PC', xPC, ' (', ve.xPC, '%)'),
         y = paste0('PC', yPC, ' (', ve.yPC, '%)')) +
    xlim(-limit, limit) + ylim(-limit, limit) +
    theme_classic() +
    theme(axis.line = element_line(size = 1, color = 'black'),
          axis.ticks = element_line(color = 'black'),
          axis.title = element_text(size = 14, face = 'bold', color = 'black'))

  # Add custom color scale if provided
  if (!is.null(customColors)) {
    pca.plot <- pca.plot + scale_color_manual(values = customColors)
  }

  # Define a function to calculate the midpoint
  midpoint <- switch(mid,
                     "middle" = mean(range(pca.df[[colorVar]], na.rm = TRUE)),  # Mean of min and max
                     "median" = median(pca.df[[colorVar]], na.rm = TRUE),       # Median
                     "mean" = mean(pca.df[[colorVar]], na.rm = TRUE),           # Mean
                     stop("'mid' must be one of 'middle', 'median', or 'mean'") # Error for invalid 'mid'
  )

  # Add custom color gradient if provided
  if (!is.null(customGradient)) {
    pca.plot <- pca.plot + scale_color_gradient2(low = paste0(customGradient$low),
                                                 mid = paste0(customGradient$mid),
                                                 high = paste0(customGradient$high),
                                                 midpoint = midpoint)
  }

  # Add optional ellipses
  if (ellipse) {
    pca.plot <- pca.plot + stat_ellipse(level = 0.95, aes_string(group = colorVar), linetype = "dashed")
  }

  # Calculate loadings
  V <- pca$rotation # Eigenvectors
  L <- diag(pca$sdev) # Diagonal matrix with square roots of eigenvalues
  loadings <- V %*% L
  colnames(loadings) <- colnames(V)  # Assign column names to loadings

  # Get loadings for specified PCs and format for plotting
  loadings.xy <- data.frame(loadings[, c(paste0('PC', xPC), paste0('PC', yPC))]) %>%
    dplyr::rename(PCx = paste0('PC', xPC), PCy = paste0('PC', yPC)) %>%
    mutate(variable = row.names(loadings),
           length = sqrt(PCx^2 + PCy^2),
           ang = atan2(PCy, PCx) * (180 / pi))

  loadings.plot <- top_n(loadings.xy, nTaxa, wt = length)

  # Adjust angles to keep labels upright
  loadings.plot <- loadings.plot %>%
    mutate(adj_ang = ifelse(ang < -90, ang + 180,
                            ifelse(ang > 90, ang - 180, ang)))

  # Rename loadings with lowest taxonomic level
  loadings.taxtab <- tax_table(ps)[row.names(loadings.plot)] %>%
    data.frame()
  loadings.taxtab <- loadings.taxtab[cbind(1:nrow(loadings.taxtab), max.col(!is.na(loadings.taxtab), ties.method = 'last'))] %>%
    data.frame()
  colnames(loadings.taxtab) <- c("name")
  loadings.taxtab$asv <- tax_table(ps)[row.names(loadings.plot)] %>%
    data.frame() %>%
    rownames()

  loadings.plot <- loadings.taxtab %>%
    dplyr::select(asv, name) %>%
    right_join(loadings.plot, by = c('asv' = 'variable'))

  # Determine the quadrant of each label
  q1 <- filter(loadings.plot, PCx > 0 & PCy > 0)
  q2 <- filter(loadings.plot, PCx < 0 & PCy > 0)
  q3 <- filter(loadings.plot, PCx < 0 & PCy < 0)
  q4 <- filter(loadings.plot, PCx > 0 & PCy < 0)

  pca.biplot <-
    pca.plot +
    geom_segment(data = loadings.plot,
                 aes(x = 0, y = 0,
                     xend = PCx, yend = PCy),
                 color = 'black',
                 arrow = arrow(angle = 15,
                               length = unit(0.1, 'inches'))) +
    labs(color = colorName)

  # Add geom_text for each quadrant with adjusted angle and justification
  if (nrow(q1) != 0) {
    pca.biplot <- pca.biplot +
      geom_text(data = q1, aes(x = PCx, y = PCy, hjust = 0, vjust = 0, angle = adj_ang,
                               label = name,
                               fontface = 'bold'),
                color = 'black', show.legend = FALSE)
  }
  if (nrow(q2) != 0) {
    pca.biplot <- pca.biplot +
      geom_text(data = q2, aes(x = PCx, y = PCy, hjust = 1, vjust = 0, angle = adj_ang,
                               label = name,
                               fontface = 'bold'),
                color = 'black', show.legend = FALSE)
  }
  if (nrow(q3) != 0) {
    pca.biplot <- pca.biplot +
      geom_text(data = q3, aes(x = PCx, y = PCy, hjust = 1, vjust = 1, angle = adj_ang,
                               label = name,
                               fontface = 'bold'),
                color = 'black', show.legend = FALSE)
  }
  if (nrow(q4) != 0) {
    pca.biplot <- pca.biplot +
      geom_text(data = q4, aes(x = PCx, y = PCy, hjust = 0, vjust = 1, angle = adj_ang,
                               label = name,
                               fontface = 'bold'),
                color = 'black', show.legend = FALSE)
  }

  return(list(pca.df = pca.df,
              pca.biplot = pca.biplot,
              loadings = loadings,
              pca.output = pca,
              scree.table = scree.table,
              scree.plot = scree.plot))
}
