--------------------
PCA plot and biplot 
Code written by Ben, adapted by Caroline 
--------------------

```{r}
pcaPlot <- function(ps, # clr transformed and filtered data
                    colorVar = NULL, # variable from samdf to color samples by
                    colorName = NULL, # what to display variable name as in legend
                    nTaxa = 10, # number of taxa to display
                    customColors = NULL, # optional named vector of colors
                    xPC = 1, # Principal Component for x-axis
                    yPC = 2,  # Principal Component for y-axis
                    ellipse = FALSE # optional add centroid ellipses 
                    ) { 
  
  if ("name" %in% colnames(data.frame(ps@sam_data))) { 
    # Prevent conflict with 'name' column
    sample_data(ps) <- ps@sam_data %>%
      data.frame() %>%
      dplyr::rename(name.x = name)
  }
  
  samdf <- data.frame(ps@sam_data) %>%
    rownames_to_column(var = 'name')
  
  # PCA
  pca <- prcomp(ps@otu_table, center = TRUE, scale = FALSE)
  
  # % variance explained
  eigs <- pca$sdev^2
  varExplained <- 100 * eigs / sum(eigs)
  names(varExplained) <- paste0('PC', seq_along(varExplained))
  
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
  
  return(list(pca.df = pca.df, pca.biplot = pca.biplot, loadings = loadings))
}
```
