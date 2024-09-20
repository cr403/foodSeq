--------------------
PCA plot and biplot 
Code written by Ben 
--------------------

```{r}
pcaPlot<-function(ps, # clr transformed and filtered data
                  colorVar, # variable from samdf to color samples by
                  colorName, # what to display variable name as in legend
                  nTaxa # number of taxa to display
                  ) { 
  if("name" %in% colnames(data.frame(ps@sam_data) )) { # there will be an error if there is a column called "name" in the ps object
  sample_data(ps)=ps@sam_data%>%
    data.frame() %>%
    dplyr::rename(name.x=name)
  }
  samdf <- data.frame(ps@sam_data) %>%
    rownames_to_column(var = 'name')
  
  # PCA
  pca <- prcomp(ps@otu_table, center = TRUE, scale = FALSE)
  
  pca.df <- data.frame(pca$x) %>% 
          rownames_to_column(var = 'name')
  
  # % variance explained
  eigs <- pca$sdev^2
  varExplained <- 100 * round(eigs/sum(eigs), 5)
  
  # pull out first 2 PC variance explained for the plot
  ve.pc1 <- as.character(round(varExplained[1], 3))
  ve.pc2 <- as.character(round(varExplained[2], 3))

  
  # Add back sample data
  pca.df <- left_join(pca.df, samdf)

  # Calculate plotting limits based on largest value observed in PC axes 1 and 2
  limit <- max(abs(pca.df[, c('PC1', 'PC2')])) +
            0.05*(max(abs(pca.df[, c('PC1', 'PC2')])))
  
  

  pca.plot <- 
       ggplot(pca.df, aes_string(x = "PC1", y = "PC2", color = colorVar)) +
       geom_point(size = 2, alpha = 0.5) +
       stat_ellipse(level = 0.95, aes_string(group = colorVar), linetype = "dashed") + # Add ellipses here
       coord_equal() +
       labs(x = paste0(' PC1 (', ve.pc1, '%)'),
            y = paste0(' PC2 (', ve.pc2, '%)')) + 
       xlim(-limit, limit) + ylim(-limit, limit)+
       theme_classic() +
       theme(axis.line = element_line(size = 1, color = 'black'),
             axis.ticks = element_line(color = 'black'),
             axis.title = element_text(size = 14, face = 'bold', color = 'black'),
             # axis.text = element_blank(),
             # legend.background = element_blank(),
             # legend.title = element_blank(),
             #legend.position = "none"
             # legend.text = element_text(size = 10, face = 'bold'),
             ) 
        # scale_color_manual(values=met.brewer('Signac',10)) ### can change color palette here 

 pca.plot
  # Biplot
  
  # Calculate loadings
  V <- pca$rotation # Eigenvectors
  L <- diag(pca$sdev) # Diag mtx w/sqrts of eigenvalues on diag.
  loadings <- V %*% L
       
  # Get loadings for first 2 PCs and format for plotting
  pythag <- function(a, b){sqrt(a^2 + b^2)}
  loadings.12 <- data.frame(loadings[, 1:2]) %>%
       dplyr::rename(PC1 = X1, PC2 = X2) %>% 
       mutate(variable = row.names(loadings)) %>% 
       mutate(length = pythag(PC1, PC2), slope = PC2/PC1, ang = atan(slope)*(180/pi))
  
  loadings.plot <- top_n(loadings.12, nTaxa, wt = length) 
  
   loadings.plot <- top_n(loadings.12, nTaxa, wt = length) 
  
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
       dplyr::select(asv,name) %>% 
       right_join(loadings.plot, by = c('asv' = 'variable'))
  
  # What quadrant of the plot is the label in?
  q1 <- filter(loadings.plot, PC1 > 0 & PC2 > 0)
  q2 <- filter(loadings.plot, PC1 < 0 & PC2 > 0)
  q3 <- filter(loadings.plot, PC1 < 0 & PC2 < 0)
  q4 <- filter(loadings.plot, PC1 > 0 & PC2 < 0)
       
  pca.biplot <- 
       pca.plot + 
       geom_segment(data = loadings.plot,
                    aes(x = 0, y = 0, 
                        xend = PC1, yend = PC2),
                    color = 'black',
                    arrow = arrow(angle = 15, 
                                  length = unit(0.1, 'inches'))) + 
    labs(color = colorName)
  
  # Then add geom_text quadrant-by-quadrant, aligning text accordingly
       if (dim(q1)[1] != 0) {
            pca.biplot <- pca.biplot +
                 geom_text(data = q1, aes(x = PC1, y = PC2, hjust = 0, angle = ang,
                                          label=paste0('   ', name),
                                          fontface = 'bold'),
                           color = 'black', show.legend = FALSE)
       }
       if (dim(q2)[1] != 0) {
            pca.biplot <- pca.biplot +
                 geom_text(data = q2, aes(x = PC1, y = PC2, hjust = 1, angle = ang,
                                          label=paste0(name, '   '),
                                          fontface = 'bold'),
                           color = 'black', show.legend = FALSE)
       }
       if (dim(q3)[1] != 0) {
            pca.biplot <- pca.biplot +
                 geom_text(data = q3, aes(x = PC1, y = PC2, hjust = 1, angle = ang,
                                          label=paste0(name, '   '),
                                          fontface = 'bold'),
                           color = 'black', show.legend = FALSE)
       }
       if (dim(q4)[1] != 0) {
            pca.biplot <- pca.biplot +
                 geom_text(data = q4, aes(x = PC1, y = PC2, hjust = 0, angle = ang,
                                          label=paste0('   ', name),
                                          fontface = 'bold'),
                           color = 'black', show.legend = FALSE)
       }
  
  #print(pca.biplot)
  
  return(list(pca.df = pca.df, pca.biplot = pca.biplot, loadings = loadings))
}
```
