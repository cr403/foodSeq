library(ggridges) 

pcDens <- function(data, # pcaPlot() output
         x = "PC1", # PC of interest 
         y = "country", # grouping variable
         title = NULL,
         scale = 2, 
         alpha = 0.7, 
         title.size = 14,
         text.size = 12
){
  data$pca.df %>%
    data.frame() %>% 
    ggplot(aes(x = .data[[x]], y = fct_reorder(.data[[y]], -.data[[x]]), fill = .data[[y]])) + 
    geom_density_ridges(scale = 2, alpha = 0.7) + 
  labs(title = title) +
    theme(axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12), 
          axis.title.y = element_blank(), 
          legend.position = "none")
}
