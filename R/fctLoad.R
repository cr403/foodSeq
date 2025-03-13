library(phyloseq) 
library(tidyverse) 

fctLoad <- function(physeq, 
                    pca_output, 
                    PCx = 1,
                    nTaxa = 10){
  
  # Extract taxa information 
  taxtab <- ps.plot@tax_table %>%
    data.frame() %>%
    rownames_to_column(var = "asv") %>%
    mutate(lowestLevel = coalesce(species, genus, family, order, class, phylum, superkingdom)) 
  
  # Extract loadings 
  load.df <- pca_output$loadings %>% 
    data.frame() %>%
    rownames_to_column(var = "asv") %>%
    left_join(taxtab, by = "asv") 
  
  # Order loadings by absolute value 
  PC_col <- sym(paste0("PC", PCx)) 
  
  result <- load.df %>% 
    slice_max(order_by = abs(-!!PC_col), n = nTaxa) %>%
    relocate(CommonName, .after = asv) %>%
    relocate(PC_col, .after = CommonName) %>%
    mutate(CommonName = case_when(
      is.na(CommonName) ~ lowestLevel, 
      is.na(CommonName) ~ paste0("NA", row_number()), 
      TRUE ~ CommonName))
  
  # Plot top nTaxa 
  plot <- result %>% 
    ggplot(aes(x = fct_reorder(CommonName, !!PC_col), y = !!PC_col)) + 
    geom_bar(stat = "identity") + 
    coord_flip() + 
    theme(axis.title = element_text(size = 16, face = "bold"), 
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14))
  
  return(list(load.df = result, plot = plot)) 
  
}
