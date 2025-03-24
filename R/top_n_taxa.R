library(phyloseq) 
library(tidyverse) 

top_n_taxa <- function(physeq,
                       n = 10,
                       name = ShortName, # dynamically pass in the name of the common name column to use 
                       title = NA
                      ){

  ps <- physeq
  
  seqtab <- ps@otu_table %>% 
    data.frame() %>%
    rownames_to_column(var = "sample") %>% 
    pivot_longer(-sample, names_to = "asv", values_to = "count") %>%
    mutate(presence = ifelse(count > 0, 1, 0)) %>%
    group_by(asv) %>%
    summarise(prevalence = sum(presence) / n_distinct(sample) * 100, .groups = "drop")
  
  top_asvs <- slice_max(seqtab, order_by = prevalence, n = n)

  taxtab <- ps@tax_table %>% 
    data.frame() %>%
    rownames_to_column(var = "asv")
  
  top_taxa <- taxtab %>% 
    filter(asv %in% top_asvs$asv) %>%
    left_join(top_asvs, by = "asv") %>% # add prevalence column 
    mutate(lowestLevel = coalesce({{name}}, species, genus, family, order, class, phylum, superkingdom),
           lowestLevel = case_when(
             is.na(lowestLevel) ~ paste0("NA", row_number()),
             TRUE ~ lowestLevel
           )) %>%
    arrange(-prevalence)
  
  top_taxa_plot <- top_taxa %>% 
    ggplot(aes(x = fct_reorder(lowestLevel, prevalence, .desc = TRUE), y = prevalence)) + 
    geom_bar(stat = "identity") + 
    labs(x = "", 
         y = "% samples w/ taxa") + 
    theme(title = element_text(size = 16, face = "bold"), 
          axis.title = element_text(size = 14), 
          axis.text = element_text(size = 12), 
          axis.text.x = element_text(angle = 45, hjust = 1))

  if(is.na(title)) { 
      top_taxa_plot <- top_taxa_plot + 
        labs(title = paste0("Top ", n, " taxa"))
    } else { 
      top_taxa_plot <- top_taxa_plot + 
        labs(title = title) 
      }
  
  return(list(top_taxa, top_taxa_plot)) 
}
