library(tidyverse) 
library(phyloseq) 

naCheck <- function(physeq){
  studies <- physeq@sam_data %>% data.frame() %>% pull(study) %>% unique() 

for(i in studies) {
  
  df <- physeq %>% 
    prune_samples(get_variable(physeq, "study") %in% i, .) %>%
    prune_taxa(taxa_sums(.) > 0, .) %>% 
    tax_table() %>% 
    data.frame() %>% 
    filter(is.na(superkingdom)) 
  
  cat(i, " has ", nrow(df), " NA's", "\n")
  }
}
