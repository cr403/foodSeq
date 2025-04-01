library(tidyverse) 
library(phyloseq) 

naCheck <- function(physeq){
  studies <- physeq@sam_data %>% data.frame() %>% pull(study) %>% unique() 

for(i in studies) {
  df <- newbn %>% 
    subset_samples(study %in% i) %>%
    prune_taxa(taxa_sums(.) > 0, .) %>% 
    tax_table() %>% 
    data.frame() %>% 
    filter(is.na(superkingdom)) 
  
  if(nrow(df) == 0) {
    cat(i, " does not have NA's", "\n") 
  } else { 
    cat(i, " has NA's", "\n") 
    }
  }
}
