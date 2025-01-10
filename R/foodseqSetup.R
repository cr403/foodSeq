```{r]
foodseqSetup <- function(ps) {
  ### Relative abundance generation 
  ps.ra <- transform_sample_counts(ps, function(x){x/sum(x)})
  
  ### CLR Transform and Filter out non-foods/NA's 
  ps.filt.clr <- ps %>% 
    prune_samples(sample_sums(.) >0, .) %>% # Remove samples that do not have any reads, will mess up PCA plot 
    microbiome::transform(transform = "compositional") %>% 
    microbiome::transform(transform = "clr") %>% # clr transform 
    subset_taxa(!is.na(superkingdom)) # Remove NA's
  
  ### Update read counts in phyloseq object
  sample_data(ps.filt.clr)$reads <- sample_sums(ps.filt.clr)
  
  ### Alpha diversity metrics 
  pMR <- ifelse(subset_taxa(ps.filt.clr, !is.na(superkingdom))@otu_table>0,1,0) %>%      # Converts OTU table values to presence/absence
  rowSums()
  ps.filt.clr@sam_data$pMR <- pMR
  
  # This is done with raw relative abundance, not clr transformed or raw abundance data 
  shan <- estimate_richness(ps.ra, measures = "Shannon") %>%
    rownames_to_column(var = "test") 
  
  sample_data(ps.ra) <- ps.ra@sam_data %>%
    data.frame() %>%
    rownames_to_column(var = "test") %>%
    left_join(shan, by = "test") %>% 
    column_to_rownames(var = "test") %>%
    sample_data()
  
  # ### Beta diversity metrics
  # bc_dist <- phyloseq::distance(ps.ra, method = "bray") # non-clr filtered
  # bc_matrix <- as.matrix(bc_dist)
  
  return(list(ps.ra = ps.ra, # Relative abundance 
              ps.filt = ps.filt, # Foods only 
              ps.filt.clr = ps.filt.clr # Foods only, CLR-transformed 
              #bc_dist = bc_dist 
              #bc_matrix = bc_matrix
              ))
}
```
