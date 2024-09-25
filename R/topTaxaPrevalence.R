```{r}
library(phyloseq)
library(dplyr)
library(ggplot2)
library(tidyr)

# Function to plot the top n most prevalent taxa
topTaxaPrevalence <- function(phyloseq_obj, n, title) {
  
  # Melt the phyloseq object to long format
  phyloseq_long <- psmelt(phyloseq_obj)
  
  # Coalesce the taxonomic levels
  phyloseq_long <- phyloseq_long %>%
    mutate(lowestLevel = coalesce(species, genus, family, order, class, phylum, superkingdom))
  
  # Calculate prevalence for each OTU
  prevalence_df <- phyloseq_long %>%
    group_by(OTU, lowestLevel) %>%
    summarise(prevalence = sum(Abundance > 0) / n_distinct(Sample)) %>%
    ungroup()
  
  # Get the top n most prevalent taxa
  top_taxa <- prevalence_df %>%
    top_n(n, prevalence) %>%
    arrange(desc(prevalence))
  
  # Reorder the factor levels of lowestLevel based on prevalence
  top_taxa <- top_taxa %>%
      mutate(lowestLevel = factor(lowestLevel, levels = unique(reorder(lowestLevel, -prevalence))))
  
  # Plot the top n most prevalent taxa with y-axis as percentage
  topTaxaPlot <- ggplot(top_taxa, aes(x = reorder(lowestLevel, -prevalence), y = prevalence * 100, fill=lowestLevel)) +
    geom_bar(stat = "identity") +
    labs(x = "Taxa", y = "Percent of Samples", title = title, fill = "lowestLevel") +
    # scale_fill_discrete(guide = guide_legend(reverse = FALSE)) +  # Reverse the legend order
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Format y-axis as percentage
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(topTaxaPlot)
  
  # example input: 
  topTaxaPrevalence(ps, 10, "Top 10 Taxa at Baseline")
}
```
