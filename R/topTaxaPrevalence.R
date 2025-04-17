#' @title Top Taxa Prevalence
#'
#' @description finds and plots the top n most prevalent taxa in a phyloseq object
#'
#' @param phyloseq_obj phyloseq of interest
#' @param n number of taxa to plot
#' @param title optional title of plot
#' @param amplicon trnl or 12s
#'
#' @return bar plot of top n taxa
#' @return top.taxa = data frame of top taxa
#'
#' @export
# Function to plot the top n most prevalent taxa
topTaxaPrevalence <- function(phyloseq_obj, n, title, amplicon) {

  # Subset the phyloseq object for the desired amplicon
  if(amplicon == "trnl") {
    phyloseq_obj <- phyloseq_obj %>%
      subset_samples(amplicon=="trnl")
  }
  if(amplicon == "12s") {
    phyloseq_obj <- phyloseq_obj %>%
      subset_samples(amplicon=="12s")
  }
  else(phyloseq_obj <- phyloseq_obj)

  # Melt the phyloseq object to long format
  phyloseq_long <- psmelt(phyloseq_obj)

  # Coalesce the taxonomic levels
  phyloseq_long <- phyloseq_long %>%
    mutate(lowestLevel = coalesce(species, genus, family, order, class, phylum, superkingdom)) %>%
    filter(!is.na(lowestLevel))

  # Calculate prevalence for each OTU (accounts for distinct OTU's assigning to the same lowestLevel)
  prevalence_df <- phyloseq_long %>%
    group_by(lowestLevel, Sample) %>%
    summarise(Abundance = sum(Abundance > 0)) %>%  # Check if any OTU for this lowestLevel is present in the sample
    ungroup() %>%
    group_by(lowestLevel) %>%
    summarise(prevalence = sum(Abundance > 0) / n_distinct(Sample)) %>%  # Recalculate prevalence
    ungroup()

  ## Get the top n most prevalent taxa
  # Use this if there are a lot of ties within prevalence -- will appear if there are more than 10 taxa on graph
  top_taxa <- prevalence_df %>%
    slice_max(prevalence, n = n, with_ties = FALSE) %>%
    arrange(desc(prevalence))
  ## Use this if you want equal prevalence values to be counted together towards n
  #top_taxa <- prevalence_df %>%
  #  top_n(n, prevalence) %>%
  #  arrange(desc(prevalence))

  # Reorder the factor levels of lowestLevel based on prevalence
  top_taxa <- top_taxa %>%
    mutate(lowestLevel = factor(lowestLevel, levels = unique(reorder(lowestLevel, -prevalence))))

  # Plot the top n most prevalent taxa with y-axis as percentage
  topTaxaPlot <- ggplot(top_taxa, aes(x = reorder(lowestLevel, -prevalence), y = prevalence * 100, fill=lowestLevel)) +
    geom_bar(stat = "identity") +
    labs(y = "% of Samples", title = title) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Format y-axis as percentage
    ylim(0,100) +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          title = element_text(size = 14))

  return(list(plot = topTaxaPlot, top.taxa = top_taxa))

  # example input:
  topTaxaPrevalence(ps, 10, "Top 10 Taxa at Baseline", "trnl")
}
