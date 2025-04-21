#' @title Top Taxa
#'
#' @description finds and plots top taxa from a phyloseq
#'
#' @param physeq phyloseq object
#' @param n number of taxa to plot
#' @param name name of common name column (e.g., ShortName)
#' @param title optional title
#' @param remNA option to keep/remove taxa that don't have common name assignment
#' @param color option to add/remove color to the graph
#'
#' @return top_taxa = data frame of top n taxonomy table
#' @return top_taxa_plot = bar graph of top n taxa
#'
#' @export

top_n_taxa <- function(physeq,
  n = 10,
  name = ShortName, # dynamically pass in the name of the common name column to use
  title = NA,
  remNA = FALSE, # option to keep/remove taxa that don't have common name assignment
  color = TRUE, # option to add/remove color to graph
  ylim100 = FALSE # option to make ylim = (0,100)
){

  ps <- physeq

  seqtab <- ps@otu_table %>%
    data.frame() %>%
    rownames_to_column(var = "sample") %>%
    pivot_longer(-sample, names_to = "asv", values_to = "count") %>%
    mutate(presence = ifelse(count > 0, 1, 0)) %>%
    group_by(asv) %>%
    summarise(prevalence = sum(presence) / n_distinct(sample) * 100, .groups = "drop")

  taxtab <- ps@tax_table %>%
    data.frame() %>%
    rownames_to_column(var = "asv")

  if(remNA) {
    taxtab <- taxtab %>%
      filter(!is.na({{name}}))
  }

  top_taxa <- taxtab %>%
    left_join(seqtab, by = "asv") %>% # add prevalence column
    slice_max(order_by = prevalence, n = n) %>%
    mutate(lowestLevel = coalesce({{name}}, species, genus, family, order, class, phylum, superkingdom),
           lowestLevel = case_when(
             is.na(lowestLevel) ~ paste0("NA", row_number()),
             TRUE ~ lowestLevel
           ),
           {{name}} := factor({{name}}, levels = unique({{name}}[order(prevalence)])),
           lowestLevel = factor(lowestLevel, levels = unique(lowestLevel[order(prevalence)])))

  top_taxa_plot <- top_taxa %>%
    ggplot(aes(x = lowestLevel, y = prevalence)) +
    # geom_bar(stat = "identity") +
    coord_flip() +
    labs(x = "",
         y = "% samples w/ taxa") +
    theme(title = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  if (ylim100) {
    top_taxa_plot <- top_taxa_plot + ylim(0,100)
  } else {
    top_taxa_plot <- top_taxa_plot
    }

  if (color) {
    top_taxa_plot <- top_taxa_plot + geom_bar(stat = "identity", aes(fill = {{name}}))
  } else {
    top_taxa_plot <- top_taxa_plot + geom_bar(stat = "identity")
  }

  top_taxa_plot <- top_taxa_plot

  if(is.na(title)) {
    top_taxa_plot <- top_taxa_plot +
      labs(title = paste0("Top ", n, " taxa"))
  } else {
    top_taxa_plot <- top_taxa_plot +
      labs(title = title)
  }

  return(list(df = top_taxa, plot = top_taxa_plot))
}
