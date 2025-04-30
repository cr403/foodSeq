#' @title Top Taxa
#'
#' @description finds and plots top taxa from a phyloseq
#'
#' @param physeq phyloseq object
#' @param n number of taxa to plot
#' @param name name of common name column (e.g., ShortName)
#' @param labWidth character length for string breakpoints (used for long labels)
#' @param title optional title
#' @param remNA option to keep/remove taxa that don't have common name assignment
#' @param color option to add/remove color to the graph
#' @param facet option to add facet_wrap()
#' @param colorGlobal option to make colors across facets the same (e.g., Poaceae always same color as Poaceae)
#'
#' @return top_taxa = data frame of top n taxonomy table
#' @return top_taxa_plot = bar graph of top n taxa
#'
#' @export

top_n_taxa <- function(physeq,
  n = 10,
  name = "ShortName", # dynamically pass in the name of the common name column to use
  labWidth = 60, # character length for string breakpoints
  title = NA,
  remNA = FALSE, # option to keep/remove taxa that don't have common name assignment
  color = TRUE, # option to add/remove color to graph
  colorGlobal = FALSE, # option to make colors across facets the same (e.g., Poaceae always same color as Poaceae)
  ylim100 = FALSE, # option to make ylim = (0,100)
  facet = NULL, # optional facet_wrap
  nrow = NULL # optional facet_wrap nrow
){

  ps <- physeq

  # Convert OTU table to presence/absence
  seqtab.pa <- ps@otu_table %>%
    data.frame() %>%
    rownames_to_column(var = "samid") %>%
    pivot_longer(-samid, names_to = "asv") %>%
    mutate(presence = ifelse(value > 0, 1, 0))

  if (!is.null(facet)) {
    # Extract facet variable
    samdf <- ps@sam_data %>%
      data.frame() %>%
      rownames_to_column(var = "samid") %>%
      select(samid, .data[[facet]])

    # Group to facet variable and calculate prevalence
    seqtab <- seqtab.pa %>%
      left_join(samdf, by = "samid") %>%
      group_by(asv, .data[[facet]]) %>%
      summarise(prevalence = sum(presence) / n_distinct(samid) * 100, .groups = "drop")

  } else {
    # Calculate prevalence globally
    seqtab <- seqtab.pa
      group_by(asv) %>%
      summarise(prevalence = sum(presence) / n_distinct(sample) * 100, .groups = "drop")
  }

  # Extract tax table
  taxtab <- ps@tax_table %>%
    data.frame() %>%
    rownames_to_column(var = "asv")

  # Remove NA's if desired
  if(remNA) {
    taxtab <- taxtab %>%
      filter(!is.na(.data[[name]]))
  }

  # Add taxonomy labels
  top_taxa <- top_taxa %>%
    left_join(seqtab, by = "asv") %>%
    arrange(across(any_of(facet)), desc(prevalence), label) %>%
    mutate(label = coalesce(.data[[name]], species, genus, family, order, class, phylum, superkingdom),
           label = case_when(
             is.na(label) ~ paste0("NA", row_number()),
             TRUE ~ label
           ),
           label = wrapLabels(label, width = labWidth),
           label = factor(label)) %>%
    group_by(across(any_of(facet))) %>%
    arrange(desc(prevalence), .by_group = TRUE) %>%
    slice_head(n = n) %>%
    mutate(rank = factor(row_number()))

  # Plot
  if (!is.null(facet)) {
    top_taxa_plot <- top_taxa %>%
      ggplot(aes(x = tidytext::reorder_within(label, prevalence, .data[[facet]]), y = prevalence)) +
      facet_wrap(~.data[[facet]], scales = "free_y", nrow = nrow) +
      tidytext::scale_x_reordered()
  } else {
    top_taxa_plot <- top_taxa %>%
      ggplot(aes(x = fct_reorder(label, prevalence), y = prevalence))
  }

  top_taxa_plot <- top_taxa_plot +
    coord_flip() +
    labs(x = "",
         y = "% samples w/ taxa") +
    theme(title = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 12),
          legend.position = "none")

  if (ylim100) {
    top_taxa_plot <- top_taxa_plot + ylim(0,100)
  } else {
    top_taxa_plot <- top_taxa_plot
    }

  if (color) {
    if (!is.null(facet)) {
      if (colorGlobal) {
        top_taxa_plot <- top_taxa_plot + geom_bar(stat = "identity", aes(fill = rank))
      } else {
        top_taxa_plot <- top_taxa_plot + geom_bar(stat = "identity", aes(fill = tidytext::reorder_within(label, prevalence, .data[[facet]])))
      }
    } else {
      top_taxa_plot <- top_taxa_plot + geom_bar(stat = "identity", aes(fill = fct_reorder(label, prevalence)))
    }
  } else {
    top_taxa_plot <- top_taxa_plot + geom_bar(stat = "identity")
  }

  if(is.na(title)) {
    top_taxa_plot <- top_taxa_plot +
      labs(title = paste0("Top ", n, " taxa"))
  } else {
    top_taxa_plot <- top_taxa_plot +
      labs(title = title)
  }

  return(list(df = top_taxa, plot = top_taxa_plot))
}
