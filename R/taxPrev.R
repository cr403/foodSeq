#' @title Specific Taxa Prevalence
#'
#' @description Gives the prevalence of a particular list of taxa (asv's) in a phyloseq object. Optional faceting.
#'
#' @param physeq phyloseq object
#' @param taxList list of asv's of interest
#' @param name name of common name column in tax table of phyloseq object
#' @param facet name of samdf column by which to facet (optional)
#' @param stripText optional labels for strip text
#' @param nrow number of rows for facet_wrap()
#' @param title optional title of plot
#' @param taxListTitle optional title for taxList (y axis title)
#' @param color option for graph to be colored
#' @param localCol option for colors to be assigned locally versus globally (e.g., color by rank vs taxa)
#' @param ylim100 option for setting ylim(0,100)
#' @param labWidth option to set character breakpoint for label wrapping
#' @param numLab option to add value labels to bar chart
#' @param titleSize set plot title size
#' @param axtitleSize set axis title size
#' @param textSize set plot text size
#' @param stripSize set strip text size
#'
#' @return df  = data frame of taxa and prevalence
#' @return plot  = bar chart of taxa and prevalence
#'
#' @export
taxPrev <- function(physeq, # phyloseq object
                    taxList, # ASV's of interest
                    name = "ShortName", # Common name column in tax table
                    facet = NULL, # optional facet variable
                    nrow = NULL, # nrow for facet_wrap()
                    title = NULL, # optional title
                    taxListTitle = NULL, # optional axis title
                    color = TRUE, # optional colored bars
                    localCol = TRUE, # option for local vs global colors
                    ylim100 = TRUE, # option to set ylim(0,100) vs letting ggplot decide
                    labWidth = 60, # can edit wrap length
                    numLab = FALSE, # option to add value labels to bar chart
                    titleSize = 16, # title text size
                    axtitleSize = 16, # axis title size
                    textSize = 12, # text size
                    stripSize = 10 # strip text size
){
  ps <- physeq

  # Extract tax table
  taxdf <- ps@tax_table %>%
    data.frame() %>%
    rownames_to_column(var = "asv") %>%
    filter(asv %in% taxList) %>%
    mutate(label = ifelse(is.na(.data[[name]]), lowestLevel, .data[[name]]), # build label for plotting
           label = wrapLabels(label, labWidth))

  # Extract sample data
  samdf <- ps@sam_data %>%
    data.frame() %>%
    rownames_to_column(var = "samid")

  # Take out facet variable from sample data
  if(!is.null(facet)) {
    samdf <- samdf %>%
      select(samid, .data[[facet]])
  } else {
    samdf <- samdf %>%
      select(samid)
  }

  # Extract otu table and join with tax/samdf tables
  seqdf <- ps@otu_table %>%
    data.frame() %>%
    rownames_to_column(var = "samid") %>%
    pivot_longer(-samid, names_to = "asv") %>%
    filter(asv %in% taxList) %>%
    left_join(samdf, by = "samid") %>%
    left_join(taxdf, by = "asv") %>%
    mutate(presence = ifelse(value > 0, 1, 0)) # Convert to presence/absence


  # Build prevalence data frame within facet
  if (!is.null(facet)) {
    taxList.prev <- seqdf %>%
      group_by(.data[[facet]], label) %>%
      summarise(prevalence = sum(presence) / n_distinct(samid) * 100) %>%
      ungroup()

    prev.plot <- taxList.prev %>%
      ggplot(aes(x = tidytext::reorder_within(label, by = prevalence, within = .data[[facet]]), y = prevalence)) +
      tidytext::scale_x_reordered() +
      facet_wrap(~.data[[facet]], scales = "free_y", nrow = nrow)

  # Build prevalence data frame within entire phyloseq
  } else {
    taxList.prev <- seqdf %>%
      group_by(label) %>%
      summarise(prevlaence = sum(presence) / n_distinct(samid) * 100) %>%
      ungroup()

    prev.plot <- taxList.prev %>%
      ggplot(aes(x = fct_reorder(label, prevalence), y = prevalence))
  }

  # Add bars with color to graph
  if (color) {

    # Local coloring
    if (localCol) {
      prev.plot <- prev.plot + geom_bar(stat = "identity", aes(fill = tidytext::reorder_within(label, by = prevalence, within = .data[[facet]]))) +
        theme(legend.position = "none")

    # Global coloring
    } else {
      prev.plot <- prev.plot + geom_bar(stat = "identity", aes(fill = label)) +
        theme(legend.position = "none")
    }

  # Add bars without color to graph
  } else {
    prev.plot <- prev.plot + geom_bar(stat = "identity")
  }

  # Build remaining parts of graph
  prev.plot <- prev.plot +
    coord_flip() +
    labs(y = "Prevalence (% Samples w/ Taxa)") +
    theme(plot.title = element_text(size = titleSize, face = "bold"),
          axis.title = element_text(size = axtitleSize, face = "bold"),
          axis.text = element_text(size = textSize),
          strip.text = element_text(size = stripSize))

  # Set axis limits
  if(ylim100) {
    prev.plot <- prev.plot +
      ylim(0,100)
  }

  # Add title (default is no title)
  if(!is.null(title)) {
    prev.plot <- prev.plot + labs(title = title)
  }

  # Add axis title (default is no title)
  if(!is.null(taxListTitle)) {
    prev.plot <- prev.plot + labs(x = taxListTitle)
  } else {
    prev.plot <- prev.plot + labs(x = "")
  }

  # Add number labels (default is no label)
  if(numLab) {
    prev.plot <- prev.plot +
      geom_text(aes(label = round(prevalence, 1)), hjust = -0.5) +
      expand_limits(y = 105)
  }

  return(list(df = taxList.prev,
              plot = prev.plot))
}
