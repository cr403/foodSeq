#' @title Specific Taxa Prevalence
#'
#' @description Gives the prevalence of a particular list of taxa (asv's) in a phyloseq object. Optional faceting.
#'
#' @param physeq phyloseq object
#' @param taxList list of asv's of interest
#' @param name name of common name column in tax table of phyloseq object
#' @param facet name of samdf column by which to facet (optional)
#' @param nrow number of rows for facet_wrap()
#' @param title optional title of plot
#' @param taxListTitle optional title for taxList (y axis title)
#' @param color option for graph to be colored
#' @param localCol option for colors to be assigned locally versus globally (e.g., color by rank vs taxa)
#' @param ylim100 option for setting ylim(0,100)
#' @param labWidth option to set character breakpoint for label wrapping
#' @param titleSize set size in element_text()
#' @param textSize set size in element_text()
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
                    titleSize = 16,
                    textSize = 12
){
  ps <- physeq

  taxdf <- ps@tax_table %>%
    data.frame() %>%
    rownames_to_column(var = "asv") %>%
    filter(asv %in% taxList) %>%
    mutate(label = ifelse(is.na(.data[[name]]), lowestLevel, .data[[name]]),
           label = wrapLabels(label, labWidth))

  samdf <- ps@sam_data %>%
    data.frame() %>%
    rownames_to_column(var = "samid")

  if(!is.null(facet)) {
    samdf <- samdf %>%
      select(samid, .data[[facet]])
  } else {
    samdf <- samdf %>%
      select(samid)
  }

  seqdf <- ps@otu_table %>%
    data.frame() %>%
    rownames_to_column(var = "samid") %>%
    pivot_longer(-samid, names_to = "asv") %>%
    filter(asv %in% taxList) %>%
    left_join(samdf, by = "samid") %>%
    left_join(taxdf, by = "asv") %>%
    mutate(presence = ifelse(value > 0, 1, 0))

  taxList.prev <- seqdf %>%
    group_by(.data[[facet]], label) %>%
    summarise(prevalence = sum(presence) / n_distinct(samid) * 100) %>%
    ungroup()

  if (!is.null(facet)) {
    prev.plot <- taxList.prev %>%
      ggplot(aes(x = tidytext::reorder_within(label, by = prevalence, within = .data[[facet]]), y = prevalence)) +
      facet_wrap(~.data[[facet]], scales = "free_y", nrow = nrow) +
      tidytext::scale_x_reordered()
  } else {
    prev.plot <- taxList.prev %>%
      ggplot(aes(x = fct_reorder(label, prevalence), y = prevalence))
  }

  if (color) {
    if (localCol) {
      prev.plot <- prev.plot + geom_bar(stat = "identity", aes(fill = tidytext::reorder_within(label, by = prevalence, within = .data[[facet]]))) +
        theme(legend.position = "none")
    } else {
      prev.plot <- prev.plot + geom_bar(stat = "identity", aes(fill = label)) +
        theme(legend.position = "none")
    }

  } else {
    prev.plot <- prev.plot + geom_bar(stat = "identity")
  }

  prev.plot <- prev.plot +
    coord_flip() +
    labs(y = "Prevalence (% Samples w/ Taxa)") +
    theme(title = element_text(size = titleSize, face = "bold"),
          axis.text = element_text(size = textSize))

  if(ylim100) {
    prev.plot <- prev.plot +
      ylim(0,100)
  }

  if(!is.null(title)) {
    prev.plot <- prev.plot + labs(title = title)
  }

  if(!is.null(taxListTitle)) {
    prev.plot <- prev.plot + labs(x = taxListTitle)
  } else {
    prev.plot <- prev.plot + labs(x = "")
  }

  return(list(df = taxList.prev,
              plot = prev.plot))
}
