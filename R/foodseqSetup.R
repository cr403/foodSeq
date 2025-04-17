#' @title FoodSeq Setup Pipeline
#'
#' @description This function runs raw phyloseq object through typical inital processing steps.
#'
#' @param physeq raw phyloseq object
#' @param amplicon trnl or 12s
#' @param collapse optional to run collapseNoMismatch()
#'
#' @return ps.raw = phyloseq object with no changes from physeq other than adding p/aMR, Shannon, and reads to samdf (no NA/human reads removed)
#' @return ps = phyloseq object wiht no changes from physeq other than changes to ps.raw + animals glommed
#' @return ps.ra = changes to ps + transformed to relative abundance
#' @return ps.filt.clr = changes to ps + NA/human reads removed + clr transformation
#'
#' @export
foodseqSetup <- function(physeq,
                         amplicon = "",
                         collapse = FALSE){
  ps <- physeq
  amplicon <- tolower(amplicon) # change casing for matching

  # Initializes sample data if you're working with raw phyloseq
  if(is.null(ps@sam_data)) {
    ids <- sample_names(ps)
    samdf <- data.frame(row.names = ids) %>%
      mutate(Sample_ID = row.names(.))
    sample_data(ps) <- sample_data(samdf)
  }

  # trnL
  if(amplicon == "trnl") {

    if(collapse) {
      # CollapseNoMismatch
      seqtab.merged <- ps@otu_table

      before_collapse <- dim(seqtab.merged) # Check dimensions before and after collapse
      print(paste("Before collapse-- ", "Rows: ", before_collapse[1], "Columns: ", before_collapse[2]))

      seqtab.merged <- collapseNoMismatch(seqtab.merged)

      after_collapse <- dim(seqtab.merged)
      print(paste("After collapse-- ", "Rows: ", after_collapse[1], "Columns: ", after_collapse[2]))
      print(paste("Taxa Collapsed: ", before_collapse[2] - after_collapse[2]))

      otu_table(ps) <- otu_table(seqtab.merged, taxa_are_rows = FALSE) # Add back OTU table to ps
    }

    # Relative Abundance
    ps.ra <- transform_sample_counts(ps, function(x){x/sum(x)})

    # CLR-transform and filter
    ps.filt.clr <- ps %>%
      prune_samples(sample_sums(.) >0, .) %>% # Remove samples that do not have any reads, will mess up PCA plot
      microbiome::transform(transform = "compositional") %>%
      microbiome::transform(transform = "clr") %>%
      subset_taxa(!is.na(superkingdom)) # Remove non-foods

    # Update read counts in phyloseq object
    sample_data(ps.filt.clr)$reads <- sample_sums(ps.filt.clr)
    sample_data(ps)$reads <- sample_sums(ps)

    # pMR
    pMR <- ifelse(subset_taxa(ps.filt.clr, !is.na(superkingdom))@otu_table>0,1,0) %>%      # Converts OTU table values to presence/absence
      rowSums()
    ps.filt.clr@sam_data$pMR <- pMR

    ps@sam_data$pMR <- 0 # initialize pMR variable in ps
    matching_samples <- rownames(ps.filt.clr@sam_data) %in% rownames(ps@sam_data)
    ps@sam_data[rownames(ps.filt.clr@sam_data)[matching_samples], "pMR"] <- ps.filt.clr@sam_data$pMR

    # Shannon Index
    shannon <- estimate_richness(ps, split = TRUE, measures = "Shannon")
    ps@sam_data$Shannon_diversity_plants <- shannon$Shannon
  }

  # 12SV5
  if(amplicon == "12s") {
    # Replace "NA" strings with NA value
    tax_table(ps) <- tax_table(ps) %>%
      as.data.frame() %>%
      mutate(across(everything(), ~ ifelse(. == "NA", NA, .))) %>%
      as.matrix() %>%
      tax_table()

    # Add human_percent to sam_data for each individual
    human_asvs <- ps@tax_table %>%
      data.frame() %>%
      rownames_to_column(var = "asv") %>%
      dplyr::filter(species == "Homo sapiens") %>%
      pull(asv)

    total_reads <- ps@otu_table %>%
      data.frame() %>%
      rownames_to_column(var = "sample") %>%
      pivot_longer(-sample, names_to = "asv") %>%
      group_by(sample) %>%
      summarise(total_reads = sum(value, na.rm = TRUE))

    human_reads <- ps@otu_table %>%
      data.frame() %>%
      rownames_to_column(var = "sample") %>%
      pivot_longer(-sample, names_to = "asv") %>%
      dplyr::filter(asv %in% human_asvs) %>%
      group_by(sample) %>%
      summarise(human_reads = sum(value, na.rm = TRUE))

    human_percent <- total_reads %>%
      left_join(human_reads, by = "sample") %>%
      mutate(human_reads_percent = human_reads/total_reads*100) %>%
      select(sample, human_reads, human_reads_percent)

    sample_data(ps) <- ps@sam_data %>%
      data.frame() %>%
      rownames_to_column(var = "sample") %>%
      left_join(human_percent, by = "sample") %>%
      column_to_rownames(var = "sample") %>%
      sample_data()

    # Calculate and print % human reads
    total_human_reads <- ps %>%
      subset_taxa(species == "Homo sapiens") %>%
      sample_sums()

    total_reads <- sample_sums(ps)

    percent_human_reads <- round((sum(total_human_reads)/sum(total_reads))*100, 2)

    print(paste0(percent_human_reads, "% reads were assigned to Homo sapien"))

    # Glom animals
    tax_table(ps) <- ps@tax_table %>%
      data.frame() %>%
      mutate(lowestLevel = coalesce(species, genus, family, order, class, phylum, kingdom)) %>%
      as.matrix() %>%
      tax_table()

    ps.glom <- tax_glom(ps, taxrank = "lowestLevel")

    # Relative Abundance
    ps.ra <- transform_sample_counts(ps.glom, function(x){x/sum(x)})

    # CLR-transform and Filter NA's
    ps.filt.clr <- ps.glom %>%
      prune_samples(sample_sums(.) >0, .) %>% # Remove samples that do not have any reads, will mess up PCA plot
      microbiome::transform(transform = "compositional") %>%
      microbiome::transform(transform = "clr") %>%
      subset_taxa(!is.na(kingdom)) %>%  # Remove non-foods
      subset_taxa(lowestLevel != "Homo sapiens") # Remove human reads

    # Remove human reads
    ps.noHuman <- ps.glom %>%
      subset_taxa(lowestLevel != "Homo sapiens")

    # aMR
    aMR <- ifelse(subset_taxa(ps.noHuman, !is.na(kingdom))@otu_table>0,1,0) %>%      # Converts OTU table values to presence/absence
      rowSums()
    ps.noHuman@sam_data$aMR <- aMR

    ps@sam_data$pMR <- 0 # initialize aMR variable in ps
    matching_samples <- rownames(ps.noHuman@sam_data) %in% rownames(ps@sam_data)
    ps@sam_data[rownames(ps.noHuman@sam_data)[matching_samples], "aMR"] <- ps.noHuman@sam_data$aMR[matching_samples]

    ps.filt.clr@sam_data$aMR <- 0 # initialize aMR variable in ps.filt.clr
    matching_samples <- rownames(ps.noHuman@sam_data) %in% rownames(ps.filt.clr@sam_data)
    ps.filt.clr@sam_data[rownames(ps.noHuman@sam_data)[matching_samples], "aMR"] <- ps.noHuman@sam_data$aMR[matching_samples]

    ps.ra@sam_data$aMR <- 0 # initialize aMR variable in ps.ra
    matching_samples <- rownames(ps.noHuman@sam_data) %in% rownames(ps.ra@sam_data)
    ps.ra@sam_data[rownames(ps.noHuman@sam_data)[matching_samples], "aMR"] <- ps.noHuman@sam_data$aMR

    physeq@sam_data$aMR <- 0 # initialize aMR variable in physeq
    matching_samples <- rownames(ps.noHuman@sam_data) %in% rownames(physeq@sam_data)
    physeq@sam_data[rownames(ps.noHuman@sam_data)[matching_samples], "aMR"] <- ps.noHuman@sam_data$aMR

    # Shannon Index
    shannon <- estimate_richness(ps.noHuman, split = TRUE, measures = "Shannon")
    ps.noHuman@sam_data$Shannon_diversity_animals <- shannon$Shannon

    # Update read counts
    sample_data(ps.filt.clr)$reads <- sample_sums(ps.filt.clr)
    sample_data(ps)$reads <- sample_sums(ps)
  }

  # ### Beta diversity metrics -- This needs to be incorporated later
  # bc_dist <- phyloseq::distance(ps.ra, method = "bray") # non-clr filtered
  # bc_matrix <- as.matrix(bc_dist)

  return(list(ps.raw = physeq, # Raw ps with no changes, but p/aMR, Shannon, and reads are updated (no NA/human filtered)
              ps = ps, # Raw ps + animals glommed
              ps.ra = ps.ra, # Relative abundance
              ps.filt.clr = ps.filt.clr # Foods only, CLR-transformed
              #bc_dist = bc_dist
              #bc_matrix = bc_matrix
  ))
}
