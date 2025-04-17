#' @title Parallel Analysis
#'
#' @description Runs parallel analysis (permutation) for principal component importance
#'
#' @param physe filtered and clr transformed phyloseq object on which pcaPlot() was run
#' @param pcaOutput object storing the output from pcaPlot()
#' @param B number of permutations
#' @param centile percentile threshold
#' @param seed random seed for reproducibility
#'
#' @return cutoff_pc = last significantly important pc
#' @return var_explained = variance explained from PC1-cutoff_pc
#' @return paran.df = data frame used to generate plot
#' @return paran.plot = parallel analysis plot (scree table with cutoff)
#' @return qc.df = data frame with prcomp() numbers from internal and external (pcaPlot()) pca generation
#' @return qc.plot = scatterplot of eigenvalues from internal and external (pcaPlot()) pca generation
#'
#' @export
paranPC <- function(physeq, # Filtered and CLR transformed phyloseq
                    pcaOutput, # object storing output of pcaPlot()
                    B = 1000,
                    centile = 95,
                    seed = 123
){
  # setup parallelization
  future::plan("multisession", workers = parallel::detectCores() - 1) # permutation step will use all computer cores available except 1
  set.seed(seed)

  # extract structure information
  seqtab <- physeq@otu_table %>%
    as.matrix

  n <- nrow(seqtab)
  p <- ncol(seqtab)

  # compute observed eigenvalues
  obs_pca <- prcomp(seqtab, center = TRUE, scale. = FALSE)
  obs_eig <- obs_pca$sdev^2

  # run B permutations in parallel
  rand_eig <- future.apply::future_sapply(1:B, function(b) {
    perm_data <- apply(seqtab, 2, sample)  # column-wise permutation
    eigen(cov(perm_data))$values
  })
  rand_eig <- t(rand_eig)  # rows = permutations, cols = PCs

  # compute quantiles
  rand_means <- apply(rand_eig, 2, mean)
  rand_percentile <- apply(rand_eig, 2, function(x) quantile(x, probs = centile/100))

  # find cutoff
  cutoff_pc <- sum(obs_eig > rand_percentile)
  variance_explained <- sum(obs_eig[1:cutoff_pc]) / sum(obs_eig)

  # plot
  paran.df <- data.frame(
    PC = paste0("PC", 1:length(obs_eig)),
    PC_index = 1:length(obs_eig),
    observed_eig = obs_eig,
    rand_mean = rand_means,
    rand_percentile = rand_percentile
  ) %>%
    pivot_longer(cols = c(observed_eig, rand_mean, rand_percentile),
                 names_to = "type", values_to = "eigenvalue")

  paran.plot <- paran.df %>%
    ggplot(aes(x = PC_index, y = eigenvalue, color = type)) +
    geom_line(size = 1) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_vline(xintercept = cutoff_pc, color = "red", linetype = "dashed") +
    annotate("text", x = cutoff_pc + 0.5, y = max(obs_eig),
             label = paste0("Retain up to PC ", cutoff_pc,
                            " (", round(variance_explained*100,1), "% variance)"),
             hjust=0, vjust=1.2, color="red") +
    labs(title="Parallel Analysis (Permutation)",
         x="Principal Component",
         y="Eigenvalue")

  # ensure that eigenvalues generated here are the same as those generated in pcaPlot()
  qc.df <- paran.df %>%
    pivot_wider(id_cols= PC, names_from = "type", values_from = "eigenvalue") %>%
    left_join(pcaOutput$scree.table, by = "PC")

  qc.plot <- qc.df %>%
    ggplot(aes(x = observed_eig, y = Eigenvalue)) +
    geom_point() +
    labs(title = "QC: pcaPlot() versus paranPC() Eigenvalues")

  # Output
  print(paran.plot)

  return(list(
    cutoff_pc = cutoff_pc,
    var_explained = variance_explained,
    paran.df = paran.df,
    paran.plot = paran.plot,
    qc.df = qc.df,
    qc.plot = qc.plot
  ))
}
