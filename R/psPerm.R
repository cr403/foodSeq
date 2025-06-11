#' @title Phyloseq PERMANOVA
#'
#' @description function to run permanova on foodseq phyloseq object
#'
#' @param physeq filtered and clr transfomed phyloseq object
#' @param equation e.g., "country + age + treatment" or "country*age + treatment"
#' @param method euclidean
#' @param by "terms" for sequential assessment/control for covariates in order, "margins" for individual assessment/control for all other variables
#' @param permutations number of permutations
#' @param strata grouping variable
#' @param seed random seed for reproducibility
#'
#' @return PERMANOVA output
#'
#' @export
psPerm <- function(physeq,
                   equation, # e.g., "country + age + treatment" or "country*age + treatment"
                   method = "euclidean",
                   by = "terms",
                   permutations = 1000,
                   strata = NULL,
                   seed = 123){

  # Set seed for reproducibility
  set.seed(seed)

  # Extract sample data
  samdf <- data.frame(sample_data(physeq))
  message("Initial number of samples: ", nsamples(physeq))

  ################### Handle Predictor Variables ###################

  # Extract predictors from equation
  predictors <- all.vars(as.formula(paste("~", equation)))

  # Ensure predictors exist
  missing_predictors <- setdiff(predictors, colnames(samdf))
  if (length(missing_predictors) > 0) {
    stop("The following variable(s) in the formula are not found in the sample data: ",
         paste(missing_predictors, collapse = ", "))
  }

  # Check for missing predictor values
  missing_samples <- samdf %>%
    dplyr::select(all_of(predictors)) %>%
    apply(1, function(row) any(is.na(row)))

  num_removed <- sum(missing_samples)

  if (num_removed > 0) {
    removed_ids <- rownames(samdf)[missing_samples]
    warning(
      paste(
        num_removed, "samples removed due to missing values in",
        paste(predictors, collapse = " or "), "\n",
        "Samples removed:", paste(removed_ids, collapse = ", "), "\n"
      )
    )

    # Prune those samples from the phyloseq object
    physeq <- prune_samples(!missing_samples, physeq)
  }

  # Update samdf after pruning
  samdf <- data.frame(sample_data(physeq))
  message("Number of samples after missing predictor pruning: ", nsamples(physeq))

  ################### Handle strata variable ###################

  if (!is.null(strata)) {
    if (!strata %in% colnames(samdf)) {
      stop("Strata variable '", strata, "' not found in sample data.")
    }

    # Remove samples with NA in strata
    missing_strata <- is.na(samdf[[strata]])
    num_removed_strata_na <- sum(missing_strata)
    if (num_removed_strata_na > 0) {
      removed_ids <- rownames(samdf)[missing_strata]
      warning(
        paste(
          num_removed_strata_na, "samples removed due to missing values in", strata, "\n",
          "Samples removed:", paste(removed_ids, collapse = ", "), "\n"
        )
      )
      physeq <- prune_samples(!missing_strata, physeq)
      samdf <- data.frame(sample_data(physeq))
    }

    # Remove strata groups with < 2 samples
    strata_counts <- table(samdf[[strata]])
    invalid_strata <- names(strata_counts[strata_counts < 2])
    if (length(invalid_strata) > 0) {
      # Which samples will be removed?
      remove_idx <- samdf[[strata]] %in% invalid_strata
      removed_ids <- rownames(samdf)[remove_idx]

      warning(
        paste(
          length(invalid_strata), "sample(s) had <2 samples within strata and will be removed:",
          paste(invalid_strata, collapse = ", "), "\n",
          "Samples removed:", paste(removed_ids, collapse = ", "), "\n"
        )
      )
      keep_idx <- !remove_idx
      physeq <- prune_samples(keep_idx, physeq)
      samdf <- data.frame(sample_data(physeq))
    }
  }

  message("Number of samples after strata pruning: ", nsamples(physeq))

  ################### Run PERMANOVA ###################

  # Build distance object
  dist_obj <- phyloseq::distance(physeq, method = method)

  # Build formula: distance ~ predictors
  form <- as.formula(paste("dist_obj ~", equation))
  if(!is.null(strata)) {
    final_counts <- table(samdf[[strata]])
    if(any(final_counts < 2)) {
      stop("At least one stratum has <2 samples after all pruning. Cannot run PERMANOVA.")
    }

    if(nrow(samdf) < 2) {
      stop("Fewer than 2 samples remain after pruning. Cannot run PERMANOVA")
    }
    vegan::adonis2(formula = form, by = by, data = samdf, method = method, permutations = permutations, strata = samdf[[strata]])
  } else {
    vegan::adonis2(formula = form, by = by, data = samdf, method = method, permutations = permutations)
  }
}
