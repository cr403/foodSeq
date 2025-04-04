library(vegan) 
library(phyloseq)
library(tidyverse) 

# Note - older versions of adonis2 (pre-August 2024) have `by="terms"` as default, but to have things separated out by term in 
# the new version, this must be input manually. 
# 
# Note - `by = "terms"` used for a sequential test of model terms where the order of terms can matter. For example `~ country + age`, 
# you could interpret as age explains x% of variation when controlled for country. `~ age + country` can be interpreted as country explains 
# x% variation when controlling for age. This is because PERMANOVA uses sequential sums of squares. So make sure to use the right order 
# depending on the question you're asking. 
#
# The other option is to use `by = "margin"` where the order of terms does not matter. For example the output of `~ country + age` could be 
# interpreted as age explains x% of variation when controlled for country while country explains y% of variance when controlled for age. In
# other words, using `by = "margin"`, predictors are tested while controlling for all other variables in the equation. 

psPerm <- function(physeq, 
                   equation, # e.g., "country + age + treatment"
                   method = "euclidean", 
                   by = "terms",
                   permutations = 1000,
                   strata = NULL,
                   seed = 123){
  # Set seed for reproducibility
  set.seed(seed)

  # Extract sam_data
  samdf <- physeq@sam_data %>%
    data.frame() 
  
  # Extract predictors from equation parameter 
  predictors <- all.vars(as.formula(paste("~", equation))) 
  
  # Ensure that all predictors exist in sample data 
  missing_predictors <- setdiff(predictors, colnames(samdf))
  
  if(length(missing_predictors) > 0) {
    stop("The following variable(s) in the formula are not found in the sample data: ",
         paste(missing_predictors, collapse = ", "))
  }
  
  # Identify rows with any missing values in predictor variables 
  missing_samples <- samdf %>% 
    select(all_of(predictors)) %>%
    apply(1, function(row) any(is.na(row)))
  
  num_removed <- sum(missing_samples) 
  
  # Drop samples with missing predictor variables 
  if(num_removed > 0) {
    warning(paste(num_removed, "samples removed due to missing values in:", 
                  paste(predictors, collapse = " or "), 
                  "\n"))
    
    # Prune those samples from the phyloseq object
    physeq <- prune_samples(!missing_samples, physeq)
  }

  # Handle strata variable  
  if(!is.null(strata)) {
    if(!strata %in% colnames(samdf)) { 
      stop("Strata variable '", strata, "' not found in sample data.")
    }
    # Ensure that all strata groups have more than 1 sample 
    strata_counts <- table(samdf[[strata]]) 
    invalid_strata <- names(strata_counts[strata_counts < 2])
    
    if(length(invalid_strata) > 0) {
      warning(length(invalid_strata), " strata group(s) had <2 samples and will be removed: ",
              paste(invalid_strata, collapse = ", "))
      
      keep_idx <- !(samdf[[strata]] %in% invalid_strata)
      
      # Prune those samples from phyloseq object 
      physeq <- prune_samples(keep_idx, physeq)
    }
  } 

  samdf <- data.frame(sample_data(physeq)) # update samdf after all pruning 
  
  # Euclidean distance object
  dist_obj <- phyloseq::distance(physeq, method = method)
  
  # Construct formula as response ~ predictors 
  form <- as.formula(paste("dist_obj ~", equation))
  
  # Run PERMANOVA
  if(!is.null(strata)) { 
    adonis2(formula = form, by = by, data = samdf, method = method, permutations = permutations, strata = samdf[[strata]])
  } else {
  adonis2(formula = form, by = by, data = samdf, method = method, permutations = permutations) 
  }
}
