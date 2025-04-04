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
                   equation, # e.g., "country + age + treatment" or "country*age + bmi"
                   method = "euclidean", 
                   by = "terms",
                   permutations = 1000){

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
    samdf <- samdf[!missing_samples, , drop = FALSE]
  }
  
  # Euclidean distance matrix
  trnL_dist_run <- as.data.frame(as.matrix(phyloseq::distance(physeq, method = method)))
  
  # Construct formula as response ~ predictors 
  form <- as.formula(paste("trnL_dist_run ~", equation)) 
  
  # Run PERMANOVA 
  adonis2(form, by = by, data = samdf, method = method, permutations = permutations) 
}
