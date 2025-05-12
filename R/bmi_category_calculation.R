# This is just being kept in here to document how i calculated BMI for the SAGIP cohort
# cdc_bmi <- read_csv('~/Library/CloudStorage/Box-Box/project_davidlab/LAD_LAB_Personnel/Caroline_R/4_Templates/CDC_age2-19_BMIpercentiles.csv') %>%
#   data.frame() %>%
#   mutate(across(everything(), as.numeric))
#
# cdc_colnames <- colnames(cdc_bmi)
#
#
#
#
#
#
#
#
# calculate_bmi_zscore <- function(bmi, L, M, S) {
#   if (any(is.na(c(bmi, L, M, S)))) {
#     return(NA_real_)
#   }
#   if (L == 0) {
#     z <- log(bmi / M) / S
#   } else {
#     z <- ((bmi / M)^L - 1) / (L * S)
#   }
#   return(z)
# }
#
# # Function to convert z-score to percentile
# z_to_percentile <- function(z) {
#   pnorm(z) * 100
# }
#
#
#
#
#
#
#
#
# sample_data(ps) <- ps@sam_data %>%
#   data.frame() %>%
#   rownames_to_column(var = "samid") %>%
#   mutate(age_months = age*12 + 0.5, # Add 0.5 to match CDC data
#          sex_bmi = ifelse(gender == "Male", 1, 2),
#          bmi = as.numeric(bmi))  %>% # NOTE: a few samples were missing gender, get automatically assigned to female
#   left_join(cdc_bmi, by = c("age_months" = "Agemos", "sex_bmi" = "Sex")) %>% # age-sex matching
#   rowwise() %>%
#   mutate(
#     bmi_z = calculate_bmi_zscore(bmi, L, M, S),
#     bmi_percentile = z_to_percentile(bmi_z),
#     bmi_category = case_when(
#       bmi_percentile < 5 ~ "underweight",
#       bmi_percentile >= 5 & bmi_percentile < 85 ~ "healthy weight",
#       bmi_percentile >= 85 & bmi_percentile < 95 ~ "overweight",
#       bmi_percentile >= 95 ~ "obesity"
#     )
#   )  %>%
#   ungroup() %>%
#   mutate(bmi_category = factor(bmi_category, levels = c("healthy weight", "underweight", "overweight", "obesity"))) %>%
#   select(-any_of(cdc_colnames)) %>%
#   select(-c("age_months", "sex_bmi", "bmi_z")) %>%
#   column_to_rownames(var = "samid") %>%
#   sample_data()
