# example_sheepgrass.R
# ------------------------------------------------------------
# Example analysis of the sheepgrass germination data using the
# probability-scale two-factor Binomial LRT.
#
# Files expected in the same folder:
#   factorial_binomial_LRT.R
#   sheepgrass.csv
#
# Required packages:
#   nloptr
#   multcompView
# ------------------------------------------------------------

# Uncomment once if the packages are not installed:
# install.packages(c("nloptr", "multcompView"))

source("factorial_binomial_LRT.R")

# Read data.
# Required columns:
#   Fac_A  : Factor A (germplasm here)
#   Fac_B  : Factor B (year here)
#   g_ijk  : number of germinated seeds
#   n_ijk  : number of tested seeds
#
# The function accepts numeric, character, or factor labels directly.
sheepgrass <- read.csv("sheepgrass.csv", stringsAsFactors = FALSE)

# Inspect the first few observations.
head(sheepgrass)

# Run the overall probability-scale Binomial LRT.
# Because the germplasm-by-year interaction is not significant in this
# dataset, the function automatically performs pairwise comparisons
# for the significant main effects under the additive probability-scale
# model, followed by Holm adjustment and compact letter displays.
fit <- germination_factorial_LRT(
  data = sheepgrass,
  alpha = 0.05,
  posthoc = TRUE,
  adjust = "holm",
  cld = TRUE
)

# Main printed output:
print(fit)

# ------------------------------------------------------------
# Extract specific results
# ------------------------------------------------------------

# Overall tests for germplasm, year, and interaction
fit$overall_tests

# Maximized log-likelihoods for the fitted models
fit$logLik

# Fitted probabilities under the additive model
fit$fitted_probabilities$M_AplusB

# Germplasm pairwise comparisons (Holm adjusted)
fit$posthoc$Factor_A

# Year pairwise comparison
fit$posthoc$Factor_B

# Compact letter displays
fit$cld$Factor_A
fit$cld$Factor_B

# Numerical diagnostics for constrained estimation
fit$additive_fit$diagnostics


# ------------------------------------------------------------
# Optional: export results
# ------------------------------------------------------------

write.csv(
  fit$overall_tests,
  "sheepgrass_overall_tests.csv",
  row.names = FALSE
)

if (!is.null(fit$posthoc$Factor_A)) {
  write.csv(
    fit$posthoc$Factor_A,
    "sheepgrass_posthoc_germplasm.csv",
    row.names = FALSE
  )
}

if (!is.null(fit$posthoc$Factor_B)) {
  write.csv(
    fit$posthoc$Factor_B,
    "sheepgrass_posthoc_year.csv",
    row.names = FALSE
  )
}

if (!is.null(fit$cld$Factor_A)) {
  write.csv(
    fit$cld$Factor_A,
    "sheepgrass_cld_germplasm.csv",
    row.names = FALSE
  )
}

if (!is.null(fit$cld$Factor_B)) {
  write.csv(
    fit$cld$Factor_B,
    "sheepgrass_cld_year.csv",
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# Expected overall test results (minor numerical differences
# may occur across platforms/optimizer versions)
# ------------------------------------------------------------
#
# Germplasm:
#   LR = 156.40, df = 7, p < 0.001
#
# Year:
#   LR = 13.45, df = 1, p < 0.001
#
# Germplasm x Year interaction:
#   LR = 8.65, df = 7, p = 0.279
#
# The fitted probability difference for 2016 - 2017 is
# approximately -0.046 under the additive model.
