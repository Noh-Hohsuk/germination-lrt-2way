# Factorial germination analysis on the probability and logit scales

This repository contains the MATLAB and R implementations accompanying the manuscript  
**"Assessing treatment interactions on the probability scale in factorial germination experiments."**

The proposed method provides likelihood-based tests of main effects and interaction directly on the germination-probability scale. The R implementation additionally provides conventional and Firth logistic regression, allowing factorial interaction to be analysed on either the probability or logit scale according to the scientific question.

## Repository structure

- `R/` : R implementation for probability- and logit-scale factorial analyses
- `example/` : Example analysis using the sheepgrass data
- `data/` : Sheepgrass germination data
- `MATLAB/` : MATLAB code used for the simulation studies

## R implementation

The main R function is:

```r
germination_factorial_LRT()
```

Required packages for the probability-scale analysis and compact letter displays are:

```r
install.packages(c("nloptr", "multcompView"))
```

Firth logistic regression additionally requires:

```r
install.packages("logistf")
```

Load the function using:

```r
source("R/germination_factorial_LRT.R")
```

### Input data

The input data frame should contain the following four columns:

- `Fac_A`: levels of Factor A
- `Fac_B`: levels of Factor B
- `g_ijk`: number of germinated seeds
- `n_ijk`: number of tested seeds

`Fac_A` and `Fac_B` may be numeric, character, or factor variables. Original factor labels are retained in the output.

## Probability-scale analysis

Use the constrained Binomial LRT when main effects and interaction are to be interpreted directly in terms of germination probabilities:

```r
dat <- read.csv("data/sheepgrass.csv")

fit <- germination_factorial_LRT(
  data = dat,
  scale = "probability",
  method = "mle",
  alpha = 0.05,
  posthoc = TRUE,
  adjust = "holm",
  cld = TRUE
)

print(fit)
```

On the probability scale, absence of interaction corresponds to additivity of the cell germination probabilities. Thus, treatment effects are interpreted as absolute differences in germination probability.

The additive probability-scale model is fitted by constrained maximum likelihood estimation so that all fitted probabilities remain between 0 and 1.

## Logit-scale analysis

Conventional binomial logistic regression is available by specifying:

```r
fit_logit <- germination_factorial_LRT(
  data = dat,
  scale = "logit",
  method = "mle",
  alpha = 0.05,
  posthoc = TRUE,
  adjust = "holm",
  cld = TRUE
)

print(fit_logit)
```

For datasets in which ordinary maximum likelihood estimation is affected by separation or extreme estimates, Firth logistic regression is available using:

```r
fit_firth <- germination_factorial_LRT(
  data = dat,
  scale = "logit",
  method = "firth",
  alpha = 0.05,
  posthoc = TRUE,
  adjust = "holm",
  cld = TRUE
)

print(fit_firth)
```

The probability- and logit-scale analyses test different definitions of interaction. Absence of interaction on the probability scale means additivity of germination probabilities, whereas absence of interaction in a conventional logistic model means additivity on the log-odds scale. The choice of scale should therefore be guided by the scientific question.

Ordinary and Firth logistic regression share the same logit-scale interaction hypothesis but differ in estimation.

## Post-hoc comparisons

The R implementation links the overall factorial tests to follow-up comparisons.

When the interaction is not significant, pairwise comparisons are available for factors with significant overall effects. When the interaction is significant, simple-effect comparisons can be performed within levels of the other factor.

P-values for multiple pairwise comparisons are adjusted using Holm's method by default. Compact letter displays are also provided.

For probability-scale analyses, contrasts are expressed as differences in fitted germination probabilities.

For logit-scale analyses, contrasts are tested on the log-odds scale and are reported together with odds ratios. Corresponding differences in fitted germination probabilities are also provided to aid interpretation.

## Sheepgrass example

The probability-scale analysis of the sheepgrass data reported in the manuscript can be reproduced by running:

```r
source("example/example_sheepgrass.R")
```

The expected overall results are approximately:

- Germplasm: LR = 156.40, df = 7, p < 0.001
- Year: LR = 13.45, df = 1, p < 0.001
- Germplasm × Year interaction: LR = 8.65, df = 7, p = 0.279

Minor numerical differences may occur across platforms or optimizer versions.

## MATLAB implementation

The `MATLAB/` directory contains the MATLAB code used for the simulation studies reported in the manuscript.

## Statistical note on boundary probabilities

The usual chi-square reference distribution for likelihood ratio statistics relies on standard regularity conditions. If a true cell germination probability is exactly 0 or 1, these regularity conditions may fail.

This situation should be distinguished from observing a germination proportion of 0 or 1 in a finite sample when the underlying probability lies strictly between 0 and 1. The boundary simulations in the manuscript examine the empirical behaviour of the proposed procedure under the configurations considered; they do not establish general validity of the usual chi-square approximation for true boundary probabilities.
