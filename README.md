# Probability-scale Binomial LRT for factorial germination experiments

This repository contains the MATLAB and R implementations accompanying the manuscript
"Assessing treatment interactions on the probability scale in factorial germination experiments."

## Repository structure

- `R/` : R implementation of the proposed method
- `example/` : Example analysis using the sheepgrass data
- `data/` : Sheepgrass germination data
- `MATLAB/` : MATLAB code used for simulation studies and analysis

## R implementation

Required packages:

```r
install.packages(c("nloptr", "multcompView"))

## Basic use

```r
source("R/factorial_binomial_LRT.R")

dat <- read.csv("data/sheepgrass.csv")

fit <- germination_factorial_LRT(
  data = dat,
  alpha = 0.05,
  posthoc = TRUE,
  adjust = "holm",
  cld = TRUE
)

print(fit)
```

The input data frame should contain the following four columns:

- `Fac_A`: levels of Factor A
- `Fac_B`: levels of Factor B
- `g_ijk`: number of germinated seeds
- `n_ijk`: number of tested seeds

`Fac_A` and `Fac_B` may be numeric, character, or factor variables. Original factor labels are retained in the output.

## Post-hoc comparisons

The R implementation automatically selects the post-hoc analysis according to the interaction test.

If the interaction is not significant, pairwise comparisons are performed for significant main effects under the additive probability-scale model.

If the interaction is significant, post-hoc comparisons are performed as simple-effect tests within levels of the other factor.

P-values for multiple pairwise comparisons are adjusted using Holm's method by default. Compact letter displays are also provided.

## Sheepgrass example

The sheepgrass analysis can be reproduced by running:

```r
source("examples/example_sheepgrass.R")
```

The expected overall results are approximately:

- Germplasm: LR = 156.40, df = 7, p < 0.001
- Year: LR = 13.45, df = 1, p < 0.001
- Germplasm × Year interaction: LR = 8.65, df = 7, p = 0.279

## MATLAB implementation

The `MATLAB/` directory contains the MATLAB code used for the simulation studies and analysis reported in the manuscript.

## Notes

The usual chi-square reference distribution for the likelihood ratio statistics relies on standard regularity conditions. If a true cell probability is exactly 0 or 1, these regularity conditions may fail. An observed germination proportion of 0 or 1 in a finite sample does not, by itself, imply that the true probability is on the boundary.
