# LegATo: Longitudinal mEtaGenomic Analysis Toolkit <img src="https://github.com/aubreyodom/Legato-docs/blob/main/legato-logo.png?raw=true" align="right" width="170" />

## What is LegATo?

LegATo is a suite of open-source software tools for longitudinal microbiome analysis. It is extendable to
several different study forms with optimal ease-of-use for researchers. Microbiome time-series data
presents distinct challenges including complex covariate dependencies and variety of longitudinal study
designs. This toolkit will allow researchers to determine which microbial taxa are affected over time by
perturbations such as onset of disease or lifestyle choices, and to predict the effects of these perturbations
over time, including changes in composition or stability of commensal bacteria. 

LegATo integrates visualization, modeling and testing procedures. It is currently in development, but it will soon be supplemented by hierarchical clustering tools and multivariate generalized estimating equations (JGEEs) to adjust for the compositional nature of microbiome data. Other tools will be implemented as needed.

# Documentation
Documentation and tutorials for MetaScope are available at our [website](https://github.com/aubreyodom/Legato).

# Installation
LegATo requires R Version 4.3.

Install the development version of the package from Github:

```
if (!requireNamespace("devtools", quietly=TRUE))
  install.packages("devtools")
devtools::install_github("aubreyodom/LegATo")
```
