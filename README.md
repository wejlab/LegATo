
# LegATo
### A Longitudinal mEtaGenomic Analysis Toolkit <img src="https://github.com/aubreyodom/Legato-docs/blob/main/legato-logo.jpg?raw=true" align="right" width="140">

LegATo is a suite of open-source software tools for longitudinal microbiome analysis. It is extendable to
several different study forms with optimal ease-of-use for researchers. Microbiome time-series data
presents distinct challenges including complex covariate dependencies and variety of longitudinal study
designs. This toolkit will allow researchers to determine which microbial taxa are affected over time by
perturbations such as onset of disease or lifestyle choices, and to predict the effects of these perturbations
over time, including changes in composition or stability of commensal bacteria. 

LegATo integrates visualization, modeling and testing procedures. It is currently in development, but it will soon be supplemented by hierarchical clustering tools and multivariate generalized estimating equations (JGEEs) to adjust for the compositional nature of microbiome data.

### The Story Behind the Name
In music, legato indicates that notes are played or sung smoothly and connected, without a noticeable break between them. The LegATo package facilitates a cohesive and interconnected understanding of the microbial communities represented by the samples, much like the smooth connection of musical notes in a legato passage. 

Therefore, LegATo metaphorically represents the smooth and connected analysis of longitudinal metagenomic data, drawing inspiration from the musical term to convey a sense of continuity and harmony in the modeling process.

## Documentation
Documentation and tutorials for LegATo are available at our [website](https://aubreyodom.github.io/LegATo-docs/).

Check out a thorough tutorial on proper usage of our package [here](https://aubreyodom.github.io/LegATo-docs/articles/LegATo_vignette.html).

## Installation
LegATo requires R Version 4.3.

Install the development version of the package from Github:

```
if (!requireNamespace("devtools", quietly=TRUE))
  install.packages("devtools")
devtools::install_github("aubreyodom/LegATo")
```
