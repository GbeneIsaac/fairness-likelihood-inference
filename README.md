# fairness-likelihood-inference
Simulation code for evaluating fairness in likelihood-based inference systems using RMEP/RMED metrics under hierarchical mixture models.

# Fairness in Likelihood-Based Inference Systems

This repository contains simulation and analysis code supporting the paper:

**“Statistical Analysis of Fairness in Likelihood-Based Inference Systems”**
submitted to *The American Statistician*.
*Authors: Isaac Gbene, Dylan Borchert, Andrew Simpson,
Christopher Saunders, and Semhar Michael∗
Department of Mathematics and Statistics, South Dakota State University


## Overview

The code implements Monte Carlo simulation studies to evaluate fairness properties of likelihood-based inference systems, with a focus on likelihood ratio (LR) performance across subpopulations. Fairness is assessed using the **Rates of Misleading Evidence in favor of the Prosecution (RMEP)** and **Defense (RMED)** under hierarchical mixture models.

The framework allows systematic investigation of:
- Subpopulation imbalance via mixing proportions (π)
- Dimensionality effects (p = 1, 2, 5)
- Within- and between-source variability
- Density-based vs parametric likelihood ratios
- Relative changes in RMEP across groups

## Repository Structure
├── 1DSIM.R # Univariate simulations (p = 1)
├── 2DSIM.R # Bivariate simulations (p = 2)
├── 5DSIM.R # Multivariate simulations (p = 5)
├── UPDATEDRMEPMTV.R # Core RMEP/RMED computation functions
├── Copper_Wire_Simulation.R # Hierarchical copper wire simulations
├── GlassDataAnalysis.R # Glass evidence analysis
├── hier_mix_mod_bivariate_scatterplots_copperwire.R
└── Results/ # Simulation outputs (generated)


## Methods

Simulations are based on hierarchical mixture models with:
- Between-source variability modeled via multivariate normal mixtures
- Within-source variability modeled using conditional normal sampling
- Likelihood ratios computed using the `comparison` package

Fairness is quantified through:
- Group-conditional RMEP and RMED
- Relative change metrics between subpopulations
- Monte Carlo uncertainty via standard deviation and standard error estimates

## Reproducibility

All simulations are fully reproducible using fixed random seeds. Parameter settings (π, p, B, m, nᵢ) correspond directly to those reported in the manuscript.

## Software Requirements

- R (≥ 4.0)
- Packages: `MixSim`, `mvtnorm`, `comparison`, `ggplot2`

## Citation

If you use this code, please cite the associated manuscript (citation to be added upon publication).

