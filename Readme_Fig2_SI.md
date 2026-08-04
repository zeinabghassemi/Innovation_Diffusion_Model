## Supplementary Figure S2: Income distribution sensitivity analysis

### Script

`Figure2_SI_IncomeDistribution.m`

### Description

This script reproduces Supplementary Figure S2 and evaluates the robustness of the model outcomes under different assumptions about income distributions.

The analysis considers two aspects of income heterogeneity:

1. Changes in the relative mean income between high-income and low-income communities.
2. Changes in the variance of income distributions.

For each case, the model evaluates the impact on:

- equilibrium disparity between communities;
- overall adoption proportion.

The simulations are repeated under different levels of intergroup interaction and learning rates.

### Parameters explored

#### Mean income ratio sensitivity

The ratio between community income parameters is varied while evaluating multiple levels of homophily:

- h = 1
- h = 0.875
- h = 0.625
- h = 0.5

#### Variance sensitivity

The variance of the income distribution is modified using the Gamma distribution parameters.

#### Learning rate sensitivity

Learning parameters:

- E = 0.16
- E = 0.18
- E = 0.20

are evaluated under different income distribution assumptions.

### Outputs

The script generates figures showing:

- normalized equilibrium disparity;
- normalized overall adoption proportion.

### Running the code

Run:
Figure2_SI_IncomeDistribution.m
