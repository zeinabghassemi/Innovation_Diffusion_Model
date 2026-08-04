# Innovation Diffusion Model

MATLAB code for reproducing the simulations presented in the manuscript on innovation diffusion, adoption dynamics, and inequality between heterogeneous communities.

The repository implements an innovation diffusion model that integrates:
- social learning and intergroup interactions,
- income-based heterogeneity between communities,
- experience-based cost reduction (learning effects),
- affordability constraints,
- adoption disparities between high-income and low-income groups.

The code reproduces the main-text figures and supplementary analyses, including the effects of:
- intergroup interaction on adoption and disparity,
- learning rate improvements,
- income distribution assumptions,
- population size,
- community composition,
- model parameter sensitivity.

---

## Model Overview

The model considers two interacting communities:

- High-income community ($c_H$)
- Low-income community ($c_L$)

Each community has its own adoption dynamics governed by:
- social influence,
- external adoption effects,
- affordability constraints,
- innovation cost evolution.

The model tracks:

- overall adoption proportion,
- community-level adoption,
- transient disparity,
- equilibrium disparity,
- cumulative disparity.

The innovation cost decreases through experience-based learning, where higher cumulative adoption accelerates future cost reductions.

## Requirements

The code was developed using MATLAB.

Required MATLAB functionality:

- MATLAB Optimization Toolbox  
  (for `fsolve`)
- MATLAB Statistics and Machine Learning Toolbox  
  (for Gamma distribution functions such as `gamcdf`)
- MATLAB ODE solvers  
  (for numerical integration using `ode45`)

