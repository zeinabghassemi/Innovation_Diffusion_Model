# Innovation Diffusion Model with Income Disparities

This repository contains MATLAB codes used to reproduce the simulation results presented in the manuscript.

The model investigates how technology learning dynamics influence innovation adoption and disparities between high-income and low-income communities.

## Repository contents

### Figure2_learning_rate.m

This script reproduces the learning-rate analysis presented in Figure 2 of the manuscript.

The script evaluates how different technology learning rates affect:

- overall innovation adoption trajectories;
- transient adoption disparity between communities;
- cumulative disparity over time.

The simulations are performed for three learning-rate scenarios:

- E = 0.16 (low learning rate)
- E = 0.18 (moderate learning rate)
- E = 0.20 (high learning rate)

while keeping the intergroup interaction parameter fixed.

## Requirements

The code was developed and tested using:

- MATLAB R2024b

Required MATLAB toolboxes:

- Optimization Toolbox (for `fsolve`)
- Statistics and Machine Learning Toolbox (for `gamcdf`)

## Model parameters

The baseline simulations use:

- Total population size: n = 8000
- Two communities:
  - high-income community
  - low-income community
- Initial adoption conditions as specified in the manuscript
- Time horizon: 120 model time units

## Running the code

To reproduce Figure 2:

1. Open MATLAB.
2. Navigate to the repository folder.
3. Run:
