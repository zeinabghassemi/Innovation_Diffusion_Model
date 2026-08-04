## Figure 1: Impact of intergroup interaction

### Script

`Figure1_intergroup_interaction.m`

### Description

This script reproduces the intergroup interaction analysis presented in Figure 1 of the manuscript.

The simulation investigates how different levels of connectivity between high-income and low-income communities affect:

- overall innovation adoption dynamics;
- transient adoption disparity between communities;
- cumulative disparity over time.

The model considers four interaction scenarios:

- isolated communities: \(I=0\)
- weak intergroup interaction: \(I=0.14\)
- moderate intergroup interaction: \(I=0.37\)
- fully connected communities: \(I=0.5\)

while keeping the learning rate fixed.

### Outputs

The script generates:

- adoption trajectories \(Z(t)\);
- transient disparity trajectories;
- subgroup adoption trajectories;
- cumulative disparity across interaction levels;
- zoomed-in equilibrium adoption dynamics.

### Running the code

Open MATLAB and run:
