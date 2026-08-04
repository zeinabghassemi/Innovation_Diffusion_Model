%% Figure 3: Impact of Seeding Policies
%
% This script reproduces the heatmaps presented in Figure 3
% of the manuscript.
%
% The script visualizes the effect of different seeding policies
% under varying:
%
% 1) intergroup interaction levels (I)
% 2) learning rates (E)
%
% Outputs:
% - Overall adoption heatmaps
% - Equilibrium disparity heatmaps
%
% MATLAB version:
% MATLAB R2024b

clear;
clc;
close all;

% Define the data matrix
data = [
    0.799, 0.801, 0.803, 0.805;
    0.735, 0.736, 0.737, 0.738;
    0.709, 0.710, 0.712, 0.713;
];

% Define row and column labels
row_labels = {'S_1', 'S_2', 'S_3'};
col_labels = {'I=0', 'I=0.14', 'I=0.37', 'I=0.5'};

% Create the heatmap
figure;
h = heatmap(col_labels, row_labels, data);

% Set color scaling to enhance contrast
h.Colormap = parula;
h.ColorLimits = [min(data(:)), max(data(:))];

% Add title
%title('Policy Effect Heatmap');

% Add labels
ylabel('Policy');
xlabel('Intergroup Interaction (I)');

% Add title
%title('Overall Adoption');

%%%%%%%%%%%%%%%%%
%%%%%%% Disparity 
%%%%%%%%%%%%%%%%%

data = [
    0.400, 0.397, 0.393, 0.391;
    0.29, 0.287, 0.284, 0.282;
    0.368, 0.364, 0.359, 0.355;
];

% Define row and column labels
row_labels = {'S_1', 'S_2', 'S_3'};
col_labels = {'I=0', 'I=0.14', 'I=0.37', 'I=0.5'};

% Create the heatmap
figure;
h = heatmap(col_labels, row_labels, data);

% Set color scaling to enhance contrast
h.Colormap = parula;
h.ColorLimits = [min(data(:)), max(data(:))];

% Add title
%title('Policy Effect Heatmap');

% Add labels
ylabel('Policy');
xlabel('Intergroup Interaction (I)');

% Add title
%title('Equilibrium Disparity');


%%%%%%%%%%%%%%%
% Learning rate
%%%%%%%%%%%%%%%
% Define the data matrix
data = [
    0.801, 0.846, 0.910;
    0.736, 0.789, 0.860;
    0.711, 0.767, 0.844;
];

% Define row and column labels
row_labels = {'S_1', 'S_2', 'S_3'};
col_labels = {'E=0.16', 'E=0.18', 'E=0.2'};

% Create the heatmap
figure;
h = heatmap(col_labels, row_labels, data);

% Set color scaling to enhance contrast
h.Colormap = parula;
h.ColorLimits = [min(data(:)), max(data(:))];

% Add labels
ylabel('Policy');
xlabel('Learning Rate (E)');

% Add title
%title('Overall Adoption');

%%%%%%%%%%%%%%%%%
%%%%%%% Disparity 
%%%%%%%%%%%%%%%%%

data = [
    0.397, 0.305, 0.167;
    0.287, 0.249, 0.182;
    0.364, 0.311, 0.223;
];

% Define row and column labels
row_labels = {'S_1', 'S_2', 'S_3'};
col_labels = {'E=0.16', 'E=0.18', 'E=0.2'};

% Create the heatmap
figure;
h = heatmap(col_labels, row_labels, data);

% Set color scaling to enhance contrast
h.Colormap = parula;
h.ColorLimits = [min(data(:)), max(data(:))];


% Add labels
ylabel('Policy');
xlabel('Learning Rate (E)');

% Add title
%title('Equilibrium Disparity');
