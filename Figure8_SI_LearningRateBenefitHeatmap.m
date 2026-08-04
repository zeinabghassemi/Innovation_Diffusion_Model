%% Supplementary Figure S8: Learning Rate, Initial Price, and Benefit Sensitivity
%
% This script reproduces Supplementary Figure S8.
%
% The analysis examines how the interaction between technology
% learning rate, initial innovation price, and adoption benefit
% influences equilibrium adoption and disparity outcomes.
%
% Parameters explored:
% - Learning rate (E)
% - Initial price (rho_0)
% - Adoption benefit (q)
%
% Outputs:
% - Equilibrium adoption heatmaps
% - Equilibrium disparity heatmaps
%
% MATLAB version:
% MATLAB R2024b

clear;
clc;
close all;

% model's parameters
global xi P eta etai b a n E rho0 alpha_H alpha_L t_initial t_final

xi = 1;
P = 0.05;
p2 = 0.045;

eta = 0.75; %h
etai = 1 - eta; 
b = 0.15; %q
a = 0;
E = 0.2; % the experience parameter E for solar panels
% E1 = 0.16, E2 = 0.18, E3 = .
rho0 = 35; 
n = 8000;

t_initial = 0;
t_final = 1200;

% parameters for the gamma distribution for community 1: G_H
alpha_H = 4;
beta_H = 4;

% parameters for the gamma distribution for community 2: G_L
alpha_L = alpha_H;
beta_L = 2;

Y0_H = 0.01;
Y0_L = 0.01;
Y0 = (Y0_H + Y0_L) / 2;

Z0_H = Y0_H * (1 - gamcdf(rho0, alpha_H, beta_H));
Z0_L = Y0_L * (1 - gamcdf(rho0, alpha_H, beta_L));
Z0 = (Z0_H + Z0_L) / 2;

% the time span
tspan = [t_initial, t_final];

% initial Z, Y and price vectors
Y = zeros(1, t_final);
Z = zeros(1, t_final);
price = rho0 * ones(1, t_final);
Y_H = zeros(1, t_final);
Y_L = zeros(1, t_final);
Z_H = zeros(1, t_final);
Z_L = zeros(1, t_final);
Disparity_OverallVector = zeros(1, t_final);

% Initial conditions
Y(1) = Y0;
Y_L(1) = Y0_L;
Y_H(1) = Y0_H;
Z(1) = Z0;
Z_L(1) = Z0_L;
Z_H(1) = Z0_H;
price(1) = rho0;

% the functions are solved for each iteration of the loop
ode = @(t, Y, Zi, Zj) (xi * P + eta * (b * (1-a) + a * (2 * Zi - 1)) * Zi + etai * (b * (1-a) + a * (2 * Zj - 1)) * Zj) * (1 - Y) - p2 * (Y - Zi);

equationZ = @(Z, price, Y, beta) Y * (1 - gamcdf(price, alpha_H, beta)) - Z;

equationPrice = @(price, Z, rho0, n, E) rho0 * (n * Z) ^ -E - price;

% the Z and price equations for each iteration of the loop
for ts = 2:t_final
    [t, YH_sol] = ode45(@(t, Y_H) ode(t, Y_H, Z_H(ts-1), Z_L(ts-1)), [ts-1, ts], Y_H(ts-1));
    Y_H(ts) = YH_sol(end) ;

    [t, YL_sol] = ode45(@(t, Y_L) ode(t, Y_L, Z_L(ts-1), Z_H(ts-1)), [ts-1, ts], Y_L(ts-1));
    Y_L(ts) = YL_sol(end) ;

    Y(ts) = (Y_H(ts) + Y_L(ts)) / 2;

    price(ts) = fsolve(@(price) equationPrice(price, Z(ts-1), rho0, n , E), price(ts-1));

    Z_H(ts) = fsolve(@(Z_H) equationZ(Z_H, price(ts), Y_H(ts), beta_H), Z_H(ts-1));
    Z_L(ts) = fsolve(@(Z_L) equationZ(Z_L, price(ts), Y_L(ts), beta_L), Z_L(ts-1));

    Z(ts) = (Z_H(ts) + Z_L(ts)) / 2;

    Dsiparity_OverallVector(ts) = abs(Z_H(ts) - Z_L(ts));
end

Dsiparity_Overall = sum(abs(Z_H - Z_L));
Disparity_SteadyState = abs(Z_H(t_final) - Z_L(t_final));
Disp_ZZ_SteadyState = abs(Z_H(t_final) + Z_L(t_final))/2;

%%% Heat maps
% Example data (Replace these with actual data for each heatmap)
data1 = [0.37 0.46 0.50; 
         0.35 0.45 0.499;
         0.34 0.44 0.498]; 

data2 = [0.72 0.59 0.48; 
         0.73 0.61 0.49;
         0.74 0.62 0.50];

data3 = [0.27 0.37 0.45; 
         0.25 0.35 0.44;
         0.24 0.34 0.43];

data4 = [0.81 0.71 0.61; 
         0.82 0.73 0.63;
         0.83 0.74 0.64];

data5 = [0.18 0.27 0.36; 
         0.17 0.26 0.34;
         0.16 0.25 0.33];

data6 = [0.87 0.80 0.72; 
         0.88 0.81 0.73;
         0.89 0.82 0.74];

% Define x and y labels
xLabels = {'\rho_0 = 25', '\rho_0 = 30', '\rho_0 = 35'}; 
yLabels = {'q = 0.1', 'q = 0.15', 'q = 0.2'}; 

% Create figure
figure;

% Create 6 heatmaps in a 2-row, 3-column grid
data_list = {data1, data2, data3, data4, data5, data6}; % Store all data in a cell array
titles = {'Disparity for E=0.16', 'Adoption for E=0.16', 'Disparity for E=0.18', ...
          'Adoption for E=0.18', 'Disparity for E=0.2', 'Adoption for E=0.2'}; % Titles for each heatmap

for i = 1:6
    subplot(3, 2, i); % Create subplot in a 2x3 grid
    h = heatmap(xLabels, yLabels, data_list{i});
    h.Title = titles{i}; % Assign title
    colormap('parula'); % Set colormap
end
 
%%%%%%%%% Z S1S7
figure(3);
plot(1:t_final, Z, 'm', 'LineWidth', 3);
title('Plot of Z(t)');
hold on
 
% plot(1:t_final, Z_H, 'm', 'LineWidth', 3);
% hold on
% plot(1:t_final, Z_L, 'r', 'LineWidth', 3);
% hold on
% legend('Z','Z_H','Z_L')


figure(4);
plot(1:t_final, Y, 'm', 'LineWidth', 3);
title('Plot of Y(t)');
hold on

% plot(1:t_final, Y_H, 'm', 'LineWidth', 3);
% hold on
% plot(1:t_final, Y_L, 'r', 'LineWidth', 3);
% hold on
% legend('Y','Y_H','Y_L')


figure(5);
plot(1:t_final, price, 'm', 'LineWidth', 3);
title('Plot of price(t)');
hold on

figure(6);
plot(1:t_final, Z_H, 'm', 'LineWidth', 3);
title('Plot of Z_H(t)');
hold on

figure(7);
plot(1:t_final, Z_L, 'm', 'LineWidth', 3);
title('Plot of Z_L(t)');
hold on

figure(8);
plot(1:t_final, Y_H, 'm', 'LineWidth', 3);
title('Plot of Y_H(t)');
hold on

figure(9);
plot(1:t_final, Y_L, 'm', 'LineWidth', 3);
title('Plot of Y_L(t)');
hold on