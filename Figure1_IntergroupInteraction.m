%% Figure 1: Impact of Intergroup Interaction
%
% This script reproduces Figure 1 of the manuscript.
%
% The script evaluates innovation diffusion under different
% levels of intergroup interaction between communities.
%
% MATLAB version: R2024b

clear;
clc;
close all;

% Model parameters

global xi P eta b a n E rho0 alpha_H alpha_L beta_L beta_H t_initial t_final

xi = 0.05;
P = 0.05;
eta = 1;
etai = 1-eta;
b = 0.2;
a = 0;
E = 0.168; % the experience parameter E for solar panels
rho0 = 30; 
n = 8000;
t_initial = 0;
t_final = 120;
p2 = 0.05;
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
Z_H1 = zeros(1, t_final);
Z_L1 = zeros(1, t_final);


% Initial conditions
Y(1) = Y0;
Y_L(1) = Y0_L;
Y_H(1) = Y0_H;
Z(1) = Z0;
Z_L(1) = Z0_L;
Z_H(1) = Z0_H;
price(1) = rho0;

% the functions are solved for each iteration of the loop
ode = @(t, Y, Zi, Zj) ((xi * P + eta * (b * (1-a) + a * (2 * Zi - 1)) * Zi + etai * (b * (1-a) + a * (2 * Zj - 1)) * Zj) * (1 - Y)- p2 * (Y - Zi));

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
Z_H1 = Z_H;
Z_L1 = Z_L;

Dsiparity_Overall = sum(abs(Z_H - Z_L));
Disparity_SteadyState = abs(Z_H(t_final) - Z_L(t_final));

% plot Z(t)
figure(1);
plot(1:t_final, Z, 'r', 'LineWidth', 3);
xlabel('Time');
ylabel('Overall Adoption');
%title('Adoption Proportion for Homophily Values (h)');
hold on
grid on;

% plot Y(t)
figure(2);
plot(1:t_final, Y, 'm', 'LineWidth', 3);
xlabel('Time (t)');
ylabel('Y(t)');
title('Plot of Y(t) for the learning parameter E=0.168');
hold on
grid on

% plot price
figure(3);
plot(1:t_final, price, 'm', 'LineWidth', 3);
xlabel('Time (t)');
ylabel('Price');
title('Plot of the Price for the learning parameter E=0.168');
grid on
hold on

% plot Dsiparity_OverallVector
figure(13);
plot(1:t_final, Dsiparity_OverallVector(1:t_final), 'r', 'LineWidth', 3);
xlabel('Time');
ylabel('Transient Disparity');
%title('Disparity for Homophily Values (h)');
hold on
grid on;

% plot ZH(t)
figure(4);
plot(1:t_final, Z_H, 'm', 'LineWidth', 3);
xlabel('Time (t)');
ylabel('Zh(t)');
title('Plot of ZH(t) for the learning parameter E=0.168');
hold on
grid on;

% plot ZL(t)
figure(5);
plot(1:t_final, Z_L, 'm', 'LineWidth', 3);
xlabel('Time (t)');
ylabel('Zl(t)');
title('Plot of ZL(t) for the learning parameter E=0.168');
hold on
grid on
Z1= Z;
%%%%%% Second etai
eta = 0.86;
etai = 0.14;

% initial Z, Y and price vectors
Y = zeros(1, t_final);
Z = zeros(1, t_final);
price = rho0 * ones(1, t_final);
Y_H = zeros(1, t_final);
Y_L = zeros(1, t_final);
Z_H = zeros(1, t_final);
Z_L = zeros(1, t_final);
Disparity_OverallVector = zeros(1, t_final);

Z_H2 = zeros(1, t_final);
Z_L2 = zeros(1, t_final);

% Initial conditions
Y(1) = Y0;
Y_L(1) = Y0_L;
Y_H(1) = Y0_H;
Z(1) = Z0;
Z_L(1) = Z0_L;
Z_H(1) = Z0_H;
price(1) = rho0;

% the functions are solved for each iteration of the loop
ode = @(t, Y, Zi, Zj) ((xi * P + eta * (b * (1-a) + a * (2 * Zi - 1)) * Zi + etai * (b * (1-a) + a * (2 * Zj - 1)) * Zj) * (1 - Y)- p2 *(Y - Zi));

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

Z_H2 = Z_H;
Z_L2 = Z_L;

Dsiparity_Overall = sum(abs(Z_H - Z_L));
Disparity_SteadyState = abs(Z_H(t_final) - Z_L(t_final));

% plot Z(t)
figure(1);
plot(1:t_final, Z, '--', 'Color', 'b', 'LineWidth', 3);
hold on

% plot Y(t)
figure(2);
plot(1:t_final, Y, 'g', 'LineWidth', 3);
hold on

% plot price
figure(3);
plot(1:t_final, price, 'g', 'LineWidth', 3);
hold on

% plot Dsiparity_OverallVector
figure(13);
plot(1:t_final, Dsiparity_OverallVector(1:t_final), '--', 'Color', 'b', 'LineWidth', 3);
hold on

% plot ZH(t)
figure(4);
plot(1:t_final, Z_H, 'g', 'LineWidth', 3);

% plot ZL(t)
figure(5);
plot(1:t_final, Z_L, 'g', 'LineWidth', 3);
Z2 = Z;
%%%%% Third etai
eta = 0.63;
etai = 0.37;

% initial Z, Y and price vectors
Y = zeros(1, t_final);
Z = zeros(1, t_final);
price = rho0 * ones(1, t_final);
Y_H = zeros(1, t_final);
Y_L = zeros(1, t_final);
Z_H = zeros(1, t_final);
Z_L = zeros(1, t_final);
Disparity_OverallVector = zeros(1, t_final);

Z_H3 = zeros(1, t_final);
Z_L3 = zeros(1, t_final);

% Initial conditions
Y(1) = Y0;
Y_L(1) = Y0_L;
Y_H(1) = Y0_H;
Z(1) = Z0;
Z_L(1) = Z0_L;
Z_H(1) = Z0_H;
price(1) = rho0;

% the functions are solved for each iteration of the loop
ode = @(t, Y, Zi, Zj) ((xi * P + eta * (b * (1-a) + a * (2 * Zi - 1)) * Zi + etai * (b * (1-a) + a * (2 * Zj - 1)) * Zj) * (1 - Y)- p2 *(Y - Zi));

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

Z_H3 = Z_H;
Z_L3 = Z_L;

Dsiparity_Overall = sum(abs(Z_H - Z_L));
Disparity_SteadyState = abs(Z_H(t_final) - Z_L(t_final));

% plot Z(t)
figure(1);
plot(1:t_final, Z, ':', 'Color', 'g' , 'LineWidth', 4);
hold on

% plot Y(t)
figure(2);
plot(1:t_final, Y, 'b', 'LineWidth', 3);
hold on

% plot price
figure(3);
plot(1:t_final, price, 'b', 'LineWidth', 3);
hold on

% plot Dsiparity_OverallVector
figure(13);
plot(1:t_final, Dsiparity_OverallVector(1:t_final), ':', 'Color',' g', 'LineWidth', 4);
hold on

% plot ZH(t)
figure(4);
plot(1:t_final, Z_H, 'b', 'LineWidth', 3);

% plot ZL(t)
figure(5);
plot(1:t_final, Z_L, 'b', 'LineWidth', 3);

Z3 = Z;
%%%%%% Forth etai
eta = 0.5;
etai = 0.5;

% initial Z, Y and price vectors
Y = zeros(1, t_final);
Z = zeros(1, t_final);
price = rho0 * ones(1, t_final);
Y_H = zeros(1, t_final);
Y_L = zeros(1, t_final);
Z_H = zeros(1, t_final);
Z_L = zeros(1, t_final);
Disparity_OverallVector = zeros(1, t_final);

Z_H4 = zeros(1, t_final);
Z_L4 = zeros(1, t_final);

% Initial conditions
Y(1) = Y0;
Y_L(1) = Y0_L;
Y_H(1) = Y0_H;
Z(1) = Z0;
Z_L(1) = Z0_L;
Z_H(1) = Z0_H;
price(1) = rho0;

% the functions are solved for each iteration of the loop
ode = @(t, Y, Zi, Zj) ((xi * P + eta * (b * (1-a) + a * (2 * Zi - 1)) * Zi + etai * (b * (1-a) + a * (2 * Zj - 1)) * Zj) * (1 - Y)- p2 * (Y - Zi));

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

Z_H4 = Z_H;
Z_L4 = Z_L;
Z4 = Z;

Dsiparity_Overall = sum(abs(Z_H - Z_L));
Disparity_SteadyState = abs(Z_H(t_final) - Z_L(t_final));

% plot Z(t)
figure(1);
plot(1:t_final, Z, '-.', 'Color', 'k', 'LineWidth', 3);
hold on
legend('I = 0', 'I = 0.14', 'I = 0.37', 'I = 0.5')
ylim([0,1]);
xlim([0,120]);

% plot Y(t)
figure(2);
plot(1:t_final, Y, '--r', 'LineWidth', 3);
hold on
legend('\eta_{-i} = 0', '\eta_{-i} = 0.125', '\eta_{-i} = 0.375', '\eta_{-i} = 0.5')


% plot price
figure(3);
plot(1:t_final, price, '--r', 'LineWidth', 3);
hold on
legend('\eta_{-i} = 0', '\eta_{-i} = 0.125', '\eta_{-i} = 0.375', '\eta_{-i} = 0.5')

% plot Dsiparity_OverallVector
figure(13);
plot(1:t_final, Dsiparity_OverallVector(1:t_final), '-.', 'Color','k', 'LineWidth', 3);
legend('I = 0', 'I = 0.14', 'I = 0.37', 'I = 0.5')
hold on
ylim([0,1]);
xlim([0,120]);

Dist12H = sum(abs(Z_H1 - Z_H2));
Dist12L = sum(abs(Z_L1 - Z_L2));

Dist23H = sum(abs(Z_H2 - Z_H3));
Dist23L = sum(abs(Z_L2 - Z_L3));

Dist34H = sum(abs(Z_H3 - Z_H4));
Dist34L = sum(abs(Z_L3 - Z_L4));

% plot ZH(t)
figure(4);
plot(1:t_final, Z_H, '--r', 'LineWidth', 3);

% plot ZL(t)
figure(5);
plot(1:t_final, Z_L, '--r', 'LineWidth', 3);

% plot zhs
figure(6);
plot(1:t_final, Z_H1, 'r', 'LineWidth', 3);
hold on
plot(1:t_final, Z_H2, 'b', 'LineWidth', 3);
hold on
plot(1:t_final, Z_H3, 'g', 'LineWidth', 3);
hold on
plot(1:t_final, Z_H4, 'm', 'LineWidth', 3);
legend('\eta_{-i} = 0', '\eta_{-i} = 0.125', '\eta_{-i} = 0.375', '\eta_{-i} = 0.5')
grid on
ylim([0,1]);
xlim([0,120]);
xlabel('Time (t)');
ylabel('Z_H(t)');
title('Plot of ZH(t) for the learning parameter E=0.168');
% plot zls
figure(7);
plot(1:t_final, Z_L1, 'r', 'LineWidth', 3);
hold on
plot(1:t_final, Z_L2, 'b', 'LineWidth', 3);
hold on
plot(1:t_final, Z_L3, 'g', 'LineWidth', 3);
hold on
plot(1:t_final, Z_L4, 'm', 'LineWidth', 3);
legend('\eta_{-i} = 0', '\eta_{-i} = 0.125', '\eta_{-i} = 0.375', '\eta_{-i} = 0.5')
grid on
ylim([0,1]);
xlim([0,120]);
xlabel('Time (t)');
ylabel('ZL(t)');
title('Plot of ZL(t) for the learning parameter E=0.168');
% Example data (2x3 matrix)
data = [2.8383, 1.0336, 0.4055, 0.0955, 0.1035; 
        1.344, 1.4432, 1.4864, 1.49, 1.4635]; % Your actual data should replace this

% Define x and y labels (Plain text for the heatmap)
%xLabels = {'0 \rightarrow 0.05', '0.05 \rightarrow 0.1', '0.1 \rightarrow 0.15', '0.15 \rightarrow 0.2', '0.2 \rightarrow 0.25', '0.25 \rightarrow 0.3', '0.3 \rightarrow 0.35', '0.35 \rightarrow 0.4', '0.4 \rightarrow 0.45', '0.45 \rightarrow 0.5'}; % Temporary plain labels for x-axis
xLabels = {'0 \rightarrow 0.1', '0.1 \rightarrow 0.2', '0.2 \rightarrow 0.3', '0.3 \rightarrow 0.4', '0.4 \rightarrow 0.5'}; % Temporary plain labels for x-axis
yLabels = {'Disparity within c_Ls', 'Disparity within c_Hs'}; % Labels for y-axis

% Ensure the figure is new and there is no conflict with 'hold on'
figure; % Create a new figure

% Create heatmap
h = heatmap(xLabels, yLabels, data);

% Customize heatmap
h.Title = 'Heatmap of disparity within communities vs Interactability';
h.XLabel = '\eta_{-i}';

% Optional: Customize colormap
colormap('parula');

% Data from the table
%a = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1]; % x-axis values (parameter h)
a = [0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0];
Dispchcl = [32.2, 33.3, 34.5, 35.8, 37, 38.5, 40, 41.8, 43.8, 46.6, 50.9]; % y-axis values

% Create a figure
figure(17);

% Create a bar chart
b = bar(a, Dispchcl, 'FaceColor', 'flat'); % Bar chart with flat coloring

% Set all bars to blue
b.FaceColor = 'b'; % Blue color for all bars

% Add labels, title, and legend
xlabel('Intergroup Interaction', 'Interpreter', 'latex');
ylabel('Cumulative Disparity', 'Interpreter', 'latex');
ylim([30, 55]);
%title('Cumulative Disparity for Homophily Values (h)', 'Interpreter', 'latex');

% Customize grid and aesthetics
grid on;
set(gca, 'FontSize', 12);

hold off;

% --------------------------------
% Create separate zoom-in figure
% --------------------------------

x_zoom = 100:120;

figure(9);
clf;
plot(x_zoom, Z1(x_zoom), 'r', 'LineWidth', 20);
hold on;
plot(x_zoom, Z2(x_zoom), '--b', 'LineWidth', 20);
plot(x_zoom, Z3(x_zoom), ':g', 'LineWidth', 20);
plot(x_zoom, Z4(x_zoom), '-.k', 'LineWidth', 20);

grid on;
box on;

xlabel('Time', 'FontName', 'Times New Roman', 'FontSize', 140);
ylabel('Adoption', 'FontName', 'Times New Roman', 'FontSize', 140);

xlim([100 120]);

y_min = min([Z1(x_zoom), Z2(x_zoom), Z3(x_zoom), Z4(x_zoom)]);
y_max = max([Z1(x_zoom), Z2(x_zoom), Z3(x_zoom), Z4(x_zoom)]);

y_margin = 0.002;   % try 0.001, 0.002, or 0.005
ylim([y_min - y_margin, y_max + y_margin]);

xticks([100 110 120]);

ytick_values = linspace(y_min, y_max, 4);
yticks(ytick_values);
ytickformat('%.3f');

set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 140);
set(gca, 'LineWidth', 0.8);

% Make the figure small and clean
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1 1 2.4 1.7]);

% Save as EPS
set(gcf, 'Renderer', 'painters');
print(gcf, 'Fig1_zoom_in.eps', '-depsc', '-r600');