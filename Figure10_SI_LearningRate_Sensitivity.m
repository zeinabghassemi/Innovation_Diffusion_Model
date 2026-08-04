%% Supplementary Figure S10: Learning Rate Sensitivity Analysis
%
% This script reproduces Supplementary Figure S10.
%
% The analysis evaluates how technology learning rates affect:
% - overall adoption dynamics;
% - transient disparity;
% - cumulative disparity.
%
% The simulations are conducted under fixed intergroup interaction
% and varying learning rates.
%
% MATLAB version:
% MATLAB R2024b

clear;
clc;
close all;

% model's parameters
global xi P eta b a n E rho0 alpha_H alpha_L beta_L beta_H t_initial t_final

xi = .05;
P = 0.05;
eta = 0.86;
etai = 1-eta;
b = 0.4;
a = 0.1;
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
%title('Adoption Proportion for Learning rates (E)');
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
%title('Disparity for Learning rates (E)');
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
%%%%%% Second E
E = 0.18;

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

%%%%% Third E
E = 0.202;

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
plot(1:t_final, Z, ':', 'Color','g', 'LineWidth', 4);
hold on
legend('E=0.16', 'E=0.18', 'E=0.2')

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
plot(1:t_final, Dsiparity_OverallVector(1:t_final),':', 'Color', ' g', 'LineWidth', 4);
hold on
legend('E=0.16', 'E=0.18', 'E=0.2')
ylim([0, 1]);
% Create a figure
% Create a figure
figure(7);
a = [0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.2]; % x-axis values (parameter h)
Dispchcl = [41.59, 39.59, 37.44, 35.2, 32.91, 30.63, 28.38, 26.19, 24.08]; % y-axis values

% Create a bar chart
b = bar(a, Dispchcl, 'FaceColor', 'flat'); % Bar chart with flat coloring

% Set all bars to blue
b.FaceColor = 'b'; % Blue color for all bars

% Add labels, title, and legend
xlabel('Learning Rate', 'Interpreter', 'latex');
ylabel('Cumulative Disparity', 'Interpreter', 'latex');
ylim([20, 55]);
%xlim([0.157, 0.203]);
%title('Cumulative Disparity for Learning rates', 'Interpreter', 'latex');

% Customize grid and aesthetics
grid on;
set(gca, 'FontSize', 12);

% Add colorbar for reference (though all bars are now blue)

hold off;
