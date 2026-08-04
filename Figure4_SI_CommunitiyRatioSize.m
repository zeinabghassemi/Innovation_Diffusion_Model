%% Supplementary Figure S4: Community Size Ratio Sensitivity Analysis
%
% This script reproduces Supplementary Figure S4.
%
% The analysis examines how the relative sizes of the high-income
% and low-income communities affect adoption and disparity outcomes.
%
% Two sensitivity analyses are performed:
%
% (1) Community size ratio under different levels of
%     intergroup interaction/homophily.
%
% (2) Community size ratio under different technology learning rates.
%
% Outputs:
% - Equilibrium disparity
% - Overall adoption
%
% MATLAB version:
% MATLAB R2024b

clear;
clc;
close all;
% model's parameters
global xi P eta b a n E rho0 alpha_H alpha_L beta_L beta_H t_initial t_final

xi = 1;
P = 0.05;
b = 0.2;
a = 0;
E = 0.16; % the experience parameter E for solar panels
rho0 = 33; 
n = 8000;
t_initial = 0;
t_final = 500;
p2 = 0.05;

% parameters for the gamma distribution for community 1: G_H
alpha_H = 4;
beta_H = 4;
beta_L = 2;

% arrays to store results
%i_values = 0.1:0.01:0.9;
i_values = 0.5:0.1:2;
disparity_h1 = zeros(size(i_values));
disparity_h2 = zeros(size(i_values));
disparity_h3 = zeros(size(i_values));
disparity_h4 = zeros(size(i_values));

z_h1 = zeros(size(i_values));
z_h2 = zeros(size(i_values));
z_h3 = zeros(size(i_values));
z_h4 = zeros(size(i_values));

nH = n/2;
nL = n/2;
% loop i
index = 1;
for i = i_values
    nH = (i/(i+1)) * n;
    nL = nH/i;
    % h values and their corresponding etai
    h_values = [1, 0.875, 0.8, 0.5];
    disparities = zeros(1, length(h_values));
    zs = zeros(1, length(h_values));

    % loop h
    for h_index = 1:length(h_values)
        h = h_values(h_index);
        etai = 1 - h;

        % initial conditions
        Y0_H = 0.01;
        Y0_L = 0.01;
        Y0 = (Y0_H + Y0_L) / 2;

        Z0_H = Y0_H * (1 - gamcdf(rho0, alpha_H, beta_H));
        Z0_L = Y0_L * (1 - gamcdf(rho0, alpha_H, beta_L));
        Z0 = Z0_H * nH/n + Z0_L * nL/n;

        % initial vectors
        Z_H = zeros(1, t_final);
        Z_L = zeros(1, t_final);
        price = rho0 * ones(1, t_final);
        Y_H = zeros(1, t_final);
        Y_L = zeros(1, t_final);

        Y_H(1) = Y0_H;
        Y_L(1) = Y0_L;
        Z_H(1) = Z0_H;
        Z_L(1) = Z0_L;
        price(1) = rho0;

        % ODE and equations
        ode = @(t, Y, Zi, Zj, ni, nj) ((xi * P + h * (b * (1-a) + a * (2 * Zi - 1)) * Zi * ni/(n/2) + etai * (b * (1-a) + a * (2 * Zj - 1)) * Zj * nj/(n/2)) * (1 - Y) - p2 * (Y - Zi));
        equationZ = @(Z, price, Y, beta) Y * (1 - gamcdf(price, alpha_H, beta)) - Z;
        equationPrice = @(price, Z, rho0, n, E) rho0 * (n * Z) ^ -E - price;

        % time span
        for ts = 2:t_final
            [t, YH_sol] = ode45(@(t, Y_H) ode(t, Y_H, Z_H(ts-1), Z_L(ts-1), nH, nL), [ts-1, ts], Y_H(ts-1));
            Y_H(ts) = YH_sol(end);

            [t, YL_sol] = ode45(@(t, Y_L) ode(t, Y_L, Z_L(ts-1), Z_H(ts-1), nL, nH), [ts-1, ts], Y_L(ts-1));
            Y_L(ts) = YL_sol(end);

            price(ts) = fsolve(@(price) equationPrice(price, (Z_H(ts-1) * nH/n + Z_L(ts-1) * nL/n), rho0, n, E), price(ts-1));
            Z_H(ts) = fsolve(@(Z_H) equationZ(Z_H, price(ts), Y_H(ts), beta_H), Z_H(ts-1));
            Z_L(ts) = fsolve(@(Z_L) equationZ(Z_L, price(ts), Y_L(ts), beta_L), Z_L(ts-1));
        end

        % steady-state disparity
        disparities(h_index) = abs(Z_H(t_final) - Z_L(t_final));
        zs(h_index) = (Z_H(ts-1) * nH/n + Z_L(ts-1) * nL/n);
    end

    % disparities for each h
    disparity_h1(index) = disparities(1);
    disparity_h2(index) = disparities(2);
    disparity_h3(index) = disparities(3);
    disparity_h4(index) = disparities(4);

    z_h1(index) = zs(1);
    z_h2(index) = zs(2);
    z_h3(index) = zs(3);
    z_h4(index) = zs(4);

    index = index + 1;
end

% Plots
figure(1);
tiledlayout(2,2); % Using tiled layout for better spacing

% First subplot
nexttile;
hold on;
plot(i_values, disparity_h1, 'r', 'LineWidth', 3, 'DisplayName', 'I = 0');
plot(i_values, disparity_h2, 'b--', 'LineWidth', 3, 'DisplayName', 'I = 0.14');
plot(i_values, disparity_h3, 'g:', 'LineWidth', 4, 'DisplayName', 'I = 0.37');
plot(i_values, disparity_h4, 'k.-', 'LineWidth', 3, 'DisplayName', 'I = 0.5');
xlabel('Ratio of Communities Size: c_H to c_L');
ylim([0.35 0.6]);
ylabel('Equilibrium Disparity');
%title('Disparity vs Sizes Ratio for h');
legend('show');
grid on;
hold off;

figure(1);
nexttile;
hold on;
plot(i_values, z_h1, 'r', 'LineWidth', 3, 'DisplayName', 'I = 0');
plot(i_values, z_h2, 'b--', 'LineWidth', 3, 'DisplayName', 'I = 0.14');
plot(i_values, z_h3, 'g:', 'LineWidth', 4, 'DisplayName', 'I = 0.37');
plot(i_values, z_h4, 'k.-', 'LineWidth', 3, 'DisplayName', 'I = 0.5');
xlabel('Ratio of Communities Size: c_H to c_L');
%xlim([1.5 9]);
ylabel('Overall Adoption');
%title('Adoption Proportion vs Sizes Ratio for h');
legend('show');
grid on;
hold off;

%%%%%%%% Learning parameter E
xi = 1;
P = 0.05;
b = 0.2;
a = 0;
h = 0.875;
rho0 = 30; 
n = 8000;
t_initial = 0;
t_final = 500;
p2 = 0.05;

% parameters for the gamma distribution for community 1: G_H
alpha_H = 4;
beta_H = 4;
beta_L = 2;

% arrays to store results
%i_values = 0.1:0.01:0.9;
i_values = 0.5:0.1:2;
disparity_e1 = zeros(size(i_values));
disparity_e2 = zeros(size(i_values));
disparity_e3 = zeros(size(i_values));

z_e1 = zeros(size(i_values));
z_e2 = zeros(size(i_values));
z_e3 = zeros(size(i_values));
z_e4 = zeros(size(i_values));
nH = n/2;
nL = n/2;
% loop i
index = 1;
for i = i_values
    nH = (i/(i+1)) * n;
    nL = nH/i;
 
    % e values and their corresponding etai
    e_values = [0.16, 0.18, 0.202];
    disparities = zeros(1, length(e_values));
    zs = zeros(1, length(e_values));

    % loop E
    for e_index = 1:length(e_values)
        e = e_values(e_index);

        % initial conditions
        Y0_H = 0.01;
        Y0_L = 0.01;
        Y0 = (Y0_H + Y0_L) / 2;

        Z0_H = Y0_H * (1 - gamcdf(rho0, alpha_H, beta_H));
        Z0_L = Y0_L * (1 - gamcdf(rho0, alpha_H, beta_L));
        Z0 = (Z0_H + Z0_L) / 2;

        % initial vectors
        Z_H = zeros(1, t_final);
        Z_L = zeros(1, t_final);
        price = rho0 * ones(1, t_final);
        Y_H = zeros(1, t_final);
        Y_L = zeros(1, t_final);

        Y_H(1) = Y0_H;
        Y_L(1) = Y0_L;
        Z_H(1) = Z0_H;
        Z_L(1) = Z0_L;
        price(1) = rho0;

        % ODE and equations
        ode = @(t, Y, Zi, Zj, ni, nj) ((xi * P + h * (b * (1-a) + a * (2 * Zi - 1)) * Zi * ni/(n/2) + etai * (b * (1-a) + a * (2 * Zj - 1)) * Zj * nj/(n/2)) * (1 - Y) - p2 * (Y - Zi));
        equationZ = @(Z, price, Y, beta) Y * (1 - gamcdf(price, alpha_H, beta)) - Z;
        equationPrice = @(price, Z, rho0, n, E) rho0 * (n * Z) ^ -e - price;

        % time span
        for ts = 2:t_final
            [t, YH_sol] = ode45(@(t, Y_H) ode(t, Y_H, Z_H(ts-1), Z_L(ts-1), nH, nL), [ts-1, ts], Y_H(ts-1));
            Y_H(ts) = YH_sol(end);

            [t, YL_sol] = ode45(@(t, Y_L) ode(t, Y_L, Z_L(ts-1), Z_H(ts-1), nL, nH), [ts-1, ts], Y_L(ts-1));
            Y_L(ts) = YL_sol(end);

            price(ts) = fsolve(@(price) equationPrice(price, (Z_H(ts-1) * nH/n + Z_L(ts-1) * nL/n), rho0, n, E), price(ts-1));
            Z_H(ts) = fsolve(@(Z_H) equationZ(Z_H, price(ts), Y_H(ts), beta_H), Z_H(ts-1));
            Z_L(ts) = fsolve(@(Z_L) equationZ(Z_L, price(ts), Y_L(ts), beta_L), Z_L(ts-1));
        end

        % steady-state disparity
        disparities(e_index) = abs(Z_H(t_final) - Z_L(t_final));
        zs(e_index) = (Z_H(ts-1) * nH/n + Z_L(ts-1) * nL/n);
    end
    % disparities for each e
    disparity_e1(index) = disparities(1);
    disparity_e2(index) = disparities(2);
    disparity_e3(index) = disparities(3);

    z_e1(index) = zs(1);
    z_e2(index) = zs(2);
    z_e3(index) = zs(3);

    index = index + 1;
end

% Plots
figure(1);
nexttile;
hold on;
plot(i_values, disparity_e1, 'r', 'LineWidth', 3, 'DisplayName', 'E = 0.16');
plot(i_values, disparity_e2, 'b--', 'LineWidth', 3, 'DisplayName', 'E = 0.18');
plot(i_values, disparity_e3, 'g:', 'LineWidth', 4, 'DisplayName', 'E = 0.2');
xlabel('Ratio of Communities Size: c_H to c_L');
ylim([0 1]);
ylabel('Equilibrium Disparity');
%title('Disparity vs Sizes Ratio for E');
legend('show');
grid on;
hold off;

figure(1);
nexttile;
hold on;
plot(i_values, z_e1, 'r', 'LineWidth', 3, 'DisplayName', 'E = 0.16');
plot(i_values, z_e2, 'b--', 'LineWidth', 3, 'DisplayName', 'E = 0.18');
plot(i_values, z_e3, 'g:', 'LineWidth', 4, 'DisplayName', 'E = 0.2');
ylim([0 1]);
xlabel('Ratio of Communities Size: c_H to c_L');
ylabel('Overall Adoption');
%title('Adoption Proportion vs Sizes Ratio for E');
legend('show');
grid on;
hold off;