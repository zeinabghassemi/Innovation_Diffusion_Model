%% Supplementary Figure S2: Income Distribution Sensitivity Analysis
%
% This script reproduces Supplementary Figure S2.
%
% The analysis investigates how assumptions regarding the income
% distribution influence innovation adoption and equilibrium disparity.
%
% Four sensitivity analyses are performed:
%
% (1) Variation in mean income ratio between communities.
% (2) Variation in income variance ratio.
% (3) Learning parameter sensitivity under mean changes.
% (4) Learning parameter sensitivity under variance changes.
%
% Outputs:
% - Normalized equilibrium disparity
% - Normalized overall adoption proportion
%
% MATLAB implementation of the innovation diffusion model.

clear;
clc;
close all;
% model's parameters
global xi P eta b a n E rho0 alpha_H alpha_L beta_L beta_H t_initial t_final

xi = 1;
P = 0.05;
b = 0.2;
a = 0;
E = 0.168; % the experience parameter E for solar panels
rho0 = 20; 
n = 8000;
t_initial = 0;
t_final = 500;
p2 = 0.045;

% parameters for the gamma distribution for community 1: G_H
alpha_H = 4;
beta_H = 4;

% arrays to store results
%i_values = 0.1:0.01:0.9;
i_values = 1:0.1:5;
disparity_h1 = zeros(size(i_values));
disparity_h2 = zeros(size(i_values));
disparity_h3 = zeros(size(i_values));
disparity_h4 = zeros(size(i_values));

z_h1 = zeros(size(i_values));
z_h2 = zeros(size(i_values));
z_h3 = zeros(size(i_values));
z_h4 = zeros(size(i_values));

% loop i
index = 1;
for i = i_values
    %beta_L = i * beta_H;
    beta_L = 2;
    %beta_H = sqrt(i) * beta_L; %for variance
    beta_H = i * beta_L; %for means and standard deviations

    % h values and their corresponding etai
    h_values = [1, 0.875, 0.625, 0.5];
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
        ode = @(t, Y, Zi, Zj) ((xi * P + h * (b * (1-a) + a * (2 * Zi - 1)) * Zi + etai * (b * (1-a) + a * (2 * Zj - 1)) * Zj) * (1 - Y) - p2 * (Y - Zi));
        equationZ = @(Z, price, Y, beta) Y * (1 - gamcdf(price, alpha_H, beta)) - Z;
        equationPrice = @(price, Z, rho0, n, E) rho0 * (n * Z) ^ -E - price;

        % time span
        for ts = 2:t_final
            [t, YH_sol] = ode45(@(t, Y_H) ode(t, Y_H, Z_H(ts-1), Z_L(ts-1)), [ts-1, ts], Y_H(ts-1));
            Y_H(ts) = YH_sol(end);

            [t, YL_sol] = ode45(@(t, Y_L) ode(t, Y_L, Z_L(ts-1), Z_H(ts-1)), [ts-1, ts], Y_L(ts-1));
            Y_L(ts) = YL_sol(end);

            price(ts) = fsolve(@(price) equationPrice(price, (Z_H(ts-1) + Z_L(ts-1)) / 2, rho0, n, E), price(ts-1));
            Z_H(ts) = fsolve(@(Z_H) equationZ(Z_H, price(ts), Y_H(ts), beta_H), Z_H(ts-1));
            Z_L(ts) = fsolve(@(Z_L) equationZ(Z_L, price(ts), Y_L(ts), beta_L), Z_L(ts-1));
        end

        % steady-state disparity
        disparities(h_index) = abs(Z_H(t_final) - Z_L(t_final));
        zs(h_index) = (Z_H(t_final) + Z_L(t_final)) / 2;
    end

    disparities = disparities - min(disparities);
    zs = zs - min(zs);

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
hold on;
plot(i_values, disparity_h1, 'r', 'LineWidth', 3, 'DisplayName', 'h = 1');
plot(i_values, disparity_h2, 'b', 'LineWidth', 3, 'DisplayName', 'h = 0.86');
plot(i_values, disparity_h3, 'g', 'LineWidth', 4, 'DisplayName', 'h = 0.63');
plot(i_values, disparity_h4, 'k', 'LineWidth', 3, 'DisplayName', 'h = 0.5');
xlabel('Ratio of Income c_H to c_L');
xlim([1 5]);
ylabel('Normalized Equilibrium Disparity');
%title('Means and Standard Deviation Ratios for Different Homophily Values (h)');
legend('show');
grid on;
hold off;

figure(2);
hold on;
plot(i_values, z_h1, 'r', 'LineWidth', 3, 'DisplayName', 'h = 1');
plot(i_values, z_h2, 'b', 'LineWidth', 3, 'DisplayName', 'h = 0.86');
plot(i_values, z_h3, 'g', 'LineWidth', 4, 'DisplayName', 'h = 0.63');
plot(i_values, z_h4, 'k', 'LineWidth', 3, 'DisplayName', 'h = 0.5');
xlabel('Ratio of Income c_H to c_L');
xlim([1 5]);
ylabel('Normalized Equilibrium Adoption Proportion');
title('Means and Standard Deviation Ratios');
legend('show');
grid on;
hold off;

figure(5);
tiledlayout(2,2); % Using tiled layout for better spacing

% First subplot
nexttile;
hold on;
plot(i_values, disparity_h1, 'r', 'LineWidth', 3, 'DisplayName', 'I = 0');
plot(i_values, disparity_h2, 'b--', 'LineWidth', 3, 'DisplayName', 'I = 0.14');
plot(i_values, disparity_h3, 'g:', 'LineWidth', 4, 'DisplayName', 'I = 0.37');
plot(i_values, disparity_h4, 'k-.', 'LineWidth', 3, 'DisplayName', 'I = 0.5');
xlabel('Ratio of Income c_H to c_L');
xlim([1 5]);
ylabel('Normalized Equilibrium Disparity');
title('Means and Std. Dev. Ratios');
legend('show');
grid on;
hold off;

% Second subplot
figure(5)
nexttile;
hold on;
plot(i_values, z_h1, 'r', 'LineWidth', 3, 'DisplayName', 'I = 0');
plot(i_values, z_h2, 'b--', 'LineWidth', 3, 'DisplayName', 'I = 0.14');
plot(i_values, z_h3, 'g:', 'LineWidth', 4, 'DisplayName', 'I = 0.37');
plot(i_values, z_h4, 'k-.', 'LineWidth', 3, 'DisplayName', 'I = 0.5');
xlabel('Ratio of Income c_H to c_L');
xlim([1 5]);
ylabel('Normalized Overall Adoption');
title('Means and Std. Dev. Ratios');
legend('show');
grid on;
hold off;

%%%%% Variance
xi = 1;
P = 0.05;
b = 0.2;
a = 0;
E = 0.168; % the experience parameter E for solar panels
rho0 = 20; 
n = 8000;
t_initial = 0;
t_final = 500;
p2 = 0.045;

% parameters for the gamma distribution for community 1: G_H
alpha_H = 4;
beta_H = 4;

% arrays to store results
%i_values = 0.1:0.01:0.9;
i_values = 1:0.1:5;
disparity_h1 = zeros(size(i_values));
disparity_h2 = zeros(size(i_values));
disparity_h3 = zeros(size(i_values));
disparity_h4 = zeros(size(i_values));

z_h1 = zeros(size(i_values));
z_h2 = zeros(size(i_values));
z_h3 = zeros(size(i_values));
z_h4 = zeros(size(i_values));

% loop i
index = 1;
for i = i_values
    %beta_L = i * beta_H;
    beta_L = 2;
    beta_H = sqrt(i) * beta_L; %for variance
    %beta_H = i * beta_L; %for means and standard deviations

    % h values and their corresponding etai
    h_values = [1, 0.875, 0.625, 0.5];
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
        ode = @(t, Y, Zi, Zj) ((xi * P + h * (b * (1-a) + a * (2 * Zi - 1)) * Zi + etai * (b * (1-a) + a * (2 * Zj - 1)) * Zj) * (1 - Y) - p2 * (Y - Zi));
        equationZ = @(Z, price, Y, beta) Y * (1 - gamcdf(price, alpha_H, beta)) - Z;
        equationPrice = @(price, Z, rho0, n, E) rho0 * (n * Z) ^ -E - price;

        % time span
        for ts = 2:t_final
            [t, YH_sol] = ode45(@(t, Y_H) ode(t, Y_H, Z_H(ts-1), Z_L(ts-1)), [ts-1, ts], Y_H(ts-1));
            Y_H(ts) = YH_sol(end);

            [t, YL_sol] = ode45(@(t, Y_L) ode(t, Y_L, Z_L(ts-1), Z_H(ts-1)), [ts-1, ts], Y_L(ts-1));
            Y_L(ts) = YL_sol(end);

            price(ts) = fsolve(@(price) equationPrice(price, (Z_H(ts-1) + Z_L(ts-1)) / 2, rho0, n, E), price(ts-1));
            Z_H(ts) = fsolve(@(Z_H) equationZ(Z_H, price(ts), Y_H(ts), beta_H), Z_H(ts-1));
            Z_L(ts) = fsolve(@(Z_L) equationZ(Z_L, price(ts), Y_L(ts), beta_L), Z_L(ts-1));
        end

        % steady-state disparity
        disparities(h_index) = abs(Z_H(t_final) - Z_L(t_final));
        zs(h_index) = (Z_H(t_final) + Z_L(t_final)) / 2;
    end

    disparities = disparities - min(disparities);
    zs = zs - min(zs);

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
figure(3);
hold on;
plot(i_values, disparity_h1, 'r', 'LineWidth', 3, 'DisplayName', 'I = 0');
plot(i_values, disparity_h2, 'b', 'LineWidth', 3, 'DisplayName', 'I = 0.14');
plot(i_values, disparity_h3, 'g', 'LineWidth', 4, 'DisplayName', 'I = 0.37');
plot(i_values, disparity_h4, 'k', 'LineWidth', 3, 'DisplayName', 'I = 0.5');
xlabel('Ratio of Income c_H to c_L');
xlim([1 5]);
ylabel('Normalized Equilibrium Disparity');
title('Variance Ratios');
legend('show');
grid on;
hold off;

figure(4);
hold on;
plot(i_values, z_h1, 'r', 'LineWidth', 3, 'DisplayName', 'I = 0');
plot(i_values, z_h2, 'b', 'LineWidth', 3, 'DisplayName', 'I = 0.14');
plot(i_values, z_h3, 'g', 'LineWidth', 4, 'DisplayName', 'I = 0.37');
plot(i_values, z_h4, 'k', 'LineWidth', 3, 'DisplayName', 'I = 0.5');
xlabel('Ratio of Income c_H to c_L');
xlim([1 5]);
ylabel('Normalized Overall Adoption');
title('Variance Ratios');
legend('show');
grid on;
hold off;

% Third subplot
figure(5)
nexttile;
hold on;
plot(i_values, disparity_h1, 'r', 'LineWidth', 3, 'DisplayName', 'I = 0');
plot(i_values, disparity_h2, 'b--', 'LineWidth', 3, 'DisplayName', 'I = 0.14');
plot(i_values, disparity_h3, 'g:', 'LineWidth', 4, 'DisplayName', 'I = 0.37');
plot(i_values, disparity_h4, 'k-.', 'LineWidth', 3, 'DisplayName', 'I = 0.5');
xlabel('Ratio of Income c_H to c_L');
xlim([1 5]);
ylabel('Normalized Equilibrium Disparity');
title('Variance Ratios');
legend('show');
grid on;
hold off;

% Fourth subplot
figure(5)
nexttile;
hold on;
plot(i_values, z_h1, 'r', 'LineWidth', 3, 'DisplayName', 'I = 0');
plot(i_values, z_h2, 'b--', 'LineWidth', 3, 'DisplayName', 'I = 0.14');
plot(i_values, z_h3, 'g:', 'LineWidth', 4, 'DisplayName', 'I = 0.37');
plot(i_values, z_h4, 'k-.', 'LineWidth', 3, 'DisplayName', 'I = 0.5');
xlabel('Ratio of Income c_H to c_L');
xlim([1 5]);
ylabel('Normalized Overall Adoption');
title('Variance Ratios');
legend('show');
grid on;
hold off;
%%%%%%%% Learning parameter E
xi = 1;
P = 0.05;
b = 0.2;
a = 0;
h = 0.875;
rho0 = 20; 
n = 8000;
t_initial = 0;
t_final = 500;
p2 = 0.045;

% parameters for the gamma distribution for community 1: G_H
alpha_H = 4;
beta_H = 4;

% arrays to store results
%i_values = 0.1:0.01:0.9;
i_values = 1:0.1:5;
disparity_e1 = zeros(size(i_values));
disparity_e2 = zeros(size(i_values));
disparity_e3 = zeros(size(i_values));

z_e1 = zeros(size(i_values));
z_e2 = zeros(size(i_values));
z_e3 = zeros(size(i_values));
z_e4 = zeros(size(i_values));

% loop i
index = 1;
for i = i_values
    beta_L = 2;
    beta_H = i * beta_L;

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
        ode = @(t, Y, Zi, Zj) ((xi * P + h * (b * (1-a) + a * (2 * Zi - 1)) * Zi + etai * (b * (1-a) + a * (2 * Zj - 1)) * Zj) * (1 - Y) - p2 * (Y - Zi));
        equationZ = @(Z, price, Y, beta) Y * (1 - gamcdf(price, alpha_H, beta)) - Z;
        equationPrice = @(price, Z, rho0, n, E) rho0 * (n * Z) ^ - e - price;

        % time span
        for ts = 2:t_final
            [t, YH_sol] = ode45(@(t, Y_H) ode(t, Y_H, Z_H(ts-1), Z_L(ts-1)), [ts-1, ts], Y_H(ts-1));
            Y_H(ts) = YH_sol(end);

            [t, YL_sol] = ode45(@(t, Y_L) ode(t, Y_L, Z_L(ts-1), Z_H(ts-1)), [ts-1, ts], Y_L(ts-1));
            Y_L(ts) = YL_sol(end);

            price(ts) = fsolve(@(price) equationPrice(price, (Z_H(ts-1) + Z_L(ts-1)) / 2, rho0, n, E), price(ts-1));
            Z_H(ts) = fsolve(@(Z_H) equationZ(Z_H, price(ts), Y_H(ts), beta_H), Z_H(ts-1));
            Z_L(ts) = fsolve(@(Z_L) equationZ(Z_L, price(ts), Y_L(ts), beta_L), Z_L(ts-1));
        end

        % steady-state disparity
        disparities(e_index) = (Z_H(t_final) - Z_L(t_final));
        zs(e_index) = (Z_H(t_final) + Z_L(t_final)) / 2;
    end

    % disparities for each h
    disparity_e1(index) = disparities(1);
    disparity_e2(index) = disparities(2);
    disparity_e3(index) = disparities(3);

    z_e1(index) = zs(1);
    z_e2(index) = zs(2);
    z_e3(index) = zs(3);

    index = index + 1;
end


figure(6);
tiledlayout(2,2); % Using tiled layout for better spacing

% First subplot
nexttile;
hold on;
plot(i_values, disparity_e1, 'r', 'LineWidth', 3, 'DisplayName', 'E = 0.16');
plot(i_values, disparity_e2, 'b--', 'LineWidth', 3, 'DisplayName', 'E = 0.18');
plot(i_values, disparity_e3, 'g:', 'LineWidth', 4, 'DisplayName', 'E = 0.2');
xlabel('Ratio of Income c_H to c_L');
xlim([1 5]);
ylabel('Equilibrium Disparity');
title('Means and Std. Dev. Ratios');
legend('show');
grid on;
hold off;

% Second subplot
figure(6)
nexttile;
hold on;
plot(i_values, z_e1, 'r', 'LineWidth', 3, 'DisplayName', 'E = 0.16');
plot(i_values, z_e2, 'b--', 'LineWidth', 3, 'DisplayName', 'E = 0.18');
plot(i_values, z_e3, 'g:', 'LineWidth', 4, 'DisplayName', 'E = 0.2');
xlabel('Ratio of Income c_H to c_L');
xlim([1 5]);
ylabel('Overall Adoption');
title('Means and Std. Dev. Ratios');
legend('show');
grid on;
hold off;
% Plots

%%%% Variance
xi = 1;
P = 0.05;
b = 0.2;
a = 0;
h = 0.875;
rho0 = 20; 
n = 8000;
t_initial = 0;
t_final = 500;
p2 = 0.045;

% parameters for the gamma distribution for community 1: G_H
alpha_H = 4;
beta_H = 4;

% arrays to store results
%i_values = 0.1:0.01:0.9;
i_values = 1:0.1:5;
disparity_e1 = zeros(size(i_values));
disparity_e2 = zeros(size(i_values));
disparity_e3 = zeros(size(i_values));

z_e1 = zeros(size(i_values));
z_e2 = zeros(size(i_values));
z_e3 = zeros(size(i_values));
z_e4 = zeros(size(i_values));

% loop i
index = 1;
for i = i_values
    beta_L = 2;
    beta_H = sqrt(i) * beta_L;

    % e values and their corresponding etai
    e_values = [0.135, 0.168, 0.202];
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
        ode = @(t, Y, Zi, Zj) ((xi * P + h * (b * (1-a) + a * (2 * Zi - 1)) * Zi + etai * (b * (1-a) + a * (2 * Zj - 1)) * Zj) * (1 - Y) - p2 * (Y - Zi));
        equationZ = @(Z, price, Y, beta) Y * (1 - gamcdf(price, alpha_H, beta)) - Z;
        equationPrice = @(price, Z, rho0, n, E) rho0 * (n * Z) ^ - e - price;

        % time span
        for ts = 2:t_final
            [t, YH_sol] = ode45(@(t, Y_H) ode(t, Y_H, Z_H(ts-1), Z_L(ts-1)), [ts-1, ts], Y_H(ts-1));
            Y_H(ts) = YH_sol(end);

            [t, YL_sol] = ode45(@(t, Y_L) ode(t, Y_L, Z_L(ts-1), Z_H(ts-1)), [ts-1, ts], Y_L(ts-1));
            Y_L(ts) = YL_sol(end);

            price(ts) = fsolve(@(price) equationPrice(price, (Z_H(ts-1) + Z_L(ts-1)) / 2, rho0, n, E), price(ts-1));
            Z_H(ts) = fsolve(@(Z_H) equationZ(Z_H, price(ts), Y_H(ts), beta_H), Z_H(ts-1));
            Z_L(ts) = fsolve(@(Z_L) equationZ(Z_L, price(ts), Y_L(ts), beta_L), Z_L(ts-1));
        end

        % steady-state disparity
        disparities(e_index) = (Z_H(t_final) - Z_L(t_final));
        zs(e_index) = (Z_H(t_final) + Z_L(t_final)) / 2;
    end

    % disparities for each h
    disparity_e1(index) = disparities(1);
    disparity_e2(index) = disparities(2);
    disparity_e3(index) = disparities(3);

    z_e1(index) = zs(1);
    z_e2(index) = zs(2);
    z_e3(index) = zs(3);

    index = index + 1;
end

figure(6)
nexttile;
hold on;
plot(i_values, disparity_e1, 'r', 'LineWidth', 3, 'DisplayName', 'E = 0.16');
plot(i_values, disparity_e2, 'b--', 'LineWidth', 3, 'DisplayName', 'E = 0.18');
plot(i_values, disparity_e3, 'g:', 'LineWidth', 4, 'DisplayName', 'E = 0.2');
xlabel('Ratio of Income c_H to c_L');
xlim([1 5]);
ylabel('Equilibrium Disparity');
title('Variance Ratios');
legend('show');
grid on;
hold off;

% Second subplot
figure(6)
nexttile;
hold on;
plot(i_values, z_e1, 'r', 'LineWidth', 3, 'DisplayName', 'E = 0.16');
plot(i_values, z_e2, 'b--', 'LineWidth', 3, 'DisplayName', 'E = 0.18');
plot(i_values, z_e3, 'g:', 'LineWidth', 4, 'DisplayName', 'E = 0.2');
xlabel('Ratio of Income c_H to c_L');
xlim([1 5]);
ylabel('Overall Adoption');
title('Variance Ratios');
legend('show');
grid on;
hold off;