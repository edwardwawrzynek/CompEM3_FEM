clc; clear; close all;

intermediatePlots = false;
mesh_h = @(lambda) lambda/80;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PEC cylinder, TM_z  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

c0 = 3*1e8; % speed of light
f0 = 10e9; % frequency
lambda0 = c0/f0; % lambda
k0 = 2*pi/lambda0; % wave vector magnitude
geo_a = 0.015; geo_b = 3*lambda0; % spatial
beta = -k0^2; % beta
alpha = 1; % alpha

% save properties
properties_a = {c0; f0; lambda0; k0; beta; alpha; beta; geo_a; geo_b};

% set up the problem
a = FEMProblem(alpha, beta, 0);

% create the geometry
a = a.makeRing(geo_a, geo_b);

% boundary conditions
% (1) Inner PEC Boundary
a = a.setDirichletBC(5:8, 0);
% (2) First Order ABC
gamma_abc = 1j*k0 + 1/(2*geo_b);
kappa_abc = @(x,y) (gamma_abc - 1j*k0*x/geo_b) .* exp(-1j*k0*x);
a = a.setNeumannBC(1:4, gamma_abc, kappa_abc);

% create the mesh
a = a.createMesh(mesh_h(lambda0));

if intermediatePlots
    fig_a = figure('Color','White', 'Units','normalized', 'Position', [0.02 0.04 0.96 0.5]);
    t_a = tiledlayout(fig_a, 1, 3, 'TileSpacing','compact', 'Padding','compact');

    % plot the geometry
    nexttile;
    a.plotGeo();
    set(gcf, 'Color', 'w');
    set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
    box on; grid on;
    axis equal tight;
    colormap inferno;
    title("Geometry")
    xlabel("x [m]")
    ylabel("y [m]")

    % plot the mesh
    nexttile;
    a.plotMesh(false);
    set(gcf, 'Color', 'w');
    set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
    box on; grid on;
    axis equal tight;
    colormap inferno;
    title("Mesh")
    xlabel("x [m]")
    ylabel("y [m]")

    % plot phi
    nexttile;
    a = a.assemble();
    set(gcf, 'Color', 'w');
    set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
    box on; grid on;
    axis equal tight;
    colormap inferno;
    title("\Phi(x,y)")
    xlabel("x [m]")
    ylabel("y [m]")

else
    % plot phi
    figure();
    a = a.assemble();
    set(gcf, 'Color', 'w');
    set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
    box on; grid on;
    axis equal tight;
    colormap inferno;
    title("|E(x,y)| [V/m]")
    xlabel("x [m]")
    ylabel("y [m]")
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PEC cylinder, TE_z  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

c0 = 3*1e8; % speed of light
f0 = 10e9; % frequency
lambda0 = c0/f0; % lambda
k0 = 2*pi/lambda0; % wave vector magnitude
geo_a = 0.015; geo_b = 3*lambda0; % spatial
beta = -k0^2; % beta
alpha = 1; % alpha

% save properties
properties_b = {c0; f0; lambda0; k0; beta; alpha; beta; geo_a; geo_b};

% set up the problem
b = FEMProblem(alpha, beta, 0);

% create the geometry
b = b.makeRing(geo_a, geo_b);

% boundary conditions
% (1) Inner PEC Boundary
b = b.setNeumannBC(5:8);
% (2) First Order ABC
gamma_abc = 1j*k0 + 1/(2*geo_b);
kappa_abc = @(x,y) (gamma_abc - 1j*k0*x/geo_b) .* exp(-1j*k0*x);
b = b.setNeumannBC(1:4, gamma_abc, kappa_abc);

% create the mesh
b = b.createMesh(mesh_h(lambda0));

if intermediatePlots
    fig_b = figure('Color','White', 'Units','normalized', 'Position', [0.02 0.04 0.96 0.5]);
    t_b = tiledlayout(fig_b, 1, 3, 'TileSpacing','compact', 'Padding','compact');

    % plot the geometry
    nexttile;
    b.plotGeo();
    set(gcf, 'Color', 'w');
    set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
    box on; grid on;
    axis equal tight;
    colormap inferno;
    title("Geometry")
    xlabel("x [m]")
    ylabel("y [m]")

    % plot the mesh
    nexttile;
    b.plotMesh(false);
    set(gcf, 'Color', 'w');
    set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
    box on; grid on;
    axis equal tight;
    colormap inferno;
    title("Mesh")
    xlabel("x [m]")
    ylabel("y [m]")

    % plot phi
    nexttile;
    b = b.assemble();
    set(gcf, 'Color', 'w');
    set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
    box on; grid on;
    axis equal tight;
    colormap inferno;
    title("\Phi(x,y)")
    xlabel("x [m]")
    ylabel("y [m]")

else
    % plot phi
    figure();
    b = b.assemble();
    set(gcf, 'Color', 'w');
    set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
    box on; grid on;
    axis equal tight;
    colormap inferno;
    title("|H(x,y)| [A/m]")
    xlabel("x [m]")
    ylabel("y [m]")
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dielectric cylinder, TM_z  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c0 = 3*1e8; % speed of light
f0 = 10e9; % frequency
lambda0 = c0/f0; % lambda
k0 = 2*pi/lambda0; % wave vector magnitude
eps_r = 9; % dielectric permittivity
geo_a = 0.015; geo_b = 3*lambda0; % spatial
beta = @(rad) -k0^2*(1 + (eps_r-1)*(rad < geo_a)); % beta
alpha = 1; % alpha

% save properties
properties_c = {c0; f0; lambda0; k0; beta; alpha; beta; geo_a; geo_b; eps_r};

% set up the problem
c = FEMProblem(alpha, beta, 0);

% create the geometry
c = c.makeDisk(geo_a, geo_b);

% boundary conditions
% (1) First Order ABC
gamma_abc = 1j*k0 + 1/(2*geo_b);
kappa_abc = @(x,y) (gamma_abc - 1j*k0*x/geo_b) .* exp(-1j*k0*x);
c = c.setNeumannBC(1:4, gamma_abc, kappa_abc);

% create the mesh
c = c.createMesh(mesh_h(lambda0));

if intermediatePlots
    fig_c = figure('Color','White', 'Units','normalized', 'Position', [0.02 0.04 0.96 0.5]);
    t_c = tiledlayout(fig_c, 1, 3, 'TileSpacing','compact', 'Padding','compact');

    % plot the geometry
    nexttile;
    c.plotGeo();
    set(gcf, 'Color', 'w');
    set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
    box on; grid on;
    axis equal tight;
    colormap inferno;
    title("Geometry")
    xlabel("x [m]")
    ylabel("y [m]")

    % plot the mesh
    nexttile;
    c.plotMesh(false);
    set(gcf, 'Color', 'w');
    set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
    box on; grid on;
    axis equal tight;
    colormap inferno;
    title("Mesh")
    xlabel("x [m]")
    ylabel("y [m]")

    % plot phi
    nexttile;
    c = c.assemble();
    set(gcf, 'Color', 'w');
    set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
    box on; grid on;
    axis equal tight;
    colormap inferno;
    title("\Phi(x,y)")
    xlabel("x [m]")
    ylabel("y [m]")

else
    % plot phi
    figure();
    c = c.assemble();
    set(gcf, 'Color', 'w');
    set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
    box on; grid on;
    axis equal tight;
    colormap inferno;
    title("|E(x,y)| [V/m]")
    xlabel("x [m]")
    ylabel("y [m]")
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dielectric cylinder, TE_z  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c0 = 3*1e8; % speed of light
f0 = 10e9; % frequency
lambda0 = c0/f0; % lambda
k0 = 2*pi/lambda0; % wave vector magnitude
eps_r = 9; % dielectric permittivity
geo_a = 0.015; geo_b = 3*lambda0; % spatial
beta = -k0^2; % beta
alpha = @(rad) 1./(1 + (eps_r-1)*(rad < geo_a)); % alpha

% save properties
properties_d = {c0; f0; lambda0; k0; beta; alpha; beta; geo_a; geo_b; eps_r};

% set up the problem
d = FEMProblem(alpha, beta, 0);

% create the geometry
d = d.makeDisk(geo_a, geo_b);

% boundary conditions
% (1) First Order ABC
gamma_abc = 1j*k0 + 1/(2*geo_b);
kappa_abc = @(x,y) (gamma_abc - 1j*k0*x/geo_b) .* exp(-1j*k0*x);
d = d.setNeumannBC(1:4, gamma_abc, kappa_abc);

% create the mesh
d = d.createMesh(mesh_h(lambda0));

if intermediatePlots
    fig_d = figure('Color','White', 'Units','normalized', 'Position', [0.02 0.04 0.96 0.5]);
    t_d = tiledlayout(fig_d, 1, 3, 'TileSpacing','compact', 'Padding','compact');

    % plot the geometry
    nexttile;
    d.plotGeo();
    set(gcf, 'Color', 'w');
    set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
    box on; grid on;
    axis equal tight;
    colormap inferno;
    title("Geometry")
    xlabel("x [m]")
    ylabel("y [m]")

    % plot the mesh
    nexttile;
    d.plotMesh(false);
    set(gcf, 'Color', 'w');
    set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
    box on; grid on;
    axis equal tight;
    colormap inferno;
    title("Mesh")
    xlabel("x [m]")
    ylabel("y [m]")

    % plot phi
    nexttile;
    d = d.assemble();
    set(gcf, 'Color', 'w');
    set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
    box on; grid on;
    axis equal tight;
    colormap inferno;
    title("\Phi(x,y)")
    xlabel("x [m]")
    ylabel("y [m]")

else
    % plot phi
    figure();
    d = d.assemble();
    set(gcf, 'Color', 'w');
    set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
    box on; grid on;
    axis equal tight;
    colormap inferno;
    title("|H(x,y)| [A/m]")
    xlabel("x [m]")
    ylabel("y [m]")
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Total, Incident, Scattered %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CASE 1: TMz PEC
plotIST_fields(a, properties_a);

% CASE 2: TEz PEC
plotIST_fields(b, properties_b);

% CASE 3: TMz Dielectric
plotIST_fields(c, properties_c);

% CASE 4: TEz Dielectric
plotIST_fields(d, properties_d);

%%
%%%%%%%%%%%%%%%%%%
%%% Echo Width %%%
%%%%%%%%%%%%%%%%%%

% sweep the angle
angle = linspace(0, 2*pi, 360);

% CASE 1: TMz PEC
k0 = properties_a{4}; R = properties_a{8};
sigma_simulated  = simulated_echo_width(a, properties_a, 0.05, 360, 50);
sigma_analytical = bistatic_echo_width_tm_pec(angle, k0, R);

% plot
figure();
plot(angle*180/pi, sigma_simulated, 'LineWidth', 2); hold on;
plot(angle*180/pi, sigma_analytical, 'LineWidth', 2);
set(gca,'FontName', 'Times', 'FontSize', 18);
set(gcf, 'Color', 'w')
legend("Simulated", "Analytic");
title("Bistatic Echowidth; TM_z PEC")
xlabel("\phi [deg]");
ylabel("\sigma(\phi)");
grid on; box on;

% CASE 2: TEz PEC
k0 = properties_b{4}; R = properties_b{8};
sigma_simulated  = simulated_echo_width(b, properties_b, 0.05, 360, 50);
sigma_analytical = bistatic_echo_width_te_pec(angle, k0, R);

% plot
figure();
plot(angle*180/pi, sigma_simulated, 'LineWidth', 2); hold on;
plot(angle*180/pi, sigma_analytical, 'LineWidth', 2);
set(gca,'FontName', 'Times', 'FontSize', 18);
set(gcf, 'Color', 'w')
legend("Simulated", "Analytic");
title("Bistatic Echowidth; TE_z PEC")
xlabel("\phi [deg]");
ylabel("\sigma(\phi)");
grid on; box on;

% CASE 3: TMz Dielectric
k0 = properties_c{4}; R = properties_c{8}; eps_r = properties_c{10};
sigma_simulated  = simulated_echo_width(c, properties_c, 0.05, 360, 50);
sigma_analytical = bistatic_echo_width_tm_eps(angle, k0, R, eps_r);

% plot
figure();
plot(angle*180/pi, sigma_simulated, 'LineWidth', 2); hold on;
plot(angle*180/pi, sigma_analytical, 'LineWidth', 2);
set(gca,'FontName', 'Times', 'FontSize', 18);
set(gcf, 'Color', 'w')
legend("Simulated", "Analytic");
title("Bistatic Echowidth; TM_z Dielectric")
xlabel("\phi [deg]");
ylabel("\sigma(\phi)");
grid on; box on;

% CASE 4: TEz Dielectric
k0 = properties_d{4}; R = properties_d{8}; eps_r = properties_d{10};
sigma_simulated  = simulated_echo_width(d, properties_d, 0.05, 360, 10);
sigma_analytical = bistatic_echo_width_te_eps(angle, k0, R, eps_r);

% plot
figure();
plot(angle*180/pi, sigma_simulated, 'LineWidth', 2); hold on;
plot(angle*180/pi, sigma_analytical, 'LineWidth', 2);
set(gca,'FontName', 'Times', 'FontSize', 18);
set(gcf, 'Color', 'w')
legend("Simulated", "Analytic");
title("Bistatic Echowidth; TE_z Dielectric")
xlabel("\phi [deg]");
ylabel("\sigma(\phi)");
grid on; box on;