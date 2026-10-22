clc; clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PEC cylinder, TM_z  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig_a = figure('Color','White', 'Units','normalized', 'Position', [0.02 0.04 0.96 0.5]);
t_a = tiledlayout(fig_a, 1, 3, 'TileSpacing','compact', 'Padding','compact');

c = 3*1e8; % speed of light
f0 = 10e9; % frequency
lambda0 = c/f0; % lambda
k0 = 2*pi/lambda0; % wave vector magnitude

% set up the problem
a = FEMProblem(1, -k0^2, 0);

% create the geometry
geo_a = 0.015; geo_b = 3*lambda0;
a = a.makeRing(geo_a, geo_b);

% boundary conditions
% (1) Inner PEC Boundary
a = a.setDirichletBC(5:8, 0);
% (2) First Order ABC
gamma_abc = 1j*k0 + 1/(2*geo_b);
kappa_abc = @(x,y) (gamma_abc - 1j*k0*x/geo_b) .* exp(-1j*k0*x);
a = a.setNeumannBC(1:4, gamma_abc, kappa_abc);

% create the mesh
a = a.createMesh(lambda0/20);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PEC cylinder, TE_z  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig_b = figure('Color','White', 'Units','normalized', 'Position', [0.02 0.04 0.96 0.5]);
t_b = tiledlayout(fig_b, 1, 3, 'TileSpacing','compact', 'Padding','compact');

c = 3*1e8; % speed of light
f0 = 10e9; % frequency
lambda0 = c/f0; % lambda
k0 = 2*pi/lambda0; % wave vector magnitude

% set up the problem
b = FEMProblem(1, -k0^2, 0);

% create the geometry
geo_a = 0.015; geo_b = 3*lambda0;
b = b.makeRing(geo_a, geo_b);

% boundary conditions
% (1) Inner PEC Boundary
b = b.setNeumannBC(5:8);
% (2) First Order ABC
gamma_abc = 1j*k0 + 1/(2*geo_b);
kappa_abc = @(x,y) (gamma_abc - 1j*k0*x/geo_b) .* exp(-1j*k0*x);
b = b.setNeumannBC(1:4, gamma_abc, kappa_abc);

% create the mesh
b = b.createMesh(lambda0/20);

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
