clc; clear; close all;

%%%%%%%%%%%%%%%%%%%
%%% Unit test 1 %%%
%%%%%%%%%%%%%%%%%%%
fig_a = figure('Color','White', 'Units','normalized', 'Position', [0.02 0.04 0.96 0.5]);
t_a = tiledlayout(fig_a, 1, 3, 'TileSpacing','compact', 'Padding','compact');

% define the FEM problem and assign properties:
% ep_r_center = 1       pde parameter (alpha)
% beta = 0              pde parameter
% gamma = 0             neumann bc term
% kappa_n = 0           neumann bc term
% f = 0                 forcing
a = FEMProblem(1, 0, 0, 0, 0);

% create the geometry
geo_a = 0.5; geo_b = 2;
a = a.makeRing(geo_a, geo_b);

% set dirichlet BC
a = a.setDirichletBC(1:4,0);
a = a.setDirichletBC(5:8,1);

% create the mesh
a = a.createMesh(0.025);

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

% Capacitance check
C_analytical = 2 * pi * 8.85 / log(geo_b/geo_a);
C_approx = a.Cap * 1e12; 
per_Err = abs(C_analytical - C_approx)/C_analytical * 100;

disp("Capcitance (Unit test 1) = " + C_approx);
disp("Percent Error = " + per_Err + "%");

%%%%%%%%%%%%%%%%%%%
%%% Unit test 2 %%%
%%%%%%%%%%%%%%%%%%%
fig_b = figure('Color','White', 'Units','normalized', 'Position', [0.02 0.04 0.96 0.5]);
t_b = tiledlayout(fig_b, 1, 3, 'TileSpacing','compact', 'Padding','compact');

% define the FEM problem and assign properties:
% ep_r_center = 1       pde parameter (alpha)
% beta = 0              pde parameter
% gamma = 0             neumann bc term
% kappa_n = 0           neumann bc term
% f = 0                 forcing
b = FEMProblem(1, 0, 0, 0, 0);

% create the geometry
geo_a = 0.25; geo_b = 1;
b = b.makeSlice(geo_a, geo_b, pi/4);

% set dirichlet BC
b = b.setDirichletBC(3,0);
b = b.setDirichletBC(4,1);

% create a mesh such that h = 0.05 is the target element size
b = b.createMesh(0.025);

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

% Capacitance check
C_analytical = 2 * pi * 8.85 / log(geo_b/geo_a);
C_approx = (b.Cap * 8) * 1e12; 
per_Err = abs(C_analytical - C_approx)/C_analytical * 100;

disp("Capcitance (Unit test 2) = " + C_approx);
disp("Percent Error = " + per_Err + "%");

%%%%%%%%%%%%%%%%%%%
%%% Unit test 3 %%%
%%%%%%%%%%%%%%%%%%%
fig_c = figure('Color','White', 'Units','normalized', 'Position', [0.02 0.04 0.96 0.5]);
t_c = tiledlayout(fig_c, 1, 3, 'TileSpacing','compact', 'Padding','compact');

% define the FEM problem and assign properties:
% ep_r_center = eps_r(rho)       pde parameter (alpha)
% beta = 0                       pde parameter
% gamma = 0                      neumann bc term
% kappa_n = 0                    neumann bc term
% f = 0                          forcing
geo_a = 0.5; geo_b = 2;
eps_m = 5;
eps = @(rho) eps_m * (geo_a/rho);
c = FEMProblem(eps, 0, 0, 0, 0);

% create the geometry
c = c.makeSlice(geo_a, geo_b, pi/4);

% set dirichlet BC
c = c.setDirichletBC(3,0);
c = c.setDirichletBC(4,1);

% create a mesh such that h = 0.05 is the target element size
c = c.createMesh(0.025);

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

% Capacitance check
C_analytical = 2 * pi * 8.85 * eps_m * geo_a / (geo_b - geo_a);
C_approx = (c.Cap * 8) * 1e12; 
per_Err = abs(C_analytical - C_approx)/C_analytical * 100;

disp("Capcitance (Unit test 3) = " + C_approx);
disp("Percent Error = " + per_Err + "%");

%%
%%%%%%%%%%%%%%%%%%%
%%% Unit test 4 %%%
%%%%%%%%%%%%%%%%%%%
fig_d = figure('Color','White', 'Units','normalized', 'Position', [0.02 0.04 0.96 0.5]);
t_d = tiledlayout(fig_d, 1, 3, 'TileSpacing','compact', 'Padding','compact');

% define the FEM problem and assign properties:
% ep_r_center = eps_r(rho)       pde parameter (alpha)
% beta = 0                       pde parameter
% gamma = 0                      neumann bc term
% kappa_n = rho_s / ep0          neumann bc term
% f = 0                          forcing
geo_a = 0.5; geo_b = 2;
eps_m = 5;
eps = @(rho) eps_m * (geo_a/rho);
d = FEMProblem(eps, 0, 0, 0, 0);

% create the geometry
d = d.makeSlice(geo_a, geo_b, pi/4);

% set dirichlet BC
d = d.setDirichletBC(3,0);

% set neumann BC
rho_s = @(x,y) eps_m / (geo_b - geo_a);
d = d.setNeumannBC(4, 0, rho_s);

% create a mesh such that h = 0.05 is the target element size
d = d.createMesh(0.025);

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
