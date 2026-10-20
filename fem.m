clc; clear; close all;

eps0 = 8.85 * 1e-12; % permittivity of vacuum

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
eps_r = 1;
a = FEMProblem(eps_r, 0, 0, 0, 0);

% create the geometry
geo_a = 0.25; geo_b = 1;
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
C_analytical = 2 * pi * eps0 / log(geo_b/geo_a) * 1e12;
C_approx = a.calculateCap() * 1e12; 
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
eps_r = 1;
b = FEMProblem(eps_r, 0, 0, 0, 0);

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
C_analytical = 2 * pi * eps0 / log(geo_b/geo_a) * 1e12;
C_approx = (b.calculateCap() * 8) * 1e12; 
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
geo_a = 0.25; geo_b = 1;
eps_m = 1;
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

% plot epsilon
figure();
c = c.plotEps();
set(gcf, 'Color', 'w');
set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
box on; grid on;
axis equal tight;
colormap inferno;
title("\epsilon_r(x,y)")
xlabel("x [m]")
ylabel("y [m]")

% Capacitance check
C_analytical = 2 * pi * eps0 * eps_m * geo_a / (geo_b - geo_a) * 1e12;
C_approx = (c.calculateCap() * 8) * 1e12; 
per_Err = abs(C_analytical - C_approx)/C_analytical * 100;

disp("Capcitance (Unit test 3) = " + C_approx);
disp("Percent Error = " + per_Err + "%");

%%%%%%%%%%%%%%%%%%%
%%% Unit test 4 %%%
%%%%%%%%%%%%%%%%%%%
fig_d = figure('Color','White', 'Units','normalized', 'Position', [0.02 0.04 0.96 0.5]);
t_d = tiledlayout(fig_d, 1, 3, 'TileSpacing','compact', 'Padding','compact');

% define the FEM problem and assign properties:
% ep_r_center = eps(rad)         pde parameter (alpha)
% beta = 0                       pde parameter
% gamma = 0                      neumann bc term
% kappa_n = rho_s/eps0           neumann bc term
% f = 0                          forcing
geo_a = 0.25; geo_b = 1;
eps_m = 1;
eps = @(rho) eps_m * (geo_a/rho);
rho_s = 1e-9;
kappa_n = rho_s/eps0;
d = FEMProblem(eps, 0, 0, kappa_n, 0);

% create the geometry
d = d.makeSlice(geo_a, geo_b, pi/4);

% set BC
d = d.setDirichletBC(3,0);
d = d.setNeumannBC(4);

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

% Phi test:
% nodes along the bottom edge
bottom_nodes = findNodes(d.mesh, "region", "Edge", 1); 
x_coordinate = d.mesh.Nodes(1,bottom_nodes);

% sort by increasing radius outward
[rad, idx] = sort(x_coordinate);

% calculate both
phi_analytical = kappa_n*(geo_b-rad)/eps_m;
phi_approx = d.phi(bottom_nodes(idx));

% plot phi_analytical and phi_approx
figure(); hold on;
plot(rad, phi_approx, 'b-', "LineWidth", 2);
plot(rad, phi_analytical, 'r--', "LineWidth", 2);
set(gcf, 'Color', 'w');
set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
box on; grid on;
axis tight;
legend(["Numerical","Analytical"]);
title("Unit Test 4")
xlabel("\rho [m]")
ylabel("\Phi [V]")

%%%%%%%%%%%%%%%%%%%
%%% Unit test 5 %%%
%%%%%%%%%%%%%%%%%%%
fig_e = figure('Color','White', 'Units','normalized', 'Position', [0.02 0.04 0.96 0.5]);
t_e = tiledlayout(fig_e, 1, 3, 'TileSpacing','compact', 'Padding','compact');

% define the FEM problem and assign properties:
% ep_r_center = 5                pde parameter (alpha)
% beta = 0                       pde parameter
% gamma = 0                      neumann bc term
% kappa_n = rho_s/eps0           neumann bc term
% f = rho_v/eps0                 forcing
geo_a = 0.5; geo_b = 8;
eps_r = 5;
A = 10; B = 1/10;
kappa_n = A/(2*pi*geo_a); % rho_s/eps0
f = B/(pi*geo_a^2); % rho_v/eps0
e = FEMProblem(eps_r, 0, 0, kappa_n, f);

% create the geometry
e = e.makeSlice(geo_a, geo_b, pi/4);

% set BC
e = e.setDirichletBC(3,0);
e = e.setNeumannBC(4);

% create a mesh such that h = 0.05 is the target element size
e = e.createMesh(0.05);

% plot the geometry
nexttile;
e.plotGeo();
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
e.plotMesh(false);
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
e = e.assemble();
set(gcf, 'Color', 'w');
set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
box on; grid on;
axis equal tight;
colormap inferno;
title("\Phi(x,y)")
xlabel("x [m]")
ylabel("y [m]")

% Phi test:
% nodes along the bottom edge
bottom_nodes = findNodes(e.mesh, "region", "Edge", 1); 
x_coordinate = e.mesh.Nodes(1,bottom_nodes);

% sort by increasing radius outward
[rad, idx] = sort(x_coordinate);

% calculate both
phi_analytical = (A - B).*log(geo_b./rad)./(2*pi*eps_r) + B.*(geo_b^2 - rad.^2)./(4*pi*geo_a^2*eps_r);
phi_approx = e.phi(bottom_nodes(idx));

% plot phi_analytical and phi_approx
figure(); hold on;
plot(rad, phi_approx, 'b-', "LineWidth", 2);
plot(rad, phi_analytical, 'r--', "LineWidth", 2);
set(gcf, 'Color', 'w');
set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
box on; grid on;
axis tight;
legend(["Numerical","Analytical"]);
title("Unit Test 5")
xlabel("\rho [m]")
ylabel("\Phi [V]")