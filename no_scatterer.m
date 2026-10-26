clc; clear; close all;

intermediatePlots = true;
mesh_h = @(lambda) lambda/80;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dielectric cylinder, TM_z  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c0 = 3*1e8; % speed of light
f0 = 10e9; % frequency
lambda0 = c0/f0; % lambda
k0 = 2*pi/lambda0; % wave vector magnitude
eps_r = 1; % dielectric permittivity
geo_a = 0.015; geo_b = 4*lambda0; % spatial
beta = @(rad) -k0^2; % beta
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


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Total, Incident, Scattered %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CASE 3: TMz Dielectric
plotIST_fields(c, properties_c);