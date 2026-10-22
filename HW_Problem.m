clc; clear; close all;

intermediatePlots = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PEC cylinder, TM_z  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

c0 = 3*1e8; % speed of light
f0 = 10e9; % frequency
lambda0 = c0/f0; % lambda
k0 = 2*pi/lambda0; % wave vector magnitude
beta = -k0^2; % beta
alpha = 1; % alpha
geo_a = 0.015; geo_b = 3*lambda0;

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
a = a.createMesh(lambda0/20);

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
    title("\Phi(x,y)")
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
beta = -k0^2; % beta
alpha = 1; % alpha
geo_a = 0.015; geo_b = 3*lambda0;

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
b = b.createMesh(lambda0/20);

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
    title("\Phi(x,y)")
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
beta = @(rad) -k0^2*(1 + (eps_r-1)*(rad < geo_a)); % beta
alpha = 1; % alpha
geo_a = 0.015; geo_b = 3*lambda0;

% save properties
properties_c = {c0; f0; lambda0; k0; beta; alpha; beta; geo_a; geo_b};

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
c = c.createMesh(lambda0/40);

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
    title("\Phi(x,y)")
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
beta = -k0^2; % beta
alpha = @(rad) 1./(1 + (eps_r-1)*(rad < geo_a)); % alpha
geo_a = 0.015; geo_b = 3*lambda0;

% save properties
properties_d = {c0; f0; lambda0; k0; beta; alpha; beta; geo_a; geo_b};

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
d = d.createMesh(lambda0/40);

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
    title("\Phi(x,y)")
    xlabel("x [m]")
    ylabel("y [m]")
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Total, Incident, Scattered %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CASE 1: TMz PEC
[Ei_a Es_a Et_a] = plotIST_fields(a, properties_a);

% CASE 2: TEz PEC
[Ei_b Es_b Et_b] = plotIST_fields(b, properties_b);

% CASE 3: TMz Dielectric
[Ei_c Es_c Et_c] = plotIST_fields(c, properties_c);

% CASE 4: TEz Dielectric
[Ei_d Es_d Et_d] = plotIST_fields(d, properties_d);

function G = evalInterp(nodeX, nodeY, E, x, y)
    F = scatteredInterpolant(nodeX, nodeY, real(E), 'linear', 'none');
    G = F(x, y);
end

function [Ei, Es, Et] = plotIST_fields(obj, properties)
    nodeX = obj.mesh.Nodes(1,:);
    nodeY = obj.mesh.Nodes(2,:);
    R = properties{8};
    f = properties{2};
    k = 2*pi*f/3e8;
    Ei = exp(-1j*k*nodeX).';
    Et = obj.phi;
    Es = Et - Ei;
    
    % interpolate onto a grid
    x = linspace(min(nodeX), max(nodeX), 1000);
    y = linspace(min(nodeY), max(nodeY), 1000);
    [X, Y] = meshgrid(x, y);
    Gi = evalInterp(nodeX.', nodeY.', real(Ei), X, Y);
    Gt = evalInterp(nodeX.', nodeY.', real(Et), X, Y);
    Gs = evalInterp(nodeX.', nodeY.', real(Es), X, Y);

    if ~any(nodeX.^2 + nodeY.^2 < (0.95*R)^2) % if PEC
        m = X.^2 + Y.^2 < R^2;
        Gi(m) = NaN; Gt(m) = NaN; Gs(m) = NaN;
    end
    
    % plot
    fig = figure('Color','White', 'Units','normalized', 'Position', [0.02 0.04 0.96 0.5]);
    t = tiledlayout(fig, 1, 3, 'TileSpacing','compact', 'Padding','compact');
    
    % incident field
    nexttile;
    imagesc(x, y, Gi); hold on;
    plot(R*cos(linspace(0,2*pi,300)), R*sin(linspace(0,2*pi,300)), 'k', 'LineWidth', 1.5);
    set(gca,'YDir','normal', 'TickLength', [0,0], 'FontName', 'Times', 'FontSize', 18);
    set(gcf, 'Color', 'w')
    axis equal tight;
    colormap inferno;
    colorbar;
    cl = 0.50*max(abs(Gi), [], 'all');
    clim([-cl cl]);
    title("Incident");
    xlabel("x [m]")
    ylabel("y [m]")
    
    % scattered field
    nexttile;
    imagesc(x, y, Gs); hold on;
    plot(R*cos(linspace(0,2*pi,300)), R*sin(linspace(0,2*pi,300)), 'k', 'LineWidth', 1.5);
    set(gca,'YDir','normal', 'TickLength', [0,0], 'FontName', 'Times', 'FontSize', 18);
    set(gcf, 'Color', 'w')
    axis equal tight;
    colormap inferno;
    colorbar;
    cl = 0.50*max(abs(Gs), [], 'all');
    clim([-cl cl]);
    title("Scattered");
    xlabel("x [m]")
    ylabel("y [m]")
    
    % total field
    nexttile;
    imagesc(x, y, Gt); hold on;
    plot(R*cos(linspace(0,2*pi,300)), R*sin(linspace(0,2*pi,300)), 'k', 'LineWidth', 1.5);
    set(gca,'YDir','normal', 'TickLength', [0,0], 'FontName', 'Times', 'FontSize', 18);
    set(gcf, 'Color', 'w')
    axis equal tight;
    colormap inferno;
    colorbar;
    cl = 0.50*max(abs(Gt), [], 'all');
    clim([-cl cl]);
    title("Total");
    xlabel("x [m]")
    ylabel("y [m]")
    hold off;
end

%%
%%%%%%%%%%%%%%%%%%
%%% Echo Width %%%
%%%%%%%%%%%%%%%%%%