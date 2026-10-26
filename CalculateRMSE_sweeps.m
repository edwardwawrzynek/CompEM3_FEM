clc; clear; close all;

% choose sweep variable
% 1 - h, spatial sampling
% 2 - truncation distance from scatterer
% 3 - order of the ABC
sweep_var = 2;

if sweep_var == 1
    % h-sweep
    sweep = [1/10];
end

if sweep_var == 2
    % trunc-sweep
    sweep = 3:1:7;
end

if sweep_var == 3
    % order of ABC
    sweep = 1:3;
end

% initialize
RMSE_TM_PEC = cell(1, numel(sweep));
RMSE_TE_PEC = cell(1, numel(sweep));
RMSE_TM_Diel = cell(1, numel(sweep));
RMSE_TE_Diel = cell(1, numel(sweep));

for k = 1:numel(sweep)
    
    disp(k);

    if sweep_var == 1
        % for h-sweep
        mesh_h = @(lambda) lambda * sweep(k);
    else
        mesh_h = @(lambda) lambda * 1/100;
    end

    if sweep_var == 2
        % for trunc-sweep
        trunc = @(lambda) lambda * sweep(k);
    else
        trunc = @(lambda) lambda * 3;
    end
    
    disp("TM");
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PEC cylinder, TM_z  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    c0 = 3*1e8; % speed of light
    f0 = 10e9; % frequency
    lambda0 = c0/f0; % lambda
    k0 = 2*pi/lambda0; % wave vector magnitude
    geo_a = 0.015; geo_b = trunc(lambda0); % spatial
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
    
    a = a.assemble();
    
    disp("TE");
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PEC cylinder, TE_z  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    c0 = 3*1e8; % speed of light
    f0 = 10e9; % frequency
    lambda0 = c0/f0; % lambda
    k0 = 2*pi/lambda0; % wave vector magnitude
    geo_a = 0.015; geo_b = trunc(lambda0); % spatial
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
    
    b = b.assemble();
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Dielectric cylinder, TM_z  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    c0 = 3*1e8; % speed of light
    f0 = 10e9; % frequency
    lambda0 = c0/f0; % lambda
    k0 = 2*pi/lambda0; % wave vector magnitude
    eps_r = 9; % dielectric permittivity
    geo_a = 0.015; geo_b = trunc(lambda0); % spatial
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
    
    c = c.assemble();
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Dielectric cylinder, TE_z  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    c0 = 3*1e8; % speed of light
    f0 = 10e9; % frequency
    lambda0 = c0/f0; % lambda
    k0 = 2*pi/lambda0; % wave vector magnitude
    eps_r = 9; % dielectric permittivity
    geo_a = 0.015; geo_b = trunc(lambda0); % spatial
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
    
    d = d.assemble();
    
    %%%%%%%%%%%%%%%%%%
    %%% Echo Width %%%
    %%%%%%%%%%%%%%%%%%
    
    % sweep the angle
    angle = linspace(0, 2*pi, 360);
    
    % CASE 1: TMz PEC
    k0 = properties_a{4}; R = properties_a{8};
    sigma_simulated  = simulated_echo_width(a, properties_a, 0.9*geo_b, 360, 10);
    sigma_analytical = bistatic_echo_width_tm_pec(angle, k0, R);
    RMSE_TM_PEC{k} = sqrt(sum((sigma_simulated(:) - sigma_analytical(:)).^2)) ./ sqrt(sum(sigma_analytical(:).^2));
    
    % CASE 2: TEz PEC
    k0 = properties_b{4}; R = properties_b{8};
    sigma_simulated  = simulated_echo_width(b, properties_b, 0.9*geo_b, 360, 10);
    sigma_analytical = bistatic_echo_width_te_pec(angle, k0, R);
    RMSE_TE_PEC{k} = sqrt(sum((sigma_simulated(:) - sigma_analytical(:)).^2)) ./ sqrt(sum(sigma_analytical(:).^2));
    
    % CASE 3: TMz Dielectric
    k0 = properties_c{4}; R = properties_c{8}; eps_r = properties_c{10};
    sigma_simulated  = simulated_echo_width(c, properties_c, 0.9*geo_b, 360, 10);
    sigma_analytical = bistatic_echo_width_tm_eps(angle, k0, R, eps_r);
    RMSE_TM_Diel{k} = sqrt(sum((sigma_simulated(:) - sigma_analytical(:)).^2)) ./ sqrt(sum(sigma_analytical(:).^2));
    
    % CASE 4: TEz Dielectric
    k0 = properties_d{4}; R = properties_d{8}; eps_r = properties_d{10};
    sigma_simulated  = simulated_echo_width(d, properties_d, 0.05, 360, 10);
    sigma_analytical = bistatic_echo_width_te_eps(angle, k0, R, eps_r);
    RMSE_TE_Diel{k} = sqrt(sum((sigma_simulated(:) - sigma_analytical(:)).^2)) ./ sqrt(sum(sigma_analytical(:).^2));
end

%% 
%%%%%%%%%%%%%%%%%
%%% Plot RMSE %%%
%%%%%%%%%%%%%%%%%

if sweep_var == 1
    figure(); hold on;
    plot(sweep, cell2mat(RMSE_TM_PEC), 'LineWidth', 2);
    plot(sweep, cell2mat(RMSE_TE_PEC), ":", 'LineWidth', 2);
    plot(sweep, cell2mat(RMSE_TM_Diel), "--", 'LineWidth', 2);
    plot(sweep, cell2mat(RMSE_TE_Diel), "-.", 'LineWidth', 2);
    set(gca,'FontName', 'Times', 'FontSize', 18);
    set(gca,'XScale', 'log', 'YScale', 'log');
    ylim([1e-2,1e1]);
    set(gcf, 'Color', 'w');
    legend(["TM_z PEC", "TE_z PEC", "TM_z Dielectric", "TE_z Dielectric"], 'Location', 'northwest');
    title("Bistatic Echowidth RMSE")
    xlabel("h/\lambda_0");
    ylabel("RMSE");
    grid on; box on;
end

if sweep_var == 2
    figure(); hold on;
    plot(sweep, cell2mat(RMSE_TM_PEC), 'LineWidth', 2);
    plot(sweep, cell2mat(RMSE_TE_PEC), ":", 'LineWidth', 2);
    plot(sweep, cell2mat(RMSE_TM_Diel), "--", 'LineWidth', 2);
    plot(sweep, cell2mat(RMSE_TE_Diel), "-.", 'LineWidth', 2);
    set(gca,'FontName', 'Times', 'FontSize', 18);
    set(gca,'XScale', 'log', 'YScale', 'log');
    set(gcf, 'Color', 'w');
    legend(["TM_z PEC", "TE_z PEC", "TM_z Dielectric", "TE_z Dielectric"], 'Location', 'southeast');
    title("Bistatic Echowidth RMSE")
    xlabel("Truncation distance [\lambda_0]");
    ylabel("RMSE");
    grid on; box on;
end

if sweep_var == 3
    figure(); hold on;
    plot(sweep, cell2mat(RMSE_TM_PEC), 'LineWidth', 2);
    plot(sweep, cell2mat(RMSE_TE_PEC), ":", 'LineWidth', 2);
    plot(sweep, cell2mat(RMSE_TM_Diel), "--", 'LineWidth', 2);
    plot(sweep, cell2mat(RMSE_TE_Diel), "-.", 'LineWidth', 2);
    set(gca,'FontName', 'Times', 'FontSize', 18);
    set(gca,'XScale', 'log', 'YScale', 'log');
    set(gcf, 'Color', 'w');
    legend(["TM_z PEC", "TE_z PEC", "TM_z Dielectric", "TE_z Dielectric"], 'Location', 'southeast');
    title("Bistatic Echowidth RMSE")
    xlabel("Order of the ABC");
    ylabel("RMSE");
    grid on; box on;
end