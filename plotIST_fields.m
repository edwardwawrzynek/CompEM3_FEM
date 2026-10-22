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