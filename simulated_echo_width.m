function sigma_simulated = simulated_echo_width(obj, properties, R_samp, Npts, Nmax)

    %%% Setup
    % node coordinates
    nodeX = obj.mesh.Nodes(1,:).';
    nodeY = obj.mesh.Nodes(2,:).';

    % wave propagation vector
    k0  = properties{4};

    Ei = exp(-1j*k0*nodeX); % incident field
    Es = obj.phi(:) - Ei; % scattered field

    %%% Sample scattered field on circle 
    % points on the circle:
    angle = linspace(0, 2*pi, Npts+1);
    angle = angle(1:end-1); % exclude 2pi
    x = R_samp*cos(angle); 
    y = R_samp*sin(angle);

    % interpolate real and imag components
    Fr = scatteredInterpolant(nodeX, nodeY, real(Es), 'linear', 'none');
    Fi = scatteredInterpolant(nodeX, nodeY, imag(Es), 'linear', 'none');

    % apply over the circle and recombine
    Es_c = Fr(x, y) + 1j*Fi(x, y);

    % azimuthal harmonics
    a_n = fftshift(fft(Es_c)) / Npts;

    n = (0:(Npts-1)) - floor(Npts/2);

    % truncate to Nmax
    %a_n = a_n(abs(n) <= Nmax);
    %n = n(abs(n) <= Nmax);

    % propagate to far field
    c_n = a_n ./ besselh(n, 2, k0*R_samp);

    % initialize
    sigma_simulated = zeros(size(angle));

    % calculate bistatic echo width
    for p = 1:numel(angle)
        sigma_simulated(p) = 4/k0 * abs(sum(c_n .* exp(1j*(2*n+1)*pi/4) .* exp(1j*n*angle(p))))^2;
    end
end