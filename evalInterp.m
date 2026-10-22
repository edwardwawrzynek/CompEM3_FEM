function G = evalInterp(nodeX, nodeY, E, x, y)
    F = scatteredInterpolant(nodeX, nodeY, real(E), 'linear', 'none');
    G = F(x, y);
end