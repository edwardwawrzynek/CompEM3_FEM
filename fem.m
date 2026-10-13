a = FEMProblem(1, 1, 1, 1, 1, 0.2, 1.0);
a = a.createMesh(0.05);
figure;
a.plotGeo();
figure;
a.plotMesh(false);
figure;
a.assemble();