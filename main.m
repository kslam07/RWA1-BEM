%% ability to pass args when constructing class looks like a PitA
%  do this instead
a = BEMSolver;
a.nBlades = 3;
a.nSegments =10;
a.TSR = 6;
a.nSegments = 50;
a.spacing = '0';
a.nBlades = 3;
a.atol = 1e-6;
a.nIter = 100;
a.bladePitch = 2;

%% Initialise some other attributes
a = a.init();

%% run the BEM
results = a.solveRotor();

%% post-process results
figure
plot(results(:,1), results(:, 2), results(:,1), results(:,3));
grid on