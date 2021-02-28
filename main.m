%% ability to pass args when constructing class looks like a PitA
%  do this instead
a = BEMSolver;
a.nBlades = 3;
a.nSegments =10;
a.TSR = 6;
a.nSegments = 40;
a.spacing = 'cosine';
a.nBlades = 3;
a.atol = 1e-6;
a.nIter = 100;

%% Initialise some other attributes
a = a.init();

%% run the BEM
results = a.solveRotor();