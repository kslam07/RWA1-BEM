%% ability to pass args when constructing class looks like a PitA
%  do this instead
a = BEMSolver;
a.nBlades = 3;
a.TSR = 6;
a.nSegments = 100;
a.spacing = "cosine";
a.atol = 1e-6;
a.nIter = 100;
a.bladePitch = -2;
a.nPsi = 50;

%% Initialise some other attributes
a = a.init();

%% run the BEM
res = a.solveRotor();

%% post-process results
figure
plot(res.rR, res.a, res.rR, res.aprime);
grid on

% figure
% plot(res.rR, res.fTot);
% grid on

% aSkew
% nPsi = linspace(0, 2*pi, a.nPsi);
% x = res.rR' .* cos(nPsi);
% y = res.rR' .* sin(nPsi);
% 
% [X, Y] = meshgrid(x, y);
% figure

% as = linspace(-0.5, 1, 50);
% CT = CTfunction(as, true);
% a2 = calcGlauertCorr(as);

% figure
% plot(a2, CT);

% function a = calcGlauertCorr(CT)
%     % CHECKED
%     % computes the Glauert correction for heavily loaded rotors
%     if CT < (2*sqrt(1.816)-1.816)
%         a = 0.5 - sqrt(1-CT)/2;
%     else
%         a = 1+(CT-1.816)/(4*sqrt(1.816)-4);
%     end
% end
% 
% function CT = CTfunction(a, glauert)
%     CT = 4.*a.*(1-a);
%     if glauert
%         CT1=1.816;
%         a1=1-sqrt(CT1)/2;
%         CT(a>a1) = CT1-4*(sqrt(CT1)-1)*(1-a(a>a1));
%     end
% end