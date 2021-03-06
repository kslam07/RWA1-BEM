%% ability to pass args when constructing class looks like a PitA
%  do this instead
a = BEMSolver;
a.nBlades = 3;
a.TSR = 8;
a.nSegments = 100;
a.spacing = "cosine";
a.atol = 1e-6;
a.nIter = 100;
a.bladePitch = -2;

%% Initialise some other attributes
a = a.init();

%% run the BEM
results = a.solveRotor();


%% post-process results
figure
plot(results(:,1), results(:, 2),results(:,1), results(:, 3));
grid on

figure
plot(results(:, 1), results(:, 6));
grid on


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