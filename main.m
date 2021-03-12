%% ability to pass args when constructing class looks like a PitA
%  do this instead
solver = BEMSolver;
solver.nBlades = 3;
solver.TSR = 6;
solver.nAnnulus = 100;
solver.spacing = "0";
solver.atol = 1e-6;
solver.nIter = 100;
solver.bladePitch = -2;
solver.nPsi = 10;
solver.yawAngle = 15;

%% Initialise some other attributes
solver = solver.init();

%% run the BEM
solver = solver.solveStreamtube();

%% post-process results
% figure
% plot(res.rR, res.alpha, res.rR, res.phi);

figure
hold on
for i = 1:solver.nPsi
    plot(solver.rR, solver.a(:, i), solver.rR, solver.aprime(:, i), ...
        '--', "linewidth", 1.3, "DisplayName", ...
        ['yaw:' num2str(rad2deg(solver.psiSegment(i)))]);
end

legend("show")
ylim([0 1])
xlim([solver.rRootRatio 1])
grid

% plot(res.rR, res.CT, res.rR, res.CQ, res.rR, res.CN); % coefficients
% legend('CT', 'CQ', 'CN')

% plot(res.rR, res.Ax, res.rR, res.Az); %absolute forces
% grid on

% [t,r]=meshgrid(res.psi,res.rR);
% x = r.*cos(t);
% y = r.*sin(t);
% contourf(x, y, res.askew)
% colorbar

% plot(sum(res.ThrustIter,2))

% figure
% plot(solver.rR, solver.fTot(:, 1)); % crrection
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