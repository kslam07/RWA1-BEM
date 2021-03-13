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

% figure(1)
% plot(res.rR, res.alpha*180/pi, res.rR, res.phi*180/pi);
% xlabel('r/R (-)')
% ylabel('angle (deg)')
% legend('\alpha', "\phi")
% grid on
% 
% figure(2)
% plot(res.rR, res.a, res.rR, res.aprime);
% xlabel('r/R (-)')
% ylabel('Induced velocity (-)')
% legend('a', "a'")
% grid on
% 
% figure(3)
% plot(res.rR, res.CT, res.rR, res.CQ, res.rR, res.CN); % coefficients
% xlabel('r/R (-)')
% ylabel('C (-)')
% legend('C_T', 'C_Q', 'C_N')
% grid on
% 
% figure(4)
% plot(res.rR, res.Ax, res.rR, res.Az); %absolute forces
% xlabel('r/R (-)')
% ylabel('F (N)')
% legend('F_x','F_z')
% grid on
% 
% figure(5)
% [t,r]=meshgrid(res.psi,res.rR);
% x = r.*cos(t);
% y = r.*sin(t);
% contourf(x, y, res.askew)
% h=colorbar;
% ylabel(h,'      \psi','Rotation',0,'FontSize',14)
% xlabel('x/R (-)')
% ylabel('y/R (-)')
% grid on
% 
% figure(6)
% plot(sum(res.ThrustIter,1)) % might need to add cutoff for plot
% xlabel('Iteration (-)')
% ylabel('T (N)')
% grid on
% 
% figure(7)
% a.nSegments = 10; % influence number of annuli
% a = a.init();
% res = a.solveRotor();
% plot(res.rR, res.CT,'o');
% hold on
% a.nSegments = 50;
% a = a.init();
% res = a.solveRotor();
% plot(res.rR, res.CT,'x');
% hold on
% a.nSegments = 100;
% a = a.init();
% res = a.solveRotor();
% plot(res.rR, res.CT);
% xlabel('r/R (-)')
% ylabel('C_T (-)')
% legend('N_{segments}=10','N_{segments}=50','N_{segments}=100','Location','south')
% grid on
% 
% figure(8)
% a.spacing='0';
% a = a.init();
% res = a.solveRotor();
% plot(res.rR, res.CT);
% hold on
% a.spacing = 'cosine';
% a = a.init();
% res = a.solveRotor();
% plot(res.rR, res.CT);
% xlabel('r/R (-)')
% ylabel('C_T (-)')
% legend('Equal spacing', 'Cosine spacing','Location','south')
% grid on
% 
% figure(9)
% plot(res.rR, res.fTot); % correction
% xlabel('r/R (-)')
% ylabel('f (-)')
% grid on
% 
% figure(10) % made nondimensional with (np.pi*Uinf**2/(NBlades*Omega)
% % plot(res.rR,res.a.*4.*(1-res.a)./(1+res.aprime))
% % hold on
% plot(res.rR,res.gamma,'x')
% xlabel('r/R (-)')
% ylabel('\Gamma (-)')
% grid on

%% aSkew
% nPsi = linspace(0, 2*pi, a.nPsi);
% x = res.rR' .* cos(nPsi);
% y = res.rR' .* sin(nPsi);
% 
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