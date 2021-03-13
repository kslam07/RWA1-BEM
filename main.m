%% ability to pass args when constructing class looks like a PitA
%  do this instead
solver = BEMSolver;
solver.nBlades = 5;
solver.TSR = 12;
solver.nAnnulus = 50;
solver.spacing = "cosine";
solver.atol = 1e-4;
solver.nIter = 100;
solver.bladePitch = -2;
solver.nPsi = 50;
solver.yawAngle = 30;

%% Initialise some other attributes
solver = solver.init();

%% run the BEM
solver = solver.solveStreamtube();

%% post-process solverults

figure
plot(solver.rR, rad2deg(solver.alpha), "linewidth", 1.3);

figure
hold on
for i = 1:solver.nPsi
    plot(solver.rR, solver.a(:, i), solver.rR, solver.aprime(:, i), ...
        '--', "linewidth", 1.3, "DisplayName", ...
        ['yaw:' num2str(rad2deg(solver.psiSegment(i)))]);
end

% legend("show")
% ylim([0 1])
% xlim([solver.rRootRatio 1])
% grid
% 
% figure(1)
% plot(solver.rR(:,1), solver.alpha(:,1)*180/pi, solver.rR(:,1), solver.phi(:,1)*180/pi);
% xlabel('r/R (-)')
% ylabel('angle (deg)')
% legend('\alpha', "\phi")
% grid on
% 
% figure(2)
% plot(solver.rR, solver.a, solver.rR, solver.aprime);
% xlabel('r/R (-)')
% ylabel('Induced velocity (-)')
% legend('a', "a'")
% grid on
% 
figure(3)
plot(solver.rR, solver.CT)%, solver.rR, solver.Cq, solver.rR, solver.CN); % coefficients
xlabel('r/R (-)')
ylabel('C (-)')
legend('C_T')%, 'C_Q', 'C_N')
grid on

% figure(4)
% plot(solver.rR(:,1), solver.Ax(:,1), solver.rR, solver.Az); %absolute forces
% xlabel('r/R (-)')
% ylabel('F (N)')
% legend('F_x','F_z')
% grid on
% 
figure(5)
[t,r]=meshgrid(solver.psiSegment,solver.rR);
x = r.*cos(t);
y = r.*sin(t);
contourf(x, y, solver.a)
h=colorbar;
ylabel(h,'\psi','Rotation',0,'FontSize',14)
xlabel('x/R (-)')
ylabel('y/R (-)')
grid on
% colormap("summer")
% 
% figure(6)
% plot(squeeze(sum(solver.thrustIter(:,2,:),1))) % might need to add cutoff for plot
% xlabel('Iteration (-)')
% ylabel('T (N)')
% grid on
% 
% figure(7)
% solver = BEMSolver;
% solver.nAnnulus = 10; % influence number of annuli
% solver = solver.init();
% solver = solver.solveStreamtube();
% plot(solver.rR, solver.CT,'o');
% hold on
% solver = BEMSolver;
% solver.nAnnulus = 50;
% solver = solver.init();
% solver = solver.solveStreamtube();
% plot(solver.rR, solver.CT,'x');
% hold on
% solver = BEMSolver;
% solver.nAnnulus = 100;
% solver = solver.init();
% solver = solver.solveStreamtube();
% plot(solver.rR, solver.CT);
% xlabel('r/R (-)')
% ylabel('C_T (-)')
% legend('N_{segments}=10','N_{segments}=50','N_{segments}=100','Location','south')
% grid on
% 
% figure(8)
% solver.spacing='0';
% solver = solver.init();
% solver = solver.solveStreamtube();
% plot(solver.rR, solver.CT);
% hold on
% solver.spacing = 'cosine';
% solver = solver.init();
% solver = solver.solveStreamtube();
% plot(solver.rR, solver.CT);
% xlabel('r/R (-)')
% ylabel('C_T (-)')
% legend('Equal spacing', 'Cosine spacing','Location','south')
% grid on
% 
% figure(9)
% plot(solver.rR, solver.fTot); % correction
% xlabel('r/R (-)')
% ylabel('f (-)')
% grid on
% 
% figure(10) % made nondimensional with (np.pi*Uinf**2/(NBlades*Omega)
% % plot(solver.rR,solver.solver.*4.*(1-solver.a)./(1+solver.aprime))
% % hold on
% plot(solver.rR,solver.gamma,'x')
% xlabel('r/R (-)')
% ylabel('\Gamma (-)')
% grid on

%% aSkew
% nPsi = linspace(0, 2*pi, solver.nPsi);
% x = solver.rR' .* cos(nPsi);
% y = solver.rR' .* sin(nPsi);
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
%     CT = 4.*solver.*(1-a);
%     if glauert
%         CT1=1.816;
%         a1=1-sqrt(CT1)/2;
%         CT(a>a1) = CT1-4*(sqrt(CT1)-1)*(1-a(a>a1));
%     end
% end