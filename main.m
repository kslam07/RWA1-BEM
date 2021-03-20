%% ability to pass args when constructing class looks like a PitA
%  do this instead
solver = BEMSolver;
solver.nBlades = 3;
solver.TSR = 8;
solver.nAnnulus = 50;
solver.spacing = "0";
solver.atol = 1e-4;
solver.nIter = 100;
solver.bladePitch = -2;
solver.nPsi = 50;
solver.yawAngle = 0;
solver.uInf = 10;

%% Initialise some other attributes
solver = solver.init();

%% run the BEM
solver = solver.solveStreamtube();

solver.TSR=6;
solver=solver.init();
solver6=solver.solveStreamtube();

solver.TSR=8;
solver=solver.init();
solver8=solver.solveStreamtube();

solver.TSR=10;
solver=solver.init();
solver10=solver.solveStreamtube();

%% post-process solverults
% figure
% plot(solver.rR, rad2deg(solver.alpha), "linewidth", 1.3);
% 
% figure
% hold on
% for i = 1:solver.nPsi
%     plot(solver.rR, solver.a(:, i), solver.rR, solver.aprime(:, i), ...
%         '--', "linewidth", 1.3, "DisplayName", ...
%         ['yaw:' num2str(rad2deg(solver.psiSegment(i)))]);
% end

% figure(50)
% [~, idx] = min(abs(solver.rR - 0.9));
% plot(rad2deg(solver.psiSegment), solver.a(idx, :)); 
% legend("show")
% ylim([0 1])
% xlim([solver.rRootRatio 1])
% grid

%% Plots for report

% figure(1)
% hold on
% plot(solver6.rR, mean(solver6.alpha,2)*180/pi, solver6.rR, mean(solver6.phi,2)*180/pi,'--');
% plot(solver8.rR, mean(solver8.alpha,2)*180/pi, solver8.rR, mean(solver8.phi,2)*180/pi,'--');
% plot(solver10.rR, mean(solver10.alpha,2)*180/pi, solver10.rR, mean(solver10.phi,2)*180/pi,'--');
% xlabel('r/R (-)')
% ylabel('angle (deg)')
% legend('\alpha_{TSR=6}', "\phi_{TSR=6}",'\alpha_{TSR=8}','\phi_{TSR=8}','\alpha_{TSR=10}','\phi_{TSR=10}')
% grid on
% export_fig 'plot_alpha.png'

% figure(2)
% hold on
% plot(solver6.rR, mean(solver6.a,2), solver6.rR, mean(solver6.aprime,2),'--');
% plot(solver8.rR, mean(solver8.a,2), solver8.rR, mean(solver8.aprime,2),'--');
% plot(solver10.rR, mean(solver10.a,2), solver10.rR, mean(solver10.aprime,2),'--');
% xlabel('r/R (-)')
% ylabel('Induced velocity (-)')
% legend('a_{TSR=6}', "a'_{TSR=6}",'a_{TSR=8}', "a'_{TSR=8}",'a_{TSR=10}', "a'_{TSR=10}")
% grid on
% export_fig 'plot_a.png'

% figure(3)
% hold on
% plot(solver6.rR, mean(solver6.CT,2), solver6.rR, mean(solver6.Cq,2),'--', solver6.rR, mean(solver6.CN,2),'-.'); % coefficients
% plot(solver8.rR, mean(solver8.CT,2), solver8.rR, mean(solver8.Cq,2),'--', solver8.rR, mean(solver8.CN,2),'-.'); % coefficients
% plot(solver10.rR, mean(solver10.CT,2), solver10.rR, mean(solver10.Cq,2),'--', solver10.rR, mean(solver10.CN,2),'-.'); % coefficients
% xlabel('r/R (-)')
% ylabel('C (-)')
% legend('C_T_{,TSR=6}', 'C_Q_{,TSR=6}', 'C_N_{,TSR=6}','C_T_{,TSR=8}', 'C_Q_{,TSR=8}', 'C_N_{,TSR=8}','C_T_{,TSR=10}', 'C_Q_{,TSR=10}', 'C_N_{,TSR=10}','Location','eastoutside')
% grid on
% export_fig 'plot_ct.png'

% figure(4)
% hold on
% plot(solver6.rR, mean(solver6.Ax(:,1)*solver6.rho,2), solver6.rR, mean(solver6.Az*solver6.rho,2),'--'); %absolute forces
% plot(solver8.rR, mean(solver8.Ax(:,1)*solver8.rho,2), solver8.rR, mean(solver8.Az*solver8.rho,2),'--'); %absolute forces
% plot(solver10.rR, mean(solver10.Ax(:,1)*solver10.rho,2), solver10.rR, mean(solver10.Az*solver10.rho,2),'--'); %absolute forces
% xlabel('r/R (-)')
% ylabel('F (N)')
% legend('F_{x,TSR=6}','F_{z,TSR=6}','F_{x,TSR=8}','F_{z,TSR=8}','F_{x,TSR=10}','F_{z,TSR=10}','Location','northwest')
% grid on
% export_fig 'plot_forces.png'

% figure(5)
% t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
% nexttile
% [t,r]=meshgrid(solver6.psiSegment,solver6.rR);
% x = r.*cos(t);
% y = r.*sin(t);
% pplot = pcolor(x, y, rad2deg(solver6.alpha));
% h=colorbar;
% ylabel(h,'\alpha','Rotation',0,'FontSize',14)
% xlabel('x/R (-)')
% ylabel('y/R (-)')
% grid on
% set(pplot, "edgeColor", "none");
% colormap('default')
% title 'TSR=6'
% % colormap("summer")
% axis square
% nexttile
% [t,r]=meshgrid(solver8.psiSegment,solver8.rR);
% x = r.*cos(t);
% y = r.*sin(t);
% pplot = pcolor(x, y, rad2deg(solver8.alpha));
% h=colorbar;
% ylabel(h,'\alpha','Rotation',0,'FontSize',14)
% xlabel('x/R (-)')
% ylabel('y/R (-)')
% grid on
% set(pplot, "edgeColor", "none");
% colormap('default')
% axis square
% title 'TSR=8'
% nexttile
% [t,r]=meshgrid(solver10.psiSegment,solver10.rR);
% x = r.*cos(t);
% y = r.*sin(t);
% pplot = pcolor(x, y, rad2deg(solver10.alpha));
% h=colorbar;
% ylabel(h,'\alpha','Rotation',0,'FontSize',14)
% xlabel('x/R (-)')
% ylabel('y/R (-)')
% grid on
% set(pplot, "edgeColor", "none");
% colormap('default')
% title 'TSR=10'
% axis square
% export_fig 'plot_azimuthal.png'

% figure(6)
% plot(squeeze(sum(solver8.thrustIter(:,2,:),1))) % might need to add cutoff for plot
% xlabel('Iteration (-)')
% ylabel('T (N)')
% xlim([0,21])
% grid on
% export_fig 'plot_convergence.png'

% figure(7)
% solver.nAnnulus = 10; % influence number of annuli
% solver = solver.init();
% solver = solver.solveStreamtube();
% plot(solver.rR, mean(solver.CT,2),'o','Color','red');
% % plot(solver.rR, mean(solver.CP,2),'o','Color','red');
% hold on
% solver.nAnnulus = 50;
% solver = solver.init();
% solver = solver.solveStreamtube();
% plot(solver.rR, mean(solver.CT,2),'x','Color','black');
% hold on
% solver.nAnnulus = 100;
% solver = solver.init();
% solver = solver.solveStreamtube();
% plot(solver.rR, mean(solver.CT,2),'Color','blue');
% xlabel('r/R (-)')
% ylabel('C_T (-)')
% legend('N_{segments}=10','N_{segments}=50','N_{segments}=100','Location','south')
% grid on
% export_fig 'plot_annuli.png'

% figure(8)
% solver.nAnnulus=15;
% solver.spacing='0';
% solver = solver.init();
% solver = solver.solveStreamtube();
% plot(solver.rR, mean(solver.CT,2));
% hold on
% solver.spacing = 'cosine';
% solver = solver.init();
% solver = solver.solveStreamtube();
% plot(solver.rR, mean(solver.CT,2),'--');
% xlabel('r/R (-)')
% ylabel('C_T (-)')
% legend('Equal spacing', 'Cosine spacing','Location','south')
% grid on
% export_fig 'plot_spacing.png'

% figure(9)
% plot(solver8.rR, solver8.fTot); % correction
% xlabel('r/R (-)')
% ylabel('f (-)')
% grid on
% export_fig 'plot_tipcorr.png'

% figure(10) % made nondimensional with (np.pi*Uinf**2/(NBlades*Omega)
% % plot(solver.rR,solver.solver.*4.*(1-solver.a)./(1+solver.aprime))
% % hold on
% plot(solver8.rR,solver8.gamma)
% xlabel('r/R (-)')
% ylabel('\Gamma (-)')
% grid on
% export_fig 'plot_circ.png'

% figure(11)
% data=enthalpy(solver8);
% hold on
% plot(solver8.rR,data(1,:),'HandleVisibility','off')
% plot(solver8.rR,data(2,:),'DisplayName','Station 1 & 2')
% plot(solver8.rR,data(3,:),'HandleVisibility','off')
% plot(solver8.rR,data(4,:),'DisplayName','Station 3')
% plot(solver8.rR,data(5,:),'DisplayName','Station 4')
% legend()
% xlabel('r/R (-)')
% ylabel('h_s (J/kg)')
% ylim([min(min(data(4,:)))*0.99,max(data(:,1))*1.01])
% export_fig 'plot_enthalpy.png'

%%
% figure(1)
% plot(solver.rR, mean(solver.alpha,2)*180/pi, solver.rR, mean(solver.phi,2)*180/pi);
% xlabel('r/R (-)')
% ylabel('angle (deg)')
% legend('\alpha', "\phi")
% grid on

% figure(2)
% plot(solver.rR, mean(solver.a,2), solver.rR, mean(solver.aprime,2));
% xlabel('r/R (-)')
% ylabel('Induced velocity (-)')
% legend('a', "a'")
% grid on

% figure(3)
% plot(solver.rR, mean(solver.CT,2), solver.rR, mean(solver.Cq,2), solver.rR, mean(solver.CN,2)); % coefficients
% xlabel('r/R (-)')
% ylabel('C (-)')
% legend('C_T', 'C_Q', 'C_N')
% grid on

% figure(4)
% plot(solver.rR, mean(solver.Ax(:,1)*solver.rho,2), solver.rR, mean(solver.Az*solver.rho,2)); %absolute forces
% xlabel('r/R (-)')
% ylabel('F (N)')
% legend('F_x','F_z','Location','northwest')
% grid on

% figure(5)
% [t,r]=meshgrid(solver.psiSegment,solver.rR);
% x = r.*cos(t);
% y = r.*sin(t);
% pplot = pcolor(x, y, rad2deg(solver.alpha));
% h=colorbar;
% ylabel(h,'\alpha','Rotation',0,'FontSize',14)
% xlabel('x/R (-)')
% ylabel('y/R (-)')
% grid on
% set(pplot, "edgeColor", "none");
% colormap('default')
% % colormap("summer")

% figure(6)
% plot(squeeze(sum(solver.thrustIter(:,2,:),1))) % might need to add cutoff for plot
% xlabel('Iteration (-)')
% ylabel('T (N)')
% grid on

% figure(7)
% solver = BEMSolver;
% solver.nAnnulus = 10; % influence number of annuli
% solver = solver.init();
% solver = solver.solveStreamtube();
% plot(solver.rR, mean(solver.CT,2),'o','Color','red');
% hold on
% solver = BEMSolver;
% solver.nAnnulus = 50;
% solver = solver.init();
% solver = solver.solveStreamtube();
% plot(solver.rR, mean(solver.CT,2),'x','Color','black');
% hold on
% solver = BEMSolver;
% solver.nAnnulus = 100;
% solver = solver.init();
% solver = solver.solveStreamtube();
% plot(solver.rR, mean(solver.CT,2),'Color','blue');
% xlabel('r/R (-)')
% ylabel('C_T (-)')
% legend('N_{segments}=10','N_{segments}=50','N_{segments}=100','Location','south')
% grid on

% figure(8)
% solver.spacing='0';
% solver = solver.init();
% solver = solver.solveStreamtube();
% plot(solver.rR, mean(solver.CT,2));
% hold on
% solver.spacing = 'cosine';
% solver = solver.init();
% solver = solver.solveStreamtube();
% plot(solver.rR, mean(solver.CT,2));
% xlabel('r/R (-)')
% ylabel('C_T (-)')
% legend('Equal spacing', 'Cosine spacing','Location','south')
% grid on

% figure(9)
% plot(solver.rR, solver.fTot); % correction
% xlabel('r/R (-)')
% ylabel('f (-)')
% grid on

% figure(10) % made nondimensional with (np.pi*Uinf**2/(NBlades*Omega)
% % plot(solver.rR,solver.solver.*4.*(1-solver.a)./(1+solver.aprime))
% % hold on
% plot(solver.rR,solver.gamma,'x')
% xlabel('r/R (-)')
% ylabel('\Gamma (-)')
% grid on

% p=101325;
% rho=1.225;
% pt=p+0.5*rho*solver.uInf^2;
% v2=solver.uInf*(1-mean(solver.a,2));
% p2=pt-0.5*rho*v2.^2;
% 
% pressureJump=sum(solver.Ax*solver.rho,2)./solver.areaAnnulus;
% v4=solver.uInf*(1-2*mean(solver.a,2));
% p3=p2-pressureJump;
% pt3=p3+0.5*rho*v2.^2;
% p4=pt3-0.5*rho*v4.^2.
% 
% enthalpy1=ones(1,solver.nAnnulus)*(p/rho+solver.uInf^2/2);
% enthalpy2=mean(ones(1,solver.nAnnulus).*(p2/rho+v2.^2/2),2);
% enthalpy3=mean(ones(1,solver.nAnnulus).*(p3/rho+v2.^2/2),2);
% enthalpy4=mean(ones(1,solver.nAnnulus).*(p4/rho+v4.^2/2),2);
% enthalpy3mean=mean(mean(enthalpy3,2));

% data=enthalpy(solver);
% 
% figure(11)
% hold on
% plot(solver.rR,data(1,:),'HandleVisibility','off')
% plot(solver.rR,data(2,:),'DisplayName','Station 1 & 2')
% plot(solver.rR,data(3,:),'DisplayName','Station 3')
% plot(solver.rR,data(4,:),'HandleVisibility','off')
% plot(solver.rR,data(5,:),'DisplayName','Station 4')
% legend()
% xlabel('r/R (-)')
% ylabel('h_s (J/kg)')
% ylim([min(min(data(4,:)))*0.99,max(data(:,1))*1.01])


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