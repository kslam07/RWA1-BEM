%% ability to pass args when constructing class looks like a PitA
%  do this instead
solver = BEMSolver;
solver.nBlades = 3;
solver.TSR = 8;
solver.nAnnulus = 100;
solver.spacing = "0";
solver.atol = 1e-4;
solver.nIter = 100;
solver.bladePitch = -2;
solver.nPsi = 100;
solver.yawAngle = 0;
solver.uInf = 10;

%% Initialise some other attributes
solver = solver.init();

%% run the BEM
solver = solver.solveStreamtube();

%% General plots
% figure(1)
% plot(squeeze(sum(solver.thrustIter(:,2,:),1))) % might need to add cutoff for plot
% xlabel('Iteration (-)')
% ylabel('T (N)')
% xlim([0,21])
% grid on
% set(gcf,'color','w')
% export_fig 'plot_convergence.png'

% figure(2)
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

% figure(3)
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
% set(gcf,'color','w')
% export_fig 'plot_spacing.png'

% solver = BEMSolver;
% solver.nBlades = 3;
% solver.TSR = 8;
% solver.nAnnulus = 50;
% solver.spacing = "0";
% solver.atol = 1e-4;
% solver.nIter = 100;
% solver.bladePitch = -2;
% solver.nPsi = 50;
% solver.yawAngle = 0;
% solver.uInf = 10;
% solver = solver.init();
% solver = solver.solveStreamtube();

% figure(4)
% plot(solver.rR, solver.fTot); % correction
% xlabel('r/R (-)')
% ylabel('f (-)')
% grid on
% set(gcf,'color','w')
% export_fig 'plot_tipcorr.png'

% figure(5) % made nondimensional with (np.pi*Uinf**2/(NBlades*Omega)
% % plot(solver.rR,solver.solver.*4.*(1-solver.a)./(1+solver.aprime))
% % hold on
% plot(solver.rR,solver.gamma)
% xlabel('r/R (-)')
% ylabel('\Gamma (-)')
% grid on
% set(gcf,'color','w')
% export_fig 'plot_circ.png'

figure("defaultAxesFontSize", 18)
[data, R1, R2, R3, R4] =enthalpy(solver);
hold on
plot(R1*solver.rR,data(1,:),'DisplayName', "Station 1", "Marker", 'o', ...
    "markerIndices", 1:5:length(solver.rR))
plot(R2*solver.rR,data(2,:),'DisplayName','Station 2')
plot(R3*solver.rR,data(3,:),'DisplayName','Station 3')
plot(R4*solver.rR,data(4,:),'DisplayName','Station 4')
% plot(R4*solver.rR,data(5,:),'DisplayName','Station 4')
legend()
xlabel('R (-)')
ylabel('h_s (J/kg)')
ylim([min(min(data(4,:)))*0.99,max(data(:,1))*1.01])
grid on
set(gcf,'color','w')
export_fig 'plot_enthalpy.png'

%% Plots for TSR

solver.TSR=6;
solver=solver.init();
solver6=solver.solveStreamtube();

solver.TSR=8;
solver=solver.init();
solver8=solver.solveStreamtube();

solver.TSR=10;
solver=solver.init();
solver10=solver.solveStreamtube();

% Total Thrust
% sum(mean(solver6.Ax,2))
% sum(mean(solver8.Ax,2))
% sum(mean(solver10.Ax,2))

% Total torque
% sum(mean(solver6.Az*solver6.rR*50,2))
% sum(mean(solver8.Az*solver8.rR*50,2))
% sum(mean(solver10.Az*solver10.rR*50,2))

% figure(1)
% hold on
% plot(solver6.rR, mean(solver6.alpha,2)*180/pi, solver6.rR, mean(solver6.phi,2)*180/pi,'--');
% plot(solver8.rR, mean(solver8.alpha,2)*180/pi, solver8.rR, mean(solver8.phi,2)*180/pi,'--');
% plot(solver10.rR, mean(solver10.alpha,2)*180/pi, solver10.rR, mean(solver10.phi,2)*180/pi,'--');
% xlabel('r/R (-)')
% ylabel('angle (deg)')
% legend('\alpha_{TSR=6}', "\phi_{TSR=6}",'\alpha_{TSR=8}','\phi_{TSR=8}','\alpha_{TSR=10}','\phi_{TSR=10}')
% grid on
% set(gcf,'color','w')
% export_fig 'TSR_Alpha.png'

% figure(2)
% hold on
% plot(solver6.rR, mean(solver6.a,2), solver6.rR, mean(solver6.aprime,2),'--');
% plot(solver8.rR, mean(solver8.a,2), solver8.rR, mean(solver8.aprime,2),'--');
% plot(solver10.rR, mean(solver10.a,2), solver10.rR, mean(solver10.aprime,2),'--');
% xlabel('r/R (-)')
% ylabel('Induced velocity factor (-)')
% legend('a_{TSR=6}', "a'_{TSR=6}",'a_{TSR=8}', "a'_{TSR=8}",'a_{TSR=10}', "a'_{TSR=10}")
% grid on
% set(gcf,'color','w')
% export_fig 'TSR_a.png'

% figure(3)
% hold on
% plot(solver6.rR, mean(solver6.CT,2), solver6.rR, mean(solver6.Cq,2),'--', solver6.rR, mean(solver6.CN,2),'-.'); % coefficients
% plot(solver8.rR, mean(solver8.CT,2), solver8.rR, mean(solver8.Cq,2),'--', solver8.rR, mean(solver8.CN,2),'-.'); % coefficients
% plot(solver10.rR, mean(solver10.CT,2), solver10.rR, mean(solver10.Cq,2),'--', solver10.rR, mean(solver10.CN,2),'-.'); % coefficients
% xlabel('r/R (-)')
% ylabel('C (-)')
% legend('C_T_{,TSR=6}', 'C_Q_{,TSR=6}', 'C_N_{,TSR=6}','C_T_{,TSR=8}', 'C_Q_{,TSR=8}', 'C_N_{,TSR=8}','C_T_{,TSR=10}', 'C_Q_{,TSR=10}', 'C_N_{,TSR=10}','Location','eastoutside')
% grid on
% set(gcf,'color','w')
% export_fig 'TSR_C.png'

% figure(4)
% hold on
% plot(solver6.rR, mean(solver6.Ax(:,1)*solver6.rho,2), solver6.rR, mean(solver6.Az*solver6.rho,2),'--'); %absolute forces
% plot(solver8.rR, mean(solver8.Ax(:,1)*solver8.rho,2), solver8.rR, mean(solver8.Az*solver8.rho,2),'--'); %absolute forces
% plot(solver10.rR, mean(solver10.Ax(:,1)*solver10.rho,2), solver10.rR, mean(solver10.Az*solver10.rho,2),'--'); %absolute forces
% xlabel('r/R (-)')
% ylabel('F (N)')
% legend('F_{ax,TSR=6}','F_{t,TSR=6}','F_{ax,TSR=8}','F_{t,TSR=8}','F_{ax,TSR=10}','F_{t,TSR=10}','Location','northwest')
% grid on
% set(gcf,'color','w')
% export_fig 'TSR_F.png'

% plot CL as function of chord distribution

% compute chord length at each rR
% chordTSR = solver8.computeChordLength(solver6.rR);
% polarData = table2array(readtable("polar_DU95W180.xlsx"));
% LDpoint = polarData(polarData(:,2) == max(polarData(:,2)), :);
% 
% figure("defaultAxesFontSize", 18)
% tiledlayout(1,2, "TileSpacing", "compact",'Padding', 'none')
% ax1=nexttile;
% hold on
% plot(chordTSR, solver6.Cl(:,1), "linewidth", 1.5)
% plot(chordTSR, solver8.Cl(:,1), "linewidth", 1.5)
% plot(chordTSR, solver10.Cl(:,1), "linewidth", 1.5)
% legend("$\lambda=6$", "$\lambda=8$", "$\lambda=10$", "interpreter", "latex")
% set(gca, 'XDir','reverse')
% xlabel("chord length [m]")
% ylabel("$C_{L}$", "interpreter", "latex")
% grid
% ax2=nexttile;
% hold on 
% plot(polarData(:,3), polarData(:,2), "linewidth", 1.5)
% scatter(LDpoint(3), LDpoint(2), 50, "Marker", 'x', "lineWidth", 2.5)
% legend("Airfoil polar", "Max. L/D")
% grid
% xlabel("$C_{D}$ [-]", "interpreter", "latex")
% ylabel("$C_{L}$", "interpreter", "latex")
% set(gcf,'color','w');
% linkaxes([ax1 ax2], "y")

%% Plots for Skew

solver = BEMSolver;
solver.nBlades = 3;
solver.TSR = 8;
solver.nAnnulus = 100;
solver.spacing = "0";
solver.atol = 1e-4;
solver.nIter = 100;
solver.bladePitch = -2;
solver.nPsi = 100;
solver.yawAngle = 0;
solver.uInf = 10;
solver = solver.init();
solverS0 = solver.solveStreamtube();

solver.yawAngle=15;
solver=solver.init();
solverS15=solver.solveStreamtube();

solver.yawAngle=30;
solver=solver.init();
solverS30=solver.solveStreamtube();

[t,r]=meshgrid(solverS0.psiSegment,solverS0.rR);
x = r.*sin(t);
y = r.*cos(t);

% Total thrust
% sum(mean(solverS0.Ax,2))
% sum(mean(solverS15.Ax,2))
% sum(mean(solverS30.Ax,2))

% Total torque
% sum(mean(solverS0.Az*solverS0.rR*50,2))
% sum(mean(solverS15.Az*solverS15.rR*50,2))
% sum(mean(solverS30.Az*solverS30.rR*50,2))

% figure("defaultAxesFontSize", 18)
% minVal = min([solverS0.alpha(:); solverS15.alpha(:); solverS30.alpha(:)]);
% maxVal = max([solverS0.alpha(:); solverS15.alpha(:); solverS30.alpha(:)]);
% t1 = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
% nexttile
% pplot = contourf(x, y, rad2deg(solverS0.alpha));
% h=colorbar;
% caxis manual
% caxis([rad2deg(minVal), rad2deg(maxVal)])
% xlabel(h,'\alpha','Rotation',0,'FontSize',18)
% xlabel('x/R (-)', 'FontSize',17)
% ylabel('y/R (-)', 'FontSize',17)
% grid on
% colormap('default')
% title('\gamma = 0', 'FontSize',17)
% axis square
% nexttile
% pplot = contourf(x, y, rad2deg(solverS15.alpha));
% h=colorbar;
% caxis manual
% caxis([rad2deg(minVal), rad2deg(maxVal)])
% ylabel(h,'\alpha','Rotation',0,'FontSize',18)
% xlabel('x/R (-)', 'FontSize',17)
% grid on
% colormap('default')
% axis square
% title('\gamma = 15', 'FontSize',17)
% nexttile
% caxis manual
% caxis([rad2deg(minVal), rad2deg(maxVal)])
% pplot = contourf(x, y, rad2deg(solverS30.alpha));
% h=colorbar;
% ylabel(h,'\alpha','Rotation',0,'FontSize',18)
% xlabel('x/R (-)', 'FontSize',17)
% grid on
% colormap('default')
% title('\gamma = 30', 'FontSize',16)
% axis square
% set(gcf, 'Position', get(0, 'Screensize'));
% set(gcf,'color','w')
% export_fig 'SKEW_angles.png'
% 
% figure("defaultAxesFontSize", 18)
% minVal = min([solverS0.phi(:); solverS15.phi(:); solverS30.phi(:)]);
% maxVal = max([solverS0.phi(:); solverS15.phi(:); solverS30.phi(:)]);
% t1 = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
% nexttile
% pplot = contourf(x, y, rad2deg(solverS0.phi));
% h=colorbar;
% caxis manual
% caxis([rad2deg(minVal), rad2deg(maxVal)])
% ylabel(h,'\phi','Rotation',0,'FontSize',14)
% xlabel('x/R (-)', 'FontSize',17)
% ylabel('y/R (-)', 'FontSize',17)
% grid on
% % set(pplot, "edgeColor", "none");
% colormap('default')
% title('\gamma = 0', 'FontSize',17)
% % colormap("summer")
% axis square
% nexttile
% pplot = contourf(x, y, rad2deg(solverS15.phi));
% h=colorbar;
% caxis manual
% caxis([rad2deg(minVal), rad2deg(maxVal)])
% ylabel(h,'\phi','Rotation',0,'FontSize',14)
% xlabel('x/R (-)', 'FontSize',17)
% % ylabel('y/R (-)')
% grid on
% % set(pplot, "edgeColor", "none");
% colormap('default')
% axis square
% title('\gamma = 15', 'FontSize',17)
% nexttile
% pplot = contourf(x, y, rad2deg(solverS30.phi));
% h=colorbar;
% caxis manual
% caxis([rad2deg(minVal), rad2deg(maxVal)])
% ylabel(h,'\phi','Rotation',0,'FontSize',14)
% xlabel('x/R (-)', 'FontSize',17)
% % ylabel('y/R (-)')
% grid on
% % set(pplot, "edgeColor", "none");
% colormap('default')
% title('\gamma = 30', 'FontSize',16)
% axis square
% set(gcf, 'Position', get(0, 'Screensize'));
% set(gcf,'color','w')
% export_fig 'SKEW_phi'
% 
% figure("defaultAxesFontSize", 18)
% minVal = 0; %min([solverS0.a(:); solverS15.a(:); solverS30.a(:)]);
% maxVal = 0.6; %max([solverS0.a(:); solverS15.a(:); solverS30.a(:)]);
% t1 = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
% nexttile
% pplot = contourf(x, y, solverS0.a);
% h=colorbar;
% caxis manual
% caxis([minVal, maxVal])
% ylabel(h,'a','Rotation',0,'FontSize',14)
% xlabel('x/R (-)', 'FontSize',17)
% ylabel('y/R (-)', 'FontSize',17)
% grid on
% % set(pplot, "edgeColor", "none");
% colormap('default')
% title('\gamma = 0', 'FontSize',17)
% % colormap("summer")
% axis square
% nexttile
% pplot = contourf(x, y, solverS15.a);
% h=colorbar;
% caxis manual
% caxis([minVal, maxVal])
% ylabel(h,'a','Rotation',0,'FontSize',14)
% xlabel('x/R (-)', 'FontSize',17)
% % ylabel('y/R (-)')
% grid on
% % set(pplot, "edgeColor", "none");
% colormap('default')
% axis square
% title('\gamma = 15', 'FontSize',17)
% nexttile
% pplot = contourf(x, y, solverS30.a);
% h=colorbar;
% caxis manual
% caxis([minVal, maxVal])
% ylabel(h,'a','Rotation',0,'FontSize',14)
% xlabel('x/R (-)', 'FontSize',17)
% % ylabel('y/R (-)')
% grid on
% % set(pplot, "edgeColor", "none");
% colormap('default')
% title('\gamma = 30', 'FontSize',17)
% axis square
% set(gcf, 'Position', get(0, 'Screensize'));
% set(gcf,'color','w')
% export_fig 'SKEW_a.png'
% 
% figure("defaultAxesFontSize", 18)
% minVal = min([solverS0.aprime(:); solverS15.aprime(:); solverS30.aprime(:)]);
% maxVal = max([solverS0.aprime(:); solverS15.aprime(:); solverS30.aprime(:)]);
% t1 = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
% nexttile
% pplot = contourf(x, y, solverS0.aprime);
% h=colorbar;
% caxis manual
% caxis([minVal, maxVal])
% ylabel(h,"a'",'Rotation',0,'FontSize',14)
% xlabel('x/R (-)', 'FontSize',17)
% ylabel('y/R (-)', 'FontSize',17)
% grid on
% % set(pplot, "edgeColor", "none");
% colormap('default')
% title('\gamma = 0', 'FontSize',17)
% % colormap("summer")
% axis square
% nexttile
% pplot = contourf(x, y, solverS15.aprime);
% h=colorbar;
% caxis manual
% caxis([minVal, maxVal])
% ylabel(h,"a'",'Rotation',0,'FontSize',14)
% xlabel('x/R (-)', 'FontSize',17)
% % ylabel('y/R (-)')
% grid on
% % set(pplot, "edgeColor", "none");
% colormap('default')
% axis square
% title('\gamma = 15', 'FontSize',17)
% nexttile
% pplot = contourf(x, y, solverS30.aprime);
% h=colorbar;
% caxis manual
% caxis([minVal, maxVal])
% ylabel(h,"a'",'Rotation',0,'FontSize',14)
% xlabel('x/R (-)', 'FontSize',17)
% % ylabel('y/R (-)')
% grid on
% % set(pplot, "edgeColor", "none");
% colormap('default')
% title('\gamma = 30', 'FontSize',17)
% axis square
% set(gcf, 'Position', get(0, 'Screensize'));
% set(gcf,'color','w')
% export_fig 'SKEW_at.png'
% 
% figure("defaultAxesFontSize", 18)
% minVal = min([solverS0.CT(:); solverS15.CT(:); solverS30.CT(:)]);
% maxVal = max([solverS0.CT(:); solverS15.CT(:); solverS30.CT(:)]);
% t1 = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
% nexttile
% pplot = contourf(x, y, solverS0.CT);
% h=colorbar;
% caxis manual
% caxis([minVal, maxVal])
% ylabel(h,'C_T','Rotation',0,'FontSize',14)
% xlabel('x/R (-)', 'FontSize',17)
% ylabel('y/R (-)', 'FontSize',17)
% grid on
% % set(pplot, "edgeColor", "none");
% colormap('default')
% title('\gamma=0', 'FontSize',17)
% % colormap("summer")
% axis square
% nexttile
% pplot = contourf(x, y, solverS15.CT);
% h=colorbar;
% caxis manual
% caxis([minVal, maxVal])
% ylabel(h,'C_T','Rotation',0,'FontSize',14)
% xlabel('x/R (-)', 'FontSize',17)
% % ylabel('y/R (-)', 'FontSize',17)
% grid on
% % set(pplot, "edgeColor", "none");
% colormap('default')
% axis square
% title('\gamma = 15', 'FontSize',17)
% nexttile
% pplot = contourf(x, y, solverS30.CT);
% h=colorbar;
% caxis manual
% caxis([minVal, maxVal])
% ylabel(h,'C_T','Rotation',0,'FontSize',14)
% xlabel('x/R (-)', 'FontSize',17)
% % ylabel('y/R (-)')
% grid on
% % set(pplot, "edgeColor", "none");
% colormap('default')
% title('\gamma = 30', 'FontSize',17)
% axis square
% set(gcf, 'Position', get(0, 'Screensize'));
% set(gcf,'color','w')
% export_fig 'SKEW_CT.png'
% 
% figure("defaultAxesFontSize", 18)
% minVal = min([solverS0.CN(:); solverS15.CN(:); solverS30.CN(:)]);
% maxVal = max([solverS0.CN(:); solverS15.CN(:); solverS30.CN(:)]);
% t1 = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
% nexttile
% pplot = contourf(x, y, solverS0.CN);
% h=colorbar;
% caxis manual
% caxis([minVal, maxVal])
% ylabel(h,'C_N','Rotation',0,'FontSize',14)
% xlabel('x/R (-)', 'FontSize',17)
% ylabel('y/R (-)', 'FontSize',17)
% grid on
% % set(pplot, "edgeColor", "none");
% colormap('default')
% title('\gamma = 0', 'FontSize',17)
% % colormap("summer")
% axis square
% nexttile
% pplot = contourf(x, y, solverS15.CN);
% h=colorbar;
% caxis manual
% caxis([minVal, maxVal])
% ylabel(h,'C_N','Rotation',0,'FontSize',14)
% xlabel('x/R (-)', 'FontSize',17)
% % ylabel('y/R (-)')
% grid on
% % set(pplot, "edgeColor", "none");
% colormap('default')
% axis square
% title('\gamma = 15', 'FontSize',17)
% nexttile
% pplot = contourf(x, y, solverS30.CN);
% h=colorbar;
% caxis manual
% caxis([minVal, maxVal])
% ylabel(h,'C_N','Rotation',0)
% xlabel('x/R (-)', 'FontSize',17)
% % ylabel('y/R (-)')
% grid on
% % set(pplot, "edgeColor", "none");
% colormap('default')
% title('\gamma = 30', 'FontSize',17)
% axis square
% set(gcf, 'Position', get(0, 'Screensize'));
% set(gcf,'color','w')
% export_fig 'SKEW_CN.png'
% 
% figure("defaultAxesFontSize", 18)
% minVal = min([solverS0.Cq(:); solverS15.Cq(:); solverS30.Cq(:)]);
% maxVal = max([solverS0.Cq(:); solverS15.Cq(:); solverS30.Cq(:)]);
% t1 = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
% nexttile
% pplot = contourf(x, y, solverS0.Cq);
% h=colorbar;
% caxis manual
% caxis([minVal, maxVal])
% ylabel(h,'C_Q','Rotation',0,'FontSize',14)
% xlabel('x/R (-)', 'FontSize',17)
% ylabel('y/R (-)', 'FontSize',17)
% grid on
% % set(pplot, "edgeColor", "none");
% colormap('default')
% title('\gamma = 0', 'FontSize',17)
% colormap("summer")
% axis square
% nexttile
% pplot = contourf(x, y, solverS15.Cq);
% h=colorbar;
% caxis manual
% caxis([minVal, maxVal])
% ylabel(h,'C_Q','Rotation',0,'FontSize',14)
% xlabel('x/R (-)', 'FontSize',17)
% % ylabel('y/R (-)',, 'FontSize',17)
% grid on
% % set(pplot, "edgeColor", "none");
% colormap('default')
% axis square
% title('\gamma = 15', 'FontSize',17)
% nexttile
% pplot = contourf(x, y, solverS30.Cq);
% h=colorbar;
% caxis manual
% caxis([minVal, maxVal])
% ylabel(h,'C_Q','Rotation',0,'FontSize',14)
% xlabel('x/R (-)', 'FontSize',17)
% % ylabel('y/R (-)')
% grid on
% % set(pplot, "edgeColor", "none");
% colormap('default')
% title('\gamma = 30', 'FontSize',17)
% axis square
% set(gcf, 'Position', get(0, 'Screensize'));
% set(gcf,'color','w')
% export_fig 'SKEW_CQ.png'

% chordS15 = solverS15.computeChordLength(solverS15.rR);
% figure("defaultAxesFontSize", 18)
% tiledlayout(1,2, "TileSpacing", "compact",'Padding', 'none')
% ax1=nexttile;
% hold on
% plot(chordS15, mean(solverS0.Cl, 2), "linewidth", 1.5)
% plot(chordS15, mean(solverS15.Cl, 2), "linewidth", 1.5)
% plot(chordS15, mean(solverS30.Cl, 2), "linewidth", 1.5)
% legend("$\gamma=0$", "$\gamma=15$", "$\gamma=30$", "interpreter", "latex")
% set(gca, 'XDir','reverse')
% xlabel("chord length [m]")
% ylabel("$C_{L}$", "interpreter", "latex")
% grid
% ax2=nexttile;
% hold on 
% plot(polarData(:,3), polarData(:,2), "linewidth", 1.5)
% scatter(LDpoint(3), LDpoint(2), 50, "Marker", 'x', "lineWidth", 2.5)
% legend("Airfoil polar", "Max. L/D")
% grid
% xlabel("$C_{D}$ [-]", "interpreter", "latex")
% ylabel("$C_{L}$", "interpreter", "latex")
% set(gcf,'color','w');
% linkaxes([ax1 ax2], "y")
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