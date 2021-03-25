data = table2array(readtable("polar_DU95W180.xlsx"));

figure("defaultAxesFontSize", 18)
tiledlayout(1, 2, "tileSpacing", "compact")

nexttile
plot(data(:,1), data(:, 2), "linewidth", 1.3)
xlabel("$\alpha$ [$^{\circ}$]", "interpreter", "latex")
ylabel("$C_{l}$ [-]", "interpreter", "latex")
grid
nexttile
plot(data(:,3), data(:, 2), "linewidth", 1.3)
xlabel("$\alpha$ [$^{\circ}$]", "interpreter", "latex")
ylabel("$C_{d}$ [-]", "interpreter", "latex")
grid
set(gcf,'color','w');