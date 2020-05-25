% plot zonal mean omega vs. omega_QG

EPS = false;

close all hidden
figsize = [10, 10, 450, 160];

fig = figure('pos', figsize);
ax1 = subplot(1, 1, 1);

colors = get(gca,'colororder');
hold on;

p1_S = plot(ax1, plot_x_S + bias, - omega_h_mean_2(indS));
p1_N = plot(ax1, plot_x_N - bias, - omega_h_mean_2(indN));

p2_S = plot(ax1, plot_x_S + bias, - omega_QG_h_mean_2(indS));
p2_N = plot(ax1, plot_x_N - bias, - omega_QG_h_mean_2(indN));

p3_S = plot(ax1, plot_x_S + bias, - omega_r_mean_2(indS));
p3_N = plot(ax1, plot_x_N - bias, - omega_r_mean_2(indN));

p4_S = plot(ax1, plot_x_S + bias, - omega_QG_r_mean_2(indS));
p4_N = plot(ax1, plot_x_N - bias, - omega_QG_r_mean_2(indN));

set([p1_S, p2_S, p3_S, p4_S, ...
     p1_N, p2_N, p3_N, p4_N], 'LineWidth', 1)
set([p1_S, p1_N], 'LineWidth', 1, 'Color', colors(1, :), 'Linestyle', '-');
set([p2_S, p2_N], 'LineWidth', 1, 'Color', colors(1, :), 'Linestyle', '-.');
set([p3_S, p3_N], 'LineWidth', 1, 'Color', colors(2, :), 'Linestyle', '-');
set([p4_S, p4_N], 'LineWidth', 1, 'Color', colors(2, :), 'Linestyle', '-.');

lgd = legend( ...
[p1_S, p2_S, p3_S, p4_S], {...
'Historical $\overline{\omega}$', 'Historical $\overline{\omega}_{QG}$', ...
'RCP8.5 $\overline{\omega}$'    , 'RCP8.5 $\overline{\omega}_{QG}$'}, ...
'location', 'southeast', ...
'interpreter', 'latex');
lgd.Position = [0.86, 0.25, 0.08, 0.5];
set(lgd,'Box','off')
set(ax1,'TickDir','out');


if ~GFDL
    axis([-73 + bias 73 - bias 0.4 1.2]);
    yticks([0:0.2:1.2])
else
    axis([-73 + bias 73 - bias 0.2 0.9]);
    yticks([0:0.1:1.0])
end
xticks([-70.0 + bias, -50 + bias, lat_1 + bias, ...
        lat_2 - bias, 50 - bias, 70.0 - bias]);
xticklabels({'70$^\circ$S', '50$^\circ$S', [num2str(-lat_1), '$^\circ$S'], [num2str(lat_2), '$^\circ$N'], '50$^\circ$N', '70$^\circ$N'})
xlabel('Latitude', 'interpreter', 'latex');

ylabel('$-\omega$ (Pa s$^{-1}$)', 'interpreter', 'latex')

set(gca, 'Position', [0.08, 0.12, 0.72, 0.76]);
set(gca, 'ticklabelinterpreter', 'latex')
break_axis('handle', gca, 'position', 0.0, 'axis', 'x', 'length', 0.20)
hold off;

if EPS
    print(fig, figname, '-depsc', '-r0', '-painters')
else
    %saveas(gca, figname, 'png');
    print(fig, figname, '-dpng', '-r300'); 
end




