
if ~isempty(strfind(filename, 'CESM_6hourly'))
    EPS = true;
else
    EPS = false;
end

% read in global temperature changes
if strfind(filename, 'CESM')
    load('/disk7/ziweili/CESM_LENS/exp/global_temperature/global_avg_T.mat');
elseif strfind(filename, 'GFDL')
    load('/disk7/ziweili/test1_GFDL/exp/global_temperature/global_avg_T.mat');
end
delta_T = mean(T_avg_r) - mean(T_avg_h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dry decomposition, zonal average
figsize_decomp = [10, 10, 500, 250];
left_margin     = 0.06;
right_margin    = 0.22;
down_margin     = 0.13;
up_margin       = 0.08;
distance        = 0.10;
width           = 1 - left_margin - right_margin;
height          = (1 - down_margin - up_margin - distance) / 2;

position_1 = [left_margin, down_margin + height + distance, width, height];
position_2 = [left_margin, down_margin                    , width, height];

omega_QG_h_mean = nanmean(omega_QG_h_rec .* NaN_matrix_m, 2);

fig = figure('pos', figsize_decomp);

lat_1 = -30;
lat_2 = 30;
indN = plot_x > 30;
indS = plot_x < -30;
plot_x_N = plot_x(indN);
plot_x_S = plot_x(indS);
colors = get(gca,'colororder');
colors(7, 1) = colors(7, 1) + 0.2;
colors(7, 3) = colors(7, 3) - 0.05;
bias = 25;

ax1 = subplot(2, 1, 1);
plot(-90:90, zeros(size(-90:90)), 'LineWidth', 0.5, 'Color', [0.3, 0.3, 0.3], 'Linestyle', '-.')
hold on;

if strfind(filename, 'GFDL')
    f = @(x) x;
else
    f = @one_two_one; % define the smoothing function
end

p0_d = plot(ax1, plot_x_S + bias, f(nanmean(d_omega             (indS, :) .* NaN_matrix_d(indS, :), 2)) ./ omega_h_mean_d(indS));
p1_d = plot(ax1, plot_x_S + bias, f(nanmean(d_omega_QG          (indS, :) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p2_d = plot(ax1, plot_x_S + bias, f(nanmean(d_sigma             (indS, :) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p3_d = plot(ax1, plot_x_S + bias, f(nanmean(d_J                 (indS, :) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p4_d = plot(ax1, plot_x_S + bias, f(nanmean(d_k2_l2             (indS, :) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p5_d = plot(ax1, plot_x_S + bias, f(nanmean(d_m2                (indS, :) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p6_d = plot(ax1, plot_x_S + bias, f(nanmean(d_Adv               (indS, :) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p7_d = plot(ax1, plot_x_S + bias, f(nanmean(d_dtheta_dp_ma_omega(indS, :) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p8_d = plot(ax1, plot_x_S + bias, f(nanmean(d_rec               (indS, :) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));

set(p0_d, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', '--');
set(p1_d, 'LineWidth', 1.5, 'Color', 'black');
set(p2_d, 'LineWidth', 0.8, 'Color', colors(1, :), 'Marker', 's', 'MarkerSize', 2);
set(p3_d, 'LineWidth', 1.5, 'Color', colors(7, :));
set(p4_d, 'LineWidth', 1.5, 'Color', colors(3, :));
set(p5_d, 'LineWidth', 1.5, 'Color', colors(5, :));
set(p6_d, 'LineWidth', 1.5, 'Color', colors(6, :));
set(p7_d, 'LineWidth', 1.5, 'Color', colors(7, :), 'Linestyle', '-.');
set(p8_d, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', ':');

p0_d = plot(ax1, plot_x_N - bias, f(nanmean(d_omega             (indN, :) .* NaN_matrix_d(indN, :), 2)) ./ omega_h_mean_d(indN));
p1_d = plot(ax1, plot_x_N - bias, f(nanmean(d_omega_QG          (indN, :) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));
p2_d = plot(ax1, plot_x_N - bias, f(nanmean(d_sigma             (indN, :) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));
p3_d = plot(ax1, plot_x_N - bias, f(nanmean(d_J                 (indN, :) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));
p4_d = plot(ax1, plot_x_N - bias, f(nanmean(d_k2_l2             (indN, :) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));
p5_d = plot(ax1, plot_x_N - bias, f(nanmean(d_m2                (indN, :) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));
p6_d = plot(ax1, plot_x_N - bias, f(nanmean(d_Adv               (indN, :) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));
p7_d = plot(ax1, plot_x_N - bias, f(nanmean(d_dtheta_dp_ma_omega(indN, :) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));
p8_d = plot(ax1, plot_x_N - bias, f(nanmean(d_rec               (indN, :) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));

set(p0_d, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', '--');
set(p1_d, 'LineWidth', 1.5, 'Color', 'black');
set(p2_d, 'LineWidth', 0.8, 'Color', colors(1, :), 'Marker', 's', 'MarkerSize', 2);
set(p3_d, 'LineWidth', 1.5, 'Color', colors(7, :));
set(p4_d, 'LineWidth', 1.5, 'Color', colors(3, :));
set(p5_d, 'LineWidth', 1.5, 'Color', colors(5, :));
set(p6_d, 'LineWidth', 1.5, 'Color', colors(6, :));
set(p7_d, 'LineWidth', 1.5, 'Color', colors(7, :), 'Linestyle', '-.');
set(p8_d, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', ':');

title(['(a) Dry Decomposition'], 'FontWeight', 'normal', 'interpreter', 'latex')
ylabel('\%K$^{-1}$', 'Interpreter', 'latex', 'FontSize', 9);
ylab1 = get(ax1, 'ylabel');
set(ylab1, 'rotation', 0)
set(ylab1,'Units','normalized');
ylab1.Position(1) = ylab1.Position(1) + 0.06;
ylab1.Position(2) = ylab1.Position(2) / height * (2.1 * height);


if strfind(filename, 'GFDL')
    axis([-73 + bias 73 - bias -0.40 0.65]);
    yticks([-0.15:0.05:0.15]*delta_T);
    yticklabels({'-15', '-10', '-5', '0', '5', '10', '15'})
else
    axis([-73 + bias 73 - bias -0.20 0.30]);
    yticks([-0.06:0.03:0.06]*delta_T);
    yticklabels({'-6', '-3', '0', '3', '6'})
end
xticks([-70.0 + bias, -50 + bias, lat_1 + bias, ...
        lat_2 - bias,  50 - bias, 70.0 - bias])
xticklabels([])
%xticklabels({'-60', num2str(lat_1), num2str(lat_2), '60'})

set([ax1], 'Position', position_1);
set(ax1, 'TickLabelInterpreter', 'latex')
break_axis('handle', ax1, 'position', 0.0, 'axis', 'x', 'length', 0.20)
set(ax1,'TickDir','out');
%ax1.XTicklabels = [];
box off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% moist decomposition, zonal average
ax2 = subplot(2, 1, 2);

d_rec_m_full_paper    = d_sigma_m + d_k2_l2_m + d_m2_m + d_Adv_m + d_J_res;

plot(-90:90, zeros(size(-90:90)), 'LineWidth', 0.5, 'Color', [0.3, 0.3, 0.3], 'Linestyle', '-.')
hold on;

p1_m   = plot(plot_x_S + bias, f(nanmean(d_omega_QG (indS, :) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS));
p2_m   = plot(plot_x_S + bias, f(nanmean(d_Adv_m    (indS, :) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS));
p3_1_m = plot(plot_x_S + bias, f(nanmean(d_sigma_m  (indS, :) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS));
p3_2_m = plot(plot_x_S + bias, f(nanmean(d_k2_l2_m  (indS, :) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS));
p3_3_m = plot(plot_x_S + bias, f(nanmean(d_m2_m     (indS, :) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS));
p4_1_m = plot(plot_x_S + bias, f(nanmean(d_J_res    (indS, :) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS));
%p4_2_m = plot(plot_x_S + bias, f(nanmean(d_excessive_2(indS, :) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS), 'o');
%p5_m   = plot(plot_x_S + bias, f(nanmean(d_epsilon    (indS, :) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS));
p6_m   = plot(plot_x_S + bias, f(nanmean(d_rec_m_full_paper(indS, :) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS));

set(p1_m, 'LineWidth', 1.5, 'Color', 'black');
set(p2_m, 'LineWidth', 1.5, 'Color', colors(6, :));
set(p3_1_m, 'LineWidth', 0.8, 'Color', colors(1, :), 'Marker', 's', 'MarkerSize', 2);
set(p3_2_m, 'LineWidth', 1.5, 'Color', colors(3, :));
set(p3_3_m, 'LineWidth', 1.5, 'Color', colors(5, :));
set(p4_1_m, 'LineWidth', 1.5, 'Color', colors(7, :), 'Linestyle', '--', 'MarkerSize', 2);
%set(p4_2_m, 'LineWidth', 0.8, 'Color', colors(2, :), 'Linestyle', '--', 'MarkerSize', 2);
%set(p5_m, 'LineWidth', 1.5, 'Color', colors(7, :));
set(p6_m, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', ':');

p1_m   = plot(plot_x_N - bias, f(nanmean(d_omega_QG(indN, :) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN));
p2_m   = plot(plot_x_N - bias, f(nanmean(d_Adv_m   (indN, :) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN));
p3_1_m = plot(plot_x_N - bias, f(nanmean(d_sigma_m (indN, :) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN));
p3_2_m = plot(plot_x_N - bias, f(nanmean(d_k2_l2_m (indN, :) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN));
p3_3_m = plot(plot_x_N - bias, f(nanmean(d_m2_m    (indN, :) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN));
p4_1_m = plot(plot_x_N - bias, f(nanmean(d_J_res   (indN, :) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN));
%p4_2_m = plot(plot_x_N - bias, f(nanmean(d_excessive_2(indN, :) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN), 'o');
%p5_m   = plot(plot_x_N - bias, f(nanmean(d_epsilon    (indN, :) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN));
p6_m   = plot(plot_x_N - bias, f(nanmean(d_rec_m_full_paper(indN, :) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN));

set(p1_m, 'LineWidth', 1.5, 'Color', 'black');
set(p2_m, 'LineWidth', 1.5, 'Color', colors(6, :));
set(p3_1_m, 'LineWidth', 0.8, 'Color', colors(1, :), 'Marker', 's', 'MarkerSize', 2);
set(p3_2_m, 'LineWidth', 1.5, 'Color', colors(3, :));
set(p3_3_m, 'LineWidth', 1.5, 'Color', colors(5, :));
set(p4_1_m, 'LineWidth', 1.5, 'Color', colors(7, :), 'Linestyle', '--');
%set(p4_2_m, 'LineWidth', 0.8, 'Color', colors(2, :), 'Linestyle', '--', 'MarkerSize', 2);
%set(p5_m, 'LineWidth', 1.5, 'Color', colors(7, :));
set(p6_m, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', ':');



if strfind(filename, 'GFDL')
    axis([-73 + bias 73 - bias -0.25 0.38]);
    yticks([-0.075:0.025:0.075]*delta_T);
    yticklabels({'-7.5', '-5', '-2.5', '0', '2.5', '5', '7.5'})
else
    axis([-73 + bias 73 - bias -0.09 0.15]);
    yticks([-0.06:0.015:0.06]*delta_T);
    yticklabels({'-6', '-4.5', '-3', '-1.5', '0', '1.5', '3', '4.5', '6'})
end
xticks([-70.0 + bias, -50 + bias, lat_1 + bias, ...
        lat_2 - bias, 50 - bias, 70.0 - bias]);
xticklabels({'$70^\circ$S', '$50^\circ$S', ['$', num2str(abs(lat_1)), '^\circ$S'], ...
             ['$', num2str(lat_2), '^\circ$N'], '$50^\circ$N', '$70^\circ$N'})
set(ax2, 'TickLabelInterpreter', 'latex')

    
lgd = legend( ...
[p7_d, p3_d, p0_d, p1_m, p6_m, p3_3_m, p2_m, p3_2_m, p3_1_m, p4_1_m], ...
{'$\frac{p}{\kappa}\sigma^*$', ...
 '${J}$'...
 ...%'${c_p{\frac{T}{\theta}\frac{d\theta}{dp}\big|_{\theta*}\omega}}$' , ...
 ...'(a) $\frac{p}{\kappa}\sigma^*$; (b) $\overline{Res}$', ...
 '$\overline{\omega}$' , ...
 '$\overline{\omega}_{QG}$' , ...
 '${\omega_{QG-rec}}$', ...
 '$m$', ...
 '$\overline{Adv}$' , ...
 '$k,k_J$', ...
 '(a) $\hat{\sigma}$, (b) $\hat{\sigma}_m$', ...
 '$\overline{J}_{Res}$'}, ...
 ...'$J''_{\frac{d\theta}{dp}\big|_{\theta*}}$' , ...
 ...'$J''_{\omega-\omega_{QG}}$'},  ...
 'location', 'best', ...
 'interpreter', 'latex');
lgd.Position = [left_margin + width + 0.07, 0.10, 0.08, 0.80];
set(lgd,'Box','off')

title(['(b) Moist Decomposition'], 'FontWeight', 'normal', 'Interpreter', 'latex');
xlabel('Latitude', 'Interpreter', 'latex');
ylabel('\%K$^{-1}$', 'Interpreter', 'latex', 'FontSize', 9);
ylab2 = get(ax2, 'ylabel');
set(ylab2, 'rotation', 0)
set(ylab2,'Units','normalized');
ylab2.Position(1) = ylab2.Position(1) + 0.05;
ylab2.Position(2) = ylab2.Position(2) / height * (2.1 * height);

set([ax2], 'Position', position_2);
set(ax2, 'TickLabelInterpreter', 'latex')
break_axis('handle', ax2, 'position', 0.0, 'axis', 'x', 'length', 0.20)
set(ax2,'TickDir','out');
hold off;
box off;


figname = [plot_path_ensemble, 'omega_decomposition_moist_paper_', plot_level_name];

if EPS
    print(fig, figname, '-depsc', '-r0', '-painters')
else
    saveas(fig, figname, 'png')
end


