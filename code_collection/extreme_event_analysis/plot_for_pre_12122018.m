
% close all previous figures
close all hidden

figsize_pre = [10, 10, 450, 160];

width_thick = 1.5;
width_thin = 0.2;

widths = [  width_thin, ...    % d_omega
            width_thin, ...     % d_omega_QG
            0.8, ...    % d_sigma
            width_thick, ...     % d_J
            width_thin, ...     % d_k2_l2
            width_thin, ...     % d_m2
            width_thin, ...     % d_Adv
            width_thick, ...     % d_dtheta_dp_ma_omega
            width_thin];     % d_rec
			


%%%%%%%%%%%%%% dry decomposition %%%%%%%%%%%%%%%%%%

fig = figure('pos', figsize_pre);

ax1 = subplot(1, 1, 1);
plot(ax1, -70:70, zeros(size(-70:70)), 'LineWidth', 0.8, 'Color', [0.3, 0.3, 0.3], 'Linestyle', '-.')
hold on;

p0_d_S = plot(ax1, plot_x_S + bias, f(nanmean(d_omega             (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2)) ./ omega_h_mean_d(indS));
p0_d_N = plot(ax1, plot_x_N - bias, f(nanmean(d_omega             (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2)) ./ omega_h_mean_d(indN));

p1_d_S = plot(ax1, plot_x_S + bias, f(nanmean(d_omega_QG          (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p1_d_N = plot(ax1, plot_x_N - bias, f(nanmean(d_omega_QG          (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));

p2_d_S = plot(ax1, plot_x_S + bias, f(nanmean(d_sigma             (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p2_d_N = plot(ax1, plot_x_N - bias, f(nanmean(d_sigma             (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));

p3_d_S = plot(ax1, plot_x_S + bias, f(nanmean(d_J                 (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p3_d_N = plot(ax1, plot_x_N - bias, f(nanmean(d_J                 (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));

p4_d_S = plot(ax1, plot_x_S + bias, f(nanmean(d_k2_l2             (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p4_d_N = plot(ax1, plot_x_N - bias, f(nanmean(d_k2_l2             (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));

p5_d_S = plot(ax1, plot_x_S + bias, f(nanmean(d_m2                (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p5_d_N = plot(ax1, plot_x_N - bias, f(nanmean(d_m2                (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));

p6_d_S = plot(ax1, plot_x_S + bias, f(nanmean(d_Adv               (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p6_d_N = plot(ax1, plot_x_N - bias, f(nanmean(d_Adv               (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));

p7_d_S = plot(ax1, plot_x_S + bias, f(nanmean(d_dtheta_dp_ma_omega(indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p7_d_N = plot(ax1, plot_x_N - bias, f(nanmean(d_dtheta_dp_ma_omega(indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));

p8_d_S = plot(ax1, plot_x_S + bias, f(nanmean(d_rec               (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p8_d_N = plot(ax1, plot_x_N - bias, f(nanmean(d_rec               (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));

set([p0_d_S, p0_d_N], 'LineWidth', widths(1), 'Color', 'black', 'Linestyle', '--');
set([p1_d_S, p1_d_N], 'LineWidth', widths(2), 'Color', 'black');
set([p2_d_S, p2_d_N], 'LineWidth', widths(3), 'Color', colors(1, :), 'Marker', 's', 'MarkerSize', 2);
set([p3_d_S, p3_d_N], 'LineWidth', widths(4), 'Color', colors(7, :));
set([p4_d_S, p4_d_N], 'LineWidth', widths(5), 'Color', colors(3, :));
set([p5_d_S, p5_d_N], 'LineWidth', widths(6), 'Color', colors(5, :));
set([p6_d_S, p6_d_N], 'LineWidth', widths(7), 'Color', colors(6, :));
set([p7_d_S, p7_d_N], 'LineWidth', widths(8), 'Color', colors(7, :), 'Linestyle', '--');
set([p8_d_S, p8_d_N], 'LineWidth', widths(9), 'Color', 'black', 'Linestyle', ':');

lgd = legend( ...
[p7_d_S, p3_d_S, p0_d_S, p1_d_S, p8_d_S, p5_d_S, p6_d_S, p4_d_S, p2_d_S], ...
{'$\frac{p}{\kappa}\sigma^*$', ...
 '${J}$'...
 '${\omega}$' , ...
 '${\omega}_{QG}$' , ...
 '${\omega_{QG-rec}}$', ...
 '$m$', ...
 '${Adv}$' , ...
 '$k,k_J$', ...
 '${\sigma}$'}, ...
 'location', 'best', ...
 'interpreter', 'latex');
lgd.Position = [0.845, -0.05, 0.08, 1.03];
set(lgd,'Box','off')

axis([-73 + bias 73 - bias -0.20 0.30]);
yticks([-0.06:0.03:0.06]*delta_T);
yticklabels({'-6', '-3', '0', '3', '6'})
xticks([-70.0 + bias, -50 + bias, lat_1 + bias, ...
        lat_2 - bias, 50 - bias, 70.0 - bias]);
xticklabels({'70S', '50S', [num2str(-lat_1), 'S'], [num2str(lat_2), 'N'], '50N', '70N'})
ylabel('(\%/K)', 'Interpreter', 'latex');
ylab = get(ax1,'ylabel');
set(ylab,'Units','normalized');
ylab.Position(1) = ylab.Position(1) + 0.1;
ylab.Position(2) = ylab.Position(2) - 1.5;

set(gca, 'Position', [0.08, 0.12, 0.72, 0.76]);
set(gca, 'TickLabelInterpreter', 'latex')
break_axis('handle', gca, 'position', 0.0, 'axis', 'x', 'length', 0.20)
hold off;
box off;

figname = [plot_path_ensemble, 'omega_decomposition_dry_pre_12122018_', plot_level_name, '_6'];
print(fig, figname, '-depsc', '-r0', '-painters')

%saveas(gca, [plot_path_ensemble, 'omega_decomposition_moist_pre_12122018_', plot_level_name], format_str);
clf;


if strfind(pwd_str, 'GFDL')
    f = @(x) x;
else
    f = @one_two_one; % define the smoothing function
end

%%%%%%%%%%%%%% moist decomposition %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% 1st plot %%%%%%%%%%%%%%%%%%

widths = [  width_thick, ...    % d_sigma_m
            width_thick, ...     % d_k2_l2_m
            width_thick, ...    % d_m2_m
            width_thick, ...     % d_Adv_m
            width_thick, ...     % d_rec_m_full
            width_thick, ...     % d_omega_QG
            width_thick];        % d_J_res
        

fig = figure('pos', figsize_pre);
plot(-70:70, zeros(size(-70:70)), 'LineWidth', 0.8, 'Color', [0.3, 0.3, 0.3], 'Linestyle', '-.')
hold on;
p3_1_m_S = plot(plot_x_S + bias, f(nanmean(d_sigma_m    (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean(indS));
p3_1_m_N = plot(plot_x_N - bias, f(nanmean(d_sigma_m    (indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean(indN));
set([p3_1_m_S, p3_1_m_N], 'LineWidth', widths(1), 'Color', colors(1, :));

axis([-73 + bias 73 - bias -0.08 0.15]);
yticks([-0.06:0.015:0.06]*delta_T);
yticklabels({'-6', '-4.5', '-3', '-1.5', '0', '1.5', '3', '4.5', '6'})
xticks([-70.0 + bias, -50 + bias, lat_1 + bias, ...
        lat_2 - bias, 50 - bias, 70.0 - bias]);
xticklabels({'70S', '50S', [num2str(-lat_1), 'S'], [num2str(lat_2), 'N'], '50N', '70N'})
ylabel('(\%/K)', 'Interpreter', 'latex');
ax = gca;
ylab = get(ax,'ylabel');
set(ylab,'Units','normalized');
ylab.Position(1) = ylab.Position(1) + 0.1;
ylab.Position(2) = ylab.Position(2) - 1.5;

set(gca, 'Position', [0.08, 0.12, 0.72, 0.76]);
set(gca, 'TickLabelInterpreter', 'latex')
break_axis('handle', gca, 'position', 0.0, 'axis', 'x', 'length', 0.20)
hold off;
box off;

%saveas(gca, [plot_path_ensemble, 'omega_decomposition_moist_pre_12122018_', plot_level_name, '_1'], format_str);



%%%%%%%%%%%%%% 2nd plot %%%%%%%%%%%%%%%%%%
hold on;
p3_2_m_S = plot(plot_x_S + bias, f(nanmean(d_k2_l2_m    (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean(indS));
p3_2_m_N = plot(plot_x_N - bias, f(nanmean(d_k2_l2_m    (indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean(indN));
set([p3_2_m_S, p3_2_m_N], 'LineWidth', widths(2), 'Color', colors(3, :));

%saveas(gca, [plot_path_ensemble, 'omega_decomposition_moist_pre_12122018_', plot_level_name, '_2'], format_str);

%%%%%%%%%%%%%% 3rd plot %%%%%%%%%%%%%%%%%%
hold on;
p3_3_m_S = plot(plot_x_S + bias, f(nanmean(d_m2_m       (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean(indS));
p3_3_m_N = plot(plot_x_N - bias, f(nanmean(d_m2_m       (indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean(indN));
set([p3_3_m_S, p3_3_m_N], 'LineWidth', widths(3), 'Color', colors(5, :));

%saveas(gca, [plot_path_ensemble, 'omega_decomposition_moist_pre_12122018_', plot_level_name, '_3'], format_str);

%%%%%%%%%%%%%% 4th plot %%%%%%%%%%%%%%%%%%
hold on;
p2_m_S   = plot(plot_x_S + bias, f(nanmean(d_Adv_m      (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean(indS));
p2_m_N   = plot(plot_x_N - bias, f(nanmean(d_Adv_m      (indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean(indN));
set([p2_m_S, p2_m_N], 'LineWidth', widths(4), 'Color', colors(6, :));

%saveas(gca, [plot_path_ensemble, 'omega_decomposition_moist_pre_12122018_', plot_level_name, '_4'], format_str);

%%%%%%%%%%%%%% 5th plot %%%%%%%%%%%%%%%%%%
hold on;
p6_m_S   = plot(plot_x_S + bias, f(nanmean(d_rec_m_full(indS, :, plot_ind).* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean(indS));
p6_m_N   = plot(plot_x_N - bias, f(nanmean(d_rec_m_full(indN, :, plot_ind).* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean(indN));
set([p6_m_S, p6_m_N], 'LineWidth', widths(5), 'Color', 'black', 'Linestyle', ':');

%saveas(gca, [plot_path_ensemble, 'omega_decomposition_moist_pre_12122018_', plot_level_name, '_5'], format_str);

%%%%%%%%%%%%%% 6th plot %%%%%%%%%%%%%%%%%%
hold on;
p1_m_S   = plot(plot_x_S + bias, f(nanmean(d_omega_QG   (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean(indS));
p1_m_N   = plot(plot_x_N - bias, f(nanmean(d_omega_QG   (indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean(indN));
set([p1_m_S, p1_m_N], 'LineWidth', widths(6), 'Color', 'black');

%%%%%%%%%%%%%% 7th plot %%%%%%%%%%%%%%%%%%
hold on;
%p7_m_S  = plot(plot_x_S + bias, nanmean((d_excessive_1(indS, :, plot_ind) + ...
%                                         d_excessive_2(indS, :, plot_ind) + ...
%                                         d_epsilon    (indS, :, plot_ind)).* NaN_matrix_m(indS, :), 2) ./ omega_QG_h_mean(indS));
%p7_m_N  = plot(plot_x_N - bias, nanmean((d_excessive_1(indN, :, plot_ind) + ...
%                                         d_excessive_2(indN, :, plot_ind) + ...
%                                         d_epsilon    (indN, :, plot_ind)).* NaN_matrix_m(indN, :), 2) ./ omega_QG_h_mean(indN));
p7_m_S  = plot(plot_x_S + bias, f(nanmean((d_J_res    (indS, :, plot_ind)).* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean(indS));
p7_m_N  = plot(plot_x_N - bias, f(nanmean((d_J_res    (indN, :, plot_ind)).* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean(indN));
set([p7_m_S, p7_m_N], 'LineWidth', widths(7), 'Color', colors(4, :));


lgd = legend( ...
[p1_m_N, p6_m_N, p3_1_m_N, p3_3_m_N, p2_m_N, p3_2_m_N, p7_m_N], ...
{...%'$\omega$' , ...
'$\omega_{QG}$' , ...
'$\omega_{QG-rec}$', ...
'$\sigma_m$', ...
'$m$', ...
'$Adv$', ...
'${k, k_J}$', ...
'$J_{res}$'}, ...
'interpreter', 'latex');
lgd.Position = [0.845, 0.08, 0.08, 0.80];
set(lgd,'Box','off')

%saveas(gca, [plot_path_ensemble, 'omega_decomposition_moist_pre_12122018_', plot_level_name, '_6'], format_str);
figname = [plot_path_ensemble, 'omega_decomposition_moist_pre_12122018_', plot_level_name, '_6'];
print(fig, figname, '-depsc', '-r0', '-painters')


