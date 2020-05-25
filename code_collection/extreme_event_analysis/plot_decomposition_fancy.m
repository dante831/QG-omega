
if ~isempty(strfind(pwd_str, 'GFDL')) || ...
   ~isempty(strfind(pwd_str, 'daily'))
    EPS = false;
else
    EPS = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dry decomposition, 2D map

figsize_multipanel = [10, 10, 600, 500];
main_fig = fig_class;
main_fig.handle = figure('pos', figsize_multipanel);
main_fig.row = 4;
main_fig.col = 2;
main_fig.axes = {};

% parameters of the plot
left_margin     = 0.05;
right_margin    = 0.09;
down_margin     = 0.09;
up_margin       = 0.05;
distance        = 30;
margins = [left_margin, right_margin, down_margin, up_margin];
lat_mask = 3;
ZERO = true;
FLIP = true;
LIGHT_COLOR = false;
ddd = 1 - down_margin - up_margin;
cb_pos = [1 - right_margin + 0.02, ddd/4 + down_margin, 0.02, ddd/4*2];
cb_string = '\%K$^{-1}$';
cb_label_pos_bias = [-0.18, -0.465];
cb_ticks = -0.12:0.06:0.12;
cb_ticklabels = {'-12', '-6', '0', '6', '12'};
climit = [-0.15, 0.15];


temp = d_omega(:, :, plot_ind) ./ Omega_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(a) $\overline{\omega}$';
main_fig = subplot_2D_map(main_fig, @jet, margins, 1, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_omega_QG(:, :, plot_ind) ./ Omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(b) $\overline{\omega}_{QG}$';
main_fig = subplot_2D_map(main_fig, @jet, margins, 2, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_rec(:, :, plot_ind) ./ Omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(c) All contributions';
main_fig = subplot_2D_map(main_fig, @jet, margins, 3, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_J(:, :, plot_ind) ./ Omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(d) Diabatic heating $\overline{J}$';
main_fig = subplot_2D_map(main_fig, @jet, margins, 4, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_sigma(:, :, plot_ind) ./ Omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(e) Static stability $\hat{\sigma}$';
main_fig = subplot_2D_map(main_fig, @jet, margins, 5, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_k2_l2(:, :, plot_ind) ./ Omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(f) Horizontal wavenumbers ${k},{k}_J$';
main_fig = subplot_2D_map(main_fig, @jet, margins, 6, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_m2(:, :, plot_ind) ./ Omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(g) Vertical wavenumber ${m}$';
main_fig = subplot_2D_map(main_fig, @jet, margins, 7, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_Adv(:, :, plot_ind) ./ Omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(h) Advective forcing $\overline{Adv}$';
main_fig = subplot_2D_map(main_fig, @jet, margins, 8, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

positions = set_axes(main_fig, figsize_multipanel, margins, distance, distance);
figname = [plot_path_ensemble, 'dry_decomposition_', plot_level_name];

if EPS
    print(main_fig.handle, figname, '-depsc', '-r0', '-painters')
else
    saveas(main_fig.handle, figname, 'png')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dry decomposition, zonal average
figsize_decomp = [10, 10, 500, 250];
left_margin     = 0.08;
right_margin    = 0.20;
down_margin     = 0.12;
up_margin       = 0.10;
distance        = 0.10;
width           = 1 - left_margin - right_margin;
height          = (1 - down_margin - up_margin - distance) / 2;

position_1 = [left_margin, down_margin + height + distance, width, height];
position_2 = [left_margin, down_margin                    , width, height];

omega_QG_h_mean_m = nanmean(omega_QG_h_rec(:, :, plot_ind) .* NaN_matrix_m, 2);

fig2 = figure('pos', figsize_decomp);

lat_1 = -30;
lat_2 = 30;
indN = plot_x > 30;
indS = plot_x < -30;
plot_x_N = plot_x(indN);
plot_x_S = plot_x(indS);
colors = get(gca,'colororder');
bias = 25;

ax1 = axes(fig2);
plot(-90:90, zeros(size(-90:90)), 'LineWidth', 0.8, 'Color', [0.3, 0.3, 0.3], 'Linestyle', '-.')
hold on;

%p1 = plot(plot_x_S + bias, nanmean(d_omega_QG   (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2) ./ omega_QG_h_mean(indS));

p0_d = plot(ax1, plot_x_S + bias, nanmean(d_omega             (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2) ./ omega_h_mean_d(indS));
p1_d = plot(ax1, plot_x_S + bias, nanmean(d_omega_QG          (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2) ./ omega_QG_h_mean_d(indS));
p2_d = plot(ax1, plot_x_S + bias, nanmean(d_sigma             (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2) ./ omega_QG_h_mean_d(indS));
p3_d = plot(ax1, plot_x_S + bias, nanmean(d_J                 (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2) ./ omega_QG_h_mean_d(indS));
p4_d = plot(ax1, plot_x_S + bias, nanmean(d_k2_l2             (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2) ./ omega_QG_h_mean_d(indS));
p5_d = plot(ax1, plot_x_S + bias, nanmean(d_m2                (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2) ./ omega_QG_h_mean_d(indS));
p6_d = plot(ax1, plot_x_S + bias, nanmean(d_Adv               (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2) ./ omega_QG_h_mean_d(indS));
p7_d = plot(ax1, plot_x_S + bias, nanmean(d_dtheta_dp_ma_omega(indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2) ./ omega_QG_h_mean_d(indS));
p8_d = plot(ax1, plot_x_S + bias, nanmean(d_rec               (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2) ./ omega_QG_h_mean_d(indS));

set(p0_d, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', '-.');
set(p1_d, 'LineWidth', 1.5, 'Color', 'black');
set(p2_d, 'LineWidth', 1.5, 'Color', colors(1, :));
set(p3_d, 'LineWidth', 1.5, 'Color', colors(2, :));
set(p4_d, 'LineWidth', 1.5, 'Color', colors(3, :));
set(p5_d, 'LineWidth', 1.5, 'Color', colors(5, :));
set(p6_d, 'LineWidth', 1.5, 'Color', colors(6, :));
set(p7_d, 'LineWidth', 1.5, 'Color', colors(7, :), 'Linestyle', '--');
set(p8_d, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', ':');

p0_d = plot(ax1, plot_x_N - bias, nanmean(d_omega             (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2) ./ omega_h_mean_d(indN));
p1_d = plot(ax1, plot_x_N - bias, nanmean(d_omega_QG          (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2) ./ omega_QG_h_mean_d(indN));
p2_d = plot(ax1, plot_x_N - bias, nanmean(d_sigma             (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2) ./ omega_QG_h_mean_d(indN));
p3_d = plot(ax1, plot_x_N - bias, nanmean(d_J                 (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2) ./ omega_QG_h_mean_d(indN));
p4_d = plot(ax1, plot_x_N - bias, nanmean(d_k2_l2             (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2) ./ omega_QG_h_mean_d(indN));
p5_d = plot(ax1, plot_x_N - bias, nanmean(d_m2                (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2) ./ omega_QG_h_mean_d(indN));
p6_d = plot(ax1, plot_x_N - bias, nanmean(d_Adv               (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2) ./ omega_QG_h_mean_d(indN));
p7_d = plot(ax1, plot_x_N - bias, nanmean(d_dtheta_dp_ma_omega(indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2) ./ omega_QG_h_mean_d(indN));
p8_d = plot(ax1, plot_x_N - bias, nanmean(d_rec               (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2) ./ omega_QG_h_mean_d(indN));

set(p0_d, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', '-.');
set(p1_d, 'LineWidth', 1.5, 'Color', 'black');
set(p2_d, 'LineWidth', 1.5, 'Color', colors(1, :));
set(p3_d, 'LineWidth', 1.5, 'Color', colors(2, :));
set(p4_d, 'LineWidth', 1.5, 'Color', colors(3, :));
set(p5_d, 'LineWidth', 1.5, 'Color', colors(5, :));
set(p6_d, 'LineWidth', 1.5, 'Color', colors(6, :));
set(p7_d, 'LineWidth', 1.5, 'Color', colors(7, :), 'Linestyle', '--');
set(p8_d, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', ':');

title(['Dry Decomposition'], 'FontWeight', 'normal')
%xlabel('Latitude', 'interpreter', 'latex');
%ylabel('$\Delta\omega/\omega$', 'interpreter', 'latex');

pwd_str = pwd;
if strfind(pwd_str, 'GFDL')
    axis([-63 + bias 63 - bias -0.20 0.50]);
    yticks([-0.2, 0, 0.2, 0.4]);
    yticklabels({'$-20\%$', '$0$', '$20\%$', '$40\%$'})
else
    axis([-63 + bias 63 - bias -0.20 0.30]);
    yticks([-0.2, 0, 0.2]);
    yticklabels({'$-20\%$', '$0$', '$20\%$'})
end
xticks([-60.0 + bias, -45 + bias, lat_1 + bias, ...
        lat_2 - bias,  45 - bias, 60.0 - bias])
xticklabels([])
%xticklabels({'-60', num2str(lat_1), num2str(lat_2), '60'})

set([ax1], 'Position', position_1);
set(ax1, 'TickLabelInterpreter', 'latex')
break_axis('handle', ax1, 'position', 0.0, 'axis', 'x', 'length', 0.20)
%ax1.XTicklabels = [];
box off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% moist decomposition, zonal average

ax2 = axes(fig2);
plot(-70:70, zeros(size(-70:70)), 'LineWidth', 0.8, 'Color', [0.3, 0.3, 0.3], 'Linestyle', '-.')
hold on;
%p0_r = plot(plot_x_S + bias, d_omega      (indS));
    
p1_m = plot(plot_x_S + bias, nanmean(d_omega_QG   (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2) ./ omega_QG_h_mean_m(indS));
p2_m = plot(plot_x_S + bias, nanmean(d_Adv_m      (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2) ./ omega_QG_h_mean_m(indS));
%p3_m = plot(plot_x_S + bias, nanmean(d_denom    (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2) ./ omega_QG_h_mean_m(indS));
p3_1_m = plot(plot_x_S + bias, nanmean(d_sigma_m  (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2) ./ omega_QG_h_mean_m(indS));
p3_2_m = plot(plot_x_S + bias, nanmean(d_k2_l2_m  (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2) ./ omega_QG_h_mean_m(indS));
p3_3_m = plot(plot_x_S + bias, nanmean(d_m2_m     (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2) ./ omega_QG_h_mean_m(indS));
if FULL
    p4_1_m = plot(plot_x_S + bias, nanmean(d_excessive_1(indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2) ./ omega_QG_h_mean_m(indS));
    p4_2_m = plot(plot_x_S + bias, nanmean(d_excessive_2(indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2) ./ omega_QG_h_mean_m(indS));
    p5_m   = plot(plot_x_S + bias, nanmean(d_epsilon    (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2) ./ omega_QG_h_mean_m(indS));
    p6_m   = plot(plot_x_S + bias, nanmean(d_rec_m_full (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2) ./ omega_QG_h_mean_m(indS));
else
    p6_m = plot(plot_x_S + bias, nanmean(d_rec_m  (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2) ./ omega_QG_h_mean_m(indS));
end

%set(p0_m, 'LineWidth', 1, 'Color', 'red');
set(p1_m, 'LineWidth', 1.5, 'Color', 'black');
set(p2_m, 'LineWidth', 1.5, 'Color', colors(6, :));
%set(p3_m, 'LineWidth', 1.5, 'Color', 'red');
set(p3_1_m, 'LineWidth', 1.5, 'Color', colors(1, :));
set(p3_2_m, 'LineWidth', 1.5, 'Color', colors(3, :));
set(p3_3_m, 'LineWidth', 1.5, 'Color', colors(5, :));
if FULL
    set(p4_1_m, 'LineWidth', 1.5, 'Color', colors(4, :));
    set(p4_2_m, 'LineWidth', 1.5, 'Color', colors(2, :));
    set(p5_m, 'LineWidth', 1.5, 'Color', colors(7, :));
end
set(p6_m, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', ':');

%p0_m = plot(plot_x_N - bias, d_omega      (indN));
p1_m   = plot(plot_x_N - bias, nanmean(d_omega_QG(indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2) ./ omega_QG_h_mean_m(indN));
p2_m   = plot(plot_x_N - bias, nanmean(d_Adv_m   (indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2) ./ omega_QG_h_mean_m(indN));
%p3_m = plot(plot_x_N - bias, nanmean(d_denom    (indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2) ./ omega_QG_h_mean_m(indN));
p3_1_m = plot(plot_x_N - bias, nanmean(d_sigma_m (indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2) ./ omega_QG_h_mean_m(indN));
p3_2_m = plot(plot_x_N - bias, nanmean(d_k2_l2_m (indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2) ./ omega_QG_h_mean_m(indN));
p3_3_m = plot(plot_x_N - bias, nanmean(d_m2_m    (indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2) ./ omega_QG_h_mean_m(indN));
if FULL
    p4_1_m = plot(plot_x_N - bias, nanmean(d_excessive_1(indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2) ./ omega_QG_h_mean_m(indN));
    p4_2_m = plot(plot_x_N - bias, nanmean(d_excessive_2(indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2) ./ omega_QG_h_mean_m(indN));
    p5_m   = plot(plot_x_N - bias, nanmean(d_epsilon    (indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2) ./ omega_QG_h_mean_m(indN));
    p6_m   = plot(plot_x_N - bias, nanmean(d_rec_m_full (indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2) ./ omega_QG_h_mean_m(indN));
else
    p6_m = plot(plot_x_N - bias, nanmean(d_rec_m  (indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2) ./ omega_QG_h_mean_m(indN));
end

set(p1_m, 'LineWidth', 1.5, 'Color', 'black');
set(p2_m, 'LineWidth', 1.5, 'Color', colors(6, :));
%set(p3_m, 'LineWidth', 1.5, 'Color', 'red');
set(p3_1_m, 'LineWidth', 1.5, 'Color', colors(1, :));
set(p3_2_m, 'LineWidth', 1.5, 'Color', colors(3, :));
set(p3_3_m, 'LineWidth', 1.5, 'Color', colors(5, :));
if FULL
    set(p4_1_m, 'LineWidth', 1.5, 'Color', colors(4, :));
    set(p4_2_m, 'LineWidth', 1.5, 'Color', colors(2, :));
    set(p5_m, 'LineWidth', 1.5, 'Color', colors(7, :));
end
set(p6_m, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', ':');



if strfind(pwd_str, 'GFDL')
    axis([-63 + bias 63 - bias -0.20 0.50]);
    yticks([-0.2, 0, 0.2, 0.4]);
    yticklabels({'$-20\%$', '$0$', '$20\%$', '$40\%$'})
else
    axis([-63 + bias 63 - bias -0.20 0.30]);
    yticks([-0.2, 0, 0.2]);
    yticklabels({'$-20\%$', '$0$', '$20\%$'})
end
xticks([-60.0 + bias, -45 + bias, lat_1 + bias, ...
        lat_2 - bias, 45 - bias, 60.0 - bias]);
xticklabels({'-60', '-45', num2str(lat_1), num2str(lat_2), '45', '60'})

if FULL
    lgd = legend( ...
    [p0_d, p1_m, p6_m, p3_d, p3_2_m, p3_1_m, p3_3_m, p2_m, p7_d, p4_1_m, p4_2_m, p5_m], ...
    {'$\overline{\omega}$' , ...
     '$\overline{\omega}_{QG}$' , ...
     '$\overline{\omega}_{QG-rec}$', ...
     '$\overline{J}$'...
     '${k}^2$ and ${k}_J^2$', ...
     '$\hat{\sigma}$ or $\hat{\sigma}_m$', ...
     '${{m}^2}$', ...
     '$\overline{Adv}$' , ...
     '$\overline{{\frac{T}{\theta}\frac{d\theta}{dp}\big|_{\theta*}\omega}}$' , ...
     'Unbal$._{\frac{d\theta}{dp}\big|_{\theta*}}$' , ...
     'Unbal$._{\omega-\omega_{QG}}$' , ...
     '$\overline{\epsilon}$'}, ...
    'location', 'best', ...
    'interpreter', 'latex');
    lgd.Position = [left_margin + width + 0.07, 0.10, 0.08, 0.80];
else
    lgd = legend( ...
    [p0_d, p1_m, p6_m, p3_d, p3_2_m, p3_1_m, p3_3_m, p2_m, p7_d], ...
    {'$\overline{\omega}$' , ...
     '$\overline{\omega}_{QG}$' , ...
     '$\overline{\omega}_{QG-rec}$', ...
     '$\overline{J}$'...
     '${k}^2$ and ${k}_J^2$', ...
     '$\hat{\sigma}$ or $\hat{\sigma}_m$', ...
     '${m}^2$', ...
     '$\overline{Adv}$' , ...
     '$\overline{{\frac{T}{\theta}\frac{d\theta}{dp}\big|_{\theta*}\omega}}$'}, ...
     'location', 'best', ...
     'interpreter', 'latex');
    lgd.Position = [left_margin + width + 0.07, 0.15, 0.08, 0.66];
end
set(lgd,'Box','off')

title(['Moist Decomposition'], 'FontWeight', 'normal');
xlabel('Latitude');
ylabel('$\Delta\omega/\omega$', 'Interpreter', 'latex');
ylab = get(ax2, 'ylabel');
%set(ylab, 'rotation', 0)
set(ylab,'Units','normalized'); 
ylab.Position(2) = ylab.Position(2) / height * (2 * height + distance);

set([ax2], 'Position', position_2);
set(ax2, 'TickLabelInterpreter', 'latex')
break_axis('handle', ax2, 'position', 0.0, 'axis', 'x', 'length', 0.20)
hold off;
box off;

figname = [plot_path_ensemble, 'omega_decomposition_moist_fancy_', plot_level_name];
if FULL
    figname = [figname, '_full'];
end

if EPS
    print(fig2, figname, '-depsc', '-r0', '-painters')
else
    saveas(fig2, figname, 'png')
end




