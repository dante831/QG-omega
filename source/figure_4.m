
if ~isempty(strfind(filename, 'CESM_6hourly'))
    EPS = true;
else
    EPS = false;
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
SMOOTH_2 = true;

temp = d_omega ./ Omega_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(a) $\overline{\omega}$';
main_fig = subplot_2D_map(main_fig, @jet, margins, 1, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_omega_QG ./ Omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(b) $\overline{\omega}_{QG}$';
main_fig = subplot_2D_map(main_fig, @jet, margins, 2, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_rec ./ Omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(c) All contributions';
main_fig = subplot_2D_map(main_fig, @jet, margins, 3, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_J ./ Omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(d) Diabatic heating $\overline{J}$';
main_fig = subplot_2D_map(main_fig, @jet, margins, 4, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_sigma ./ Omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(e) Static stability $\hat{\sigma}$';
main_fig = subplot_2D_map(main_fig, @jet, margins, 5, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_k2_l2 ./ Omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(f) Horizontal wavenumbers ${k},{k}_J$';
main_fig = subplot_2D_map(main_fig, @jet, margins, 6, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_m2 ./ Omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(g) Vertical wavenumber ${m}$';
main_fig = subplot_2D_map(main_fig, @jet, margins, 7, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_Adv ./ Omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
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


