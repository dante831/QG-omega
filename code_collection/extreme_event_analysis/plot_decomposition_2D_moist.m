
if ~isempty(strfind(pwd_str, 'GFDL')) || ...
   ~isempty(strfind(pwd_str, 'daily'))
    EPS = false;
else
    EPS = true;
end

if GFDL
    SMOOTH_2_backup = SMOOTH_2;
    SMOOTH_2 = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% terms of the QG-omega equation, 2D map

figsize_multipanel = [10, 10, 600, 500];
main_fig = fig_class;
main_fig.handle = figure('pos', figsize_multipanel);
main_fig.row = 4;
main_fig.col = 2;
main_fig.axes = {};

% parameters of the plot
left_margin     = 0.05;
right_margin    = 0.085;
down_margin     = 0.05;
up_margin       = 0.05;
distance        = 30;
factor          = 1.0;
margins = [left_margin, right_margin, down_margin, up_margin];
lat_mask = 3;
ZERO = true;
FLIP = true;
LIGHT_COLOR = false;
canvas_h = 1 - down_margin - up_margin;
panel_height = ((1 - down_margin - up_margin) - (main_fig.row - 1) * distance / figsize_multipanel(4)) / 4;
panel_width = ((1 - left_margin - right_margin) - (main_fig.col - 1) * distance / figsize_multipanel(3)) / 4;

colorbar_height = panel_height;
colorbar_width = 0.02;

f_cb_bias = @(cb_pos, climit) [- cb_pos(1) + 0.1, - cb_pos(2) + climit(1) - 0.05*(climit(2) - climit(1))];

% colorbar
cb_pos = [1 - right_margin + 0.02, canvas_h/4 + down_margin, 0.02, canvas_h/4*2];
cb_string = '\%K$^{-1}$';
cb_label_pos_bias = [-0.18, -0.43];
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

temp = d_rec_m_full_paper(:, :, plot_ind) ./ Omega_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(c) All contributions';
main_fig = subplot_2D_map(main_fig, @jet, margins, 3, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR, 'colorbar', true);

temp = d_J_res(:, :, plot_ind) ./ Omega_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(d) Residual heating $\overline{J}_{res}$';
main_fig = subplot_2D_map(main_fig, @jet, margins, 4, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_sigma_m(:, :, plot_ind) ./ Omega_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(e) Moist stability $\hat{\sigma}_m$';
main_fig = subplot_2D_map(main_fig, @jet, margins, 5, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_k2_l2_m(:, :, plot_ind) ./ Omega_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(f) Horizontal wavenumbers $k$, $k_J$';
main_fig = subplot_2D_map(main_fig, @jet, margins, 6, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_m2_m(:, :, plot_ind) ./ Omega_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(g) Vertical wavenumber $m$';
main_fig = subplot_2D_map(main_fig, @jet, margins, 7, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_Adv_m(:, :, plot_ind) ./ Omega_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(h) Advective forcing $\overline{Adv}$';
main_fig = subplot_2D_map(main_fig, @jet, margins, 8, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

%{
%temp = (d_J_res(:, :, plot_ind) - d_epsilon(:, :, plot_ind)) ./ Omega_h_mean_d / delta_T .* NaN_matrix_d;
temp = d_excessive(:, :, plot_ind) ./ Omega_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(g) Residual heating omega part';
main_fig = subplot_2D_map(main_fig, @jet, margins, 7, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_epsilon(:, :, plot_ind) ./ Omega_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(h) Residual heating $\epsilon$ part';
main_fig = subplot_2D_map(main_fig, @jet, margins, 8, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);
%}
positions = set_axes(main_fig, figsize_multipanel, margins, distance, distance);
figname = [plot_path_ensemble, 'moist_decomposition_', plot_level_name];

if EPS
    print(main_fig.handle, figname, '-depsc', '-r0', '-painters')
else
    saveas(main_fig.handle, figname, 'png')
end

if GFDL
    SMOOTH_2 = SMOOTH_2_backup;
end
