

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% contribution of k and kJ, 2D map

figsize_multipanel = [10, 10, 350, 150];
main_fig = fig_class;
main_fig.handle = figure('pos', figsize_multipanel);
main_fig.row = 1;
main_fig.col = 1;
main_fig.axes = {};

% parameters of the plot
left_margin     = 0.09;
right_margin    = 0.11;
down_margin     = 0.13;
up_margin       = 0.12;
distance        = 30;
margins = [left_margin, right_margin, down_margin, up_margin];
lat_mask = 0;
ZERO = true;
FLIP = true;
LIGHT_COLOR = false;
ddd = 1 - down_margin - up_margin;
cb_pos = [1 - right_margin + 0.02, ddd/10 + down_margin, 0.025, ddd/10*8];
cb_string = '\%K$^{-1}$';
cb_label_pos_bias = [- cb_pos(1) + 1.5, - cb_pos(2) - 0.03];
%cb_label_pos_bias = [];
cb_ticks = -0.02:0.01:0.02;
cb_ticklabels = {'-2', '-1', '0', '1', '2'};
climit = [-0.025, 0.025];

temp = d_k2_l2(:, :, plot_ind) ./ Omega_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = 'Combined contribution of $k^2$ and $k_J^2$';
main_fig = subplot_2D_map(main_fig, @jet, margins, 1, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

positions = set_axes(main_fig, figsize_multipanel, margins, distance, distance);
figname = [plot_path_ensemble, 'dry_d_k2_l2'];

saveas(main_fig.handle, figname, 'png')

