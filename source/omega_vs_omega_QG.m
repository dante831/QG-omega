
if ~isempty(strfind(pwd_str, 'GFDL')) || ...
   ~isempty(strfind(pwd_str, 'daily')) 
    EPS = false;    
else
    EPS = true;    
end

figsize_multipanel = [10, 10, 300, 250];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot omega vs omega_QG

main_fig_1 = fig_class;
main_fig_1.handle = figure('pos', figsize_multipanel);
main_fig_1.row = 2;
main_fig_1.col = 1;
main_fig_1.axes = {};

% parameters of the plot
left_margin     = 0.08;
right_margin    = 0.10;
down_margin     = 0.09;
up_margin       = 0.09;
distance        = 30;
margins = [left_margin, right_margin, down_margin, up_margin];
lat_mask = 3;
ZERO = false;
FLIP = false;
LIGHT_COLOR = false;
ddd = 1 - down_margin - up_margin;
cb_pos = [1 - right_margin + 0.02, ddd/6 + down_margin, 0.025, ddd/6*4];
cb_string = 'Pa s$^{-1}$';
cb_label_pos_bias = [0.50, -1.90];
cb_ticks = -1.5:0.3:0;
cb_ticklabels = {};
climit = [-1.6, 0.1];
[ax1, ax1_1, ax1_2, ax2, ax3] = deal({});

temp = omega_h .* NaN_matrix_d;
plot_title = '(a) $\overline{\omega}$';
main_fig_1 = subplot_2D_map(...
    main_fig_1, @parula, margins, 1, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = omega_QG_h .* NaN_matrix_d;
plot_title = '(b) $\overline{\omega}_{QG}$';
main_fig_1 = subplot_2D_map(...
    main_fig_1, @parula, margins, 2, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

positions = set_axes(main_fig_1, figsize_multipanel, margins, distance, distance);
figname = [plot_path_ensemble, 'omega_vs_omega_QG'];

if EPS
    print(main_fig_1.handle, figname, '-depsc', '-r0', '-painters')
else
    print(main_fig_1.handle, figname, '-dpng', '-r300');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot contributions of k2 and l2

main_fig_2 = fig_class;
main_fig_2.handle = figure('pos', figsize_multipanel);
main_fig_2.row = 2;
main_fig_2.col = 1;
main_fig_2.axes = {};

% parameters of the plot
left_margin     = 0.08;
right_margin    = 0.10;
down_margin     = 0.09;
up_margin       = 0.09;
distance        = 30;
margins = [left_margin, right_margin, down_margin, up_margin];
lat_mask = 3;
ZERO = true;
FLIP = true;
LIGHT_COLOR = false;
ddd = 1 - down_margin - up_margin;
cb_pos = [1 - right_margin + 0.02, ddd/6 + down_margin, 0.025, ddd/6*4];
cb_string = '\%K$^{-1}$';
cb_label_pos_bias = [0.50, -0.390];
cb_ticks = -0.12:0.06:0.12;
cb_ticklabels = {'-12', '-6', '0', '6', '12'};
climit = [-0.15, 0.15];
[ax1, ax1_1, ax1_2, ax2, ax3] = deal({});

temp = d_k2 ./ omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(a) $k$';
main_fig_2 = subplot_2D_map(...
    main_fig_2, @jet, margins, 1, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_l2 ./ omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(b) $k_J$';
main_fig_2 = subplot_2D_map(...
    main_fig_2, @jet, margins, 2, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

positions = set_axes(main_fig_2, figsize_multipanel, margins, distance, distance);
figname = [plot_path_ensemble, 'dry_d_k2_vs_dry_d_kJ2'];

if EPS
    print(main_fig_2.handle, figname, '-depsc', '-r0', '-painters')
else
    print(main_fig_2.handle, figname, '-dpng', '-r300');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the approximated dry decomposition

main_fig_3 = fig_class;
main_fig_3.handle = figure('pos', figsize_multipanel);
main_fig_3.row = 2;
main_fig_3.col = 1;
main_fig_3.axes = {};

% parameters of the plot
left_margin     = 0.08;
right_margin    = 0.10;
down_margin     = 0.09;
up_margin       = 0.09;
distance        = 30;
margins = [left_margin, right_margin, down_margin, up_margin];
lat_mask = 3;
ZERO = true;
FLIP = true;
LIGHT_COLOR = false;
ddd = 1 - down_margin - up_margin;
cb_pos = [1 - right_margin + 0.02, ddd/6 + down_margin, 0.025, ddd/6*4];
cb_string = '\%K$^{-1}$';
cb_label_pos_bias = [0.50, -0.390];
cb_ticks = -0.12:0.06:0.12;
cb_ticklabels = {'-12', '-6', '0', '6', '12'};
climit = [-0.15, 0.15];
[ax1, ax1_1, ax1_2, ax2, ax3] = deal({});

temp = (omega_QG_r - omega_QG_h) ./ omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(a) $\overline{\omega}_{QG}$';
main_fig_3 = subplot_2D_map(...
    main_fig_3, @jet, margins, 1, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = - (k2_r - k2_h) ./ k2_h + (l2_r - l2_h) ./ l2_h - (sigma_r - sigma_h) ./ sigma_h + (J_r - J_h) ./ J_h;
temp = temp / delta_T .* NaN_matrix_d;
plot_title = '(b) Approximate dry decomposition';
main_fig_3 = subplot_2D_map(...
    main_fig_3, @jet, margins, 2, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

positions = set_axes(main_fig_3, figsize_multipanel, margins, distance, distance);
figname = [plot_path_ensemble, 'approx_drydecomp'];

print(main_fig_3.handle, figname, '-dpng', '-r300');





