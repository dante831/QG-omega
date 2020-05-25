
if ~isempty(strfind(pwd_str, 'GFDL')) || ...
   ~isempty(strfind(pwd_str, 'daily'))
    EPS = false;
else
    EPS = true;
end

term2_h = sigma_h .* (- k2_h) .* omega_QG_h;
term2_r = sigma_r .* (- k2_r) .* omega_QG_r;
term3_h = C_h;
term3_r = C_r;
term4_h = Adv_h;
term4_r = Adv_r;
term5_h = F0.^2 .* (- m2_h) .* omega_QG_h;
term5_r = F0.^2 .* (- m2_r) .* omega_QG_r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% terms of the QG-omega equation, 2D map

figsize_multipanel = [10, 10, 600, 260];
main_fig = fig_class;
main_fig.handle = figure('pos', figsize_multipanel);
main_fig.row = 2;
main_fig.col = 2;
main_fig.axes = {};

% parameters of the plot
left_margin     = 0.05;
right_margin    = 0.085;
down_margin     = 0.08;
up_margin       = 0.08;
distance        = 30;
factor          = 1e16;
margins = [left_margin, right_margin, down_margin, up_margin];
lat_mask = 3;
ZERO = true;
FLIP = false;
LIGHT_COLOR = false;
canvas_h = 1 - down_margin - up_margin;
panel_height = ((1 - down_margin - up_margin) - (main_fig.row - 1) * distance / figsize_multipanel(4)) / main_fig.row;
panel_width = ((1 - left_margin - right_margin) - (main_fig.col - 1) * distance / figsize_multipanel(3)) / main_fig.col;

colorbar_height = panel_height;
colorbar_width = 0.02;
cmap = @jet;

f_cb_bias = @(cb_pos, climit) [- cb_pos(1) - 1.0, - cb_pos(2) + climit(2) + 0.15*(climit(2) - climit(1))];

% colorbar
cb_string = '$\times10^{-16}$m kg$^{-1}$s$^{-1}$';
cb_ticks = (-4e-16:2e-16:4e-16) * factor;
cb_ticklabels = {'-4', '-2', '0', '2', '4'};
climit = [-5e-16, 5e-16] * factor;
cb_pos = [1 - right_margin + 0.02, down_margin + panel_height * 1 + distance/figsize_multipanel(4) * 1, colorbar_width, colorbar_height];
cb_label_pos_bias = f_cb_bias(cb_pos, climit);

temp = term2_h * factor;
plot_title = '(a) $\overline{\nabla^2(\sigma\omega)}$';
main_fig = subplot_2D_map(main_fig, cmap, margins, 1, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR, 'colorbar', true);

temp = - term3_h * factor;
plot_title = '(b) $\frac{\kappa}{p}\overline{\nabla^2J}$';
main_fig = subplot_2D_map(main_fig, cmap, margins, 2, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

% colorbar
cb_string = '$\times10^{-17}$m kg$^{-1}$s$^{-1}$';
cb_ticks = (-8e-17:4e-17:8e-17) * factor;
cb_ticklabels = {'-8', '-4', '0', '4', '8'};
climit = [-10e-17, 10e-17] * factor;
cb_pos = [1 - right_margin + 0.02, down_margin, colorbar_width, colorbar_height];
cb_label_pos_bias = f_cb_bias(cb_pos, climit);

temp = term5_h * factor;
plot_title = '(c) $f_0^2\overline{\partial_p^2\omega_{QG}}$';
main_fig = subplot_2D_map(main_fig, cmap, margins, 3, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR, 'colorbar', true);

cb_pos = [1 - right_margin + 0.02, down_margin, colorbar_width, colorbar_height];
cb_label_pos_bias = f_cb_bias(cb_pos, climit);
temp = - term4_h * factor;
plot_title = '(d) $-\overline{Adv}$';
main_fig = subplot_2D_map(main_fig, cmap, margins, 4, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

positions = set_axes(main_fig, figsize_multipanel, margins, distance, distance);
figname = [plot_path_ensemble, 'terms_', plot_level_name, '_2D'];

if EPS
    print(main_fig.handle, figname, '-depsc', '-r0', '-painters')
else
    saveas(main_fig.handle, figname, 'png')
end


