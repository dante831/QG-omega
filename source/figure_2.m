
if ~isempty(strfind(filename, 'CESM_6hourly'))
    EPS = true;
else
    EPS = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dry decomposition, 2D map

figsize_multipanel = [10, 10, 600, 500];
fig = figure('pos', figsize_multipanel);
hold on
row = 4;
col = 2;

% parameters of the plot
left_margin     = 0.05;
right_margin    = 0.09;
down_margin     = 0.09;
up_margin       = 0.05;
margins = [left_margin, right_margin, down_margin, up_margin];
lat_mask = 7;
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

ax = gca;
temp = d_omega ./ omega_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(a) $\overline{\omega}$';
subplot_2D_map(@jet, margins, ax, row, col, 1, figsize_multipanel, ...
    temp, X, Y, plot_title, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_omega_QG ./ omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(b) $\overline{\omega}_{QG}$';
subplot_2D_map(@jet, margins, ax, row, col, 2, figsize_multipanel, ...
    temp, X, Y, plot_title, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_rec ./ omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(c) All contributions';
subplot_2D_map(@jet, margins, ax, row, col, 3, figsize_multipanel, ...
    temp, X, Y, plot_title, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_J ./ omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(d) Diabatic heating $\overline{J}$';
subplot_2D_map(@jet, margins, ax, row, col, 4, figsize_multipanel, ...
    temp, X, Y, plot_title, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_sigma ./ omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(e) Static stability $\hat{\sigma}$';
subplot_2D_map(@jet, margins, ax, row, col, 5, figsize_multipanel, ...
    temp, X, Y, plot_title, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_k2_l2 ./ omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(f) Horizontal wavenumbers ${k},{k}_J$';
subplot_2D_map(@jet, margins, ax, row, col, 6, figsize_multipanel, ...
    temp, X, Y, plot_title, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_m2 ./ omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(g) Vertical wavenumber ${m}$';
subplot_2D_map(@jet, margins, ax, row, col, 7, figsize_multipanel, ...
    temp, X, Y, plot_title, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_Adv ./ omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(h) Advective forcing $\overline{Adv}$';
subplot_2D_map(@jet, margins, ax, row, col, 8, figsize_multipanel, ...
    temp, X, Y, plot_title, SMOOTH_2, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

set(gca,'TickLabelInterpreter','latex')

figname = [plot_path_ensemble, 'dry_decomposition_', plot_level_name];

addpath('/disk7/ziweili/CESM_LENS/altmany-export_fig-9676767/')

if EPS
    print(fig, figname, '-depsc', '-r0', '-painters')
else
    saveas(fig, figname, 'png')
end


