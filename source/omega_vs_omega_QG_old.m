
EPS = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot omega vs omega_QG

figsize_multipanel = [10, 10, 400, 300];
fig = figure('pos', figsize_multipanel);
hold on
row = 2;
col = 1;
ax = gca;

% parameters of the plot
left_margin     = 0.08;
right_margin    = 0.10;
down_margin     = 0.09;
up_margin       = 0.09;
margins = [left_margin, right_margin, down_margin, up_margin];
lat_mask = 7;
ZERO = false;
FLIP = false;
LIGHT_COLOR = false;
ddd = 1 - down_margin - up_margin;
cb_pos = [1 - right_margin + 0.02, ddd/6 + down_margin, 0.025, ddd/6*4];
cb_string = 'Pa s$^{-1}$';
cb_label_pos_bias = [0.15, -1.70];
cb_ticks = -1.4:0.2:0;
cb_ticklabels = {};
climit = [-1.4, 0];
[ax1, ax1_1, ax1_2, ax2, ax3] = deal({});

temp = omega_h .* NaN_matrix_d;
plot_title = '(a) $\overline{\omega}$';
%[ax1, ax1_1, ax1_2, ax2, ax3] = subplot_2D_map(ax1, ax1_1, ax1_2, ax2, ax3, ...
subplot_2D_map(...
    @parula, margins, ax, row, col, 1, figsize_multipanel, ...
    temp, X, Y, plot_title, SMOOTH, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = omega_QG_h .* NaN_matrix_d;
plot_title = '(b) $\overline{\omega}_{QG}$';
%[ax1, ax1_1, ax1_2, ax2, ax3] = subplot_2D_map(ax1, ax1_1, ax1_2, ax2, ax3, ...
subplot_2D_map(...
    @parula, margins, ax, row, col, 2, figsize_multipanel, ...
    temp, X, Y, plot_title, SMOOTH, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

set(gca,'TickLabelInterpreter','latex')
figname = [plot_path_ensemble, 'omega_vs_omega_QG'];
if EPS
    print(fig, figname, '-depsc', '-r0', '-painters')
else
    print(fig, figname, '-dpng', '-r300');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot contributions of k2 and l2
figsize_multipanel = [10, 10, 400, 300];
fig = figure('pos', figsize_multipanel);
hold on
row = 2;
col = 1;
ax = gca;

% parameters of the plot
left_margin     = 0.08;
right_margin    = 0.10;
down_margin     = 0.09;
up_margin       = 0.09;
margins = [left_margin, right_margin, down_margin, up_margin];
lat_mask = 7;
ZERO = true;
FLIP = true;
LIGHT_COLOR = false;
ddd = 1 - down_margin - up_margin;
cb_pos = [1 - right_margin + 0.02, ddd/6 + down_margin, 0.025, ddd/6*4];
cb_string = '\%K$^{-1}$';
cb_label_pos_bias = [0.15, -0.395];
cb_ticks = -0.12:0.06:0.12;
cb_ticklabels = {'-12', '-6', '0', '6', '12'};
climit = [-0.15, 0.15];
[ax1, ax1_1, ax1_2, ax2, ax3] = deal({});

temp = d_k2 ./ omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(a) $k$';
%[ax1, ax1_1, ax1_2, ax2, ax3] = subplot_2D_map(ax1, ax1_1, ax1_2, ax2, ax3, ...
subplot_2D_map(...
    @jet, margins, ax, row, col, 1, figsize_multipanel, ...
    temp, X, Y, plot_title, SMOOTH, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

temp = d_l2 ./ omega_QG_h_mean_d / delta_T .* NaN_matrix_d;
plot_title = '(b) $k_J$';
%[ax1, ax1_1, ax1_2, ax2, ax3] = subplot_2D_map(ax1, ax1_1, ax1_2, ax2, ax3, ...
subplot_2D_map(...
    @jet, margins, ax, row, col, 2, figsize_multipanel, ...
    temp, X, Y, plot_title, SMOOTH, climit, num_event_d, num_threshold, ...
    cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
    lat_mask, ZERO, FLIP, LIGHT_COLOR);

set(gca,'TickLabelInterpreter','latex')
figname = [plot_path_ensemble, 'dry_d_k2_vs_dry_d_kJ2'];
if EPS
    print(fig, figname, '-depsc', '-r0', '-painters')
else
    print(fig, figname, '-dpng', '-r300');
end

