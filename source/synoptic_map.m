
%%%%%%%%%%%%%%%%%%%%%%%% plot a 'synoptic' map around the event %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get 3D fields

nc_filename = './data/example_event_synoptic_map.nc';
[T_s, ua_s, va_s, ug_s, vg_s, omega_s, omega_b_s, lon_s, lat_s, level_s, event_timespan] = read_ncfile(nc_filename);
box_max_y_2 = size(T_s, 1); box_max_x_2 = size(T_s, 2);

lon_indices_plot = event_lonspan((1+end)/2) - (box_max_x_2 - 1)/2 : event_lonspan((1+end)/2) + (box_max_x_2 - 1)/2;
lat_indices_plot = event_latspan((1+end)/2) - (box_max_y_2 - 1)/2 : event_latspan((1+end)/2) + (box_max_y_2 - 1)/2;

% smooth T field
T_smoothed = T_s;
for k = 1 : length(level)
    for t = 1 : length(event_timespan)
        T_smoothed(:, :, k, t) = smooth2a(T_s(:, :, k, t), window_x, window_y);
        % renormalize so that the smoothed field has the same variance as the previous field
        T_smoothed(:, :, k, t) = smooth_normalize(T_s(:, :, k, t), T_smoothed(:, :, k, t));
    end
end
T_s = T_smoothed;

lonmax_plot = lon_series_original(lon_indices_plot(end));
lonmin_plot = lon_series_original(lon_indices_plot(1));
latmax_plot = lat_series_original(lat_indices_plot(end));
latmin_plot = lat_series_original(lat_indices_plot(1));
dlat = lat_series_original(2) - lat_series_original(1);
dlon = lon_series_original(2) - lon_series_original(1);
lat_plot = lat_series_original(lat_indices_plot);
lon_plot = lon_series_original(lon_indices_plot);

% define location matrices
[X, Y] = meshgrid(lon_plot, lat_plot);
[X_inversion, Y_inversion] = meshgrid(lon_plot(((box_max_x_2 - computation_x) / 2 + 1) : (box_max_x_2 + computation_x) / 2), ...
                                      lat_plot(((box_max_y_2 - computation_y) / 2 + 1) : (box_max_y_2 + computation_y) / 2));
[Lat2, Lon2] = meshgrid(lat_plot, lon_plot);

% get Q_vector
phi_plot = lat_plot / 180.0 * 3.1415926;
lambda_plot = lon_plot / 180.0 * 3.1415926;
[~, Phi_s] = meshgrid(lambda_plot, phi_plot);
[A_s, B_s] = Q_vector(level, ug_s, vg_s, T_s, Phi_s, event_timespan, dphi, dlambda, f0, beta, GEO_T);
A_s = -A_s; % get positive divergence

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start to plot the figure

figsize_multipanel = [10, 10, 600, 450];
main_fig = fig_class;
main_fig.row = 2;
main_fig.col = 2;
main_fig.axes = {};
main_fig.handle = figure('pos', figsize_multipanel);
plot_level = 50000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw precipitation time series on the left, and 3D contour of omega_QG on the right

%main_fig = threeD_event_plot(main_fig, 0, historical_quantile, ...
%                    lat, lon, level, omega_QG, event_precip, plot_level, event1, ...
%                    lon_series_original(event_lon), lat_series_original(event_lat), true);
main_fig = threeD_event_plot(main_fig, 0, historical_quantile, ...
                    lat, lon, level, omega_QG, event_precip(x_ind, y_ind, :), plot_level, event_timespan, ...
                    lon(omega_x_ind), lat(omega_y_ind), true);

%%%%%%%%%%%%%%%%%%%%%%%% synoptic map #1, surface pressure and omega %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_ind = 3;
ax1 = axes(main_fig.handle);
main_fig.axes{fig_ind} = [ax1];

hold on; box on;

% plot 500hPa omega
n_level = 20;
climit = [-2.0, 2.0];
inc = (climit(2) - climit(1)) / n_level;
levs = climit(1) : inc : climit(2);

[c1_1, h1_1] = contourf(ax1, X, Y, omega_s(:, :, level_s == plot_level, 2), levs, 'Linestyle', 'none');
cb = colorbar(ax1);
cb.Label.Interpreter = 'latex';
cb.Label.String = 'Pa s$^{-1}$';
cb.Label.Rotation = 0;
cpos = get(cb,'Position');
cb.Label.Position = cpos(1:2) + [-0.30, climit(1) - 0.25];
cb.Ticks = -1.8:0.6:1.8;
cb.TickLabelInterpreter = 'latex';
set(cb, 'FontSize', 10);

% get two opposite colors from a specified symmetric colormap
addpath('./source/cbrewer/')
colors_1 = colormap(ax1, cbrewer('div', 'RdBu', 26, 'pchip'));
colors_1 = colors_1([4:13, 14:23], :);
colors_1(n_level / 2 : n_level / 2 + 1, :) = 0.99;
colormap(ax1, flipud(colors_1))
caxis(ax1, climit)

set(ax1,'TickDir','out');
set(ax1, 'TickLabelInterpreter', 'latex')
set(ax1, 'xtick', [170:5:185])
set(ax1, 'ytick', [35:5:50])
set(ax1, 'Xticklabel', {'170$^\circ$E', '175$^\circ$E', '180$^\circ$', '5$^\circ$W'})
set(ax1, 'Yticklabel', {'35$^\circ$N', '40$^\circ$N', '45$^\circ$N', '50$^\circ$N'})
title(ax1, '(c) $\omega$ at 500hPa and surface pressure anomaly', 'interpreter', 'latex', 'FontSize', 11)

% plot inversion domain
% left
plot(ax1, ones(computation_x, 1) * lon_series_original(lon_indices_plot((box_max_x_2 - computation_x)/2 + 1)), ...
          lat_series_original(lat_indices_plot((box_max_y_2 - computation_y)/2 + 1: (box_max_y_2 + computation_y)/2)), ...
          'color', [0.5, 0.5, 0.5], 'linewidth', 1.5)
% right
plot(ax1, ones(computation_x, 1) * lon_series_original(lon_indices_plot((box_max_x_2 + computation_x)/2)), ...
          lat_series_original(lat_indices_plot((box_max_y_2 - computation_y)/2 + 1: (box_max_y_2 + computation_y)/2)), ...
          'color', [0.5, 0.5, 0.5], 'linewidth', 1.5)
% up
plot(ax1, lon_series_original(lon_indices_plot((box_max_x_2 - computation_x)/2 + 1: (box_max_x_2 + computation_x)/2)), ...
          ones(computation_y, 1) * lat_series_original(lat_indices_plot((box_max_y_2 + computation_y)/2)), ...
          'color', [0.5, 0.5, 0.5], 'linewidth', 1.5)
% down
plot(ax1, lon_series_original(lon_indices_plot((box_max_x_2 - computation_x)/2 + 1: (box_max_x_2 + computation_x)/2)), ...
          ones(computation_y, 1) * lat_series_original(lat_indices_plot((box_max_y_2 - computation_y)/2 + 1)), ...
          'color', [0.5, 0.5, 0.5], 'linewidth', 1.5)

% plot coastlines
load coastlines;
plot(ax1, [coastlon; NaN; coastlon + 360; NaN; coastlon + 720 ], [coastlat; NaN; coastlat; NaN; coastlat], 'Color', 'black');

ax2 = axes(main_fig.handle);
main_fig.axes{fig_ind} = [main_fig.axes{fig_ind}, ax2];

levels_2 = [-3 : 3] * 10; % in hPa
[c2, h2] = contour(ax2, Lon2, Lat2, PS(:, :, 2)/100, levels_2, 'ShowText', 'on', 'linewidth', 1.0); % in hPa
clabel(c2, h2, levels_2, 'interpreter', 'latex', 'LabelSpacing', 100, 'FontSize', 8);
colors_2 = colormap(ax2, cbrewer('seq', 'Greens', 26, 'pchip'));
colormap(ax2, flipud(colors_2));

axis(ax1, [lonmin_plot - 0.2, lonmax_plot + 0.2, latmin_plot - 0.2, latmax_plot + 0.2]);
axis(ax2, [lonmin_plot - 0.2, lonmax_plot + 0.2, latmin_plot - 0.2, latmax_plot + 0.2]);

ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];

%%%%%%%%%%%%%%%%%%%%%%%% synoptic map #2, 500hPa omega_QG and Q-vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig_ind = 4;
ax3 = axes(main_fig.handle);
main_fig.axes{fig_ind} = ax3;
hold on; box on;

% plot 500hPa omega_QG
n_level = 20;
C_TERM = false;
if C_TERM
    climit = [-3.0, 3.0]*1e-16;
    inc = (climit(2) - climit(1)) / n_level;
    levs = climit(1) : inc : climit(2);
    [c1_1, h1_1] = contourf(ax3, X_inversion, Y_inversion, ...
            C(:, :, level == plot_level, 2), levs, 'Linestyle', 'none');
else
    climit = [-2.0, 2.0];
    inc = (climit(2) - climit(1)) / n_level;
    levs = climit(1) : inc : climit(2);
    [c1_1, h1_1] = contourf(ax3, X_inversion, Y_inversion, ...
            omega_QG(:, :, level == plot_level, 2), levs, 'Linestyle', 'none');
end
cb = colorbar(ax3);
cb.Label.Interpreter = 'latex';
cb.Label.Rotation = 0;
cpos = get(cb,'Position');
cb.Label.Position = cpos(1:2) + [-0.30, climit(1) - 0.25];
cb.Label.String = 'Pa s$^{-1}$';
cb.Ticks = -1.8:0.6:1.8;
if C_TERM
    cb.Label.String = 'Pa s$^{-1}$';
    cb.Ticks = (-2.7:0.9:2.7)*1e-16;
    cb.Label.Position = cpos(1:2) + [-0.30, climit(1) - 0];
end
cb.TickLabelInterpreter = 'latex';
set(cb, 'FontSize', 10);

colormap(ax3, flipud(colors_1))
caxis(ax3, climit)
set(ax3,'TickDir','out');
set(ax3, 'TickLabelInterpreter', 'latex')
%set(ax3, 'xtick', [160:10:200])
%set(ax3, 'ytick', [20:10:70])
%set(ax3, 'Xticklabel', {'160$^\circ$E', '170$^\circ$E', '180$^\circ$', '10$^\circ$W', '20$^\circ$W'})
%set(ax3, 'Yticklabel', {'20$^\circ$N', '30$^\circ$N', '40$^\circ$N', '50$^\circ$N', '60$^\circ$N', '70$^\circ$N'})
set(ax3, 'xtick', [170:5:185])
set(ax3, 'ytick', [35:5:50])
set(ax3, 'Xticklabel', {'170$^\circ$E', '175$^\circ$E', '180$^\circ$', '5$^\circ$W'})
set(ax3, 'Yticklabel', {'35$^\circ$N', '40$^\circ$N', '45$^\circ$N', '50$^\circ$N'})
if C_TERM
    title(ax3, '(d) $-\frac{\kappa}{p}\nabla^2J$ and $Adv$ at 500hPa', 'interpreter', 'latex', 'FontSize', 11)
else
    title(ax3, '(d) $\omega_{QG}$ and $2\nabla\cdot\mathbf{Q}$ at 500hPa', 'interpreter', 'latex', 'FontSize', 11)
end
% plot inversion domain
% left
plot(ax3, ones(computation_x, 1) * lon_series_original(lon_indices_plot((box_max_x_2 - computation_x)/2 + 1)), ...
          lat_series_original(lat_indices_plot((box_max_y_2 - computation_y)/2 + 1: (box_max_y_2 + computation_y)/2)), ...
          'color', [0.5, 0.5, 0.5], 'linewidth', 1.5)
% right
plot(ax3, ones(computation_x, 1) * lon_series_original(lon_indices_plot((box_max_x_2 + computation_x)/2)), ...
          lat_series_original(lat_indices_plot((box_max_y_2 - computation_y)/2 + 1: (box_max_y_2 + computation_y)/2)), ...
          'color', [0.5, 0.5, 0.5], 'linewidth', 1.5)
% up
plot(ax3, lon_series_original(lon_indices_plot((box_max_x_2 - computation_x)/2 + 1: (box_max_x_2 + computation_x)/2)), ...
          ones(computation_y, 1) * lat_series_original(lat_indices_plot((box_max_y_2 + computation_y)/2)), ...
          'color', [0.5, 0.5, 0.5], 'linewidth', 1.5)
% down
plot(ax3, lon_series_original(lon_indices_plot((box_max_x_2 - computation_x)/2 + 1: (box_max_x_2 + computation_x)/2)), ...
          ones(computation_y, 1) * lat_series_original(lat_indices_plot((box_max_y_2 - computation_y)/2 + 1)), ...
          'color', [0.5, 0.5, 0.5], 'linewidth', 1.5)

% plot coastlines
load coastlines;
plot(ax3, [coastlon; NaN; coastlon + 360; NaN; coastlon + 720 ], [coastlat; NaN; coastlat; NaN; coastlat], 'Color', 'black');

ax4 = axes(main_fig.handle);
hold on;
main_fig.axes{fig_ind} = [main_fig.axes{fig_ind}, ax4];

% divergence of the Q-vector
levels_3 = [-6, -4, -2, 2, 4, 6]; % in 1e-17, unit of s^-3*Pa^-1
if C_TERM
    temp = -A_s + B_s;
    [c2, h2] = contour(ax4, X_inversion, Y_inversion, ...
        temp(((box_max_x_2 - computation_x) / 2 + 1) : (box_max_x_2 + computation_x) / 2, ...
             ((box_max_y_2 - computation_y) / 2 + 1) : (box_max_y_2 + computation_y) / 2, ...
             level_s == plot_level, 2) * 1e17, ...
        levels_3, 'linewidth', 1.0);
else
    [c2, h2] = contour(ax4, X_inversion, Y_inversion, ...
        A_s(((box_max_x_2 - computation_x) / 2 + 1) : (box_max_x_2 + computation_x) / 2, ...
            ((box_max_y_2 - computation_y) / 2 + 1) : (box_max_y_2 + computation_y) / 2, ...
            level_s == plot_level, 2) * 1e17, ...
        levels_3, 'linewidth', 1.0);
end
clabel(c2, h2, levels_3, 'interpreter', 'latex', 'LabelSpacing', 100, 'FontSize', 8);

colors_3 = colormap(ax4, cbrewer('div', 'PRGn', 26, 'pchip'));
colormap(ax4, flipud(colors_3))

ax4.Visible = 'off';
ax4.XTick = [];
ax4.YTick = [];

axis(ax3, [lonmin_plot - 0.2, lonmax_plot + 0.2, latmin_plot - 0.2, latmax_plot + 0.2]);
axis(ax4, [lonmin_plot - 0.2, lonmax_plot + 0.2, latmin_plot - 0.2, latmax_plot + 0.2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set location of the axes

% parameters of the plot
left_margin     = 0.08;
right_margin    = 0.12;
down_margin     = 0.05;
up_margin       = 0.05;
margins = [left_margin, right_margin, down_margin, up_margin];
distance_x = 105; % in unit of pixels
distance_y = 70;

positions = set_axes(main_fig, figsize_multipanel, margins, distance_x, distance_y);

% minor adjustments
positions{1} = positions{1} + [0, 0, 0.035, 0];
set(main_fig.axes{1}, 'Position', positions{1});
positions{2} = positions{2} + [0.005, -0.03, -0.01, 0.03];
set(main_fig.axes{2}, 'Position', positions{2});
main_fig.axes{2}(2).Colorbar.Position = main_fig.axes{2}(2).Colorbar.Position + [-0.027, 0, 0, 0];

% print the figure
print(main_fig.handle, [plot_path, 'synoptic_map'], '-depsc', '-r0', '-painters')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


