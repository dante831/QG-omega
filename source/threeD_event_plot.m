
function main_fig = threeD_event_plot(main_fig, t, quantile_p, ... 
							lat, lon, level, omega_QG, event_precip, plot_level, event_timespan, ...
                            omega_x, omega_y, PLOT_OMEGA)

% reverse the omega_QG field
omega_QG = flip(omega_QG, 3);
level2 = 1e5 - level; % level2 is the absolute z position on the plot

%[Lat, Lon] = meshgrid(lat, lon);
[Lon, Lat] = meshgrid(lon, lat);
[Lat_1, Lon_1] = meshgrid(lat, lon);
[Level_y, Lat_z] = meshgrid(level2, lat);
[Level_x, Lon_z] = meshgrid(level2, lon);

ind_x = (length(lon) + 1) / 2 + 4;
ind_y = (length(lat) + 1) / 2 + 4;
lower_level = 50000; % the lower limit of the plot

plot_level = 50000; % vertical z location
plot_level_end = 0;

omega_QG_xy = omega_QG(:, :, length(level) - find(level == plot_level) + 1, 2);
omega_QG_yz = squeeze(omega_QG(:, ind_x, length(level):-1:1, 2));
omega_QG_xz = squeeze(omega_QG(ind_y, :, length(level):-1:1, 2));
omega_QG_xy_top = omega_QG(:, :, 1, 2);
omega_QG_xy_top(1:ind_y-1, 1:ind_x-1) = NaN;
omega_QG_yz_top = squeeze(omega_QG(1, :, length(level):-1:1, 2));
omega_QG_yz_top(1:ind_y-1, level2>(1e5 - plot_level)) = NaN; % level2 > (1e5 - plot_level) because the entire graph is reversed
omega_QG_xz_top = squeeze(omega_QG(:, 1, length(level):-1:1, 2));
omega_QG_xz_top(1:ind_x-1, level2>(1e5 - plot_level)) = NaN;

% see more of the precip field
omega_QG_xz_top(:, level2<(1e5 - lower_level)) = NaN;
omega_QG_xz    (:, level2<(1e5 - lower_level)) = NaN;
omega_QG_yz_top(:, level2<(1e5 - lower_level)) = NaN;
omega_QG_yz    (:, level2<(1e5 - lower_level)) = NaN;


event_day = mod(event_timespan(2), 365);
day_series = event_day + ([-30+t:30+t])*0.25;
event_year = 1990 + floor(event_timespan(2)/365);

colors = colormap(jet(140));

clevels = 80;
fig = main_fig.handle;
position_1 = [0.08, 0.15, 0.35, 0.75];
position_2 = [0.55, 0.12, 0.30, 0.80];


% Left figure: 1D plot precip at the center
ax1 = axes(fig);
main_fig.axes{1} = [ax1];
temp = squeeze(event_precip((1+end)/2, (1+end)/2, (day_series + (event_year - 1990)*365) / 0.25));
plot(day_series, temp, '-s', 'Linewidth', 0.8, 'markersize', 1.5)
hold on;
scatter(event_day + t*0.25, temp((end+1)/2), 'o', 'MarkerFaceColor', 'red')
plot(day_series, ones(size(day_series))*quantile_p, 'Linewidth', 0.8, 'Linestyle', '--')
text(event_day + t*0.25 + 2, quantile_p - 4, '99.9th percentile', 'interpreter', 'latex')
title(ax1, '(a) Precipitation time series', 'interpreter', 'latex', 'FontSize', 11)
xlabel(['days since 01/01/', num2str(event_year)], 'interpreter', 'latex', 'FontSize', 11)
ylabel('mm day$^{-1}$', 'interpreter', 'latex')
set(ax1, 'TickLabelInterpreter', 'latex')
axis([day_series(1), day_series(end), 0.0, 140.0])
set([ax1], 'Position', position_1);
set(ax1, 'TickDir', 'out');
% remove the upper and right ticks on the box
box off

% Right figure: plot 2D precip around the event
ax2 = axes(fig);
main_fig.axes{2} = [ax2];
hold on;

expand_factor = 1;
levels = [10, 30, 50, 70, 90];
[C, hContour] = contour3(ax2, Lon_1, Lat_1, event_precip(:, :, event_timespan(2)/0.25 + t)*expand_factor, levels, 'ShowText','on');
clabel(C, hContour, [10, 50, 90], 'interpreter', 'latex', 'LabelSpacing', 100, 'FontSize', 8);
scatter3(ax2, lon((end+1)/2), lat((end+1)/2), ...
    event_precip((end+1)/2, (end+1)/2, event_day/0.25 + t)*expand_factor, 25, 'o', ...
    'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'none')
caxis(ax2, [0, 100]*expand_factor)

if ~PLOT_OMEGA
    cb = colorbar;
end
title(ax2, '(b) $\omega_{QG}$ and precipitation', 'interpreter', 'latex', 'FontSize', 11)
xlabel('lon', 'interpreter', 'latex')
xticks([170, 175, 180])
xticklabels({'$170^\circ$E', '$175^\circ$E', '$180^\circ$'})
yticks([40, 45, 50])
yticklabels({'$40^\circ$N', '$45^\circ$N', '$50^\circ$N'})
ylabel('lat', 'interpreter', 'latex')
zlabel('hPa', 'interpreter', 'latex')
zticklabs = [0, 25000, 50000, 75000, 100000];
zticks(zticklabs);
labels = {};
for i = 1 : length(zticklabs)
    labels{i} = num2str(1e3 - zticklabs(i)/100);
end
zticklabels(labels)
set(ax2, 'TickLabelInterpreter', 'latex')
colormap('parula')
axis([lon(1) - 0.1, lon(end), lat(1) - 0.1, lat(end), 0, 100000])
view(-40, 30)
grid on

% plot omega field
if PLOT_OMEGA

    ax3 = axes(fig);
    main_fig.axes{2} = [main_fig.axes{2}, ax3];
    set(ax3, 'color', 'none');
    set(ax3, 'xtick', [])
    set(ax3, 'ytick', [])
    set(ax3, 'ztick', [])

    hold on;
    s1 = surf(ax3, Lon, Lat, ones(size(Lon))*(1e5 - plot_level), clevels, 'Linestyle', 'none');
    s1.CData = omega_QG_xy;
    s1.FaceColor = 'interp';
    % rainfall center
    plot3([lon((end+1)/2), lon((end+1)/2)], [lat((end+1)/2), lat((end+1)/2)], [1e5 - plot_level, plot_level_end], ...
            'k', 'markersize', 5, 'MarkerFaceColor', 'red', 'markeredgecolor', 'none', 'linestyle', 'none')
    % omega_QG center
    scatter3(lon((end+1)/2), lat((end+1)/2), 1e5 - plot_level, 25, 'o', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'none')
    scatter3(omega_x, omega_y, 1e5 - plot_level, 70, 'p', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'none')

    s2 = surf(ax3, ones(size(Lat_z))*lon(ind_x), Lat_z, Level_y, clevels, 'Linestyle', 'none');
    s2.CData = omega_QG_yz;
    s2.EdgeColor = 'none';
    s2.FaceColor = 'interp';
    s3 = surf(ax3, Lon_z, ones(size(Lon_z))*lat(ind_y), Level_x, clevels, 'Linestyle', 'none');
    s3.CData = omega_QG_xz;
    s3.EdgeColor = 'none';
    s3.FaceColor = 'interp';
    s4 = surf(ax3, Lon, Lat, ones(size(Lon))*level2(end), clevels, 'Linestyle', 'none');
    s4.CData = omega_QG_xy_top;
    %s4.EdgeColor = [0, 0, 0];
    %s4.LineStyle = '-'
    s4.FaceColor = 'interp';
    s5 = surf(ax3, ones(size(Lat_z))*lon(1), Lat_z, Level_y, clevels, 'Linestyle', 'none');
    s5.CData = omega_QG_yz_top;
    %s5.EdgeColor = 'none';
    s5.FaceColor = 'interp';
    s6 = surf(ax3, Lon_z, ones(size(Lon_z))*lat(1), Level_x, clevels, 'Linestyle', 'none');
    s6.CData = omega_QG_xz_top;
    %s6.EdgeColor = 'none';
    s6.FaceColor = 'interp';

    d = 0.1;
    line_width = 0.8;
    plot3(ax3, [lon(1),     lon(ind_x)], [lat(1)    , lat(1)] - d, [1e5 - plot_level, 1e5 - plot_level], 'k', 'LineWidth', line_width)
    plot3(ax3, [lon(ind_x), lon(end)]  , [lat(1)    , lat(1)] - d, [level2(end), level2(end)] + 200, 'k', 'LineWidth', line_width)
    plot3(ax3, [lon(1),     lon(ind_x)], [lat(ind_y), lat(ind_y)] - d, [level2(end), level2(end)] + 200, 'k', 'LineWidth', line_width)
    plot3(ax3, [lon(1),     lon(end)]  , [lat(end)  , lat(end)]        , [level2(end), level2(end)] + 200, 'k', 'LineWidth', line_width)
    plot3(ax3, [lon(1),     lon(ind_x)], [lat(ind_y), lat(ind_y)] - d, [1e5 - plot_level, 1e5 - plot_level], 'k', 'LineWidth', line_width)

    plot3(ax3, [lon(1)    , lon(1)] - d    , [lat(1)  , lat(ind_y)], [1e5 - plot_level, 1e5 - plot_level], 'k', 'LineWidth', line_width)
    plot3(ax3, [lon(1)    , lon(1)] - d    , [lat(end), lat(ind_y)], [level2(end), level2(end)] + 200, 'k', 'LineWidth', line_width)
    plot3(ax3, [lon(ind_x), lon(ind_x)] - d, [lat(1)  , lat(ind_y)], [level2(end), level2(end)] + 200, 'k', 'LineWidth', line_width)
    plot3(ax3, [lon(end) , lon(end)]         , [lat(1)  , lat(end)]  , [level2(end), level2(end)] + 200, 'k', 'LineWidth', line_width)
    plot3(ax3, [lon(ind_x), lon(ind_x)] - d, [lat(1)  , lat(ind_y)], [1e5 - plot_level, 1e5 - plot_level], 'k', 'LineWidth', line_width)

    plot3(ax3, [lon(1)    , lon(1)] - d    , [lat(ind_y), lat(ind_y)] - d,[1e5 - plot_level, level2(end)], 'k', 'LineWidth', line_width)
    plot3(ax3, [lon(ind_x), lon(ind_x)] - d, [lat(1)    , lat(1)] - d    ,[1e5 - plot_level, level2(end)], 'k', 'LineWidth', line_width)
    plot3(ax3, [lon(ind_x), lon(ind_x)] - d, [lat(ind_y), lat(ind_y)] - d,[1e5 - plot_level, level2(end)], 'k', 'LineWidth', line_width)
    plot3(ax3, [lon(1)    , lon(1)] - d    , [lat(end)  , lat(end)]        ,[plot_level_end    , level2(end)], 'k', 'LineWidth', line_width)
    plot3(ax3, [lon(end)  , lon(end)]        , [lat(1)    , lat(1)] - d    ,[plot_level_end    , level2(end)], 'k', 'LineWidth', line_width)
    plot3(ax3, [lon(1)    , lon(1)] - d    , [lat(1)    , lat(1)] - d    ,[plot_level_end    , 1e5 - plot_level], 'k', 'LineWidth', line_width)

    plot3(ax3, [lon(ind_x), lon(end)]      , [lat(1)    , lat(1)]     - d, [1e5 - plot_level, 1e5 - plot_level], '--k', 'LineWidth', line_width)
    plot3(ax3, [lon(1)    , lon(1)]     - d, [lat(ind_y), lat(end)]   - d, [1e5 - plot_level, 1e5 - plot_level], '--k', 'LineWidth', line_width)
    plot3(ax3, [lon(1)    , lon(1)]     - d, [lat(ind_y), lat(ind_y)] - d, [plot_level_end, 1e5 - plot_level]  , '--k', 'LineWidth', line_width)
    plot3(ax3, [lon(ind_x), lon(ind_x)] - d, [lat(1)    , lat(1)]     - d, [plot_level_end, 1e5 - plot_level]  , '--k', 'LineWidth', line_width)

    caxis(ax3, [-1.5, 0.5])
    colormap(ax3, colors)
    cb = colorbar;
    axis([lon(1) - 0.1, lon(end), lat(1) - 0.1, lat(end), 0, 100000])
    %set(gca, 'Zdir', 'reverse')
    view(-40, 30)
    grid on
    set([ax3], 'Position', position_2);

    linkprop([ax2, ax3], 'CameraPosition');

end

set([ax2], 'Position', position_2);

% deal with colorbar
cpos = get(cb,'Position');
if PLOT_OMEGA
    cb.Label.String = 'Pa s$^{-1}$';
    cb.Label.Rotation = 0;
    cb.Label.Position = cpos(1:2) + [-0.30, -1.7];
    cb.Ticks = [-1.2:0.3:0.3];
else
    cb.Label.String = 'mm day$^{-1}$';
    cb.Ticks = [0:20:100];
end
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';
set(cb, 'FontSize', 10);


