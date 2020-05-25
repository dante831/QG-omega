function event_distribution(climate_string, events_sta, lat_series, lon_series, lat_indices, lon_indices, plot_path)

%ps_filename = '/net/chickpea/volume1/scratch/ziweili/test1_GFDL/ps/ps_avg_historical.nc';
ps_filename = ['/net/aimsir/archive1/ziweili/CESM_LENS/output/ps_avg_historical.nc'];
ps_avg = ncread(ps_filename, 'ps_avg');
lon = lon_series(lon_indices);
p_max = 75000;

% create a 5-by-5 box around each point in ind to make sure that it shows excluded events
% instead of just places over p_max
ind = ps_avg < p_max;
[x_, y_] = find(ind);
box_min = 9;
xx = (box_min - 1) / 2;
for k = 1 : length(x_)
    ind(mod([-xx : xx] + x_(k) - 1, length(ind(:, 1))) + 1, mod([-xx : xx] + y_(k) - 1, length(ind(1, :))) + 1) = true;
end

%if isequal(climate, 'historical')
%    events_sta = events_sta_historical;
%elseif isequal(climate, 'rcp85')
%    events_sta = events_sta_rcp85;
%end
events_map = zeros(length(lat_series), length(lon_series));
for k = 1 : length(events_sta.events_sta)
    lons_ = [events_sta.events_sta{k}.event_lon];
    lats_ = [events_sta.events_sta{k}.event_lat];
    j = lat_indices(k);
    for l = 1 : length(lons_)
        i = find(lon == lons_(l));
        events_map(j, i) = events_map(j, i) + 1;
    end
end

figure('pos',[10 10 600 350])
ax1 = axes;
box on;
hold on;
[X, Y] = meshgrid(lon_series, lat_series);
%temp = sum(events_map, 2) + double(sum(events_map, 2) == 0);
%events_map = events_map ./ repmat(temp, 1, length(lon_series));
[c1, h1] = contourf(X, Y, events_map, 20, 'Linestyle', 'none');
cbar = colorbar;
cbar.Label.String = 'count';
colormap(ax1, 'bone');
colormap(ax1, flipud(colormap(ax1)));
ax2 = axes;
[X2, Y2] = meshgrid(lon_series, lat_series);
%[c2, h2] = contour(X2, Y2, orog', 20);
[c2, h2] = contour(X2, Y2, ps_avg', ':');
%[c2, h2] = contour(X2, Y2, ind', 20);
cbar = colorbar;
cbar.Label.String = 'pressure';
colormap(ax2, 'jet')
colormap(ax2, flipud(colormap(ax2)));
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
load coastlines;
ax3 = axes;
plot([coastlon; NaN; coastlon + 360; NaN; coastlon + 720 ], [coastlat; NaN; coastlat; NaN; coastlat], 'Color', 'black')
ax3.Visible = 'off';
ax3.XTick = [];
ax3.YTick = [];
xlim([0 360]);
ylim([-90 90]);
ax4 = axes;
[c4, h4] = contour(X2, Y2, ind', [1, 1], 'r');
ax4.Visible = 'off';
ax4.XTick = [];
ax4.YTick = [];
if isequal(climate_string, 'historical')
    title(ax1, 'Historical Event Distribution');
elseif isequal(climate_string, 'rcp85')
    title(ax1, 'Rcp85 Event Distribution');
end
%xlabel(ax1, 'Longitude');
%ylabel(ax1, 'Latitude');
set([ax1, ax2, ax3, ax4],'Position',[.17 .11 .685 .815]);
cb1 = colorbar(ax1,'Position',[.09 .11 .0275 .815]);
cb1.Label.String = 'Number of events';
cb2 = colorbar(ax2,'Position',[.88 .11 .0275 .815]);
cb2.Label.String = 'Surface pressure';
linkaxes([ax1, ax2, ax3, ax4], 'xy');
saveas(gca, [plot_path, climate_string, '_event_distribution'], 'png');
hold off;




