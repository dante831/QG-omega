function plot_2D_map(figsize, temp, X, Y, plot_filename, plot_title, plot_path, smooth_window_x, smooth_window_y, SMOOTH, climit);

    figure('pos', figsize)
    box on;
    hold on;
    ax1 = gca;
    if SMOOTH
        temp = smooth2a(temp, smooth_window_x, smooth_window_y);
    end
    inc = (climit(2) - climit(1)) / 50;
    levs = climit(1) : inc : climit(2);
    [c1, h1] = contourf(ax1, X, Y, temp, levs, 'LineColor', 'none');
    cbar = colorbar;
    colormap(ax1, 'jet')
    axis([X(1, 1) - 0.2, X(1, end) + 0.2, Y(1, 1), Y(end, 1)])
    caxis(climit)
    title(plot_title, 'interpreter', 'latex')
    xlabel('lon')
    ylabel('lat')
    ax2 = axes;
    load coastlines;
    plot(ax2, [coastlon; NaN; coastlon + 360; NaN; coastlon + 720 ], [coastlat; NaN; coastlat; NaN; coastlat], 'Color', 'black')
    axis([X(1, 1), X(1, end), Y(1, 1), Y(end, 1)])
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    set([ax1, ax2], 'Position', [.05 .11 .88 .815]);
    %cb1 = colorbar(ax1,'Position',[.09 .11 .0275 .815]);
    %cb1.Label.String = '';
    linkaxes([ax1, ax2], 'xy');
    saveas(gca, [plot_path, plot_filename], 'png')
    clf;

