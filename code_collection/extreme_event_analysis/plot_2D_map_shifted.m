function plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event, num_threshold, ZERO, FLIP);

    figure('pos', figsize)
    hold on;
    
    ax1 = gca;
    if SMOOTH
        temp = smooth2a(temp, smooth_window_x, smooth_window_y);
    end
    if exist('num_event') && exist('num_threshold')
        temp(num_event < num_threshold) = NaN; % mask out grid points whose events are fewer than num_threshold
    end
    n_level = 20;
    inc = (climit(2) - climit(1)) / n_level;
    levs = climit(1) : inc : climit(2);
    half_range = max(max(temp(:)) - sum(climit)/2, - min(temp(:)) + sum(climit)/2);
    levs(1)     = min(levs(1),   sum(climit)/2 - half_range);
    levs(end)   = max(levs(end), sum(climit)/2 + half_range);
    
    %[c1, h1] = contourf(ax1, [X; X - 360], [Y; Y], [temp; temp], levs, 'LineColor', 'none');
    [c1, h1] = contourf(ax1, [X - 360, X], [Y, Y], [temp, temp], levs, 'LineColor', 'none');

    set(ax1, 'xtick', [-120:60:120])
    set(ax1, 'ytick', [-60:30:60])
    set(ax1, 'Xticklabel', {'120W', '60W', '0', '60E', '120E'})
    set(ax1, 'Yticklabel', {'60S', '30S', '0', '30N', '60N'})
    set(ax1,'TickLabelInterpreter','latex')

    if exist('FLIP')&&FLIP
        colors = flipud(colormap(jet(28))); % 23 and 28 are symmetrical
    else
        colors = colormap(jet(28));
    end
    if ZERO
        %colors = colors([1 : floor(end * 1 / 2.5), ceil(end * (1 - 1 / 2.5)) : end], :);
        colors = colors([1:10, 19:28], :);
        colors(n_level / 2 : n_level / 2 + 1, :) = 0.99;
        %colors(floor(end / 2) - 1 : ceil(end / 2 + 1) + 1, :) = 0.99;
        %colors(1 : floor(end / 2) - 1) = interp()
    end
    colormap(ax1, colors);
    %axis([X(1, 1) - 0.2, X(1, end) + 0.2, Y(1, 1), Y(end, 1)])
    axis([X(1, 1) - 180, X(1, end) - 180, Y(1, 1), Y(end, 1)])
    caxis(climit)
    title(plot_title, 'interpreter', 'latex')
    
    %% deal with the box
    box on;
    ax1_pos = ax1.Position; % position of first axes
    ax1_1 = axes('Position',ax1_pos,...
            'XAxisLocation','bottom',...
            'YAxisLocation','left',...
            'Color','none');
    set(ax1_1, 'xtick', [])
    set(ax1_1, 'ytick', [])
    ax1_2 = axes('Position',ax1_pos,...
            'XAxisLocation','top',...
            'YAxisLocation','right',...
            'Color','none');
    set(ax1_2, 'xtick', [])
    set(ax1_2, 'ytick', [])


    %ax1.LineWidth = 1.5;
    %xlabel('lon')
    %ylabel('lat')
    ax2 = axes;
   
    % plot mask
    if exist('num_event') && exist('num_threshold')
        temp2 = temp;
        s = size(temp2);
        [ind_y, ind_x] = find(isnan(temp));
        temp2(~isnan(temp)) = NaN;
        aaa = 0.1;
        temp2(sub2ind(s, ind_y, ind_x)) = aaa;
        % map one point further around the NaN points
        temp2(sub2ind(s, min(ind_y + 1, size(temp2, 1)), ind_x)) = aaa; 
        temp2(sub2ind(s, max(ind_y - 1, 1), ind_x)) = aaa;
        temp2(sub2ind(s, ind_y, mod(ind_x + 1 - 1, size(temp2, 2)) + 1)) = aaa; 
        temp2(sub2ind(s, ind_y, mod(ind_x - 1 - 1, size(temp2, 2)) + 1)) = aaa;
        % one choice:
        [c2, h2] = contourf(ax2, [X - 360, X], [Y, Y], [temp2, temp2], [0:0.1:5], 'LineColor', 'none');
        colormap(ax2, 'bone')
        % another choice:
        %ch = pcolor(ax2, [X - 360, X] - 0.5, [Y, Y] - 0.5, [temp2, temp2])
        %patches = findobj(ch,'-property','FaceAlpha');
        %for ph = patches
        %    set(ph, 'FaceAlpha', 0.8, 'EdgeAlpha',0)
        %end
    end
    axis([X(1, 1) - 180, X(1, end) - 180, Y(1, 1), Y(end, 1)])
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];



    ax3 = axes;
    load coastlines;
    plot(ax3, [coastlon; NaN; coastlon + 360; NaN; coastlon + 720 ], [coastlat; NaN; coastlat; NaN; coastlat], 'Color', 'black')
     
    axis([X(1, 1) - 180, X(1, end) - 180, Y(1, 1), Y(end, 1)])
    ax3.Visible = 'off';
    ax3.XTick = [];
    ax3.YTick = [];
     
    set([ax1, ax1_1, ax1_2, ax2, ax3], 'Position', [.07 .12 .82 .76]);
    cb = colorbar(ax1,'Position',[.92 .12 .0275 .76]);
    cpos = get(cb,'Position');
    cb.Label.Position = [cpos(1) - 0.33, cpos(2) - 0.83];
    cb.Label.Rotation = 0;
    cb.TickLabelInterpreter = 'latex';
    %cb.Ticks = [-0.12:0.06:0.12];
    set(cb, 'FontSize', 10);

    linkaxes([ax1, ax1_1, ax1_2, ax2, ax3], 'xy');
    set(ax1, 'TickDir', 'out');
    %savefig([plot_path, plot_filename])
    %saveas(gca, [plot_path, plot_filename], 'pdf')
    saveas(gca, [plot_path, plot_filename], 'png')
    clf;

