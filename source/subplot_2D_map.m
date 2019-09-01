function subplot_2D_map(cmp, margins, ax, row, col, fig_index, figsize, ...
        temp, X, Y, plot_title, SMOOTH_2, climit, num_event, num_threshold, ...
        cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
        lat_mask, ZERO, FLIP, LIGHT_COLOR);


    distance = 30;
    x_distance  = distance / figsize(3);
    y_distance  = distance / figsize(4);
        
    left_margin     = margins(1);
    right_margin    = margins(2);
    down_margin     = margins(3);
    up_margin       = margins(4);


    width_fig   = (1 - left_margin - right_margin - x_distance * (col - 1)) / col;
    height_fig  = (1 - down_margin - up_margin    - y_distance * (row - 1)) / row;

    position = [left_margin + (mod(fig_index - 1, col)) * (width_fig + x_distance), ...
                1 - (up_margin + ceil(fig_index / col) * (height_fig + y_distance)) + y_distance, ...
                width_fig, height_fig];

    ax1{fig_index} = subplot(row, col, fig_index);
    %ax1{fig_index} = axes;
    
    if SMOOTH_2
        N_s = 3; % have CESM be the same as GFDL, according to Paul
        for i = 1 : N_s
            temp = one_two_one_2D(temp, 2);
        end
    end
    if exist('num_event') && exist('num_threshold')
        temp(num_event < num_threshold) = NaN; % mask out grid points whose events are fewer than num_threshold
    end
    
    % mask out the tropics
    temp(Y < lat_mask & Y > -lat_mask) = NaN;

    % deal with colormaps
    n_level = 20;  
    inc = (climit(2) - climit(1)) / n_level;
    levs = climit(1) : inc : climit(2);
    half_range = max(max(temp(:)) - sum(climit)/2, - min(temp(:)) + sum(climit)/2);
    % let the maximum and minimum level incorporate the whole range of plotted data
    levs(1)     = min(levs(1),   sum(climit)/2 - half_range); 
    levs(end)   = max(levs(end), sum(climit)/2 + half_range);
    
    [c1, h1] = contourf(ax1{fig_index}, [X - 360, X], [Y, Y], [temp, temp], levs, 'LineColor', 'none');
    
    if FLIP
        % 23 and 28 are symmetrical, flip to make it consistent with Pfahl et al, 2017
        colors = flipud(colormap(ax1{fig_index}, cmp(28)));
    else
        colors = colormap(ax1{fig_index}, cmp(28)); 
    end

    if ZERO
        % get two opposite colors from a specified symmetric colormap
        if LIGHT_COLOR
            colors = colors([3:12, 17:26], :);
        else
            colors = colors([1:10, 19:28], :);
        end
        colors(n_level / 2 : n_level / 2 + 1, :) = 0.99;
    end
    colormap(ax1{fig_index}, colors);
    %axis([X(1, 1) - 0.2, X(1, end) + 0.2, Y(1, 1), Y(end, 1)])
    axis([X(1, 1) - 180, X(1, end) - 180, Y(1, 1), Y(end, 1)])
    caxis(climit)
    title(ax1{fig_index}, plot_title, 'interpreter', 'latex')
    
    set(ax1{fig_index}, 'xtick', [-120:60:120])
    if ~(ceil(fig_index / col) == row)
        set(ax1{fig_index}, 'Xticklabel', [])
    else
        set(ax1{fig_index}, 'Xticklabel', {'120W', '60W', '0', '60E', '120E'})
    end
    set(ax1{fig_index}, 'ytick', [-60:30:60])
    if ~(mod(fig_index - 1, col) == 0)
        set(ax1{fig_index}, 'Yticklabel', [])
    else
        set(ax1{fig_index}, 'Yticklabel', {'60S', '30S', '0', '30N', '60N'})
    end
    set(ax1{fig_index}, 'TickLabelInterpreter', 'latex')
    set(ax1{fig_index}, 'linewidth', 0.8)


    %% deal with the box
    box on;
    ax1_pos = ax1{fig_index}.Position; % position of first axes
    ax1_1{fig_index} = axes('Position',ax1_pos,...
            'XAxisLocation','bottom',...
            'YAxisLocation','left',...
            'Color','none');
    set(ax1_1{fig_index}, 'xtick', [])
    set(ax1_1{fig_index}, 'ytick', [])
    ax1_2{fig_index} = axes('Position',ax1_pos,...
            'XAxisLocation','top',...
            'YAxisLocation','right',...
            'Color','none');
    set(ax1_2{fig_index}, 'xtick', [])
    set(ax1_2{fig_index}, 'ytick', [])

    % plot mask
    ax2{fig_index} = axes('Position', ax1_pos);
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
        [c2, h2] = contourf(ax2{fig_index}, [X - 360, X], [Y, Y], [temp2, temp2], [0:0.1:5], 'LineColor', 'none');
        colormap(ax2{fig_index}, 'bone')
    end
    axis([X(1, 1) - 180, X(1, end) - 180, Y(1, 1), Y(end, 1)])
    ax2{fig_index}.Visible = 'off';
    ax2{fig_index}.XTick = [];
    ax2{fig_index}.YTick = [];

    
    ax3{fig_index} = axes('Position', ax1_pos);
    load coastlines;
    plot(ax3{fig_index}, [coastlon; NaN; coastlon + 360; NaN; coastlon + 720 ], [coastlat; NaN; coastlat; NaN; coastlat], 'Color', 'black')
     
    axis([X(1, 1) - 180, X(1, end) - 180, Y(1, 1), Y(end, 1)])
    ax3{fig_index}.Visible = 'off';
    ax3{fig_index}.XTick = [];
    ax3{fig_index}.YTick = [];
     
    set([ax1{fig_index}, ax1_1{fig_index}, ax1_2{fig_index}, ax2{fig_index}, ax3{fig_index}], 'Position', position);
    set(ax1{fig_index},'TickDir','out');
    if fig_index == 1
        cb = colorbar(ax1{fig_index}, 'Position', cb_pos);
        cb.Label.String = cb_string;
        cpos = get(cb,'Position');
        cb.Label.Position = cpos(1:2) + cb_label_pos_bias;
        cb.Label.Rotation = 0;
        cb.Label.Interpreter = 'latex';
        cb.Ticks = cb_ticks;
        if ~isempty(cb_ticklabels)
            cb.TickLabels = cb_ticklabels;
        end
        cb.TickLabelInterpreter = 'latex';
        set(cb, 'FontSize', 10);
    else
        colorbar(ax1{fig_index}, 'off')
    end
    linkaxes([ax1{fig_index}, ax1_1{fig_index}, ax1_2{fig_index}, ax2{fig_index}, ax3{fig_index}], 'xy');
  
