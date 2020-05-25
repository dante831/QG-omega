function main_fig = subplot_2D_map(main_fig, cmp, margins, fig_index, figsize, ...
        temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event, num_threshold, ...
        cb_pos, cb_string, cb_label_pos_bias, cb_ticks, cb_ticklabels, ...
        lat_mask, ZERO, FLIP, LIGHT_COLOR, varargin);
  
    % suppress warning for constant contour when plotting the grey masking
    MSGID = 'MATLAB:contour:ConstantData';
    warning('off', MSGID)

    while ~isempty(varargin)
        switch lower(varargin{1})  
            case 'colorbar'
                COLORBAR = varargin{2};
            otherwise
                error(['Unexpected option: ' varargin{1}])
        end
        varargin(1:2) = [];  
    end
    if exist('COLORBAR') ~= 1
        COLORBAR = false;
    end

    ax1 = axes(main_fig.handle);
    hold on;
    
    if SMOOTH
        pwd_str = pwd;
        if ~isempty(strfind(pwd_str, 'CESM')) 
            N_s = 3; % have CESM be the same as GFDL, according to Paul
        else
            N_s = 3;
        end
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
    
    [c1, h1] = contourf(ax1, [X - 360, X], [Y, Y], [temp, temp], levs, 'LineColor', 'none');
    
    if FLIP
        % 23 and 28 are symmetrical, flip to make it consistent with Pfahl et al, 2017
        colors = flipud(colormap(ax1, cmp(28)));
    else
        colors = colormap(ax1, cmp(28)); 
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
    colormap(ax1, colors);
    %axis([X(1, 1) - 0.2, X(1, end) + 0.2, Y(1, 1), Y(end, 1)])
    axis([X(1, 1) - 180, X(1, end) - 180, Y(1, 1), Y(end, 1)])
    caxis(climit)
    title(ax1, plot_title, 'interpreter', 'latex')
    
    set(ax1, 'xtick', [-120:60:120])
    if ~(ceil(fig_index / main_fig.col) == main_fig.row)
        set(ax1, 'Xticklabel', [])
    else
        set(ax1, 'Xticklabel', {'120W', '60W', '0', '60E', '120E'})
    end
    set(ax1, 'ytick', [-60:30:60])
    if ~(mod(fig_index - 1, main_fig.col) == 0)
        set(ax1, 'Yticklabel', [])
    else
        set(ax1, 'Yticklabel', {'60S', '30S', '0', '30N', '60N'})
    end
    set(ax1, 'TickLabelInterpreter', 'latex')
    set(ax1, 'linewidth', 0.8)


    %% deal with the box
    box on;
    ax1_pos = ax1.Position; % position of first axes
    ax1_1 = axes(main_fig.handle, 'Position',ax1_pos,...
            'XAxisLocation','bottom',...
            'YAxisLocation','left',...
            'Color','none');
    set(ax1_1, 'xtick', [])
    set(ax1_1, 'ytick', [])
    ax1_2 = axes(main_fig.handle, 'Position',ax1_pos,...
            'XAxisLocation','top',...
            'YAxisLocation','right',...
            'Color','none');
    set(ax1_2, 'xtick', [])
    set(ax1_2, 'ytick', [])

    % plot mask
    ax2 = axes(main_fig.handle, 'Position', ax1_pos);
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
        [c2, h2] = contourf(ax2, [X - 360, X], [Y, Y], [temp2, temp2], [0:0.1:5], 'LineColor', 'none');
        colormap(ax2, 'bone')
    end
    axis([X(1, 1) - 180, X(1, end) - 180, Y(1, 1), Y(end, 1)])
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];

    
    ax3 = axes(main_fig.handle, 'Position', ax1_pos);
    load coastlines;
    plot(ax3, [coastlon; NaN; coastlon + 360; NaN; coastlon + 720 ], [coastlat; NaN; coastlat; NaN; coastlat], 'Color', 'black')
     
    axis([X(1, 1) - 180, X(1, end) - 180, Y(1, 1), Y(end, 1)])
    ax3.Visible = 'off';
    ax3.XTick = [];
    ax3.YTick = [];
       
    set(ax1,'TickDir','out');
    if fig_index == 1 || COLORBAR
        cb = colorbar(ax1, 'Position', cb_pos);
        cb.Label.Interpreter = 'latex';
        cb.Label.String = cb_string;
        cpos = get(cb,'Position');
        cb.Label.Position(1:2) = cpos(1:2) + cb_label_pos_bias;
        cb.Label.Rotation = 0;
        cb.Ticks = cb_ticks;
        if ~isempty(cb_ticklabels)
            cb.TickLabels = cb_ticklabels;
        end
        cb.TickLabelInterpreter = 'latex';
        set(cb, 'FontSize', 10);
    else
        colorbar(ax1, 'off')
    end
    
    main_fig.axes{fig_index} = [ax1, ax1_1, ax1_2, ax2, ax3];
  
