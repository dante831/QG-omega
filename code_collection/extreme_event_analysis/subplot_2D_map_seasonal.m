function subplot_2D_map_seasonal(ax, row, col, fig_index, figsize, ...
        temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event, num_threshold, ZERO, tags);

    distance = 30;
    x_distance  = distance / figsize(3);
    y_distance  = distance / figsize(4);
        
    left_margin     = 0.05;
    right_margin    = 0.09;
    down_margin     = 0.10;
    up_margin       = 0.10;

    y_min = 0;
    y_max = 80;

    width_fig   = (1 - left_margin - right_margin - x_distance * (col - 1)) / col;
    height_fig  = (1 - down_margin - up_margin    - y_distance * (row - 1)) / row;

    position = [left_margin + (mod(fig_index - 1, col)) * (width_fig + x_distance), ...
                1 - (up_margin + ceil(fig_index / col) * (height_fig + y_distance)) + y_distance, ...
                width_fig, height_fig];

    ax1 = subplot(row, col, fig_index);

    if SMOOTH
        %temp = smooth2a(temp, smooth_window_x, smooth_window_y);
        pwd_str = pwd;
        if ~isempty(strfind(pwd_str, 'CESM'))
            N_s = 10;
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

    % deal with colormaps
    n_level = 20;  
    inc = (climit(2) - climit(1)) / n_level;
    levs = climit(1) : inc : climit(2);
    half_range = max(max(temp(:)) - sum(climit)/2, - min(temp(:)) + sum(climit)/2);
    % let the maximum and minimum level incorporate the whoe range of plotted data
    levs(1)     = min(levs(1),   sum(climit)/2 - half_range); 
    levs(end)   = max(levs(end), sum(climit)/2 + half_range);
    
    %[c1, h1] = contourf(ax1, [X; X - 360], [Y; Y], [temp; temp], levs, 'LineColor', 'none');
    [c1, h1] = contourf(ax1, [X - 360, X], [Y, Y], [temp, temp], levs, 'LineColor', 'none');
    

    colors = flipud(colormap(jet(28))); % 23 and 28 are symmetrical, flip to make it consistent with Pfahl et al, 2017

    if ZERO
        %colors = colors([1 : floor(end / 2.5), (end + 1) / 2, end - ceil(end / 2.5) + 1 : end], :); 
        colors = colors([1:10, 19:28], :);
                % get the blue-ish and red-ish color from jet colormap
        colors(n_level / 2 : n_level / 2 + 1, :) = 0.99;
    end
    colormap(ax1, colors);
    axis([X(1, 1) - 180, X(1, end) - 180, y_min, y_max])
    caxis(climit)
    title(ax1, plot_title, 'interpreter', 'latex', 'FontWeight', 'normal')
    
    set(ax1, 'xtick', [-120:60:120])
    if ~(ceil(fig_index / col) == row)
        set(ax1, 'Xticklabel', [])
    else
        set(ax1, 'Xticklabel', {'$120^\circ$W', '$60^\circ$W', '$0^\circ$', '$60^\circ$E', '$120^\circ$E'})
    end
    set(ax1, 'ytick', [10, 40, 70])
    if ~(mod(fig_index - 1, col) == 0)
        set(ax1, 'Yticklabel', [])
    else
        set(ax1, 'Yticklabel', {'$10^\circ$N', '$40^\circ$N', '$70^\circ$N'})
    end
    set(ax1, 'TickLabelInterpreter', 'latex')
    set(ax1, 'linewidth', 0.8)


    %% deal with the box
    box on;
    ax1_pos = ax1.Position; % position of first axes
    ax1_1 = axes('Position',ax1_pos,...
            'XAxisLocation','bottom',...
            'YAxisLocation','left',...
            'Color','none');
    axis([X(1, 1) - 180, X(1, end) - 180, y_min, y_max])
    set(ax1_1, 'xtick', [])
    set(ax1_1, 'ytick', [])
    ax1_2 = axes('Position',ax1_pos,...
            'XAxisLocation','top',...
            'YAxisLocation','right',...
            'Color','none');
    axis([X(1, 1) - 180, X(1, end) - 180, y_min, y_max])
    set(ax1_2, 'xtick', [])
    set(ax1_2, 'ytick', [])


    
    %ax1.LineWidth = 1.5;
    %xlabel('lon')
    %ylabel('lat')
    ax2 = axes('Position', ax1_pos);
   
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
    axis([X(1, 1) - 180, X(1, end) - 180, y_min, y_max])
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];

    
    ax3 = axes('Position', ax1_pos);
    load coastlines;
    plot(ax3, [coastlon; NaN; coastlon + 360; NaN; coastlon + 720 ], [coastlat; NaN; coastlat; NaN; coastlat], 'Color', 'black')
     
    axis([X(1, 1) - 180, X(1, end) - 180, y_min, y_max])
    ax3.Visible = 'off';
    ax3.XTick = [];
    ax3.YTick = [];



    %% plot land masses that are selected for seasonal analysis
    ax4_1 = axes('Position', ax1_pos);
    ax4_2 = axes('Position', ax1_pos);
    
    % Paul's idea: use two straght lines instead of contours to indicate land masses
    line_colors = get(gca,'colororder');
    h4_1 = plot(ax4_1, [-10, 65, 65, 105, 105, 120, 145, -10, -10], ...
                       [ 35, 35, 45,  45,  35,  35,  60,  60,  35], 'color', line_colors(7, :), 'LineWidth', 1.5);
    h4_2 = plot(ax4_2, [-120, -75, -60, -135, -120], ...
                       [  35,  35,  60,   60,   35], 'color', line_colors(7, :), 'LineWidth', 1.5);

    
    %tags(isnan(tags(:))) = 0;
    %temp = tags;
    %[c4, h4] = contour(ax4, [X - 360, X], [Y, Y], [temp, temp], [0.9, 0.9], 'linestyle', '-', 'linecolor', line_colors(7, :));
    %h4.LineWidth = 1.0;
    
    axis([ax4_1, ax4_2], [X(1, 1) - 180, X(1, end) - 180, y_min, y_max])
    ax4_1.Visible = 'off';
    ax4_1.XTick = [];
    ax4_1.YTick = [];
    ax4_2.Visible = 'off';
    ax4_2.XTick = [];
    ax4_2.YTick = [];

    
    % calibrate the positions of axes
    set([ax1, ax1_1, ax1_2, ax2, ax3, ax4_1, ax4_2], 'Position', position);
    set(ax1,'TickDir','out');
    if fig_index == 1
        ddd = 1 - down_margin - up_margin;
        cb = colorbar(ax1,'Position',[1 - right_margin + 0.02, ddd/8 + down_margin, 0.015, ddd/4*3]);
        cb.Label.String = '(\%/K)';
        cpos = get(cb,'Position');
        cb.Label.Position = [cpos(1) - 0.18, cpos(2) - 0.355];
        cb.Label.Rotation = 0;
        cb.Label.Interpreter = 'latex';
        cb.Ticks = [-0.12:0.06:0.12];
        cb.TickLabels = {'-12', '-6', '0', '6', '12'};
        cb.TickLabelInterpreter = 'latex';
        set(cb, 'FontSize', 10);
    else
        colorbar(ax1, 'off')
    end
    linkaxes([ax1, ax1_1, ax1_2, ax2, ax3, ax4_1, ax4_2], 'xy');
  
