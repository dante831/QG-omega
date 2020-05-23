function positions = set_axes(main_fig, figsize, margins, x_distance, y_distance)

    x_distance  = x_distance / figsize(3);
	y_distance  = y_distance / figsize(4);

	left_margin     = margins(1);
	right_margin    = margins(2);
	down_margin     = margins(3);
	up_margin       = margins(4);

	width_fig   = (1 - left_margin - right_margin - x_distance * (main_fig.col - 1)) / main_fig.col;
	height_fig  = (1 - down_margin - up_margin    - y_distance * (main_fig.row - 1)) / main_fig.row;

    positions = {};
	n_fig = min(main_fig.row*main_fig.col, length(main_fig.axes));
	for fig_ind = 1 : n_fig
		positions{fig_ind} = [left_margin + (mod(fig_ind - 1, main_fig.col)) * (width_fig + x_distance), ...
				1 - (up_margin + ceil(fig_ind / main_fig.col) * (height_fig + y_distance)) + y_distance, ...
				width_fig, height_fig];
	end

	for fig_ind = 1 : n_fig
		set(main_fig.axes{fig_ind}, 'Position', positions{fig_ind});
	end

end
