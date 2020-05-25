function multifile_quantile_month(files, years, percentage, output_name, start_year, months)

    mon_days = [31,  59,  90, 120, 151, 181, 212, 243, 273, 304, 334, 365];

    time = ncread(files{1}, 'time');
    lat = ncread(files{1}, 'lat');
    lon = ncread(files{1}, 'lon');
    Q = zeros([length(lon), length(lat)]);
	[day_ind, t_ind] = deal([]);
	for m = 1 : length(months)
    	if months(m) == 1
        	lower_bnd = 1;
	    else
    	    lower_bnd = mon_days(months(m) - 1) * 4 + 1;
	    end
    	day_ind = [day_ind, lower_bnd : mon_days(months(m)) * 4];
	end
	for y = 1 : length(years)
    	t_ind = [t_ind, (years(y) - start_year) * 4 * 365 + day_ind];
	end
    
    %N = (years(end) - years(1) + 1) * 4 * 365;
    N = length(t_ind);
    P = zeros([length(lon), length(lat), length(t_ind) * length(files)]);

    for f = 1 : length(files)
        disp(['f = ', num2str(f)])
        temp = ncread(files{f}, 'PRECT');
        temp_2 = zeros(size(temp));
        temp_2(:, :, 1:end-1) = 0.5 * temp(:, :, 1:end-1) + 0.5 * temp(:, :, 2:end);
        temp_2(:, :, end) = temp(:, :, end);
        P(:, :, ((f - 1) * N + 1 : f * N)) = temp_2(:, :, t_ind);
    end
    clear('temp', 'temp_2');

    for j = 1 : length(lat)
        disp(['j = ', num2str(j)])
        for i = 1 : length(lon)
            temp_p = squeeze(P(i, j, :));
            %Q(i, j) = quantile(temp_p(temp_p ~= 0), percentage) * 86400 * 1000;
            % Paul's correction, use the full pdf instead of only precipitating days, according to DOI:10.1007/s10584-016-1669-2
            Q(i, j) = quantile(temp_p, percentage) * 86400 * 1000;

        end
    end

    save(output_name, 'Q', '-v7.3')

end

