function multifile_quantile(files, years, percentage, output_name, start_year)

    time = ncread(files{1}, 'time');
    lat = ncread(files{1}, 'lat');
    lon = ncread(files{1}, 'lon');
    Q = zeros([length(lon), length(lat)]);
    N = (years(end) - years(1) + 1) * 4 * 365;
    t_ind = (years(1) - start_year) * 4 * 365 + 1 : (years(end) + 1 - start_year) * 4 * 365;
    P = zeros([length(lon), length(lat), length(t_ind) * length(files)]);

    for f = 1 : length(files)
        disp(['f = ', num2str(f)])
        temp = ncread(files{f}, 'PRECT');
        %{
        temp_2 = zeros(size(temp));
        temp_2(:, :, 1:end-1) = 0.5 * temp(:, :, 1:end-1) + 0.5 * temp(:, :, 2:end);
        temp_2(:, :, end) = temp(:, :, end);
        P(:, :, ((f - 1) * N + 1 : f * N)) = temp_2(:, :, t_ind);
        %}
        P(:, :, ((f - 1) * N + 1 : f * N)) = temp(:, :, t_ind);
    end
    clear('temp');

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

