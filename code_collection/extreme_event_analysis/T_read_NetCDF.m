function T_out = ...
    T_read_NetCDF(input_year, input_path, event_latspan, event_lonspan, time_indices, level_indices, string_1)

    g = 9.8;

    if ~isempty(strfind(input_path, 'CESM'))
        if input_year < 1990 || (input_year > 2005 && input_year < 2071) || input_year > 2080
            disp(['year ', num2str(input_year), ' does not exist in the data, abandoned']);
            return;   
        end
    elseif ~isempty(strfind(input_path, 'GFDL'))
        if input_year >= 2000 && input_year < 2080
            disp(['year ', num2str(input_year), ' does not exist in the data, abandoned']);
            return;
        end
    end

    start = [1, event_latspan(1), level_indices(1), time_indices(1)];
    count = [Inf, length(event_latspan), length(level_indices), length(time_indices)];

    if ~isempty(string_1)
        T = ncread([input_path, 'ta_', string_1, '_', num2str(input_year), '.nc'], 'ta_plevel', start, count);
    else
        T = ncread([input_path, 'ta_', num2str(input_year), '.nc'], 'ta_plevel', start, count);
    end
    T = T(event_lonspan, :, :, :);

    
    [T_out] = ...
        deal(zeros(length(event_latspan), length(event_lonspan), length(level_indices), length(time_indices)));
    for k = 1 : length(level_indices)
        for t = 1 : length(time_indices)
            T_out(:, :, k, t) = T(:, :, k, t)';
        end
    end


end
