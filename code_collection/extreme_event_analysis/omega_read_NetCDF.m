function omega_out = ...
    omega_read_NetCDF(input_year, input_path, event_latspan, event_lonspan, time_indices, level_indices, string_1)

    g = 9.8;

    if input_year < 1990 || (input_year > 2005 && input_year < 2071) || input_year > 2080
        disp(['year ', num2str(input_year), ' does not exist in the data, abandoned']);
        return;   
    end

    start = [1, event_latspan(1), level_indices(1), time_indices(1)];
    count = [Inf, length(event_latspan), length(level_indices), length(time_indices)];

    omega = ncread([input_path, 'omega_', string_1, '_', num2str(input_year), '.nc'], 'omega_plevel', start, count);
    omega = omega(event_lonspan, :, :, :);

    
    [omega_out] = ...
        deal(zeros(length(event_latspan), length(event_lonspan), length(level_indices), length(time_indices)));
    for k = 1 : length(level_indices)
        for t = 1 : length(time_indices)
            omega_out(:, :, k, t) = omega(:, :, k, t)';
        end
    end


end
