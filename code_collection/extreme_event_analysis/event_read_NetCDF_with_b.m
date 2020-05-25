function [T_out, ug_out, vg_out, omega_out, omega_b_out] = ...
    event_read_NetCDF_with_b(input_year, input_path, event_latspan, event_lonspan, time_indices, level_indices, string_1, JJA_DJF)

    g = 9.8;

    if input_year < 1990 || (input_year > 2005 && input_year < 2071) || input_year > 2080
        disp(['year ', num2str(input_year), ' does not exist in the data, abandoned']);
        return;   
    end
    if input_year >= 1990 && input_year <= 2005
        string_2 = 'h';
    elseif input_year >= 2071 && input_year <= 2080
        string_2 = 'r';
    end


    start = [1, event_latspan(1), level_indices(1), time_indices(1)];
    count = [Inf, length(event_latspan), length(level_indices), length(time_indices)];
    %start = [1, 1, 1, time_indices(1)];
    %count = [Inf, Inf, Inf, length(time_indices)];

    T = ncread([input_path, 'ta_', string_1, '_', num2str(input_year), '.nc'], 'ta_plevel', start, count);
    T = T(event_lonspan, :, :, :);
    ug = ncread([input_path, 'ug_', string_1, '_', num2str(input_year), '.nc'], 'ug', start, count);
    ug = ug(event_lonspan, :, :, :);
    vg = ncread([input_path, 'vg_', string_1, '_', num2str(input_year), '.nc'], 'vg', start, count);
    vg = vg(event_lonspan, :, :, :);
    omega = ncread([input_path, 'omega_', string_1, '_', num2str(input_year), '.nc'], 'omega_plevel', start, count);
    omega = omega(event_lonspan, :, :, :);

    % add climatologically averaged omega as lateral boundary condition
    start = [1, event_latspan(1), level_indices(1)];
    count = [Inf, length(event_latspan), length(level_indices)];
    if ~exist('JJA_DJF') || JJA_DJF == 0
        omega_b = ncread([input_path, 'omega_avg_', string_1, '_', string_2, '.nc'], 'omega_avg', start, count);
    elseif JJA_DJF == 1
        omega_b = ncread([input_path, 'omega_avg_', string_2, '_JJA', '.nc'], 'omega_avg', start, count);
    elseif JJA_DJF == 2
        omega_b = ncread([input_path, 'omega_avg_', string_2, '_DJF', '.nc'], 'omega_avg', start, count);
    end    
    omega_b = omega_b(event_lonspan, :, :);
    
    [T_out, ug_out, vg_out, omega_out] = ...
        deal(zeros(length(event_latspan), length(event_lonspan), length(level_indices), length(time_indices)));
    omega_b_out = zeros(length(event_latspan), length(event_lonspan), length(level_indices));
    for k = 1 : length(level_indices)
        for t = 1 : length(time_indices)
            T_out(:, :, k, t) = T(:, :, k, t)';
            ug_out(:, :, k, t) = ug(:, :, k, t)';
            vg_out(:, :, k, t) = vg(:, :, k, t)';
            omega_out(:, :, k, t) = omega(:, :, k, t)';
        end
        omega_b_out(:, :, k) = omega_b(:, :, k)';
    end


end
