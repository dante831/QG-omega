function [T, ua, va, ug, vg, omega, omega_b, q, tag] = event_read_Wrapper_uava(Ny, Nx, Nz, ...
        event_year, input_path, event_latspan, event_lonspan, event_day, level_indices, string_1, JJA_DJF)

    [T, ua, va, ug, vg, omega, q] = deal(nan(Ny, Nx, Nz, 3));
    omega_b = nan(Ny, Nx, Nz);
    if mod(event_day, 365) == 0.25
    
        time_indices = 365 / 0.25;
        input_year = event_year - 1;
        [T(:, :, :, 1), ua(:, :, :, 1), va(:, :, :, 1), ug(:, :, :, 1), vg(:, :, :, 1), omega(:, :, :, 1), q(:, :, :, 1), omega_b] =...
                event_read_NetCDF_with_b_uava_q(input_year, input_path, event_latspan, event_lonspan, time_indices, level_indices, string_1, JJA_DJF);
        time_indices = [1, 2];
        input_year = event_year;
        [T(:, :, :, 2:3), ua(:, :, :, 2:3), va(:, :, :, 2:3), ug(:, :, :, 2:3), vg(:, :, :, 2:3), omega(:, :, :, 2:3), q(:, :, :, 2:3), ~] =...
                event_read_NetCDF_with_b_uava_q(input_year, input_path, event_latspan, event_lonspan, time_indices, level_indices, string_1, JJA_DJF);
    
    elseif mod(event_day, 365) == 0
    
        time_indices = 365 / 0.25 + [-1, 0];
        input_year = event_year;
        [T(:, :, :, 1:2), ua(:, :, :, 1:2), va(:, :, :, 1:2), ug(:, :, :, 1:2), vg(:, :, :, 1:2), omega(:, :, :, 1:2), q(:, :, :, 1:2), omega_b] =...
                event_read_NetCDF_with_b_uava_q(input_year, input_path, event_latspan, event_lonspan, time_indices, level_indices, string_1, JJA_DJF);
        time_indices = 1;
        input_year = event_year + 1;
        if input_year == 2006 || input_year == 2081
            disp(['year ', num2str(input_year), ' does not exist in the data, abandoned']);
            tag = false;
            return 
        end
        [T(:, :, :, 3), ua(:, :, :, 3), va(:, :, :, 3), ug(:, :, :, 3), vg(:, :, :, 3), omega(:, :, :, 3), q(:, :, :, 3), ~] =...
        event_read_NetCDF_with_b_uava_q(input_year, input_path, event_latspan, event_lonspan, time_indices, level_indices, string_1, JJA_DJF);
    
    else
    
        time_indices = mod(event_day, 365) / 0.25 + [-1, 0, 1];
        input_year = event_year;
        [T, ua, va, ug, vg, omega, q, omega_b] = ...
                event_read_NetCDF_with_b_uava_q(input_year, input_path, event_latspan, event_lonspan, time_indices, level_indices, string_1, JJA_DJF);
    
    end
    tag = true;

end
