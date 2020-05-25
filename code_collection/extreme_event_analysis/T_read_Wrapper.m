function [T, tag] = T_read_Wrapper(Ny, Nx, Nz, ...
        event_year, input_path, event_latspan, event_lonspan, event_day, level_indices, string_1)

    T = deal(nan(Ny, Nx, Nz, 3));
    if mod(event_day, 365) == 0.25
    
        time_indices = 365 / 0.25;
        input_year = event_year - 1;
        T(:, :, :, 1) = ...
                T_read_NetCDF(input_year, input_path, event_latspan, event_lonspan, time_indices, level_indices, string_1);
        time_indices = [1, 2];
        input_year = event_year;
        T(:, :, :, 2:3) = ...
                T_read_NetCDF(input_year, input_path, event_latspan, event_lonspan, time_indices, level_indices, string_1);
    
    elseif mod(event_day, 365) == 0
    
        time_indices = 365 / 0.25 + [-1, 0];
        input_year = event_year;
        T(:, :, :, 1:2) = ...
                T_read_NetCDF(input_year, input_path, event_latspan, event_lonspan, time_indices, level_indices, string_1);
        time_indices = 1;
        input_year = event_year + 1;
        if input_year == 2006 || input_year == 2081
            disp(['year ', num2str(input_year), ' does not exist in the data, abandoned']);
            tag = false;
            return 
        end
        T(:, :, :, 3) = ...
                T_read_NetCDF(input_year, input_path, event_latspan, event_lonspan, time_indices, level_indices, string_1);
    
    else
    
        time_indices = mod(event_day, 365) / 0.25 + [-1, 0, 1];
        input_year = event_year;
        T = ...
                T_read_NetCDF(input_year, input_path, event_latspan, event_lonspan, time_indices, level_indices, string_1);
    
    end
    tag = true;

end
