function [T_out, ua_out, va_out, ug_out, vg_out, omega_out, omega_b_out, ...
          lon, lat, level, event_timespan] = read_ncfile(nc_filename)

    var_names = {'T', 'ua', 'va', 'ug', 'vg', 'omega', 'omega_b'};

    lon = ncread(nc_filename, 'longitude');
    lat = ncread(nc_filename, 'latitude');
    level = ncread(nc_filename, 'level');
    event_timespan = ncread(nc_filename, 'time');

    for n = 1 : length(var_names)
        eval([var_names{n}, ' = ncread(nc_filename, ''', var_names{n}, ''');'])
        eval([var_names{n}, '_out = zeros(length(lat), length(lon), length(level), length(event_timespan));'])
    end
    for n = 1 : length(var_names)
        for k = 1 : length(level)
            for t = 1 : length(event_timespan)
                eval([var_names{n}, '_out(:, :, k, t) = ', var_names{n}, '(:, :, k, t)'';'])
            end
        end
    end




