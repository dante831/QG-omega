
% calculate precipitation quantile

function precip = precip_with_mask_month(n, lon_indices, lat_indices, files, MASK, h_max, ocean_only, p_max, months)

    mon_days = [31,  59,  90, 120, 151, 181, 212, 243, 273, 304, 334, 365];    

    % interpolate from pr data field to variable data field

    time = ncread(files{1}, 'time');
    nt = length(time) / n * length(files);
    precip = nan([length(lon_indices), length(lat_indices), nt]);
    temp_weight = [0.5, ones(1, n - 1), 0.5];
    weights = repmat(reshape(temp_weight / sum(temp_weight), 1, 1, n + 1), length(lon_indices), length(lat_indices), 1);
    temp_weight_end = [0.5, ones(1, n - 1)];
    weights_end = repmat(reshape(temp_weight_end / sum(temp_weight_end), 1, 1, n), length(lon_indices), length(lat_indices), 1);

    for f = 1 : length(files)

        filename = files{f};

        % compute 6-hourly mean precipitation centered on the available data points
        % ---0---|---0---|---0---|---0---|  ---0---|---0---|---0---|---0---| 6 hourly mean precip data
        % -------0-------0-------0-------0  -------0-------0-------0-------0 6 hourly variable data
        % ---|---0---|---0---|---0---|---0  ---|---0---|---0---|---0---|---0 6 hourly precip
        % ---|---------------0------------  ---|---------------0------------ daily precip

        for t = 1 : length(time) / n % every file contains five years' precipitation data
            day = (mod((t - 1) * n, 365 * 4) + n) * 6 / 24; % 6 is 6-hourly, 24 is 24 hours per day
            mon_ind = min(find(day <= mon_days)); % find out which month it is
            if ~any(mon_ind == months)

                precip(:, :, (f - 1) * length(time) / n + t) = 0;

            elseif t == length(time) / n && f == length(files)

                start = [1, 1, (t - 1) * n + 1];
                count = [Inf, Inf, n];
                temp = ncread(filename, 'PRECT', start, count);
                precip(:, :, (f - 1) * length(time) / n + t) = sum(temp(lon_indices, lat_indices, :) .* weights_end, 3);

            elseif t == length(time) / n && f ~= length(files)

                start = [1, 1, (t - 1) * n + 1];
                count = [Inf, Inf, n];
                temp1 = ncread(filename, 'PRECT', start, count);

                start = [1, 1, 1];
                count = [Inf, Inf, 1];
                temp2 = ncread(files{f + 1}, 'PRECT', start, count);
                precip(:, :, (f - 1) * length(time) / n + t) = ...
                        sum(temp1(lon_indices, lat_indices, :) .* weights(:, :, 1 : n) + ...
                            temp2(lon_indices, lat_indices, :) .* weights(:, :, n + 1), 3);

            else

                start = [1, 1, (t - 1) * n + 1];
                count = [Inf, Inf, n + 1];
                temp = ncread(filename, 'PRECT', start, count);
                precip(:, :, (f - 1) * length(time) / n + t) = sum(temp(lon_indices, lat_indices, :) .* weights, 3);

            end

        end

    end

    % do masking
    if MASK
        small_value = 0.0;
        precip = pressure_masking(precip, small_value, lon_indices, lat_indices, p_max, files{1});
    end

    % convert from m/s to mm/day
    precip = precip * 86400 * 1000;

end
