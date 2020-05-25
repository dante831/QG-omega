
% calculate precipitation quantile

function [precip, precip_quantile] = percentile_grid_with_mask(n, lon_indices, lat_indices, files, percentage, MASK, h_max, ocean_only, p_mask, months)

    % interpolate from pr data field to variable data field

    time = ncread(files{1}, 'time');
    nt = length(time) / n * length(files);
    precip = nan([length(lon_indices), length(lat_indices), nt]);
    precip_zonal = nan([length(lat_indices), length(lon_indices) * nt]);
    temp_weight = [0.5, ones(1, n - 1), 0.5];
    weights = repmat(reshape(temp_weight / sum(temp_weight), 1, 1, n + 1), length(lon_indices), length(lat_indices), 1);
    temp_weight_end = [0.5, ones(1, n - 1)];
    weights_end = repmat(reshape(temp_weight_end / sum(temp_weight_end), 1, 1, n), length(lon_indices), length(lat_indices), 1);
    
    for f = 1 : length(files)
        
        filename = files{f};
        
        % compute 6-hourly mean precipitation
        % ---0---|---0---|---0---|---0---|  ---0---|---0---|---0---|---0---| 6 hourly mean precip data
        % -------0-------0-------0-------0  -------0-------0-------0-------0 6 hourly variable data
        % ---|---0---|---0---|---0---|---0  ---|---0---|---0---|---0---|---0 6 hourly precip
        % ---|---------------0------------  ---|---------------0------------ daily precip
        
        for t = 1 : length(time) / n % every file contains five years' precipitation data
            if t == length(time) / n && f == length(files)
                
                start = [1, 1, (t - 1) * n + 1];
                count = [Inf, Inf, n];
                temp = ncread(filename, 'PRECT', start, count);
                %precip(:, :, (f - 1) * length(time) / n + t) = mean(temp(lon_indices, lat_indices, :), 3);
                precip(:, :, (f - 1) * length(time) / n + t) = sum(temp(lon_indices, lat_indices, :) .* weights_end, 3);

            elseif t == length(time) / n && f ~= length(files)
                
                start = [1, 1, (t - 1) * n + 1];
                count = [Inf, Inf, n];
                temp1 = ncread(filename, 'PRECT', start, count);
                
                start = [1, 1, 1];
                count = [Inf, Inf, 1];
                temp2 = ncread(files{f + 1}, 'PRECT', start, count);
                %precip(:, :, (f - 1) * length(time) / n + t) = ...
                %        (temp2(lon_indices, lat_indices, :) + mean(temp1(lon_indices, lat_indices, :), 3) * (n - 1)) / n;
                precip(:, :, (f - 1) * length(time) / n + t) = ...
                        sum(temp1(lon_indices, lat_indices, :) .* weights(:, :, 1 : n) + ...
                            temp2(lon_indices, lat_indices, :) .* weights(:, :, n + 1), 3);
                
            else
                
                start = [1, 1, (t - 1) * n + 1];
                count = [Inf, Inf, n + 1];
                temp = ncread(filename, 'PRECT', start, count);
                %precip(:, :, (f - 1) * length(time) / n + t) = mean(temp(lon_indices, lat_indices, :), 3);
                precip(:, :, (f - 1) * length(time) / n + t) = sum(temp(lon_indices, lat_indices, :) .* weights, 3);
                
            end
            
        end

    end

    % do masking
    if MASK
        small_value = 0.0;
        precip = pressure_masking(precip, small_value, lon_indices, lat_indices, p_mask, files{1});
    end

    % take zonal average and convert to mm/day
    precip = precip * 86400 * 1000;
    precip_quantile = zeros(length(lon_indices), length(lat_indices));

    % calculate quantile

    for j = 1 : length(lat_indices)
        for i = 1 : length(lon_indices)
            precip_t = squeeze(precip(i, j, :));
            %precip_quantile(i, j) = quantile(precip_t(precip_t ~= 0), percentage); % remove non-precipitating points
            % Paul's correction, use the full pdf instead of only precipitating days, according to DOI:10.1007/s10584-016-1669-2
            precip_quantile(i, j) = quantile(precip_t, percentage);
        end
    end
    
end

