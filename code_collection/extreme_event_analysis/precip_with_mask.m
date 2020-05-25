function precip = precip_with_mask(n, lon_indices, lat_indices, files, MASK, h_max, ocean_only, p_mask)

    % interpolate from pr data field to variable data field

    time = ncread(files{1}, 'time'); % get the original time series of one file
    nt = length(time) / n * length(files); % calculate the total length of the precip time series
    precip = nan([length(lon_indices), length(lat_indices), nt]);
    
    % define unnormalized weights of the window
    temp_weight = ones(1, n); 
    temp_weight_end = ones(1, n); % boundary at the end of data
    
    % populate over the grid and normalize
    weights = repmat(reshape(temp_weight / sum(temp_weight), 1, 1, n), length(lon_indices), length(lat_indices), 1); 
    weights_end = repmat(reshape(temp_weight_end / sum(temp_weight_end), 1, 1, n), length(lon_indices), length(lat_indices), 1);

    for f = 1 : length(files)

        filename = files{f};

        % compute 6-hourly mean precipitation centered on the available data points
        % ---0---|---0---|---0---|---0---|  ---0---|---0---|---0---|---0---| 6 hourly instantaneous precip data
        % ---0---|---0---|---0---|---0---|  ---0---|---0---|---0---|---0---| 6 hourly variable data
        % ---0---|---0---|---0---|---0---|  ---0---|---0---|---0---|---0---| 6 hourly precip
        % ---------------0---------------|  ---------------0---------------| daily precip

        for t = 1 : length(time) / n % every file contains five years' precipitation data

            start = [1, 1, (t - 1) * n + 1];
            count = [Inf, Inf, n];
            temp = ncread(filename, 'PRECT', start, count);
            precip(:, :, (f - 1) * length(time) / n + t) = sum(temp(lon_indices, lat_indices, :) .* weights, 3);

        end

    end

    % do masking
    if MASK
        small_value = 0.0;
        % set grid points with climatological pressure values that are below p_mask to small_value, 
        % so that they will not exceed the 99.9th percentile, and therefore excluded from analysis
        precip = pressure_masking(precip, small_value, lon_indices, lat_indices, p_mask, files{1});
    end

    % convert from m/s to mm/day
    precip = precip * 86400 * 1000;

end
