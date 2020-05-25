function omega = omega_from_continuity_v1(ua_filename, va_filename, ps_filename, lon_indices, lat_indices, time_indices)
   
    % This function calculates vertical velocity omega from continuity equation 
    % (2.5) in Simmons & Burridge 1981 paper. 

    % some constants:

    R = 6371000.0;

    % latitudinal and longitudinal indices
    
    lat_series = ncread(ua_filename, 'lat');
    lon_series = ncread(ua_filename, 'lon');

    if lon_indices(1) - 1 < 1
        lon_min = length(lon_series);
    else
        lon_min = lon_indices(1) - 1;
    end
    if lon_min - 1 < 1
        lon_min_2 = length(lon_series);
    else
        lon_min_2 = lon_min - 1;
    end
    if lon_indices(end) + 1 > length(lon_series)
        lon_max = 1;
    else
        lon_max = lon_indices(end) + 1;
    end
    if lon_max + 1 > length(lon_series)
        lon_max_2 = 1;
    else
        lon_max_2 = lon_max + 1;
    end
    lon_indices2 = [lon_min_2; lon_min; lon_indices; lon_max; lon_max_2];

    lat = lat_series(lat_indices);
    lon = lon_series(lon_indices2);

    phi = lat / 180.0 * 3.1415926;
    lambda = lon / 180.0 * 3.1415926;
    dphi = (lat_series(end) - lat_series(1)) / (length(lat_series) - 1) / 180.0 * 3.1415926; 
    dlambda = (lon_series(end) - lon_series(1)) / (length(lon_series) - 1) / 180.0 * 3.1415926;

    % selecting data
    start = [1, lat_indices(1), 1, time_indices(1)];
    count = [Inf, length(lat_indices), Inf, length(time_indices)];
    ua = ncread(ua_filename, 'U', start, count);
    va = ncread(va_filename, 'V', start, count);
    ua = ua(lon_indices2, :, :, :);
    va = va(lon_indices2, :, :, :);
    start = [1, lat_indices(1), time_indices(1)];
    count = [Inf, length(lat_indices), length(time_indices)];
    ps = ncread(ps_filename, 'PS', start, count);
    ps = ps(lon_indices2, :, :);
    p0 = ncread(ua_filename, 'P0');

    % flip the vertical coordinates from top-down to bottom-up
    a = fliplr(ncread(ua_filename, 'hyam')')';
    b = fliplr(ncread(ua_filename, 'hybm')')';
    a_bnds = fliplr(ncread(ua_filename, 'hyai')')';
    b_bnds = fliplr(ncread(ua_filename, 'hybi')')';
    
    % revert ua and va structure
    temp_ua = ua;
    temp_va = va;
    for k = 1 : length(a)
        ua(:, :, k, :) = temp_ua(:, :, length(a) + 1 - k, :);
        va(:, :, k, :) = temp_va(:, :, length(a) + 1 - k, :);
    end
    clear('temp_ua', 'temp_va')

    % calculate pressure in Pa
    [p, p_lower_bnds, p_upper_bnds] = deal(zeros([length(lon_indices2), length(lat_indices), length(a), length(time_indices)]));
    for t = 1 : length(time_indices)
        for k = 1 : length(a)
            for j = 1 : length(lat_indices)
                for i = 1 : length(lon_indices2)
                    p(i, j, k, t) = a(k) * p0 + b(k) * ps(i, j, t);
                    % in this context, p_bnds[:, :, k, :] = p_{k + 1/2}
                    p_lower_bnds(i, j, k, t) = a_bnds(k)     * p0 + b_bnds(k)     * ps(i, j, t);
                    p_upper_bnds(i, j, k, t) = a_bnds(k + 1) * p0 + b_bnds(k + 1) * ps(i, j, t);
                end
            end
        end
    end
    
    % calculate omega by (2.5, 3.12, 3.13) in Simmons and Burrige, 1981

    omega = zeros([length(lon_indices), length(lat_indices), length(a), length(time_indices)]);

    p_lower_bnds  = permute(p_lower_bnds, [3, 1, 2, 4]);
    p_upper_bnds  = permute(p_upper_bnds, [3, 1, 2, 4]);
    p             = permute(p           , [3, 1, 2, 4]);
    ua            = permute(ua          , [3, 1, 2, 4]);
    va            = permute(va          , [3, 1, 2, 4]);
    omega         = permute(omega       , [3, 1, 2, 4]);
    dim = size(omega);
    Phi = repmat(reshape(phi, [1, 1, dim(3), 1]), [dim(1), dim(2), 1, dim(4)]);
    i_ind = 3 : length(lon_indices) + 2;
    j_ind = 2 : length(lat_indices) - 1;
    j_ind_2 = 3 : length(lat_indices) - 2;
    delta_p = p_lower_bnds - p_upper_bnds;
    [temp_omega] = deal(zeros([length(a), 1]));
    temp_k = length(a) : -1 : 1;

    alpha = 1 - p_upper_bnds ./ delta_p .* ...
            log(p_lower_bnds ./ p_upper_bnds);
    coef = log(p_lower_bnds ./ p_upper_bnds);
    [temp_2, temp] = deal(zeros(size(ua)));

    % boundaries
    for t = 1 : length(time_indices)
        for i = 3 : length(lon_indices) + 2
            j = 1; % second order one-sided finite difference
            temp(:, i, j, t) = - (1 / (R * cos(phi(j))) * ...
                            (- 1 / 2 * cos(phi(j + 2)) * va(:, i, j + 2, t) .* delta_p(:, i, j + 2, t) ...
                             + 2     * cos(phi(j + 1)) * va(:, i, j + 1, t) .* delta_p(:, i, j + 1, t) ...
                             - 3 / 2 * cos(phi(j))     * va(:, i, j,     t) .* delta_p(:, i, j    , t)) / dphi);
            for j = [2, length(lat_indices) - 1]
                temp(:, i, j, t) = - (1 / (R * cos(phi(j))) .* ...
                            (  1 / 2 * cos(phi(j + 1)) .* va(:, i, j + 1, t) .* delta_p(:, i, j + 1, t) ...
                             - 1 / 2 * cos(phi(j - 1)) .* va(:, i, j - 1, t) .* delta_p(:, i, j - 1, t)) / dphi);
            end
            j = length(lat_indices);
                temp(:, i, j, t) = - (1 / (R * cos(phi(j))) .* ...
                            (- 1 / 2 * cos(phi(j))     .* va(:, i, j,     t) .* delta_p(:, i, j    , t) ...
                             + 2     * cos(phi(j - 1)) .* va(:, i, j - 1, t) .* delta_p(:, i, j - 1, t) ...
                             - 3 / 2 * cos(phi(j - 2)) .* va(:, i, j - 2, t) .* delta_p(:, i, j - 2, t)) / dphi);
        end
    end
    % 4th order accuracy with momentum conservation:
    % interior
    temp(:, i_ind, j_ind_2, :) = - (1 ./ (R * cos(Phi(:, :, j_ind_2, :))) .* ...
            (- 1 / 12 * cos(Phi(:, :, j_ind_2 + 2, :)) .* va(:, i_ind, j_ind_2 + 2, :) .* delta_p(:, i_ind, j_ind_2 + 2, :) ...
             + 2 / 3  * cos(Phi(:, :, j_ind_2 + 1, :)) .* va(:, i_ind, j_ind_2 + 1, :) .* delta_p(:, i_ind, j_ind_2 + 1, :) ...
             - 2 / 3  * cos(Phi(:, :, j_ind_2 - 1, :)) .* va(:, i_ind, j_ind_2 - 1, :) .* delta_p(:, i_ind, j_ind_2 - 1, :) ...
             + 1 / 12 * cos(Phi(:, :, j_ind_2 - 2, :)) .* va(:, i_ind, j_ind_2 - 2, :) .* delta_p(:, i_ind, j_ind_2 - 2, :)) / dphi);
    temp(:, i_ind, :, :) = temp(:, i_ind, :, :) - (1 ./ (R * cos(Phi)) .* (- 1 / 12 * ua(:, i_ind + 2, :, :) .* delta_p(:, i_ind + 2, :, :) ...  
                                                                           + 2 / 3  * ua(:, i_ind + 1, :, :) .* delta_p(:, i_ind + 1, :, :) ...
                                                                           - 2 / 3  * ua(:, i_ind - 1, :, :) .* delta_p(:, i_ind - 1, :, :) ...
                                                                           + 1 / 12 * ua(:, i_ind - 2, :, :) .* delta_p(:, i_ind - 2, :, :)) / dlambda);

    for t = 1 : length(time_indices)
        for j = 1 : length(lat_indices)
            for i = 1 : length(lon_indices)
                summation = 0;
                for k = length(a) - 1 : -1 : 1 % integrate from top-down (in a bottom-up vertical cooridinate)
                    summation = summation + temp(k + 1, i + 2, j, t);
                    temp_2(k, i + 2, j, t) = coef(k, i + 2, j, t) * summation + alpha(k, i + 2, j, t) * temp(k, i + 2, j, t);
                end
                temp_2(end, i + 2, j, t) = alpha(end, i + 2, j, t) * temp(end, i + 2, j, t);
            end
        end
    end
    
    % interior
    % temp_k transforms the coordinate from bottom-up to top-down, consistent with the original pressure storage 
    % (the pressure is top-down on the outside in main_sigma_to_p.m)
    omega(temp_k, :, j_ind, :) = 1 ./ delta_p(:, i_ind, j_ind, :) .* p(:, i_ind, j_ind, :) .* temp_2(:, i_ind, j_ind, :) + ...
                                p(:, i_ind, j_ind, :) .* (va(:, i_ind, j_ind, :) / (R * 2 * dphi) .* ...
                                    (- 1 ./ delta_p(:, i_ind, j_ind + 1, :) .* (p_upper_bnds(:, i_ind, j_ind + 1, :) .* log(p_upper_bnds(:, i_ind, j_ind + 1, :)) - ...
                                                                                p_lower_bnds(:, i_ind, j_ind + 1, :) .* log(p_lower_bnds(:, i_ind, j_ind + 1, :))) + ...
                                       1 ./ delta_p(:, i_ind, j_ind - 1, :) .* (p_upper_bnds(:, i_ind, j_ind - 1, :) .* log(p_upper_bnds(:, i_ind, j_ind - 1, :)) - ...
                                                                                p_lower_bnds(:, i_ind, j_ind - 1, :) .* log(p_lower_bnds(:, i_ind, j_ind - 1, :)))));  ...
    
    % latitudinal boundaries
    for t = 1 : length(time_indices)
        for i = 1 : length(lon_indices)
            j = 1;
            temp_omega(temp_k) = 1 ./ delta_p(:, i + 2, j, t) .* p(:, i + 2, j, t) .* temp_2(:, i + 2, j, t);
            omega(temp_k, i, j, t) = temp_omega(temp_k) + p(:, i + 2, j, t) .* (va(:, i + 2, j, t) / (R * dphi) .* ...
                                (- 1 ./ delta_p(:, i + 2, j + 1, t) .* (p_upper_bnds(:, i + 2, j + 1, t) .* log(p_upper_bnds(:, i + 2, j + 1, t)) - ...
                                                                        p_lower_bnds(:, i + 2, j + 1, t) .* log(p_lower_bnds(:, i + 2, j + 1, t))) + ...
                                   1 ./ delta_p(:, i + 2, j, t)     .* (p_upper_bnds(:, i + 2, j, t)     .* log(p_upper_bnds(:, i + 2, j, t)) - ...
                                                                        p_lower_bnds(:, i + 2, j, t)     .* log(p_lower_bnds(:, i + 2, j, t)))));
            j = length(lat_indices);
            temp_omega(temp_k) = 1 ./ delta_p(:, i + 2, j, t) .* p(:, i + 2, j, t) .* temp_2(:, i + 2, j, t);
            omega(temp_k, i, j, t) = temp_omega(temp_k) + p(:, i + 2, j, t) .* (va(:, i + 2, j, t) / (R * dphi) .* ...
                                (- 1 ./ delta_p(:, i + 2, j, t)     .* (p_upper_bnds(:, i + 2, j, t)     .* log(p_upper_bnds(:, i + 2, j, t)) - ...
                                                                        p_lower_bnds(:, i + 2, j, t)     .* log(p_lower_bnds(:, i + 2, j, t))) + ...
                                   1 ./ delta_p(:, i + 2, j - 1, t) .* (p_upper_bnds(:, i + 2, j - 1, t) .* log(p_upper_bnds(:, i + 2, j - 1, t)) - ...
                                                                        p_lower_bnds(:, i + 2, j - 1, t) .* log(p_lower_bnds(:, i + 2, j - 1, t)))));
        end
    end

    omega(temp_k, :, :, :) = omega(temp_k, :, :, :) + p(:, i_ind, :, :) .* (...
                                ua(:, i_ind, :, :) ./ (R * cos(Phi) * 2 * dlambda) .* ...
                                (- 1 ./ delta_p(:, i_ind + 1, :, :) .* (p_upper_bnds(:, i_ind + 1, :, :) .* log(p_upper_bnds(:, i_ind + 1, :, :)) - ...
                                                                        p_lower_bnds(:, i_ind + 1, :, :) .* log(p_lower_bnds(:, i_ind + 1, :, :))) + ...
                                   1 ./ delta_p(:, i_ind - 1, :, :) .* (p_upper_bnds(:, i_ind - 1, :, :) .* log(p_upper_bnds(:, i_ind - 1, :, :)) - ...
                                                                        p_lower_bnds(:, i_ind - 1, :, :) .* log(p_lower_bnds(:, i_ind - 1, :, :)))));
    
    omega = permute(omega, [2, 3, 1, 4]);

                        


