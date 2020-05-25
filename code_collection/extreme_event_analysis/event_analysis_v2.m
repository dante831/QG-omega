
% compute event-wise QG_omega for precipitation extremes
% final version v2.0, edited on May 17th, 2019 by Ziwei Li

%{
precip_ = event_precip_historical;
days_ = days_historical;
jj = 1;
ii = 1;
m = 1;
%}

function event_analysis_v2(box_min, box_max, p_min, p_max, lat_series_original, lon_series_original, lat_indices, lon_indices, ...
                                    days_, precip_, matfilename, input_path, output_path, NETCDF, NETCDF_interval, n_ensemble, ...
                                    FIND_CENTER, JJA_DJF, ACCU_SIGMA, TRADITIONAL_ADV, DAILY, ACCU_B, SMOOTH, FIND_CENTER_ADV)
    
    % definitions of useful constants and arrays

    Omega = 7.2921e-5;  % Earth's angular speed in rad/s 
    R = 6371000.0;      % average radius from https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html, in meters
    Ra = 287.04;        % specific gas constant, in J/(kg*K)
    cp = 1005.0;        % heat capacity at 300K, in J/(kg*K)
    kappa = Ra/cp; 
    dt = 3600.0 * 6.0;  % time difference, in seconds
    window_x = 1;       % smoothing window. This translates to a (2*window_x+1)-by-(2*window_y+1) smoothing filter
    window_y = 1;
    NetCDF_label = 0;   % NetCDF file counter
            
    lon_series = lon_series_original(lon_indices); % longitude series that covers the regions in the analysis
    lat_series = lat_series_original(lat_indices); % latitude series that covers the regions in the analysis
    dlambda    = (lon_series(2) - lon_series(1)) / 180.0 * 3.1415926; % x angular increment
    dphi       = (lat_series(2) - lat_series(1)) / 180.0 * 3.1415926; % y angular increment
    
    [num_event, num_event_discarded, num_event_discarded_misalign] = ...
            deal(zeros(length(lat_indices), length(lon_indices))); % define array for number of events
    events_sta = deal(cell(length(lat_indices), length(lon_indices))); % define cell array for "event statistics"
    pwd_str = pwd; % get current path

    if ~exist('DAILY')
        DAILY = false;
    end

    if ~exist('FIND_CENTER_ADV')
        FIND_CENTER_ADV = false;
    end
    
    %[computation_x, computation_y] = deal(box_max);
    box_max_x = box_max;
    box_max_y = box_max;

    for jj = 1 : length(lat_indices)
        for ii = 1 : length(lon_indices)
            disp(['lat: ', num2str(lat_series(jj)), 'lon: ', num2str(lon_series(ii))]);
            if ~isempty(days_{jj, ii})

                for m = 1 : length(days_{jj, ii})
                    disp(['m = ', num2str(m)])
                    event_day = days_{jj, ii}(m); 
                    event_lon = lon_indices(ii);
                    event_lat = lat_indices(jj);    
                    computation_x = box_max_x;
                    computation_y = box_max_y;

                    if strfind(pwd_str, 'CESM')
                        string_1 = num2str(n_ensemble, '%.3d'); % string for the ensemble number 
                        if DAILY
                            event_year = floor((event_day - 1.0) / 365) + 1990; % daily
                        else
                            event_year = floor((event_day - 0.25) / 365) + 1990; % 6-hourly
                        end

                        if event_year >= 1990 && event_year <= 2005
                            string_2 = 'historical';
                            tag = true;
                        elseif event_year >= 2071 && event_year <= 2080
                            string_2 = 'rcp85';
                            tag = true;
                        else
                            disp(['year ', num2str(event_year), ' does not exist in the data, abandoned']);
                            tag = false;
                        end
                    elseif strfind(pwd_str, 'GFDL')
                        if DAILY
                            event_year = floor((event_day - 1.0) / 365) + 1980; % daily
                        else
                            event_year = floor((event_day - 0.25) / 365) + 1980; % 6-hourly
                        end
                        if event_year < 2000
                            string_1 = 'historical';
                            tag = true;
                        elseif event_year >= 2080
                            string_1 = 'rcp85';
                            tag = true;
                        else
                            disp(['year ', num2str(event_year), ' does not exist in the data, abandoned']);
                            tag = false;    
                        end
                    end
                    if ~tag; continue; end

                    % calculate lat and lon spans for the event
                    [event_latspan, event_lonspan, center_y, computation_y] = ...
                            event_span_Wrapper(lat_series_original, lon_series_original, event_lat, event_lon, ...
                            box_max_y, computation_x, computation_y);

                    % reset the level indices for each event
                    plevels = double(ncread([input_path, 'ta_', string_1, '_', num2str(event_year), '.nc'], 'level'));
                    p_start_0 = find(plevels == p_max); % lower boundary start point
                    p_end = find(plevels == p_min); % highest level
                    level_indices = 1 : p_end;
                    
                    if DAILY
                        event_timespan = event_day - 1 + [0.25, 0.50, 0.75, 1.00];
                        time_indices = mod(event_timespan - 0.25,  365) / 0.25 + 1;
                    else
                        event_timespan = event_day + [- 0.25, 0.00, 0.25]; 
                    end
                    
                    % read in event data
                    if DAILY
                        if strfind(pwd_str, 'CESM')
                            [T, ua, va, ug, vg, omega, omega_b_clim] = event_read_NetCDF_with_b_uava...
                                (event_year, input_path, event_latspan, event_lonspan, ...
                                time_indices, level_indices, string_1, JJA_DJF);
                        elseif strfind(pwd_str, 'GFDL')
                            [T, ua, va, ug, vg, omega, omega_b_clim] = event_read_NetCDF_with_b_uava_GFDL...
                                (event_year, input_path, event_latspan, event_lonspan, ...
                                time_indices, level_indices);
                        end
                    else
                        if strfind(pwd_str, 'CESM')
                            if strfind(pwd_str, '005')
                                [T, ua, va, ug, vg, omega, omega_b_clim, q, tag] = ...
                                    event_read_Wrapper_uava_q(computation_y, computation_x, length(level_indices), ...
                                        event_year, input_path, event_latspan, event_lonspan, event_day, ...
                                        level_indices, string_1, JJA_DJF);
                            else
                                [T, ua, va, ug, vg, omega, omega_b_clim, tag] = ...
                                    event_read_Wrapper_uava(computation_y, computation_x, length(level_indices), ...
                                        event_year, input_path, event_latspan, event_lonspan, event_day, ...
                                        level_indices, string_1, JJA_DJF);
                            end
                        elseif strfind(pwd_str, 'GFDL')
                            [T, ua, va, ug, vg, omega, omega_b_clim, tag] = ...
                                event_read_Wrapper_uava_GFDL(computation_y, computation_x, length(level_indices), ...
                                    event_year, input_path, event_latspan, event_lonspan, event_day, ...
                                    level_indices, string_1);
                        end
                    end

                    % shrink_box, if range as large as box_max has NaN value, then shrink until box_min
                    [x_ind, y_ind, computation_x, computation_y, p_start, NaN_tag] = ...
                            shrink_box_v1(omega, box_min, p_max, computation_x, computation_y, p_start_0, p_end, plevels);
                    if NaN_tag; disp('unavailable data in omega, continue');continue; end;

                    % subselect all fields and indices according to the result from shrink_box
                    level_indices = p_start : p_end;
                    level = plevels(level_indices);
                    T       = T      (y_ind, x_ind, level_indices, :);
                    ug      = ug     (y_ind, x_ind, level_indices, :);
                    vg      = vg     (y_ind, x_ind, level_indices, :);
                    ua      = ua     (y_ind, x_ind, level_indices, :);
                    va      = va     (y_ind, x_ind, level_indices, :);
                    omega   = omega  (y_ind, x_ind, level_indices, :);
                    if ~DAILY && ~isempty(strfind(pwd_str, '005'))
                        q       = q      (y_ind, x_ind, level_indices, :);
                    end

                    % omega boundaries
                    omega_b_clim = omega_b_clim(y_ind, x_ind, level_indices, :);
                    omega_b_full = omega;
                    
                    event_latspan = event_latspan(y_ind);
                    event_lonspan = event_lonspan(x_ind);

                    % define some other useful arrays
                    lat = lat_series_original(event_latspan);
                    lon = lon_series_original(event_lonspan);
                    phi = lat / 180.0 * 3.1415926;
                    lambda = lon / 180.0 * 3.1415926;
                    [Lambda, Phi] = meshgrid(lambda, phi);

                    % smooth T field
                    T_smoothed = T;
                    for k = 1 : length(level)
                        for t = 1 : length(event_timespan)
                            if SMOOTH
                                T_smoothed(:, :, k, t) = smooth2a(T(:, :, k, t), window_x, window_y);
                                % renormalize so that the smoothed field has the same variance as the previous field
                                T_smoothed(:, :, k, t) = smooth_normalize(T(:, :, k, t), T_smoothed(:, :, k, t));
                            end
                        end
                    end
                    T = T_smoothed;
                    clear('T_smoothed');

                    % define Corriolis parameters
                    f0 = 2 * Omega * sin(lat_series_original(lat_indices(jj)) / 180 * 3.1415926);
                    beta = 2 * Omega / R * cos(lat_series_original(lat_indices(jj)) / 180 * 3.1415926);
                    
                    % define place-holders for the inversion process
                    [A, B, C, rhs, J, Qx, Qy, ...
                     dug_dlambda, dvg_dlambda, dug_dphi, dvg_dphi, ...
                     dT_dlambda, dT_dphi, sigma_accu, SSigma] = ...
                                deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
                    if DAILY
                        array_t = 1;
                    else
                        array_t = length(event_timespan);
                    end
                    [sigma, T_avg, dtheta_dp_ma_avg] = deal(zeros(length(level), array_t));
                    [omega_QG, omega_QG_full, omega_QG_Adv, omega_QG_full_b, omega_QG_interior, ...
                     omega_QG_lateral_b, omega_QG_lateral_clim_b] = ...
                            deal(zeros(length(event_latspan), length(event_lonspan), length(level), array_t));
                    
                    if DAILY
                        % zonal- and time-averaged temperature
                        T_avg = reshape(mean(mean(mean(T, 1), 2), 4), size(level));
                        [dtheta_dp_ma_avg, ~, theta_avg] = eff_stat_stab(level, T_avg);
 
                        % static stability: sigma
                        sigma =  - Ra * T_avg ./ (level .* theta_avg) .* gradient(theta_avg, level);
    
                        % spatially varying static stability: sigma_accu
                        % this quantity is mainly used to calculate J term in this daily case
                        Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
                        theta = mean(T, 4) .* (1e5 ./ Level) .^ kappa;
                        sigma_accu = - Ra * mean(T, 4) ./ (Level .* theta) .* d_dp(theta, level);
                        
                        % according to Paul, J should be calculated using 4 time slices
                        sigma_accu_t = zeros(length(event_latspan), length(event_lonspan), length(level), 4);
                        for t = 1 : length(event_timespan)
                            theta_t = T(:, :, :, t) .* (1e5 ./ Level) .^ kappa;
                            sigma_accu_t(:, :, :, t) = - Ra * T(:, :, :, t) ./ (Level .* theta_t) .* d_dp(theta_t, level);
                        end
                    else
                        for t = 1 : length(event_timespan)

                            T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));
                            [dtheta_dp_ma_avg(:, t), ~, theta_avg] = eff_stat_stab(level, T_avg(:, t));

                            % static stability: sigma
                            sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg) .* gradient(theta_avg, level);

                            % spatially varying static stability: sigma_accu
                            Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
                            theta = T(:, :, :, t) .* (1e5 ./ Level) .^ kappa;
                            sigma_accu(:, :, :, t) = - Ra * T(:, :, :, t) ./ (Level .* theta) .* d_dp(theta, level);

                        end
                    end
                                       
                    % Adv forcing
                    if TRADITIONAL_ADV % traditional form
                        % get relative vorticity
                        zeta = zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan));
                        for t = 1 : length(event_timespan)
                            for k = 1 : length(level)
                                zeta(:, :, k, t) = curl_spherical(phi, ug(:, :, k, t), vg(:, :, k, t), dphi, dlambda);
                            end
                        end
                        [A, A1, A2, A3] = traditional_A(phi, level, ug, vg, zeta, event_timespan, dphi, dlambda, f0, beta);
                        B = traditional_B(phi, level, ug, vg, T, event_timespan, dphi, dlambda, f0, true);
                    else % Q vector
                        [A, B] = Q_vector_v1(level, ug, vg, T, Phi, event_timespan, dphi, dlambda, f0, beta, true);
                    end

                    if DAILY
                        A = mean(A, 4);
                        B = mean(B, 4);
                    end
                  
                    % term C (diabatic heating)
                    if DAILY
                        %[~, ~, J1, ~] = diabatic_heating(phi, T, level, ua, va, sigma_accu, omega, ...
                        %            event_timespan, dphi, dlambda, dt);
                        %J2 = - sigma_accu .* mean(omega, 4) / Ra .* Level * cp;
                        %J = mean(J1, 4) + J2;
                        [C1, C2, J1, J2] = diabatic_heating(phi, T, level, ua, va, sigma_accu_t, omega, ...
                                    event_timespan, dphi, dlambda, dt);
                        J = mean(J1 + J2, 4);
                        C = mean(C1 + C2, 4); 
                            %this should be exactly the same as a level-by-level calculation of C
                    else
                        [C1, C2, J1, J2] = diabatic_heating(phi, T, level, ua, va, sigma_accu, omega, ...
                                                    event_timespan, dphi, dlambda, dt);
                        J = J1 + J2;
                        C = C1 + C2;
                    end

                    % right hand side
                    rhs = A + B + C;
                     
                    % remove negative values in sigma_accu for inversion. 
                    % Note: we use unsmoothed sigma_accu for diabatic heating, and smoothed sigma_accu for inversion
                    [sigma_accu, tag_sigma] = sigma_remove_negative(sigma_accu, sigma);
                    if ~tag_sigma; disp('too many points with negative stability!');continue; end

                    % smooth sigma_accu field while preserving first and second moment
                    sigma_accu_smoothed = sigma_accu;
                    for k = 1 : length(level)
                        for t = 1 : array_t
                            if SMOOTH
                                sigma_accu_smoothed(:, :, k, t) = smooth2a(sigma_accu(:, :, k, t), window_x, window_y);
                                % renormalize so that the smoothed field has the same variance as the previous field
                                sigma_accu_smoothed(:, :, k, t) = smooth_normalize(sigma_accu(:, :, k, t), sigma_accu_smoothed(:, :, k, t));
                            end
                        end
                    end
                    sigma_accu = sigma_accu_smoothed;
                    clear('sigma_accu_smoothed');

                    % do the inversion
                    Zeros =  zeros(size(omega(:, :, :, t)));
                    for t = 1 : array_t
                        % default inversion, with zero upper boundary, zero lower boundary, and climatological lateral boundaries
                        omega_QG(:, :, :, t)          = SIP_diagnostic_v1_sigma(R, f0, length(lon), length(phi), ...
                                    length(level), phi, lambda, level, ...
                                    sigma_accu(:, :, :, t), dphi, dlambda, rhs(:, :, :, t), omega_b_clim, false);
                        % inversion with zero upper boundary, full lower and lateral boundaries
                        omega_QG_full(:, :, :, t)     = SIP_diagnostic_v1_sigma(R, f0, length(lon), length(phi), ...
                                    length(level), phi, lambda, level, ...
                                    sigma_accu(:, :, :, t), dphi, dlambda, rhs(:, :, :, t), omega_b_full(:, :, :, t), true);
                        % inversion with only Adv in the rhs
                        omega_QG_Adv(:, :, :, t)      = SIP_diagnostic_v1_sigma(R, f0, length(lon), length(phi), ...
                                    length(level), phi, lambda, level, ...
                                    sigma_accu(:, :, :, t), dphi, dlambda, A(:, :, :, t) + B(:, :, :, t), omega_b_full(:, :, :, t), true);
                        %{
                        % interior-only, with zero on all six boundaries
                        omega_QG_interior(:, :, :, t) = SIP_diagnostic_v1_sigma(R, f0, length(lon), length(phi), ...
                                    length(level), phi, lambda, level, ...
                                    sigma_accu(:, :, :, t), dphi, dlambda, rhs(:, :, :, t), Zeros, false);
                        % boundary inversion, with full omega as lateral and lower boundaries
                        omega_QG_full_b(:, :, :, t)   = SIP_diagnostic_v1_sigma(R, f0, length(lon), length(phi), ...
                                    length(level), phi, lambda, level, ...
                                    sigma_accu(:, :, :, t), dphi, dlambda, Zeros, omega_b_full(:, :, :, t), true);
                        % boundary inversion, with full omega as lateral boundaries and zero lower boundary
                        omega_QG_lateral_b(:, :, :, t) = SIP_diagnostic_v1_sigma(R, f0, length(lon), length(phi), ...
                                    length(level), phi, lambda, level, ...
                                    sigma_accu(:, :, :, t), dphi, dlambda, Zeros, omega_b_full(:, :, :, t), false);
                        % boundary inversion, with climatological omega as lateral boundaries and zero lower boundary and no interior
                        omega_QG_lateral_clim_b(:, :, :, t) = SIP_diagnostic_v1_sigma(R, f0, length(lon), length(phi), ...
                                    length(level), phi, lambda, level, ... 
                                    sigma_accu(:, :, :, t), dphi, dlambda, Zeros, omega_b_clim, false);
                        %}
                    end
                    if ACCU_B
                        omega_QG = omega_QG_full;
                    end

                    % write data statistics
                    t0 = (array_t + 1) / 2;
                    aaa = omega_QG(:, :, level == 50000, t0);
                    if ~any(isnan(aaa(:)))
                        
                        % select the center of omega_QG, so that k2 is well-behaved
                        temp_level = 50000;
                        tag_find_center = 1;
                        if FIND_CENTER
                            I = computation_x - 5;
                        	[event_lon, event_lat, tag_find_center] = find_center_omega(event_lon, event_lat, ...
                            	                omega_QG(:, :, level == temp_level, t0), length(lon_series_original), I);
                        elseif FIND_CENTER_ADV
                            I = computation_x - 5;
                            [event_lon, event_lat, tag_find_center] = find_center_omega(event_lon, event_lat, ...
                                -(A(:, :, level == temp_level, t0) + B(:, :, level == temp_level, t0)), length(lon_series_original), I);
                        end
	                    if ~tag_find_center
                            if FIND_CENTER
                                disp('local omega_QG maximum too far from precip maximum, continue');
                            elseif FIND_CENTER_ADV
                                disp('local Adv maximum too far from precip maximum, continue');
                            end
                            num_event_discarded_misalign(jj, ii) = num_event_discarded_misalign(jj, ii) + 1;
                            continue;
                        end
                        indx = find(lon == lon_series_original(event_lon));
                        indy = find(lat == lat_series_original(event_lat));
                        [Ny, Nx] = size(omega_QG(:, :, 1, 1));
                        y_indices = max(indy - (Ny - indy), 1) : min(2 * indy - 1, Ny);
                        x_indices = max(indx - (Nx - indx), 1) : min(2 * indx - 1, Nx);
                        temp_1 = omega_QG(indy, indx, :, t0);
                        
                        % information about the event
                        events_sta{jj, ii}(end + 1) = event;
                        events_sta{jj, ii}(end).event_timespan = event_timespan;
                        events_sta{jj, ii}(end).event_year = event_year;
                        events_sta{jj, ii}(end).event_latspan = event_latspan;
                        events_sta{jj, ii}(end).event_lonspan = event_lonspan;
                        events_sta{jj, ii}(end).computation_x = computation_x;
                        events_sta{jj, ii}(end).computation_y = computation_y;
                        events_sta{jj, ii}(end).event_lat = lat_series_original(lat_indices(jj));
                        events_sta{jj, ii}(end).event_lon = lon_series_original(lon_indices(ii));
                        if FIND_CENTER || FIND_CENTER_ADV
                            events_sta{jj, ii}(end).event_lat_center = lat_series_original(event_lat);
                            events_sta{jj, ii}(end).event_lon_center = lon_series_original(event_lon);
                        end
                        events_sta{jj, ii}(end).event_level = level;
                        events_sta{jj, ii}(end).precip = precip_{jj, ii}(m);
                        events_sta{jj, ii}(end).Adv_500hPa_full = A(y_indices, x_indices, level == 50000, t0) + ...
                                                                  B(y_indices, x_indices, level == 50000, t0);

                        % write in different contributions to omega_QG
                        events_sta{jj, ii}(end).omega_QG                = squeeze(omega_QG(indy, indx, :, t0));
                        events_sta{jj, ii}(end).omega_QG_full           = squeeze(omega_QG_full(indy, indx, :, t0));
                        events_sta{jj, ii}(end).omega_QG_Adv            = squeeze(omega_QG_Adv(indy, indx, :, t0));
                        %{
                        events_sta{jj, ii}(end).omega_QG_interior       = squeeze(omega_QG_interior(indy, indx, :, t0));
                        events_sta{jj, ii}(end).omega_QG_full_b         = squeeze(omega_QG_full_b(indy, indx, :, t0));
                        events_sta{jj, ii}(end).omega_QG_lateral_b      = squeeze(omega_QG_lateral_b(indy, indx, :, t0));
                        events_sta{jj, ii}(end).omega_QG_lateral_clim_b = squeeze(omega_QG_lateral_clim_b(indy, indx, :, t0));
                        %}

                        % write in omega data
                        if DAILY
                            events_sta{jj, ii}(end).omega_500hPa_full = ...
                                    squeeze(mean(omega(y_indices, x_indices, level == 50000, :), 4));
                            events_sta{jj, ii}(end).omega = ...
                                    squeeze(mean(omega(indy, indx, :, :), 4));
                        else
                            events_sta{jj, ii}(end).omega_500hPa_full = omega(y_indices, x_indices, level == 50000, t0);
                            events_sta{jj, ii}(end).omega = squeeze(omega(indy, indx, :, t0));
                        end
                        events_sta{jj, ii}(end).omega_QG_500hPa_full = omega_QG(y_indices, x_indices, level == 50000, t0);
                        events_sta{jj, ii}(end).omega_QG_max = min(temp_1);
                        events_sta{jj, ii}(end).omega_QG_max_p = level(temp_1 == min(temp_1));

                        % write in diagnostics

                        [del_omega_QG, del_J, del_sigma_omega_QG, del_sigma_accu] = deal(zeros(length(level), 1));
                        for k = 1 : length(level)
                            temp_J = spherical_laplacian(J(:, :, k, t0), phi, dphi, dlambda);
                            del_J(k) = temp_J(indy, indx);

                            temp_del_omega_QG = spherical_laplacian(omega_QG(:, :, k, t0), phi, dphi, dlambda);
                            del_omega_QG(k) = temp_del_omega_QG(indy, indx);
                            
                            if ACCU_SIGMA
                                temp_del_sigma_omega_QG = spherical_laplacian...
                                        (omega_QG(:, :, k, t0).*sigma_accu(:, :, k, t0), phi, dphi, dlambda);
                                del_sigma_omega_QG(k) = temp_del_sigma_omega_QG(indy, indx);
                            end
                        end
                        clear('temp_J', 'temp_omega_QG');
                        
                        if ACCU_SIGMA
                            temp_sigma_omega_QG = reshape(...
                                    omega_QG(indy, indx, :, t0).*sigma_accu(indy, indx, :, t0), size(del_sigma_omega_QG));
                            k2 = - del_sigma_omega_QG ./ temp_sigma_omega_QG;
                        else
                            temp_omega_QG = reshape(omega_QG(indy, indx, :, t0), size(del_omega_QG));
                            k2 = - del_omega_QG ./ temp_omega_QG;
                        end
                        
                        % horizontal structure of sigma_star*omega_QG
                        temp_L = 5;
                        %temp_L = min(5, floor(((computation_x + 1)/2 - 1)/2)*2 + 1);
                        dL = (temp_L - 1)/2;
                        [sigma_star, del_sigma_star_omega_QG] = deal(zeros(temp_L, temp_L, length(level)));
                        if DAILY
                            T2 = mean(T, 4);
                        else
                            T2 = T;
                        end
                        for ji = 1 : temp_L^2
                            [j, i] = ind2sub([temp_L, temp_L], ji);
                            temp_T = squeeze(T2(indy - dL - 1 + j, indx - dL - 1 + i, :, t0));
                            [temp, ~, temp_theta] = eff_stat_stab(level, temp_T);
                            temp_dtheta_dp_ma = temp_T ./ temp_theta .* temp;
                            sigma_star(j, i, :) = - temp_dtheta_dp_ma * cp * kappa ./ level;
                        end
                        for k = 1 : length(level)
                            del_sigma_star_omega_QG(:, :, k) = spherical_laplacian...
                                    (omega_QG(indy+[-dL:dL], indx+[-dL:dL], k, t0).*sigma_star(:, :, k), ...
                                    phi(indy+[-dL:dL]), dphi, dlambda);
                        end
                        sigma_star_omega_QG = reshape(...
                                omega_QG(indy+[-dL:dL], indx+[-dL:dL], :, t0).*sigma_star, ...
                                size(del_sigma_star_omega_QG));
                        k2_star = - squeeze(del_sigma_star_omega_QG(dL+1, dL+1, :) ./ ...
                                    sigma_star_omega_QG(dL+1, dL+1, :));

                        % horizontal structure of J
                        temp_J = reshape(J(indy, indx, :, t0), size(del_J));
                        l2 = - del_J ./ temp_J;

                        % vertical structure of omega
                        d2_dp2_omega_QG = squeeze(d2_dp2(omega_QG(indy, indx, :, t0), level));
                        m2 = - squeeze(d2_dp2_omega_QG) ./ squeeze(omega_QG(indy, indx, :, t0));

                        events_sta{jj, ii}(end).k2 = k2;
                        events_sta{jj, ii}(end).k2_star = k2_star;
                        events_sta{jj, ii}(end).sigma_star = squeeze(sigma_star(dL+1, dL+1, :));
                        events_sta{jj, ii}(end).m2 = m2;
                        events_sta{jj, ii}(end).l2 = l2;
                        events_sta{jj, ii}(end).A = squeeze(A(indy, indx, :, t0));
                        events_sta{jj, ii}(end).B = squeeze(B(indy, indx, :, t0));
                        events_sta{jj, ii}(end).C = squeeze(C(indy, indx, :, t0));
                        events_sta{jj, ii}(end).J_center = squeeze(J(indy, indx, :, t0));
                        events_sta{jj, ii}(end).J_500hPa = squeeze(J(y_indices, x_indices, level == 50000, t0));
                        events_sta{jj, ii}(end).sigma = sigma(:, t0);
                        events_sta{jj, ii}(end).sigma_accu = squeeze(sigma_accu(indy, indx, :, t0));
                        events_sta{jj, ii}(end).dtheta_dp_ma_avg = dtheta_dp_ma_avg(:, t0);
                        events_sta{jj, ii}(end).T = squeeze(T2(indy, indx, :, t0)); % for calculation of dtheta_dp_ma*omega
                        events_sta{jj, ii}(end).T_avg = T_avg(:, t0);
                        events_sta{jj, ii}(end).f0 = f0;
                        num_event(jj, ii) = num_event(jj, ii) + 1;


                        % test if events are recorded accurately
                        if rand < 0.01
                            if ACCU_SIGMA
                                temp = (events_sta{jj, ii}(end).k2 .* events_sta{jj, ii}(end).sigma_accu + ...
                                        events_sta{jj, ii}(end).m2 .* f0^2) .* events_sta{jj, ii}(end).omega_QG + ...
                                        events_sta{jj, ii}(end).A + events_sta{jj, ii}(end).B + events_sta{jj, ii}(end).C;
                            else
                                temp = (events_sta{jj, ii}(end).k2 .* events_sta{jj, ii}(end).sigma + ...
                                        events_sta{jj, ii}(end).m2 .* f0^2) .* events_sta{jj, ii}(end).omega_QG + ...
                                        events_sta{jj, ii}(end).A + events_sta{jj, ii}(end).B + events_sta{jj, ii}(end).C;
                            end
                            temp2 = temp - events_sta{jj, ii}(end).C + events_sta{jj, ii}(end).l2 .* kappa ./ level .* events_sta{jj, ii}(end).J_center;
                            disp(['If this number: ', num2str(max(abs([temp(floor(end/2)), temp2(floor(end/2))]))), ' is larger than 1e-30, ', ...
                                  'something is wrong about event recording']);
                        end

                        if NETCDF == true && rand < 1 / NETCDF_interval

                            disp('writing NetCDF file');
                            NetCDF_label = NetCDF_label + 1;

                            % write in QG omega field
                            [omega_QG_in, omega_in] = deal(zeros(length(lon(x_indices)), length(lat(y_indices)), ...
                                                        length(level), array_t));

                            for k = 1 : length(level)
                                for t = 1 : array_t
                                    omega_QG_in(:, :, k, t) = omega_QG(y_indices, x_indices, k, t)';
                                end
                            end

                            nc_filename = [output_path, string_1, '_jj=', num2str(jj, '%10.3d'), '_ii=', num2str(ii, '%10.3d'), ...
                                            '_m=', num2str(m, '%10.3d'), '_omega_QG', '.nc'];
                            if DAILY
                                NetCDF_time = event_timespan(1);
                            else
                                NetCDF_time = event_timespan;
                            end
                            writeNetCDF(nc_filename, 'omega_QG', omega_QG_in, lat(y_indices), lon(x_indices) + 10000, NetCDF_time, level);
                                % 10000 here is to deal with a bug in ncview: it shows fillvalues mysteriously when coordinates are correct

                            % write in event omega field
                            for k = 1 : length(level)
                                if DAILY
                                    omega_in(:, :, k) = mean(omega(y_indices, x_indices, k, :), 4)';
                                else
                                    for t = 1 : array_t
                                        omega_in(:, :, k, t) = omega(y_indices, x_indices, k, t)';
                                    end
                                end
                            end
                            nc_filename = [output_path, string_1, '_jj=', num2str(jj, '%10.3d'), '_ii=', num2str(ii, '%10.3d'), ...
                                           	'_m=', num2str(m, '%10.3d'), '_omega', '.nc'];
                            writeNetCDF(nc_filename, 'omega', omega_in, lat(y_indices), lon(x_indices) + 10000, NetCDF_time, level);

                        end
                    else
                        disp(['NaN found in inverted omega_QG, not saved, m = ', num2str(m)]);
                        continue;
                    end
                end
                num_event_discarded(jj, ii) = length(days_{jj, ii}) - num_event(jj, ii);
            else
                disp(['jj = ', num2str(lat_indices(jj)), ...
                      'ii = ', num2str(lon_indices(ii)), ' statistics does not exist']);
            end
        end
    end

    disp(num_event);
    
    save([output_path, matfilename], 'num_event', 'num_event_discarded', 'num_event_discarded_misalign', ...
            'plevels', 'lon_series', 'lat_series', 'events_sta', '-v7.3');


