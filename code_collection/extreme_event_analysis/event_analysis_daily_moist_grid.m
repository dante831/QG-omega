
% compute event-wise QG_omega for precipitation extremes
% final version, edited on Feb. 6th, 2019 by Ziwei Li

%{
precip_ = event_precip_historical;
days_ = days_historical;
jj = 50;
ii = 120;
m = 1;
%}

function event_analysis_daily_moist_grid(box_min, box_max, p_min, p_max, lat_series_original, lon_series_original, lat_indices, lon_indices, ...
                                    days_, precip_, matfilename, input_path, output_path, NETCDF, NETCDF_interval, n_ensemble, ...
                                    FIND_CENTER, JJA_DJF, ACCU_SIGMA, SIGMA_LOCAL, SIGMA_2)
    
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
    
    [num_event, num_event_discarded, num_event_discarded_misalign] = deal(zeros(length(lat_indices), length(lon_indices))); % define array for number of events
    events_sta = deal(cell(length(lat_indices), length(lon_indices))); % define cell array for "event statistics"
    pwd_str = pwd; % get current path

    
    for jj = 1 : length(lat_indices)
        for ii = 1 : length(lon_indices)
            disp(['lat: ', num2str(lat_series(jj)), 'lon: ', num2str(lon_series(ii))]);
            if ~isempty(days_{jj, ii})

                for m = 1 : length(days_{jj, ii})
                    disp(['m = ', num2str(m)])
                    event_day = days_{jj, ii}(m); 
                    event_lon = lon_indices(ii);
                    event_lat = lat_indices(jj);
                    
                    if strfind(pwd_str, 'CESM')

                        %[computation_x, computation_y] = deal(box_max);
                        box_max_x = box_max;
                        box_max_y = box_max;
                        computation_y = box_max_y;
                        computation_x = box_max_x;

                        %event_year = floor((event_day - 0.25) / 365) + 1990; % 6-hourly
                        event_year = floor((event_day - 1.0) / 365) + 1990; % daily
                        string_1 = num2str(n_ensemble, '%.3d'); % string for the ensemble number 

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

                        box_max_x = box_max;
                        box_max_y = box_max;
                        computation_y = box_max_y;
                        computation_x = box_max_x;

                        %event_year = floor((event_day - 0.25) / 365) + 1980; % 6-hourly
                        event_year = floor((event_day - 1.0) / 365) + 1980; % daily
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
                    
                    % compute indices to invert omega equation in
                    % lon_series_original and lat_series_original

                    plevels = double(ncread([input_path, 'ta_', string_1, '_', num2str(event_year), '.nc'], 'level'));
                    p_start_0 = 1; % lowest level
                    p_start = find(plevels == p_max); % lower boundary start point
                    p_end = find(plevels == p_min); % highest level
                    level_indices = p_start_0 : p_end;

                    % read in omega field
                    % daily:
                    event_timespan = event_day - 1 + [0.25, 0.50, 0.75, 1.00];
                    time_indices = mod(event_timespan - 0.25,  365) / 0.25 + 1;
                    if strfind(pwd_str, 'CESM')
                        omega = omega_read_NetCDF(event_year, input_path, event_latspan, event_lonspan, ...
                                                    time_indices, level_indices, string_1);
                    elseif strfind(pwd_str, 'GFDL')
                        omega = omega_read_NetCDF_GFDL(event_year, input_path, event_latspan, event_lonspan, ...
                                                    time_indices, level_indices);
                    end
                    
                    % read in event data
                    % daily:
                    if strfind(pwd_str, 'CESM')
                        [T, ua, va, ug, vg, omega, omega_b] = event_read_NetCDF_with_b_uava...
                                (event_year, input_path, event_latspan, event_lonspan, time_indices, level_indices, string_1, JJA_DJF);
                    elseif strfind(pwd_str, 'GFDL')
                        [T, ua, va, ug, vg, omega, omega_b] = event_read_NetCDF_with_b_uava_GFDL...
                                (event_year, input_path, event_latspan, event_lonspan, time_indices, level_indices);
                    end

                    % get rid of events which have NaN values between p_start
                    % and p_end, then expand the selected data to the 
                    % lowest level where it has no NaN values

                    % add shrink_box, if range as large as box_max has NaN value, then shrink until box_min
                    [x_ind, y_ind, computation_x, computation_y, p_start, NaN_tag] = ...
                            shrink_box_v1(omega, box_min, p_max, computation_x, computation_y, p_start, p_end, plevels);
                    if NaN_tag; disp('NaN found in omega, continue');continue; end;

                    % subselect all fields and indices according to the result from shrink_box
                    level_indices = p_start : p_end;
                    level = plevels(level_indices);
                    T       = T      (y_ind, x_ind, level_indices, :);
                    ug      = ug     (y_ind, x_ind, level_indices, :);
                    vg      = vg     (y_ind, x_ind, level_indices, :);
                    ua      = ua     (y_ind, x_ind, level_indices, :);
                    va      = va     (y_ind, x_ind, level_indices, :);
                    omega   = omega  (y_ind, x_ind, level_indices, :);
                    omega_b = omega_b(y_ind, x_ind, level_indices, :);
                    event_latspan = event_latspan(y_ind);
                    event_lonspan = event_lonspan(x_ind);

                    % define some other useful arrays
                    lat = lat_series_original(event_latspan);
                    lon = lon_series_original(event_lonspan);
                    phi = lat / 180.0 * 3.1415926;
                    lambda = lon / 180.0 * 3.1415926;
                    [~, Phi] = meshgrid(lambda, phi);

                    % smooth T field
                    T_smoothed = T;
                    for k = 1 : length(level)
                        for t = 1 : length(event_timespan)
                            T_smoothed(:, :, k, t) = smooth2a(T(:, :, k, t), window_x, window_y);
                            % renormalize so that the smoothed field has the same variance as the previous field
                            T_smoothed(:, :, k, t) = smooth_normalize(T(:, :, k, t), T_smoothed(:, :, k, t));
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
                    [sigma, sigma_eff, dtheta_dp_ma, T_avg, dtheta_dp_ma_avg] = deal(zeros(length(level)));
                    omega_QG = zeros(length(event_latspan), length(event_lonspan), length(level));
                    
                    % zonal- and time-averaged temperature    
                    T_avg = reshape(mean(mean(mean(T, 1), 2), 4), size(level));
                        
                    % effective static stability: sigma_eff
                    temp_omega = squeeze(squeeze(mean(omega, 4)));
                    omega_prime = temp_omega - repmat(mean(mean(temp_omega, 1), 2), ...
                            [size(temp_omega(:, :, 1)), 1]);
                    omega_up = temp_omega; omega_up(temp_omega > 0) = 0;
                    omega_up_prime = omega_up - repmat(mean(mean(omega_up, 1), 2), ...
                            [size(omega_up(:, :, 1)), 1]);
                    % equation (5) in O'Gorman, 2011
                    lambda_eff = squeeze(mean(mean(omega_up .* omega_up_prime, 1), 2) ./ ...
                                     mean(mean(omega_up.^2, 1), 2));
                    Lambda_eff = repmat(reshape(lambda_eff, [1, 1, length(lambda_eff)]), [size(omega_up(:, :, 1)), 1]);
                    epsilon = omega_up_prime - Lambda_eff .* omega_prime;
                    e2 = mean(mean(epsilon.^2)) ./ mean(mean(omega_up_prime.^2)); % this is fairly small between 925hPa and 300hPa
                    [dtheta_dp_ma_avg, dtheta_dp_eff, theta_avg] = eff_stat_stab(level, T_avg, lambda_eff);
                    sigma_eff = - Ra * T_avg ./ (level .* theta_avg) .* dtheta_dp_eff;
                    
                    % static stability: sigma
                    sigma =  - Ra * T_avg ./ (level .* theta_avg) .* gradient(theta_avg, level);

                    % Q vector
                    [A, B] = Q_vector_v1(level, ug, vg, T, Phi, event_timespan, dphi, dlambda, f0, beta, true);
                    A = mean(A, 4);
                    B = mean(B, 4);
                    
                    % term C (latent heating)
                    
                    [~, ~, J1, ~] = diabatic_heating(phi, T, level, ua, va, sigma_accu, omega, ...
                                event_timespan, dphi, dlambda, dt);
                                            
                    % spatially varying static stability: sigma_accu
                    % this quantity is mainly used to calculate J term            
                    Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
                    theta = mean(T, 4) .* (1e5 ./ Level) .^ kappa;
                    sigma_accu = - Ra * mean(T, 4) ./ (Level .* theta) .* d_dp(theta, level);
                   
                    % calculate terms in the right-hand-side
                    J2 = - sigma_accu .* mean(omega, 4) / Ra .* Level * cp;
                    J = mean(J1, 4) + J2;
                    C = deal(zeros(length(event_latspan), length(event_lonspan), length(level)));
                    for k = 1 : length(level)
                        C(:, :, k) = - kappa / level(k) * spherical_laplacian(J(:, :, k), phi, dphi, dlambda);
                    end

                    % right hand side
                    rhs = A + B + C;
                    
                    % remove negative values in sigma_accu
                    [sigma_accu, tag_sigma] = sigma_remove_negative(sigma_accu, sigma);
                    if ~tag_sigma; disp('too many points with negative stability!');continue; end

                    % smooth sigma_accu field
                    sigma_accu_smoothed = sigma_accu;
                    for k = 1 : length(level)
                        sigma_accu_smoothed(:, :, k) = smooth2a(sigma_accu(:, :, k), 1, 1);
                        sigma_accu_smoothed(:, :, k) = smooth_normalize(sigma_accu(:, :, k), sigma_accu_smoothed(:, :, k));
                    end
                    sigma_accu = sigma_accu_smoothed;
                    clear('sigma_accu_smoothed');

                    % do the inversion
                    if ACCU_SIGMA
                        omega_QG(:, :, :, 1) = SIP_diagnostic_v1_sigma(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
                                                                sigma_accu, dphi, dlambda, rhs, omega_b, false);
                    elseif exist('SIGMA_LOCAL') && SIGMA_LOCAL % Paul's idea, use local sigma_accu as zonal averaged sigma
                        sigma = squeeze(sigma_accu(ceil(end/2), ceil(end/2), :, 1));
                        omega_QG(:, :, :, 1) = SIP_diagnostic_v1(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
                                                                sigma, dphi, dlambda, rhs, omega_b, false);
                    else
                        omega_QG(:, :, :, 1) = SIP_diagnostic_v1(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
                                                                sigma, dphi, dlambda, rhs, omega_b, false);
                    end

                    % write data statistics
                    aaa = omega_QG(:, :, level == 50000);
                    if ~any(isnan(aaa(:)))
                        
                        % select the center of omega_QG, so that k2 is well-behaved
                        temp_level = 50000;
                        if FIND_CENTER
                        	[event_lon, event_lat, tag] = find_center_omega(event_lon, event_lat, ...
                            	                            omega_QG(:, :, level == temp_level), length(lon_series_original));
                        end
	                    if ~tag
                            disp('local omega_QG maximum too far from precip maximum, continue');
                            num_event_discarded_misalign(jj, ii) = num_event_discarded_misalign(jj, ii) + 1;
                            continue;
                        end
                        indx = find(lon == lon_series_original(event_lon));
                        indy = find(lat == lat_series_original(event_lat));
                        [Ny, Nx] = size(omega_QG(:, :, 1, 1));
                        y_indices = max(indy - (Ny - indy), 1) : min(2 * indy - 1, Ny);
                        x_indices = max(indx - (Nx - indx), 1) : min(2 * indx - 1, Nx);
                        temp_1 = omega_QG(indy, indx, :);
                        
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
                        if FIND_CENTER
                            events_sta{jj, ii}(end).event_lat_center = lat_series_original(event_lat);
                            events_sta{jj, ii}(end).event_lon_center = lon_series_original(event_lon);
                        end
                        events_sta{jj, ii}(end).event_level = level;
                        events_sta{jj, ii}(end).precip = precip_{jj, ii}(m);
                        events_sta{jj, ii}(end).Adv_500hPa_full = A(y_indices, x_indices, level == 50000) + ...
                                                                  B(y_indices, x_indices, level == 50000);

                        % write in omega data

                        events_sta{jj, ii}(end).omega = squeeze(mean(omega(indy, indx, :, :), 4));
                        events_sta{jj, ii}(end).omega_QG = squeeze(omega_QG(indy, indx, :));
                        events_sta{jj, ii}(end).omega_500hPa_full = squeeze(mean(omega(y_indices, x_indices, level == 50000, :), 4));
                        %events_sta{jj, ii}(end).omega_QG_400hPa_full = omega_QG(y_indices, x_indices, level == 40000, :);
                        events_sta{jj, ii}(end).omega_QG_500hPa_full = omega_QG(y_indices, x_indices, level == 50000, :);
                        %events_sta{jj, ii}(end).omega_QG_600hPa_full = omega_QG(y_indices, x_indices, level == 60000, :);
                        %events_sta{jj, ii}(end).omega_QG_max = min(temp_1);
                        %events_sta{jj, ii}(end).omega_QG_max_p = level(temp_1 == min(temp_1));

                        % write in diagnostics

                        [del_omega_QG, del_J, del_sigma_omega_QG, del_sigma_accu] = deal(zeros(length(level), 1));
                        for k = 1 : length(level)
                            temp_J = spherical_laplacian(J(:, :, k), phi, dphi, dlambda);
                            del_J(k) = temp_J(indy, indx);

                            temp_del_omega_QG = spherical_laplacian(omega_QG(:, :, k), phi, dphi, dlambda);
                            del_omega_QG(k) = temp_del_omega_QG(indy, indx);
                            
                            if ACCU_SIGMA
                                temp_del_sigma_omega_QG = spherical_laplacian...
                                        (omega_QG(:, :, k).*sigma_accu(:, :, k), phi, dphi, dlambda);
                                del_sigma_omega_QG(k) = temp_del_sigma_omega_QG(indy, indx);
                            elseif SIGMA_2
                                temp_del_sigma_omega_QG = spherical_laplacian...
                                        (omega_QG(:, :, k).*SSigma(:, :, k), phi, dphi, dlambda);
                                del_sigma_omega_QG(k) = temp_del_sigma_omega_QG(indy, indx);
                            end

                        end
                        clear('temp_J', 'temp_omega_QG');
                        if ACCU_SIGMA
                            temp_sigma_omega_QG = reshape(...
                                    omega_QG(indy, indx, :).*sigma_accu(indy, indx, :), size(del_sigma_omega_QG));
                            k2 = - del_sigma_omega_QG ./ temp_sigma_omega_QG;
                        elseif SIGMA_2
                            temp_sigma_omega_QG = reshape(...
                                    omega_QG(indy, indx, :).*SSigma(indy, indx, :), size(del_sigma_omega_QG));
                            k2 = - del_sigma_omega_QG ./ temp_sigma_omega_QG;
                        else
                            temp_omega_QG = reshape(omega_QG(indy, indx, :), size(del_omega_QG));
                            k2 = - del_omega_QG ./ temp_omega_QG;
                        end
                        
                        temp_J = reshape(J(indy, indx, :), size(del_J));
                        l2 = - del_J ./ temp_J;

                        d2_dp2_omega_QG = squeeze(d2_dp2(omega_QG(indy, indx, :), level));
                        m2 = - squeeze(d2_dp2_omega_QG) ./ squeeze(omega_QG(indy, indx, :));

                        events_sta{jj, ii}(end).k2 = k2;
                        events_sta{jj, ii}(end).m2 = m2;
                        events_sta{jj, ii}(end).l2 = l2;
                        events_sta{jj, ii}(end).A = squeeze(A(indy, indx, :));
                        events_sta{jj, ii}(end).B = squeeze(B(indy, indx, :));
                        events_sta{jj, ii}(end).C = squeeze(C(indy, indx, :));
                        events_sta{jj, ii}(end).J_center = squeeze(J(indy, indx, :));
                        %events_sta{jj, ii}(end).J_400hPa = squeeze(J(y_indices, x_indices, level == 40000));
                        events_sta{jj, ii}(end).J_500hPa = squeeze(J(y_indices, x_indices, level == 50000));
                        %events_sta{jj, ii}(end).J_600hPa = squeeze(J(y_indices, x_indices, level == 60000));
                        events_sta{jj, ii}(end).sigma = sigma;
                        events_sta{jj, ii}(end).sigma_eff = sigma_eff;
                        if exist('SIGMA_2') && SIGMA_2
                            events_sta{jj, ii}(end).sigma_accu = squeeze(SSigma(indy, indx, :));
                        else
                            events_sta{jj, ii}(end).sigma_accu = squeeze(sigma_accu(indy, indx, :));
                        end
                        events_sta{jj, ii}(end).dtheta_dp_ma_avg = dtheta_dp_ma_avg;
                        events_sta{jj, ii}(end).T = squeeze(mean(T(indy, indx, :, :), 4)); % for calculation of dtheta_dp_ma*omega
                        events_sta{jj, ii}(end).T_avg = T_avg;
                        events_sta{jj, ii}(end).f0 = f0;
                        num_event(jj, ii) = num_event(jj, ii) + 1;

                        % update pointer

                        %%%%%%%%%%%%%% DON'T DELETE %%%%%%%%%%%%%%%%%%
                        % test if events are recorded accurately
                        if rand < 0.01
                            if exist('SIGMA_2') && SIGMA_2
                                temp = (events_sta{jj, ii}(end).k2 .* events_sta{jj, ii}(end).sigma_accu + ...
                                        events_sta{jj, ii}(end).m2 .* f0^2) .* events_sta{jj, ii}(end).omega_QG + ...
                                        events_sta{jj, ii}(end).A + events_sta{jj, ii}(end).B;
                            elseif ACCU_SIGMA
                                temp = (events_sta{jj, ii}(end).k2 .* events_sta{jj, ii}(end).sigma_accu + ...
                                        events_sta{jj, ii}(end).m2 .* f0^2) .* events_sta{jj, ii}(end).omega_QG + ...
                                        events_sta{jj, ii}(end).A + events_sta{jj, ii}(end).B + events_sta{jj, ii}(end).C;
                            else
                                temp = (events_sta{jj, ii}(end).k2 .* events_sta{jj, ii}(end).sigma + ...
                                        events_sta{jj, ii}(end).m2 .* f0^2) .* events_sta{jj, ii}(end).omega_QG + ...
                                        events_sta{jj, ii}(end).A + events_sta{jj, ii}(end).B + events_sta{jj, ii}(end).C;
                            end
                            disp(['If this number: ', num2str(temp(floor(end/2))), ' is larger than 1e-30, ', ...
                                  'something is wrong about event recording']);
                        end

                        if exist('SIGMA_2') && SIGMA_2
                            disp('event saved')
                        end
                    
                        if NETCDF == true && rand < 1 / NETCDF_interval

                            disp('writing NetCDF file');
                            NetCDF_label = NetCDF_label + 1;

                            % write in QG omega field
                            [omega_QG_in, omega_in] = deal(zeros(length(lon(x_indices)), length(lat(y_indices)), length(level)));

                            for k = 1 : length(level)
                                omega_QG_in(:, :, k) = omega_QG(y_indices, x_indices, k)';
                            end

                            nc_filename = [output_path, string_1, '_jj=', num2str(jj, '%10.3d'), '_ii=', num2str(ii, '%10.3d'), ...
                                            '_m=', num2str(m, '%10.3d'), '_omega_QG', '.nc'];
                            writeNetCDF(nc_filename, 'omega_QG', omega_QG_in, lat(y_indices), lon(x_indices) + 10000, event_timespan(1), level);
                                % 10000 here is to deal with a bug in ncview: it shows fillvalues mysteriously when coordinates are correct

                            % write in event omega field
                            for k = 1 : length(level)
                                omega_in(:, :, k) = mean(omega(y_indices, x_indices, k, :), 4)';
                            end
                            nc_filename = [output_path, string_1, '_jj=', num2str(jj, '%10.3d'), '_ii=', num2str(ii, '%10.3d'), ...
                                           	'_m=', num2str(m, '%10.3d'), '_omega', '.nc'];
                            writeNetCDF(nc_filename, 'omega', omega_in, lat(y_indices), lon(x_indices) + 10000, event_timespan(1), level);

                        end
                    else
                        disp(['NaN found in omega_QG, not saved, m = ', num2str(m)]);
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




