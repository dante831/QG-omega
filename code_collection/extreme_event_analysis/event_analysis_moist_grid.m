
% compute event-wise QG_omega for precipitation extremes

function event_analysis_moist_grid(box_min, box_max, p_min, p_max, lat_series_original, lon_series_original, lat_indices, lon_indices, ...
                                    days_, precip_, matfilename, input_path, output_path, NETCDF, NETCDF_interval, n_ensemble, ...
                                    FIND_CENTER, JJA_DJF, ACCU_SIGMA, SIGMA_LOCAL, SIGMA_2, SMOOTH)
    
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

                        event_year = floor((event_day - 0.25) / 365) + 1990; 
                        string_1 = num2str(n_ensemble, '%.3d'); % string for the ensemble number 
    
                        % get whether it is historical or rcp85
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

                        event_year = floor((event_day - 0.25) / 365) + 1980;
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
                    if strfind(pwd_str, 'CESM')
                        [omega, tag] = omega_read_Wrapper(computation_y, computation_x, length(level_indices), ...
                                event_year, input_path, event_latspan, event_lonspan, event_day, level_indices, string_1);
                    elseif strfind(pwd_str, 'GFDL')
                        [omega, tag] = omega_read_Wrapper_GFDL(computation_y, computation_x, length(level_indices), ...
                                event_year, input_path, event_latspan, event_lonspan, event_day, level_indices, string_1);
                    end
                    if ~tag; continue; end;
                    
                    % read in event data
                    if strfind(pwd_str, 'CESM')
                        if strfind(pwd_str, '005')
                            [T, ua, va, ug, vg, omega, omega_b, q, tag] = event_read_Wrapper_uava_q(computation_y, computation_x, length(level_indices), ...
                                    event_year, input_path, event_latspan, event_lonspan, event_day, level_indices, string_1, JJA_DJF);
                        else
                            [T, ua, va, ug, vg, omega, omega_b, tag] = event_read_Wrapper_uava(computation_y, computation_x, length(level_indices), ...
                                    event_year, input_path, event_latspan, event_lonspan, event_day, level_indices, string_1, JJA_DJF);
                        end
                    elseif strfind(pwd_str, 'GFDL')
                        [T, ua, va, ug, vg, omega, omega_b, tag] = event_read_Wrapper_uava_GFDL(computation_y, computation_x, length(level_indices), ...
                                event_year, input_path, event_latspan, event_lonspan, event_day, level_indices, string_1);
                    end
                    if ~tag; continue; end;

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
                    if strfind(pwd_str, '005')
                        q       = q      (y_ind, x_ind, level_indices, :);
                    end
                    omega_b = omega_b(y_ind, x_ind, level_indices, :);
                    event_latspan = event_latspan(y_ind);
                    event_lonspan = event_lonspan(x_ind);
                    
                    % define some other useful arrays
                    lat = lat_series_original(event_latspan);
                    lon = lon_series_original(event_lonspan);
                    phi = lat / 180.0 * 3.1415926;
                    lambda = lon / 180.0 * 3.1415926;
                    event_timespan = event_day + [- 0.25, 0.00, 0.25]; % convert to 'days since 1990'
                    [~, Phi] = meshgrid(lambda, phi);

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
                    [sigma, sigma_eff, T_avg, dtheta_dp_ma_avg] = deal(zeros(length(level), length(event_timespan)));
                    omega_QG = zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan));
                    
                    for t = 1 : length(event_timespan)
                        
                        T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));
                        
                        % effective static stability: sigma_eff
                        temp_omega = squeeze(omega(:, :, :, t));
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
                        [dtheta_dp_ma_avg(:, t), dtheta_dp_eff, theta_avg] = eff_stat_stab(level, T_avg(:, t), lambda_eff);
                        sigma_eff(:, t) = - Ra * T_avg(:, t) ./ (level .* theta_avg) .* dtheta_dp_eff;

                        % static stability: sigma
                        % the following formula is more accurate since it takes into account the moist part in the exponent 
                        % when calculating the potential temperature
                        sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg) .* gradient(theta_avg, level);

                        % spatially varying static stability: sigma_accu
                        % this quantity is mainly used to calculate J term
                        Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
                        theta = T(:, :, :, t) .* (1e5 ./ Level) .^ kappa;
                        sigma_accu(:, :, :, t) = - Ra * T(:, :, :, t) ./ (Level .* theta) .* d_dp(theta, level);
                    
                    end

                    % Q vector
                    [A, B] = Q_vector(level, ug, vg, T, Phi, event_timespan, dphi, dlambda, f0, beta, true);

                    % traditional form
                    if strfind(pwd_str, 'traditional')
                        zeta = zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan));
                        % get relative vorticity
                        for t = 1 : length(event_timespan)
                            for k = 1 : length(level)
                                %zeta(:, :, k, t) = curl(phi, ua(:, :, k, t), va(:, :, k, t), dphi, dlambda);
                                zeta(:, :, k, t) = curl(phi, ug(:, :, k, t), vg(:, :, k, t), dphi, dlambda);
                            end
                        end
                        % term A
                        [A, ~, ~, ~] = traditional_A(phi, level, ug, vg, zeta, event_timespan, dphi, dlambda, f0, beta);
                        % term B
                        B = traditional_B(phi, level, ug, vg, T, event_timespan, dphi, dlambda);
                    end

                    % term C (diabatic heating)
                    [C1, C2, J1, J2] = diabatic_heating(phi, T, level, ua, va, sigma_accu, omega, ...
                                                event_timespan, dphi, dlambda, dt);
                    J = J1 + J2;
                    C = C1 + C2;
                        
                    % right hand side
                    rhs = A + B + C;
                    
                    % remove negative values in sigma_accu
                    [sigma_accu, tag_sigma] = sigma_remove_negative(sigma_accu, sigma);
                    if ~tag_sigma; disp('too many points with negative stability!');continue; end
       
                    % smooth sigma_accu field
                    sigma_accu_smoothed = sigma_accu;
                    for k = 1 : length(level)
                        for t = 1 : length(event_timespan)
                            if SMOOTH
                                sigma_accu_smoothed(:, :, k, t) = smooth2a(sigma_accu(:, :, k, t), window_x, window_y);
                                sigma_accu_smoothed(:, :, k, t) = smooth_normalize(sigma_accu(:, :, k, t), sigma_accu_smoothed(:, :, k, t));
                            end
                        end
                    end
                    sigma_accu = sigma_accu_smoothed;
                    clear('sigma_accu_smoothed');

                    if exist('SIGMA_2') && SIGMA_2
                        C_temp = C(:, :, level == 50000, 2);
                        AB_temp = A(:, :, level == 50000, 2) + B(:, :, level == 50000, 2);
                        if mean(AB_temp(:)) < 0.5*mean(C_temp(:)) || mean(AB_temp(:)) < 0
                            disp('QG forcing too weak, abandoned')
                            continue;
                        end
                    end

                    % do the inversion
                    for t = 1 : length(event_timespan)
                        if ACCU_SIGMA
                            omega_QG(:, :, :, t) = SIP_diagnostic_v1_sigma(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
                                                                    sigma_accu(:, :, :, t), dphi, dlambda, rhs(:, :, :, t), omega_b, false);
                        elseif exist('SIGMA_LOCAL') && SIGMA_LOCAL % Paul's idea, use local sigma_accu as zonal averaged sigma
                            sigma(:, t) = squeeze(sigma_accu(ceil(end/2), ceil(end/2), :, t));
                            omega_QG(:, :, :, t) = SIP_diagnostic_v1(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
                                                                    sigma(:, t), dphi, dlambda, rhs(:, :, :, t), omega_b, false);
                        elseif exist('SIGMA_2') && SIGMA_2
                            tag_sigma_small = false;
                            [SSigma_t, SSigma_t_0] = deal(repmat(reshape(sigma(:, t), [1, 1, length(level)]), [length(phi), length(lon), 1]));
                            Plevels = repmat(reshape(level, [1, 1, length(level)]), [length(phi), length(lon), 1]);
                            Dtheta_dp_ma_avg = repmat(reshape(dtheta_dp_ma_avg(:, t), [1, 1, length(level)]), [length(phi), length(lon), 1]);
                            Sigma_m = max(SSigma_t_0 + Ra ./ Plevels .* Dtheta_dp_ma_avg, 0.3 * SSigma_t_0);

                            SSigma_t = Sigma_m    .* double(omega(:, :, :, t) < 0) + ...
                                       SSigma_t_0 .* double(omega(:, :, :, t) >= 0);
                            if tag_sigma_small
                                disp('warning: some sigma_m are set to reduced sigma*0.2')
                            end
                            rhs(:, :, :, t) = A(:, :, :, t) + B(:, :, :, t);
                            r = 1;
                            temp_omega_QG_2 = omega(:, :, :, t);
                            temp_omega_QG_2 = zeros(size(omega(:, :, :, t)));
                            iterations = 0;
                            temp_omega_QG = SIP_diagnostic_v1_sigma(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
                                                            SSigma_t, dphi, dlambda, rhs(:, :, :, t),  omega_b, false);
                            omega_QG(:, :, :, t) = temp_omega_QG;
                            SSigma(:, :, :, t) = SSigma_t;
                        else
                            omega_QG(:, :, :, t) = SIP_diagnostic_v1(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
                                                                    sigma(:, t), dphi, dlambda, rhs(:, :, :, t),  omega_b, false);
                        end
                    end
                    
                    % write data statistics
                    aaa = omega_QG(:, :, level == 50000, 2);
                    if ~any(isnan(aaa(:)))

                        % select the center of omega_QG, so that k2 is well-behaved
                        temp_level = 50000;
                        if FIND_CENTER
                        	[event_lon, event_lat, tag] = find_center_omega(event_lon, event_lat, ...
                            	                            omega_QG(:, :, level == temp_level, 2), length(lon_series_original));
                        end
                        if ~tag; 
                            disp('local omega_QG maximum too far from precip maximum, continue');
                            num_event_discarded_misalign(jj, ii) = num_event_discarded_misalign(jj, ii) + 1;
                            continue; 
                        end;
                        
                        indx = find(lon == lon_series_original(event_lon));
                        indy = find(lat == lat_series_original(event_lat));
                        [Ny, Nx] = size(omega_QG(:, :, 1, 1));
                        y_indices = max(indy - (Ny - indy), 1) : min(2 * indy - 1, Ny);
                        x_indices = max(indx - (Nx - indx), 1) : min(2 * indx - 1, Nx);
                        temp_1 = omega_QG(indy, indx, :, 2);
                        
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
                        events_sta{jj, ii}(end).Adv_500hPa_full = A(y_indices, x_indices, level == 50000, 2) + ...
                                                                  B(y_indices, x_indices, level == 50000, 2);

                        % write in omega data
                        events_sta{jj, ii}(end).omega = squeeze(omega(indy, indx, :, 2));
                        events_sta{jj, ii}(end).omega_QG = squeeze(omega_QG(indy, indx, :, 2));
                        events_sta{jj, ii}(end).omega_500hPa_full = omega(y_indices, x_indices, level == 50000, 2);
                        %events_sta{jj, ii}(end).omega_QG_400hPa_full = omega_QG(y_indices, x_indices, level == 40000, 2);
                        events_sta{jj, ii}(end).omega_QG_500hPa_full = omega_QG(y_indices, x_indices, level == 50000, 2);
                        %events_sta{jj, ii}(end).omega_QG_600hPa_full = omega_QG(y_indices, x_indices, level == 60000, 2);
                        %events_sta{jj, ii}(end).omega_QG_max = min(temp_1);
                        %events_sta{jj, ii}(end).omega_QG_max_p = level(temp_1 == min(temp_1));

                        % write in diagnostics
                        [del_omega_QG, del_J, del_sigma_omega_QG, del_sigma_accu] = deal(zeros(length(level), 1));
                        for k = 1 : length(level)
                            temp_J = spherical_laplacian(J(:, :, k, 2), phi, dphi, dlambda);
                            del_J(k) = temp_J(indy, indx);
                            
                            temp_del_omega_QG = spherical_laplacian(omega_QG(:, :, k, 2), phi, dphi, dlambda);
                            del_omega_QG(k) = temp_del_omega_QG(indy, indx);

                            if ACCU_SIGMA
                                temp_del_sigma_omega_QG = spherical_laplacian...
                                        (omega_QG(:, :, k, 2).*sigma_accu(:, :, k, 2), phi, dphi, dlambda);
                                del_sigma_omega_QG(k) = temp_del_sigma_omega_QG(indy, indx);
                            elseif SIGMA_2
                                temp_del_sigma_omega_QG = spherical_laplacian...
                                        (omega_QG(:, :, k, 2).*SSigma(:, :, k, 2), phi, dphi, dlambda);
                                del_sigma_omega_QG(k) = temp_del_sigma_omega_QG(indy, indx);
                            end

                        end
                        clear('temp_J', 'temp_omega_QG');
                        if ACCU_SIGMA
                            temp_sigma_omega_QG = reshape(...
                                    omega_QG(indy, indx, :, 2).*sigma_accu(indy, indx, :, 2), size(del_sigma_omega_QG));
                            k2 = - del_sigma_omega_QG ./ temp_sigma_omega_QG;
                        elseif SIGMA_2
                            temp_sigma_omega_QG = reshape(...
                                    omega_QG(indy, indx, :, 2).*SSigma(indy, indx, :, 2), size(del_sigma_omega_QG));
                            k2 = - del_sigma_omega_QG ./ temp_sigma_omega_QG;
                        else
                            temp_omega_QG = reshape(omega_QG(indy, indx, :, 2), size(del_omega_QG));
                            k2 = - del_omega_QG ./ temp_omega_QG;
                        end

                        temp_J = reshape(J(indy, indx, :, 2), size(del_J));
                        l2 = - del_J ./ temp_J;

                        d2_dp2_omega_QG = squeeze(d2_dp2(omega_QG(indy, indx, :, 2), level));
                        m2 = - squeeze(d2_dp2_omega_QG) ./ squeeze(omega_QG(indy, indx, :, 2));

                        events_sta{jj, ii}(end).k2 = k2;
                        events_sta{jj, ii}(end).m2 = m2;
                        events_sta{jj, ii}(end).l2 = l2;
                        events_sta{jj, ii}(end).A = squeeze(A(indy, indx, :, 2));
                        events_sta{jj, ii}(end).B = squeeze(B(indy, indx, :, 2));
                        events_sta{jj, ii}(end).C = squeeze(C(indy, indx, :, 2));
                        events_sta{jj, ii}(end).J_center = squeeze(J(indy, indx, :, 2));
                        %events_sta{jj, ii}(end).J_400hPa = squeeze(J(y_indices, x_indices, level == 40000, 2));
                        events_sta{jj, ii}(end).J_500hPa = squeeze(J(y_indices, x_indices, level == 50000, 2));
                        %events_sta{jj, ii}(end).J_600hPa = squeeze(J(y_indices, x_indices, level == 60000, 2));
                        events_sta{jj, ii}(end).sigma = sigma(:, 2);
                        events_sta{jj, ii}(end).sigma_eff = sigma_eff(:, 2);
                        if exist('SIGMA_2') && SIGMA_2
                            events_sta{jj, ii}(end).sigma_accu = squeeze(SSigma(indy, indx, :, 2));
                        else
                            events_sta{jj, ii}(end).sigma_accu = squeeze(sigma_accu(indy, indx, :, 2));
                        end
                        events_sta{jj, ii}(end).dtheta_dp_ma_avg = dtheta_dp_ma_avg(:, 2);
                        events_sta{jj, ii}(end).T = squeeze(T(indy, indx, :, 2));
                        events_sta{jj, ii}(end).T_avg = T_avg(:, 2);
                        events_sta{jj, ii}(end).f0 = f0;
                        if strfind(pwd_str, '005')
                            events_sta{jj, ii}(end).q = squeeze(q(indy, indx, :, 2));
                        end
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
                            [omega_QG_in, omega_in] = deal(zeros(length(lon(x_indices)), length(lat(y_indices)), ...
                                        length(level), length(event_timespan)));

                            for k = 1 : length(level)
                                for t = 1 : length(event_timespan)
                                    omega_QG_in(:, :, k, t) = omega_QG(y_indices, x_indices, k, t)';
                                end
                            end
                            nc_filename = [output_path, string_1, '_jj=', num2str(jj, '%10.3d'), '_ii=', num2str(ii, '%10.3d'), ...
                                            '_m=', num2str(m, '%10.3d'), '_omega_QG', '.nc'];
                            writeNetCDF(nc_filename, 'omega_QG', omega_QG_in, lat(y_indices), lon(x_indices) + 10000, event_timespan, level);
                                % 10000 here is to deal with a bug in ncview: it shows fillvalues mysteriously when coordinates are correct

                            % write in event omega field
                            for k = 1 : length(level)
                                for t = 1 : length(event_timespan)
                                    omega_in(:, :, k, t) = omega(y_indices, x_indices, k, t)';
                                end
                            end
                            nc_filename = [output_path, string_1, '_jj=', num2str(jj, '%10.3d'), '_ii=', num2str(ii, '%10.3d'), ...
                                            '_m=', num2str(m, '%10.3d'), '_omega', '.nc'];
                            writeNetCDF(nc_filename, 'omega', omega_in, lat(y_indices), lon(x_indices) + 10000, event_timespan, level);

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




