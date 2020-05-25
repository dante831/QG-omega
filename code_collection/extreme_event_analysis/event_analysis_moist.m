% compute event-wise QG_omega for precipitation extremes

function event_analysis_moist(box_min, box_max, p_min, p_max, lat_series_original, lon_series_original, lat_indices, lon_indices, ...
                lons_, days_, precip_, matfilename, input_path, output_path, NETCDF, NETCDF_interval, sigma_tag, sigma_smooth_window, n_ensemble)
    
    Omega = 7.2921e-5;
    R = 6371000.0;
    Ra = 287.04;
    cp = 1005.0;
    kappa = Ra/cp;
    dt = 3600.0 * 6.0;
    window_x = 1;
    window_y = 1;
            
    lon_series = lon_series_original(lon_indices);
    lat_series = lat_series_original(lat_indices);
    dphi = (lat_series(2) - lat_series(1)) / 180.0 * 3.1415926;
    dlambda = (lon_series(2) - lon_series(1)) / 180.0 * 3.1415926;
    
    num_event = deal(zeros(length(lat_indices), 1));
    events_sta = cell(length(lat_indices), 1);

    
    for jj = 1 : length(lat_indices)
        disp(['lat: ', num2str(lat_series(jj))]);
        if ~isempty(lons_{jj})
            
            for m = 1 : length(lons_{jj})
                
                event_lon = find(lon_series_original == lon_series(lons_{jj}(m))); % this is event longitude index in lon_series_original
                event_lat = lat_indices(jj);
                event_day = days_{jj}(m); 
                
                %[computation_x, computation_y] = deal(box_max);
                box_max_x = 29;
                box_max_y = 27;
                computation_y = box_max_y;
                computation_x = box_max_x;

                event_year = floor((event_day - 0.25) / 365) + 1990;
                string_1 = num2str(n_ensemble, '%.3d');
                if event_year >= 1990 && event_year <= 2005
                    string_2 = 'historical';
                elseif event_year >= 2071 && event_year <= 2080
                    string_2 = 'rcp85';
                else
                    disp(['year ', num2str(event_year), ' does not exist in the data, abandoned']);
                    continue;
                end
                            
  				[event_latspan, event_lonspan, center_y, computation_y] = ...
                        event_span_Wrapper(lat_series_original, lon_series_original, event_lat, event_lon, ...
                        box_max_y, computation_x, computation_y);
                
                % compute indices to invert omega equation in
                % lon_series_original and lat_series_original
                
                plevels = double(ncread([input_path, 'ta_', string_1, '_', num2str(event_year), '.nc'], 'level'));
                p_start_0 = 1; % lowest level
                p_start = find(plevels == p_max);
                p_end = find(plevels == p_min); % highest level
                level_indices = p_start_0 : p_end;
                
                % find center of omega at 500hPa
                [~, ~, ~, omega, ~, ~] = event_read_Wrapper(computation_y, computation_x, length(level_indices), ...
                        event_year, input_path, event_latspan, event_lonspan, event_day, level_indices, string_1);
                temp_level = 50000;
                [event_lon, event_lat, tag] = find_center_omega(event_lon, event_lat, ...
                                                omega(:, :, level_indices(plevels == temp_level), 2), length(lon_series_original));
                if ~tag
                    disp('local omega maxima too far from precip maxima, continue');
                    continue;
                end
                
                [event_latspan, event_lonspan, center_y, computation_y] = ...
                        event_span_Wrapper(lat_series_original, lon_series_original, event_lat, event_lon, ...
                        box_max_y, computation_x, computation_y);

                % read in data
                
                [T, ug, vg, omega, omega_b, tag] = event_read_Wrapper(computation_y, computation_x, length(level_indices), ...
                        event_year, input_path, event_latspan, event_lonspan, event_day, level_indices, string_1);
                if ~tag; continue; end;
                
                % get rid of events which have NaN values between p_start
                % and p_end, then expand the selected data to the 
                % lowest level where it has no NaN values
                
                % add shrink_box, if range as large as box_max has NaN value, then shrink until box_min
                                
                shrink_box = true;
                NaN_tag = false;
                while shrink_box
                    x_ind = (box_max_x - computation_x) / 2 + 1 : (box_max_x + computation_x) / 2;
                    y_ind = center_y - (computation_y - 1) / 2 : center_y + (computation_y - 1) / 2;
                    temp = omega(y_ind, x_ind, p_start : p_end, :);
                    if any(isnan(temp(:)))
                        if computation_x > box_min && computation_y > box_min
                            computation_x = computation_x - 2;
                            computation_y = computation_y - 2;
                        else
                            shrink_box = false;
                            NaN_tag = true;
                        end
                    else 
                        while ~any(isnan(temp(:))) && p_start >= 1
                            p_start = p_start - 1;
                            if p_start == 0
                                break;
                            end
                            temp = omega(:, :, p_start : p_end, :);
                        end
                        p_start = p_start + 1;
                        shrink_box = false;
                    end
                end
                if NaN_tag
                    disp('NaN found in omega, continue');
                    continue;
                end
                level_indices = p_start : p_end;
                level = plevels(level_indices);
                T = T(y_ind, x_ind, level_indices, :);
                ug = ug(y_ind, x_ind, level_indices, :);
                vg = vg(y_ind, x_ind, level_indices, :);
                omega = omega(y_ind, x_ind, level_indices, :);
                omega_b = omega_b(y_ind, x_ind, level_indices, :);
                event_latspan = event_latspan(y_ind);
                event_lonspan = event_lonspan(x_ind);
                lat = lat_series_original(event_latspan);
                lon = lon_series_original(event_lonspan);
                phi = lat / 180.0 * 3.1415926;
                lambda = lon / 180.0 * 3.1415926;
                %event_timespan = 43800 + event_day + [0.25, 0.50, 0.75, 1.00];
                event_timespan = event_day + [- 0.25, 0.00, 0.25]; % convert to 'days since 1990'
                [~, Phi] = meshgrid(lambda, phi);

                % smooth T field
                T_smoothed = T;
                for k = 1 : length(level)
                    for t = 1 : length(event_timespan)
                        T_smoothed(:, :, k, t) = smooth2a(T(:, :, k, t), window_x, window_y);
                        T_smoothed(:, :, k, t) = smooth_normalize(T(:, :, k, t), T_smoothed(:, :, k, t));
                    end
                end
                T = T_smoothed;
                clear('T_smoothed');

                f0 = 2 * Omega * sin(lat_series_original(lat_indices(jj)) / 180 * 3.1415926);
                f = 2 * Omega * sin(lat_series_original(lat_indices) / 180 * 3.1415926);
                beta = 2 * Omega / R * cos(lat_series_original(lat_indices(jj)) / 180 * 3.1415926);
                [A, B, C, J, J1, J2, Qx, Qy, ...
                 dug_dlambda, dvg_dlambda, dug_dphi, dvg_dphi, ...
                 dT_dlambda, dT_dphi, sigma_accu] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
                [sigma, sigma_eff, dtheta_dp_ma, T_avg] = deal(zeros(length(level), length(event_timespan)));

                omega_QG = zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan));
                
                if sigma_tag
                    sigma_accu_smoothed = zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan));
                end
                
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
                    [dtheta_dp_ma(:, t), dtheta_dp_eff, theta_avg] = eff_stat_stab(level, T_avg(:, t), lambda_eff);
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
                    
                    % deal with points where sigma_accu is not smooth
                    if sigma_tag
                        if sigma_smooth_window ~= 0
                            for k = 1 : length(level)
                                sigma_accu_smoothed(:, :, k, t) = smooth2a(sigma_accu(:, :, k, t), sigma_smooth_window, sigma_smooth_window);
                            end
                        else
                            sigma_accu_smoothed(:, :, :, t) = sigma_accu(:, :, :, t);
                        end
                    end

                    [ug_smoothed, vg_smoothed] = deal(zeros(size(A(:, :, :, t))));

                    % term A
                    
                    for k = 1 : length(level)
                        
                        ug_smoothed(:, :, k) = smooth_normalize(ug(:, :, k, t), smooth2a(ug(:, :, k, t), 1, 1));
                        vg_smoothed(:, :, k) = smooth_normalize(vg(:, :, k, t), smooth2a(vg(:, :, k, t), 1, 1));

                        % probably the ug and vg field need smoothing as well...
                        dug_dlambda(:, :, k, t) = d_dlambda(ug(:, :, k, t), dlambda);
                        dvg_dlambda(:, :, k, t) = d_dlambda(vg(:, :, k, t), dlambda);
                        dug_dphi(:, :, k, t) = d_dphi(ug(:, :, k, t), dphi);
                        dvg_dphi(:, :, k, t) = d_dphi(vg(:, :, k, t), dphi);

                        dT_dlambda(:, :, k, t) = d_dlambda(T(:, :, k, t), dlambda);
                        dT_dphi(:, :, k, t) = d_dphi(T(:, :, k, t), dphi);
                        Qx(:, :, k, t) = 1 / R^2 * (...
                            (1 ./ cos(Phi).^2) .* dug_dlambda(:, :, k, t) .* dT_dlambda(:, :, k, t) + ...
                            (1 ./ cos(Phi))    .* dvg_dlambda(:, :, k, t) .* dT_dphi(:, :, k, t));
                        Qy(:, :, k, t) = 1 / R^2 * (...
                            (1 ./ cos(Phi))    .* dug_dphi(:, :, k, t)    .* dT_dlambda(:, :, k, t) + ...
                                                  dvg_dphi(:, :, k, t)    .* dT_dphi(:, :, k, t));
                        div_temp = div(Qx(:, :, k, t), Qy(:, :, k, t), Phi, dphi, dlambda);
                        A(:, :, k, t) = 2 * Ra / level(k) * div_temp;
                    end
                    clear('temp_x', 'temp_y', 'div_temp');
                                        
                    % term B
                    
                    %temp = d_dp(vg(:, :, :, t), level);
                    temp = d_dp(vg_smoothed, level);
                    for j = 1 : length(lat)
                        B(j, :, :, t) = f0 * beta * temp(j, :, :);
                    end
                    clear('temp');
                    
                    % term C
        
                    if t == 1
                        T1 = T(:, :, :, t);
                    end
                    if t < length(event_timespan)
                        T2 = T(:, :, :, t + 1);
                    end
                    for k = 1 : length(level)
                        
                        % calculate diabatic heating using temperature field and geostrophic wind
                        if t == length(event_timespan)
                            temp = (T1(:, :, k) - T0(:, :, k)) / dt;
                        elseif t == 1
                            temp = (T2(:, :, k) - T1(:, :, k)) / dt;
                        else
                            temp = (T2(:, :, k) - T0(:, :, k)) / (2 * dt);
                        end
                        
                        J1(:, :, k, t) = (temp + v_del(phi, ug(:, :, k, t), vg(:, :, k, t), T(:, :, k, t), dphi, dlambda)) * cp;
                        J2(:, :, k, t) = - sigma_accu(:, :, k, t) .* omega(:, :, k, t) / Ra * level(k) * cp;
                        J(:, :, k, t) = J1(:, :, k, t) + J2(:, :, k, t);
                        
                        % calculate the radiative and latent heating portion of the diabatic heating
                        % note that this calculation is only valid if the field is close to saturation, as discussed
                        % in O'Gorman, 2011
                        %J_latent(:, :, k, t) = omega(:, :, k, t) .* heaviside( - omega(:, :, k, t)) .* ...
                        %        repmat(dtheta_dp_ma(k, t) * T_avg(k, t) / theta_avg(k), size(omega(:, :, k, t))) .* cp;
                        %J_rad(:, :, k, t) = J(:, :, k, t) - J_latent(:, :, k, t);

                        C(:, :, k, t) = - kappa / level(k) * spherical_laplacian(J(:, :, k, t), phi, dphi, dlambda);

                    end

                    T0 = T1;
                    T1 = T2;

                    % right hand side
                    
                    rhs = A(:, :, :, t) + B(:, :, :, t) + C(:, :, :, t);
                    
                    % smoothing the rhs (this turns out to be redundant and it wipes out useful information)
                    %for k = 1 : length(level);
                    %    rhs(:, :, k) = smooth2a(rhs(:, :, k), window_x, window_y);
                    %end

        
                                         % SIP_diagnostic(R, f0, Ni,                   Nj,            Nk, phi, lambda, level, ...
                                                           %sigma,       dphi, dlambda, rhs,  omega)
                    if sigma_tag
                        omega_QG(:, :, :, t) = SIP_sigma(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
                                                            sigma_accu_smoothed(:, :, :, t), dphi, dlambda, rhs,  omega_b);
                    else
                        omega_QG(:, :, :, t) = SIP_diagnostic(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
                                                            sigma(:, t), dphi, dlambda, rhs,  omega_b);
                    end
                end
                
                

                if ~any(reshape(isnan(omega_QG(:, :, :, 2)), 1, []))
                    
                    % select the center of omega_QG
                    [event_lon, event_lat, tag] = find_center_omega(event_lon, event_lat, ...
                                                    omega_QG(:, :, level == temp_level, 2), length(lon_series_original));
                    if ~tag
                        disp('local omega_QG maxima too far from omega maxima, continue');
                        continue;
                    end
                    indx = find(lon == lon_series_original(event_lon));
                    indy = find(lat == lat_series_original(event_lat));
                    [Ny, Nx] = size(omega_QG(:, :, 1, 1));
                    y_indices = max(indy - (Ny - indy), 1) : min(2 * indy - 1, Ny);
                    x_indices = max(indx - (Nx - indx), 1) : min(2 * indx - 1, Nx);
                    %indy = find(event_latspan == lat_indices(jj));
                    temp_1 = omega_QG(indy, indx, :, 2);
                
                    events_sta{jj}(end + 1) = event;
                    
                    % information about the event

                    events_sta{jj}(end).event_timespan = event_timespan;
                    events_sta{jj}(end).event_year = event_year;
                    events_sta{jj}(end).event_latspan = event_latspan;
                    events_sta{jj}(end).event_lonspan = event_lonspan;
                    events_sta{jj}(end).computation_x = computation_x;
                    events_sta{jj}(end).computation_y = computation_y;
                    events_sta{jj}(end).event_lat = lat_series_original(lat_indices(jj));
                    events_sta{jj}(end).event_lon = lon_series(lons_{jj}(m));
                    events_sta{jj}(end).event_level = level;
                    events_sta{jj}(end).precip = precip_{jj}(m);
                    
                    % write in omega data

                    events_sta{jj}(end).omega = squeeze(omega(indy, indx, :, 2));
                    events_sta{jj}(end).omega_QG = squeeze(omega_QG(indy, indx, :, 2));
                    %events_sta{jj}(end).omega_QG_400hPa_full = omega_QG(y_indices, x_indices, level == 40000, 2);
                    events_sta{jj}(end).omega_QG_500hPa_full = omega_QG(y_indices, x_indices, level == 50000, 2);
                    %events_sta{jj}(end).omega_QG_600hPa_full = omega_QG(y_indices, x_indices, level == 60000, 2);
                    %events_sta{jj}(end).omega_QG_max = min(temp_1);
                    %events_sta{jj}(end).omega_QG_max_p = level(temp_1 == min(temp_1));
                    
                    % write in diagnostics
                    del_J = zeros(length(level), 1);
                    for k = 1 : length(level)
                        temp_J = spherical_laplacian(J(:, :, k, 2), phi, dphi, dlambda);
                        del_J(k) = temp_J(indy, indx);
                    end
                    clear('temp_J');
                    temp_J = reshape(J(indy, indx, :, 2), size(del_J));
                    l2 = - del_J ./ temp_J;
                    if sigma_tag
                        [del_omega_sigma_QG] = deal(zeros(length(level), 1));
                        for k = 1 : length(level)
                            temp_del_omega_sigma_QG = spherical_laplacian(omega_QG(:, :, k, 2) .* ...
                                sigma_accu(:, :, k, 2), phi, dphi, dlambda);
                            del_omega_sigma_QG(k) = temp_del_omega_sigma_QG(indy, indx);
                        end
                        temp_omega_QG = reshape(omega_QG(indy, indx, :, 2), size(del_omega_sigma_QG));
                        temp_sigma_accu = reshape(sigma_accu(indy, indx, :, 2), size(del_omega_sigma_QG));
                        k2 = - del_omega_sigma_QG ./ temp_omega_QG ./ temp_sigma_accu;
                    else
                        [del_omega_QG] = deal(zeros(length(level), 1));
                        for k = 1 : length(level)
                            temp_del_omega_QG = spherical_laplacian(omega_QG(:, :, k, 2), phi, dphi, dlambda);
                            del_omega_QG(k) = temp_del_omega_QG(indy, indx);
                        end
                        temp_omega_QG = reshape(omega_QG(indy, indx, :, 2), size(del_omega_QG));
                        k2 = - del_omega_QG ./ temp_omega_QG;
                    end

                    d2_dp2_omega_QG = d_dp(d_dp(omega_QG(indy, indx, :, 2), level), level);
                    m2 = - squeeze(d2_dp2_omega_QG) ./ squeeze(omega_QG(indy, indx, :, 2));

                    events_sta{jj}(end).k2 = k2;
                    events_sta{jj}(end).m2 = m2;
                    events_sta{jj}(end).l2 = l2;
                    events_sta{jj}(end).A = squeeze(A(indy, indx, :, 2));
                    events_sta{jj}(end).B = squeeze(B(indy, indx, :, 2));
                    events_sta{jj}(end).C = squeeze(C(indy, indx, :, 2));
                    %events_sta{jj}(end).C1 = C1(indy, indx, :, 2);
                    %events_sta{jj}(end).C2 = C2(indy, indx, :, 2);
                    %events_sta{jj}(end).Qx = Qx(indy-1:indy+1, indx-1:indx+1, :, 2);
                    %events_sta{jj}(end).Qy = Qy(indy-1:indy+1, indx-1:indx+1, :, 2);
                    %events_sta{jj}(end).dug_dlambda = dug_dlambda(indy, indx, :, 2);
                    %events_sta{jj}(end).dvg_dlambda = dvg_dlambda(indy, indx, :, 2);
                    %events_sta{jj}(end).dug_dphi = dug_dphi(indy, indx, :, 2);
                    %events_sta{jj}(end).dvg_dphi = dvg_dphi(indy, indx, :, 2);
                    %events_sta{jj}(end).dT_dlambda = dT_dlambda(indy, indx, :, 2);
                    %events_sta{jj}(end).dT_dphi = dT_dphi(indy, indx, :, 2);
                    events_sta{jj}(end).J_center = squeeze(J(indy, indx, :, 2));
                    %events_sta{jj}(end).J1 = squeeze(J1(indy, indx, :, 2));
                    %events_sta{jj}(end).J2 = squeeze(J2(indy, indx, :, 2));
                    %events_sta{jj}(end).J_400hPa = squeeze(J(y_indices, x_indices, level == 40000, 2));
                    events_sta{jj}(end).J_500hPa = squeeze(J(y_indices, x_indices, level == 50000, 2));
                    %events_sta{jj}(end).J_600hPa = squeeze(J(y_indices, x_indices, level == 60000, 2));
                    if sigma_tag
                        events_sta{jj}(end).sigma = squeeze(sigma_accu_smoothed(indy, indx, :, 2));
                    else
                        events_sta{jj}(end).sigma = sigma(:, 2);
                    end
                    events_sta{jj}(end).sigma_eff = sigma_eff(:, 2);
                    events_sta{jj}(end).sigma_accu = squeeze(sigma_accu(indy, indx, :, 2));
                    events_sta{jj}(end).T = squeeze(T(indy, indx, :, 2));
                    events_sta{jj}(end).T_avg = squeeze(T_avg(:, 2));

                    % update pointer

                    num_event(jj) = num_event(jj) + 1;
                    
                    %disp(num2str(length(events_sta{jj})));
                    
                    if NETCDF == true && mod(length(events_sta{jj}), NETCDF_interval) == 0
                        
                        disp('writing NetCDF file');
                        label = length(events_sta{jj});
                        
                        % write in QG omega field
                        [omega_QG_in, omega_in] = deal(zeros(length(lon(x_indices)), length(lat(y_indices)), ...
                                    length(level), length(event_timespan)));
                                                
                        for k = 1 : length(level)
                            for t = 1 : length(event_timespan)
                                omega_QG_in(:, :, k, t) = omega_QG(y_indices, x_indices, k, t)';
                            end
                        end

                        nc_filename = [output_path, string_2, '_', string_1, '_lat=', num2str(lat_series_original(lat_indices(jj))), ...
                                       '_', num2str(label, '%10.4d'), '_omega_QG1', '.nc'];
                        writeNetCDF(nc_filename, 'omega_QG', omega_QG_in, lat(y_indices), lon(x_indices), event_timespan, level);
					    
                        % write in event omega field
                        for k = 1 : length(level)
                            for t = 1 : length(event_timespan)
                                omega_in(:, :, k, t) = omega(y_indices, x_indices, k, t)';
                            end
                        end
                        nc_filename = [output_path, string_2, '_', string_1, '_lat=', num2str(lat_series_original(lat_indices(jj))), ...
                                       '_', num2str(label, '%10.4d'), '_omega', '.nc'];
                        writeNetCDF(nc_filename, 'omega', omega_in, lat(y_indices), lon(x_indices), event_timespan, level);
                        
                    end
                else
                    disp(['NaN found in omega_QG, not saved, m = ', num2str(m)]);
                    continue;
                end
            end
        else
            disp([num2str(lat_indices(jj)), ' statistics does not exist']);
        end

    end
    

    disp(num_event);
    
    save([output_path, matfilename], 'num_event', 'plevels', 'lon_series', 'lat_series', 'events_sta', '-v7.3');




