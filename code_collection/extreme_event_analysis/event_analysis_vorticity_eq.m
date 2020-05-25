
% compute event-wise QG_omega for precipitation extremes

function event_analysis_vorticity_eq(box_min, box_max, p_min, p_max, lat_series_original, lon_series_original, lat_indices, lon_indices, ...
                                    days_, precip_, matfilename, input_path, output_path, NETCDF, NETCDF_interval, n_ensemble, ...
                                    FIND_CENTER, JJA_DJF, ACCU_SIGMA, SIGMA_LOCAL, SIGMA_2, TRENBERTH)
    
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
    
    num_event = deal(zeros(length(lat_indices), length(lon_indices)));
    events_sta = deal(cell(length(lat_indices), length(lon_indices)));
    pwd_str = pwd;

    
    for jj = 1 : length(lat_indices)
        NetCDF_label = 0;
        for ii = 1 : length(lon_indices)
            disp(['lat: ', num2str(lat_series(jj)), 'lon: ', num2str(lon_series(ii))]);
            if ~isempty(days_{jj, ii})

                for m = 1 : length(days_{jj, ii})
                    disp(['m = ', num2str(m)])
                    event_day = days_{jj, ii}(m); 
                    event_lon = lon_indices(ii);
                    event_lat = lat_indices(jj);

                    if strfind(pwd_str, 'CESM')

                        box_max_x = 29;
                        box_max_y = 29;
                        computation_x = box_max_x;
                        computation_y = box_max_y;

                        event_year = floor((event_day - 0.25) / 365) + 1990;
                        string_1 = num2str(n_ensemble, '%.3d');
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
                        
                        box_max_x = 19;
                        box_max_y = 19;
                        computation_x = box_max_x;
                        computation_y = box_max_y;

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

                    [event_latspan, event_lonspan, center_y, computation_y] = ...
                            event_span_Wrapper(lat_series_original, lon_series_original, event_lat, event_lon, ...
                            box_max_y, box_max_x, computation_y);
                    
                    % compute indices to invert omega equation in
                    % lon_series_original and lat_series_original

                    plevels = double(ncread([input_path, 'ta_', string_1, '_', num2str(event_year), '.nc'], 'level'));
                    p_start_0 = 1; % lowest level
                    p_start = find(plevels == p_max);
                    p_end = find(plevels == p_min); % highest level
                    level_indices = p_start_0 : p_end;

                    % find center of omega at 500hPa
                    if strfind(pwd_str, 'CESM')
                        [omega, tag] = omega_read_Wrapper(computation_y, computation_x, length(level_indices), ...
                                event_year, input_path, event_latspan, event_lonspan, event_day, level_indices, string_1);
                    elseif strfind(pwd_str, 'GFDL')
                        [omega, tag] = omega_read_Wrapper_GFDL(computation_y, computation_x, length(level_indices), ...
                                event_year, input_path, event_latspan, event_lonspan, event_day, level_indices, string_1);
                    end

                    if ~tag; continue; end;
                    temp_level = 50000;
                    if FIND_CENTER
                        [event_lon, event_lat, tag] = find_center_omega(event_lon, event_lat, ...
                                                        omega(:, :, level_indices(plevels == temp_level), 2), length(lon_series_original));
                    end
                    if ~tag
                        disp('local omega maxima too far from precip maxima, continue');
                        continue;
                    end
                    

                    [event_latspan, event_lonspan, center_y, computation_y] = ...
                            event_span_Wrapper(lat_series_original, lon_series_original, event_lat, event_lon, ...
                            box_max_y, computation_x, computation_y);

                    % read in data
                    if strfind(pwd_str, 'CESM')
                        %[T, ug, vg, omega, omega_b, tag] = event_read_Wrapper(computation_y, computation_x, length(level_indices), ...
                        %        event_year, input_path, event_latspan, event_lonspan, event_day, level_indices, string_1, JJA_DJF);
                        [T, ua, va, ug, vg, omega, omega_b, tag] = event_read_Wrapper_uava(computation_y, computation_x, length(level_indices), ...
                                event_year, input_path, event_latspan, event_lonspan, event_day, level_indices, string_1, JJA_DJF);
                    elseif strfind(pwd_str, 'GFDL')
                        %[T, ug, vg, omega, omega_b, tag] = event_read_Wrapper_GFDL(computation_y, computation_x, length(level_indices), ...
                        %        event_year, input_path, event_latspan, event_lonspan, event_day, level_indices, string_1);
                        [T, ua, va, ug, vg, omega, omega_b, tag] = event_read_Wrapper_uava_GFDL(computation_y, computation_x, length(level_indices), ...
                                event_year, input_path, event_latspan, event_lonspan, event_day, level_indices, string_1);
                    end
                    if ~tag; continue; end;

                    % get rid of events which have NaN values between p_start
                    % and p_end, then expand the selected data to the 
                    % lowest level where it has no NaN values

                    % add shrink_box, if range as large as box_max has NaN value, then shrink until box_min

                    [x_ind, y_ind, computation_x, computation_y, p_start, NaN_tag] = ...
                            shrink_box_v1(omega, box_max_x, box_max_y, box_min, p_max, computation_x, computation_y, p_start, p_end, plevels);
                    %{  
                    shrink_box = true;
                    NaN_tag = false;
                    p_start = 1;
                    while shrink_box
                        x_ind = (box_max_x - computation_x) / 2 + 1 : (box_max_x + computation_x) / 2;
                        y_ind = (box_max_y - computation_y) / 2 + 1 : (box_max_y + computation_y) / 2;
                        temp = omega(y_ind, x_ind, p_start : p_end, :);
                        if any(isnan(temp(:))) && (computation_x > box_min && computation_y > box_min)
                            computation_x = computation_x - 2;
                            computation_y = computation_y - 2;
                        elseif any(isnan(temp(:)))
                            while any(isnan(temp(:))) && p_start <= find(plevels == p_max)
                                p_start = p_start + 1;
                                if p_start == (find(plevels == p_max) + 1)
                                    break;
                                    NaN_tag = true;
                                    shrink_box = false;
                                end
                                temp = omega(y_ind, x_ind, p_start : p_end, :);
                            end
                            shrink_box = false;
                        else
                            shrink_box = false;
                        end
                    end
                    %}
                    if NaN_tag
                        disp('NaN found in omega, continue');
                        continue;
                    end

                    level_indices = p_start : p_end;
                    level = plevels(level_indices);
                    T = T(y_ind, x_ind, level_indices, :);
                    ug = ug(y_ind, x_ind, level_indices, :);
                    vg = vg(y_ind, x_ind, level_indices, :);
                    ua = ua(y_ind, x_ind, level_indices, :);
                    va = va(y_ind, x_ind, level_indices, :);
                    omega = omega(y_ind, x_ind, level_indices, :);
                    omega_b = omega_b(y_ind, x_ind, level_indices, :);
                    event_latspan = event_latspan(y_ind);
                    event_lonspan = event_lonspan(x_ind);
                    lat = lat_series_original(event_latspan);
                    lon = lon_series_original(event_lonspan);
                    phi = lat / 180.0 * 3.1415926;
                    lambda = lon / 180.0 * 3.1415926;
                    event_timespan = event_day + [- 0.25, 0.00, 0.25]; % convert to 'days since 1990'
                    %event_timespan = 43800 + event_day + [0.25, 0.50, 0.75, 1.00];
                    [~, Phi] = meshgrid(lambda, phi);

                    % smooth T field
                    T_smoothed = T;
                    for k = 1 : length(level)
                        for t = 1 : length(event_timespan)
                            T_smoothed(:, :, k, t) = smooth2a(T(:, :, k, t), window_x, window_y);
                            % the smooth_normalize here could be problematic when box_min is very small
                            T_smoothed(:, :, k, t) = smooth_normalize(T(:, :, k, t), T_smoothed(:, :, k, t));
                        end
                    end
                    T = T_smoothed;
                    clear('T_smoothed');

                    f0 = 2 * Omega * sin(lat_series_original(lat_indices(jj)) / 180 * 3.1415926);
                    beta = 2 * Omega / R * cos(lat_series_original(lat_indices(jj)) / 180 * 3.1415926);
                    [A, B, C, D, rhs, J, Qx, Qy, ...
                     dug_dlambda, dvg_dlambda, dug_dphi, dvg_dphi, ...
                     dT_dlambda, dT_dphi, sigma_accu, SSigma, ...
                     zeta, ug_del_zeta, ua_del_zeta, theta] = ...
                                deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
                    [sigma, sigma_eff, dtheta_dp_ma, T_avg, dtheta_dp_ma_avg] = deal(zeros(length(level), length(event_timespan)));

                    omega_QG = zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan));

                    %% calculate stabilities

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
                        Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
                        theta(:, :, :, t) = T(:, :, :, t) .* (1e5 ./ Level) .^ kappa;
                        sigma_accu(:, :, :, t) = - Ra * T(:, :, :, t) ./ (Level .* theta(:, :, :, t)) .* d_dp(theta(:, :, :, t), level);
                    end

                    % calculate vorticity

                    for t = 1 : length(event_timespan)
                        for k = 1 : length(level)
                            %zeta(:, :, k, t) = curl(phi, ug(:, :, k, t), vg(:, :, k, t), dphi, dlambda);
                            zeta(:, :, k, t) = curl(phi, ua(:, :, k, t), va(:, :, k, t), dphi, dlambda);
                        end
                    end

                    % term A
                    %[A1, A2] = traditional_A(phi, level, ug, vg, zeta, event_timespan, dphi, dlambda, f0, beta);
                    [A, A1, A2, A3] = traditional_A(phi, level, ua, va, zeta, event_timespan, dphi, dlambda, f0, beta);

                    % term B
                    % geostrophic version:
                    %B = traditional_B(phi, level, ug, vg, T, event_timespan, dphi, dlambda);
                    % non-approximated version:
                    %B = traditional_B(phi, level, ua, va, T, event_timespan, dphi, dlambda);
                    %B_backup = B;

                    % Q vector
                    %[A, B] = Q_vector(level, ug, vg, T, Phi, event_timespan, dphi, dlambda, f0, beta);

                    % term C (latent heating)
                    [C1, C2, J1, J2] = diabatic_heating(phi, T, level, ua, va, sigma_accu, omega, ...
                                              event_timespan, dphi, dlambda, dt);
                    J = J1 + J2;
                    C = C1 + C2;
                    
                    % term B (change of vorticity)
                    B = vorticity_tendency(zeta, level, event_timespan, f0, dt);

                    %{
                    % compare Q-vector and traditional form
                    % two terms that cancel each other:
                    k = find(level == 50000); t = 2;
                    dug_dp = d_dp(ug(:, :, :, t), level);
                    dvg_dp = d_dp(vg(:, :, :, t), level);
                    dzeta_dp = d_dp(zeta(:, :, :, t), level);
                    A1 = f0 * v_del(phi, dug_dp(:, :, k), dvg_dp(:, :, k), zeta(:, :, k, t), dphi, dlambda);
                    B1 = Ra .* T(:, :, k, t) ./ (theta(:, :, k, t) .* level(k)) .* ...
                            v_del(phi, spherical_laplacian(ug(:, :, k, t), phi, dphi, dlambda), ...
                                       spherical_laplacian(vg(:, :, k, t), phi, dphi, dlambda), theta(:, :, k, t), dphi, dlambda);
                    A2 = f0 * v_del(phi, ug(:, :, k, t), vg(:, :, k, t), dzeta_dp(:, :, k), dphi, dlambda);
                    B2 = Ra .* T(:, :, k, t) ./ (theta(:, :, k, t) .* level(k)) .* ...
                            v_del(phi, ug(:, :, k, t), vg(:, :, k, t), ...
                            spherical_laplacian(theta(:, :, k, t), phi, dphi, dlambda), dphi, dlambda);
                    %A3 = f0 * beta * dvg_dp(:, :, k, t);
                    %B3 = 2 * Ra .* T(:, :, k, t) ./ (theta(:, :, k, t) .* level(k)) .* ...
                    %}      
                    
                    % smooth rhs
                    A1_smoothed = A1;
                    A2_smoothed = A2;
                    A3_smoothed = A3;
                    B_smoothed = B;
                    for k = 1 : length(level)
                        for t = 1 : length(event_timespan)
                            A1_smoothed(:, :, k, t) = smooth2a(A1(:, :, k, t), window_x, window_y);
                            A2_smoothed(:, :, k, t) = smooth2a(A2(:, :, k, t), window_x, window_y);
                            A3_smoothed(:, :, k, t) = smooth2a(A3(:, :, k, t), window_x, window_y);
                            B_smoothed(:, :, k, t) = smooth2a(B(:, :, k, t), window_x, window_y);
                            %rhs_smoothed(:, :, k, t) = smooth2a(rhs(:, :, k, t), window_x, window_y);
                            %rhs_smoothed(:, :, k, t) = smooth_normalize(rhs(:, :, k, t), rhs_smoothed(:, :, k, t));
                        end
                    end
                    A1 = A1_smoothed;
                    A2 = A2_smoothed;
                    A3 = A3_smoothed;
                    B = B_smoothed;
                    
                    % Trenberth version
                    if TRENBERTH
                        A = A1 + A3;
                        rhs = A;
                    else
                        A = A1 + A2 + A3;
                        rhs = A + B;
                    end
                    
                    %rhs_vort = A1 + A2 + A3 + D;
                    %rhs_thermo = -C1 + C;
                    %rhs_Q_vector = A + B + C;

                    [omega_QG_thermo, omega_QG_Q_vector, omega_QG_vort] = deal(zeros(size(omega_QG)));
                    for t = 1 : length(event_timespan)
                        %======= thermodynamic equation with original boundary condition =======%
                        % this gives you exactly the original omega field, but may be unstable
                        %omega_QG_thermo(:, :, :, t) = SIP_diagnostic_v1_sigma(R, 0, length(lon), length(phi), length(level), phi, lambda, level, ...
                        %                                    sigma_accu(:, :, :, t), dphi, dlambda, rhs_thermo(:, :, :, t),  omega(:, :, :, t), false);
                        %======= thermodynamic equation with climatological boundary condition =======% 
                        %omega_QG_thermo(:, :, :, t) = SIP_diagnostic_v1_sigma(R, 0, length(lon), length(phi), length(level), phi, lambda, level, ...
                        %                                    sigma_accu(:, :, :, t), dphi, dlambda, rhs_thermo(:, :, :, t),  omega_b, false);
                        %======= vorticity equation =======%
                        omega_QG(:, :, :, t)   = SIP_diagnostic_v1(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
                                                            zeros(size(sigma, 1), 1), dphi, dlambda, rhs(:, :, :, t),  omega_b, false);
                        %======= full QG equation =======%
                        %omega_QG_Q_vector(:, :, :, t) = SIP_diagnostic_v1_sigma(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
                        %                                    sigma_accu(:, :, :, t), dphi, dlambda, rhs_Q_vector(:, :, :, t),  omega_b, false);

                    end
                    

                    
                    if ~any(isnan(omega_QG(:, :, level == 50000, 2)))

                        % select the center of omega_QG
                        if FIND_CENTER
                        	[event_lon, event_lat, tag] = find_center_omega(event_lon, event_lat, ...
                            	                            omega_QG(:, :, level == temp_level, 2), length(lon_series_original));
                        end
	                    
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
                        
                        events_sta{jj, ii}(end + 1) = event;

                        % information about the event

                        events_sta{jj, ii}(end).event_timespan = event_timespan;
                        events_sta{jj, ii}(end).event_year = event_year;
                        events_sta{jj, ii}(end).event_latspan = event_latspan;
                        events_sta{jj, ii}(end).event_lonspan = event_lonspan;
                        events_sta{jj, ii}(end).computation_x = computation_x;
                        events_sta{jj, ii}(end).computation_y = computation_y;
                        events_sta{jj, ii}(end).event_lat = lat_series_original(lat_indices(jj));
                        events_sta{jj, ii}(end).event_lon = lon_series(ii);
                        events_sta{jj, ii}(end).event_level = level;
                        events_sta{jj, ii}(end).precip = precip_{jj, ii}(m);
                        events_sta{jj, ii}(end).omega_500hPa_full = omega(y_indices, x_indices, level == 50000, 2);

                        % write in omega data

                        events_sta{jj, ii}(end).omega = squeeze(omega(indy, indx, :, 2));
                        events_sta{jj, ii}(end).omega_QG = squeeze(omega_QG(indy, indx, :, 2));
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
                        events_sta{jj, ii}(end).A1 = squeeze(A1(indy, indx, :, 2));
                        events_sta{jj, ii}(end).A2 = squeeze(A2(indy, indx, :, 2));
                        events_sta{jj, ii}(end).A3 = squeeze(A3(indy, indx, :, 2));
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

                        % update pointer

                        %%%%%%%%%%%%%% DON'T DELETE %%%%%%%%%%%%%%%%%%
                        % test if events are recorded accurately
                        if rand < 0.01
                            if TRENBERTH
                                temp = (events_sta{jj, ii}(end).m2 .* f0^2) .* events_sta{jj, ii}(end).omega_QG + ...
                                        events_sta{jj, ii}(end).A;
                            else
                                temp = (events_sta{jj, ii}(end).m2 .* f0^2) .* events_sta{jj, ii}(end).omega_QG + ...
                                        events_sta{jj, ii}(end).A + events_sta{jj, ii}(end).B;
                            end
                            disp(['If this number: ', num2str(temp(floor(end/2))), ' is larger than 1e-30, ', ...
                                  'something is wrong about event recording']);
                        end


                        num_event(jj, ii) = num_event(jj, ii) + 1;
                        
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
            else
                disp(['jj = ', num2str(lat_indices(jj)), ...
                      'ii = ', num2str(lon_indices(ii)), ' statistics does not exist']);
            end
        end
    end
    

    disp(num_event);
    
    save([output_path, matfilename], 'num_event', 'plevels', 'lon_series', 'lat_series', 'events_sta', '-v7.3');




