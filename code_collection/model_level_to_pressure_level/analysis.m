function analysis(n_ensemble, year, output_path, input_path, latmax, latmin, lonmax, lonmin, plevels)    

    % GFDL data description:

    % the data is archived as: [day, level, lat, lon]
    % but matlab reads the data as: [lon, lat, level, day], [144, 90, 48, 1460]

    addpath('source/')
    % some constants:

    Ra = 287.04;
    g = 9.8;
    SMOOTH = true;
    days_normal = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    %days_leap = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

    % define several selection parameters:
    if year <= 2005
        filename_part_1 = ['b.e11.B20TRC5CNBDRD.f09_g16.', num2str(n_ensemble, '%.3d'), '.cam.h2.'];
        filename_part_2 = ['.1990010100Z-2005123118Z.nc'];
        time_start = (year - 1990) * 365 * 4;
    else
        filename_part_1 = ['b.e11.BRCP85C5CNBDRD.f09_g16.', num2str(n_ensemble, '%.3d'), '.cam.h2.'];
        filename_part_2 = ['.2071010100Z-2080123118Z.nc'];
        time_start = (year - 2071) * 365 * 4;
    end

    ta_filename = [input_path, filename_part_1, 'T', filename_part_2];
    lat_series = ncread(ta_filename, 'lat');
    lon_series = ncread(ta_filename, 'lon');
    [lat_indices, lon_indices] = latlonindices(lat_series, lon_series, latmin, latmax, lonmin, lonmax);
    lat = lat_series(lat_indices);
    lon = lon_series(lon_indices);
    phi = lat / 180.0 * 3.1415926;
    lambda = lon / 180.0 * 3.1415926;
    [Phi_3D, Lambda_3D, level_3D] = meshgrid(phi, lambda, plevels);

    days = days_normal;
    disp(['Year of ', num2str(year)]);
    ua_filename = [input_path, filename_part_1,  'U', filename_part_2];
    va_filename = [input_path, filename_part_1,  'V', filename_part_2];
    ps_filename = [input_path, filename_part_1, 'PS', filename_part_2];
    q_filename  = [input_path, filename_part_1,  'Q', filename_part_2];
    time_original = ncread(ta_filename, 'time');
    time_indices_original = linspace(1, length(time_original), length(time_original));
    %time_indices_original = linspace(1, 124, 124);

    ndays = floor(length(time_indices_original) / 4.0);
    nmonths = 12;

    for month = 1 : nmonths
        for day = 1 : min(days(month), ndays - sum(days(1 : month - 1)))
            disp(['month: ', num2str(month), ' day: ', num2str(day), ' ', datestr(now)]);
            time_indices = linspace(time_start + (sum(days(1 : month - 1)) + day - 1) * 4 + 1, ...
                                    time_start + (sum(days(1 : month - 1)) + day - 1) * 4 + 4, 4);
            time = time_original(time_indices);
            
            % calculate pressure
    
            start = [1, lat_indices(1), time_indices(1)];
            count = [Inf, length(lat_indices), length(time_indices)];
            a = ncread(ua_filename, 'hyam');
            b = ncread(ua_filename, 'hybm');
            p0 = ncread(ua_filename, 'P0');
            ps = ncread(ps_filename, 'PS', start, count);
            p = zeros([length(lon_indices), length(lat_indices), length(a), length(time_indices)]);
            for t = 1 : length(time_indices)
                for k = 1 : length(a)
                    for j = 1 : length(lat_indices)
                        for i = 1 : length(lon_indices)
                            p(i, j, k, t) = a(k) * p0 + b(k) * ps(i, j, t);
                        end
                    end
                end
            end

            start = [1, lat_indices(1), 1, time_indices(1)];
            count = [Inf, length(lat_indices), Inf, length(time_indices)];
            
            %%%%%%%%% DON'T DELETE! %%%%%%%%%%%
            %{
            % compute omega from continuity equation
  
            disp('compute omega...');
            omega_plevel = zeros([length(lon_indices), length(lat_indices), length(plevels), length(time_indices)]);
            omega_eta = omega_from_continuity_v1(ua_filename, va_filename, ps_filename, lon_indices, ...
                                lat_indices, time_indices);
            for t = 1 : length(time_indices)
                omega_plevel(:, :, :, t) = ND_interp_v1(omega_eta(:, :, :, t), p(:, :, :, t), plevels, ps(:, :, t), 'linear');
            end
            %clear('omega_eta');
            if SMOOTH == true
                window_x = 1;
                window_y = 1;
                for k = 1 : length(plevels)
                    for t = 1 : length(time)
                        omega_plevel(:, :, k, t) = smooth2a_periodic(omega_plevel(:, :, k, t), window_x, window_y, 1);
                    end
                end
            end

            % output to NetCDF file
    
            nc_filename = [output_path, 'omega_', num2str(n_ensemble, '%.3d'), '_', num2str(year, '%10.4d'), ...
                            '_', num2str(month, '%10.2d'), '_', num2str(day, '%10.2d'), '.nc'];
            varname = 'omega_plevel';
            writeNetCDF(nc_filename, varname, omega_plevel, lat, lon, time, plevels);       
            clear('omega_plevel');

            %% interpolate other variables 
            % temperature

            disp('interpolate other variables...');
                
            disp('interpolate ta...');
            ta = ncread(ta_filename, 'T', start, count);
            ta = ta(lon_indices, :, :, :); % this treatment is to deal with areas across lon = 0 line 
            ta_plevel = zeros([length(lon_indices), length(lat_indices), length(plevels), length(time_indices)]);
            for t = 1 : length(time_indices)
                ta_plevel(:, :, :, t) = ND_interp_v1(ta(:, :, :, t), p(:, :, :, t), plevels, ps(:, :, t), 'linear');
            end
            clear('ta');
            
            nc_filename = [output_path, 'ta_', num2str(n_ensemble, '%.3d'), '_', num2str(year, '%10.4d'), ...
                            '_', num2str(month, '%10.2d'), '_', num2str(day, '%10.2d'), '.nc'];
            varname = 'ta_plevel';
            writeNetCDF(nc_filename, varname, ta_plevel, lat, lon, time, plevels);
            clear('ta_plevel');
            
            % x-direction wind

            disp('interpolate ua...');
            ua = ncread(ua_filename, 'U', start, count);
            ua = ua(lon_indices, :, :, :);
            ua_plevel = zeros([length(lon_indices), length(lat_indices), length(plevels), length(time_indices)]);
            for t = 1 : length(time_indices)
                ua_plevel(:, :, :, t) = ND_interp_v1(ua(:, :, :, t), p(:, :, :, t), plevels, ps(:, :, t), 'linear');
            end
            clear('ua');
            nc_filename = [output_path, 'ua_', num2str(n_ensemble, '%.3d'), '_', num2str(year, '%10.4d'), ...
                            '_', num2str(month, '%10.2d'), '_', num2str(day, '%10.2d'), '.nc'];
            varname = 'ua_plevel';
            writeNetCDF(nc_filename, varname, ua_plevel, lat, lon, time, plevels);
            clear('ua_plevel');

            % y-direction wind

            disp('interpolate va...');
            va = ncread(va_filename, 'V', start, count);
            va = va(lon_indices, :, :, :);
            va_plevel = zeros([length(lon_indices), length(lat_indices), length(plevels), length(time_indices)]);
            for t = 1 : length(time_indices)
                va_plevel(:, :, :, t) = ND_interp_v1(va(:, :, :, t), p(:, :, :, t), plevels, ps(:, :, t), 'linear');
            end
            clear('va');
            nc_filename = [output_path, 'va_', num2str(n_ensemble, '%.3d'), '_', num2str(year, '%10.4d'), ... 
                            '_', num2str(month, '%10.2d'), '_', num2str(day, '%10.2d'), '.nc'];
            varname = 'va_plevel';
            writeNetCDF(nc_filename, varname, va_plevel, lat, lon, time, plevels);
            clear('va_plevel');
            

            % specific humidity

            disp('interpolate Q...');
            q = ncread(q_filename, 'Q', start, count);
            q = q(lon_indices, :, :, :);
            q_plevel = zeros([length(lon_indices), length(lat_indices), length(plevels), length(time_indices)]);
            for t = 1 : length(time_indices)
                q_plevel(:, :, :, t) = ND_interp_v1(q(:, :, :, t), p(:, :, :, t), plevels, ps(:, :, t), 'linear');
            end
            clear('q')
            nc_filename = [output_path, 'q_', num2str(n_ensemble, '%.3d'), '_', num2str(year, '%10.4d'), ...
                            '_', num2str(month, '%10.2d'), '_', num2str(day, '%10.2d'), '.nc'];
            varname = 'q_plevel';
            writeNetCDF(nc_filename, varname, q_plevel, lat, lon, time, plevels);
            clear('q_plevel');
            %}

            clear('p');
        end
    end
    
end
    

