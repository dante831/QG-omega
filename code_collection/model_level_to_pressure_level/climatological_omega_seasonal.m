    
    % get the averaged surface pressure
    
    addpath('/net/aimsir/archive1/ziweili/CESM_LENS/source');

    mon_days = [31,  59,  90, 120, 151, 181, 212, 243, 273, 304, 334, 365];
    months = [12, 1, 2];
    latmax = 90; % degrees
    latmin = -90;
    lonmax = 360; % degrees
    lonmin = 0;
    n_ensembles = [1, 2, 3, 4, 5, 35];

    plevels = [100000, 97500, 95000, 92500, 90000, 87500, 85000, 82500, 80000, 77500, 75000, 70000, 65000, ...
                60000, 55000, 50000, 45000, 40000, 35000, 30000, 25000, 22500, 20000, 17500, 15000, 12500, ...
                10000,  7000,  5000,  3000,  2000,  1000,   700,   500,   300,   200,   100];

    data_path = '/net/aimsir/archive1/ziweili/CESM_LENS/data/';
    input_path = '/net/aimsir/archive1/ziweili/CESM_LENS/output/';
    output_path = input_path;

    day_ind = [];
    for m = 1 : length(months)
        if months(m) == 1
            day_ind = [day_ind, 1:(31*4)];
        else
            day_ind = [day_ind, mon_days(months(m) - 1)*4 + 1 : mon_days(months(m))*4];
        end
    end

    ps_filename_0 = [data_path, 'b.e11.B20TRC5CNBDRD.f09_g16.001.cam.h2.PS.1990010100Z-2005123118Z.nc'];
    lat_series = ncread(ps_filename_0, 'lat');
    lon_series = ncread(ps_filename_0, 'lon');
    [lat_indices, lon_indices] = latlonindices(lat_series, lon_series, latmin, latmax, lonmin, lonmax);
    lon = lon_series(lon_indices);
    lat = lat_series(lat_indices);

    for climate = 1 : 2
    
        if climate == 1
            string_2 = 'h';
            years = linspace(1991, 2000, 10);
        else
            string_2 = 'r';
            years = linspace(2071, 2080, 10);
        end   

        omega_avg_ensemble = {};
        for e = 1 : length(n_ensembles)

            string_1 = num2str(n_ensembles(e), '%.3d');
            disp(['ensemble ', string_1, ' for ', string_2]);

            omega_avg_ensemble{e} = zeros(length(lon), length(lat), length(plevels));
            for k = 1 : length(plevels)
                temp_omega = zeros(length(lon), length(lat), length(years));
                for f = 1 : length(years)
                    omega_filename = [input_path, 'omega_', string_1, '_', num2str(years(f)), '.nc'];
                    start = [lon_indices(1), lat_indices(1), k, 1];
                    count = [length(lon_indices), length(lat_indices), 1, Inf];
                    temp = squeeze(ncread(omega_filename, 'omega_plevel', start, count));
                    temp_omega(:, :, f) = mean(temp(:, :, day_ind), 3, 'omitnan');
                end
                omega_avg_ensemble{e}(:, :, k) = mean(temp_omega, 3, 'omitnan');
            end
        end
        omega_avg = (omega_avg_ensemble{1} + ...
                     omega_avg_ensemble{2} + ...
                     omega_avg_ensemble{3} + ...
                     omega_avg_ensemble{4} + ...
                     omega_avg_ensemble{5} + ...
                     omega_avg_ensemble{6}) / 6;
        nc_filename = [output_path, 'omega_avg_', string_2, '_DJF', '.nc'];
        varname = 'omega_avg';
        writeNetCDF_omega(nc_filename, varname, omega_avg, lat, lon, plevels);
    end


