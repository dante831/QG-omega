
    % convert small daily .nc files into large yearly .nc files

    output_path = '/archive1/ziweili/CESM_LENS/output/';
    %years = [1991 : 2000]; 
    years = [2071 : 2080];
    n_ensemble = 5;
    
    %varnames = {'omega'; 'ua'; 'va'; 'ta'; 'ug'; 'vg'};
    %varnames_file = {'omega_plevel'; 'ua_plevel'; 'va_plevel'; 'ta_plevel'; 'ug'; 'vg'};
    varnames = {'q'};
    varnames_file = {'q_plevel'};

    for y = 1 : length(years)
    
        year = years(y);
        disp(['year of ', num2str(year)]);
    
        for v = 1 : length(varnames)
        
            disp(['var: ', varnames{v}]);
            files = dir([output_path, varnames{v}, '_', num2str(n_ensemble, '%.3d'), '_', num2str(year), '_*.nc']);
            for f = 1 : length(files)
                files(f).name = [output_path, files(f).name];
            end
            lat = ncread(files(1).name, 'latitude');
            lon = ncread(files(1).name, 'longitude');
            plevels = ncread(files(1).name, 'level');
            time_start = double(ncread(files(1).name, 'time', 1, 1));
            time = linspace(time_start, time_start + length(files) - 0.25, length(files) * 4);
        
            var = zeros([length(lon), length(lat), length(plevels), length(time)]);
            for f = 1 : length(files)
                var(:, :, :, (f - 1) * 4 + 1 : (f - 1) * 4 + 4) = ncread(files(f).name, varnames_file{v});
            end
        
            nc_filename = [output_path, varnames{v}, '_', num2str(n_ensemble, '%.3d'), '_', num2str(years(y), '%10.4d'), '.nc'];
            writeNetCDF(nc_filename, varnames_file{v}, var, lat, lon, time, plevels);
    
        end
    
    end
