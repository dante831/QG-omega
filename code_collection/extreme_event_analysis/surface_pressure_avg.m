    
    % get the averaged surface pressure
    addpath('/net/aimsir/archive1/ziweili/CESM_LENS/source/');

    latmax = 90; % degrees
    latmin = -90;
    lonmax = 360; % degrees
    lonmin = 0;
    n_ensemble = 1;

    %rcp85path = strcat('/net/chickpea/volume1/scratch/jgdwyer/cmip5/nomads.gfdl.noaa.gov/dods-data/CMIP5/output1/', ...
    %            'NOAA-GFDL/GFDL-CM3/rcp85/6hr/atmos/6hrLev/r1i1p1/v20110601/ps/');
    %historicalpath = strcat('/net/chickpea/volume1/scratch/jgdwyer/cmip5/nomads.gfdl.noaa.gov/dods-data/CMIP5/output1/', ...
    %                'NOAA-GFDL/GFDL-CM3/historical/6hr/atmos/6hrLev/r1i1p1/v20110601/ps/');
    path = '/net/aimsir/archive1/ziweili/CESM_LENS/data/';
    output_path = '/net/aimsir/archive1/ziweili/CESM_LENS/output/';
    %output_path = '/net/chickpea/volume1/scratch/ziweili/test1_GFDL/ps/';

    ps_filename_0 = [path, 'b.e11.B20TRC5CNBDRD.f09_g16.001.cam.h2.PS.1990010100Z-2005123118Z.nc'];
    lat_series = ncread(ps_filename_0, 'lat');
    lon_series = ncread(ps_filename_0, 'lon');
    [lat_indices, lon_indices] = latlonindices(lat_series, lon_series, latmin, latmax, lonmin, lonmax);
    lon = lon_series(lon_indices);
    lat = lat_series(lat_indices);

    for climate = 1 : 2
        if climate == 1
            years = linspace(1990, 2005, 16);
            string = 'historical';
            ps_filename = [path, 'b.e11.B20TRC5CNBDRD.f09_g16.', num2str(n_ensemble, '%.3d'), '.cam.h2.PS.1990010100Z-2005123118Z.nc'];
        else
            years = linspace(2071, 2080, 10);
            string = 'rcp85';
            ps_filename = [path, 'b.e11.BRCP85C5CNBDRD.f09_g16.', num2str(n_ensemble, '%.3d'), '.cam.h2.PS.2071010100Z-2080123118Z.nc'];
        end
        ps = zeros(length(lon), length(lat), length(years));
        start = [lon_indices(1), lat_indices(1), 1];
        count = [length(lon_indices), length(lat_indices), Inf];
        ps = ncread(ps_filename, 'PS', start, count);
        ps_avg = mean(ps, 3);
        %{
        for f = 1 : length(years)
            ps_filename = [path, 'ps_6hrLev_GFDL-CM3_', string, '_r1i1p1_', num2str(years(f)), '010100-', num2str(years(f)), '123123.nc'];
            start = [lon_indices(1), lat_indices(1), 1];
            count = [length(lon_indices), length(lat_indices), Inf];
            temp_ps = ncread(ps_filename, 'ps', start, count);
            ps(:, :, f) = mean(temp_ps, 3);
        
        end
        %}
        nc_filename = [output_path, 'ps_avg_', string, '.nc'];
        varname = 'ps_avg';
        writeNetCDF_ps(nc_filename, varname, ps_avg, lat, lon);
    end


