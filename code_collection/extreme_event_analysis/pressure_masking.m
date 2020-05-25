function precip = pressure_masking(precip, small_value, lon_indices, lat_indices, p_max, filename)

    pwd_str = pwd;
    if strfind(pwd_str, 'CESM')
        if strfind(filename, 'B20TRC5CNBDRD')
            string = 'historical';
        elseif strfind(filename, 'BRCP85C5CNBDRD')
            string = 'rcp85';
        end
        %ps_filename = ['/net/aimsir/archive1/ziweili/CESM_LENS/output/ps_avg_', string, '.nc'];
        ps_filename = ['/archive1/ziweili/CESM_LENS/output/ps_avg_', string, '.nc'];
    elseif strfind(pwd_str, 'GFDL')
        if strfind(filename, 'historical')
            string = 'historical';
        elseif strfind(filename, 'rcp85')
            string = 'rcp85';
        end
        ps_filename = ['/net/chickpea/volume1/scratch/ziweili/test1_GFDL/ps/ps_avg_', string, '.nc'];
    end
    
    lat = ncread(ps_filename, 'latitude');
    lon = ncread(ps_filename, 'longitude');
    ps_avg = ncread(ps_filename, 'ps_avg');
    ind = ps_avg(lon_indices, lat_indices) >= p_max;
    %{
    if ocean_only
        [Lat, Lon] = meshgrid(lat(lat_indices), lon(lon_indices));
        % mask out land areas where tag is equal to 0
        Tag = reshape(event_tag(Lon(:), Lat(:)), size(Lat)) ~= 0;
        ind = p_ind & Tag;
    else
        ind = p_ind;
    end
    %}

    for t = 1 : length(precip(1, 1, :))
        temp_precip = precip(:, :, t);
        temp_precip(~ind) = small_value;
        precip(:, :, t) = temp_precip;
    end

end
