function precip = geology_masking(precip, small_value, lon_indices, lat_indices, ocean_only, land_only)

    ps_filename = ['/net/chickpea/volume1/scratch/ziweili/test1_GFDL/ps/ps_avg_rcp85.nc'];

    lat = ncread(ps_filename, 'latitude');
    lon = ncread(ps_filename, 'longitude');
    [Lat, Lon] = meshgrid(lat(lat_indices), lon(lon_indices));
    ind = true(size(Lat));
    if ocean_only
        % select ocean areas where tag is not equal to 0
        Tag1 = reshape(event_tag(Lon(:), Lat(:)), size(Lat)) ~= 0;
        ind = ind & Tag1;
    end

    if land_only
        % select land areas where tag is equal to 0
        Tag2 = reshape(event_tag(Lon(:), Lat(:)), size(Lat)) == 0;
        ind = ind & Tag2;
    end

    for t = 1 : length(precip(1, 1, :))
        temp_precip = precip(:, :, t);
        % set the precip of the unselected grid points to a small value (presumably zero)
        temp_precip(~ind) = small_value;
        precip(:, :, t) = temp_precip;
    end



