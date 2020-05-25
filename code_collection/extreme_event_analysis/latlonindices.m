function [lat_indices, lon_indices] = latlonindices(lat_series, lon_series, latmin, latmax, lonmin, lonmax)
    
    lat_indices = find(lat_series >= latmin & lat_series <= latmax);

    if lonmax > lonmin
        lon_indices = find(lon_series >= lonmin & lon_series <= lonmax);
    else
        lon_indices = [find(lon_series >= lonmin)', find(lon_series <= lonmax)']';
    end

end
