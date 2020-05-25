function [event_latspan, event_lonspan, center_y, computation_y] = ...
                        event_span_Wrapper(lat_series_original, lon_series_original, event_lat, event_lon, ...
                        box_max_y, computation_x, computation_y)

	lon1 = lon_series_original(mod(event_lon - (computation_x - 1) / 2 - 1, length(lon_series_original)) + 1); % left longitude
    lon2 = lon_series_original(mod(event_lon + (computation_x - 1) / 2 - 1, length(lon_series_original)) + 1); % right longitude
    if event_lat - (box_max_y - 1) / 2 < 1
        lat1 = lat_series_original(1);
        lat2 = lat_series_original(event_lat + (box_max_y - 1) / 2);
        %center_y = (box_max_y + 1) / 2 - (1 - event_lat + (box_max_y - 1) / 2); % calculate the y index of event center in the event domain
        center_y = event_lat;
        %computation_y = (event_lat + (box_max_y - 1) / 2 - 1 + 1 - center_y) * 2 + 1;
        computation_y = center_y * 2 - 1;   % computation_y is shrinked so that the vertical extent of the event does not exceed the 
                                            % boundaries of the original domain
    elseif event_lat + (box_max_y - 1) / 2 > length(lat_series_original)
        lat1 = lat_series_original(event_lat - (box_max_y - 1) / 2);
        lat2 = lat_series_original(length(lat_series_original));
        center_y = (box_max_y + 1) / 2;
        computation_y = (length(lat_series_original) - event_lat) * 2 + 1;
        %computation_y = (length(lat_series_original) - event_lat + (box_max_y - 1) / 2 + 1 - center_y) * 2 + 1;
    else
        lat1 = lat_series_original(event_lat - (box_max_y - 1) / 2);
        lat2 = lat_series_original(event_lat + (box_max_y - 1) / 2);
        center_y = (box_max_y + 1) / 2;
    end
    [event_latspan, event_lonspan] = latlonindices(lat_series_original, lon_series_original, lat1, lat2, lon1, lon2);

end
