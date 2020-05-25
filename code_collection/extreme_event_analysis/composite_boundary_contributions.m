function [omega_QG_interior, omega_QG_full, omega_QG_lateral_b, omega_QG_lateral_clim_b] = ...
            composite_boundary_contributions(nt, events_sta, num_event, ...
        lat_series, lon_series, box_max_x, box_max_y, plevels, string_1, input_path)
   
[omega_QG_interior, omega_QG_full, omega_QG_lateral_b, omega_QG_lateral_clim_b] = ...
    deal(zeros([size(events_sta.events_sta), length(plevels)]));
[N_lat, N_lon] = size(events_sta.events_sta);

first_event_index = min(find(~cellfun(@isempty, {events_sta.events_sta{:}})));

tags = true(4, 1);
if isempty(events_sta.events_sta{first_event_index}(1).omega_QG_interior)
    tags(1) = false;
end
if isempty(events_sta.events_sta{first_event_index}(1).omega_QG_full)
    tags(2) = false;
end
if isempty(events_sta.events_sta{first_event_index}(1).omega_QG_lateral_b)
    tags(3) = false;
end
if isempty(events_sta.events_sta{first_event_index}(1).omega_QG_lateral_clim_b)
    tags(4) = false;
end


for j = 1 : N_lat        
    for i = 1 : N_lon
        N = length(events_sta.events_sta{j, i});
        if N == 0
            [omega_QG_interior(j, i, :), omega_QG_full(j, i, :), ...
             omega_QG_lateral_b(j, i, :), omega_QG_lateral_clim_b(j, i, :)] = deal(NaN);
            continue;
        end
        ind = zeros(size(plevels));
        [omega_QG_interior_center, omega_QG_full_center, omega_QG_lateral_b_center, ...
         omega_QG_lateral_clim_b_center] = deal(zeros(size(plevels)));
        for n = 1 : N
            level = events_sta.events_sta{j, i}(n).event_level;
            lat = lat_series(events_sta.events_sta{j, i}(n).event_latspan);
            lon = lon_series(events_sta.events_sta{j, i}(n).event_lonspan);
            nlat = length(lat);
            nlon = length(lon);
            [~, map, map2] = intersect(level, plevels);
            map = sort(map);
            map2 = sort(map2);
            [temp_omega_QG_interior, temp_omega_QG_full, temp_omega_QG_lateral_b, ...
             temp_omega_QG_lateral_clim_b] = deal(nan(size(plevels)));
            if tags(1)
                temp_omega_QG_interior(map2)  = mean(events_sta.events_sta{j, i}(n).omega_QG_interior (map, :), 2);
            end
            if tags(2)
                temp_omega_QG_full(map2)      = mean(events_sta.events_sta{j, i}(n).omega_QG_full   (map, :), 2);
            end
            if tags(3)
                temp_omega_QG_lateral_b(map2) = mean(events_sta.events_sta{j, i}(n).omega_QG_lateral_b(map, :), 2);
            end
            if tags(4)
                temp_omega_QG_lateral_clim_b(map2)    = mean(events_sta.events_sta{j, i}(n).omega_QG_lateral_clim_b   (map, :), 2);
            end

            omega_QG_interior_center(map2)  = omega_QG_interior_center(map2)  + squeeze(temp_omega_QG_interior(map2));
            omega_QG_full_center(map2)      = omega_QG_full_center(map2)      + squeeze(temp_omega_QG_full(map2));
            omega_QG_lateral_b_center(map2) = omega_QG_lateral_b_center(map2) + squeeze(temp_omega_QG_lateral_b(map2));
            omega_QG_lateral_clim_b_center(map2)    = omega_QG_lateral_clim_b_center(map2)    + squeeze(temp_omega_QG_lateral_clim_b(map2));

            ind(map2) = ind(map2) + 1;
        end

        %if isnan(omega_QG_interior_center(plevels == 50000))
        %    disp(['Warning: NaN in omega_QG_interior_center(j=', num2str(j), ', i=', num2str(i), ') at 500hPa'])
        %end

        omega_QG_interior(j, i, :)  = omega_QG_interior_center ./ ind;
        omega_QG_full(j, i, :)      = omega_QG_full_center ./ ind;
        omega_QG_lateral_b(j, i, :) = omega_QG_lateral_b_center ./ ind;
        omega_QG_lateral_clim_b(j, i, :)    = omega_QG_lateral_clim_b_center ./ ind;
    
    end
end

end

