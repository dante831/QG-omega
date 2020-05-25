function [k2, l2, m2, sigma] = k2_l2_m2_grid( ...
        events_sta, ...
        lat_series, lon_series, box_max_x, box_max_y, ...
        sigma_tag, plevels, input_path, SIGMA_2)

Omega = 7.2921e-5;
Ra = 287.04;
cp = 1005.0;
kappa = Ra/cp;

lat = [];
for j = 1 : length(events_sta.events_sta(:, 1))
    max_i = length(events_sta.events_sta(1, :));
    i = 1;
    while i <= max_i && isempty(events_sta.events_sta{j, i})
        i = i + 1;
    end
    if i <= max_i
        lat(j) = events_sta.events_sta{j, i}(1).event_lat;
    else
        lat(j) = 2 * lat(j - 1) - lat(j - 2);
    end
end
[~, ~, lat_indices] = intersect(lat, lat_series);

[k2, l2, m2] = deal(zeros([size(events_sta.events_sta), length(plevels)]));
[k2_center, sigma, m2_center, l2_center, sigma_center] = deal(zeros([size(events_sta.events_sta), length(plevels)]));

[N_lat, N_lon] = size(events_sta.events_sta);
count_skip = 0;

for j = 1 : N_lat

    for i = 1 : N_lon

        N = length(events_sta.events_sta{j, i});
        if N == 0
            [k2(j, i, :), l2(j, i, :), m2(j, i, :), sigma(j, i, :), ...
             k2_center(j, i, :), ...
             m2_center(j, i, :), sigma_center(j, i, :)...
             l2_center(j, i, :)] = deal(NaN);
            continue;
        end
        event_latind = find(lat_series == events_sta.events_sta{j, i}(1).event_lat);
        lat0 = lat_series(event_latind - (box_max_y - 1) / 2 : event_latind + (box_max_y - 1) / 2);
        lon0 = lon_series(1 : box_max_x);
        [omega_comp, J_term, sig_om_comp] = deal(zeros(box_max_y, box_max_x));
        ind = zeros(size(plevels));
        f0_0 = 2 * Omega * sin(lat_series(lat_indices(j)) / 180 * 3.1415926);

        for n = 1 : N

            level = events_sta.events_sta{j, i}(n).event_level;
            lat = lat_series(events_sta.events_sta{j, i}(n).event_latspan);
            lon = lon_series(events_sta.events_sta{j, i}(n).event_lonspan);
            nlat = length(lat);
            nlon = length(lon);
            [~, map, map2] = intersect(level, plevels);
            map = sort(map);
            map2 = sort(map2);
            
            [sigma_event, k2_event, l2_event, m2_event] = deal(nan([size(plevels), 1]));
            if sigma_tag
                sigma_event(map2) = events_sta.events_sta{j, i}(n).sigma_accu(map);
            elseif SIGMA_2
                temp_sigma = events_sta.events_sta{j, i}(n).sigma(map) + ...
                             Ra ./ level(map) .* events_sta.events_sta{j, i}(n).dtheta_dp_ma_avg(map);
                temp_sigma(temp_sigma < 0) = events_sta.events_sta{j, i}(n).sigma(map(temp_sigma < 0)) * 0.05;
                sigma_event(map2) = temp_sigma;
            else
                sigma_event(map2) = events_sta.events_sta{j, i}(n).sigma(map);
            end

            k2_event(map2)        = events_sta.events_sta{j, i}(n).k2(map);
            l2_event(map2)        = events_sta.events_sta{j, i}(n).l2(map);
            m2_event(map2)        = events_sta.events_sta{j, i}(n).m2(map);
            
            l2_center(j, i, map2) = squeeze(l2_center(j, i, map2)) + l2_event(map2);
            m2_center(j, i, map2) = squeeze(m2_center(j, i, map2)) + m2_event(map2);
            k2_center(j, i, map2) = squeeze(k2_center(j, i, map2)) + k2_event(map2);
            sigma_center(j, i, map2) = squeeze(sigma_center(j, i, map2)) + mean(sigma_event(map2), 2);

            ind(map2) = ind(map2) + 1;
        
        end
        l2_center(j, i, :)              = squeeze(l2_center(j, i, :)) ./ ind;
        k2_center(j, i, :)              = squeeze(k2_center(j, i, :)) ./ ind;
        sigma_center(j, i, :)           = squeeze(sigma_center(j, i, :)) ./ ind;
        m2_center(j, i, :)              = squeeze(m2_center(j, i, :)) ./ ind;

        k2(j, i, :)         = k2_center(j, i, :);
        l2(j, i, :)         = l2_center(j, i, :);
        m2(j, i, :)         = squeeze(m2_center(j, i, :));
        sigma(j, i, :)      = squeeze(sigma_center(j, i, :));

    end
end


