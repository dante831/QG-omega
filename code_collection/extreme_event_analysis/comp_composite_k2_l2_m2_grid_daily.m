function [k2, l2, m2, dtheta_dp_ma, dtheta_dp_ma_omega, ...
        omega_sigma_center, omega_k2_sigma_center, ...
        J_center, J_l2_center, d2_dp2_omega_QG_grid, omega_QG_grid, omega_grid] = comp_composite_k2_l2_m2_grid_daily( ...
        nt, events_sta, num_event, plot_level, ...
        lat_series, lon_series, box_max_x, box_max_y, ...
        sigma_tag)

lat = [];
for j = 1 : length(events_sta.events_sta(:, 1))
    i = 1;
    while isempty(events_sta.events_sta{j, i})
        i = i + 1;
    end
    lat(j) = events_sta.events_sta{j, i}(1).event_lat;
end
[~, ~, lat_indices] = intersect(lat, lat_series);
%[lat_indices, lon_indices] = latlonindices(lat_series, lon_series, latmin, latmax, lonmin, lonmax);
%lat = lat_series(lat_indices);
%lon = lon_series(lon_indices);
%dphi = (lat_series(2) - lat_series(1)) / 180.0 * 3.1415926;
%dlambda = (lon_series(2) - lon_series(1)) / 180.0 * 3.1415926;
% compute composite k^2 and l^2

plevels = [85000, 82500, 80000, 77500, 75000, 70000, 65000, ...
           60000, 55000, 50000, 45000, 40000, 35000, 30000, 25000, 22500, 20000, 17500, 15000]';

kappa = 2./7.;

[k2, l2, m2, d2_dp2_omega_QG_grid, ...
    omega_QG, omega_QG_grid, omega_grid] = deal(zeros(size(events_sta.events_sta)));
[dtheta_dp_ma, dtheta_dp_ma_omega] = deal(zeros([size(events_sta.events_sta), 4]));
[omega_k2_sigma_center, omega_sigma_center, J_l2_center, J_center] = deal(zeros(size(events_sta.events_sta)));

[N_lat, N_lon] = size(events_sta.events_sta);

for j = 1 : N_lat

    for i = 1 : N_lon

        %N = num_event.num_event(j, i);
        N = length(events_sta.events_sta{j, i});
        if N == 0
            [k2(j, i), l2(j, i), m2(j, i), dtheta_dp_ma(j, i, :), ...
             dtheta_dp_ma_omega(j, i, :), d2_dp2_omega_QG_grid(j, i), ...
             omega_QG(j, i), omega_QG_grid(j, i), omega_grid(j, i), ...
             omega_k2_sigma_center(j, i), omega_sigma_center(j, i), ...
             J_l2_center(j, i), J_center(j, i)] = deal(NaN);
            continue;
        end
        event_latind = find(lat_series == events_sta.events_sta{j, i}(1).event_lat);
        lat0 = lat_series(event_latind - (box_max_y - 1) / 2 : event_latind + (box_max_y - 1) / 2);
        lon0 = lon_series(1 : box_max_x);
        [omega_comp, J_term, sig_om_comp] = deal(zeros(box_max_y, box_max_x));
        ind = zeros(length(plevels), 1);
        ind_2 = zeros(box_max_y, box_max_x);
        [T, omega_QG_center, omega_center] = deal(zeros(size(plevels)));
        [dtheta_dp_ma_event, dtheta_dp_ma_omega_event] = deal(zeros(length(plevels), 4));
        for n = 1 : N
            %if events_sta.events_sta{j, i}(n).computation_x ~= 9
            %    continue
            %end
            level = events_sta.events_sta{j, i}(n).event_level;
            %if level(1) > 85000
            %    continue
            %end
            lat = lat_series(events_sta.events_sta{j, i}(n).event_latspan);
            lon = lon_series(events_sta.events_sta{j, i}(n).event_lonspan);
            nlat = length(lat);
            nlon = length(lon);
            %if nlat > nlat0
            %    lat0 = lat;
            %    lon0 = lon;
            %end
            %if plot_level == 40000
            %    omega_QG = reshape(events_sta.events_sta{j, i}(n).omega_QG_400hPa_full, nlat, nlon, nt);
            %    temp2 = mean(events_sta.events_sta{j, i}(n).J_400hPa, 3);
            %elseif plot_level == 50000
            %    omega_QG = reshape(events_sta.events_sta{j, i}(n).omega_QG_500hPa_full, nlat, nlon, nt);
            %    temp2 = mean(events_sta.events_sta{j, i}(n).J_500hPa, 3);
            %elseif plot_level == 60000
            %    omega_QG = reshape(events_sta.events_sta{j, i}(n).omega_QG_600hPa_full, nlat, nlon, nt);
            %    temp2 = mean(events_sta.events_sta{j, i}(n).J_600hPa, 3);
            %else
            %    disp(['Error: omega_QG and J field on level ', num2str(plot_level), ' does not exist']);
            %    return;
            %end
            sigma = events_sta.events_sta{j, i}(n).sigma(level == plot_level, :);
            k2_event = events_sta.events_sta{j, i}(n).k2(level == plot_level, :);
            l2_event = events_sta.events_sta{j, i}(n).l2(level == plot_level, :);
            temp_omega_1 = events_sta.events_sta{j, i}(n).omega_QG(level == plot_level, :);
            temp_J_1 = events_sta.events_sta{j, i}(n).J_center(level == plot_level, :);
            [~, map, map2] = intersect(level, plevels);
            map = sort(map);
            map2 = sort(map2);
            temp_omega = mean(events_sta.events_sta{j, i}(n).omega, 2);
            temp_omega_2 = events_sta.events_sta{j, i}(n).omega;
            temp_omega_QG = mean(events_sta.events_sta{j, i}(n).omega_QG, 2);
            temp_dtheta_dp_ma = zeros(length(level), 4);
            for t = 1 : 4
                temp_T = events_sta.events_sta{j, i}(n).T(:, t);
                [~, temp, temp_theta] = compute_moist_adia_sigma(temp_T, level);
                temp_dtheta_dp_ma(:, t) = temp_T ./ temp_theta .* temp;
            end
            %temp_T = mean(events_sta.events_sta{j, i}(n).T, 2);
            %[~, temp, temp_theta] = compute_moist_adia_sigma(temp_T, level);
            %temp_dtheta_dp_ma = temp_T ./ temp_theta .* temp;
            dtheta_dp_ma_omega_event(map2, :) = dtheta_dp_ma_omega_event(map2, :) + ...
                    temp_dtheta_dp_ma(map, :) .* temp_omega_2(map, :);
            dtheta_dp_ma_event(map2, :) = dtheta_dp_ma_event(map2, :) + squeeze(temp_dtheta_dp_ma(map, :));
            T                 (map2) = T                 (map2) + squeeze(temp_T(map));
            omega_center      (map2) = omega_center      (map2) + squeeze(temp_omega(map));
            omega_QG_center   (map2) = omega_QG_center   (map2) + squeeze(temp_omega_QG(map));
            J_l2_center(j, i)        = J_l2_center(j, i)        + mean(temp_J_1 .* l2_event);
            J_center(j, i)           = J_center(j, i)           + mean(temp_J_1);
            if sigma_tag
                omega_k2_sigma_center(j, i) = omega_k2_sigma_center(j, i) + mean(temp_omega_1 .* sigma .* k2_event);
                omega_sigma_center(j, i)    = omega_sigma_center(j, i)    + mean(temp_omega_1 .* sigma);
            else
                omega_k2_sigma_center(j, i) = omega_k2_sigma_center(j, i) + mean(temp_omega_1 .* k2_event);
                omega_sigma_center(j, i)    = omega_sigma_center(j, i)    + mean(temp_omega_1);
            end
            ind(map2) = ind(map2) + 1;
            select_x = (box_max_x - nlon) / 2 + 1: (box_max_x + nlon) / 2;
            select_y = (box_max_y - nlat) / 2 + 1: (box_max_y + nlat) / 2;
            %omega_comp (select_y, select_x) = omega_comp (select_y, select_x) + mean(omega_QG, 3);
            %J_term     (select_y, select_x) = J_term     (select_y, select_x) + temp2;
            %sig_om_comp(select_y, select_x) = sig_om_comp(select_y, select_x) + mean(omega_QG, 3) * sigma;    
            ind_2      (select_y, select_x) = ind_2      (select_y, select_x) + 1;
        
            if isnan(omega_QG_center)
                break
            end
        end 
        dtheta_dp_ma_event = dtheta_dp_ma_event ./ repmat(ind, 1, 4);
        dtheta_dp_ma_omega_event = dtheta_dp_ma_omega_event ./ repmat(ind, 1, 4);
        omega_center    = omega_center    ./ ind;
        omega_QG_center = omega_QG_center ./ ind;
        T               = T               ./ ind;
        %omega_comp      = omega_comp      ./ ind_2;
        %J_term          = J_term          ./ ind_2;
        %sig_om_comp     = sig_om_comp     ./ ind_2;
        J_l2_center(j, i) = J_l2_center(j, i) / N;
        J_center(j, i) = J_center(j, i) / N;
        omega_k2_sigma_center(j, i) = omega_k2_sigma_center(j, i) / N;
        omega_sigma_center(j, i) = omega_sigma_center(j, i) / N;

        lon0 = mod(lon0 - lon0(1), 360);
        phi = lat0 / 180.0 * 3.1415926;
        lambda = lon0 / 180.0 * 3.1415926;
        dphi    = mod(phi(2)    - phi(1)   , 360);
        dlambda = mod(lambda(2) - lambda(1), 360);

        %temp3 = spherical_laplacian(omega_comp, phi, dphi, dlambda);
        %k2(j, i) = - temp3(floor((box_max_y + 1) / 2), floor((box_max_x + 1) / 2)) / omega_comp(floor((box_max_y + 1) / 2), floor((box_max_x + 1) / 2));
        k2(j, i) = omega_k2_sigma_center(j, i) / omega_sigma_center(j, i);
        %temp4 = spherical_laplacian(J_term, phi, dphi, dlambda);
        %l2(j, i) = - temp4(floor((box_max_y + 1) / 2), floor((box_max_x + 1) / 2)) / J_term(floor((box_max_y + 1) / 2), floor((box_max_x + 1) / 2));
        l2(j, i) = J_l2_center(j, i) / J_center(j, i);
 
        d2_dp2_omega_QG = d_dp(d_dp(omega_QG_center, plevels), plevels);
        m2(j, i) = - d2_dp2_omega_QG(plevels == plot_level) ./ omega_QG_center(plevels == plot_level);
        %m2(j, i) = abs(d2_dp2_omega_QG(plevels == plot_level)) ./ abs(omega_QG_center(plevels == plot_level));
        d2_dp2_omega_QG_grid(j, i) = d2_dp2_omega_QG(plevels == plot_level);
        omega_QG_grid(j, i) = omega_QG_center(plevels == plot_level);
        omega_grid(j, i) = omega_center(plevels == plot_level);

        [~, temp, theta] = compute_moist_adia_sigma(T, plevels);

        % calculating moist adiabatic stability using the composite field
        %dtheta_dp_ma(j, i) = T(plevels == plot_level) / theta(plevels == plot_level)...
        %        * temp(plevels == plot_level);
        % A better alternative: calculate the same thing event-wise
        dtheta_dp_ma(j, i, :) = dtheta_dp_ma_event(plevels == plot_level, :);
        dtheta_dp_ma_omega(j, i, :) = dtheta_dp_ma_omega_event(plevels == plot_level, :);


    end
end


