function [k2, k2_star, l2, m2, dtheta_dp_ma, dtheta_dp_ma_avg, ...
        sigma, sigma_star, omega_sigma_center, omega_k2_sigma_center, ...
        J_center, J_l2_center, d2_dp2_omega_QG, omega_QG, omega, ...
        Adv, C, A1, A2, A3, B, num_event, k2_m, T_adv, precip, ...
        omega_500hPa_std, omega_QG_500hPa_std] = comp_composite_k2_l2_m2_grid_v1( ...
        nt, events_sta, num_event, ...
        lat_series, lon_series, box_max_x, box_max_y, ...
        sigma_tag, plevels, string_1, input_path, SIGMA_2, VORTICITY, TRENBERTH)

Omega = 7.2921e-5;
Ra = 287.04;
cp = 1005.0;
kappa = Ra/cp;

[latmax, latmin, lonmax, lonmin] = get_lat_lon_limits(pwd);
[lat_indices, lon_indices] = latlonindices(lat_series, lon_series, latmin, latmax, lonmin, lonmax);

[k2, k2_m, k2_star, l2, m2, dtheta_dp_ma, dtheta_dp_ma_omega, d2_dp2_omega_QG, ...
    omega_QG, omega] = deal(zeros([size(events_sta.events_sta), length(plevels)]));
[sigma, sigma_star, omega_k2_sigma_center, omega_sigma_center, omega_QG_m2_center, ...
    J_l2_center, J_center, omegaQG_sigmastar_center, omegaQG_k2_sigmastar_center] = ...
    deal(zeros([size(events_sta.events_sta), length(plevels)]));
[Adv, C, A1, A2, A3, B, T_adv] = deal(zeros([size(events_sta.events_sta), length(plevels)]));
[precip, omega_500hPa_std, omega_QG_500hPa_std] = deal(zeros(size(events_sta.events_sta)));

[N_lat, N_lon] = size(events_sta.events_sta);
count_skip = 0;

for j = 1 : N_lat

    for i = 1 : N_lon

        N = length(events_sta.events_sta{j, i});
        
        if N == 0

            [k2(j, i, :), k2_m(j, i, :), k2_star(j, i, :), l2(j, i, :), m2(j, i, :), ...
             dtheta_dp_ma(j, i, :), dtheta_dp_ma_omega(j, i, :), ...
             sigma(j, i, :), sigma_star(j, i, :), d2_dp2_omega_QG(j, i, :), ...
             omega_QG(j, i, :), omega(j, i, :), precip(j, i), ...
             omega_QG_500hPa_std(j, i), omega_500hPa_std(j, i), ...
             omega_k2_sigma_center(j, i, :), omega_sigma_center(j, i, :), ...
             omega_QG_m2_center(j, i, :), ...
             J_l2_center(j, i, :), J_center(j, i, :)] = deal(NaN);

            [Adv(j, i, :), C(j, i, :), ...
             A1(j, i, :), A2(j, i, :), A3(j, i, :), B(j, i, :), ...
             T_adv(j, i, :)] = deal(NaN);

            continue;
        end

        event_latind = find(lat_series == events_sta.events_sta{j, i}(1).event_lat);
        lat0 = lat_series(event_latind - (box_max_y - 1) / 2 : event_latind + (box_max_y - 1) / 2);
        lon0 = lon_series(1 : box_max_x);
        [omega_comp, J_term, sig_om_comp] = deal(zeros(box_max_y, box_max_x));
        ind = zeros(size(plevels));
        [T, omega_QG_center, omega_center, ...
         dtheta_dp_ma_event, dtheta_dp_ma_omega_event, dtheta_dp_ma_avg_event, ...
         Adv_center, C_center, A1_center, A2_center, A3_center, B_center, ...
         k2_m_sigma_m_center, sigma_m_center, ...
         T_adv_center] = deal(zeros(size(plevels)));
        [omega_QG_500hPa_temp, omega_500hPa_temp] = deal(zeros(N, 1));
        precip_center = 0;
        f0_0 = 2 * Omega * sin(lat_series(lat_indices(j)) / 180 * 3.1415926);

        for n = 1 : N

            level = events_sta.events_sta{j, i}(n).event_level;
            
            if VORTICITY
                if exist('TRENBERTH') && TRENBERTH
                    aaa = events_sta.events_sta{j, i}(n).A;
                else
                    aaa = max([events_sta.events_sta{j, i}(n).A1, ...
                               events_sta.events_sta{j, i}(n).A2, ...
                               events_sta.events_sta{j, i}(n).A3, ...
                               events_sta.events_sta{j, i}(n).B]')';
                end
                m2_temp = events_sta.events_sta{j, i}(n).m2;
                if abs(aaa(level == 50000)) > 1e-16 || m2_temp(level == 50000) < 0
                    %disp(['skipped event jj = ', num2str(j), ' ii = ', num2str(i), ' n = ', num2str(n)])
                    count_skip = count_skip + 1;
                    num_event.num_event(j, i) = num_event.num_event(j, i) - 1;
                    continue;
                end
            end
            lat = lat_series(events_sta.events_sta{j, i}(n).event_latspan);
            lon = lon_series(events_sta.events_sta{j, i}(n).event_lonspan);
            nlat = length(lat);
            nlon = length(lon);
            [~, map, map2] = intersect(level, plevels);
            map = sort(map);
            map2 = sort(map2);
            
            [sigma_event, sigma_star_event, k2_event, k2_m_event, l2_event, k2_star_event, m2_event, temp_omega_1, J_event, ...
             temp_omega, temp_omega_QG, temp_Adv, temp_T, temp_T_avg, ...
             temp_dtheta_dp_ma_avg, temp_C, ...
             temp_A1, temp_A2, temp_A3, temp_B, ...
             temp_k2_m_sigma_m, temp_sigma_m, ...
             temp_T_adv] = deal(nan([size(plevels), 1]));
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
            if ~any(~isnan(l2_event)) 
                % if all elements in l2_event is nan
                % this is a result of inversion without J
                l2_event(map2) = 0;
            end
            m2_event(map2)        = events_sta.events_sta{j, i}(n).m2(map);
            J_event(map2)         = events_sta.events_sta{j, i}(n).J_center(map);
            temp_omega(map2)      = events_sta.events_sta{j, i}(n).omega(map);
            temp_omega_QG(map2)   = events_sta.events_sta{j, i}(n).omega_QG(map);
            temp_sigma_m(map2) = sigma_event(map2) .* temp_omega_QG(map2) + kappa ./ plevels(map2) .* J_event(map2);
            
            temp_k2_m_sigma_m(map2) = k2_event(map2) .* sigma_event(map2) .* temp_omega_QG(map2) + l2_event(map2) .* kappa ./ plevels(map2) .* J_event(map2);
            temp_T_adv(map2) = J_event(map2) + sigma_event(map2) .* temp_omega(map2) .* plevels(map2) ./ kappa;
            if exist('TRENBERTH') && TRENBERTH
                temp_Adv(map2)        = mean(events_sta.events_sta{j, i}(n).A       (map, :), 2);
            else
                temp_Adv(map2)        = mean(events_sta.events_sta{j, i}(n).A       (map, :), 2) + ...
                                        mean(events_sta.events_sta{j, i}(n).B       (map, :), 2);
            end
            if ~isempty(events_sta.events_sta{j, i}(n).A1)
                temp_A1(map2)         = mean(events_sta.events_sta{j, i}(n).A1      (map, :), 2);
                temp_A2(map2)         = mean(events_sta.events_sta{j, i}(n).A2      (map, :), 2);
                temp_A3(map2)         = mean(events_sta.events_sta{j, i}(n).A3      (map, :), 2);
            end
            temp_B(map2)          = mean(events_sta.events_sta{j, i}(n).B       (map, :), 2);
            temp_C(map2)          = mean(events_sta.events_sta{j, i}(n).C       (map, :), 2);
            temp_T(map2)          = mean(events_sta.events_sta{j, i}(n).T       (map, :), 2);
            if ~isempty(events_sta.events_sta{j, i}(n).dtheta_dp_ma_avg)
                temp_dtheta_dp_ma_avg(map2) = mean(events_sta.events_sta{j, i}(n).dtheta_dp_ma_avg(map, :), 2);
            end
            
            % check inversion accuracy
            ind_p = plevels == 50000;
            tag = false;
            if SIGMA_2
                temp = (k2_event(ind_p)*sigma_event(ind_p) + f0_0^2*m2_event(ind_p))*temp_omega_QG(ind_p) + temp_Adv(ind_p);
            elseif VORTICITY
                temp = f0_0^2*m2_event(ind_p)*temp_omega_QG(ind_p) + temp_Adv(ind_p);
            else
                %temp = (k2_event(ind_p)*sigma_event(ind_p) + f0_0^2*m2_event(ind_p))*temp_omega_QG(ind_p) + temp_Adv(ind_p) + temp_C(ind_p);
                temp = (k2_event(ind_p)*sigma_event(ind_p) + f0_0^2*m2_event(ind_p))*temp_omega_QG(ind_p) + ...
                        temp_Adv(ind_p) + l2_event(ind_p).*kappa./plevels(ind_p).*J_event(ind_p);
            end
            if abs(temp) < 1e-28
                tag = true;
            end
            
            if ~tag
                disp(['warning, equation accuracy not satisfied, jj = ', num2str(j), '; ii = ', num2str(i), ...
                        ' temp = ', num2str(temp)]);
            end

            [~, temp, temp_theta] = compute_moist_adia_sigma(temp_T, plevels);
            temp_dtheta_dp_ma = temp_T ./ temp_theta .* temp;
            if ~isempty(events_sta.events_sta{j, i}(n).k2_star)
                K_STAR = true;
                k2_star_event(map2) = events_sta.events_sta{j, i}(n).k2_star(map);
            else
                K_STAR = false;
            end
            
            % calculate sigma_star if not saved in event_analysis_v1.m
            %if isempty(events_sta.events_sta{j, i}(n).sigma_star)
                sigma_star_event(:) = - temp_dtheta_dp_ma * cp * kappa ./ plevels;
            %else
            %    sigma_star_event(map2) = events_sta.events_sta{j, i}(n).sigma_star(map);
            %end

            % get omega and omega_QG at 500hPa for std
            omega_QG_500hPa_temp(n) = events_sta.events_sta{j, i}(n).omega_QG(level == 50000);
            omega_500hPa_temp   (n) = events_sta.events_sta{j, i}(n).omega   (level == 50000);

            dtheta_dp_ma_event(map2) = dtheta_dp_ma_event(map2) + squeeze(temp_dtheta_dp_ma(map2));
            dtheta_dp_ma_omega_event(map2) = dtheta_dp_ma_omega_event(map2) + temp_dtheta_dp_ma(map2) .* temp_omega(map2);
            dtheta_dp_ma_avg_event(map2) = dtheta_dp_ma_avg_event(map2) + squeeze(temp_dtheta_dp_ma_avg(map2));
            T              (map2) = T              (map2) + squeeze(temp_T(map2));
            omega_center   (map2) = omega_center   (map2) + squeeze(temp_omega(map2));
            omega_QG_center(map2) = omega_QG_center(map2) + squeeze(temp_omega_QG(map2));
            precip_center         = precip_center + events_sta.events_sta{j, i}(n).precip;
            Adv_center     (map2) = Adv_center     (map2) + squeeze(temp_Adv(map2));
            A1_center      (map2) = A1_center      (map2) + squeeze(temp_A1(map2));
            A2_center      (map2) = A2_center      (map2) + squeeze(temp_A2(map2));
            A3_center      (map2) = A3_center      (map2) + squeeze(temp_A3(map2));
            B_center       (map2) = B_center       (map2) + squeeze(temp_B(map2));
            C_center       (map2) = C_center       (map2) + squeeze(temp_C(map2));
            sigma_m_center (map2) = sigma_m_center (map2) + squeeze(temp_sigma_m(map2));
            k2_m_sigma_m_center(map2) = k2_m_sigma_m_center(map2) + squeeze(temp_k2_m_sigma_m(map2));
            T_adv_center   (map2) = T_adv_center   (map2) + squeeze(temp_T_adv(map2));
            J_l2_center       (j, i, map2) = squeeze(J_l2_center       (j, i, map2)) + J_event(map2) .* l2_event(map2);
            J_center          (j, i, map2) = squeeze(J_center          (j, i, map2)) + J_event(map2);
            omega_QG_m2_center(j, i, map2) = squeeze(omega_QG_m2_center(j, i, map2)) + temp_omega_QG(map2) .* m2_event(map2);
            
            omega_k2_sigma_center(j, i, map2) = squeeze(omega_k2_sigma_center(j, i, map2)) + ...
                        temp_omega_QG(map2) .* sigma_event(map2) .* k2_event(map2);
            omega_sigma_center(j, i, map2) = squeeze(omega_sigma_center(j, i, map2))    + mean(temp_omega_QG(map2) .* sigma_event(map2), 2);

            omegaQG_sigmastar_center(j, i, map2) = squeeze(omegaQG_sigmastar_center(j, i, map2)) + ...
                            temp_omega_QG(map2) .* sigma_star_event(map2);
            omegaQG_k2_sigmastar_center(j, i, map2) = squeeze(omegaQG_k2_sigmastar_center(j, i, map2)) + ...
                            temp_omega_QG(map2) .* sigma_star_event(map2) .* k2_star_event(map2);
            
            ind(map2) = ind(map2) + 1;
        
        end
        if isnan(omega_QG_center(plevels == 50000))
            disp('warning: NaN in omega_QG_center(j=', num2tr(j), ', i=', num2str(i), ') at 500hPa')
        end
        dtheta_dp_ma_event = dtheta_dp_ma_event ./ ind;
        dtheta_dp_ma_omega_event = dtheta_dp_ma_omega_event ./ ind;
        dtheta_dp_ma_avg_event = dtheta_dp_ma_avg_event ./ ind;
        omega_center    = omega_center    ./ ind;
        omega_QG_center = omega_QG_center ./ ind;
        Adv_center      = Adv_center      ./ ind;
        T               = T               ./ ind;
        A1_center       = A1_center       ./ ind;
        A2_center       = A2_center       ./ ind;
        A3_center       = A3_center       ./ ind;
        B_center        = B_center        ./ ind;
        C_center        = C_center        ./ ind;
        sigma_m_center  = sigma_m_center  ./ ind;
        precip_center   = precip_center    / ind(plevels == 50000);
        k2_m_sigma_m_center = k2_m_sigma_m_center ./ ind;
        omegaQG_sigmastar_center(j, i, :) = squeeze(omegaQG_sigmastar_center(j, i, :)) ./ ind;
        omegaQG_k2_sigmastar_center(j, i, :) = squeeze(omegaQG_k2_sigmastar_center(j, i, :)) ./ ind;
        T_adv_center    = T_adv_center    ./ ind;
        J_l2_center(j, i, :)            = squeeze(J_l2_center(j, i, :)) ./ ind;
        J_center(j, i, :)               = squeeze(J_center(j, i, :)) ./ ind;
        omega_k2_sigma_center(j, i, :)  = squeeze(omega_k2_sigma_center(j, i, :)) ./ ind;
        omega_sigma_center(j, i, :)     = squeeze(omega_sigma_center(j, i, :)) ./ ind;
        omega_QG_m2_center(j, i, :)     = squeeze(omega_QG_m2_center(j, i, :)) ./ ind;

        k2(j, i, :) = omega_k2_sigma_center(j, i, :) ./ omega_sigma_center(j, i, :);
        k2_m(j, i, :) = k2_m_sigma_m_center ./ sigma_m_center;
        k2_star(j, i, :) = omegaQG_k2_sigmastar_center(j, i, :) ./ omegaQG_sigmastar_center(j, i, :);
        l2(j, i, :) = J_l2_center(j, i, :) ./ J_center(j, i, :);
 
        d2_dp2_omega_QG_temp = d2_dp2(omega_QG_center, plevels);
        m2(j, i, :) = squeeze(omega_QG_m2_center(j, i, :)) ./ omega_QG_center;
        d2_dp2_omega_QG(j, i, :) = d2_dp2_omega_QG_temp;

        omega_QG(j, i, :)   = omega_QG_center;
        omega(j, i, :)      = omega_center;
        precip(j, i)        = precip_center;
        sigma(j, i, :)      = squeeze(omega_sigma_center(j, i, :) ./ omega_QG(j, i, :));
        sigma_star(j, i, :) = squeeze(omegaQG_sigmastar_center(j, i, :) ./ omega_QG(j, i, :));

        Adv(j, i, :) = Adv_center;
        A1(j, i, :) = A1_center;
        A2(j, i, :) = A2_center;
        A3(j, i, :) = A3_center;
        B(j, i, :) = B_center;
        C(j, i, :) = C_center;
        T_adv(j, i, :) = T_adv_center;
        [~, temp, theta] = compute_moist_adia_sigma(T, plevels);

        % calculating moist adiabatic stability using the composite field
        dtheta_dp_ma(j, i, :) = dtheta_dp_ma_omega_event ./ omega_center;
        dtheta_dp_ma_avg(j, i, :) = dtheta_dp_ma_avg_event;

        % calculate standard deviation of omega and omega_QG
        omega_500hPa_std   (j, i) = std(omega_500hPa_temp);
        omega_QG_500hPa_std(j, i) = std(omega_QG_500hPa_temp);

    end

end

if VORTICITY
    disp(['a total of ', num2str(count_skip), ' events were skipped']);
end
