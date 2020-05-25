function [events_sta, num_event] = remove_by_Adv(events_sta, num_event, plot_level, threshold, string_1, input_path, SIGMA_2)

    N = 0;
    Ra = 287.04;
    cpd= 1005.7;
    kappa = Ra / cpd;


    for j = 1 : length(events_sta.events_sta(:))
        
        ind = true(size(events_sta.events_sta{j}));
        for n = 1 : length(events_sta.events_sta{j})
            Adv = events_sta.events_sta{j}(n).A + events_sta.events_sta{j}(n).B;
            level = events_sta.events_sta{j}(n).event_level;
            temp_T = mean(events_sta.events_sta{j}(n).T, 2);
            J = events_sta.events_sta{j}(n).J_center;
            sigma = mean(events_sta.events_sta{j}(n).sigma, 2);
            [~, temp, temp_theta] = compute_moist_adia_sigma(temp_T, level);
            temp_dtheta_dp_ma = temp_T ./ temp_theta .* temp;
            temp_omega = mean(events_sta.events_sta{j}(n).omega, 2);
            temp_omega_QG   = mean(events_sta.events_sta{j}(n).omega_QG, 2);
            l2 = events_sta.events_sta{j}(n).l2;
            k2 = events_sta.events_sta{j}(n).k2;

            if SIGMA_2 && isempty(events_sta.events_sta{j}(n).dtheta_dp_ma_avg)
                disp('error: missing dtheta_dp_ma_avg, set to zeros');
                dtheta_dp_ma_avg = zeros(size(l2));
            else
                dtheta_dp_ma_avg = events_sta.events_sta{j}(n).dtheta_dp_ma_avg;
            end


            %this sigma_2 should be very close to events_sta.events_sta{j, i}(n).sigma
            %sigma_2 = - Ra * temp_T_avg ./ (level .* theta_avg) .* gradient(theta_avg, level); 
            

 
            higher_order = l2 .* kappa ./ level .* J - Ra ./ level .* temp_dtheta_dp_ma .* l2 .* temp_omega_QG;
            sigma_m = sigma + Ra ./ level .* temp_dtheta_dp_ma .* l2 ./ k2;
 
            if SIGMA_2
                higher_order = l2 .* kappa ./ level .* J - Ra ./ level .* dtheta_dp_ma_avg .* l2 .* temp_omega_QG;
                sigma_m = sigma + Ra ./ level .* dtheta_dp_ma_avg .* l2 ./ k2;
            end

            %higher_order = l2 .* kappa ./ level .* J - Ra ./ level .* temp_dtheta_dp_ma .* k2 .* temp_omega_QG;
            %sigma_m = sigma + Ra ./ level .* temp_dtheta_dp_ma;

            iii = find(level == plot_level);

            if ~SIGMA_2 && ...
               (sum(Adv) < threshold || ...
               Adv(level == plot_level) + higher_order(level == plot_level) < threshold || ...
               Adv(level == plot_level) < threshold)
                ind(n) = false;
                disp(['j = ', num2str(j), ', n = ', num2str(n), ' was removed']);
                N = N + 1;
            else
                events_sta.events_sta{j}(n).dtheta_dp_ma_avg = dtheta_dp_ma_avg;
            end
        end

        events_sta.events_sta{j} = events_sta.events_sta{j}(ind);
        num_event.num_event(j) = sum(double(ind));
    
    end

    disp(['A total of ', num2str(N), ' events were/was removed.']);

end

