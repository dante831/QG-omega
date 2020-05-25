function [A, B] = Adv_traditional(level, ug, vg, theta, Phi, event_timespan, dphi, dlambda, f0, beta)
                    
%{
    R = 6371000.0;
    Ra = 287.04;

    [A, B, Qx, Qy, ...
     dug_dlambda, dvg_dlambda, dug_dphi, dvg_dphi, ...
     dT_dlambda, dT_dphi] = deal(zeros(size(ug)));

    for t = 1 : length(event_timespan)
    
        for k = 1 : length(level)

            dug_dlambda(:, :, k, t) = d_dlambda(ug(:, :, k, t), dlambda);
            dvg_dlambda(:, :, k, t) = d_dlambda(vg(:, :, k, t), dlambda);
            dug_dphi(:, :, k, t) = d_dphi(ug(:, :, k, t), dphi);
            dvg_dphi(:, :, k, t) = d_dphi(vg(:, :, k, t), dphi);

            dT_dlambda(:, :, k, t) = d_dlambda(T(:, :, k, t), dlambda);
            dT_dphi(:, :, k, t) = d_dphi(T(:, :, k, t), dphi);
            Qx(:, :, k, t) = 1 / R^2 * (...
                (1 ./ cos(Phi).^2) .* dug_dlambda(:, :, k, t) .* dT_dlambda(:, :, k, t) + ...
                (1 ./ cos(Phi))    .* dvg_dlambda(:, :, k, t) .* dT_dphi(:, :, k, t));
            Qy(:, :, k, t) = 1 / R^2 * (...
                (1 ./ cos(Phi))    .* dug_dphi(:, :, k, t)    .* dT_dlambda(:, :, k, t) + ...
                                      dvg_dphi(:, :, k, t)    .* dT_dphi(:, :, k, t));
                            
            div_temp = div(Qx(:, :, k, t), Qy(:, :, k, t), Phi, dphi, dlambda);
            A(:, :, k, t) = 2 * Ra / level(k) * div_temp;
        
        end
        clear('temp_x', 'temp_y', 'div_temp');

        % term B
    
        temp = d_dp(vg(:, :, :, t), level);
        for j = 1 : size(ug, 1)
            B(j, :, :, t) = f0 * beta * temp(j, :, :);
        end 

        clear('temp');
    end
%}

    R = 6371000.0;
    Ra = 287.04;
    
    [A, B] = deal(zeros(size(ug)));

    for t = 1 : length(event_timespan)

        % term A

        temp = zeros(size(A(:, :, :, t)));
        for k = 1 : length(level)
             temp(:, :, k) = v_del(phi, ug(:, :, k, t), vg(:, :, k, t), ...
                    (1 / f0 * spherical_laplacian(z(:, :, k, t), phi, dphi, dlambda) + 2 * Omega * sin(Phi)), ...
                    dphi, dlambda);
        end
                    
        d_dptemp = d_dp(temp(:, :, :), level);
        for k = 1 : length(level)
            A(:, :, k, t) = f0 * d_dptemp(:, :, k);
        end
        clear('temp', 'd_dptemp');

                    
        % term B

        
        temp2 = d_dp(z(:, :, :, t), level);
        for k = 1 : length(level)
            B(:, :, k, t) = spherical_laplacian(...
                    v_del(phi, ug(:, :, k, t), vg(:, :, k, t), - temp2(:, :, k), dphi, dlambda), ...
                    phi, dphi, dlambda);
        end
    end


end 


