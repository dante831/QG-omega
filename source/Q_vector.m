function [A, B] = Q_vector(level, ug, vg, T, Phi, event_timespan, dphi, dlambda, f0, beta, GEO_T)

    % the Q-vector form following Dostalek, Schubert, and DeMaria, 2017
                    
    R = 6371000.0;
    Ra = 287.04;

    [A, B, Qx, Qy] = deal(zeros(size(ug)));

    for t = 1 : length(event_timespan)
 
        % term B
        dvg_dp = d_dp(vg(:, :, :, t), level);
        for j = 1 : size(ug, 1)
            B(j, :, :, t) = f0 * beta * dvg_dp(j, :, :);
        end 
       
        if GEO_T
            dug_dp = d_dp(ug(:, :, :, t), level);
        end
    
        for k = 1 : length(level)

            dug_dlambda = d_dlambda(ug(:, :, k, t), dlambda);
            dvg_dlambda = d_dlambda(vg(:, :, k, t), dlambda);
            dug_dphi = d_dphi(ug(:, :, k, t), dphi);
            dvg_dphi = d_dphi(vg(:, :, k, t), dphi);
           
            if GEO_T
                % geostrophically balanced temperature field
                dT_dlambda = f0 * level(k) ./ Ra * (-dvg_dp(:, :, k)) * R .* cos(Phi);
                dT_dphi    = f0 * level(k) ./ Ra *   dug_dp(:, :, k)  * R;
            else
                % original formulation that uses actual T field
                dT_dlambda = d_dlambda(T(:, :, k, t), dlambda);
                dT_dphi    = d_dphi(T(:, :, k, t), dphi);
            end

            Qx(:, :, k, t) = 1 / R^2 * (...
                dT_dphi .* (  1./cos(Phi) .* dvg_dlambda + ug(:, :, k, t).*tan(Phi)) ...
                - 1./cos(Phi) .* dT_dlambda .* dvg_dphi);
            Qy(:, :, k, t) = 1 / R^2 * (...
                dT_dphi .* (- 1./cos(Phi) .* dug_dlambda + vg(:, :, k, t).*tan(Phi)) ...
                + 1./cos(Phi) .* dT_dlambda .* dug_dphi);
                            
            div_temp = div(Qx(:, :, k, t), Qy(:, :, k, t), Phi, dphi, dlambda);
            A(:, :, k, t) = 2 * Ra / level(k) * div_temp;
        
        end
        clear('div_temp');
    end
        
end 
 
