function [A, B] = Q_vector(level, ug, vg, T, Phi, event_timespan, dphi, dlambda, f0, beta, GEO_T)
                    
    R = 6371000.0;
    Ra = 287.04;

    [A, B, Qx, Qy, ...
     dug_dlambda, dvg_dlambda, dug_dphi, dvg_dphi, ...
     dT_dlambda, dT_dphi] = deal(zeros(size(ug)));

    for t = 1 : length(event_timespan)
        
        if GEO_T
            dvg_dp = d_dp(vg(:, :, :, t), level);
            dug_dp = d_dp(ug(:, :, :, t), level);
        end
    
        for k = 1 : length(level)

            dug_dlambda(:, :, k, t) = d_dlambda(ug(:, :, k, t), dlambda);
            dvg_dlambda(:, :, k, t) = d_dlambda(vg(:, :, k, t), dlambda);
            dug_dphi(:, :, k, t) = d_dphi(ug(:, :, k, t), dphi);
            dvg_dphi(:, :, k, t) = d_dphi(vg(:, :, k, t), dphi);
           
            if GEO_T
                % geostrophically balanced temperature field
                dT_dlambda(:, :, k, t) = f0 * level(k) ./ Ra * (-dvg_dp(:, :, k)) * R .* cos(Phi);
                dT_dphi(:, :, k, t)    = f0 * level(k) ./ Ra *   dug_dp(:, :, k)  * R;
            else
                % original formulation that uses actual T field
                dT_dlambda(:, :, k, t) = d_dlambda(T(:, :, k, t), dlambda);
                dT_dphi(:, :, k, t) = d_dphi(T(:, :, k, t), dphi);
            end


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
        
end 
%dT_dlambda(:, :, k, t) = d_dlambda(T(:, :, k, t), dlambda);
%dT_dphi(:, :, k, t) = d_dphi(T(:, :, k, t), dphi);
 
