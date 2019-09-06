function B = traditional_B(phi, level, u, v, T, event_timespan, dphi, dlambda, f0, GEO_T)

    B = zeros(size(u));
    Ra = 287.04;
    
    % use the geostrophically balanced temperature
    for t = 1 : length(event_timespan)
        if GEO_T
            dvg_dp = d_dp(v(:, :, :, t), level);
            dug_dp = d_dp(u(:, :, :, t), level);
        end
        for k = 1 : length(level)
            if GEO_T
                u_del_T = u(:, :, k, t) .* f0 * level(k) ./ Ra .* (-dvg_dp(:, :, k)) + ...
                          v(:, :, k, t) .* f0 * level(k) ./ Ra .*   dug_dp(:, :, k);
            else
                u_del_T = v_del(phi, u(:, :, k, t), v(:, :, k, t), T(:, :, k, t), dphi, dlambda);
            end
            B(:, :, k, t) = Ra ./ level(k) .* spherical_laplacian(u_del_T, phi, dphi, dlambda);
        end
    end

   

