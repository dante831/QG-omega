function [A, A1, A2, A3] = traditional_A(phi, level, u, v, zeta, event_timespan, dphi, dlambda, f0, beta)

    [A, A1, A2, A3, u_del_dzeta_dp, du_dp_del_zeta, u_del_zeta] = deal(zeros(size(u)));

    for t = 1 : length(event_timespan)
        
        du_dp = d_dp(u(:, :, :, t), level);
        dv_dp = d_dp(v(:, :, :, t), level);
        dzeta_dp = d_dp(zeta(:, :, :, t), level);

        for k = 1 : length(level)
            
            u_del_dzeta_dp(:, :, k, t) = v_del(phi, u(:, :, k, t), v(:, :, k, t), dzeta_dp(:, :, k), dphi, dlambda);
            du_dp_del_zeta(:, :, k, t) = v_del(phi, du_dp(:, :, k), dv_dp(:, :, k), zeta(:, :, k, t), dphi, dlambda);
            u_del_zeta    (:, :, k, t) = v_del(phi, u(:, :, k, t), v(:, :, k, t), zeta(:, :, k, t), dphi, dlambda); 
        
        end
        
        A(:, :, :, t) = f0 * d_dp(u_del_zeta(:, :, :, t), level) + f0 * beta * d_dp(v(:, :, :, t), level);
        A1(:, :, :, t) = f0 * du_dp_del_zeta(:, :, :, t);
        A2(:, :, :, t) = f0 * u_del_dzeta_dp(:, :, :, t);
        A3(:, :, :, t) = f0 * beta * d_dp(v(:, :, :, t), level);
    
    end


