function [C1, C2, J1, J2] = diabatic_heating(phi, T, level, u, v, sigma_accu, omega, event_timespan, dphi, dlambda, dt)

    cp = 1005.0;
    Ra = 287.04;
    kappa = Ra / cp;
    [J1, J2, C1, C2] = deal(zeros(size(u)));
                    
    for t = 1 : length(event_timespan)
    
        if t == 1
            T1 = T(:, :, :, t);
        end
        if t < length(event_timespan)
            T2 = T(:, :, :, t + 1); 
        end

        for k = 1 : length(level)

            % calculate diabatic heating using temperature and wind fields
            if t == length(event_timespan)
                temp = (T1(:, :, k) - T0(:, :, k)) / dt;
            elseif t == 1
                temp = (T2(:, :, k) - T1(:, :, k)) / dt;
            else
                temp = (T2(:, :, k) - T0(:, :, k)) / (2 * dt);
            end
            J1(:, :, k, t) = (temp + v_del(phi, u(:, :, k, t), v(:, :, k, t), T(:, :, k, t), dphi, dlambda)) * cp;
            if size(sigma_accu, 4) == 1
                J2(:, :, k, t) = - sigma_accu(:, :, k) .* omega(:, :, k, t) / Ra * level(k) * cp;    
            else
                J2(:, :, k, t) = - sigma_accu(:, :, k, t) .* omega(:, :, k, t) / Ra * level(k) * cp;
            end    
            C1(:, :, k, t) = - kappa / level(k) * spherical_laplacian(J1(:, :, k, t), phi, dphi, dlambda);
            C2(:, :, k, t) = - kappa / level(k) * spherical_laplacian(J2(:, :, k, t), phi, dphi, dlambda);
            
        end

        T0 = T1;
        T1 = T2;    

    end




