function D = vorticity_tendency(zeta, level, event_timespan, f0, dt)

    D = zeros(size(zeta));
    for t = 1 : length(event_timespan)
    
        if t == 1
            zeta1 = zeta(:, :, :, t);
        end
        if t < length(event_timespan)
            zeta2 = zeta(:, :, :, t + 1);
        end
        temp = zeros(size(zeta1));
        if t == length(event_timespan)
            temp = f0 * d_dp((zeta1 - zeta0) / dt, level);
        elseif t == 1
            temp = f0 * d_dp((zeta2 - zeta1) / dt, level);
        else
            temp = f0 * d_dp((zeta2 - zeta0) / (2 * dt), level);
        end

        D(:, :, :, t) = temp;
        zeta0 = zeta1;
        zeta1 = zeta2;
        
    end

