function Le = e_fold_v1(omega_QG_500hPa_composite, event_omega_avg, dphi, dlambda, phi0, R)

    s = size(omega_QG_500hPa_composite);
    indx = (s(2) + 1) / 2;
    indy = (s(1) + 1) / 2;
    h_center = omega_QG_500hPa_composite(indy, indx);
    Delta = 6;

    h_x = @(x) interp1q((1-indx:indx-1)', omega_QG_500hPa_composite(indy, :)' - event_omega_avg(indy, :)' - exp(-1)*h_center, x);
    if(isnan(h_x(-Delta)))
        disp('Warning: omega field too small to find zero point in x')
        Lx = 0;
    else
        x1 = fzero(h_x, max(-Delta, 1-indx));
        x2 = fzero(h_x, min( Delta, indx-1));
        if x2 - x1 < 3 
            disp('Warning: roots too close in x direction!')
            Lx = abs((x2 + x1) / 2) * dlambda * cos(phi0) * R;
        else
            Lx = (x2 - x1) / 2 * dlambda * cos(phi0) * R;
        end
    end

    h_y = @(y) interp1q((1-indy:indy-1)', omega_QG_500hPa_composite(:, indx) - event_omega_avg(:, indx) - exp(-1)*h_center, y);
    if(isnan(h_y(-Delta)))
        disp('Warning: omega field too small to find zero point in y')
        Ly = 0;
    else
        y1 = fzero(h_y, max(-Delta, 1-indy));
        y2 = fzero(h_y, min( Delta, indy-1));
        if y2 - y1 < 3 
            disp('Warning: roots too close in y direction!')
            Ly = abs((y2 + y1) / 2) * dphi * R;
        else
            Ly = (y2 - y1) / 2 * dphi * R;
        end
    end

    Le = sqrt(Lx^2 + Ly^2) / (0.19 * 2 * pi);

end
