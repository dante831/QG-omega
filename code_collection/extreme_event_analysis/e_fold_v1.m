function Le = e_fold_v1(omega_QG_500hPa_composite, event_omega_avg, dphi, dlambda, phi0, R)
   
    nan_y = find(~isnan(omega_QG_500hPa_composite(:, (end+1)/2)));
    nan_x = find(~isnan(omega_QG_500hPa_composite((end+1)/2, :)));
    omega_QG_500hPa_composite = omega_QG_500hPa_composite(nan_y, nan_x);
    event_omega_avg = event_omega_avg(nan_y, nan_x);
    s = size(omega_QG_500hPa_composite);
    indx = (s(2) + 1) / 2;
    indy = (s(1) + 1) / 2;
    h_center = omega_QG_500hPa_composite(indy, indx);
    if strfind(pwd, 'CESM')
        Delta = 3.5;
    elseif strfind(pwd, 'GFDL')
        Delta = 2.0;
    end

    h_x = @(x) interp1q((1-indx:indx-1)', omega_QG_500hPa_composite(indy, :)' - event_omega_avg(indy, :)' - exp(-1)*h_center, x);
    x1 = fzero(h_x, max(-Delta, 1-indx));
    x2 = fzero(h_x, min( Delta, indx-1));
    if isnan(x1) && isnan(x2)
        disp('Warning: cannot find zero point in x')
        Lx = 0;
    elseif isnan(x1)
        Lx = abs(x2) / 2 * dphi * R;
    elseif isnan(x2)
        Lx = abs(x1) / 2 * dphi * R;
    else
        Lx = (abs(x1) + abs(x2)) / 2 * dlambda * cos(phi0) * R;
    end

    h_y = @(y) interp1q((1-indy:indy-1)', omega_QG_500hPa_composite(:, indx) - event_omega_avg(:, indx) - exp(-1)*h_center, y);
    y1 = fzero(h_y, max(-Delta, 1-indy));
    y2 = fzero(h_y, min( Delta, indy-1));
    if isnan(y1) && isnan(y2)
        disp('Warning: cannot find zero point in y')
        Ly = 0;
    elseif isnan(y1)
        Ly = abs(y2) / 2 * dphi * R;
    elseif isnan(y2)
        Ly = abs(y1) / 2 * dphi * R;
    else
        Ly = (abs(y1) + abs(y2)) / 2 * dphi * R;
    end

    if Lx == 0 && Ly ~= 0
        Lx = Ly;
    elseif Lx ~= 0 && Ly == 0
        Ly = Lx;
    elseif Lx == 0 && Ly == 0
        Lx = NaN;
        Ly = NaN;
    end
        
    Le = sqrt(Lx^2 + Ly^2) / (0.19 * 2 * pi);

end
