function [Le, Lx, Ly] = e_fold(omega_500hPa_composite, dphi, dlambda, phi0, R)

    s = size(omega_500hPa_composite);
    indx = (s(2) + 1) / 2;
    indy = (s(1) + 1) / 2;
    h_center = omega_500hPa_composite(indy, indx);

    h_x = @(x) interp1q((1-indx:indx-1)', omega_500hPa_composite(indy, :)' - exp(-1)*h_center, x);
    x1 = fzero(h_x, -6.0);
    x2 = fzero(h_x, 6.0);
    if x2 - x1 < 3 
        disp('Warning: roots too close in x direction!')
        Lx = abs((x2 + x1) / 2) * dlambda * cos(phi0) * R;
    else
        Lx = (x2 - x1) / 2 * dlambda * cos(phi0) * R;
    end

    h_y = @(y) interp1q((1-indy:indy-1)', omega_500hPa_composite(:, indx)  - exp(-1)*h_center, y);
    y1 = fzero(h_y, -6.0);
    y2 = fzero(h_y, 6.0);
    if y2 - y1 < 3 
        disp('Warning: roots too close in y direction!')
        Ly = abs((y2 + y1) / 2) * dphi * R;
    else
        Ly = (y2 - y1) / 2 * dphi * R;
    end

    Lx = Lx / 0.19;
    Ly = Ly / 0.19;

    Le = sqrt(Lx^2 + Ly^2) / (2 * pi);

end
